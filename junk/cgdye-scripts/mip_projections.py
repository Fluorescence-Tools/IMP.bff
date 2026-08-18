#!/usr/bin/env python
"""Integrate occupancy MRC density maps and save as translucent RGBA PNGs.

For each region produces two transparent PNGs (stackable in PowerPoint):
  proj_top_<region>.png   – density summed along the z-axis  (xy-plane, top view)
  proj_side_<region>.png  – density summed along the x-axis  (zy-plane, side view)

The density maps are centered on a fixed global grid to ensure all regions align.
"""
import struct
import sys
from pathlib import Path

import click
import numpy as np

REGION_ORDER = ("linker", "top", "middle", "bottom")
REGION_HEX = {
    "linker":  "#e07b39",
    "top":     "#3b9de0",
    "middle":  "#5bc46a",
    "bottom":  "#b84fd8",
}

def hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

def load_mrc_binary(path: Path):
    """Load MRC file using numpy binary parsing for speed and axis control.
    
    Returns (data, voxel_size, origin) where data is (nz, ny, nx).
    """
    with open(path, "rb") as f:
        header = f.read(1024)
        # Header: nx, ny, nz, mode, nxstart, nystart, nzstart, mx, my, mz
        # Indices 0, 1, 2 are nx, ny, nz
        nx, ny, nz = struct.unpack("3i", header[:12])
        # Origin is at byte 196 (3 floats)
        origin = struct.unpack("3f", header[196:208])
        # Spacing: xlen/mx, ylen/my, zlen/mz. 
        # mrc header floats at offset 40: xlen, ylen, zlen
        xlen, ylen, zlen = struct.unpack("3f", header[40:52])
        # Assuming isotropic or at least mapping and mx==nx etc.
        vs_x, vs_y, vs_z = xlen/nx, ylen/ny, zlen/nz
        
        # Read the data (floats, mode 2)
        f.seek(1024)
        count = nx * ny * nz
        data = np.fromfile(f, dtype=np.float32, count=count)
        
    data = data.reshape((nz, ny, nx))
    data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    data = np.clip(data, 0.0, None)
    return data, vs_x, np.array(origin)

def embed_in_global_grid(data, vs, origin):
    """Embed local MRC data into a fixed global grid centered at (0,0,0)."""
    nz, ny, nx = data.shape
    
    # Standardize global grid: [-50, 50] x [-50, 50] x [-60, 60]
    g_xmin, g_xmax = -50.0, 50.0
    g_ymin, g_ymax = -50.0, 50.0
    g_zmin, g_zmax = -60.0, 60.0
    
    gnx = int(round((g_xmax - g_xmin) / vs))
    gny = int(round((g_ymax - g_ymin) / vs))
    gnz = int(round((g_zmax - g_zmin) / vs))
    
    global_data = np.zeros((gnz, gny, gnx), dtype=np.float32)
    
    # Local data indices to global grid indices mapping
    # Local origin maps to these global indices:
    ix_g = int(round((origin[0] - g_xmin) / vs))
    iy_g = int(round((origin[1] - g_ymin) / vs))
    iz_g = int(round((origin[2] - g_zmin) / vs))
    
    # Calculate slice intersection
    lz_start = max(0, -iz_g)
    lz_end   = min(nz, gnz - iz_g)
    gz_start = iz_g + lz_start
    gz_end   = iz_g + lz_end
    
    ly_start = max(0, -iy_g)
    ly_end   = min(ny, gny - iy_g)
    gy_start = iy_g + ly_start
    gy_end   = iy_g + ly_end
    
    lx_start = max(0, -ix_g)
    lx_end   = min(nx, gnx - ix_g)
    gx_start = ix_g + lx_start
    gx_end   = ix_g + lx_end
    
    if (gz_start < gz_end) and (gy_start < gy_end) and (gx_start < gx_end):
        global_data[gz_start:gz_end, gy_start:gy_end, gx_start:gx_end] = \
            data[lz_start:lz_end, ly_start:ly_end, lx_start:lx_end]
            
    return global_data

def to_rgba(proj2d: np.ndarray, rgb: tuple, gamma: float) -> np.ndarray:
    """Normalise projection to [0,1], apply gamma, return (H,W,4) uint8 RGBA."""
    vmax = proj2d.max()
    if vmax == 0:
        return np.zeros((*proj2d.shape, 4), dtype=np.uint8)
    norm = np.clip(proj2d / vmax, 0.0, 1.0)
    alpha = np.power(norm, gamma)
    rgba = np.zeros((*proj2d.shape, 4), dtype=np.uint8)
    rgba[..., 0] = rgb[0]
    rgba[..., 1] = rgb[1]
    rgba[..., 2] = rgb[2]
    rgba[..., 3] = np.round(alpha * 255).astype(np.uint8)
    return rgba

def _write_png(rgba: np.ndarray, path: Path) -> None:
    """Write an 8-bit RGBA array as a PNG, using only the standard library.

    PNG is a short format when the image is already RGBA and unfiltered: an
    8-byte signature, IHDR, one zlib stream of scanlines each prefixed with
    filter byte 0, and IEND. Writing it here keeps Pillow out of IMP.bff for
    the sake of one nearest-neighbour upscale.
    """
    import struct
    import zlib

    height, width = rgba.shape[:2]
    raw = b"".join(
        b"\x00" + rgba[row].astype(np.uint8).tobytes() for row in range(height)
    )

    def chunk(tag: bytes, payload: bytes) -> bytes:
        body = tag + payload
        return struct.pack(">I", len(payload)) + body + struct.pack(">I", zlib.crc32(body))

    png = (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )
    path.write_bytes(png)


def save_png(rgba: np.ndarray, path: Path, scale: int):
    H, W = rgba.shape[:2]
    # Nearest-neighbour upscale: repeating rows and columns is exactly what
    # Image.NEAREST did, without the dependency.
    upscaled = np.repeat(np.repeat(rgba, scale, axis=0), scale, axis=1)
    _write_png(upscaled, path)
    print(f"  {path.name}  ({W*scale}×{H*scale} px)")

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--analysis-dir", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--regions", default=",".join(REGION_ORDER), show_default=True)
@click.option("--gamma", type=float, default=0.5, show_default=True)
@click.option("--scale", type=int, default=8, show_default=True)
@click.option("--output-dir", type=click.Path(path_type=Path), default=None)
@click.option("--mode", type=click.Choice(["sum", "max"]), default="sum", show_default=True)
def main(analysis_dir, regions, gamma, scale, output_dir, mode):
    """Integrate MRC density maps to transparent per-region PNGs."""
    region_list = [r.strip() for r in regions.split(",")]
    out_dir = output_dir or analysis_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Source : {analysis_dir}")
    print(f"Output : {out_dir}")
    print(f"Mode   : {mode}\n")

    project_fn = np.sum if mode == "sum" else np.max

    for region in region_list:
        mrc = analysis_dir / f"occupancy_{region}.mrc"
        if not mrc.exists():
            continue
        print(f"[{region}]", flush=True)

        data, vs, origin = load_mrc_binary(mrc)
        # Embed in centered global grid
        g_data = embed_in_global_grid(data, vs, origin)
        rgb  = hex_to_rgb(REGION_HEX.get(region, "#ffffff"))

        # Top view: project along z (axis 0) -> (ny, nx), flip y for traditional viewing
        top = project_fn(g_data, axis=0)[::-1, :]
        save_png(to_rgba(top, rgb, gamma), out_dir / f"proj_top_{region}.png", scale)

        # Side view: project along x (axis 2) -> (nz, ny), flip z for z-up
        side = project_fn(g_data, axis=2)[::-1, :]
        save_png(to_rgba(side, rgb, gamma), out_dir / f"proj_side_{region}.png", scale)

    print("\nDone.")

if __name__ == "__main__":
    main()
