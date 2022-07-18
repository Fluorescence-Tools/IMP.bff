import typing
import json
import pathlib

import IMP
import IMP.core
import IMP.atom
import IMP.pmi.restraints
import IMP.rmf
import IMP.bff

import numpy as np


class XLinkScore(typing.TypedDict):
    """
    Attributes
    ----------

    individual: list
        The scores to all potential cross linking partners. The score is the
        fraction of distances of all path that are shorter or equal to the
        cross-linker length

    total: np.ndarray
        The sum of scores to all potential cross linking partners. The score is the
        fraction of distances of all path that are shorter or equal to the
        cross-linker length

    """
    individual: list
    total: np.ndarray


class ScoreXlinkSurfaceDistance(IMP.pmi.restraints.RestraintBase):
    """
    :param linker_length:
    :param linker_width:
    :param radius:
    :param simulation_grid_spacing:
    :param min_points: the minimum number of points the the computed volume. If
        there are less points the site is not accessible and the score will be
        zero
    :param verbose:
    :return:
    """

    model: IMP.Model
    obstacles: np.ndarray
    xlinks: typing.Dict[int, typing.Dict]
    linker_length: float = 20.0
    linker_width: float = 2.0
    radius: float = 1.0
    simulation_grid_spacing: float = 3.5
    min_points: int = 100
    verbose: bool = False
    name_map: dict = None

    def compute_scores(self, **kwargs) -> XLinkScore:
        linker_length = kwargs.get('linker_length', self.linker_length)
        linker_width = kwargs.get('linker_width', self.linker_width)
        radius = kwargs.get('radius', self.radius)
        simulation_grid_spacing = kwargs.get('simulation_grid_spacing', self.simulation_grid_spacing)
        verbose = kwargs.get('verbose', self.verbose)
        min_points = kwargs.get('min_points', self.min_points)
        obstacles = self.obstacles
        scores = list()
        total_scores = list()
        for xlink_key in self.xlinks:
            xlink = self.xlinks[xlink_key]
            protein_1 = xlink['protein_1']
            protein_2 = xlink['protein_2']
            residue_1 = xlink['residue_1']
            residue_2 = xlink['residue_2']
            attachment_idx_1 = np.where((obstacles['res_id'] == residue_1) & (obstacles['protein_name'] == protein_1))[0]
            attachment_idx_2 = np.where((obstacles['res_id'] == residue_2) & (obstacles['protein_name'] == protein_2))[0]
            x_links_scores = []
            for idx_1 in attachment_idx_1:
                origin = obstacles['xyz'][idx_1]
                for idx_2 in attachment_idx_2:
                    target = obstacles['xyz'][idx_2]
                    if np.linalg.norm(origin - target) > linker_length:
                        if verbose:
                            print("Eucledian distance > linker length")
                        x_links_scores.append(0.0)
                    else:
                        av = get_path_length(
                            obstacles=self.obstacles,
                            idx=idx_1,
                            linker_length=linker_length,
                            linker_width=linker_width,
                            radius=radius,
                            simulation_grid_spacing=simulation_grid_spacing
                        )
                        if len(np.array(av.points()).T) < min_points:
                            if verbose:
                                print(protein_1, residue_1, "is not accessible")
                            x_links_scores.append(0.0)
                        else:
                            if verbose:
                                print("Eucledian distance < linker length.. Testing surface distance")
                            points = av.points()
                            xyz = points[0:3]
                            dist = points[3]
                            dist_eq = np.linalg.norm(xyz.T - target, axis=1) + dist
                            shorter = np.sum(dist_eq < (linker_length + radius))
                            score = shorter / len(dist_eq)
                            x_links_scores.append(score)
            indiviudal_scores = np.array(x_links_scores)
            scores.append(indiviudal_scores)
            total_score = indiviudal_scores.sum()
            total_scores.append(total_score)
            if verbose:
                print(xlink_key, ":", xlink, "score:", total_score)
        re: XLinkScore = {
            'total': np.array(total_scores, dtype=np.float),
            'individual': scores
        }
        return re

    def __init__(
            self,
            root_hier: IMP.atom.Hierarchy = None,
            verbose: bool = False,
            xlink_settings_file: str = '',
            rmf_file: str = '',
            weight: float = 1.0,
            label: str = 'XLinkSurfaceDistance'
    ):
        # create a new ScoreXlinkRMF object with predefined settings
        base_dir = pathlib.Path(xlink_settings_file).parent
        with open(xlink_settings_file, 'r') as fp:
            xlink_settings = json.load(fp)
        if root_hier is None:
            model = IMP.Model()
        else:
            model = root_hier.get_model()
        super().__init__(
            model,
            weight=weight,
            label=label
        )
        self.model = model
        keys, formats = list(zip(*OBSTACLES_KEYS_FORMATS))
        obstacles = np.zeros(0, dtype={
                'names': keys,
                'formats': formats
            }
        )
        self.obstacles = obstacles
        self.xlinks = dict()
        self.linker_length = xlink_settings['linker_length']
        self.linker_width = xlink_settings['linker_width']
        self.radius = xlink_settings['radius']
        self.simulation_grid_spacing = xlink_settings['simulation_grid_spacing']
        self.min_points = xlink_settings['min_points']
        self.verbose = verbose
        xlink_file = xlink_settings['xlink_file']
        if pathlib.Path(xlink_file).is_file():
            self.xlinks = IMP.bff.tools.read_xlink_table(
                fn=xlink_file
            )
        elif (base_dir / xlink_file).is_file():
            self.xlinks = IMP.bff.tools.read_xlink_table(
                fn=str(base_dir / xlink_file)
            )
        else:
            raise FileNotFoundError("Could not find XL file %s" % xlink_file)
        if pathlib.Path(rmf_file).is_file():
            self.obstacles = get_obstacles(
                rmf_file=rmf_file,
                model=self.model,
                name_map=self.name_map
            )
        # Add custom metadata (will be saved in RMF output)
        self.rs.filename = xlink_settings_file

    def __call__(
            self,
            rmf_file: str = '',
            xlink_file: str = '',
            frame_index: int = 0,
            hier: IMP.atom.Hierarchy = None,
            *args, **kwargs
    ):
        if pathlib.Path(xlink_file).is_file():
            self.xlinks = IMP.bff.tools.read_xlink_table(
                fn=xlink_file
            )
        if pathlib.Path(rmf_file).is_file():
            self.obstacles = get_obstacles(
                rmf_file=rmf_file,
                hier=hier,
                model=self.model,
                frame_index=frame_index,
                name_map=self.name_map
            )
        return self.compute_scores(**kwargs)

    def get_score(self) -> float:
        total_score = float(np.sum(self()['total']))
        return total_score


