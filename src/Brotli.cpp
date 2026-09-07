/**
 *  \file Brotli.cpp
 *  \brief The vendored brotli codec, as one translation unit.
 *
 *  `src/brotli/` holds google/brotli's decoder (MIT, see its `LICENSE`)
 *  kept in its upstream shape: the `common/` and `dec/` implementation files
 *  carry the `.inc` suffix and are pulled in here, because IMP compiles a
 *  module from the `.cpp` files sitting directly in `src/` -- `bff_all.cpp`
 *  globs that one level and no deeper -- and because IMP's second,
 *  per-file compilation mode would otherwise define every brotli symbol
 *  twice. The only edits to upstream are the renamed suffix
 *  and `<brotli/x.h>` rewritten to the relative header next to it -- so the
 *  vendored headers are used, never a brotli that happens to sit in the
 *  build environment. `DrotReader.cpp` is the sole consumer.
 *
 *  The decoder is compiled hidden: nothing named `Brotli*` leaves
 *  `libimp_bff`, so a process that also loads a real libbrotlidec (through
 *  fonttools, pillow, ...) keeps the two apart.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#pragma GCC visibility push(hidden)

#include "brotli/common/constants.inc"
#include "brotli/common/context.inc"
#include "brotli/common/dictionary.inc"
#include "brotli/common/platform.inc"
#include "brotli/common/shared_dictionary.inc"
#include "brotli/common/transform.inc"

#include "brotli/dec/bit_reader.inc"
#include "brotli/dec/decode.inc"
#include "brotli/dec/huffman.inc"
#include "brotli/dec/prefix.inc"
#include "brotli/dec/state.inc"
#include "brotli/dec/static_init.inc"


#include "brotli/enc/backward_references.inc"
#include "brotli/enc/backward_references_hq.inc"
#include "brotli/enc/bit_cost.inc"
#include "brotli/enc/block_splitter.inc"
#include "brotli/enc/brotli_bit_stream.inc"
#include "brotli/enc/cluster.inc"
#include "brotli/enc/command.inc"
#include "brotli/enc/compound_dictionary.inc"
#include "brotli/enc/compress_fragment.inc"
/* [imp.bff] Upstream compiles the two fragment compressors apart, so both
 * may define these file-local helpers (and their own MIN_RATIO). In one
 * translation unit the second copy is renamed here rather than in the
 * vendored source. */
#undef MIN_RATIO
#define BrotliStoreMetaBlockHeader BrotliStoreMetaBlockHeader_two_pass
#define RewindBitPosition RewindBitPosition_two_pass
#include "brotli/enc/compress_fragment_two_pass.inc"
#undef RewindBitPosition
#undef BrotliStoreMetaBlockHeader
#include "brotli/enc/dictionary_hash.inc"
#include "brotli/enc/encode.inc"
#include "brotli/enc/encoder_dict.inc"
/* [imp.bff] same story: brotli_bit_stream.inc has its own SortHuffmanTree. */
#define SortHuffmanTree SortHuffmanTree_entropy_encode
#include "brotli/enc/entropy_encode.inc"
#undef SortHuffmanTree
#include "brotli/enc/fast_log.inc"
#include "brotli/enc/histogram.inc"
#include "brotli/enc/literal_cost.inc"
#include "brotli/enc/memory.inc"
#include "brotli/enc/metablock.inc"
#include "brotli/enc/static_dict.inc"
#include "brotli/enc/static_dict_lut.inc"
#include "brotli/enc/static_init.inc"
#include "brotli/enc/utf8_util.inc"

#pragma GCC visibility pop
