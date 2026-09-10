/*
 * Where a probe sits on a structure, and how it gets there -- over
 * IMP::atom::Hierarchy, so connection layer. `ProbeAttachment.h` carries the
 * frame, the site, the strip and the attachment -- in flrCIF's vocabulary:
 * **probe**, which the dictionary uses throughout (`_flr_probe_list`,
 * `_flr_poly_probe_position`, `_flr_sample_probe_details`) and which does not
 * collide with mmCIF's `label_asym_id` / `label_seq_id` naming. From a
 * structure's point of view a fluorophore, a spin label and a quencher are
 * one thing: a molecule on a linker, anchored at a residue's backbone.
 *
 * **One backbone frame** serves both halves of the package: `backbone_rotation`
 * in `Rotamer.h` (core), which places rotamer libraries, is the flat-matrix
 * view of `backbone_frame`, which places explicit probes.
 *
 * **Stripping is in place.** A default that clones the structure on every call
 * hides the cost of doing so; a caller that wants the original clones it and
 * can see what that costs. This file holds no Python.
 */

// `StripSelection` is declared with the AV builder, which reads the same
// masks; declaring it twice defines SWIG's type test twice.
IMP_SWIG_VALUE(IMP::bff, ProbePosition, ProbePositions);
IMP_SWIG_VALUE(IMP::bff, ProbeAttachment, ProbeAttachments);

%feature("kwargs") IMP::bff::attach_probes;
%feature("kwargs") IMP::bff::strip_sidechain_at_site;
%feature("kwargs") IMP::bff::ProbePosition::ProbePosition;
%feature("kwargs") IMP::bff::probe_position_from_source_info;

%include "IMP/bff/ProbeAttachment.h"

%template(ProbeAttachmentList) std::vector<IMP::bff::ProbeAttachment>;

%attribute_py(IMP::bff::ProbeAttachment, Hierarchy, probe, get_probe);
%attributestring(IMP::bff::ProbeAttachment, std::string, chain, get_chain);
%attribute(IMP::bff::ProbeAttachment, int, residue, get_residue);
%attribute(IMP::bff::ProbeAttachment, int, n_stripped, get_n_stripped);

%attributestring(IMP::bff::ProbePosition, std::string, key, get_key);
%attributestring(IMP::bff::ProbePosition, std::string, source_info, get_source_info);
