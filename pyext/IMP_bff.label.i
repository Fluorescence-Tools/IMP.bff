/*
 * Where a probe sits on a structure, and how it gets there.
 * This file holds no Python.
 *
 * `ProbeAttachment.h` carries the frame, the site, the strip and the attachment --
 * in flrCIF's vocabulary: **probe**, which the dictionary uses throughout
 * (`_flr_probe_list`, `_flr_poly_probe_position`, `_flr_sample_probe_details`)
 * and which does not collide with mmCIF's `label_asym_id` / `label_seq_id`
 * naming -- the very items a position is made of. From a structure's point of
 * view a fluorophore, a spin label and a quencher are one thing: a molecule
 * on a linker, anchored at a residue's backbone. What differs is what is
 * measured (spectra for FRET, a dipole for DEER), and that lives with the
 * measurement.
 *
 * `SequenceAlignment.h` carries the fluorescent-protein detection and the
 * pLDDT segmentation, beside the Smith-Waterman they are built on;
 * `Quenching.h` has carried the quencher and its parameters since an
 * earlier batch.
 *
 * **One backbone frame** serves both halves of the package: `backbone_rotation`
 * in `Rotamer.h`, which places rotamer libraries, is the flat-matrix view
 * of `backbone_frame`, which places explicit probes. One atom-name accessor
 * too: `atom_name` in `HierarchyFrame.h`.
 *
 * **Stripping is in place.** A default that clones the structure on every call
 * hides the cost of doing so; a caller that wants the original clones it and
 * can see what that costs.
 */

// `StripSelection` is declared with the AV builder, which reads the same
// masks; declaring it twice defines SWIG's type test twice.
IMP_SWIG_VALUE(IMP::bff, ProbePosition, ProbePositions);
IMP_SWIG_VALUE(IMP::bff, ProbeAttachment, ProbeAttachments);
IMP_SWIG_VALUE(IMP::bff, SequenceSegment, SequenceSegments);

%feature("kwargs") IMP::bff::attach_probes;
%feature("kwargs") IMP::bff::strip_sidechain_at_site;
%feature("kwargs") IMP::bff::ProbePosition::ProbePosition;
%feature("kwargs") IMP::bff::probe_position_from_source_info;
%feature("kwargs") IMP::bff::get_fp_domains;
%feature("kwargs") IMP::bff::segments_from_plddt;
%feature("kwargs") IMP::bff::parse_plddt_from_pdb;

%include "IMP/bff/ProbeAttachment.h"

%template(ProbeAttachmentList) std::vector<IMP::bff::ProbeAttachment>;
%template(SequenceSegmentList) std::vector<IMP::bff::SequenceSegment>;
%template(MapIntDouble) std::map<int, double>;

%attribute_py(IMP::bff::ProbeAttachment, Hierarchy, probe, get_probe);
%attributestring(IMP::bff::ProbeAttachment, std::string, chain, get_chain);
%attribute(IMP::bff::ProbeAttachment, int, residue, get_residue);
%attribute(IMP::bff::ProbeAttachment, int, n_stripped, get_n_stripped);

%attributestring(IMP::bff::ProbePosition, std::string, key, get_key);
%attributestring(IMP::bff::ProbePosition, std::string, source_info, get_source_info);
