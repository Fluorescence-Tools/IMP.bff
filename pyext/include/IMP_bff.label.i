/*
 * The sequence side of labelling: fluorescent-protein detection and pLDDT
 * segmentation (`SequenceAlignment.h`, wrapped earlier in core.i beside the
 * Smith-Waterman they are built on) -- the value type and the keyword
 * arguments. The attachment itself (`ProbeAttachment.h`) is the connection
 * layer's, in probeattachment.i; `PhotophysicsQuenching.h` has carried the quencher and
 * its parameters since an earlier batch. This file holds no Python.
 */

IMP_SWIG_VALUE(IMP::bff, SequenceSegment, SequenceSegments);
%feature("kwargs") IMP::bff::get_fp_domains;
%feature("kwargs") IMP::bff::segments_from_plddt;
%feature("kwargs") IMP::bff::parse_plddt_from_pdb;
%template(SequenceSegmentList) std::vector<IMP::bff::SequenceSegment>;
%template(MapIntDouble) std::map<int, double>;
