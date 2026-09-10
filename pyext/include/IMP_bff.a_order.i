/* [bff] The topic files, in the order the two true entry files impose:
 * core.i's internal order for the core topics (pathmap before states
 * before av-model; VdwRadii before av; types before everything), then
 * layer.i's order for the connection layer. IMP's module tooling
 * %includes this directory's files alphabetically, but SWIG includes a
 * file once -- the later alphabetical %includes are no-ops. The
 * standalone entry reaches the same order through core.i and the
 * IMPBFF_WITH_IMP block of IMP_bff_standalone.i. */

%include "IMP_bff.core.i"
%include "IMP_bff.layer.i"
