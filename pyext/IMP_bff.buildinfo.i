/*
 * Which IMP.bff this is, for the IMP module build. The standalone build
 * defines the same name in its own entry (standalone/pyext) -- its answer
 * is "core"/"core+imp"; here, where IMP's own Python is importable beside
 * the module, the answer is "imp".
 */
%pythoncode %{
def get_build():
    """Which IMP.bff this is: the IMP module build answers "imp"."""
    return "imp"
%}
