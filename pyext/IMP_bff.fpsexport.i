/*
 * FPS's export formats.
 *
 * Everything here is C++: the file layout, the number formatting and the
 * geometry. A program in `bin/` picks the files and the flags; nothing about
 * *what a line looks like* is decided in Python, because the byte layout is
 * the interface a decade of downstream readers depend on.
 *
 * The same `std::vector` member trap as `IMP_bff.docking.i` applies to
 * `FPSResultRow::pairs`, `FPSResultTable::rows` and `FPSExportOptions::labels`:
 * bind the owning object to a name before reading a member off it.
 */

IMP_SWIG_VALUE(IMP::bff, FPSMolecule, FPSMolecules);
IMP_SWIG_VALUE(IMP::bff, FPSLabelPosition, FPSLabelPositions);
IMP_SWIG_VALUE(IMP::bff, FPSResultRow, FPSResultRows);
IMP_SWIG_VALUE(IMP::bff, FPSResultTable, FPSResultTables);
IMP_SWIG_VALUE(IMP::bff, FPSExportOptions, FPSExportOptionsList);

%feature("kwargs") IMP::bff::fps_result_row;
%feature("kwargs") IMP::bff::fps_bootstrap_table;
%feature("kwargs") IMP::bff::add_rmsd_columns;
%feature("kwargs") IMP::bff::add_best_fit;
%feature("kwargs") IMP::bff::write_fps_pymol_script;
%feature("kwargs") IMP::bff::write_fps_pymol_scripts;
%feature("kwargs") IMP::bff::write_fps_overlay;
%feature("kwargs") IMP::bff::write_fps_overlay_states;
%feature("kwargs") IMP::bff::write_fps_r_table;
%feature("kwargs") IMP::bff::write_fps_chi2_table;
%feature("kwargs") IMP::bff::write_fps_exports;
%feature("kwargs") IMP::bff::write_fps_screening_r_table;
%feature("kwargs") IMP::bff::write_fps_screening_chi2_table;

%include "IMP/bff/FPSExport.h"

%template(FPSMoleculeList) std::vector<IMP::bff::FPSMolecule>;
%template(FPSLabelPositionList) std::vector<IMP::bff::FPSLabelPosition>;
%template(FPSResultRowList) std::vector<IMP::bff::FPSResultRow>;
