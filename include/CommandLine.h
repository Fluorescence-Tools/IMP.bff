#ifndef IMPBFF_COMMANDLINE_H
#define IMPBFF_COMMANDLINE_H

#include <string>
#include <utility>
#include <vector>

#include "IMPCompatibility.h"

IMPBFF_BEGIN_NAMESPACE

//! The command line, compiled: one dispatcher, one function per sub.
/*!
    The programs in `bin/` were Python because a click command is (the
    repository's language rule). The owner's ruling of 2026-09-09 moves the
    command line itself into the library: `bin/` scripts that need nothing
    beyond the core become \c subs of one compiled dispatcher, so they ship
    with the wheel and run with no IMP at all -- the wheel's \c imp_bff
    console script is a Python shell of a few lines around \c command_line_main.
    Scripts whose subject is the heavy stack (the IMP/pmi/RMF modelling
    commands of `bin/imp_bff`) stay Python; no dispatcher makes them lighter.

    \code
    IMP.bff.command_line_main(["traj2drot", "traj.bcif", "lib.drot.pto"])
    \endcode

    \param[in] args the words after the program name; \c args[0] names the
               sub
    \return the process exit code: 0 success, 1 failure, 2 bad usage
*/
IMPBFFEXPORT int command_line_main(const std::vector<std::string>& args);

//! The subs and their one-line briefs, in dispatch order.
/*! What \c imp_bff help prints, and what a caller offers in its own UI. */
IMPBFFEXPORT std::vector<std::pair<std::string, std::string> > command_line_subcommands();

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_COMMANDLINE_H
