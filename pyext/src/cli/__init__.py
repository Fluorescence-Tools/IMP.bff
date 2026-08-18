"""Support for the ``imp_bff`` console script.

Not the whole command line. This package holds what ``bin/imp_bff`` imports --
the optimizer states that write frames during a flexible fit, and the angle-file
reader it starts from. The commands themselves (``flexfit``, ``rmsd``,
``auto_model``) live in ``bin/imp_bff``, because that is the file the console
script points at.

There is a second, unconnected command tree: ``python -m IMP.bff.cgdye.cli``
gives a ``dye`` group with its own subcommands, including ``rotamer``. It is not
reachable through ``imp_bff``. That split is a leftover of the PRD-107
migration, not a design; joining them is worth doing and has not been done.
"""

from IMP.bff.cli.flexfit import (  # noqa: F401
    WritePDBFrame,
    WriteRMFFrame,
    read_angle_file,
)

__all__ = ["WritePDBFrame", "WriteRMFFrame", "read_angle_file"]
