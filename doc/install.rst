.. _installing-openmm:

Installing IMP.bff
******************

``IMP.bff`` is installed from prebuilt binaries using the Conda package manager (https://docs.conda.io).
Conda is included as part of the Anaconda Python distribution, which you can
download from https://docs.continuum.io/anaconda/install.  This is a Python
distribution specifically designed for scientific applications, with many of the
most popular mathematical and scientific packages preinstalled.  Alternatively
you can use Miniconda (available from https://docs.conda.io/en/latest/miniconda.html),
which includes only Python itself, plus the Conda package manager.  That offers
a much smaller initial download, with the ability to then install only the
packages you want.

\1. Begin by installing a recent 64 bit, Python 3.x version of either
Anaconda or Miniconda.

2. Open a command line terminal and type the following command
.. code-block:: bash

    conda install -c conda-forge imp.bff

This will install a version of IMP.bff compiled with the latest version of IMP.
Alternatively you can request a version that is compiled for a specific IMP
version with the command
.. code-block:: bash

    conda install -c tpeulen imp.bff imp=2.17

where :code:`2.17` should be replaced with the particular IMP version
you want to target.  We build packages for IMP 2.17 and above for Linux,
Windows, and macOS (x86).  Because different IMP releases are
not binary compatible with each other, IMP.bff can only work with the particular
IMP version it was compiled with.

.. note::

    Future IMP.bff releases may be distributed together with IMP.


3. Optionally verify your installation using the unit-test that ship with the IMP.bff
git repository by typing the following commands:
.. code-block:: bash

    conda install nose git
    git clone https://gitlab.peulen.xyz/imp/bff
    cd bff/test
    nosetests

This command verifies that IMP.bff is installed and that all platforms produce consistent
results.

Compilation
-----------
A second option is to compile IMP.bff from source. This provides more flexibility,
but it is much more work, and there is rarely a need for anyone but advanced users
to compile from source. ``IMP.bff`` can be compiled and installed using the source
code provided in the git repository. ``IMP.bff`` is compiled in Docker environments.
The continuous integration pipeline in ``.gitlab-ci.yml`` can serve as a reference
for compilation. As an IMP module IMP.bff can be compiled either together with IMP
or out of tree (https://integrativemodeling.org/2.17.0/doc/manual/outoftree.html).

On windows ``IMP.bff`` is compiled with the `Visual Studio 2017 <https://visualstudio.microsoft.com/>`_. For
compilation the Visual Studio Community edition is sufficient. In addition to
Visual Studio the libraries and the include files as listed above need to be
installed. The prebuilt binaries are compiled on Windows 10 with using 64bit anaconda
Python environments `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
using the conda build recipe that is provided with the source code in the ``conda-recipe``
folder.

For MacOS the prebuilt binaries are compiled on MacOS 10.14 with the Apple clang
compiler using a anaconda distribution and the provided ``conda-recipe``.

The Linux prebuilt binaries are compiled on Ubuntu 18.04 in an anaconda distribution
and the provided ``conda-recipe``.
