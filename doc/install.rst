************
Installation
************
``IMP.bff`` can either be installed from prebuilt binaries or from the source code.
For most users it is recommended to install prebuilt ``IMP.bff`` binaries.
Below the installation using prebuilt binaries and the prerequisites to compile
``IMP.bff`` are briefly outlined.

Prebuilt binaries
=================
It is recommended to install ``IMP.bff`` for Python environments via ``conda`` by.

.. code-block:: bash

    conda install -c tpeulen imp.bff


Compilation
===========
``IMP.bff`` can be compiled and installed using the source code provided in the
git repository. ``IMP.bff`` is compiled in Docker environments. The continuous
integration pipeline in the ``.gitlab-ci.yml`` reference to setup a
compilation machine. As ``IMP`` module ``IMP.bff`` requires in addition to
the ``IMP`` build dependencies the MongoDB C and the RTTR library.

Windows
-------
On windows ``IMP.bff`` is best compiled with the `Visual Studio 2017 <https://visualstudio.microsoft.com/>`_. For
compilation the Visual Studio Community edition is sufficient. In addition to
Visual Studio the libraries and the include files as listed above need to be
installed. The prebuilt binaries are compiled on Windows 10 with using 64bit anaconda
Python environments `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
using the conda build recipe that is provided with the source code in the ``conda-recipe``
folder.

macOS
-----
For MacOS the prebuilt binaries are compiled on MacOS 10.13 with the Apple clang
compiler using a anaconda distribution and the provided ``conda-recipe``.

Linux
-----
The Linux prebuilt binaries are compiled on Ubuntu 18.04 in an anaconda distribution
and the provided ``conda-recipe``.

Conda
-----
A conda recipe is provided in the folder 'conda-recipe' to build ``IMP.bff`` with the
`conda build <https://docs.conda.io/projects/conda-build/en/latest/>`_ environment.

To build the library download the recipe, install the conda build package and use
the provided recipe to build the library.

.. code-block:: bash

    conda install conda-build
    conda build conda-recipe

