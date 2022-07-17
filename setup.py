#! /usr/bin/env python
import os
import sys
import platform
import subprocess
import multiprocessing

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


def read_version(header_file):
    version = "0.0.0"
    with open(header_file, "r") as fp:
        for line in fp.readlines():
            if "#define" in line and "IMPBFF_VERSION" in line:
                version = line.split()[-1]
    return version.replace('"', '')


VERSION = read_version(os.path.dirname(os.path.abspath(__file__)) + '/include/info.h')
LICENSE = 'Mozilla Public License 2.0 (MPL 2.0)'
NAME = "IMP.bff"  # name of the module
DESCRIPTION = "Python bindings for IMP.bff"
LONG_DESCRIPTION = """"IMP.bff is a C++ library with Python wrapper to construct \
networks of computation node with associated ports that can be deposited in a \
mongoDB. IMP.bff is based on chinet (the data management backend of chisurf)."""


class CMakeExtension(Extension):

    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


def build_swig_documentation():
    # build the documentation.i file using doxygen and doxy2swig
    if not os.path.isfile("./pyext/documentation.i"):
        print("-- building documentation.i using doxygen and doxy2swig")
        path = os.path.dirname(os.path.abspath(__file__)) + "/doc"
        env = os.environ.copy()
        subprocess.check_call(["doxygen"], cwd=path, env=env)
        subprocess.check_call(
            [
                "python", "../tools/doxy2swig.py", "./_build/xml/index.xml", "../pyext/documentation.i"],
            cwd=path,
            env=env
        )


class CMakeBuild(build_ext):

    def run(self):
        build_swig_documentation()
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        print(NAME, " VERSION:", VERSION)
        extdir = os.path.abspath(
            os.path.dirname(
                self.get_ext_fullpath(ext.name)
            )
        ).replace('\\', '/')

        cmake_args = [
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DCMAKE_SWIG_OUTDIR=' + extdir
        ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        if platform.system() == "Windows":
            cmake_args += [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                '-GVisual Studio 14 2015 Win64'
            ]
        else:
            build_args += [
                '--',
                '-j%s' % int(multiprocessing.cpu_count() * 1.5)
            ]
        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        print("BUILDING: " + " ".join(cmake_args))
        subprocess.check_call(
            ['cmake', ext.sourcedir] + cmake_args,
            cwd=self.build_temp,
            env=env
        )
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
        subprocess.check_call(['make install'] + build_args, cwd=self.build_temp)

setup(
    name=NAME,
    version=VERSION,
    license=LICENSE,
    author='Thomas-Otavio Peulen',
    author_email='thomas@peulen.xyz',
    ext_modules=[
        CMakeExtension('IMP.bff')
    ],
    cmdclass={
        'build_ext': CMakeBuild
    },
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    install_requires=[
        'numpy',
        'numba',
        'scipy',
        'click',
        'pandas',
        'tqdm'
    ],
    setup_requires=[
        'setuptools'
    ],
    zip_safe=False,
    classifiers=[
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ]
)
