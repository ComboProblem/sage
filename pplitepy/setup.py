#!/usr/bin/env python

import os
import sys

from setuptools import setup, Command
from setuptools.extension import Extension

# NOTE: setuptools build_ext does not work properly with Cython code
from distutils.command.build_ext import build_ext as _build_ext

# Adapted from Cython's new_build_ext
class build_ext(_build_ext):
    def run(self):
        # Check dependencies
        try:
            from Cython.Build.Dependencies import cythonize
        except ImportError as E:
            sys.stderr.write("Error: {0}\n".format(E))
            sys.stderr.write("The installation of ppl requires Cython\n")
            sys.exit(1)

        try:
            # We need the header files for cysignals at compile-time
            import cysignals
        except ImportError as E:
            sys.stderr.write("Error: {0}\n".format(E))
            sys.stderr.write("The installation of ppl requires cysignals\n")
            sys.exit(1)

        try:
            # We need the header files for gmpy2 at compile-time
            import gmpy2
        except ImportError as E:
            sys.stderr.write("Error: {0}\n".format(E))
            sys.stderr.write("The installation of ppl requires gmpy2\n")
            sys.exit(1)

        self.extensions[:] = cythonize(
            self.extensions,
            include_path=sys.path,
            compiler_directives={'embedsignature': True,
                                 'language_level': '3'})

        _build_ext.run(self)

class TestCommand(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess, os, tempfile, shutil

        old_path = os.getcwd()
        tempdir_path = tempfile.mkdtemp()
        try:
            shutil.copytree('./tests', tempdir_path, dirs_exist_ok=True)
            os.chdir(tempdir_path)

            if subprocess.call([sys.executable, 'runtests.py']):
                raise SystemExit("Doctest failures")

            if subprocess.call([sys.executable, 'setup.py', 'build_ext', '--inplace']) or \
                    subprocess.call([sys.executable, '-c', "import testpplpy; testpplpy.test(); testpplpy.example()"]):
                raise SystemExit("Cython test 1 failure")

            if subprocess.call([sys.executable, 'setup2.py', 'build_ext', '--inplace']) or \
                    subprocess.call([sys.executable, '-c', "import testpplpy2; testpplpy2.test(); testpplpy2.example()"]):
                raise SystemExit("Cython test 2 failure")
        finally:
            os.chdir(old_path)
            shutil.rmtree(tempdir_path)

extensions = [
    Extension('pplite.integer_conversions', sources=['pplite/integer_conversions.pyx']),
    Extension('pplite.linear_algebra', sources=['pplite/linear_algebra.pyx']),
    Extension('pplite.constraint', sources=['pplite/constraint.pyx']),
    Extension('pplite.generators', sources=['pplite/generators.pyx']),
    Extension('pplite.intervals', sources=['pplite/intervals.pyx'])
    ]
#finish writing setup at some point
setup(
    ext_modules = extensions,
    cmdclass = {'build_ext': build_ext, 'test': TestCommand},
)