#!/usr/bin/env python2

from distutils.core import setup

setup(name = 'MDAextensions',
      version = '0.1',
      packages = [
        'MDAextensions',
        'MDAextensions.analysis',
        'MDAextensions.datatools',
        'MDAextensions.scripts',
        'MDAextensions.tests'
      ],
      install_requires = [
        'numpy',
        'scipy',
        'MDAnalysis'
      ],
      )
