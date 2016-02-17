#!/usr/bin/env python

import sys

from setuptools import find_packages, setup

# get requirements from `requirements.txt`
with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(name='hori',
      description='hori',
      author='SUNCAT Center for Interface Science and Catalysis',
      url='git+ssh://git@https://github.com/SUNCAT-Center/hori.git',
      install_requires=required,
      package_data={'': ['RELEASE-VERSION']},
      packages=find_packages())
