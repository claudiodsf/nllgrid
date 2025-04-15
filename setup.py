# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
"""
A minimal setup script for SourceSpec.

This script ensures compatibility with versioneer
and dynamically generates a README for correctly resolving paths on GitHub.

All the remaining configuration is in pyproject.toml.
"""
from setuptools import setup
import versioneer

# Dynamically generate the README text for PyPI, replacing local paths
# with the corresponding URLs on GitHub.
with open('README.md', 'rb') as f:
    long_description = f.read().decode('utf-8').replace(
        ': CHANGELOG.md',
        ': https://github.com/claudiodsf/nllgrid/blob/main/CHANGELOG.md'
    )

setup(
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass()
)
