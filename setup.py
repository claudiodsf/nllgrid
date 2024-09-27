# -*- coding: utf-8 -*-
# SPDX-License-Identifier: CECILL-2.1
"""setup.py: setuptools control."""
from setuptools import setup
import versioneer

with open('README.md', 'rb') as f:
    long_descr = f.read().decode('utf-8').replace(
        '(CHANGELOG.md)',
        '(https://cdn.jsdelivr.net/gh/claudiodsf/nllgrid/CHANGELOG.md)'
    )

setup(
    name='nllgrid',
    packages=['nllgrid', ],
    include_package_data=True,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Python class for reading and writing NonLinLoc grid files',
    long_description=long_descr,
    long_description_content_type='text/markdown',
    author='Claudio Satriano',
    author_email='satriano@ipgp.fr',
    url='https://github.com/claudiodsf/nllgrid',
    license='CeCILL Free Software License Agreement, Version 2.1',
    platforms='OS Independent',
    classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre '
                'License, version 2.1 (CeCILL-2.1)',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3.12',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics'],
    python_requires='>=3.7',
    install_requires=['numpy>=1.0', 'scipy>=0.16', 'pyproj']
    )
