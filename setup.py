# -*- coding: utf-8 -*-
"""setup.py: setuptools control."""
from setuptools import setup

# import inspect
# import os
# import sys

# Import the version string.
# path = os.path.join(os.path.abspath(os.path.dirname(inspect.getfile(
#     inspect.currentframe()))), 'nllgrid')
# sys.path.insert(0, path)
# from ssp_version import get_git_version


with open('README.md', 'rb') as f:
    long_descr = f.read().decode('utf-8')


setup(
    name='nllgrid',
    packages=['nllgrid', ],
    include_package_data=True,
    version=0.1,
    description='Python class for reading and writing NonLinLoc grid files',
    long_description=long_descr,
    author='Claudio Satriano',
    author_email='satriano@ipgp.fr',
    url='http://www.ipgp.fr/~satriano',
    license='CeCILL Free Software License Agreement, Version 2.1',
    platforms='OS Independent',
    classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: CeCILL Free Software License '
                'Agreement, Version 2.1',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics'],
    install_requires=['numpy>=1.0', 'scipy>=0.16']
    )
