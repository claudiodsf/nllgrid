# SPDX-License-Identifier: CECILL-2.1
"""
Reading and writing of NonLinLoc grid files.

:copyright:
    2013-2024 Claudio Satriano <satriano@ipgp.fr>,
              Natalia Poiata <poiata@ipgp.fr>,
              Robert Pickle <rpickle@gmail.com>
:license:
    CeCILL Free Software License Agreement v2.1
    (http://www.cecill.info/licences.en.html)
"""
from .NLLGrid import NLLGrid  # noqa
from . import _version
__version__ = _version.get_versions()['version']
