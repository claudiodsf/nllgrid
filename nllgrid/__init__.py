# SPDX-License-Identifier: GPL-3.0-or-later
"""
Reading and writing of NonLinLoc grid files.

:copyright:
    2013-2025 Claudio Satriano <satriano@ipgp.fr>,
              Natalia Poiata <poiata@ipgp.fr>,
              Robert Pickle <rpickle@gmail.com>
:license:
    GNU General Public License v3.0 or later
    (https://www.gnu.org/licenses/gpl-3.0-standalone.html)
"""
from .NLLGrid import NLLGrid  # noqa
from . import _version
__version__ = _version.get_versions()['version']
