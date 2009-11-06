#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Controls which tests are run.
"""

import os

FAST, NORMAL, SLOW, EXHAUSTIVE = 0, 10, 20, 30

###############################################################################
## TESTING LEVELS

class TestLevel:

    def level_name(i):
        if i <= FAST:
            return "FAST"
        if i <= NORMAL:
            return "NORMAL"
        if i <= SLOW:
            return "SLOW"
        return "EXHAUSTIVE"

    level_name = staticmethod(level_name)

    def name_to_int(l):
        try:
            return int(l)
        except:
            pass
        l = l.upper()
        if l == "FAST":
            return FAST
        if l == "NORMAL":
            return NORMAL
        if l == "SLOW":
            return SLOW
        if l == "EXHAUSTIVE":
            return TEXHAUSTIVE
        raise ValueError("Test run level %s unrecognized" % l)

    name_to_int = staticmethod(name_to_int)

def fast_testing_notification(logger, test_object_name, message=None, level=FAST):
    if message is None:
        message = ""
    else:
        message = ": %s" % message
    logger.warning('%s Testing Level: skipping %s tests in %s%s' \
        % (TestLevel.level_name(get_current_testing_level()),
           TestLevel.level_name(level),
           test_object_name,
           message))

def get_current_testing_level():
    l = os.environ.get("DENDROPY_TESTING_LEVEL")
    if l is None:
        if "DENDROPY_FAST_TESTS" in os.environ:
            return FAST
        return NORMAL
    try:
        return name_to_int(l)
    except:
        _LOG.warn("the value %s for DENDROPY_TESTING_LEVEL is not recognized.  Using NORMAL level" % l)
    return NORMAL

def is_test_enabled(level, logger=None, module_name="", message=None):
    tl = get_current_testing_level()
    if level > tl:
        if logger:
            fast_testing_notification(logger, module_name, message, level)
        return False
    return True
