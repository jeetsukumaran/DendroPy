#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Controls which tests are run.
"""

import os
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

FAST, NORMAL, SLOW, EXHAUSTIVE = 0, 10, 20, 30

###############################################################################
## TESTING LEVELS

class TestLevel(object):

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
            return EXHAUSTIVE
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
        return TestLevel.name_to_int(l)
    except ValueError:
        _LOG.warn("the value %s for DENDROPY_TESTING_LEVEL is not recognized.  Using NORMAL level" % l)
    return NORMAL

def is_test_enabled(level, logger=None, module_name="", message=None):
    tl = get_current_testing_level()
    if level > tl:
        if logger:
            fast_testing_notification(logger, module_name, message, level)
        return False
    return True
