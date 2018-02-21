#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Handling deprecation warnings and messages correctly.
"""

import os
import warnings
from dendropy.utility import metavar

DEPRECATION_WARNING_FILTER = None
_DEPRECATION_WARNINGS_CONFIGURED = False

class CriticalDeprecationWarning(UserWarning):
    pass

def configure_deprecation_warning_behavior(warning_filter=None):
    global DEPRECATION_WARNING_FILTER
    global _DEPRECATION_WARNINGS_CONFIGURED
    if warning_filter is None:
        warning_filter = os.environ.get(metavar.DEPRECATION_WARNING_FILTER, "default")
    DEPRECATION_WARNING_FILTER = warning_filter
    warnings.simplefilter(DEPRECATION_WARNING_FILTER,
            CriticalDeprecationWarning)
    _DEPRECATION_WARNINGS_CONFIGURED = True

def _initialize_deprecation_warnings():
    global _DEPRECATION_WARNINGS_CONFIGURED
    if not _DEPRECATION_WARNINGS_CONFIGURED:
        configure_deprecation_warning_behavior()

def dendropy_deprecation_warning(**kwargs):
    _initialize_deprecation_warnings()
    leader = "  # "
    stacklevel = kwargs.pop("stacklevel", 3)
    if "message" in kwargs:
        message = kwargs["message"]
    elif "old_construct" in kwargs or "new_construct" in kwargs:
        message = []
        message.append("")
        if "preamble" in kwargs:
            message.append(leader + kwargs["preamble"])
        message.append(leader + "Instead of:")
        for construct in kwargs["old_construct"].split("\n"):
            message.append(leader + "    {}".format(construct))
        message.append(leader + "Use:")
        for construct in kwargs["new_construct"].split("\n"):
            message.append(leader + "    {}".format(construct))
        if "epilog" in kwargs:
            message.append(leader + kwargs["epilog"])
        message = "\n".join(message)
    _initialize_deprecation_warnings()
    old_formatwarning = warnings.formatwarning
    warnings.warn(
            message=message,
            category=CriticalDeprecationWarning,
            stacklevel=stacklevel,
            )
