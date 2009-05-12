#! /usr/bin/env python

############################################################################
##  setup.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
############################################################################


"""
Package setup and installation.
"""

import os
import sys
import random

__all__ = ["base.py",
        "characters.py",
        "chargen.py",
        "charmodels.py",
        "coalescent.py",
        "dataio.py",
        "datasets.py",
        "distributions.py",
        "fasta.py",
        "nexml.py",
        "nexus.py",
        "phylip.py",
        "splits.py",
        "taxa.py",
        "treedists.py",
        "treecalc.py",
        "treegen.py",
        "trees.py",
        "treesum.py",
        "utils.py",
        "xmlparser.py",
        "high_level.py",
]
__version__ = "2.4.0"
PACKAGE_NAME = "DendroPy"
PACKAGE_VERSION = __version__
PACKAGE_AUTHOR = "Jeet Sukumaran and Mark T. Holder"
PACKAGE_COPYRIGHT = "Copyright 2008 Jeet Sukumaran and Mark T. Holder."
PACKAGE_LICENSE = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
"""

# Global random number generator
GLOBAL_RNG = random.Random()

def python_version():
    "Returns Python version as float."
    major_ver = sys.version_info[0]
    minor_ver = sys.version_info[1]
    return major_ver + (float(minor_ver)/10)

def is_python_at_least(version):
    """
    Returns True if Python version is at least as high as the argument
    (a numeric value).
    """
    return bool(python_version() >= version)

###############################################################################
## LOGGING

import logging

_LOGGING_LEVEL_ENVAR="DENDROPY_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR="DENDROPY_LOGGING_FORMAT"

# global debugging flag
if "DENDROPY_DEBUG" in os.environ:
    if os.environ["DENDROPY_DEBUG"] \
        and os.environ["DENDROPY_DEBUG"].lower()[0] in ["1", "t", "y", "d"]:
        GLOBAL_DEBUG = True
    else:
        GLOBAL_DEBUG = False
else:
    GLOBAL_DEBUG = False

def get_logging_level():
    if _LOGGING_LEVEL_ENVAR in os.environ:
        if os.environ[_LOGGING_LEVEL_ENVAR].upper() == "NOTSET":
            level = logging.NOTSET
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "DEBUG":
            level = logging.DEBUG 
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "INFO":
            level = logging.INFO
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "WARNING":
            level = logging.WARNING
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "ERROR":
            level = logging.ERROR
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
    else:
        level = logging.NOTSET  
    return level        

def get_logger(name="dendropy"):
    """
    Returns a logger with name set as given, and configured 
    to the level given by the environment variable _LOGGING_LEVEL_ENVAR.
    """
    logger_set = False    
#     package_dir = os.path.dirname(module_path)
#     config_filepath = os.path.join(package_dir, _LOGGING_CONFIG_FILE)    
#     if os.path.exists(config_filepath):
#         try:
#             logging.config.fileConfig(config_filepath)
#             logger_set = True
#         except:
#             logger_set = False
    logger = logging.getLogger(name)            
    if not logger_set:    
        level = get_logging_level()
        rich_formatter = logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)d): %(levelname) 8s: %(message)s")
        simple_formatter = logging.Formatter("%(levelname) 8s: %(message)s")
        raw_formatter = logging.Formatter("%(message)s")
        default_formatter = None
        logging_formatter = default_formatter
        if _LOGGING_FORMAT_ENVAR in os.environ:
            if os.environ[_LOGGING_FORMAT_ENVAR].upper() == "RICH":
                logging_formatter = rich_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "SIMPLE":
                logging_formatter = simple_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "NONE":
                logging_formatter = None
            else:
                logging_formatter = default_formatter
        else:
            logging_formatter = default_formatter
        if logging_formatter is not None:            
            logging_formatter.datefmt='%H:%M:%S'
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
    return logger

_user_ini_checked = False
if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.dendropy/startup.py")
    if os.path.exists(p):
        execfile(p)
    del p

if __name__ == "__main__":
    sys.stdout.write("%s %s\n" % (PACKAGE_NAME, PACKAGE_VERSION))

def deprecation(message, logger_obj=None, stacklevel=3):
    try:
        import warnings
        warnings.warn(message, DeprecationWarning, stacklevel=stacklevel)
    except:
        if logger_obj:
            logger_obj.warning(message)
    
