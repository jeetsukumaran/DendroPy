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
Low-level data structures, classes, constants, and services that are used
by the DendroPy library but generally not directly by client code.
"""

import os
import sys
import random
import subprocess
from cStringIO import StringIO

###############################################################################
## VERSION INFORMATION UTILITIES

def python_version():
    "Returns Python version as float."
    major_ver = sys.version_info[0]
    minor_ver = sys.version_info[1]
    return major_ver + (float(minor_ver)/10)

def get_current_git_head(dirpath=os.path.curdir):
    # git rev-parse HEAD
    # git show --quiet --pretty=format:%H
    # git log -1 --pretty=format:%H
    # git rev-list -1 HEAD
    p = subprocess.Popen(['git show --quiet --pretty=format:%H'],
        shell=True,
        cwd=os.path.abspath(os.path.expandvars(dirpath)),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout = p.stdout.read()
    if len(stdout) == 0:
        return "[NOT A GIT REPOSITORY]"
    else:
        return stdout.replace("\n","")

def get_current_git_branch(dirpath=os.path.curdir):
    p = subprocess.Popen(['git branch'],
        shell=True,
        cwd=os.path.abspath(os.path.expandvars(dirpath)),
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    branches = [t for t in p.stdout.read().split("\n") if t]
    if len(branches) == 0:
        return "[NOT A GIT REPOSITORY]"
    else:
        for b in branches:
            if b.startswith('* '):
                return b[2:]
        return "[UNIDENTIFIABLE]"

def get_current_git_tag(dirpath=os.path.curdir):
    p = subprocess.Popen(['git describe --tag'],
            shell=True,
        cwd=os.path.abspath(os.path.expandvars(dirpath)),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    if p.stderr.read() == 'fatal: Not a git repository\n':
        return '[NOT A GIT REPOSITORY]'
    tags = [t for t in p.stdout.read().split("\n") if t]
    if len(tags) == 0:
        return 'UNSPECIFIED'
    else:
        return tags[-1]

###############################################################################
## USER-SPECIFIC

# global debugging flag
if "DENDROPY_DEBUG" in os.environ:
    if os.environ["DENDROPY_DEBUG"] \
        and os.environ["DENDROPY_DEBUG"].lower()[0] in ["1", "t", "y", "d"]:
        GLOBAL_DEBUG = True
    else:
        GLOBAL_DEBUG = False
else:
    GLOBAL_DEBUG = False

_user_ini_checked = False
if not _user_ini_checked:
    import os
    _user_ini_checked = True
    p = os.path.expanduser("~/.dendropy/startup.py")
    if os.path.exists(p):
        execfile(p)
    del p

###############################################################################
## GLOBAL RANDOM NUMBER GENERATOR

GLOBAL_RNG = random.Random()
