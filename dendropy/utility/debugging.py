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
Debugging support.
"""

import linecache

# http://www.dalkescientific.com/writings/diary/archive/2005/04/20/tracing_python_code.html
def trace_run(frame, event, arg):
    """
    Prints module name and line number of each line as it is executed.
    To activate::
        # /usr/bin/env python
        import sys
        from dendropy.utility.debugging import trace_run
        sys.settrace(trace_run)
        .
        (code)
        .
    """
    if event == "line":
        lineno = frame.f_lineno
        filename = frame.f_globals["__file__"]
        if filename == "<stdin>":
            filename = "debugging.py"
        if (filename.endswith(".pyc") or
            filename.endswith(".pyo")):
            filename = filename[:-1]
        name = frame.f_globals["__name__"]
        line = linecache.getline(filename, lineno)
        print "%s:%s: %s" % (name, lineno, line.rstrip())
    return trace_run
