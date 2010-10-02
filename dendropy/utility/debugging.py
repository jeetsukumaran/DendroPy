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
