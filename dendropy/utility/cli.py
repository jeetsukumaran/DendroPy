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
Support for CLI operations.
"""
import os
import sys

import dendropy

def confirm_overwrite(filepath,
                      replace_without_asking=False,
                      file_desc="Output",
                      out=sys.stdout):
    if os.path.exists(filepath):
        if replace_without_asking:
            overwrite = 'y'
        else:
            out.write('%s file already exists: "%s"\n' % (file_desc, filepath))
            overwrite = raw_input("Overwrite (y/N)? ")
        if not overwrite.lower().startswith("y"):
            return False
        else:
            return True
    else:
        return True

def show_splash(prog_name,
        prog_subtitle,
        prog_version,
        prog_author,
        prog_copyright,
        dest=sys.stderr,
        extended=False,
        width=70):

    lines = []
    lines.append("%s - %s" % (prog_name, prog_subtitle))
    lines.append("%s" % prog_version)
    lines.append("By %s" % prog_author)
    lines.append("(using the DendroPy Phylogenetic Computing Library Version %s)" % (dendropy.__version__))
    if extended:
        lines.append('')
        lines.extend(prog_copyright.split('\n'))
    if width is None or width <= 0:
        width = max([len(i) for i in lines]) + 1
    sbars = '=' * width
    dest.write("%s\n" % sbars)
    dest.write("%s\n" % ('\n'.join(lines)))
    dest.write("%s\n\n" % sbars)


