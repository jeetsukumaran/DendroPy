#! /usr/bin/env python

############################################################################
##  prob-coal-tree.py
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
Calculates the probability of trees under a neutral unstructured fixed effective
population size coalescent.
"""

import sys
import os
from optparse import OptionGroup
from optparse import OptionParser
from dendropy import datasets
from dendropy import coalescent

_prog_usage = '%prog [options] <tree-files>'
_prog_version = 'PROB-COAL-TREE Version 1.0'
_prog_description = 'This program returns the log probability of a tree under a coalescent model.'
_prog_author = 'Jeet Sukumaran and Mark T. Holder'
_prog_copyright = 'Copyright (C) 2009 Jeet Sukumaran and Mark T. Holder.'

def main():
    """
    Main CLI handler.
    """
    
    parser = OptionParser(usage=_prog_usage, 
        add_help_option=True, 
        version=_prog_version, 
        description=_prog_description)    
       
    parser.add_option('-n', '--pop-size', '-N',
        action='store',
        dest='pop_size',
        type='int',
        default=1,
        metavar='Ne',
        help='effective HAPLOID population size (default=%default [assumes edge lengths are in units of Ne])')

    (opts, args) = parser.parse_args()
    
    if len(args) == 0:
        sys.stderr.write("%s" % parser.get_usage())
        sys.exit(1)
        
    for a in args:
        fpath = os.path.expandvars(os.path.expanduser(a))
        if not os.path.exists(fpath):
            sys.stderr.write('File not found: "%s"\n' % fpath)
        else:
            sys.stderr.write('Reading: "%s"\n' % fpath)
            d = datasets.Dataset()
            ctrees = d.read_trees(open(fpath, "rU"), "NEXUS")
            for t in ctrees:
                p = coalescent.log_probability_of_coalescent_tree(t, opts.pop_size)
                sys.stdout.write("%s\n" % p)

if __name__ == '__main__':
    main()

    
