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
Extracts coalescent frames (# alleles : waiting time for coalescence) from a
set of trees.
"""

import sys
import os
from optparse import OptionGroup
from optparse import OptionParser
from dendropy import datasets
from dendropy import coalescent
from dendropy import distributions

_prog_usage = '%prog [options] <tree-files>'
_prog_version = 'WTD Version 1.0'
_prog_description = 'returns distribution of waiting times of coalescent frames of a set of trees'
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
        
    parser.add_option('-s', '--summarize-means',
        action='store',
        dest='summarize_means',
        default=None,
        metavar='FILENAME',
        help='summarize means to this file (default="%default")')
        
    parser.add_option('-n', '--pop-size', '-N',
        action='store',
        dest='pop_size',
        type='int',
        default=1,
        metavar='Ne',
        help='effective HAPLOID population size (for calculation of expected distribution means; default=%default [assumes edge lengths are in units of Ne])')
        
#     parser.add_option('-o', '--output-prefix',
#         action='store',
#         dest='output_prefix',
#         default="wt",
#         metavar='OUTPUT-PREFIX',
#         help='prefix for output file names (default="%default")')

    (opts, args) = parser.parse_args()
    
    if len(args) == 0:
        sys.stderr.write("%s" % parser.get_usage())
        sys.exit(1)
       
    output = sys.stdout
    output.write("k\twaiting_time\n")
    
    coal_frames = {}

    for a in args:
        fpath = os.path.expandvars(os.path.expanduser(a))
        if not os.path.exists(fpath):
            sys.stderr.write('File not found: "%s"\n' % fpath)
        else:
            sys.stderr.write('Reading: "%s"\n' % fpath)
            d = datasets.Dataset()
            ctrees = d.read_trees(open(fpath, "rU"), "NEXUS")
            for t in ctrees:
                cf = coalescent.extract_extract_coalescent_frames(t)
                for k, wt in cf:
                    output.write("%d\t%s\n" % (k, wt))
                    if k not in coal_frames:
                        coal_frames[k] = []
                    coal_frames[k].append(wt)
                    
    if opts.summarize_means is not None:
        smfile = open(os.path.expandvars(os.path.expanduser(opts.summarize_means)), "w")
        smfile.write("k\tmean_wt\texpected_wt\n")
        for k, wt in coal_frames.items():
            actual_mean = float(sum(wt))/len(wt)        
            expected_mean = float(opts.pop_size) / distributions.binomial_coefficient(k, 2)
            smfile.write("%d\t%s\t%s\n" % (k, actual_mean, expected_mean))
                
if __name__ == "__main__":
    main()