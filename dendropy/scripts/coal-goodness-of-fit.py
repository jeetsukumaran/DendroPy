#! /usr/bin/env python

############################################################################
##  coal-goodness-of-fit.py
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
Applies a goodness-of-fit test to evaluate the fit of a set of trees to 
Kingman's neutral coalescent model
"""

import sys
import os
from optparse import OptionGroup
from optparse import OptionParser
from dendropy import datasets
from dendropy import coalescent
from dendropy import distributions

_prog_usage = '%prog [options] <tree-files>'
_prog_version = 'COAL-GOODNESS-OF-FIT Version 1.0'
_prog_description = "applies a goodness-of-fit test to evaluate the fit of a set of trees to Kingman's neutral coalescent model"
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
       
    coal_frames = {}
    num_trees = 0
    for a in args:
        fpath = os.path.expandvars(os.path.expanduser(a))
        if not os.path.exists(fpath):
            sys.stderr.write('File not found: "%s"\n' % fpath)
        else:
            sys.stderr.write('Reading: "%s"\n' % fpath)
            d = datasets.Dataset()
            ctrees = d.read_trees(open(fpath, "rU"), "NEXUS")
            for t in ctrees:
                num_trees += 1
                cf = coalescent.coalescent_frames(t)
                for k, wt in cf:
                    if k not in coal_frames:
                        coal_frames[k] = []
                    coal_frames[k].append(wt)

#     sum_x = 0
#     n = 0
#     for k, wts in coal_frames.items():
#         exp_mean = float(opts.pop_size) / distributions.binomial_coefficient(k, 2)
#         for wt in wts:
#             x = pow((wt - exp_mean), 2) / exp_mean
#             sum_x += x
#             n += 1
#     sys.stdout.write("X=%s, n=%s, df=%s\n" % (sum_x, n, n-1))
    sys.stderr.write("k\tmean_wt\texpected_wt\tx\n")
    sum_x = 0
    
    for k, wt in coal_frames.items():
        obs_mean = float(sum(wt))/len(wt)        
        exp_mean = float(opts.pop_size) / distributions.binomial_coefficient(k, 2)
        x = pow((obs_mean - exp_mean), 2) / exp_mean
        sum_x += x       
        sys.stderr.write("%d\t%s\t%s\t%s\n" % (k, obs_mean, exp_mean, x))
    sys.stdout.write("X^ = %s (num. trees = %s)\n" % (sum_x, num_trees))
if __name__ == "__main__":
    main()