#! /usr/bin/env python

############################################################################
##  simulate-tree.py
##
##  command-line interface for tree simulation
##
##  Copyright 2009 Mark Holder
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


import time, random, sys, os
from optparse import OptionParser
from dendropy import TaxonSet
from dendropy.treesim import uniform_pure_birth
if __name__ == '__main__':

    parser = OptionParser(usage='%prog [options]',
        add_help_option=True,
        version='0.0.1',
        description="""Provides a command-line interface to some of dendropy's tree simulation functions.
Default is Yule tree shape with random U[0,1] branch lengths.
""")

    parser.add_option('-l', '--leaves',
        action='store',
        dest='num_leaves',
        type=int,
        default=4,
        help='The number of leaves to generate.')

    parser.add_option('-s', '--seed',
        action='store',
        dest='seed',
        type=long,
        default=0,
        help='The seed of the random number generator')

    (opts, args) = parser.parse_args()
    num_leaves = opts.num_leaves
    if num_leaves < 3:
        sys.exit("The number of leaves must be greater than 2")
    seed = opts.seed
    if seed < 0:
        sys.exit("Seed must be positive")
    if seed < 1:
        seed = time.time()
    rng = random.Random()
    sys.stderr.write("seed = %ld\n" % seed)
    rng.seed(seed)
    outstr = sys.stdout
    
    
    taxa = TaxonSet(["t%d" % i for i in xrange(1, 1 + num_leaves)])
    tree = uniform_pure_birth(taxa, rng=rng)
    for e in tree.preorder_edge_iter():
        e.length = rng.random()
    outstr.write("%s;\n" % tree.as_newick_string(preserve_spaces=True))
        

