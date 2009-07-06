#! /usr/bin/env python

############################################################################
##  incr_garli.py
##
##  Copyright 2009 Mark T. Holder
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


import sys
import re

fn = sys.argv[1]
f = open(fn, 'rU')
wfn = sys.argv[2]
wf = open(wfn, 'rU')
wt_lines = wf.readlines()

#garli_dna_model_pattern = r'\[!GarliModel  r ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) e ([.0-9]*) ([.0-9]*) ([.0-9]*) [.0-9]* a ([.0-9]*) p ([.0-9]*) \]'
garli_dna_model_pattern = r'\[!GarliModel\s*(r [.0-9]* [.0-9]* [.0-9]* [.0-9]* [.0-9]* e [.0-9]* [.0-9]* [.0-9]* [.0-9]* a [.0-9]* p [.0-9]*)\s*\]'
garli_after_name_tm_pattern = r'\s*=\s*\[&U\]\[!GarliScore ([-.0-9]*)\]\s*%s\s*(\(.*\))\s*;' % garli_dna_model_pattern
garli_after_name_tree_pattern = r'\s*=\s*\[.*\]\s*(\(.*\))\s*;'
tree_prefix = r''
garli_tree_model_pat = re.compile(tree_prefix + garli_after_name_tm_pattern)
garli_tree_pat = re.compile(tree_prefix + garli_after_name_tree_pattern)
#sys.stderr.write("garli_tree_pat = %s\n" % tree_prefix + garli_after_name_tree_pattern)

first = True
bootrep = -1
for line in f:
    m = garli_tree_model_pat.search(line)
    if m:
        has_model = True
    else:
        m = garli_tree_pat.search(line)
        has_model = False
    if m:
        g = m.groups()
        if has_model:
            score = float(g[0])
            sys.stdout.write('model = %s\n' % g[1])
        sys.stdout.write('tree = %s\n' % g[-1])
        if first: # this is a way to make 
            sys.stdout.write('treeNum = 1\n')
            first = False
        if bootrep >= 0:
            sys.stdout.write('setcharintwts = %s' % wt_lines[bootrep])
        sys.stdout.write('run\n')
        bootrep += 1
    elif True:
        sys.stderr.write("nomatch: %s\n" %line)
sys.exit(0)

