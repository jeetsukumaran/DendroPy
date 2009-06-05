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


#garli_dna_model_pattern = r'\[!GarliModel  r ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) e ([.0-9]*) ([.0-9]*) ([.0-9]*) [.0-9]* a ([.0-9]*) p ([.0-9]*) \]'
garli_dna_model_pattern = r'\[!GarliModel\s*(r [.0-9]* [.0-9]* [.0-9]* [.0-9]* [.0-9]* e [.0-9]* [.0-9]* [.0-9]* [.0-9]* a [.0-9]* p [.0-9]*)\s*\]'
garli_after_name_tm_pattern = r'\s*=\s*\[&U\]\[!GarliScore ([-.0-9]*)\]\s*%s\s*(\(.*\))\s*;' % garli_dna_model_pattern
garli_after_name_tree_pattern = r'\s*=\s*\[&U\]\[!GarliScore ([-.0-9]*)\]\s*(\(.*\))\s*;'
tree_prefix = r'\s*Tree\s*tree(\d+)'
garli_tree_model_pat = re.compile(tree_prefix + garli_after_name_tm_pattern)
garli_tree_pat = re.compile(tree_prefix + garli_after_name_tree_pattern)


first = True
for line in f:
	m = garli_tree_model_pat.match(line)
	if m:
		has_model = True
	else:
		m = garli_tree_pat.match(line)
		has_model = False
	if m:
		g = m.groups()
		result_number = g[0]
		score = float(g[1])
		if has_model:
			sys.stdout.write('model = %s\n' % g[2])
		sys.stdout.write('tree = %s\n' % g[-1])
		if first: # this is a way to make 
			sys.stdout.write('treeNum = 1\n')
			first = False
		sys.stdout.write('run\n')
	elif False:
		sys.stderr.write("nomatch: %s\n" %line)
sys.exit(0)

