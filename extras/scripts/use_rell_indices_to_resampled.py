#!/usr/bin/env python

import sys
import re

fn = sys.argv[1]
f = open(fn, 'rU')


garli_dna_model_pattern = r'\[!GarliModel  r ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) e ([.0-9]*) ([.0-9]*) ([.0-9]*) [.0-9]* a ([.0-9]*) p ([.0-9]*) \]'
#garli_dna_model_pattern = r'\[!GarliModel\s*(r [.0-9]* [.0-9]* [.0-9]* [.0-9]* [.0-9]* e [.0-9]* [.0-9]* [.0-9]* [.0-9]* a [.0-9]* p [.0-9]*)\s*\]'
garli_after_name_tm_pattern = r'\s*=\s*\[&U\]\[!GarliScore ([-.0-9]*)\]\s*%s\s*(\(.*\))\s*;' % garli_dna_model_pattern
garli_after_name_tree_pattern = r'\s*=\s*\[.*\]\s*(\(.*\))\s*;'
tree_prefix = r'\[iGarli (\d+) \] tree best'
tree_with_model_pat_str = tree_prefix + garli_after_name_tm_pattern
#sys.stderr.write(tree_with_model_pat_str)
garli_tree_model_pat = re.compile(tree_with_model_pat_str)
garli_tree_pat = re.compile(tree_prefix + garli_after_name_tree_pattern)
#sys.stderr.write("garli_tree_pat = %s\n" % tree_prefix + garli_after_name_tree_pattern)

first = True
sys.stdout.write("""#NEXUS
begin trees;
""")
tree_list = []
for line in f:
	m = garli_tree_model_pat.search(line)
	if m:
		tree_list.append(line)		
	elif garli_tree_pat.search(line):
		sys.exit("Tree without model found")
	elif False:
		sys.stderr.write("No match: %s" % line)


rell_ind_filename = sys.argv[2]
rell_ind_file = open(rell_ind_filename, 'rU')
for line in rell_ind_file:
	ls = line.strip()
	if ls:
		i = int(ls)
		sys.stdout.write(tree_list[i])


sys.stdout.write("end;\n")
