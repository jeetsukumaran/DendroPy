#!/usr/bin/env python
import sys, re
pat= re.compile(r'\[(\d+) rep (\d+)\s*\] .*\[!GarliScore ([-.0-9]+)\]')
inp = open(sys.argv[1], 'rU')
max_t = int(sys.argv[2])
stop_gen_arg= int(sys.argv[3])
outp = sys.stdout
d = {}

for line in inp:
	m = pat.match(line)
	if not m:
		sys.exit("No match %s" % line)
	g = m.groups()
	tree_n = int(g[0])
	rep_n = int(g[1])
	score = float(g[2])
	sc_list = d.get(tree_n, [])
	sc_list.append(score)
	d[tree_n] = sc_list


sorted_keys = d.keys()
sorted_keys.sort()

def calc_iter(n_taxa, stop_gen=stop_gen_arg, n_indiv=4):
	swap_iter = stop_gen*n_indiv
	per_t = 5
	sadd_iter = 0
	for i in xrange(4,n_taxa):
		sadd_iter += per_t
		per_t += 2
	return sadd_iter + swap_iter
	
def write_row(outp, ntaxa, sc_list):
	outp.write("%d" % calc_iter(ntaxa, stop_gen_arg*(ntaxa-3)))
	m = max(sc_list)
	outp.write("\t%f" % m)
	n = 0
	for n, i in enumerate(sc_list):
		outp.write("\t%f" % i)
	for i in range(n, 9):
		outp.write("\tNA")
	outp.write("\n")
	
outp.write("iter\tbest\tzer\tone\ttwo\tthr\tfou\tfiv\tsix\tsev\teig\tnin\n")
for i in sorted_keys:
	if i > max_t:
		break
	write_row(outp, i, d[i])
	
