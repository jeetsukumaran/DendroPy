#!/usr/bin/env python
import sys, re
pat= re.compile(r'(\d+)\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)')
inp = open(sys.argv[1], 'rU')
n_tax = int(sys.argv[2])

def calc_iter(n_taxa, stop_gen=200, n_indiv=4):
	swap_iter = stop_gen*n_indiv
	per_t = 5
	sadd_iter = 0
	for i in xrange(3,n_taxa):
		sadd_iter += per_t
		per_t += 2
	return sadd_iter + swap_iter

outp = sys.stdout
master_sc_list = []
parsing = False
sc_list = []
for line in inp:
    if parsing:
        m = pat.match(line)
        if m:
            g = m.groups()
            gen_n = int(g[0])
            score = float(g[1])
            n_score = calc_iter(n_tax, stop_gen=gen_n)
            sc_list.append((n_score, score))
        else:
            sys.stderr.write("break parse at %s" % line)
            parsing = False
            master_sc_list.append(sc_list)
            sc_list = []
    elif line.startswith("gen"):
        parsing = True
	
def write_row(outp, ntaxa, sc_list):
	outp.write("%d" % calc_iter(ntaxa))
	m = max(sc_list)
	outp.write("\t%f" % m)
	n = 0
	for n, i in enumerate(sc_list):
		outp.write("\t%f" % i)
	for i in range(n, 9):
		outp.write("\tNA")
	outp.write("\n")
	
outp.write("iter\tbest\tzer\tone\ttwo\tthr\tfou\tfiv\tsix\tsev\teig\tnin\n")
n_points = len(master_sc_list[0])

it_list = [iter(i) for i in master_sc_list]
for ind in xrange(n_points):
    try:
        v = [i.next() for i in it_list]
    except:
        break
    s = [i[1] for i in v]
    niter = v[0][0]
    for i in v:
        assert i[0] == niter
    b = max(s)
    outp.write("%d\t%f" % (niter, b))
    for i in s:
        outp.write("\t%f" % (i))
    outp.write("\n")

