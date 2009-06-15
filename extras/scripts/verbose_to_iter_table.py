#!/usr/bin/env python
import sys,re
pref = 'stepAddFitness'
next_pref = 'calcAvFitness'
len_pref = len(pref)
inp = sys.stdin
outp = sys.stdout

is_step_add = True

outp.write("iter\tlnL\tbestLnL\n")
curr_iter = 1
curr_best = -sys.float_info.max
for line in inp:
	if line.startswith(pref):
		lnL_list = line[len_pref:].strip().split()
		if is_step_add:
			curr_best = -sys.float_info.max
		for i in lnL_list:
			f = float(i)
			if f > curr_best:
				curr_best = f
			outp.write("%d\t%s\t%s\n" % (curr_iter, i, curr_best))
			curr_iter += 1	
			
	elif next_pref and  line.startswith(next_pref):
		pref, next_pref = next_pref, None
		len_pref = len(pref)
		sys.stderr.write("step_add_end = %d\n" % curr_iter)
		is_step_add = False
		curr_best = -sys.float_info.max

outp.write("%d\t%s\t%f\n" % (curr_iter, i, curr_best))
