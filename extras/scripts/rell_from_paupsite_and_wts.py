#!/usr/bin/env python
"""Takes a paup sitelikes score file and one or more weight files produced by the bootigarli scripts
    returns the index (starting at 0 of the tree with the best likelihood in each rell bootstrap rep."""
import sys, itertools
paupsitelikefn = sys.argv[1]
inp = open(paupsitelikefn, 'rU')
n_char = None
curr_sl = []
expected_site_ind = 1
expected_tree_ind = 1
all_sl = []

for line_number, line in enumerate(inp):
    if line.startswith('Tree\t-lnL\tSite\t-lnL'):
        if n_char is None:
            if len(curr_sl) == 0:
                continue # this will happen on the first line of the file
            n_char = len(curr_sl)
    elif line.strip():
        s = line.split()
        assert(len(s) == 2)
        n = int(s[0])
        f = float(s[1])
        if n == expected_site_ind:
            curr_sl.append(f)
            expected_site_ind += 1
        else:
            if n != expected_tree_ind:
                sys.exit("Could not parse line %d (%s) . Expecting site %d or tree %d" % (line_number, line[:-1], expected_site_ind, expected_tree_ind))
            if n_char is None:
                n_char = len(curr_sl)
            elif n_char != len(curr_sl):
                sys.exit("Expecting %d sites but found %d at line %d" % (n_char, len(curr_sl), line_number))
            all_sl.append(curr_sl)
            if sum(curr_sl) - f > 0.0001:
                sys.exit("Expecting a total likelihood of %f (based on sum of site likes) but got %f at line %d" % (sum(curr_sl), f, line_number))
            curr_sl = []
            expected_tree_ind += 1
            expected_site_ind = 1

if n_char is None:
    sys.exit("Parse failed to find any sites")
if curr_sl:
    sys.exit("Partial sitelikes found for tree %d" % expected_tree_ind)

if len(sys.argv) == 2:
    sys.stdout.write("all_site_likes = ")
    sys.stdout.write(repr(all_sl))
    sys.stdout.write("\n")
    sys.exit()

for fn in sys.argv[2:]:
    sys.stderr.write('fn = %s\n' % fn)
    f = open(fn, 'rU')
    for line_number, line in enumerate(f):
        wts = [float(i) for i in line.strip().split()]
        if len(wts) != n_char:
            sys.exit("Found only %d pattern weights (instead of %d) at line %d of %s" % (len(wts), n_char, line_number, fn))
        best = -sum([lnL*wt for lnL, wt in itertools.izip(all_sl[0], wts)])
        bestInd = 0
        #sys.stderr.write("best initial = %f\n" % best)
        for offset, site_like_vec in enumerate(all_sl[1:]):
            wtLnL = -sum([lnL*wt for lnL, wt in itertools.izip(site_like_vec, wts)])
            if wtLnL > best:
                best = wtLnL
                bestInd = offset + 1
        #sys.stderr.write("best final = %f\n" % best)
        sys.stdout.write("%d\n" % bestInd)
