#!/usr/bin/env python
import sys
import shutil
import os
import re

from dendropy import get_logger
_LOG = get_logger("igarli_master.py")



script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
if len(sys.argv) < 3:
    sys.exit("Expecting a data file (NEXUS) and a tree file as arguments")
dataf = open(sys.argv[1], 'rU')

dim_pat = re.compile(r'\s*dimensions\s*ntax\s*=\s*(\d+)\s*nchar\s*=\s*(\d+)', re.I)
rev_dim_pat = re.compile(r'\s*dimensions\s*nchar\s*=\s*(\d+)\s*ntax\s*=\s*(\d+)', re.I)
n_taxa = None
n_char = None
for line in dataf:
    m = dim_pat.match(line)
    if m:
        n_taxa, n_char = [int(i) for i in m.groups()]
        break
    else:
        m = rev_dim_pat.match(line)
        if m:
            n_char, n_taxa = [int(i) for i in m.groups()]
            break
if n_char is None:
    sys.exit("Dimensions command not found")

fmt_cmd_pat = re.compile('\s*format', re.I)
mat_cmd_pat = re.compile('\s*matrix', re.I)
fmt_cmd = None
for line in dataf:
    print line
    if fmt_cmd_pat.match(line):
        fmt_cmd = line[:-1]
    elif mat_cmd_pat.match(line):
        break
assert fmt_cmd

data_lines = []
taxa = []
for line in dataf:
    if len(line) > n_char:
        data_lines.append(line)
        taxon = line.split()[0]
        taxa.append(taxon)
    if len(data_lines) == n_taxa:
        break

tree_file = os.path.abspath(sys.argv[2])

try:
    starting_n_tax = int(sys.argv[3])
except IndexError:
    starting_n_tax = 4

orig_dir = os.path.abspath(os.curdir)
for nt in range(starting_n_tax, n_taxa + 1):
    subdir = 'incrSubdir%d' % nt
    prev_dir = 'incrSubdir%d' % (nt - 1)
    if os.path.exists(subdir):
        sys.exit("the directory %s already exists. Exiting..." % subdir)
    if nt >= 5 and not os.path.exists(prev_dir):
        sys.exit("the directory %s does not exists. Exiting..." % prev_dir)
    os.mkdir(subdir)
    os.chdir(subdir)
    _LOG.debug("Changing directory into %s" % subdir)
    try:
        if nt == 4:
            shutil.copy(tree_file, 'incrgarli.tre')
        else:
            df = open('data.nex', 'w')
            df.write("""#NEXUS
Begin Data;
    Dimensions ntax = %d nchar = %d ;
    %s
Matrix
%s
;
end;
""" % (nt, n_char, fmt_cmd, "".join(data_lines[:nt])))

            df.close()
            translate_content = ",\n  ".join(["%d %s" % (n + 1, s) for n, s in enumerate(taxa[:nt])])
            star_tree_content = """#NEXUS
begin trees;
    translate %s;
tree o = [&U] (%s);
end;
""" % (translate_content, ",".join([str(1+i) for i in range(nt)]))
            tf = open('.tmp.tre', 'w')
            tf.write(star_tree_content)
            tf.close()
            script = os.path.join(script_dir, 'igarli_one_round.sh')
            cmd = '%s ../incrSubdir%d/incrgarli.tre %d' % (script, nt - 1, nt)
            _LOG.debug("Invoking:\n  %s\nfrom %s" % (cmd, os.path.abspath(os.curdir)))
            rc = os.system(cmd)
            if not rc == 0:
                sys.exit("Exiting because %s failed with an exit code of %d" % (cmd, rc))
    finally:
        os.chdir(orig_dir)
