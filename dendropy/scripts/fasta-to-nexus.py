#! /usr/bin/env python

import os
import sys

if ("--help" in sys.argv) or ("-?" in sys.argv):
    sys.stderr.write("usage: fasta-to-nexus.py [<fasta-file-path>] [<nexus-file-path>]\n")
    sys.exit(1)
   
if len(sys.argv) < 2:
    src = sys.stdin
else:
    src_fpath = os.path.expanduser(os.path.expandvars(sys.argv[1]))
    if not os.path.exists(src_fpath):
        sys.stderr.write('Not found: "%s"' % src_fpath)
    src = open(src_fpath)        
       
if len(sys.argv) < 3:
    dest = sys.stdout # os.path.splitext(src_fpath)[0] + ".nex"
else:
    dest_fpath = os.path.expanduser(os.path.expandvars(sys.argv[2]))
    dest = open(dest_fpath, "w")

seqs = {}
cur_seq = None
lines = src.readlines()
for i in lines:
    i = i.replace("\n", "").replace("\r", "")
    if i:
        if i.startswith(">"):
            label = i[1:]
            seqs[label] = []
            cur_seq = seqs[label]
        else:
            if cur_seq is None:
                raise Exception("Sequence data found before label")
            cur_seq.extend(i.replace(" ", ""))
       
taxlabels = seqs.keys()

dest.write("#NEXUS\n\n")

dest.write("Begin Taxa;\n")
dest.write("  dimensions ntax=%d;\n" % len(seqs))
dest.write("  taxlabels\n")
for taxlabel in taxlabels:
    dest.write("    %s\n" % taxlabel)
dest.write("  ;\n")
dest.write("End;\n\n")

nchar = max([len(s) for s in seqs.values()])
dest.write("Begin Characters;\n")
dest.write("   dimensions nchar=%d;\n" % nchar)
dest.write("   format datatype=dna missing=? gap=-;\n")
dest.write("   matrix\n")
for taxlabel in taxlabels:
    dest.write("    %s      %s\n" % (taxlabel, "".join(seqs[taxlabel])))
dest.write("  ;\n")
dest.write("end;\n")
