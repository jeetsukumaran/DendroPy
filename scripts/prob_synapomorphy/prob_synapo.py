#! /usr/bin/env python
############################################################################
##  prob-synapo.py
##
##  Uses GARLI to calculate the probabilities of different state changes across
##      each edge in a tree
##
##  Copyright 2010 Mark Holder
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

"""Uses a special version of GARLI to calculate the probabilities of different
state changes across each edge in a tree
"""
import sys
import os
import itertools
import math
import shutil
from cStringIO import StringIO
from optparse import OptionGroup, OptionParser

from dendropy.utility.cli import confirm_overwrite
from dendropy import DataSet
from dendropy.treecalc import symmetric_difference
import subprocess
import dendropy

BOGUS_TAXON_LABELS = ['bogus', 'bogus2']

def write_split_file(edge, taxon_set):
    o = open('split.tre', 'w')
    contained = set()
    o.write('((')
    for n, nd in enumerate(dendropy.Node.leaf_iter(edge.head_node)):
        if n > 0:
            o.write(',')
        o.write(nd.taxon.label)
        contained.add(nd.taxon)
    o.write(')')
    for leaf in taxon_set:
        if (leaf not in contained) and (leaf.label not in BOGUS_TAXON_LABELS):
            o.write(',%s' % leaf.label)
    o.write(');\n')

def convert_data_vec(v, s):
    nc = len(v)
    for i in xrange(nc):
        v[i].value = s

def garli_safe_add_child(par, child):
    "Adds an extra node so that the tree stays binary"
    attacher_nd = dendropy.Node()
    par_children = par.child_nodes()
    len(par_children) == 2
    r_anc = par_children[-1]
    r_anc_e_len = r_anc.edge.length
    par.add_child(attacher_nd, edge_length=0.0)
    attacher_nd.add_child(child, edge_length=0.0)
    par.remove_child(r_anc)
    attacher_nd.add_child(r_anc, edge_length=r_anc_e_len)
    return attacher_nd, r_anc, r_anc_e_len


def write_garli_tree_file(tree, e, anc_taxon, des_taxon, filename, model_string):
    anc_attach = e.tail_node
    assert anc_attach
    anc_nd = dendropy.Node(taxon=anc_taxon)
    anc_extra, old_a_child, old_a_child_e_len = garli_safe_add_child(anc_attach, anc_nd)


    if e.is_terminal():
        des_attach = anc_attach
    else:
        des_attach = e.head_node
        assert des_attach
    des_nd = dendropy.Node(taxon=des_taxon)
    des_extra, old_d_child, old_d_child_e_len = garli_safe_add_child(des_attach, des_nd)

    try:
        o = open(filename, 'w')
        o.write('''#NEXUS
Begin trees;
tree inp = [&U] ''')
        tree.write(stream=o, schema="newick")
        o.write(''';
END;
BEGIN GARLI;
        %s ;
END;
''' % model_string)
        o.close()
    finally:
        des_attach.add_child(old_d_child, edge_length=old_d_child_e_len)
        des_attach.remove_child(des_extra)
        anc_attach.add_child(old_a_child, edge_length=old_a_child_e_len)
        anc_attach.remove_child(anc_extra)


GARLI_CONF_TEMPLATE = '''
[general]
datafname = %(datafname)s
constraintfile = none
streefname = %(streefname)s
attachmentspertaxon = 50
ofprefix = %(ofprefix)s
randseed = -1
availablememory = 2000
logevery = 1
saveevery = 100
refinestart = 0
refineend = 0
outputeachbettertopology = 0
outputsitelikelihoods = 1
outputcurrentbesttree = 0
enforcetermconditions = 1
genthreshfortopoterm = 10000
scorethreshforterm = 0.05
significanttopochange = 0.01
outputphyliptree = 0
outputmostlyuselessfiles = 0
writecheckpoints = 0
restart = 0
searchreps = 1
runmode = 6

datatype = %(datatype)s
ratematrix = %(ratematrix)s
statefrequencies = %(statefrequencies)s
invariantsites = %(invariantsites)s
ratehetmodel = %(ratehetmodel)s
numratecats = %(numratecats)d

[master]
nindivs = 4
holdover = 1
selectionintensity = .5
holdoverpenalty = 0
stopgen = 1
stoptime = 5000000

startoptprec = .00001
minoptprec = .00001
numberofprecreductions = 5
treerejectionthreshold = 5
topoweight = 0
modweight = .0
brlenweight = 0.
randnniweight = 0.0
randsprweight = 0.0
limsprweight =  0.0
intervallength = 100
intervalstostore = 5
uniqueswapbias = 0.1
distanceswapbias = 1.0

limsprrange = 6
meanbrlenmuts = 5
gammashapebrlen = 1000
gammashapemodel = 1000

bootstrapreps = 0
inferinternalstateprobs = 0
'''
MY_PREFIX = 'with_bogus_'

def write_garli_confs(data_files):
    tree_files = []
    o_files = []
    c_files = []
    for fn in data_files:
        assert fn.startswith(MY_PREFIX)
        #tree_files.append(MY_PREFIX + fn[len(MY_PREFIX):-4] + '.tre')
        # tree file does not need to change!
        tree_files.append(MY_PREFIX + '.tre')
        o_files.append(MY_PREFIX + fn[len(MY_PREFIX):-4] + '_out')
        c_files.append(MY_PREFIX + fn[len(MY_PREFIX):-4] + '.conf')
    # TEMPORARY -- these should be inferred from the input!
    fmt_d = { 'datatype' : 'dna',
              'ratematrix' : 'fixed',
              'statefrequencies' : 'fixed',
              'invariantsites' : 'none', #'fixed',
              'ratehetmodel' : 'gammafixed',
              'numratecats' : 4,
             }
    for d, t, o, cfn in itertools.izip(data_files, tree_files, o_files, c_files):
        fmt_d['datafname'] = d
        fmt_d['streefname'] = t
        fmt_d['ofprefix'] = o
        cfo = open(cfn, 'w')
        cfo.write(GARLI_CONF_TEMPLATE % fmt_d)
        cfo.close()
    return c_files, [i + '.sitelikes.log' for i in o_files]

def move_files(fn_list, d):
    for fn in fn_list:
        fp = os.path.join(d, fn)
        print fp, fn
        try:
            os.rename(fp, fn)
        except:
            print "failed on ", fp, fn

def invoke_garli(cf):
    invok = ['scoreOnlyGarli', cf]
    g = subprocess.Popen(invok)
    g.wait()
    if g.returncode != 0:
        sys.exit("""Invocation:\n%s\nfrom %s failed!

Note that "scoreOnlyGarli" is a special form of GARLI(Zwickl) available from Mark Holder on request.
This program calculates likelihoods without changing any parameters of the model or branch lengths)
""" % (" ".join(invok), os.path.abspath(os.curdir)))

def parse_site_likes(file_obj):
    first_line = file_obj.next().strip()
    assert first_line == 'Tree	-lnL	Site	-lnL'
    lnl_list = []
    for n, line in enumerate(file_obj):
        numb = 1 + n
        s = line.strip().split()
        if len(s) == 2:
            assert int(s[0]) == numb
            lnl_list.append(-float(s[1]))
    return lnl_list

def write_ti_probs(outf, per_ti, full_site_lnL):
    for el in per_ti:
        assert len(el) == len(full_site_lnL)
    for n, partial_el_list in enumerate(itertools.izip(*per_ti)):
        full_sll = full_site_lnL[n]
        outf.write('%d\t' % (n + 1))
        print "partial_el_list =", partial_el_list
        mx_ll = max(partial_el_list)
        print "mx_ll =", mx_ll
        scaled_l = [math.exp(i - mx_ll) for i in partial_el_list]
        f_l = sum(scaled_l)
        post_prob = [i/f_l for i in scaled_l]
        recalc_lnl = math.log(f_l) + mx_ll
        print "recalc_lnl = ", recalc_lnl
        print "full_sll = ", full_sll
        assert abs(recalc_lnl - full_sll) < 0.00001
        outf.write('%s\n' % '\t'.join(["%.6f" % i for i in post_prob]))

if __name__ == "__main__":
    """
    Main CLI handler.
    """

    parser = OptionParser(add_help_option=True,
        version="0.1",
        description="Uses a special version of GARLI to calculate the probabilities of different state changes across each edge in a tree")

    parser.add_option('-s','--separator',
        dest='separator',
        default='\t',
        help="character to use to separate/delimit columns (default=<TAB>)")

    parser.add_option('-o','--output',
        dest='output_filepath',
        default=None,
        help="path to output file (if not given, will print to standard output)")

    parser.add_option('-m', '--model',
        dest='model',
        default=None,
        help="GARLI model string")

    parser.add_option('-r', '--replace',
        action='store_true',
        dest='replace',
        default=False,
        help="replace/overwrite output file without asking if it already exists ")

    parser.add_option('-t', '--tree-file',
        dest='tree',
        default="",
        help="Tree file")

    parser.add_option('-d', '--data-file',
        dest='data',
        default="",
        help="Data file (NEXUS)")

    (opts, args) = parser.parse_args()

    scratch_dir = 'tmp_dir_prob_syn'
    if not os.path.exists(scratch_dir):
        os.mkdir(scratch_dir)
    scratch_dir = os.path.abspath(scratch_dir)
    tree_filepaths = []
    if not opts.model:
        sys.exit('Expecting a GARLI model string')
    if not opts.tree:
        sys.exit('Expecting a tree file as an argument (use -h to see all options)')
    if not opts.data:
        sys.exit('Expecting a tree file as an argument (use -h to see all options)')
    if not os.path.exists(opts.data):
        sys.exit('Data file not found: "%s"' % opts.data)
    if not os.path.exists(opts.tree):
        sys.exit('Tree file not found: "%s"' % opts.tree)

    tree_file_objs = [open(f, "rU") for f in tree_filepaths]

    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if not confirm_overwrite(output_fpath, opts.replace):
            sys.exit(1)
        output_dest = open(output_fpath, "w")


    dataset = DataSet()
    ts = dendropy.TaxonSet()
    dataset.read(stream=open(opts.data, 'rU'), schema='NEXUS', taxon_set=ts)
    if len(dataset.char_matrices) != 1:
        sys.exit("Currently the script only supports data files with a single character matrix")
    if len(dataset.tree_lists) > 0:
        sys.exit("Currently the script does not support trees in the data file")
    dataset.read(stream=open(opts.tree, 'rU'), schema='NEXUS', taxon_set=ts)
    if (len(dataset.tree_lists) != 1) or len(dataset.tree_lists[0]) != 1:
        sys.exit("Currently the script only not tree files with a single tree")

    if ts.get_taxon(label=BOGUS_TAXON_LABELS[0]) or ts.get_taxon(label=BOGUS_TAXON_LABELS[1]):
        sys.exit('Give me a break. You really have a taxon named "%s" or "%s" in your data!?\nI refuse to deal with this file.\n' % (BOGUS_TAXON_LABELS[0], BOGUS_TAXON_LABELS[1]))
    tree = dataset.tree_lists[0][0]

    matrix = dataset.char_matrices[0]
    if len(matrix.state_alphabets) != 1:
        sys.exit('Expecting the character matrix to have exactly one "state alphabet". Found %d' % len(matrix.state_alphabets))
    if os.path.exists('results'):
        sys.exit('The filepath "results" already exists. Move this file or directory out of the way before running this script.')
    state_alphabet = matrix.state_alphabets[0]

    a_taxon = ts[0]
    data_for_a_taxon = matrix.get(a_taxon)
    num_char = len(data_for_a_taxon)

    all_missing = '?' * num_char # TEMP! should get this from the data type
    data_type = 'dna' # TEMP! should get this from the datafile


    fake_data = StringIO('2 %d\n%-9s %s\n%-9s %s\n' % (num_char, BOGUS_TAXON_LABELS[0], all_missing, BOGUS_TAXON_LABELS[1], all_missing))
    fake_dataset = DataSet()
    fake_dataset.read(stream=fake_data, schema='phylip', data_type=data_type)
    fd = fake_dataset.char_matrices[0]

    matrix.extend(fd)

    fn = MY_PREFIX + 'missing.nex'
    const_data_files = [fn]
    data_files = []
    o = open(os.path.join(scratch_dir, fn), 'w')
    matrix.write_to_stream(o, schema='nexus')
    o.close()

    fundamental_states = [i for i in state_alphabet.fundamental_states() if i is not state_alphabet.gap]
    anc_taxon = ts.get_taxon(label=BOGUS_TAXON_LABELS[0])
    assert(anc_taxon)
    anc_taxon_data = matrix.get(anc_taxon)
    assert(anc_taxon_data)

    des_taxon = ts.get_taxon(label=BOGUS_TAXON_LABELS[1])
    assert(des_taxon)
    des_taxon_data = matrix.get(des_taxon)
    assert(des_taxon_data)

    for from_state in fundamental_states:
        from_state_label = from_state.symbol
        convert_data_vec(anc_taxon_data, from_state)
        for to_state in fundamental_states:
            to_state_label = to_state.symbol
            convert_data_vec(des_taxon_data, to_state)
            fn = "%s%s_to_%s.nex" % (MY_PREFIX, from_state_label, to_state_label)
            if from_state == to_state:
                const_data_files.append(fn)
            else:
                data_files.append(fn)
            o = open(os.path.join(scratch_dir, fn), 'w')
            matrix.write_to_stream(o, schema='nexus')
            o.close()

    to_move = list(data_files)
    to_move.extend(const_data_files)
    orig_dir = os.path.abspath(os.curdir)
    os.chdir(scratch_dir)
    try:
        conf_files, out_files = write_garli_confs(data_files)
        to_move.extend(conf_files)
        const_conf_files, const_out_files = write_garli_confs(const_data_files)
        to_move.extend(const_conf_files)
    finally:
        os.chdir(orig_dir)


    SL_SUFFIX = '_out.sitelikes.log'
    ti_order = []
    for from_state in fundamental_states:
        from_state_label = from_state.symbol
        for to_state in fundamental_states:
            to_state_label = to_state.symbol
            ti_order.append((from_state_label, to_state_label))


    edge_list = [i for i in tree.preorder_edge_iter()]
    for nd in tree.preorder_node_iter():
        if (nd.parent_node and len(nd.child_nodes()) > 2) or ((nd.parent_node is None) and len(nd.child_nodes()) > 3):
            sys.exit('oops! polytomy in tree. GARLI will not behave correctly!')

    for n, e in enumerate(edge_list):
        if not e.tail_node:
            continue
        sub_dir = os.path.join(scratch_dir, 'split%d' %n)
        os.mkdir(sub_dir)

        os.chdir(sub_dir)
        move_files(to_move, scratch_dir)
        try:
            write_split_file(e, ts)
            write_garli_tree_file(tree, e, anc_taxon=anc_taxon, des_taxon=des_taxon, filename=MY_PREFIX + '.tre', model_string=opts.model)
            if e.is_terminal():
                cf_list = const_conf_files
            else:
                cf_list = conf_files + const_conf_files
            for cf in cf_list:
                invoke_garli(cf)
            if not e.is_terminal():
                full_site_lnL = parse_site_likes(open(MY_PREFIX + 'missing' + SL_SUFFIX, 'rU'))
                s_ti_lnL_list = []
                for from_state_label, to_state_label in ti_order:
                    fn = "%s%s_to_%s%s" % (MY_PREFIX, from_state_label, to_state_label, SL_SUFFIX)
                    s_ti_lnL_list.append(parse_site_likes(open(fn, 'rU')))
                outf = open('change_probs.txt', 'w')
                outf.write('pos\t%s\n' % '\t'.join(['%s_to_%s' % i for i in ti_order]))
                write_ti_probs(outf, s_ti_lnL_list, full_site_lnL)
                outf.close()
        finally:
            os.chdir(orig_dir)
        os.chdir(scratch_dir)
        move_files(to_move, sub_dir)
        os.chdir(orig_dir)

    os.mkdir('results')
    for n, e in enumerate(edge_list):
        if not e.tail_node:
            continue
        src_dir = os.path.join(scratch_dir, 'split%d' %n)
        dest_dir = os.path.join('results', 'split%d' %n)
        posteriors = os.path.join(src_dir, 'change_probs.txt')
        if os.path.exists(posteriors):
            os.mkdir(dest_dir)
            shutil.copyfile(posteriors, os.path.join(dest_dir, 'change_probs.txt'))
            split_tree = os.path.join(src_dir, 'split.tre')
            shutil.copyfile(split_tree, os.path.join(dest_dir, 'split.tre'))
