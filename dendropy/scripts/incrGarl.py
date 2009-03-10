#! /usr/bin/env python

############################################################################
##  sumtrees.py
##
##  Copyright 2008 Jeet Sukumaran.
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

"""
Symmetric difference between collection of trees and a reference tree displayed
    in tabular display as a function of the split cutoff.
"""

__DEBUG__ = True

import os
import sys
import itertools
import copy
from optparse import OptionParser

from dendropy import nexus
from dendropy.splits import encode_splits, lowest_bit_only
from dendropy import treesum
from dendropy import datasets
from dendropy.trees import format_split
from dendropy import treegen
from dendropy import get_logger
from dendropy.datasets import Dataset
_LOG = get_logger("incrGarl.py")

_program_name = 'incrGarl.py'
_program_version = "0.01"
_program_author = 'Mark T. Holder'
_program_contact = 'mtholder@ku.edu'
_program_copyright = "Copyright (C) 2008 Mark T. Holder.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \
                 "There is NO WARRANTY,\nto the extent permitted by law."

TAXON_TO_TRANSLATE = {}

GARLI_GENERAL = (
			"datafname",
			"constraintfile",
			"xstreefname",
			"streefname",
			"runmode",
			"incompletetreefname",
			"attachmentspertaxon",
			"ofprefix",
			"randseed",
			"availablememory",
			"logevery",
			"saveevery",
			"refinestart",
			"outputeachbettertopology",
			"outputcurrentbesttopology",
			"enforcetermconditions",
			"genthreshfortopoterm",
			"scorethreshforterm",
			"significanttopochange",
			"outputphyliptree",
			"outputmostlyuselessfiles",
			"writecheckpoints",
			"restart",
			"outgroup",
			"searchreps",
			"datatype",
			"ratematrix",
			"statefrequencies",
			"ratehetmodel",
			"numratecats",
			"invariantsites",
			)
GARLI_MASTER = (
			"nindivs",
			"holdover",
			"selectionintensity",
			"holdoverpenalty",
			"stopgen",
			"stoptime",
			"startoptprec",
			"minoptprec",
			"numberofprecreductions",
			"treerejectionthreshold",
			"topoweight",
			"modweight",
			"brlenweight",
			"randnniweight",
			"randsprweight",
			"limsprweight",
			"intervallength",
			"intervalstostore",
			"limsprrange",
			"meanbrlenmuts",
			"gammashapebrlen",
			"gammashapemodel",
			"uniqueswapbias",
			"distanceswapbias",
			"bootstrapreps",
			"resampleproportion",
			"inferinternalstateprobs",
			)


DEFAULT_GARLI_CONF = {"datafname" : "rana.nex",
        "constraintfile" : "none",
        "xstreefname" : "rana.tre",
        "streefname" : "incomplete",
        "runmode" : "10",
        "incompletetreefname" : "rana.tre",
        "attachmentspertaxon" : "50",
        "ofprefix" : "rana2.nuc.GTRIG",
        "randseed" : "71836",
        "availablememory" : "512",
        "logevery" : "10",
        "saveevery" : "100",
        "refinestart" : "1",
        "outputeachbettertopology" : "0",
        "outputcurrentbesttopology" : "0",
        "enforcetermconditions" : "1",
        "genthreshfortopoterm" : "20000",
        "scorethreshforterm" : "0.05",
        "significanttopochange" : "0.01",
        "outputphyliptree" : "0",
        "outputmostlyuselessfiles" : "0",
        "writecheckpoints" : "0",
        "restart" : "0",
        "outgroup" : "1",
        "searchreps" : "1",
        "datatype" : "aminoacid",
        "ratematrix" : "jones",
        "statefrequencies" : "empirical",
        "ratehetmodel" : "gamma",
        "numratecats" : "4",
        "invariantsites" : "estimate",
        "nindivs" : "4",
        "holdover" : "1",
        "selectionintensity" : "0.5",
        "holdoverpenalty" : "0",
        "stopgen" : "10",
        "stoptime" : "30",
        "startoptprec" : "0.5",
        "minoptprec" : "0.01",
        "numberofprecreductions" : "10",
        "treerejectionthreshold" : "50.0",
        "topoweight" : "0.0",
        "modweight" : "0.05",
        "brlenweight" : "0.2",
        "randnniweight" : "0.1",
        "randsprweight" : "0.3",
        "limsprweight" : "0.6",
        "intervallength" : "100",
        "intervalstostore" : "5",
        "limsprrange" : "1",
        "meanbrlenmuts" : "5",
        "gammashapebrlen" : "1000",
        "gammashapemodel" : "1000",
        "uniqueswapbias" : "0.1",
        "distanceswapbias" : "1.0",
        "bootstrapreps" : "0",
        "resampleproportion" : "1.0",
        "inferinternalstateprobs" : "0",
        }

def write_garli_conf(out, conf):
    out.write("[general]\n")
    for k in GARLI_GENERAL:
        out.write("%s = %s\n" % (k, conf[k]))
    out.write("[master]\n")
    for k in GARLI_MASTER:
        out.write("%s = %s\n" % (k, conf[k]))

def rev_trans_func(t):
    global TAXON_TO_TRANSLATE
    return TAXON_TO_TRANSLATE[t]

def write_tree_file(outstream, tree, dataset):
    outstream.write("#NEXUS\nBegin Trees;\n  Translate")
    sep = ""
    for n, taxon in enumerate(dataset.taxa_blocks[0]):
        outstream.write(sep)
        sep = ',\n '
        outstream.write(" %d %s " % ((n + 1), nexus.NexusWriter.escape_token(taxon.label)))
    
    outstream.write(";\n Tree a = [&U] %s ;\nEnd;\n" % tree.compose_newick(reverse_translate=rev_trans_func))
    
def addToTree(tree, conf, dataset):
    tmp_tree_filename = ".tmp.tre"
    f = open(tmp_tree_filename, "w")
    write_tree_file(f, tree, dataset)
    f.close()

    conf["incompletetreefname"] = tmp_tree_filename

    tmp_conf_file = ".garli.conf"
    f = open(tmp_conf_file, "w")
    write_garli_conf(f, conf)
    f.close()

    
def read_garli_conf(f):
    default_conf = copy.copy(DEFAULT_GARLI_CONF)
    for line in f:
        s = line.split("=")
        if len(s) > 1:
            k = s[0].strip().lower()
            v = "=".join(s[1:]).strip()
            if not k in default_conf:
                raise RuntimeError("Key %s is not understood" % k)
            default_conf[k] = v
    return default_conf
    
if __name__ == '__main__':
    description =  '%s %s ' % (_program_name, _program_version)    
    usage = "%prog [options] <TREES FILE> [<TREES FILE> [<TREES FILE> [...]]"
    
    parser = OptionParser(usage=usage, add_help_option=True, version = _program_version, description=description)
    parser.add_option('-d','--data',  
                  dest='data_filepath',
                  default=None,
                  help="path to file containing the reference data (must be readable by GARLI and Dendropy)")  
    parser.add_option('-t','--tree',  
                  dest='intree_filepath',
                  default=None,
                  help="path to file containing the starting collection of trees (must be readable by GARLI and Dendropy)")  
    parser.add_option('-c','--conf',  
                  dest='conf',
                  default=None,
                  help="path to a garli conf file")  
    parser.add_option('-v', '--verbose', 
                      action='store_false', 
                      dest='quiet',
                      default=True,
                      help="Verbose mode") 
  
    (opts, args) = parser.parse_args()
    conf_file = opts.conf
    if conf_file is None:
        sys.exit("Expecting a conf file template for GARLI")
    data_file = opts.data_filepath
    intree_file = opts.intree_filepath
    if data_file is None:
        sys.exit("Data file must be specified")
    if intree_file is None:
        sys.exit("Input tree file must be specified")
    for f in [data_file, intree_file, conf_file]:
        if not os.path.exists(f):
            sys.exit("%s does not exist" % f)

    conf = read_garli_conf(open(conf_file, "rU"))
    write_garli_conf(sys.stdout, conf)
    d = Dataset()
    d.read(open(data_file, "rU"), format="NEXUS")
    taxa = d.taxa_blocks[0]
    full_taxa_mask = taxa.all_taxa_bitmask()
    for n, taxon in enumerate(taxa):
        TAXON_TO_TRANSLATE[taxon] = str(n + 1)
    _LOG.debug("%s = full_taxa_mask" % bin(full_taxa_mask))
    assert(len(d.taxa_blocks) == 1)
    characters = d.char_blocks[0]
    assert(len(d.char_blocks) == 1)
    assert(len(characters) == len(taxa))
    inp_trees = d.read_trees(open(intree_file, "rU"), format="NEXUS")
    assert(inp_trees)
    current_taxon_mask = None
    for tree in inp_trees:
        assert tree.taxa_block is taxa
        encode_splits(tree)
        if current_taxon_mask is None:
            current_taxon_mask = tree.seed_node.edge.clade_mask
            _LOG.debug("%s = current_taxon_mask" % bin(current_taxon_mask))
            assert( (current_taxon_mask | full_taxa_mask) == full_taxa_mask)
            toadd_taxon_mask = current_taxon_mask ^ full_taxa_mask
        else:
            assert(current_taxon_mask == tree.seed_node.edge.clade_mask)
    next_toadd = lowest_bit_only(current_taxon_mask^full_taxa_mask)
    if (next_toadd - 1) != current_taxon_mask:
        _LOG.debug("%s = next_toadd" % format_split(next_toadd, taxa=taxa))
        _LOG.debug("%s = current_taxon_mask\n(next_toadd - 1) != current_taxon_mask" % format_split(current_taxon_mask, taxa=taxa))
        sys.exit("In this version, taxa must be added to the tree in the order that they appear in the matrix")


    conf["datafname"] = data_file
    for tree in inp_trees:
        trees = addToTree(tree, conf, d)
    sys.exit(0)
    
    
    
    
    ###################################################
    # Support file idiot checking
        
    sampled_filepaths = []        
    missing = False 
    for fpath in args:
        fpath = os.path.expanduser(os.path.expandvars(fpath))        
        if not os.path.exists(fpath):
            sys.exit('Sampled trees file not found: "%s"' % fpath)
        sampled_filepaths.append(fpath)
    if not sampled_filepaths:
        sys.exit("Expecting arguments indicating files that contain sampled trees")
        
    sampled_file_objs = [open(f, "rU") for f in sampled_filepaths]

    ###################################################
    # Lots of other idiot-checking ...
    
    # target tree
    if opts.reference_tree_filepath is None:
        sys.exit("A reference tree must be specified (use -h to see all options)")
    reference_tree_filepath = os.path.expanduser(os.path.expandvars(opts.reference_tree_filepath))
    if not os.path.exists(reference_tree_filepath):
        sys.exit('Reference tree file not found: "%s"\n' % reference_tree_filepath)

    d = Dataset()
    ref_trees  = d.read_trees(open(reference_tree_filepath, 'ru'), format="NEXUS")
    
    if len(ref_trees) != 1:
        sys.exit("Expecting one reference tree")
    ref_tree = ref_trees[0]
    splits.encode_splits(ref_tree)
    assert(len(d.taxa_blocks) == 1)
    taxa = d.taxa_blocks[0]
    
    
    ###################################################
    # Main work begins here: Count the splits
    
    start_time = datetime.datetime.now()
    
    comments = []
    tsum = treesum.TreeSummarizer()
    tsum.burnin = 0
    if opts.quiet:
        tsum.verbose = False
        tsum.write_message = None
    else:
        tsum.verbose = True
        tsum.write_message = sys.stderr.write




    _LOG.debug("### COUNTING SPLITS ###\n")                
    split_distribution = tsum.count_splits_on_trees(tree_files=sampled_filepaths, tree_iterator=nexus.iterate_over_trees, taxa_block=taxa) 
        
    report = []
    report.append("%d trees read from %d files." % (tsum.total_trees_read, len(sampled_filepaths)))
    report.append("%d trees ignored in total." % (tsum.total_trees_ignored))    
    report.append("%d trees considered in total for split support assessment." % (tsum.total_trees_counted))
    report.append("%d unique taxa across all trees." % len(split_distribution.taxa_block))
    num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits = split_distribution.splits_considered()
    report.append("%d unique splits out of %d total splits counted." % (num_unique_splits, num_splits))
    report.append("%d unique non-trivial splits out of %d total non-trivial splits counted." % (num_nt_unique_splits, num_nt_splits))
        
    _LOG.debug("\n".join(report))


    con_tree = treegen.star_tree(taxa)
    taxa_mask = taxa.all_taxa_bitmask()
    splits.encode_splits(con_tree)
    leaves = con_tree.leaf_nodes()
    
    to_leaf_dict = {}
    for leaf in leaves:
        to_leaf_dict[leaf.edge.clade_mask] = leaf
    unrooted = True
    n_read = float(tsum.total_trees_read)
    sp_list = []
    for split, count in split_distribution.split_counts.iteritems():
        freq = count/n_read
        if not splits.is_trivial_split(split, taxa_mask):
            m = split & taxa_mask
            if (m != taxa_mask) and ((m-1) & m): # if not root (i.e., all "1's") and not singleton (i.e., one "1")
                if unrooted:
                    c = (~m) & taxa_mask
                    if (c-1) & c: # not singleton (i.e., one "0")
                        if 1 & m:
                            k = c
                        else:
                            k = m
                        sp_list.append((freq, k, m))
                else:
                    sp_list.append((freq, m, m))
    sp_list.sort(reverse=True)
    
    root = con_tree.seed_node
    root_edge = root.edge

    curr_freq = 1.1
    curr_all_splits_list = []
    curr_compat_splits_list = []
    all_splits_by_freq = []
    compat_splits_by_freq = []

    # Now when we add splits in order, we will do a greedy, extended majority-rule consensus tree
    for freq, split_to_add, split_in_dict in sp_list:
        if abs(curr_freq-freq) > 0.000001:
            # dropping down to the next lowest freq
            curr_l = [freq, []]
            curr_all_splits_list = curr_l[1]
            all_splits_by_freq.append(curr_l)
            curr_l = [freq, []]
            curr_compat_splits_list = curr_l[1]
            compat_splits_by_freq.append(curr_l)
            curr_freq = freq
        
        curr_all_splits_list.append(split_to_add)

        if (split_to_add & root_edge.clade_mask) != split_to_add:
            continue
        lb = splits.lowest_bit_only(split_to_add)
        one_leaf = to_leaf_dict[lb]
        parent_node = one_leaf
        while (split_to_add & parent_node.edge.clade_mask) != split_to_add:
            parent_node = parent_node.parent_node
        if parent_node is None or parent_node.edge.clade_mask == split_to_add:
            continue # split is not in tree, or already in tree.
        
        new_node = trees.Node()
        new_node_children = []
        new_edge = new_node.edge
        new_edge.clade_mask = 0
        for child in parent_node.child_nodes():
            # might need to modify the following if rooted splits
            # are used
            cecm = child.edge.clade_mask
            if (cecm & split_to_add ):
                assert cecm != split_to_add
                new_edge.clade_mask |= cecm
                new_node_children.append(child)
        # Check to see if we have accumulated all of the bits that we
        #   needed, but none that we don't need.
        if new_edge.clade_mask == split_to_add:
            for child in new_node_children:
                parent_node.remove_child(child)
                new_node.add_child(child)
            parent_node.add_child(new_node)
            con_tree.split_edges[split_to_add] = new_edge
            curr_compat_splits_list.append(split_to_add)
    ref_set = set()
    for s in ref_tree.split_edges.iterkeys():
        m = s & taxa_mask
        if 1 & m:
            k = (~m) & taxa_mask
        else:
            k = m
        if not splits.is_trivial_split(k, taxa_mask):
            ref_set.add(k)
        
    all_set = set()
    compat_set = set()
    
    _LOG.debug("%d edges is the reference tree" % (len(ref_set)))

    print "freq\tcompatFP\tcompatFN\tcompatSD\tallFP\tallFN\tallSD"
    for all_el, compat_el in itertools.izip(all_splits_by_freq, compat_splits_by_freq):
        freq = all_el[0]
        all_sp = all_el[1]
        all_set.update(all_sp)
        all_fn = len(ref_set - all_set)
        all_fp = len(all_set - ref_set)
        compat_sp = compat_el[1]
        compat_set.update(compat_sp)
        compat_fn = len(ref_set - compat_set)
        compat_fp = len(compat_set - ref_set)
        
        print "%f\t%d\t%d\t%d\t%d\t%d\t%d" % (freq, compat_fp, compat_fn, compat_fp + compat_fn, all_fp, all_fn, all_fp + all_fn )
        
