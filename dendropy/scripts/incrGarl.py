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
import re
from optparse import OptionParser
from subprocess import Popen, PIPE

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


def run_garli(conf, commands):
    tmp_conf_file = ".garli.conf"
    f = open(tmp_conf_file, "w")
    write_garli_conf(f, conf)
    f.close()
    garli_instance = Popen(["iGarli", tmp_conf_file], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    gstdout, gstderr = garli_instance.communicate("\n".join(commands))
    rc = garli_instance.wait()
    if rc != 0:
        sys.exit(gstderr)


GARLI_SCORE_PATTERN = re.compile(r"\[!GarliScore ([-0-9.]+)\]")
def read_garli_scores(inp):
    sc = []
    for line in inp:
        m = GARLI_SCORE_PATTERN.search(line)
        if m:
            g = m.group(1)
            sc.append(float(g))
    return sc

def add_to_tree(tree, conf, dataset, tree_ind):
    ofprefix = "from%d" % tree_ind
    conf["ofprefix"] = ofprefix

    tmp_tree_filename = ".tmp.tre"
    f = open(tmp_tree_filename, "w")
    write_tree_file(f, tree, dataset)
    f.close()

    conf["incompletetreefname"] = tmp_tree_filename

    
    run_garli(conf, ["run", "quit"])
    
    output_tree = ofprefix + ".best.tre"
    t = dataset.read_trees(open(output_tree, "rU"), format="NEXUS")
    del dataset.trees_blocks[-1]
    sc = read_garli_scores(open(output_tree, "rU"))
    if len(t) != len(sc):
        sys.exit("Did not read the same number of trees (%d) as scores (%d) from %s" % (len(t), len(sc), output_tree))
    for otree, osc in itertools.izip(t, sc):
        otree.score = osc
    return t
    
    
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
    dataset = Dataset()
    dataset.read(open(data_file, "rU"), format="NEXUS")
    taxa = dataset.taxa_blocks[0]
    full_taxa_mask = taxa.all_taxa_bitmask()
    for n, taxon in enumerate(taxa):
        TAXON_TO_TRANSLATE[taxon] = str(n + 1)
    _LOG.debug("%s = full_taxa_mask" % bin(full_taxa_mask))
    assert(len(dataset.taxa_blocks) == 1)
    characters = dataset.char_blocks[0]
    assert(len(dataset.char_blocks) == 1)
    assert(len(characters) == len(taxa))
    inp_trees = dataset.read_trees(open(intree_file, "rU"), format="NEXUS")
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


    conf["datafname"] = os.path.join("..", data_file)


    n = len(dataset.taxa_blocks[0])
    dirn = "t%d" % n
    if not os.path.exists(dirn):
        os.makedirs(dirn)
        if not os.path.exists(dirn):
            sys.exit("Could not make %s" % dirn)
    if not os.path.isdir(dirn):
        sys.exit("%s is not a directory" % dirn)
    orig_dir = os.getcwd()
    os.chdir(dirn)

    try:
        for tree_ind, tree in enumerate(inp_trees):
            trees = add_to_tree(tree, conf, dataset, tree_ind)
            for t in trees:
                print t.score
    finally:
        os.chdir(orig_dir)
        
    sys.exit(0)
    
