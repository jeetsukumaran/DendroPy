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
VERBOSE = True
import os
import sys
import itertools
import copy
import re
import logging
from optparse import OptionParser
from subprocess import Popen, PIPE

from dendropy import nexus
from dendropy.splits import encode_splits, lowest_bit_only, iter_split_indices, find_edge_from_split
from dendropy.characters import CharactersBlock
from dendropy.taxa import TaxaBlock
from dendropy.treedists import symmetric_difference
from dendropy import treesum
from dendropy import datasets
from dendropy.trees import format_split, TreesBlock
from dendropy import treegen
from dendropy import get_logger
from dendropy.datasets import Dataset
from dendropy.utils import LineReadingThread
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

class GARLI_ENUM:
    NORMAL_RUNMODE, INCR_RUNMODE = "0", "10"

class GarliConf(object):
    def __init__(self):
        self.datafname = "rana.nex",
        self.constraintfile =  "none"
        self.xstreefname =  "rana.tre"
        self.streefname =  "incomplete"
        self.runmode =  "10"
        self.incompletetreefname =  "rana.tre"
        self.attachmentspertaxon =  "50"
        self.ofprefix =  "rana2.nuc.GTRIG"
        self.randseed =  "71836"
        self.availablememory =  "512"
        self.logevery =  "10"
        self.saveevery =  "100"
        self.refinestart =  "1"
        self.outputeachbettertopology =  "0"
        self.outputcurrentbesttopology =  "0"
        self.enforcetermconditions =  "1"
        self.genthreshfortopoterm =  "20000"
        self.scorethreshforterm =  "0.05"
        self.significanttopochange =  "0.01"
        self.outputphyliptree =  "0"
        self.outputmostlyuselessfiles =  "0"
        self.writecheckpoints =  "0"
        self.restart =  "0"
        self.outgroup =  "1"
        self.searchreps =  "1"
        self.datatype =  "aminoacid"
        self.ratematrix =  "jones"
        self.statefrequencies =  "empirical"
        self.ratehetmodel =  "gamma"
        self.numratecats =  "4"
        self.invariantsites =  "estimate"
        self.nindivs =  "4"
        self.holdover =  "1"
        self.selectionintensity =  "0.5"
        self.holdoverpenalty =  "0"
        self.stopgen =  "10"
        self.stoptime =  "30"
        self.startoptprec =  "0.5"
        self.minoptprec =  "0.01"
        self.numberofprecreductions =  "10"
        self.treerejectionthreshold =  "50.0"
        self.topoweight =  "0.0"
        self.modweight =  "0.05"
        self.brlenweight =  "0.2"
        self.randnniweight =  "0.1"
        self.randsprweight =  "0.3"
        self.limsprweight =  "0.6"
        self.intervallength =  "100"
        self.intervalstostore =  "5"
        self.limsprrange =  "1"
        self.meanbrlenmuts =  "5"
        self.gammashapebrlen =  "1000"
        self.gammashapemodel =  "1000"
        self.uniqueswapbias =  "0.1"
        self.distanceswapbias =  "1.0"
        self.bootstrapreps =  "0"
        self.resampleproportion =  "1.0"
        self.inferinternalstateprobs =  "0"
        self.garli_instance = None

    def write_garli_conf(self, out):
        conf = self.__dict__
        out.write("[general]\n")
        for k in GARLI_GENERAL:
            out.write("%s = %s\n" % (k, conf[k]))
        out.write("[master]\n")
        for k in GARLI_MASTER:
            out.write("%s = %s\n" % (k, conf[k]))

    def run(self, commands, terminate_run=True):
        if self.garli_instance is None:
            tmp_conf_file = ".garli.conf"
            f = open(tmp_conf_file, "w")
            self.write_garli_conf(f)
            f.close()
            invoc = ["iGarli", tmp_conf_file]
            _LOG.debug("Running:\n  %s\nfrom\n  %s" % (" ".join(invoc), os.path.abspath(os.curdir)))
            _LOG.debug("### BEGIN COMMANDS ###:\n%s\n###END COMMANDS ###" % ("\n".join(commands)))
            if terminate_run:
                commands.append("quit")
            s = Popen(invoc, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            self.garli_instance = s
            self.stderrThread = LineReadingThread(stream=s.stderr, store_lines=True, subproc=s)
            self.stderrThread.start()
            self.stdoutThread = LineReadingThread(stream=s.stdout, store_lines=True, subproc=s)
            self.stdoutThread.start()
        if terminate_run:
            gstdout, gstderr = garli_instance.communicate("\n".join(commands))
            if VERBOSE:
                print gstdout
            rc = garli_instance.wait()
            if rc != 0:
                sys.exit(gstderr)
            self.garli_instance = None
        else:
            garli_prompt = "iGarli>"
            for command in commands:
                self.garli_instance.stdin.write(command)
                stderr_line = ""
                while not stderr_line.startswith(garli_prompt):
                    stderr_line = self.garli_instance.stderr.readline()
                    print "garli stderr =", stderr_line
                    stdout_lines = self.garli_instance.stdout.read()
                    print "garli stdout =", stdout_lines
            

    def check_neighborhood_after_addition(self, tree, nd, edge_dist, dataset, tree_ind):
        ofprefix = "nbhood%dfromtree%d" % (edge_dist, tree_ind)
        
        self.ofprefix = ofprefix
        
        tmp_tree_filename = ".tmp.tre"
        f = open(tmp_tree_filename, "w")
        write_newick_file(f, [tree], dataset)
        f.close()

        tmp_tree_filename = ".tmpconstrain.tre"
        f = open(tmp_tree_filename, "w")
        mapper = {}
        c = copy.deepcopy(tree, mapper)
        new_nd = mapper[id(nd)]
        nd .collapse_neighborhood(edge_dist)
        write_constraint_file(f, [tree], dataset)
        f.close()
    
        self.runmode = GARLI_ENUM.NORMAL_RUNMODE
        self.streefname = ".tmp.tre"
        self.constraintfile = tmp_tree_filename

    
        
        self.run(["run"], terminate_run=True)
        
        output_tree = ofprefix + ".best.tre"
        t = dataset.read_trees(open(output_tree, "rU"), format="NEXUS")
        del dataset.trees_blocks[-1]
        sc = read_garli_scores(open(output_tree, "rU"))
        if len(t) != len(sc):
            sys.exit("Did not read the same number of trees (%d) as scores (%d) from %s" % (len(t), len(sc), output_tree))
        for otree, osc in itertools.izip(t, sc):
            otree.score = osc
        t.sort(cmp=cmp_score)
        return t
    
    def add_to_tree(self, tree, dataset, tree_ind):
        ofprefix = "from%d" % tree_ind
    
        self.ofprefix = ofprefix
        self.streefname = "incomplete"
        self.constraintfile = "none"
        
        tmp_tree_filename = ".tmp.tre"
        f = open(tmp_tree_filename, "w")
        write_tree_file(f, [tree], dataset)
        f.close()
    
        self.incompletetreefname = tmp_tree_filename
        self.runmode = GARLI_ENUM.INCR_RUNMODE
        
        self.run(["run"], terminate_run=False)
        
        output_tree = ofprefix + ".best.tre"
        t = dataset.read_trees(open(output_tree, "rU"), format="NEXUS")
        del dataset.trees_blocks[-1]
        sc = read_garli_scores(open(output_tree, "rU"))
        if len(t) != len(sc):
            sys.exit("Did not read the same number of trees (%d) as scores (%d) from %s" % (len(t), len(sc), output_tree))
        for otree, osc in itertools.izip(t, sc):
            otree.score = osc
        t.sort(cmp=cmp_score)
        return t
        

def rev_trans_func(t):
    global TAXON_TO_TRANSLATE
    return TAXON_TO_TRANSLATE[t]

def write_tree_file(outstream, trees_block, dataset):
    outstream.write("#NEXUS\nBegin Trees;\n  Translate")
    sep = ""
    for n, taxon in enumerate(dataset.taxa_blocks[0]):
        outstream.write(sep)
        sep = ',\n '
        outstream.write(" %d %s " % ((n + 1), nexus.NexusWriter.escape_token(taxon.label)))
    for tree in trees_block:
        outstream.write(";\n Tree a = [&U] %s ;\n" % tree.compose_newick(reverse_translate=rev_trans_func))
    outstream.write("End;\n")

def write_constraint_file(outstream, trees_block, dataset):
    write_newick_file(outstream, trees_block, dataset, '+')

def write_newick_file(outstream, trees_block, dataset, pref=''):
    tree = trees_block[0]
    assert(len(trees_block) == 1)
    outstream.write("%s%s;\n" % (pref, tree.compose_newick()))
    



GARLI_SCORE_PATTERN = re.compile(r"\[!GarliScore ([-0-9.]+)\]")
def read_garli_scores(inp):
    sc = []
    for line in inp:
        m = GARLI_SCORE_PATTERN.search(line)
        if m:
            g = m.group(1)
            sc.append(float(g))
    return sc


    
def read_garli_conf(f):
    default_conf = GarliConf()
    for line in f:
        s = line.split("=")
        if len(s) > 1:
            k = s[0].strip().lower()
            v = "=".join(s[1:]).strip()
            if not k in default_conf.__dict__:
                raise RuntimeError("Key %s is not understood" % k)
            setattr(default_conf, k, v)
    return default_conf


def cmp_score(x, y):
    return cmp(x.score, y.score)
    
if __name__ == '__main__':
    description =  '%s %s ' % (_program_name, _program_version)    
    usage = "%prog --conf=template.conf --data=data.nex --tree=start.tre"
    
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
                      action='store_true', 
                      dest='verbose',
                      default=False,
                      help="Verbose mode") 
    (opts, args) = parser.parse_args()
    conf_file = opts.conf
    if conf_file is None:
        sys.exit("Expecting a conf file template for GARLI")
    if opts.verbose:
        VERBOSE = True
        _LOG.setLevel(logging.DEBUG)
    data_file = opts.data_filepath
    intree_file = opts.intree_filepath
    if data_file is None:
        sys.exit("Data file must be specified")
    if intree_file is None:
        sys.exit("Input tree file must be specified")
    for f in [data_file, intree_file, conf_file]:
        if not os.path.exists(f):
            sys.exit("%s does not exist" % f)

    garli = read_garli_conf(open(conf_file, "rU"))
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


    datafname = "data.nex"
    garli.datafname = os.path.join(datafname)


    
    nexusWriter = nexus.NexusWriter()
    
    while True:
        if intree_file:
            inp_trees = dataset.read_trees(open(intree_file, "rU"), format="NEXUS")
            intree_file = ""
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

        current_taxon_mask = inp_trees[0].seed_node.edge.clade_mask
        inds = [i for i in iter_split_indices(current_taxon_mask+1)]
        assert(len(inds) == 1)
        curr_n_taxa = inds[0]
    
        curr_n_taxa += 1
        if curr_n_taxa > len(taxa):
            break

        _LOG.debug("Adding taxon %d" % curr_n_taxa)
        dirn = "t%d" % curr_n_taxa
        if not os.path.exists(dirn):
            os.makedirs(dirn)
            if not os.path.exists(dirn):
                sys.exit("Could not make %s" % dirn)
        if not os.path.isdir(dirn):
            sys.exit("%s is not a directory" % dirn)
        orig_dir = os.getcwd()
        os.chdir(dirn)
        

        culled_taxa = TaxaBlock(taxa[:curr_n_taxa])
        culled_chars = copy.copy(characters)
        culled_chars.taxa_block = culled_taxa
        culled_chars.matrix = copy.copy(characters.matrix)
        culled_chars.matrix.clear()
        #culled_chars = characters.__class__(taxa_block=culled_taxa)
        #culled_chars.column_types = characters.column_types
        #culled_chars.markup_as_sequences = characters.markup_as_sequences
        template_matrix = characters.matrix
        for taxon in culled_taxa:
            culled_chars.matrix[taxon] = template_matrix[taxon]

        culled = Dataset()
        culled.taxa_blocks.append(culled_taxa)
        culled.char_blocks.append(culled_chars)
        
        o = open(datafname, "w")
        nexusWriter.write_dataset(culled, o);
        o.close()
        
        try:
            next_round_trees = TreesBlock(taxa_block=culled_taxa)
            
            for tree_ind, tree in enumerate(inp_trees):
                trees = garli.add_to_tree(tree, culled, tree_ind)
                to_save = []
                for t in trees:
                    print t.score
                    encode_splits(t)
                    if False:
                        split = 1 << (curr_n_taxa - 1)
                        e = find_edge_from_split(t.seed_node, split)
                        if e is None:
                            sys.exit("Could not find split %s" % (bin(split)[2:]))
                            assert e is not None
                        alt_t = garli.check_neighborhood_after_addition(t, e.head_node, 2, culled, tree_ind)
                        encode_splits(alt_t[0])
                        if symmetric_difference(alt_t[0], t) != 0:
                            e = find_edge_from_split(tree.seed_node, split)
                            further_t = garli.check_neighborhood_after_addition(alt_t, e.head_node, 3, culled, tree_ind)
                            to_save.extend(further_t)
                        else:
                            to_save.append(t)
                        
                # this is where we should evaluate which trees need to be maintained for the next round.
                next_round_trees.extend(trees)
            del dataset.trees_blocks[:]
            dataset.trees_blocks.append(next_round_trees)
            o = open("incrgarli.tre", "w")
            write_tree_file(o, next_round_trees, culled)
            o.close()
            inp_trees = next_round_trees
            
        finally:
            os.chdir(orig_dir)
            
    sys.exit(0)
    
