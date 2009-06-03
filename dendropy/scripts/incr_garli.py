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
import time
import cStringIO
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
from dendropy.utils import LineReadingThread, Event
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

SLEEP_INTERVAL = 0.05

class GarliStdErrThread(LineReadingThread):
    def __init__(self, prompt, *args, **kwargs):
        LineReadingThread.__init__(self, *args, **kwargs)
        self.prompt = prompt
        self.n_prompts_read = 0
        self.n_lines_to_skip = 0
        self.igarli_comment = (0,0)
        self.prompt_pattern = re.compile(r'iGarli\[(\d+)\s*-\s*(\d+)\]')
        self.start_of_prev_read = 0
    def wait_for_prompt(self):
        self.start_of_prev_read = self.n_lines_to_skip
        while True:
            self.line_list_lock.acquire()
            try:
                to_parse = self.lines[self.n_lines_to_skip:]
            finally:
                self.line_list_lock.release()

            if to_parse:
                for n, line in enumerate(to_parse):
                    if line.startswith(self.prompt):
                        _LOG.debug("prompt found")
                        m = self.prompt_pattern.match(line)
                        if m:
                            self.n_lines_to_skip += n + 1
                            self.igarli_comment = tuple([int(i) for i in m.groups()])
                            self.n_prompts_read += 1
                            _LOG.debug("comment from prompt = %s" % str(self.igarli_comment))
                            return True
                    else:
                        _LOG.debug("%s... does not start with %s" % (line[:10], self.prompt))
                self.n_lines_to_skip += len(to_parse)
            if self.stop_event is not None:
                if self.stop_event.isSet():
                    return False
                else:
                    _LOG.debug("stop_event not triggered")
            else:
                _LOG.debug("no stop_event registered")
            _LOG.debug("Sleeping unparsed stderr = %s" % "\nline: ".join(to_parse))
            time.sleep(SLEEP_INTERVAL)
            # check to see if iGarli is still running after we slept
            self.subproc.poll()
            if self.subproc.returncode is not None:
                self.stop_event.set()
    def lines_between_prompt(self):
        self.line_list_lock.acquire()
        try:
            return self.lines[self.start_of_prev_read:self.n_lines_to_skip]
        finally:
            self.line_list_lock.release()

garli_dna_model_pattern = r'\[!GarliModel  r ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) e ([.0-9]*) ([.0-9]*) ([.0-9]*) [.0-9]* a ([.0-9]*) p ([.0-9]*) \]'

class DNAModel(object):
    def __init__(self, params):
        assert(len(params) == 10)
        self.r_mat = params[:5]
        self.r_mat.append(1.0)
        self.freqs = params[5:8]
        self.freqs.append(1 - sum(self.freqs))
        self.shape = params[8]
        self.p_invar = params[9]
    def _normalize_r_mat(self):
        assert(len(self.r_mat) == 6)
        gt_rate = self.r_mat[-1]
        nr = [i/gt_rate for i in self.r_mat[:6]]
        nr.append(1.0)
        self.r_mat = nr
            
    def __str__(self):
        self._normalize_r_mat()
        t = (self.r_mat[0], self.r_mat[1], self.r_mat[2], self.r_mat[3], self.r_mat[4],
            self.freqs[0], self.freqs[1], self.freqs[2], 1.0 - self.freqs[0] - self.freqs[1] - self.freqs[2], 
            self.shape,
            self.p_invar)
        return r'r %.5f %.5f %.5f %.5f %.5f e %.5f %.5f %.5f %.5f a %.5f p %.5f' % t

class TreeModel(object):
    def __init__(self, name="unnamed", score=None, model=None, tree=None):
        self.model = model
        self.score = score
        self.name = name
        self.tree = tree
    def __cmp__(self, other):
        return cmp(self.score, other.score)
    def __str__(self):
        return '[iGarli %s ] tree %s = [&U][!GarliScore %f][!GarliModel %s ] %s' % (self.name, self.name, self.score, str(self.model), str(self.tree))


def tree_has_structure(t):
    "returns True if there are internal nodes that are not the root"
    root = t.seed_node
    for n in root.child_nodes():
        if not n.is_leaf():
            return True
    return False

class GarliConf(object):
    def __init__(self):
        self.line_pattern = garli_line_pattern = re.compile(r'\[iGarli (\d+) \] tree best = \[&U\]\[!GarliScore ([-.0-9]*)\]%s (\(.*\))' % garli_dna_model_pattern)
        self.model_class = DNAModel
        self.active_taxa = None
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
        self.garli_stopped_event = None

    def write_garli_conf(self, out):
        conf = self.__dict__
        out.write("[general]\n")
        for k in GARLI_GENERAL:
            out.write("%s = %s\n" % (k, conf[k]))
        out.write("[master]\n")
        for k in GARLI_MASTER:
            out.write("%s = %s\n" % (k, conf[k]))

    def set_active_taxa(self, active_taxa):
        self.active_taxa = active_taxa

    def quit(self):
        return self.run(["quit"], terminate_run=True);

    def run(self, commands, terminate_run=True):
        if self.garli_instance is None:
            assert(self.garli_stopped_event is None)
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
            e = Event()
            self.garli_stopped_event = e
            self.stderrThread = GarliStdErrThread(prompt="iGarli[", stream=s.stderr, store_lines=True, subproc=s, stop_event=e)
            self.stderrThread.start()
            self.stdoutThread = LineReadingThread(stream=s.stdout, store_lines=True, subproc=s, stop_event=e)
            self.stdoutThread.start()

        assert(self.stderrThread is not None)
        for command in commands:
            self.garli_instance.poll()
            if self.garli_instance.returncode is not None:
                _LOG.debug("triggering stop_event in command loop")
                self.garli_stopped_event.set()
            if self.stderrThread.wait_for_prompt():
                _LOG.debug("***ISSUING command ***\n%s\n" % command) 
                self.garli_instance.stdin.write(command + '\n')
            else:
                assert False, 'iGarli exited before I even had a chance to tell it "%s".  How frustrating for me!\nThe full err stream is:\n%s' % (command, "".join(self.stderrThread.lines))

        if terminate_run:
            rc = self.garli_instance.wait()
            if rc != 0:
                sys.exit(gstderr)
        else:
            self.garli_instance.poll()
        if self.garli_instance.returncode is not None:
            _LOG.debug("triggering stop_event after commands")
            self.garli_stopped_event.set()
            self.garli_stopped_event = None
            self.garli_instance = None

    def check_neighborhood_after_addition(self, tree_model, nd, edge_dist, dataset, tree_ind):
        tree = tree_model.tree
        ofprefix = "nbhood%dfromtree%d" % (edge_dist, tree_ind)
        
        self.ofprefix = ofprefix
        
        tmp_tree_filename = ".tmp.tre"
        f = open(tmp_tree_filename, "w")
        write_tree_file(f, [tree], dataset)
        f.close()
        # it seems a little odd to call this incompletetreefname rather than streefname
        #   but we'd like to trigger the interactive mode, and this is one way of doing that.
        self.incompletetreefname = tmp_tree_filename

        tmp_tree_filename = ".tmpconstrain.tre"
        f = open(tmp_tree_filename, "w")
        mapper = {}
        c = copy.deepcopy(tree, mapper)
        new_nd = mapper[id(nd)]
        nd .collapse_neighborhood(edge_dist)
        _LOG.debug("Checking neighborhood by constraining %s" % str(tree))
        write_constraint_file(f, [tree], dataset)
        f.close()
    
        self.runmode = GARLI_ENUM.INCR_RUNMODE # NORMAL_RUNMODE

        #garli does not like stars as constraints
        do_constraint = tree_has_structure(tree)
        if do_constraint:
            # cache settings
            pcf = self.constraintfile
            self.constraintfile = tmp_tree_filename
        
        self.run(["model = %s" % str(tree_model.model),
                  "run"], terminate_run=True)
        if do_constraint:
            # restore cached settings
            self.constraintfile = pcf 
        
        err_lines = self.stderrThread.lines_between_prompt()
        r = self.parse_igarli_lines(err_lines, dataset)
        r.sort(reverse=True)
        return r


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

        # store the settings that we must override
        ptw = self.topoweight
        psg = self.stopgen
        self.topoweight = 0.0
        self.stopgen = 1000
        
        self.run(["run"], terminate_run=True)

        # restore the settings that we had to override
        self.topoweight = ptw
        self.stopgen = psg
        
        err_lines = self.stderrThread.lines_between_prompt()
        r = self.parse_igarli_lines(err_lines, dataset)
        r.sort(reverse=True)
        return r
    
    def parse_igarli_lines(self, line_list, dataset):
        tm_list = []
        for line in line_list:
            if line.startswith('[iGarli '):
                m = self.line_pattern.match(line)
                assert(m)
                g = m.groups()
                result_number = g[0]
                score = float(g[1])
                model_params = [float(i) for i in g[2:-1]]
                model = self.model_class(model_params)
                tree_string = g[-1]
                newick_stream = cStringIO.StringIO(tree_string)
                tree_list = self.read_trees(dataset, newick_stream, format="newick")
                assert(len(tree_list) == 1)
                tm = TreeModel(name=result_number, score=score, model=model, tree=tree_list[0])
                tm_list.append(tm)
        return tm_list

    def read_trees(self, dataset, stream, format):
        if (format.upper() == "NEWICK") and (self.active_taxa):
            # we will add a translate dictionary using self.active_taxa
            td = {}
            assert len(dataset.taxa_blocks) == 1
            assert dataset.taxa_blocks[0] is self.active_taxa, "%s != %s" % (str(dataset.taxa_blocks[0]), str(self.active_taxa))
            tb = self.active_taxa
            for n, taxon in enumerate(self.active_taxa):
                td[str(1 + n)] = taxon
            t = dataset.read_trees(stream, format=format, translate_dict=td)
        else:
            t = dataset.read_trees(stream, format=format)
        del dataset.trees_blocks[-1]
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
    full_dataset = Dataset()
    full_dataset.read(open(data_file, "rU"), format="NEXUS")
    taxa = full_dataset.taxa_blocks[0]
    full_taxa_mask = taxa.all_taxa_bitmask()
    for n, taxon in enumerate(taxa):
        TAXON_TO_TRANSLATE[taxon] = str(n + 1)
    _LOG.debug("%s = full_taxa_mask" % bin(full_taxa_mask))
    assert(len(full_dataset.taxa_blocks) == 1)
    characters = full_dataset.char_blocks[0]
    assert(len(full_dataset.char_blocks) == 1)
    assert(len(characters) == len(taxa))


    datafname = "data.nex"
    garli.datafname = os.path.join(datafname)


    
    nexusWriter = nexus.NexusWriter()
    
    while True:
        if intree_file:
            inp_trees = full_dataset.read_trees(open(intree_file, "rU"), format="NEXUS")
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

        garli.set_active_taxa(culled_taxa)

        try:
            next_round_trees = [TreeModel(tree=i) for i in TreesBlock(taxa_block=culled_taxa)]
            
            for tree_ind, tree in enumerate(inp_trees):
                tree_model_list = garli.add_to_tree(tree, culled, tree_ind)
                to_save = []
                for tm in tree_model_list:
                    print tm.score
                    step_add_tree = tm.tree
                    encode_splits(step_add_tree)
                    split = 1 << (curr_n_taxa - 1)
                    e = find_edge_from_split(step_add_tree.seed_node, split)
                    assert e is not None, "Could not find split %s.  Root mask is %s" % (bin(split)[2:], bin(step_add_tree.seed_node.edge.clade_mask)[2:])
                    
                    garli.topoweight = 1.0
                    garli.stopgen = 1000

                    nt_list = garli.check_neighborhood_after_addition(tm, e.head_node, 2, culled, tree_ind)
                    deeper_search_start = []
                    better_tm = tm
                    for nt in nt_list:
                        encode_splits(nt.tree)
                        if symmetric_difference(nt.tree, step_add_tree) != 0:
                            deeper_search_start.append(nt)
                        elif nt.score > better_tm.score:
                            better_tm = nt


                    if deeper_search_start:
                        entire_neighborhood = [better_tm] + deeper_search_start
                        for alt_tm in deeper_search_start:
                            e = find_edge_from_split(alt_tm.tree.seed_node, split)

                            assert e is not None, "Could not find split %s.  Root mask is %s" % (bin(split)[2:], bin(alt_tm.tree.seed_node.edge.clade_mask)[2:])

                            nt_list = garli.check_neighborhood_after_addition(alt_tm, e.head_node, 3, culled, tree_ind)
                            for nt in nt_list:
                                encode_splits(nt.tree)
                                entire_neighborhood.append(nt)
                        entire_neighborhood.sort(reverse=True)
                        to_add = []
                        for nt in entire_neighborhood:
                            found = False
                            for x in to_add:
                                if symmetric_difference(x.tree, nt.tree) == 0:
                                    found = True
                                    break
                            if not found:
                                to_add.append(nt)
                        to_save.extend(to_add)
                    else:
                        to_save.append(better_tm)
                        
                # this is where we should evaluate which trees need to be maintained for the next round.
                next_round_trees.extend(to_save)
            del full_dataset.trees_blocks[:]
            inp_trees = [i.tree for i in next_round_trees]
            full_dataset.trees_blocks.append(inp_trees)
            o = open("incrgarli.tre", "w")
            write_tree_file(o, [i.tree for i in next_round_trees], culled)
            o.close()
            
        finally:
            os.chdir(orig_dir)
            
    sys.exit(0)
    
