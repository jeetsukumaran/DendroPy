#! /usr/bin/env python

############################################################################
##  incr_garli.py
##
##  Copyright 2009 Mark T. Holder
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
from dendropy.splits import encode_splits, lowest_bit_only, iter_split_indices, find_edge_from_split, collapse_conflicting, is_trivial_split, split_as_string_rev
from dendropy.characters import CharactersBlock
from dendropy.taxa import TaxaBlock
from dendropy.treedists import symmetric_difference
from dendropy import treesum
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
GARLI_ALL = tuple( list(GARLI_GENERAL) + list(GARLI_MASTER) )

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
        assert len(self.r_mat) == 6, "%s is an invalid r_mat" % str(self.r_mat)
        gt_rate = self.r_mat[-1]
        nr = [i/gt_rate for i in self.r_mat[:5]]
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
            _LOG.debug("n children = %s ; newick = %s " % (str(n.child_nodes()), t.compose_newick(reverse_translate=rev_trans_func)))
            return True
    return False

class GarliConf(object):
    def __init__(self):
        # the following are settings for the incremental behavior
        self.last_model_read = None
        self.init_tree_scoring_stopgen = 200
        self.init_stopgen = 100
        
        self.tree_scoring_modweight = 0.05
        
        self.add_tree_modweight = 0.0
        
        self.neighborhood_topoweight = 0.5
        self.neighborhood_modweight = 0.0
        self.neighborhood_stopgen = 100
        
        self.negcon_topoweight = 0.5
        self.negcon_modweight = 0.01
        self.negcon_stopgen = 100

        self.first_neighborhood = 3
        self.neighborhood_incr = 1
        self.curr_n_taxa = None
        
        self.split_diversity_multiplier = 1.0 # the number multiplied by the "split_diversity_score" when evaluating which trees to carry over
                                              # higher values mean that trees with lower likelihoods may score better.
        self._max_trees_carried_over = 1000 # the maximum number of trees to carry over between iterations
        
        # The next section are GARLI settings
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
        self._modweight =  "0.05"
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
        # finally we have attributes that are impl details
        self.active_taxa = None
        self.garli_instance = None
        self.garli_stopped_event = None
        self._cached_settings = []
        self.line_pattern = re.compile(r'\[iGarli (\d+) \] tree best = \[&U\]\[!GarliScore ([-.0-9]*)\]%s (\(.*\))' % garli_dna_model_pattern)
        self.no_model_line_pattern = re.compile(r'\[iGarli (\d+) \] tree best = \[&U\]\[!GarliScore ([-.0-9]*)\] (\(.*\))')
        self.model_class = DNAModel

    def get_max_trees_carried_over(self):
        return self._max_trees_carried_over

    max_trees_carried_over = property(get_max_trees_carried_over)
    
    def cache_settings(self):
        d = {}
        for k in GARLI_ALL:
            d[k] = getattr(self, k)
        self._cached_settings.append(d)

    def restore_settings(self):
        d = self._cached_settings.pop()
        for k in GARLI_ALL:
            setattr(self, k, d[k])

    def write_garli_conf(self, out):
        out.write("[general]\n")
        for k in GARLI_GENERAL:
            out.write("%s = %s\n" % (k, getattr(self, k)))
        out.write("[master]\n")
        for k in GARLI_MASTER:
            out.write("%s = %s\n" % (k, getattr(self, k)))

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

    def set_modweight(self, x):
        if x <= 0.0:
            self.statefrequencies = "fixed"
            self.ratematrix = "fixed"
            self.ratehetmodel = "gammafixed"
            self.invariantsites = "fixed"
            self._modweight = 0.0
        else:
            self.statefrequencies = "estimate"
            self.ratematrix = "6rate"
            self.ratehetmodel = "gamma"
            self.invariantsites = "estimate"
            self._modweight = x
    def get_modweight(self):
        return self._modweight
    
    modweight = property(get_modweight, set_modweight)
    
    def find_best_conflicting(self, starting_tree, split, dataset):
        new_starting = TreeModel(model=starting_tree.model)
        new_starting.tree = copy.deepcopy(starting_tree.tree)
        root = new_starting.tree.seed_node
        e = find_edge_from_split(root, split, root.edge.clade_mask)
        if e:
            e.collapse()

        tmp_tree_filename = ".tmp.tre"
        write_trees_to_filepath([new_starting], dataset, tmp_tree_filename)
        self.cache_settings()

        try:
    
            tmp_constrain_filename = ".tmpconstrain.tre"
            self.constraintfile = tmp_constrain_filename
            f = open(tmp_constrain_filename, "w")
            f.write("-%s\n" % split_as_string_rev(split, self.curr_n_taxa, '.', '*'))
            f.close()
    
            self.ofprefix = "negconst%d" % (split)
            
            # it seems a little odd to call this incompletetreefname rather 
            #   than streefname but we'd like to trigger the interactive mode, 
            #   and this is one way of doing that.
            self.incompletetreefname = tmp_tree_filename
            self.runmode = GARLI_ENUM.INCR_RUNMODE # NORMAL_RUNMODE
            self.topoweight = self.negcon_topoweight
            self.modweight = self.negcon_modweight
            self.stopgen = self.negcon_stopgen
            invoc = []
            if new_starting.model:
                invoc.append("model = %s" % str(new_starting.model))
            invoc.append("run")
            self.run(invoc, terminate_run=True)
        finally:
            self.restore_settings()

        err_lines = self.stderrThread.lines_between_prompt()
        r = self.parse_igarli_lines(err_lines, dataset)
        r.sort(reverse=True)
        for tm in r:
            encode_splits(tm.tree)
            assert split not in tm.tree.split_edges
        return r
       
    def check_neighborhood_after_addition(self, tree_model, nd, edge_dist, dataset, tree_ind):
        tmp_tree_filename = ".tmp.tre"
        write_trees_to_filepath([tree_model], dataset, tmp_tree_filename)
        self.cache_settings()
        tree = tree_model.tree

        # create a less resolved version to serve as a constraint
        mapper = {}
        c = copy.deepcopy(tree, mapper)
        new_nd = mapper[id(nd)]
        new_nd.collapse_neighborhood(edge_dist)
        _LOG.debug("Checking neighborhood by constraining:\n%s\nwhich was based on collapse_neighborhood %d around %s from:\n%s" % (str(c), edge_dist, new_nd.taxon, str(tree)))

        try:
    
            #garli does not like stars as constraints
            if tree_has_structure(c):
                tmp_constrain_filename = ".tmpconstrain.tre"
                self.constraintfile = tmp_constrain_filename
                f = open(tmp_constrain_filename, "w")
                write_constraint_file(f, [c], dataset)
                f.close()
    
            self.ofprefix = "nbhood%dfromtree%d" % (edge_dist, tree_ind)
            
            # it seems a little odd to call this incompletetreefname rather 
            #   than streefname but we'd like to trigger the interactive mode, 
            #   and this is one way of doing that.
            self.incompletetreefname = tmp_tree_filename
            self.runmode = GARLI_ENUM.INCR_RUNMODE # NORMAL_RUNMODE
            self.topoweight = self.neighborhood_topoweight
            self.modweight = self.neighborhood_modweight
            self.stopgen = self.neighborhood_stopgen
            invoc = []
            if tree_model.model:
                invoc.append("model = %s" % str(tree_model.model))
            elif self.modweight <= 0.0:
                if self.last_model_read:
                    invoc.append("model = %s" % str(self.last_model_read))
                else:
                    self.modweight = 0.01
            invoc.append("run")
            self.run(invoc, terminate_run=True)
        finally:
            self.restore_settings()

        err_lines = self.stderrThread.lines_between_prompt()
        r = self.parse_igarli_lines(err_lines, dataset)
        r.sort(reverse=True)
        return r


    def add_to_tree(self, tree_model, dataset, tree_ind, stop_gen):
        tmp_tree_filename = ".tmp.tre"
        write_trees_to_filepath([tree_model], dataset, tmp_tree_filename)
        self.cache_settings()
        try:
            self.ofprefix = "from%d" % tree_ind
            self.streefname = "incomplete"
            self.constraintfile = "none"
            self.topoweight = 0.0
            self.modweight = self.add_tree_modweight
            self.stopgen = stop_gen
            self.incompletetreefname = tmp_tree_filename
            self.runmode = GARLI_ENUM.INCR_RUNMODE
            invoc = ["keep = 5"]
            if tree_model.model:
                invoc.append("model = %s" % str(tree_model.model))
            invoc.append("run")
            self.run(invoc, terminate_run=True)
        finally:
            self.restore_settings()

        err_lines = self.stderrThread.lines_between_prompt()
        r = self.parse_igarli_lines(err_lines, dataset)
        r.sort(reverse=True)
        return r


    def score_tree(self, tree, dataset, tree_ind, stop_gen):

        tmp_tree_filename = ".tmp.tre"
        write_trees_to_filepath([tree], dataset, tmp_tree_filename)
        self.cache_settings()
        try:
            self.ofprefix = "from%d" % tree_ind
            self.streefname = "incomplete"
            self.constraintfile = "none"
            self.topoweight = 0.0
            self.modweight = self.tree_scoring_modweight
            self.stopgen = stop_gen
            self.incompletetreefname = tmp_tree_filename
            self.runmode = GARLI_ENUM.INCR_RUNMODE
    
            self.run(["run"], terminate_run=True)
        finally:
            self.restore_settings()

        err_lines = self.stderrThread.lines_between_prompt()
        r = self.parse_igarli_lines(err_lines, dataset)
        r.sort(reverse=True)
        return r[0]

    def parse_igarli_lines(self, line_list, dataset):
        tm_list = []
        for line in line_list:
            if line.startswith('[iGarli '):
                has_model = True
                m = self.line_pattern.match(line)
                if not m:
                    m = self.no_model_line_pattern.match(line)
                    has_model = False
                    assert m, "%s doesn't fit the pattern" % line
                g = m.groups()
                result_number = g[0]
                score = float(g[1])
                if has_model:
                    model_params = [float(i) for i in g[2:-1]]
                    model = self.model_class(model_params)
                    self.last_model_read = model
                else:
                    model = None
                tree_string = g[-1]
                newick_stream = cStringIO.StringIO(tree_string)
                tree_list = self.read_trees(dataset, newick_stream, format="newick")
                assert(len(tree_list) == 1)
                tm = TreeModel(name="tree" + result_number, score=score, model=model, tree=tree_list[0])
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

    def incrementally_build_trees(self, full_dataset, inp_trees):
        assert len(full_dataset.taxa_blocks) == 1
        taxa = full_dataset.taxa_blocks[0]
        assert len(inp_trees) > 0

        current_taxon_mask = inp_trees[0].tree.seed_node.edge.clade_mask
        inds = [i for i in iter_split_indices(current_taxon_mask + 1)]
        assert(len(inds) == 1)
        self.curr_n_taxa = inds[0]
        # reset the add_tree_stopgen to its initial setting
        self.add_tree_stopgen = self.init_stopgen
        self.tree_scoring_stop_gen = self.init_tree_scoring_stopgen
        # Score the original trees
        orig_dir = self._cd_for_work_dir()
        try:
            culled = self._write_garli_input(full_dataset)
            culled_taxa = culled.taxa_blocks[0]
            self.set_active_taxa(culled_taxa)
            inp_trees = self.score_tree_list(full_dataset, inp_trees, stop_gen=self.tree_scoring_stop_gen)
        finally:
            os.chdir(orig_dir)
        self.curr_n_taxa += 1
        
        
        while self.curr_n_taxa <= len(taxa):
            _LOG.debug("Adding taxon %d" % self.curr_n_taxa)            
            orig_dir = self._cd_for_work_dir()
            try:
                inp_trees = self._do_add_taxon_incremental_step(full_dataset, inp_trees)
            finally:
                os.chdir(orig_dir)
            self.curr_n_taxa += 1

    def _cd_for_work_dir(self):
        """changes the current directory to a working directory for the specified
        number of taxa and returns the absolute path to the previous directory."""
        dirn = "t%d" % self.curr_n_taxa
        if not os.path.exists(dirn):
            os.makedirs(dirn)
            if not os.path.exists(dirn):
                sys.exit("Could not make %s" % dirn)
        if not os.path.isdir(dirn):
            sys.exit("%s is not a directory" % dirn)

        orig_dir = os.getcwd()
        os.chdir(dirn)
        return orig_dir

    def _write_garli_input(self, full_dataset):
        assert len(full_dataset.taxa_blocks) == 1
        taxa = full_dataset.taxa_blocks[0]

        assert(len(full_dataset.char_blocks) == 1)
        characters = full_dataset.char_blocks[0]
        assert(len(characters) == len(taxa))
        
        culled_taxa = TaxaBlock(taxa[:self.curr_n_taxa])
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

        o = open(self.datafname, "w")
        nexusWriter = nexus.NexusWriter()
        nexusWriter.write_dataset(culled, o);
        o.close()
        return culled


    def select_trees_for_next_round(self, culled, curr_results):
        all_taxa_bitmask = curr_results[0].tree.seed_node.edge.clade_mask
        assert all_taxa_bitmask == ((1 << self.curr_n_taxa) - 1)
        ########################################
        # First, we make sure that there are not duplicate topologies
        # Because we reverse sort, we'll be retaining the tree with the
        #   best score
        #####
        curr_results.sort(reverse=True)
        set_of_split_sets = set()
        unique_topos = []
        for tm in curr_results:
            add_nontriv_splits_attr(tm, all_taxa_bitmask)
            if tm.splits not in set_of_split_sets:
                unique_topos.append(tm)
                set_of_split_sets.add(tm.splits)
        curr_results = unique_topos
        set_of_split_sets.clear()
        _LOG.info('There were %d unique result topologiesfor ntax = %d ' % (len(curr_results), self.curr_n_taxa))

        ########################################
        # the trees can be hefty, so lets eliminate unneeded references
        #####
        del unique_topos


        ########################################
        # Make sure to keep the next_round_trees list (and a set that helps us
        #   identify unique trees.
        #####
        ml_est = curr_results[0]
        curr_results.pop(0)
        next_round_trees = [ml_est]
        
        ########################################
        # Now we identify best trees that LACK the splits in the ML tree
        #####
        ml_split_dict = ml_est.tree.split_edges
        unanimous_splits = []
        best_disagreeing_index_set = set()
        for split in ml_split_dict.iterkeys():
            if not is_trivial_split(split, all_taxa_bitmask):
                found = False
                for n, tm in enumerate(curr_results):
                    if split not in tm.tree.split_edges:
                        best_disagreeing_index_set.add(n)
                        found = True
                        break
                if not found:
                    unanimous_splits.append(split)

        ########################################
        # now we add the trees that "must" be included because they do NOT have
        #   a split that is in the current ML tree.
        # We do this with a reverse sorted list so that we can pop them off of
        #   the curr_results list without invalidating the list of indices to move
        #####
        bdis_list = list(best_disagreeing_index_set)
        bdis_list.sort(reverse=True)
        for tree_ind in bdis_list:
            tm = curr_results.pop(tree_ind)
            next_round_trees.append(tm)
        
        ########################################
        # Now we have to augment our list of trees such that we have exemplar trees
        #   that conflict with every split in the ML tree
        # We'll do this by starting from a version of the ML tree that has been
        #   collapsed so that it does not conflict with the split
        #####
        _LOG.info('There were %d unanimous splits in the curr_results for ntax = %d ' % (len(unanimous_splits), self.curr_n_taxa))
        for split in unanimous_splits:
            best_conflicting = self.find_best_conflicting(starting_tree=ml_est, split=split, dataset=culled)
            for b in best_conflicting:
                add_nontriv_splits_attr(b, all_taxa_bitmask)

            best_conflicting.sort(reverse=True)
            tm = best_conflicting[0]
            next_round_trees.append(tm)
            curr_results.extend(best_conflicting[1:])
        
        
        ########################################
        # To keep the remaining trees in next_round_trees diverse we will
        #   try to add trees that maximize a score which is:
        #       lambda*tree_split_rarity + lnL
        #   where lambda is a tuning parameter and tree_split_rarity is:
        #       n_tree_times_splits = num_trees_in_next_round_trees * num_splits_per_tree
        #       split_occurrence = num_trees_in_next_round_trees_that_have_split
        #       tree_split_rarity = n_tree_times_splits - SUM split_occurrence
        #   in which the summation is taken over all splits in the tree
        #####
        max_len = self.max_trees_carried_over
        num_trees_to_add = max_len - len(next_round_trees)
        if num_trees_to_add > len(curr_results):
            split_count = {}
            n_tree_times_splits = 0
            for tm in next_round_trees:
                for k in tm.splits:
                    split_count[k] = split_count.get(k, 0) + 1
                    n_tree_times_splits += 1
            for tm in curr_results:
                tm.tree_split_rarity = n_tree_times_splits
                for split in tm.splits:
                    tm.tree_split_rarity -= split_count.get(split, 0)
            def split_diversity_cmp(x, y, lambda_mult=self.split_diversity_multiplier):
                x.retention_score = lambda_mult*x.tree_split_rarity + x.score
                y.retention_score = lambda_mult*y.tree_split_rarity + y.score
                return cmp(x.retention_score, y.retention_score)
            curr_results.sort(cmp=split_diversity_cmp, reverse=True)

        ########################################
        # We are now going to try to add elements (in order) from curr_results
        #   until we run out of trees to add or we reach max_len
        #####
        n_added = 0
        try:
            set_of_split_sets.clear()
            for tm in next_round_trees:
                set_of_split_sets.add(tm.splits)

            cri = iter(curr_results)
            while n_added < num_trees_to_add:
                tm = cri.next()
                if tm.splits not in set_of_split_sets:
                    set_of_split_sets.add(tm.splits)
                    next_round_trees.append(tm)
                    n_added += 1
        except StopIteration:
            pass
        _LOG.info('Added %d trees that were not "required" to guarantee that no splits were unanimous for ntax = %d' % (n_added, self.curr_n_taxa))

        ########################################
        # the trees can be hefty, so lets free unneeded memory
        #####
        del curr_results[:]
        
        ########################################
        # Finally, lets get a decent score for each tree before moving to the next round
        #   because the trees are big, we'll replace each element rather
        #   than allowing a duplicate list to be created.
        #####
        for i in range(len(next_round_trees)):
            next_round_trees[i] = self.score_tree(next_round_trees[i], culled, n, self.tree_scoring_stop_gen)

        _LOG.info('A total of %d trees were retained for ntax = %d lnL range from %f to %f' % (len(next_round_trees), self.curr_n_taxa, next_round_trees[0].score, next_round_trees[-1].score))

        return next_round_trees

    def score_tree_list(self, full_dataset, inp_trees, stop_gen):
        culled = self._write_garli_input(full_dataset)
        culled_taxa = culled.taxa_blocks[0]
        self.set_active_taxa(culled_taxa)
        rescored = []
        for tree_ind, tree in enumerate(inp_trees):
            tm = self.score_tree(tree, culled, tree_ind, stop_gen=stop_gen)
            encode_splits(tm.tree)
            rescored.append(tm)
        rescored.sort(reverse=True)
        del full_dataset.trees_blocks[:]
        full_dataset.trees_blocks.append([i.tree for i in rescored])
        o = open("incrgarli.tre", "w")
        write_tree_file(o, rescored, culled)
        o.close()
        return rescored

    def _do_add_taxon_incremental_step(self, full_dataset, inp_trees):
        culled = self._write_garli_input(full_dataset)
        culled_taxa = culled.taxa_blocks[0]
        self.set_active_taxa(culled_taxa)
        next_round_trees = []

        for tree_ind, tree in enumerate(inp_trees):
            tree_model_list = self.add_to_tree(tree, culled, tree_ind, self.add_tree_stopgen)
            to_save = []
            for tm in tree_model_list:
                print "Tree %d for %d taxa: %f" % (tree_ind, self.curr_n_taxa, tm.score)
                step_add_tree = tm.tree
                encode_splits(step_add_tree)
                split = 1 << (self.curr_n_taxa - 1)
                e = find_edge_from_split(step_add_tree.seed_node, split)
                assert e is not None, "Could not find split %s.  Root mask is %s" % (bin(split)[2:], bin(step_add_tree.seed_node.edge.clade_mask)[2:])


                nt_list = self.check_neighborhood_after_addition(tm, e.head_node, self.first_neighborhood, culled, tree_ind)
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

                        nt_list = self.check_neighborhood_after_addition(alt_tm, e.head_node, self.first_neighborhood + self.neighborhood_incr, culled, tree_ind)
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
        
        next_round_trees = self.select_trees_for_next_round(culled, next_round_trees)
        
        del full_dataset.trees_blocks[:]
        full_dataset.trees_blocks.append([i.tree for i in next_round_trees])
        o = open("incrgarli.tre", "w")
        write_tree_file(o, next_round_trees, culled)
        o.close()
        return next_round_trees

    def read_garli_conf(self, stream):
        for line in stream:
            s = line.split("=")
            if len(s) > 1:
                k = s[0].strip().lower()
                v = "=".join(s[1:]).strip()
                try:
                    getattr(self, k)
                except:
                    raise RuntimeError("Key %s is not understood" % k)
                setattr(self, k, v)

def rev_trans_func(t):
    global TAXON_TO_TRANSLATE
    return TAXON_TO_TRANSLATE[t]

def write_trees_to_filepath(tree_list, dataset, filepath):
        tmp_tree_filename = ".tmp.tre"
        f = open(filepath, "w")
        write_tree_file(f, tree_list, dataset)
        f.close()

_tm_template = ";\n Tree %s = [&U][!GarliScore %f][!GarliModel %s ] %s ;\n"
_tm_unscored_template = ";\n Tree %s = [&U][!GarliModel %s ] %s ;\n"
_tree_template = ";\n Tree a = [&U] %s ;\n"

def write_tree_file(outstream, tree_model_list, dataset):
    outstream.write("#NEXUS\nBegin Trees;\n  Translate")
    sep = ""
    for n, taxon in enumerate(dataset.taxa_blocks[0]):
        outstream.write(sep)
        sep = ',\n '
        outstream.write(" %d %s " % ((n + 1), nexus.NexusWriter.escape_token(taxon.label)))
    for tm in tree_model_list:
        try:
            tree = tm.tree
        except:
            newick = tm.compose_newick(reverse_translate=rev_trans_func)
            msg = _tree_template % newick
        else:
            newick = tree.compose_newick(reverse_translate=rev_trans_func)
            if tm.model:
                if tm.score:
                    msg = _tm_template % (tm.name, tm.score, str(tm.model), newick)
                else:
                    msg = _tm_unscored_template % (tm.name, str(tm.model), newick)
            else:
                msg = _tree_template % newick
        outstream.write(msg)

    outstream.write("End;\n")

def write_constraint_file(outstream, trees_block, dataset):
    write_newick_file(outstream, trees_block, dataset, '+')

def write_newick_file(outstream, trees_block, dataset, pref=''):
    tree = trees_block[0]
    assert(len(trees_block) == 1)
    outstream.write("%s%s;\n" % (pref, tree.compose_newick()))

def add_nontriv_splits_attr(tm, all_taxa_bitmask):
    all_spl = tm.tree.split_edges.keys()
    non_triv = [i for i in all_spl if not is_trivial_split(i, all_taxa_bitmask)]
    non_triv.sort()
    tm.splits = tuple(non_triv)
    


GARLI_SCORE_PATTERN = re.compile(r"\[!GarliScore ([-0-9.]+)\]")
def read_garli_scores(inp):
    sc = []
    for line in inp:
        m = GARLI_SCORE_PATTERN.search(line)
        if m:
            g = m.group(1)
            sc.append(float(g))
    return sc


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

    garli = GarliConf()
    garli.read_garli_conf(open(conf_file, "rU"))

    full_dataset = Dataset()
    full_dataset.read(open(data_file, "rU"), format="NEXUS")
    assert(len(full_dataset.taxa_blocks) == 1)
    taxa = full_dataset.taxa_blocks[0]
    full_taxa_mask = taxa.all_taxa_bitmask()
    for n, taxon in enumerate(taxa):
        TAXON_TO_TRANSLATE[taxon] = str(n + 1)
    _LOG.debug("%s = full_taxa_mask" % bin(full_taxa_mask))


    garli.datafname = os.path.join("data.nex")

    raw_trees = full_dataset.read_trees(open(intree_file, "rU"), format="NEXUS")
    assert(raw_trees)
    current_taxon_mask = None

    # read initial trees and verify that they have the correct set of taxa
    for tree in raw_trees:
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

    inp_trees = [TreeModel(tree=i) for i in raw_trees]

    garli.incrementally_build_trees(full_dataset, inp_trees)

    sys.exit(0)

