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


import sys
import re
import copy
from cStringIO import StringIO
from dendropy.taxa import TaxaBlock
from dendropy.datasets import Dataset
from dendropy.splits import encode_splits, is_trivial_split, find_edge_from_split, SplitDistribution
from dendropy.treesum import TreeSummarizer
from dendropy import get_logger
_LOG = get_logger("igarli_neighborhood.py")



#garli_dna_model_pattern = r'\[!GarliModel  r ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) e ([.0-9]*) ([.0-9]*) ([.0-9]*) [.0-9]* a ([.0-9]*) p ([.0-9]*) \]'
garli_dna_model_pattern = r'\[!GarliModel\s*(r [.0-9]* [.0-9]* [.0-9]* [.0-9]* [.0-9]* e [.0-9]* [.0-9]* [.0-9]* [.0-9]* a [.0-9]* p [.0-9]*)\s*\]'
garli_after_name_tm_pattern = r'\s*=\s*\[&U\]\[!GarliScore ([-.0-9]*)\]\s*%s\s*(\(.*\))\s*;' % garli_dna_model_pattern
garli_after_name_tree_pattern = r'\s*=\s*\[&U\]\[!GarliScore ([-.0-9]*)\]\s*(\(.*\))\s*;'



tree_prefix = r'\[iGarli\s*(\d+)\s*\] tree [a-zA-Z0-9_]+'
tm_pat_string = tree_prefix + garli_after_name_tm_pattern
tree_pat_string = tree_prefix + garli_after_name_tree_pattern
garli_tree_model_pat = re.compile(tm_pat_string)
garli_tree_pat = re.compile(tree_pat_string)


class ParsedTree(object):
    def __init__(self, score=None, model=None, tree_string=None, tree=None):
        self.score = score
        self.model = model
        self.tree_string = tree_string
        self.tree = tree
    def __cmp__(self, x):
        return cmp(self.score, x.score)

def add_nontriv_splits_attr(tree, all_taxa_bitmask):
    all_spl = tree.split_edges.keys()
    non_triv = []
    for i in all_spl:
        if not is_trivial_split(i, all_taxa_bitmask):
            if i & 1:
                non_triv.append(i)
            else:
                non_triv.append((~i)&all_taxa_bitmask)
    non_triv.sort()
    tree.splits = tuple(non_triv)
    tree.split_set = set(non_triv)

    
def read_add_tree_groups(f):
    first = True
    curr_tree_group = []
    all_tree_groups = []
    # this is hacky -- we use the presence of a model string to tell us 
    #   when we have hit a new group of trees that were found by adding to 
    #   the same incomplete tree (the other trees will come from the suboptimal
    #   trees that are spit out without any model info because the hook
    #   in the Finish...Stepwise...() function does not expose the model)
    curr_model = None
    for line in f:
        m = garli_tree_model_pat.match(line)
        if m:
            if first:
                first = False
                assert not curr_tree_group
            else:
                assert curr_tree_group
                for el in curr_tree_group:
                    if el.model is None:
                        el.model = curr_model
                all_tree_groups.append(curr_tree_group)
                curr_tree_group = []
                curr_model = None
            has_model = True
        else:
            m = garli_tree_pat.match(line)
            has_model = False
        if m:
            g = m.groups()
            result_number = g[0]
            score = float(g[1])
            if has_model:
                curr_model =  g[2]
            curr_tree_group.append(ParsedTree(score=score, model=None, tree_string=g[-1]))
        else:
            _LOG.debug("No match: %s" % line)
    
    
    if curr_model:
        for el in curr_tree_group:
            if el.model is None:
                el.model = curr_model
        all_tree_groups.append(curr_tree_group)
    return all_tree_groups

def get_norm_nontrivial_split_set(tree):
    add_nontriv_splits_attr(tree, all_taxa_bitmask)
    return tree.split_set

def connected_at(dist_mat, max_dist):
    _LOG.debug("dist_mat =%s\n" % "\n".join([str(row) for row in dist_mat]))

    dim = len(dist_mat)
    if dim == 0:
        return []
    assert dim == len(dist_mat[0])
    svec = [None]*dim
    svec[0] = set([0])
    for row_n, row in enumerate(dist_mat[:-1]):
        offset = 1 + row_n
        for col_disp, d in enumerate(row[offset:]):
            col_n = offset + col_disp
            if d < max_dist:
                rs = svec[row_n]
                cs = svec[col_n]
                if rs:
                    if cs:
                        if rs is not cs:
                            #merge rs and cs
                            for ind in cs:
                                svec[ind] = rs
                                rs.add(ind)
                    else:
                        rs.add(col_n)
                        svec[col_n] = rs
                else:
                    if cs is None:
                        ns = set([col_n, row_n])
                        svec[row_n] = ns
                        svec[col_n] = ns
                    else:
                        svec[row_n] = cs
                        cs.add(row_n)
    written = []
    all_connected_inds = []
    for n, el in enumerate(svec):
        if el is None:
            all_connected_inds.append([n])
        elif el not in written:
            all_connected_inds.append(list(el))
            written.append(el)
    return all_connected_inds
    

def write_unique_commands(stream, sc_tr_commands_list):
    # here we suppress repeated searches over the same positive constraint (or
    #   a more constrained search.
    to_write = [] 
    for n, stc_el in enumerate(sc_tr_commands_list):
        needed = True
        splits = stc_el.constr_splits
        score = stc_el.score
        if splits is not None:
            for other_ind, other in enumerate(sc_tr_commands_list):
                if other_ind == n:
                    continue
                _LOG.debug("n=%d, other_ind=%d, stc_el(splits = %s, score=%f), other(constr_splits = %s, score=%f)" % (n, other_ind, splits, score, other.constr_splits, other.score))
                if other_ind > n or ((other_ind < n) and (to_write[other_ind] is not None)):
                    osplits = other.constr_splits
                    #############################################################
                    # the commented out conditional reduces the number of searches, but may not generate enough diverse sets of trees
                    # if (osplits is None or len(osplits.difference(splits)) == 0) and (splits != osplits or score < other.score):
                    #############################################################
                    if osplits == splits and score < other.score:
                        needed = False
                        if score > other.score:
                            other.score = score
                            other.commands[0] = stc_el.commands[0]
                        break
        if needed:
            to_write.append(stc_el.commands)
        else:
            to_write.append(None)

    # finally we write out everything that remains in the list
    first = True
    for cmds in to_write:
        if cmds is not None:
            if first:
                cmds.insert(1, "treenum = 1\n")
                first = False
            stream.write("".join(cmds))

class ScoreConstraintCommands(object):
    def __init__(self, score, constr_splits, commands):
        self.score, self.constr_splits, self.commands = score, constr_splits, commands

def gather_neighborhood_commands(tree_list):
    dim = len(tree_list)
    mat = [[0]*dim for i in range(dim)]
    _LOG.debug("tree_list = %s\n" % str(tree_list))
    for row_n, i in enumerate(tree_list[:-1]):
        offset = 1 + row_n
        row_splits = tree_list[row_n].splits
        _LOG.debug("row_splits = %s\n" % str(row_splits))
        row = mat[row_n]
        for col_disp, j in enumerate(tree_list[offset:]):
            col_n = offset + col_disp
            col_splits = tree_list[col_n].splits
            _LOG.debug("col_splits = %s\n" % str(col_splits))
            d = len(row_splits.symmetric_difference(col_splits))
            _LOG.debug("d = %s\n" % str(d))
            row[col_n] = d
            mat[col_n][row_n] = d
    connected_indices = connected_at(mat, 4)
    _LOG.debug("connected_indices = %s\n" % str(connected_indices))
    sc_tr_commands_list = []
    for group_indices in connected_indices:
        trees = [tree_list[i] for i in group_indices]
        first_tree = trees[0]
        cmd_list = ["model = %s\ntree = %s\n" % (first_tree.model, first_tree.tree_string)]
        sc_tr_commands = ScoreConstraintCommands(first_tree.score, None, cmd_list)
        cmd_list.append("clearconstraints = 1\n")

        
        if len(trees) > 1:
            trees.sort(reverse=True, cmp=lambda x,y: cmp(x.score,y.score))
            c = copy.deepcopy(first_tree.tree)
            encode_splits(c)
            e = find_edge_from_split(c.seed_node, last_split)
            edge_dist = 3
            e.head_node.collapse_neighborhood(edge_dist)
            c.splits = get_norm_nontrivial_split_set(c)
            si = set(c.splits)
            for t in trees[1:]:
                si.intersection_update(t.splits)
            if len(si):
                sd = SplitDistribution(taxa_block=taxa_block, split_set=si)
                ts = TreeSummarizer()
                sc = ts.tree_from_splits(sd, min_freq=None, include_edge_lengths=False)
                encode_splits(sc)
                sc.splits = si
                sc_tr_commands.constr_splits = si
                cmd_list.append("posconstraint = %s\n" % sc.compose_newick(edge_lengths=False))
        else:
            c = copy.deepcopy(first_tree.tree)
            e = find_edge_from_split(c.seed_node, last_split)
            edge_dist = 3
            e.head_node.collapse_neighborhood(edge_dist)
            encode_splits(c)
            sc_tr_commands.constr_splits = get_norm_nontrivial_split_set(c)
            if len(sc_tr_commands.constr_splits) > 0:
                cmd_list.append("posconstraint = %s\n" % c.compose_newick(edge_lengths=False))
        cmd_list.append("run\n")
        sc_tr_commands_list.append(sc_tr_commands)
    return sc_tr_commands_list

def find_best_conflicting(self, starting_tree, split, dataset):
    new_starting = TreeModel(model=starting_tree.model)
    new_starting.tree = copy.deepcopy(starting_tree.tree)
    root = new_starting.tree.seed_node
    e = find_edge_from_split(root, split, root.edge.clade_mask)
    if e:
        e.collapse()

    tmp_tree_filename = ".tmp.tre"
    write_trees_to_filepath([new_starting], dataset, tmp_tree_filename)

    try:

        tmp_constrain_filename = ".tmpconstrain.tre"
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


n_tax = int(sys.argv[1])
last_split = 1 << (n_tax - 1)
all_taxa_bitmask = (1 << n_tax) - 1

add_trees_fn = sys.argv[2]
if len(sys.argv) > 3:
    nbhd_tree_groups = []
    for nbhd_tree_fn in sys.argv[3:]:
        nbhd_tree_f = open(nbhd_tree_fn, 'rU')
        nbhd_tree_groups.extend(read_add_tree_groups(nbhd_tree_f))
else:
    nbhd_tree_groups = None
    
add_trees_f = open(add_trees_fn, 'rU')
all_tree_groups = read_add_tree_groups(add_trees_f)

taxa_block = TaxaBlock([str(i+1) for i in range(n_tax)])
taxa_blocks = [taxa_block]
dataset = Dataset(taxa_blocks=taxa_blocks)

#setting this > 1.0 means that more trees are retained to the neighborhood search stage
score_diff_multiplier = 1.0
    
commands = []
if nbhd_tree_groups is None:
    _LOG.debug("Invocation of igarli_neighborhood.py with only one tree file -- need to set up initial neighborhood searches") 
    for g in all_tree_groups:
        for el in g:
            newick_string = el.tree_string
            newick_stream = StringIO(newick_string)
            t = dataset.read_trees(newick_stream, format="newick")[0]
            encode_splits(t)
            el.tree = t
        _LOG.debug("len(g) = %d" % len(g))
        opt_tree_el = g[0]
        opt_tree = opt_tree_el.tree
        opt_tree_el.splits = get_norm_nontrivial_split_set(opt_tree)
        _LOG.debug("opt_tree_el.splits = %s" % str(opt_tree_el.splits))
        unopt_score = None
        to_preserve = [opt_tree_el]
        for el in g[1:]:
            other_tree = el.tree
            el.splits = get_norm_nontrivial_split_set(other_tree)
            if unopt_score is None and el.splits == opt_tree_el.splits:
                unopt_score = el.score
            else:
                to_preserve.append(el)
    
        if unopt_score is not None:
            to_consider = to_preserve[1:]
            to_preserve = to_preserve[:1]
            min_score = unopt_score - score_diff_multiplier*(opt_tree_el.score - unopt_score)
            for el in to_consider:
                if el.score > min_score:
                    to_preserve.append(el)
        stc = gather_neighborhood_commands(to_preserve)
        commands.extend(stc)
    write_unique_commands(sys.stdout, commands)
else:
    # first we collect all of the ParsedTree objects into all_parsed_trees and we
    #   call encode_splits so that we can look up split info on each tree
    all_tree_groups.extend(nbhd_tree_groups)
    all_parsed_trees = []
    for g in all_tree_groups:
        for el in g:
            newick_string = el.tree_string
            newick_stream = StringIO(newick_string)
            t = dataset.read_trees(newick_stream, format="newick")[0]
            encode_splits(t)
            el.tree = t
            all_parsed_trees.append(el)
    assert(all_taxa_bitmask == all_parsed_trees[0].tree.seed_node.edge.clade_mask)

    ########################################
    # First, we make sure that there are not duplicate topologies
    # Because we reverse sort, we'll be retaining the tree with the
    #   best score
    #####
    all_parsed_trees.sort(reverse=True)
    set_of_split_sets = set()
    unique_topos = []
    for tm in all_parsed_trees:
        add_nontriv_splits_attr(tm.tree, all_taxa_bitmask)
        tm.splits, tm.split_set = tm.tree.splits, tm.tree.split_set
        if tm.tree.splits not in set_of_split_sets:
            unique_topos.append(tm)
            set_of_split_sets.add(tm.tree.splits)
    curr_results = unique_topos
    set_of_split_sets.clear()
    _LOG.info('There were %d unique result topologies for ntax = %d ' % (len(curr_results), n_tax))

    ########################################
    # the trees can be hefty, so lets eliminate unneeded references
    #####
    del unique_topos


    ########################################
    # Make sure to keep the current ML estimate in the next_round_trees list
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
    for split in ml_est.splits:
        found = False
        for n, tm in enumerate(curr_results):
            if split not in tm.split_set:
                best_disagreeing_index_set.add(n)
                found = True
                break
        if not found:
            unanimous_splits.append(split)

    if not unanimous_splits:
        # no more neighborhood searches are needed - exit without writing anything
        _LOG.debug("Every tree in the ML estimate is contradicted by at least one tree")
        sys.exit(0)
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
            add_nontriv_splits_attr(b.tree, all_taxa_bitmask)
            b.splits, b.split_set = tm.tree.splits, tm.tree.split_set

        best_conflicting.sort(reverse=True)
        tm = best_conflicting[0]
        next_round_trees.append(tm)
        curr_results.extend(best_conflicting[1:])
