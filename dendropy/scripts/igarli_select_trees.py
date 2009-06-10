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
_LOG = get_logger("igarli_select.py")

# this should not be hard coded
MAX_TREES_CARRIED_OVER = 1000
SPLIT_DIVERSITY_MULTIPLIER = 1.0


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
    def __str__(self):
        if self.model:
            if self.score:
                return "[!GarliScore %f][!GarliModel  %s ] %s" % (self.score, self.model, self.tree_string)
            return "[!GarliModel  %s ] %s" % (self.model, self.tree_string)
        return "[!GarliModel  %s ] %s" % (self.tree_string)
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



n_tax = int(sys.argv[1])
last_split = 1 << (n_tax - 1)
all_taxa_bitmask = (1 << n_tax) - 1

add_trees_fn = sys.argv[2]
assert len(sys.argv) > 3
nbhd_tree_groups = []
for nbhd_tree_fn in sys.argv[3:]:
    nbhd_tree_f = open(nbhd_tree_fn, 'rU')
    nbhd_tree_groups.extend(read_add_tree_groups(nbhd_tree_f))
    
add_trees_f = open(add_trees_fn, 'rU')
all_tree_groups = read_add_tree_groups(add_trees_f)

taxa_block = TaxaBlock([str(i+1) for i in range(n_tax)])
taxa_blocks = [taxa_block]
dataset = Dataset(taxa_blocks=taxa_blocks)

#setting this > 1.0 means that more trees are retained to the neighborhood search stage
score_diff_multiplier = 1.0
    
commands = []
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
# We should have already run the igarli_neighborhood script enough so that 
#   there are not unanimous_splits
#####
assert len(unanimous_splits) == 0


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
max_len = MAX_TREES_CARRIED_OVER
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
    def split_diversity_cmp(x, y, lambda_mult=SPLIT_DIVERSITY_MULTIPLIER):
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
_LOG.info('Added %d trees that were not "required" to guarantee that no splits were unanimous for ntax = %d' % (n_added, n_tax))

for n, nt in enumerate(next_round_trees):
    sys.stdout.write("Tree tree%d = [&U]%s ;\n" % (n, str(nt)))
_LOG.info('A total of %d trees were retained for ntax = %d lnL range from %f to %f' % (len(next_round_trees), n_tax, next_round_trees[0].score, next_round_trees[-1].score))

sys.exit(0)


"""     sys.stdout.write('tree = %s\n' % i)
        if first: # this is a way to make 
            sys.stdout.write('treeNum = 1\n')
            first = False
        sys.stdout.write('run\n')
    elif False:
        _LOG.debug("nomatch: %s\n" %line)
sys.exit(0)
"""
