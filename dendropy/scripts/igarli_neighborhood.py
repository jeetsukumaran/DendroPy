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

# bail out if you are given multiple files, this means that you are in the second round
#   of neighborhood searching
if len(sys.argv) > 3: 
    sys.exit(0)

#garli_dna_model_pattern = r'\[!GarliModel  r ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) ([.0-9]*) e ([.0-9]*) ([.0-9]*) ([.0-9]*) [.0-9]* a ([.0-9]*) p ([.0-9]*) \]'
garli_dna_model_pattern = r'\[!GarliModel\s*(r [.0-9]* [.0-9]* [.0-9]* [.0-9]* [.0-9]* e [.0-9]* [.0-9]* [.0-9]* [.0-9]* a [.0-9]* p [.0-9]*)\s*\]'
garli_after_name_tm_pattern = r'\s*=\s*\[&U\]\[!GarliScore ([-.0-9]*)\]\s*%s\s*(\(.*\))\s*;' % garli_dna_model_pattern
garli_after_name_tree_pattern = r'\s*=\s*\[&U\]\[!GarliScore ([-.0-9]*)\]\s*(\(.*\))\s*;'



tree_prefix = r'\[iGarli\s*(\d+)\s*\] tree best'
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
    
def read_add_tree_groups(f):
    first = True
    curr_tree_group = []
    all_tree_groups = []
    # this is hacky -- we use the presence of a model string to tell us 
    #   when we have hit a new group of trees that were found by adding to 
    #   the same incomplete tree (the other trees will come from the suboptimal
    #   trees that are spit out without any model info because the hook
    #   in the Finish...Stepwise...() function does not expose the model)
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
    
    
    if curr_model:
        for el in curr_tree_group:
            if el.model is None:
                el.model = curr_model
        all_tree_groups.append(curr_tree_group)
    return all_tree_groups

def get_norm_nontrivial_split_set(tree):
    norm_non_triv = set()
    for split in tree.split_edges.keys():
        if not is_trivial_split(split, mask):
            if split & 1:
                norm_non_triv.add(split)
            else:
                comp_split = (~split) & mask
                norm_non_triv.add(comp_split)
    return norm_non_triv

def connected_at(dist_mat, max_dist):
    sys.stderr.write("dist_mat =\n%s\n" % "\n".join([str(row) for row in dist_mat]))

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
    

first_write = True

def write_neighborhood_commands(stream, tree_list):
    global first_write
    dim = len(tree_list)
    mat = [[0]*dim for i in range(dim)]
    sys.stderr.write("tree_list = %s\n" % str(tree_list))
    for row_n, i in enumerate(tree_list[:-1]):
        offset = 1 + row_n
        row_splits = tree_list[row_n].splits
        sys.stderr.write("row_splits = %s\n" % str(row_splits))
        row = mat[row_n]
        for col_disp, j in enumerate(tree_list[offset:]):
            col_n = offset + col_disp
            col_splits = tree_list[col_n].splits
            sys.stderr.write("col_splits = %s\n" % str(col_splits))
            d = len(row_splits.symmetric_difference(col_splits))
            sys.stderr.write("d = %s\n" % str(d))
            row[col_n] = d
            mat[col_n][row_n] = d
    connected_indices = connected_at(mat, 4)
    sys.stderr.write("connected_indices = %s\n" % str(connected_indices))
    for group_indices in connected_indices:
        trees = [tree_list[i] for i in group_indices]
        first_tree = trees[0]
        stream.write("model = %s\n" % first_tree.model)
        stream.write("tree = %s\n" % first_tree.tree_string)
        if first_write:
            stream.write("treenum = 1\n")
            first_write = False
        stream.write("clearconstraints = 1\n")

        
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
                stream.write("posconstraint = %s\n" % sc.compose_newick(edge_lengths=False))
        else:
            c = copy.deepcopy(first_tree.tree)
            e = find_edge_from_split(c.seed_node, last_split)
            edge_dist = 3
            e.head_node.collapse_neighborhood(edge_dist)
            stream.write("posconstraint = %s\n" % c.compose_newick(edge_lengths=False))
        stream.write("run\n")
    
n_tax = int(sys.argv[1])
last_split = 1 << (n_tax - 1)
mask = (1 << n_tax) - 1

add_trees_fn = sys.argv[2]
add_trees_f = open(add_trees_fn, 'rU')
all_tree_groups = read_add_tree_groups(add_trees_f)

taxa_block = TaxaBlock([str(i+1) for i in range(n_tax)])
taxa_blocks = [taxa_block]
dataset = Dataset(taxa_blocks=taxa_blocks)

#setting this > 1.0 means that more trees are retained to the neighborhood search stage
score_diff_multiplier = 1.0
    

for g in all_tree_groups:
    for el in g:
        newick_string = el.tree_string
        newick_stream = StringIO(newick_string)
        t = dataset.read_trees(newick_stream, format="newick")[0]
        encode_splits(t)
        el.tree = t
    opt_tree_el = g[0]
    opt_tree = opt_tree_el.tree
    opt_tree_el.splits = get_norm_nontrivial_split_set(opt_tree)
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
    write_neighborhood_commands(sys.stdout, to_preserve)

sys.exit(0)


"""     sys.stdout.write('tree = %s\n' % i)
        if first: # this is a way to make 
            sys.stdout.write('treeNum = 1\n')
            first = False
        sys.stdout.write('run\n')
    elif False:
        sys.stderr.write("nomatch: %s\n" %line)
sys.exit(0)
"""
