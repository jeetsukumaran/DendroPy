#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Generate trees with incomplete leaf sets.
"""

import dendropy
import random
from dendropy.interop import paup
from dendropy.test.support import pathmap

def generate_pruned_trees(
        src_trees_fname,
        num_reps,
        num_trees_per_rep):
    rng = random.Random()
    trees = dendropy.TreeList.get_from_path(
            src=pathmap.tree_source_path(src_trees_fname),
            schema='nexus')
    taxa = trees.taxon_set
    # print "1 >>>>", id(taxa), ":", len(taxa)
    # for t in taxa:
    #     print repr(t)
    # input_trees = open(output_prepruned_tree_file_path, "w")
    # output_trees = open(output_postpruned_tree_file_path, "w")
    input_dataset = dendropy.DataSet(attached_taxon_set=taxa)
    output_dataset = dendropy.DataSet(attached_taxon_set=taxa)
    pruned_taxa = []
    retained_taxa = []
    for rep in range(num_reps):
        sub_trees = [dendropy.Tree(t, taxon_set=taxa) for t in rng.sample(trees, num_trees_per_rep)]
        sub_trees = dendropy.TreeList(sub_trees, taxon_set=taxa)
        sub_size = rng.randint(5, len(taxa)-5)
        assert sub_size > 0
        assert sub_size < len(taxa)
        sub_taxa = rng.sample(taxa, sub_size)
        assert len(sub_taxa) > 4
        assert len(sub_taxa) < len(taxa)
        # if retain_taxa_in_list:
        #     taxa_to_prune = [t for t in taxa if t not in sub_taxa]
        #     taxa_to_retain = sub_taxa
        # else:
        #     taxa_to_prune = sub_taxa
        #     taxa_to_retain = [t for t in taxa if t not in sub_taxa]
        taxa_to_prune = sub_taxa
        taxa_to_retain = [t for t in taxa if t not in sub_taxa]
        pruned_trees = paup.prune_taxa_from_trees(sub_trees, taxa_to_prune)
        pruned_taxa.append(taxa_to_prune)
        retained_taxa.append(taxa_to_retain)
        assert sub_trees.taxon_set is taxa
        input_dataset.add_tree_list(sub_trees)
        assert pruned_trees.taxon_set is taxa
        output_dataset.add_tree_list(pruned_trees)
    # print "2 >>>>", id(taxa), ":", len(taxa)
    # for t in taxa:
    #     print repr(t)
    for trees in input_dataset.tree_lists:
        assert trees.taxon_set is taxa
        for tree in trees:
            assert tree.taxon_set is taxa
            count = 0
            for nd in tree.postorder_node_iter():
                if nd.taxon is not None:
                    count += 1
            assert count == len(taxa)
    for trees in output_dataset.tree_lists:
        assert trees.taxon_set is taxa
        for tree in trees:
            assert tree.taxon_set is taxa
    return taxa, pruned_taxa, retained_taxa, input_dataset, output_dataset

if __name__ == "__main__":
    import sys
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("src_tree_filename",
            type=str,
            help="name of source of trees (should be located under 'data/trees' of test directory)")
    parser.add_argument("num_reps",
            type=int,
            help="number of replicates")
    parser.add_argument("num_trees_per_rep",
            type=int,
            help="number of trees subsampled from source of trees per replicate")
    parser.add_argument("output_prefix",
            type=str,
            help="prefix of output files")
    args = parser.parse_args()
    all_taxa, pruned_taxa, retained_taxa, input_dataset, output_dataset = generate_pruned_trees(args.src_tree_filename,
            args.num_reps,
            args.num_trees_per_rep)

    # print "3 >>>>", id(all_taxa), ":", len(all_taxa)
    # for t in all_taxa:
    #     print repr(t)
    input_dataset.write_to_path(("%s.pre-pruned.nex" % args.output_prefix), "nexus")
    output_dataset.write_to_path(("%s.paup-pruned.nex" % args.output_prefix), "nexus")
    f = open(("%s.retained_taxa.txt" % args.output_prefix), "w")
    for rt in retained_taxa:
        f.write("%s\n" % (" ".join(str(all_taxa.index(t)) for t in rt)))
    f = open(("%s.pruned_taxa.txt" % args.output_prefix), "w")
    for pt in pruned_taxa:
        f.write("%s\n" % (" ".join(str(all_taxa.index(t)) for t in pt)))
    for set_idx, pt in enumerate(pruned_taxa):
        print "set {}".format(set_idx+1)
        src_trees = input_dataset.tree_lists[set_idx]
        ref_trees = output_dataset.tree_lists[set_idx]
        assert src_trees.taxon_set is all_taxa
        assert ref_trees.taxon_set is all_taxa
        for tree_idx, src_tree in enumerate(src_trees):
            ref_tree = ref_trees[tree_idx]
            assert src_tree.taxon_set is all_taxa
            assert ref_tree.taxon_set is all_taxa
            count = 0
            for nd in src_tree.postorder_node_iter():
                if nd.taxon is not None:
                    count += 1
            src_tree.prune_taxa(pt)
            tree_dist = paup.symmetric_difference(src_tree, ref_tree)
            if tree_dist != 0:
                print "\n##################################################\n"
                print "Post-pruning conflict: Symmetric Difference = {}\n".format(tree_dist)
                print "Pruned taxa ({}/{}):".format(len(pt), len(all_taxa))
                for t in pt:
                    print "   {}".format(t.label)
                print
                print src_tree.as_string("newick")
                print ref_tree.as_string("newick")
                print
                print len(src_tree.leaf_nodes())
                seen_taxa = []
                for nd in src_tree.leaf_nodes():
                    if nd.taxon is None:
                        print "null leaf node!"
                    else:
                        try:
                            taxon_idx = all_taxa.index(nd.taxon)
                        except ValueError:
                            print "4 >>>>", id(all_taxa), ":", len(all_taxa)
                            for t in all_taxa:
                                print repr(t)
                            raise
                        if taxon_idx in pt:
                            print "not removed: {}".format(nd.taxon.label)
                        else:
                            seen_taxa.append(nd.taxon)
                for taxon in retained_taxa[set_idx]:
                    if taxon not in seen_taxa:
                        print "wrongly removed: {}".format(taxon.label)
                print len(ref_tree.leaf_nodes())
                print "##################################################\n"
                sys.exit(1)
    input_dataset.write_to_path(("%s.dendropy-pruned.nex" % args.output_prefix), "nexus")
