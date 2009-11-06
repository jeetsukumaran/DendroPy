#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
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
###############################################################################

"""
Functions and classes in support of comparing data.
"""

def compare_datasets(ds1, ds2, tester, distinct_taxa=True, distinct_oids=False):
    tester.logger.info("Comparing dataset taxon sets ...")
    compare_dataset_taxon_sets(ds1, ds2, tester, distinct_taxa, distinct_oids)
    tester.logger.info("Comparing dataset tree lists ...")
    compare_dataset_tree_lists(ds1, ds2, tester, distinct_taxa, distinct_oids)
    tester.logger.info("Comparing dataset character arrays ...")
    compare_dataset_char_arrays(ds1, ds2, tester, distinct_taxa, distinct_oids)

def compare_dataset_taxon_sets(ds1, ds2, tester, distinct_taxa=True, distinct_oids=False):
    tester.assertEqual(len(ds1.taxon_sets), len(ds2.taxon_sets))
    if distinct_taxa:
        tester.assertTrue(ds1.taxon_sets is not ds2.taxon_sets)
    for ts_idx, ts1 in enumerate(ds1.taxon_sets):
        ts2 = ds2.taxon_sets[ts_idx]
        tester.logger.info("Comparing taxa of taxon set %d: %d taxa vs. %d taxa" \
            % (ts_idx, len(ts1), len(ts2)))
        compare_individual_taxon_sets(ts1, ts2, tester, distinct_taxa, distinct_oids)

def compare_individual_taxon_sets(ts1, ts2, tester, distinct_taxa=True, distinct_oids=False):
    if distinct_taxa:
        tester.assertTrue(ts1 is not ts2)
    tester.assertEqual(len(ts1), len(ts2))
    if distinct_oids:
        tester.assertNotEqual(ts1.oid, ts2.oid)
    tester.assertEqual(ts1.label, ts2.label)
    for taxon_idx, taxon1 in enumerate(ts1):
        tester.logger.debug("Taxon %d: '%s' == '%s'" % (taxon_idx, taxon1.label, ts2[taxon_idx].label))
        taxon2 = ts2[taxon_idx]
        if distinct_taxa:
            tester.assertTrue(taxon1 is not taxon2)
        tester.assertEqual(taxon1.label, taxon2.label)
        if distinct_oids:
            tester.assertNotEqual(taxon1.oid, taxon2.oid)

def compare_dataset_tree_lists(ds1, ds2, tester, distinct_taxa=True, distinct_oids=False):
    tester.assertTrue(ds1.tree_lists is not ds2.tree_lists)
    tester.assertEqual(len(ds1.tree_lists), len(ds2.tree_lists))
    for tree_list_idx, tree_list1 in enumerate(ds1.tree_lists):
        tree_list2 = ds2.tree_lists[tree_list_idx]
        if distinct_taxa:
            tester.assertTrue(tree_list1.taxon_set is not tree_list2.taxon_set)
            tester.assertTrue(tree_list1.taxon_set in ds1.taxon_sets)
            tester.assertTrue(tree_list2.taxon_set in ds2.taxon_sets)
            tester.assertTrue(tree_list1.taxon_set not in ds2.taxon_sets)
            tester.assertTrue(tree_list2.taxon_set not in ds1.taxon_sets)
        compare_individual_tree_lists(tree_list1, tree_list2, tester, distinct_taxa, distinct_oids)

def compare_individual_tree_lists(tree_list1, tree_list2, tester, distinct_taxa=True, distinct_oids=False):
    tester.assertTrue(tree_list1 is not tree_list2)
    tester.assertEqual(len(tree_list1), len(tree_list2))
    if distinct_taxa:
        tester.assertTrue(tree_list1.taxon_set is not tree_list2.taxon_set)
    compare_individual_taxon_sets(tree_list1.taxon_set, tree_list2.taxon_set, tester, distinct_taxa, distinct_oids)
    tester.assertEqual(tree_list1.label, tree_list2.label)
    if distinct_oids:
        tester.assertNotEqual(tree_list1.oid, tree_list2.oid)
    for tree_idx, tree1 in enumerate(tree_list1):
        tree2 = tree_list2[tree_idx]
        compare_individual_trees(tree1, tree2, tester, distinct_taxa, distinct_oids)

def compare_individual_trees(tree1, tree2, tester, distinct_taxa=True, distinct_oids=False):
        tester.logger.debug(tree1.to_newick_str())
        tree1.debug_check_tree(logger=tester.logger)
        tester.logger.debug(tree2.to_newick_str())
        tree2.debug_check_tree(logger=tester.logger)

        tester.assertTrue(tree1 is not tree2)
        if distinct_taxa:
            tester.assertTrue(tree1.taxon_set is not tree2.taxon_set)
        tester.assertTrue(tree1.taxon_set is tree_list1.taxon_set)
        tester.assertTrue(tree2.taxon_set is tree_list2.taxon_set)

        tree1_nodes = [nd for nd in tree1.postorder_node_iter()]
        tree2_nodes = [nd for nd in tree2.postorder_node_iter()]
        tester.assertEqual(len(tree1_nodes), len(tree2_nodes))
        for nd_idx, node1 in enumerate(tree1_nodes):
            node2 = tree2_nodes[nd_idx]
            if node1.taxon is not None:
                tester.assert_(node2.taxon is not None)
                if distinct_taxa:
                    tester.assertTrue(node1.taxon is not node2.taxon)
                else:
                    tester.assertTrue(node1.taxon is node2.taxon)
                if distinct_oids:
                    tester.assertNotEqual(node1.oid, node2.oid)
                tester.assertEqual(node1.taxon.label, node2.taxon.label)
                tester.assertTrue(node1.taxon in tree1.taxon_set)
                tester.assertTrue(node2.taxon in tree2.taxon_set)
            else:
                tester.assert_(node2.taxon is None)
            if node1.edge.length is not None:
                tester.assert_(node2.edge.length is not None)
                tester.assertAlmostEqual(node1.edge.length, node2.edge.length, 3)
            else:
                tester.assert_(node2.edge.length is None)
            tester.assertEqual(len(node1.child_nodes()), len(node2.child_nodes()))

def compare_dataset_char_arrays(ds1, ds2, tester, distinct_taxa=True, distinct_oids=False):
    tester.assertEqual(len(ds1.char_arrays), len(ds2.char_arrays))
    for char_array_idx, char_array1 in enumerate(ds1.char_arrays):
        char_array2 = ds2.char_arrays[char_array_idx]
        compare_individual_char_arrays(char_array1, char_array2, tester, distinct_taxa, distinct_oids)

def compare_individual_char_arrays(char_array1, char_array2, tester, distinct_taxa=True, distinct_oids=False):
    tester.assertEqual(len(char_array1), len(char_array2))
    tester.assertEqual(len(char_array1.taxon_set), len(char_array2.taxon_set))
    for taxon_idx, taxon1 in enumerate(char_array1.taxon_set):
        tester.assertEqual(char_array1.taxon_set[taxon_idx].label,
                char_array2.taxon_set[taxon_idx].label)
        seq1 = char_array1[taxon_idx]
        seq2 = char_array2[taxon_idx]
        tester.assertEqual(len(seq1), len(seq2))
        for cell_idx, cell1 in enumerate(seq1):
            cell2 = seq2[cell_idx]
            state1 = cell1.value
            state2 = cell2.value
            tester.assertEqual(state1.symbol, state2.symbol)
            tester.assertEqual(state1.token, state2.token)
            tester.assertEqual(state1.multistate, state2.multistate)
            tester.assertEqual(state1.fundamental_symbols, state2.fundamental_symbols)
