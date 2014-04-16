#! /usr/bin/env python

import sys
import json
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
from dendropy.test.support import pathmap
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

tree_file_titles = [
    'standard-test-trees-n14-unrooted-treeshapes',
    'standard-test-trees-n10-rooted-treeshapes',
    'standard-test-trees-n12-x2',
    'standard-test-trees-n33-x100a',
    'standard-test-trees-n33-x10a',
    'standard-test-trees-n33-x10b',
    'standard-test-trees-n33-annotated',
]

schema_extension_map = {
    "newick" : "newick",
    "nexus" : "nexus",
    "json" : "json",
}

tree_filepaths = {}
for schema in schema_extension_map:
    tree_filepaths[schema] = {}
    for tree_file_title in tree_file_titles:
        tf = "{}.{}".format(tree_file_title, schema_extension_map[schema])
        tree_filepaths[schema][tree_file_title] = pathmap.tree_source_path(tf)
tree_references = {}
for tree_file_title in tree_file_titles:
    with open(tree_filepaths["json"][tree_file_title]) as src:
        tree_references[tree_file_title] = json.load(src)

class StandardTestTreeChecker(object):

    def compare_comments(self,
            item,
            check,
            metadata_extracted=False):
        check_comments = list(check["comments"])
        item_comments = list(item.comments)
        for comment in item.comments:
            try:
                check_comments.remove(comment)
            except ValueError:
                pass
            else:
                item_comments.remove(comment)
        self.assertEqual(check_comments, [])
        if metadata_extracted:
            self.assertEqual(item_comments, [])
        else:
            for idx, c in enumerate(item_comments):
                if c.startswith("&"):
                    item_comments[idx] = c[1:]
            item_metadata_comments = ",".join(item_comments)
            check_metadata_comments = ",".join(check["metadata_comments"])
            self.maxDiff = None
            self.assertEqual(item_metadata_comments, check_metadata_comments)

    def label_nodes(self, tree):
        for node_idx, node in enumerate(tree):
            if node.taxon is not None:
                node.canonical_label = node.taxon.label
            else:
                node.canonical_label = node.label

    def compare_to_check_tree(self,
            tree,
            tree_file_title,
            check_tree_idx,
            suppress_internal_node_taxa=True,
            suppress_external_node_taxa=False,
            metadata_extracted=False,
            distinct_nodes_and_edges=True):
        check_tree = tree_references[tree_file_title][str(check_tree_idx)]
        self.assertIs(tree.is_rooted, check_tree["is_rooted"])
        self.compare_comments(tree, check_tree, metadata_extracted)
        seen_taxa = []
        node_labels = []
        edge_labels = []
        num_visited_nodes = 0
        self.label_nodes(tree)
        for node_idx, node in enumerate(tree):
            num_visited_nodes += 1
            check_node = check_tree["nodes"][node.canonical_label]
            check_node_label = check_node["label"]
            self.assertEqual(node.canonical_label, check_node_label)
            # node_labels.append(node.canonical_label)
            _LOG.debug("{}: {}: {}".format(tree_file_title, check_tree_idx, node.canonical_label))

            check_node_children = check_node["children"]
            if check_node_children:
                self.assertTrue(node.is_internal())
                self.assertFalse(node.is_leaf())
                self.assertEqual(len(node._child_nodes), len(check_node_children))
                if suppress_internal_node_taxa:
                    self.assertEqual(node.label, check_node_label)
                    self.assertIs(node.taxon, None)
                    node_labels.append(node.label)
                else:
                    self.assertIsNot(node.taxon, None)
                    self.assertEqual(node.taxon.label, check_node_label)
                    self.assertIs(node.label, None)
                    seen_taxa.append(node.taxon)
            else:
                self.assertFalse(node.is_internal())
                self.assertTrue(node.is_leaf())
                self.assertEqual(len(node._child_nodes), len(check_node_children))
                if suppress_external_node_taxa:
                    self.assertEqual(node.label, check_node_label)
                    self.assertIs(node.taxon, None)
                    node_labels.append(node.label)
                else:
                    self.assertIsNot(node.taxon, None)
                    self.assertEqual(node.taxon.label, check_node_label)
                    self.assertIs(node.label, None)
                    seen_taxa.append(node.taxon)

            if node.parent_node is not None:
                if node.parent_node.is_internal:
                    if suppress_internal_node_taxa:
                        self.assertEqual(node.parent_node.label, check_node["parent"])
                        self.assertIs(node.parent_node.taxon, None)
                    else:
                        self.assertEqual(node.parent_node.taxon.label, check_node["parent"])
                        self.assertIs(node.parent_node.label, None)
                else:
                    if suppress_external_node_taxa:
                        self.assertEqual(node.parent_node.label, check_node["parent"])
                        self.assertIs(node.parent_node.taxon, None)
                    else:
                        self.assertEqual(node.parent_node.taxon.label, check_node["parent"])
                        self.assertIs(node.parent_node.label, None)
            else:
                self.assertEqual(check_node["parent"], "None")

            child_labels = []
            for ch in node.child_node_iter():
                if ch.is_internal():
                    if suppress_internal_node_taxa:
                        self.assertIs(ch.taxon, None)
                        child_labels.append(ch.label)
                    else:
                        self.assertIsNot(ch.taxon, None)
                        child_labels.append(ch.taxon.label)
                        self.assertIs(ch.label, None)
                else:
                    if suppress_external_node_taxa:
                        self.assertIs(ch.taxon, None)
                        child_labels.append(ch.label)
                    else:
                        self.assertIsNot(ch.taxon, None)
                        child_labels.append(ch.taxon.label)
                        self.assertIs(ch.label, None)
            self.assertEqual(len(child_labels), len(check_node["children"]))
            self.assertEqual(set(child_labels), set(check_node["children"]))

            edge = node.edge
            check_edge = check_tree["edges"][node.canonical_label]
            if edge.tail_node is None:
                self.assertEqual(check_edge["tail_node"], "None")
            else:
                self.assertEqual(edge.tail_node.canonical_label, check_edge["tail_node"])
            self.assertEqual(edge.head_node.canonical_label, check_edge["head_node"])
            self.assertAlmostEqual(node.edge.length, float(check_edge["length"]))

            # This hackery because NEWICK/NEXUS cannot distinguish between
            # node and edge comments, and everything gets lumped in as a
            # node comment
            if not distinct_nodes_and_edges:
                node.comments += edge.comments
                d = {
                        "comments": check_node["comments"] + check_edge["comments"],
                        "metadata_comments": check_node["metadata_comments"] + check_edge["metadata_comments"],
                        }
                self.compare_comments(node, d, metadata_extracted)
            else:
                self.compare_comments(node, check_node, metadata_extracted)
                self.compare_comments(edge, check_node, metadata_extracted)

        self.assertEqual(num_visited_nodes, len(check_tree["nodeset"]))
        self.assertEqual(len(seen_taxa), len(tree.taxon_namespace))
        self.assertEqual(set(seen_taxa), set(tree.taxon_namespace))
        node_labels.extend([t.label for t in tree.taxon_namespace])
        self.assertEqual(len(node_labels), len(check_tree["nodeset"]))
        self.assertEqual(set(node_labels), set(check_tree["nodeset"]))

    def verify_standard_trees(self,
            tree_list,
            tree_file_title,
            tree_offset=0,
            suppress_internal_node_taxa=True,
            suppress_external_node_taxa=False,
            metadata_extracted=False,
            distinct_nodes_and_edges=True):
        tree_reference = tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        if tree_offset < 0:
            if abs(tree_offset) > expected_number_of_trees:
                tree_offset = 0
            else:
                tree_offset = expected_number_of_trees + tree_offset
        self.assertEqual(len(tree_list), expected_number_of_trees-tree_offset)
        # for tree_idx, (tree, check_tree) in enumerate(zip(tree_list, tree_directory[tree_file_title])):
        for tree_idx, tree in enumerate(tree_list):
            _LOG.debug("{}: {}".format(tree_file_title, tree_idx))
            self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
            self.compare_to_check_tree(
                    tree=tree,
                    tree_file_title=tree_file_title,
                    check_tree_idx=tree_idx + tree_offset,
                    suppress_internal_node_taxa=suppress_internal_node_taxa,
                    suppress_external_node_taxa=suppress_external_node_taxa,
                    metadata_extracted=metadata_extracted,
                    distinct_nodes_and_edges=distinct_nodes_and_edges)
