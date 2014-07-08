#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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

import sys
import json
import random
import dendropy
import copy
from dendropy.test.support import pathmap
if sys.hexversion < 0x03040000:
    from dendropy.utility.filesys import pre_py34_open as open
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)

tree_file_titles = [
    "dendropy-test-trees-multifurcating-rooted-annotated",
    "dendropy-test-trees-multifurcating-rooted",
    "dendropy-test-trees-multifurcating-unrooted",
    "dendropy-test-trees-n10-rooted-treeshapes",
    "dendropy-test-trees-n12-x2",
    "dendropy-test-trees-n14-unrooted-treeshapes",
    "dendropy-test-trees-n33-unrooted-annotated-x10a",
    "dendropy-test-trees-n33-unrooted-x10a",
    "dendropy-test-trees-n33-unrooted-x10b",
    "dendropy-test-trees-n33-unrooted-x100a",
]

schema_extension_map = {
    "nexml" : "nexml",
    "newick" : "newick",
    "nexus" : "nexus",
    "json" : "json",
    "nexus-metadata-comments" : "nexus-metadata-comments.json",
}

_TREE_FILEPATHS = {}
_TREE_REFERENCES = {}
_NEXUS_METADATA_COMMENTS = {}

def setup_module():
    for schema in schema_extension_map:
        _TREE_FILEPATHS[schema] = {}
        for tree_file_title in tree_file_titles:
            tf = "{}.{}".format(tree_file_title, schema_extension_map[schema])
            _TREE_FILEPATHS[schema][tree_file_title] = pathmap.tree_source_path(tf)
    for tree_file_title in tree_file_titles:
        with open(_TREE_FILEPATHS["json"][tree_file_title]) as src:
            _TREE_REFERENCES[tree_file_title] = json.load(src)
        if "annotated" in tree_file_title:
            with open(_TREE_FILEPATHS["nexus-metadata-comments"][tree_file_title]) as src:
                _NEXUS_METADATA_COMMENTS[tree_file_title] = json.load(src)
setup_module()

class StandardTestTreesChecker(object):

    def preprocess_tree_to_be_checked(self, tree):
        pass

    def compare_annotations_to_json_metadata_dict(self,
            item,
            expected_metadata):
        item_annotations_as_dict = item.annotations.values_as_dict()
        if self.__class__.is_coerce_metadata_values_to_string:
            for k in expected_metadata:
                v = expected_metadata[k]
                if isinstance(v, list):
                    v = [str(i) for i in v]
                elif isinstance(v, tuple):
                    v = (str(i) for i in v)
                else:
                    v = str(v)
                expected_metadata[k] = v

        # for annote in item.annotations:
        #     print("{}: {}".format(annote.name, annote.value))
        # k1 = sorted(list(item_annotations_as_dict.keys()))
        # k2 = sorted(list(expected_metadata.keys()))
        # print("--")
        # for k in k1:
        #     print("'{}':'{}'".format(k, item_annotations_as_dict[k]))
        # print("--")
        # for k in k2:
        #     print("'{}':'{}'".format(k, expected_metadata[k]))
        # self.assertEqual(len(k1), len(k2))
        # self.assertEqual(set(k1), set(k2))
        # for key in set(item_annotations_as_dict.keys()):
        #     if item_annotations_as_dict[key] != expected_metadata[key]:
        #         v = expected_metadata[key]
        #         # if isinstance(v, list):
        #         #     print("{}: {}".format(v, [type(i) for i in v]))
        #         # elif isinstance(v, tuple):
        #         #     print("{}: {}".format(v, (type(i) for i in v)))
        #         # else:
        #         #     print("{}: {}".format(v, type(v)))
        #         print("**** {}:\t\t{} ({}) \t\t{} ({})".format(
        #             key,
        #             item_annotations_as_dict[key],
        #             type(item_annotations_as_dict[key]),
        #             expected_metadata[key],
        #             type(expected_metadata[key]),
        #             ))

        self.assertEqual(item_annotations_as_dict, expected_metadata)

    def check_metadata_annotations(self,
            item,
            reference):
        expected_annotations = reference["metadata"]
        self.compare_annotations_to_json_metadata_dict(item, expected_annotations)

    def check_comments(self, item, reference):
        reference_comments = list(reference["comments"])
        item_comments = list(item.comments)
        self.assertEqualUnorderedSequences(item_comments, reference_comments)

    def compare_to_reference_tree(self,
            tree,
            tree_file_title,
            reference_tree_idx):
        ref_tree = self.tree_references[tree_file_title][str(reference_tree_idx)]
        self.assertIs(tree.is_rooted, ref_tree["is_rooted"])
        if self.__class__.is_check_comments:
            self.check_comments(tree, ref_tree)
        self.check_metadata_annotations(
                item=tree,
                reference=ref_tree)
        obs_taxon_labels = []
        obs_node_labels = []
        obs_edge_labels = []
        visited_nodes = []
        self.preprocess_tree_to_be_checked(tree)
        for node_idx, node in enumerate(tree):
            visited_nodes.append(node)
            ref_node = ref_tree["nodes"][node.label]
            ref_node_label = ref_node["label"]
            ref_node_taxon_label = ref_node["taxon_label"]
            self.assertEqual(node.label, ref_node_label)
            ref_node_children = ref_node["children"]
            self.assertEqual(len(node._child_nodes), len(ref_node_children))
            if node.taxon:
                self.assertEqual(node.taxon.label, ref_node_taxon_label)
                obs_taxon_labels.append(node.taxon.label)
            else:
                self.assertEqual(ref_node_taxon_label, None)
            obs_node_labels.append(node.label)
            if ref_node_children:
                self.assertTrue(node.is_internal())
                self.assertFalse(node.is_leaf())
            else:
                self.assertFalse(node.is_internal())
                self.assertTrue(node.is_leaf())

            if node.parent_node is None:
                self.assertEqual(ref_node["parent"], None)
            else:
                self.assertEqual(node.parent_node.label, ref_node["parent"])
                if node.parent_node.taxon:
                    self.assertEqual(node.parent_node.taxon.label, ref_node["parent"])
            child_labels = [ch.label for ch in node.child_node_iter()]
            self.assertEqual(len(child_labels), len(ref_node["children"]))
            self.assertEqual(set(child_labels), set(ref_node["children"]))

            edge = node.edge
            ref_edge = ref_tree["edges"][edge.label]
            if edge.tail_node is None:
                self.assertEqual(ref_edge["tail_node"], None)
            else:
                self.assertEqual(edge.tail_node.label, ref_edge["tail_node"])
            self.assertEqual(edge.head_node.label, ref_edge["head_node"])
            self.assertAlmostEqual(node.edge.length, float(ref_edge["length"]), 3)

    def compare_to_reference_tree2(self,
            tree,
            tree_file_title,
            reference_tree_idx,
            suppress_internal_node_taxa,
            suppress_leaf_node_taxa,
            is_metadata_extracted,
            is_coerce_metadata_values_to_string,
            is_distinct_nodes_and_edges_representation,
            is_taxa_managed_separately_from_tree):
        ref_tree = self.tree_references[tree_file_title][str(reference_tree_idx)]
        self.assertIs(tree.is_rooted, ref_tree["is_rooted"])
        self.check_comments(
                tree,
                ref_tree)
        self.check_metadata_annotations(
                item=tree,
                reference=ref_tree,
                is_coerce_metadata_values_to_string=False)
        obs_taxa = []
        obs_node_labels = []
        obs_edge_labels = []
        visited_nodes = []
        for node in tree:
            if node.taxon is not None:
                node.canonical_label = node.taxon.label
            else:
                node.canonical_label = node.label
        for node_idx, node in enumerate(tree):
            visited_nodes.append(node)
            ref_node = ref_tree["nodes"][node.label]
            ref_node_label = ref_node["label"]
            self.assertEqual(node.label, ref_node_label)
            ref_edge = ref_tree["edges"][node.edge.label]
            self.assertEqual(node.edge.label, ref_edge["label"])
            ref_node_children = ref_node["children"]
            if ref_node_children:
                self.assertTrue(node.is_internal())
                self.assertFalse(node.is_leaf())
                self.assertEqual(len(node._child_nodes), len(ref_node_children))
                if suppress_internal_node_taxa:
                    self.assertEqual(node.label, ref_node_label)
                    self.assertIs(node.taxon, None)
                    obs_node_labels.append(node.label)
                else:
                    if not self.__class__.is_distinct_taxa_and_labels_on_tree:
                        self.assertIsNot(node.taxon, None)
                        self.assertEqual(node.taxon.label, ref_node_label)
                        self.assertIs(ch.label, None)
                        obs_taxa.append(node.taxon)
                    else:
                        obs_taxa.append(node.label)
                        obs_node_labels.append(node.label)
            else:
                self.assertFalse(node.is_internal())
                self.assertTrue(node.is_leaf())
                self.assertEqual(len(node._child_nodes), len(ref_node_children))
                if suppress_leaf_node_taxa:
                    self.assertEqual(node.label, ref_node_label)
                    self.assertIs(node.taxon, None)
                    obs_node_labels.append(node.label)
                else:
                    self.assertIsNot(node.taxon, None)
                    self.assertEqual(node.taxon.label, ref_node_label)
                    if not self.__class__.is_distinct_taxa_and_labels_on_tree:
                        self.assertIs(ch.label, None)
                    obs_taxa.append(node.taxon)
            if node.parent_node is not None:
                if node.parent_node.is_internal:
                    if suppress_internal_node_taxa:
                        self.assertEqual(node.parent_node.label, ref_node["parent"])
                        self.assertIs(node.parent_node.taxon, None)
                    else:
                        if not self.__class__.is_distinct_taxa_and_labels_on_tree:
                            self.assertEqual(node.parent_node.taxon.label, ref_node["parent"])
                            self.assertIs(node.parent_node.label, None)
                        else:
                            self.assertEqual(node.parent_node.label, ref_node["parent"])
                else:
                    if suppress_leaf_node_taxa:
                        self.assertEqual(node.parent_node.label, ref_node["parent"])
                        self.assertIs(node.parent_node.taxon, None)
                    else:
                        self.assertEqual(node.parent_node.taxon.label, ref_node["parent"])
                        self.assertIs(node.parent_node.label, None)
            else:
                self.assertEqual(ref_node["parent"], "None")

            child_labels = []
            for ch in node.child_node_iter():
                if ch.is_internal():
                    if suppress_internal_node_taxa:
                        self.assertIs(ch.taxon, None)
                        child_labels.append(ch.label)
                    else:
                        if not self.__class__.is_distinct_taxa_and_labels_on_tree:
                            self.assertIsNot(ch.taxon, None)
                            child_labels.append(ch.taxon.label)
                            self.assertIs(ch.label, None)
                        else:
                            self.assertEqual(node.label, ref_node_label)
                            child_labels.append(ch.label)
                else:
                    if suppress_leaf_node_taxa:
                        self.assertIs(ch.taxon, None)
                        child_labels.append(ch.label)
                    else:
                        self.assertIsNot(ch.taxon, None)
                        child_labels.append(ch.taxon.label)
                        if not self.__class__.is_distinct_taxa_and_labels_on_tree:
                            self.assertIs(ch.label, None)
            self.assertEqual(len(child_labels), len(ref_node["children"]))
            self.assertEqual(set(child_labels), set(ref_node["children"]))
            edge = node.edge
            ref_edge = ref_tree["edges"][node.canonical_label]
            if edge.tail_node is None:
                self.assertEqual(ref_edge["tail_node"], "None")
            else:
                self.assertEqual(edge.tail_node.canonical_label, ref_edge["tail_node"])
            self.assertEqual(edge.head_node.canonical_label, ref_edge["head_node"])
            self.assertAlmostEqual(node.edge.length, float(ref_edge["length"]), 3)

            # This hackery because NEWICK/NEXUS cannot distinguish between
            # node and edge comments, and everything gets lumped in as a
            # node comment
            if not is_distinct_nodes_and_edges_representation:
                node.comments += edge.comments
                d = {
                        "comments": ref_node["comments"] + ref_edge["comments"],
                        "metadata_comments": ref_node["metadata_comments"] + ref_edge["metadata_comments"],
                        }
                self.check_comments(node, d, is_metadata_extracted)
                if is_metadata_extracted:
                    obs_tuples = []
                    for o in (node, edge):
                        for a in o.annotations:
                            # print("++ {}: {} = {} ({})".format(type(o), a.name, a.value, type(a.value)))
                            v = a.value
                            if isinstance(v, list):
                                v = tuple(v)
                            obs_tuples.append( (a.name, v) )
                    exp_tuples = []
                    for idx, o in enumerate((ref_node["metadata"], ref_edge["metadata"])):
                        for k in o:
                            v = o[k]
                            # print("-- {}{}: {} = {}".format(type(o), idx+1, k, v))
                            if isinstance(v, list):
                                if is_coerce_metadata_values_to_string:
                                    v = tuple(str(vx) for vx in v)
                                else:
                                    v = tuple(v)
                            elif is_coerce_metadata_values_to_string:
                                v = str(v)
                            # print("-- {}{}: {} = {} ({})".format(type(o), idx+1, k, v, type(v)))
                            exp_tuples.append( (k, v) )
                    self.assertEqualUnorderedSequences(tuple(obs_tuples), tuple(exp_tuples))
            else:
                if self.__class__.is_check_comments:
                    self.check_comments(node, ref_node, is_metadata_extracted)
                    self.check_comments(edge, ref_edge, is_metadata_extracted)
                if self.__class__.is_metadata_extracted:
                    self.check_metadata_annotations(
                            item=node,
                            reference=ref_node,
                            is_coerce_metadata_values_to_string=is_coerce_metadata_values_to_string)
                    self.check_metadata_annotations(
                            item=edge,
                            reference=ref_edge,
                            is_coerce_metadata_values_to_string=is_coerce_metadata_values_to_string)
        self.assertEqual(len(visited_nodes), len(ref_tree["nodeset"]))
        if self.__class__.is_taxa_managed_separately_from_tree:
            self.assertEqual(len(obs_taxa), len(tree.taxon_namespace))
            self.assertEqual(set(obs_taxa), set(tree.taxon_namespace))
            obs_node_labels.extend([t.label for t in tree.taxon_namespace])
        elif not self.__class__.is_distinct_taxa_and_labels_on_tree:
            # node labels may have been interpreted as taxa depending on read mode
            obs_node_labels.extend([t.label for t in tree.taxon_namespace if t.label not in obs_node_labels])
        self.assertEqual(len(obs_node_labels), len(ref_tree["nodeset"]))
        self.assertEqual(set(obs_node_labels), set(ref_tree["nodeset"]))

    def verify_standard_trees(self,
            tree_list,
            tree_file_title,
            tree_offset=0):
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        if tree_offset < 0:
            if abs(tree_offset) > expected_number_of_trees:
                tree_offset = 0
            else:
                tree_offset = expected_number_of_trees + tree_offset
        self.assertEqual(len(tree_list), expected_number_of_trees-tree_offset)
        # for tree_idx, (tree, ref_tree) in enumerate(zip(tree_list, tree_directory[tree_file_title])):
        for tree_idx, tree in enumerate(tree_list):
            _LOG.debug("{}: {}".format(tree_file_title, tree_idx))
            self.assertIs(tree.taxon_namespace, tree_list.taxon_namespace)
            self.compare_to_reference_tree(
                    tree=tree,
                    tree_file_title=tree_file_title,
                    reference_tree_idx=tree_idx + tree_offset)

class NewickTestTreesChecker(StandardTestTreesChecker):

    @staticmethod
    def create_class_fixtures(cls,
            suppress_internal_node_taxa=True,
            suppress_leaf_node_taxa=False,
            is_metadata_extracted=False,
            is_coerce_metadata_values_to_string=True,
            is_taxa_managed_separately_from_tree=False,
            is_check_comments=True):
        cls.schema = "newick"
        cls.schema_tree_filepaths = copy.deepcopy(_TREE_FILEPATHS[cls.schema])
        cls.tree_references = copy.deepcopy(_TREE_REFERENCES)
        for tree_file_title in cls.tree_references:
            for reference_tree_idx in range(cls.tree_references[tree_file_title]["num_trees"]):
                ref_tree = cls.tree_references[tree_file_title][str(reference_tree_idx)]
                for ref_node_label in ref_tree["nodes"]:
                    ref_node = ref_tree["nodes"][ref_node_label]
                    ref_node_taxon_label = ref_node["taxon_label"]
                    # print(">>> {} <<<".format(ref_node_taxon_label))
                    # if ref_node["children"] and not suppress_internal_node_taxa:
                    #     ref_node["label"] = ref_node_taxon_label
                    #     ref_node["taxon_label"] = "None"
                    # if not ref_node["children"] and not suppress_leaf_node_taxa:
                    #     ref_node["label"] = ref_node_taxon_label
                    #     ref_node["taxon_label"] = "None"
                for ref_edge_label in ref_tree["edges"]:
                    ref_edge = ref_tree["edges"][ref_edge_label]
                    ref_edge["label"] = "None"
        cls.suppress_internal_node_taxa = suppress_internal_node_taxa
        cls.suppress_leaf_node_taxa = suppress_leaf_node_taxa
        cls.is_metadata_extracted = is_metadata_extracted
        cls.is_coerce_metadata_values_to_string = is_coerce_metadata_values_to_string
        cls.is_taxa_managed_separately_from_tree = is_taxa_managed_separately_from_tree
        cls.is_check_comments = is_check_comments

    def preprocess_tree_to_be_checked(self, tree):
        for nd in tree:
            if nd.is_internal():
                if self.__class__.suppress_internal_node_taxa:
                    self.assertIs(nd.taxon, None)
                else:
                    self.assertIsNot(nd.taxon, None)
            else:
                if self.__class__.suppress_leaf_node_taxa:
                    self.assertIs(nd.taxon, None)
                else:
                    self.assertIsNot(nd.taxon, None)
            if nd.taxon is not None:
                nd.label = nd.taxon.label
            nd.edge.label = nd.label

class NexmlTestTreesChecker(StandardTestTreesChecker):

    @staticmethod
    def create_class_fixtures(cls):
        cls.schema = "nexml"
        cls.schema_tree_filepaths = copy.deepcopy(_TREE_FILEPATHS[cls.schema])
        cls.tree_references = copy.deepcopy(_TREE_REFERENCES)
        cls.is_coerce_metadata_values_to_string = False
        cls.is_taxa_managed_separately_from_tree = True
        cls.is_check_comments = False

class StandardTreeParsingTestCase(StandardTestTreesChecker):

    def test_default_get(self):
        for tree_file_title in [
            "dendropy-test-trees-multifurcating-rooted",
            "dendropy-test-trees-multifurcating-unrooted",
            "dendropy-test-trees-n10-rooted-treeshapes",
            "dendropy-test-trees-n14-unrooted-treeshapes",
                ]:
            tree_filepath = self.schema_tree_filepaths[tree_file_title]
            with open(tree_filepath, "r") as src:
                tree_string = src.read()
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        (dendropy.TreeList.get_from_path, tree_filepath),
                        (dendropy.TreeList.get_from_stream, tree_stream),
                        (dendropy.TreeList.get_from_string, tree_string),
                        )
                for method, src in approaches:
                    tree_list = method(src, self.__class__.schema)
                    # for taxon in tree_list.taxon_namespace:
                    #     print(taxon)
                    self.verify_standard_trees(tree_list=tree_list,
                            tree_file_title=tree_file_title)

    def test_default_read(self):
        preloaded_tree_file_title = "dendropy-test-trees-n33-unrooted-x10a"
        preloaded_tree_reference = self.tree_references[preloaded_tree_file_title]
        tree_file_title = "dendropy-test-trees-n33-unrooted-x10a"
        tree_reference = self.tree_references[tree_file_title]
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                # prepopulate
                tree_list = dendropy.TreeList.get_from_path(
                        self.schema_tree_filepaths[preloaded_tree_file_title],
                        self.__class__.schema)
                # check to make sure trees were loaded
                old_len = len(tree_list)
                self.assertEqual(old_len, len(tree_list._trees))
                self.assertEqual(old_len, preloaded_tree_reference["num_trees"])
                self.verify_standard_trees(tree_list, preloaded_tree_file_title)

                # load
                old_id = id(tree_list)
                f = getattr(tree_list, method)
                trees_read = f(src, self.__class__.schema)
                new_id = id(tree_list)
                self.assertEqual(old_id, new_id)

                # make sure new trees added
                new_len = len(tree_list)
                self.assertEqual(new_len, len(tree_list._trees))
                expected_number_of_trees = tree_reference["num_trees"]
                self.assertEqual(old_len + expected_number_of_trees, new_len)
                self.assertEqual(trees_read, expected_number_of_trees)

                # check new trees
                for tree_idx, tree in enumerate(tree_list[old_len:]):
                    self.compare_to_reference_tree(
                            tree=tree,
                            tree_file_title=tree_file_title,
                            reference_tree_idx=tree_idx)

                # make sure old ones still intact
                for tree_idx, tree in enumerate(tree_list[:old_len]):
                    self.compare_to_reference_tree(
                            tree=tree,
                            tree_file_title=preloaded_tree_file_title,
                            reference_tree_idx=tree_idx)

    def test_tree_offset_get(self):
        tree_file_title = "dendropy-test-trees-n33-unrooted-x100a"
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        tree_offsets = set([0, expected_number_of_trees-1, -1, -expected_number_of_trees])
        while len(tree_offsets) < 8:
            tree_offsets.add(random.randint(1, expected_number_of_trees-2))
        while len(tree_offsets) < 12:
            tree_offsets.add(random.randint(-expected_number_of_trees-2, -2))
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        for tree_offset in tree_offsets:
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        (dendropy.TreeList.get_from_path, tree_filepath),
                        (dendropy.TreeList.get_from_stream, tree_stream),
                        (dendropy.TreeList.get_from_string, tree_string),
                        )
                for method, src in approaches:
                    tree_list = method(
                            src,
                            self.__class__.schema,
                            collection_offset=0,
                            tree_offset=tree_offset)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title,
                            tree_offset=tree_offset)

    def test_tree_offset_read(self):
        tree_file_title = "dendropy-test-trees-n33-unrooted-x100a"
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        tree_offsets = set([0, expected_number_of_trees-1, -1, -expected_number_of_trees])
        while len(tree_offsets) < 8:
            tree_offsets.add(random.randint(1, expected_number_of_trees-2))
        while len(tree_offsets) < 12:
            tree_offsets.add(random.randint(-expected_number_of_trees-2, -2))
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        for tree_offset in tree_offsets:
            with open(tree_filepath, "r") as tree_stream:
                approaches = (
                        ("read_from_path", tree_filepath),
                        ("read_from_stream", tree_stream),
                        ("read_from_string", tree_string),
                        )
                for method, src in approaches:
                    tree_list = dendropy.TreeList()
                    f = getattr(tree_list, method)
                    trees_read = f(src,
                            self.__class__.schema,
                            collection_offset=0,
                            tree_offset=tree_offset)
                    self.verify_standard_trees(
                            tree_list=tree_list,
                            tree_file_title=tree_file_title,
                            tree_offset=tree_offset)

    def test_tree_offset_without_collection_offset_get(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        approaches = (
                dendropy.TreeList.get_from_path,
                dendropy.TreeList.get_from_stream,
                dendropy.TreeList.get_from_string,
                )
        for approach in approaches:
            with self.assertRaises(TypeError):
                approach(tree_filepath, self.__class__.schema, collection_offset=None, tree_offset=0)

    def test_tree_offset_without_collection_offset_read(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        approaches = (
                "read_from_path",
                "read_from_stream",
                "read_from_string",
                )
        for approach in approaches:
            tree_list = dendropy.TreeList()
            f = getattr(tree_list, approach)
            with self.assertRaises(TypeError):
                f(tree_filepath, self.__class__.schema, collection_offset=None, tree_offset=0)

    def test_out_of_range_tree_offset_get(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(IndexError):
                    method(src, self.__class__.schema, collection_offset=0, tree_offset=expected_number_of_trees)

    def test_out_of_range_tree_offset_read(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        tree_reference = self.tree_references[tree_file_title]
        expected_number_of_trees = tree_reference["num_trees"]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(IndexError):
                    f(src, self.__class__.schema, collection_offset=0, tree_offset=expected_number_of_trees)

    def test_out_of_range_collection_offset_get(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(IndexError):
                    method(src, self.__class__.schema, collection_offset=1, tree_offset=0)

    def test_out_of_range_collection_offset_read(self):
        tree_file_title = 'dendropy-test-trees-n33-unrooted-x10a'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(IndexError):
                    f(src, self.__class__.schema, collection_offset=1, tree_offset=0)

    def test_unsupported_keyword_arguments_get(self):
        tree_file_title = 'dendropy-test-trees-n12-x2'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    (dendropy.TreeList.get_from_path, tree_filepath),
                    (dendropy.TreeList.get_from_stream, tree_stream),
                    (dendropy.TreeList.get_from_string, tree_string),
                    )
            for method, src in approaches:
                with self.assertRaises(TypeError):
                    method(src,
                            self.__class__.schema,
                            suppress_internal_taxa=True,  # should be suppress_internal_node_taxa
                            gobbledegook=False,
                            )

    def test_unsupported_keyword_arguments_read(self):
        tree_file_title = 'dendropy-test-trees-n12-x2'
        tree_filepath = self.schema_tree_filepaths[tree_file_title]
        with open(tree_filepath, "r") as src:
            tree_string = src.read()
        with open(tree_filepath, "r") as tree_stream:
            approaches = (
                    ("read_from_path", tree_filepath),
                    ("read_from_stream", tree_stream),
                    ("read_from_string", tree_string),
                    )
            for method, src in approaches:
                tree_list = dendropy.TreeList()
                f = getattr(tree_list, method)
                with self.assertRaises(TypeError):
                    f(src,
                      self.__class__.schema,
                      suppress_internal_taxa=True,  # should be suppress_internal_node_taxa
                      gobbledegook=False,
                    )

