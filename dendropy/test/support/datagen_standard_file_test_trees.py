#! /usr/bin/env python

import json
from dendropy.test.support import pathmap

tree_file_titles = [
    'standard-test-trees-n14-unrooted-treeshapes',
    'standard-test-trees-n10-rooted-treeshapes',
    'standard-test-trees-n12-x2',
    'standard-test-trees-n33-x100a',
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

