#! /usr/bin/env python

import dendropy

tree0 = dendropy.Tree.get(
        path="pythonidae.mle.nex",
        schema="nexus")

node_filter_fn = lambda nd: nd.taxon is None or nd.taxon.label.startswith("Morelia")
tree1 = tree0.extract_tree(node_filter_fn=node_filter_fn)
print(tree1.as_string("newick"))
