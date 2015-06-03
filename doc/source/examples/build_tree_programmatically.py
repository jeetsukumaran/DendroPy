import dendropy

taxon_namespace = dendropy.TaxonNamespace(["A", "B", "C", "D",])
tree = dendropy.Tree(taxon_namespace=taxon_namespace)


# Create and add a new child node to the seed node,
# assigning it an edge length:
#
#     (seed)
#      /
#     /
#    ch1
#
ch1 = tree.seed_node.new_child()
ch1.edge.length = 1

# Can also assign edge length on construction:
#
#     (seed)
#      / \
#     /   \
#   ch1   ch2
#
ch2 = tree.seed_node.new_child(edge_length=1)

# Can also add an existing node as child
#
#       (seed)
#       /   \
#      /     \
#    ch1     ch2
#   /  \     /  \
#  ch3 ch4  ch5 ch6
ch3 = dendropy.Node(edge_length=1)
ch4 = dendropy.Node(edge_length=2)
ch1.add_child(ch3)
ch1.add_child(ch4)
ch5 = dendropy.Node(edge_length=1)
ch6 = dendropy.Node(edge_length=2)
# Note: this clears/deletes existing child nodes before adding the new ones;
ch2.set_child_nodes([ch5, ch6])

# Assign taxa
ch3.taxon = taxon_namespace.get_taxon("A")
ch4.taxon = taxon_namespace.get_taxon("B")
ch5.taxon = taxon_namespace.get_taxon("C")
ch6.taxon = taxon_namespace.get_taxon("D")

print(tree.as_string("newick"))
print(tree.as_ascii_plot())
