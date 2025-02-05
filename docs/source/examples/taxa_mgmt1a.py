import warnings
import dendropy

warnings.warn(
    "This example is known to be broken! "
    "It will be fixed or removed in the future. "
    "See https://github.com/jeetsukumaran/DendroPy/issues/160 for details. "
    "Patch contributions are welcome.",
)

tree_str1 = "((A,B),C);"

tree1 = dendropy.Tree.get(data=tree_str1, schema="newick")
tree2 = dendropy.Tree.get(data=tree_str1, schema="newick")
print(tree1.taxon_namespace is  tree2.taxon_namespace) # False
for nd1, nd2 in zip(tree1, tree2):
    assert nd1.taxon is nd2.taxon # Assertion Error

