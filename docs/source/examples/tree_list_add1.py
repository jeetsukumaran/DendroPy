import warnings
import dendropy

warnings.warn(
    "This example is known to be broken! "
    "It will be fixed or removed in the future. "
    "See https://github.com/jeetsukumaran/DendroPy/issues/160 for details. "
    "Patch contributions are welcome.",
)

trees = dendropy.TreeList()
trees.read(path="sometrees.nex", schema="nexus", tree_offset=10)
trees.read(data="(A,(B,C));((A,B),C);", schema="newick")

