#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Split calculation and management.
DEPRECATED IN DENDROPY 4.
"""

import dendropy
from dendropy.utility import bitprocessing
from dendropy.utility import deprecate

def tree_from_splits(splits, taxon_set, is_rooted=False):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.tree_from_splits()'.",
            old_construct="from dendropy import treesplit\ntree = treesplit.tree_from_splits(...)",
            new_construct="import dendropy\ntree = dendropy.Tree.from_split_bitmasks(...)")
    return dendropy.Tree.from_split_bitmasks(
            split_bitmasks=splits,
            taxon_namespace=taxon_set,
            is_rooted=is_rooted)

def split_to_list(s, mask=-1, one_based=False, ordination_in_mask=False):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.split_to_list()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.split_to_list(...)",
            new_construct="from dendropy.utility import bitprocessing\nd = bitprocessing.indexes_of_set_bits(...)")
    return bitprocessing.indexes_of_set_bits(
            s=s,
            fill_bitmask=mask,
            one_based=one_based,
            ordination_in_mask=ordination_in_mask)

def is_trivial_split(split, mask):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.is_trivial_split()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.is_trivial_split(...)",
            new_construct="import dendropy\nd = dendropy.Bipartition.is_trivial_bitmask(...)")
    return dendropy.Bipartition.is_trivial_bitmask(bitmask=split, fill_bitmask=mask)

def is_non_singleton_split(split, mask):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.is_non_singleton_split()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.is_non_singleton_split(...)",
            new_construct="import dendropy\nd = not dendropy.Bipartition.is_trivial_bitmask(...)")
    return dendropy.Bipartition.is_trivial_bitmask(bitmask=split, fill_bitmask=mask)

def split_as_string(split_mask, width, symbol1=None, symbol2=None):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.split_as_string()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.split_as_string(...)",
            new_construct="""\
# if using a bipartition
d = bipartition.split_as_bitstring(...)
d = bipartition.leafset_as_bitstring(...)
# if a "raw" bitmask
from dendropy.utility import bitprocessing
d = bitprocessing.int_as_bitstring(...)""")
    return bitprocessing.int_as_bitstring(
            n=split_mask,
            length=width,
            symbol0=symbol1,
            symbol1=symbol2,
            reverse=False)

def split_as_string_rev(split_mask, width, symbol1='.', symbol2='*'):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.split_as_string_rev()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.split_as_string(...)",
            new_construct="""\
# if using a bipartition
d = bipartition.split_as_bitstring(...)[::-1]
d = bipartition.leafset_as_bitstring(...)[::-1]
# if a "raw" bitmask
from dendropy.utility import bitprocessing
d = bitprocessing.int_as_bitstring(..., reverse=True)""")
    return bitprocessing.int_as_bitstring(
            mask=split_mask,
            length=width,
            symbol0=symbol0,
            symbol1=symbol1,
            reverse=True)

def find_edge_from_split(root, split_to_find, mask=-1):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.find_edge_from_split()'",
            old_construct="from dendropy import treesplit\nd = treesplit.find_edge_from_split(...)",
            new_construct="""\
# if using a bipartition
d = bipartition.edge
# if a "raw" bitmask
d = tree.find_edge_for_split_bitmask(...)""")
    return root.find_edge_for_split_bitmask(split_to_find, fill_bitmask=mask)

def encode_splits(tree, create_dict=True, suppress_unifurcations=True):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.encode_splits()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.encode_splits(tree, ...)",
            new_construct="bipartitions = tree.encode_bipartitions(...)\nsplit_bitmasks = tree.split_bitmask_edge_map.keys()")
    return tree.encode_bipartitions(suppress_unifurcations=suppress_unifurcations)

def is_compatible(split1, split2, mask):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.is_compatible()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.is_compatible(...)",
            new_construct="""\
# if using a bipartition
d = bipartition.is_compatible_with(other_bipartition)
# if a "raw" bitmask
d = dendropy.Bipartition.is_compatible_bitmasks(m1, m2, fill_bitmask=mask)""")
    """
    Mask should have 1 for every leaf in the leaf_set
    """
    # m1 = mask & split1
    # m2 = mask & split2
    # if 0 == (m1 & m2):
    #     return True
    # c2 = mask ^ split2
    # if 0 == (m1 & c2):
    #     return True
    # c1 = mask ^ split1
    # if 0 == (c1 & m2):
    #     return True
    # if 0 == (c1 & c2):
    #     return True
    # return False
    return dendropy.Bipartition.is_compatible_bitmasks(split1, split2, fill_bitmask=mask)

def delete_outdegree_one(tree):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.delete_outdegree_one()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.delete_outdegree_one(tree)",
            new_construct="tree.suppress_unifurcations()")
    return tree.suppress_unifurcations()

def lowest_bit_only(s):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.lowest_bit_only()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.lowest_bit_only(...)",
            new_construct="from dendropy.utility import bitprocessing\nd = bitprocessing.least_significant_set_bit(...)")
    return bitprocessing.least_significant_set_bit(s)

def count_bits(a):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: 'dendropy.treesplit.count_bits()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.count_bits(...)",
            new_construct="from dendropy.utility import bitprocessing\nd = bitprocessing.num_set_bits(...)")
    return bitprocessing.num_set_bits(a)

class SplitDistribution(dendropy.SplitDistribution):
    def __init__(self, taxon_set=None, split_set=None):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.SplitDistribution' class has moved to 'dendropy.calculate.treesplit.SplitDistribution'.",
                old_construct="from dendropy import treesplit\nm = treesplit.SplitDistribution(...)",
                new_construct="import dendropy\nm = dendropy.SplitDistribution(...)")
        dendropy.SplitDistribution.__init__(self,
                taxon_namespace=taxon_set)
        if split_set:
            for split in split_set:
                self.add_split_count(split, count=1)




