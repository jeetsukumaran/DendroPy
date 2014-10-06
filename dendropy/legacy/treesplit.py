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

"""
Split calculation and management.
DEPRECATED IN DENDROPY 4: USE `dendropy.calculate.treesplit` instead.
"""

from dendropy.calculate import treesplit
from dendropy.utility import deprecate

def tree_from_splits(splits, taxon_set, is_rooted=False):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.tree_from_splits()' function has moved to 'dendropy.calculate.treesplit.weighted_tree_from_splits()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.tree_from_splits(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.tree_from_splits(...)")
    return treesplit.tree_from_splits(splits=splits,
            taxon_set=taxon_namespace,
            is_rooted=is_rooted)

def split_to_list(s, mask=-1, one_based=False, ordination_in_mask=False):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.split_to_list()' function has moved to 'dendropy.calculate.treesplit.split_to_list()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.split_to_list(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.split_to_list(...)")
    return treesplit.split_to_list(
            s=s,
            mask=mask,
            one_based=one_based,
            ordination_in_mask=ordination_in_mask)

def is_trivial_split(split, mask):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.is_trivial_split()' function has moved to 'dendropy.calculate.treesplit.is_trivial_split()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.is_trivial_split(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.is_trivial_split(...)")
    return treesplit.is_trivial_split(split=split, mask=mask)

def is_non_singleton_split(split, mask):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.is_non_singleton_split()' function has moved to 'dendropy.calculate.treesplit.is_non_singleton_split()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.is_non_singleton_split(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.is_non_singleton_split(...)")
    return treesplit.is_non_singleton_split(split=split, mask=mask)

def split_as_string(split_mask, width, symbol1=None, symbol2=None):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.split_as_string()' function has moved to 'dendropy.calculate.treesplit.split_as_string()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.split_as_string(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.split_as_string(...)")
    return treesplit.split_as_string(
            split_mask=split_mask,
            width=width,
            symbol1=symbol1,
            symbol2=symbol2)

def split_as_string_rev(split_mask, width, symbol1='.', symbol2='*'):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.split_as_string_rev()' function has moved to 'dendropy.calculate.treesplit.split_as_string_rev()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.split_as_string_rev(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.split_as_string_rev(...)")
    return treesplit.split_as_string_rev(
            split_mask=split_mask,
            width=width,
            symbol1=symbol1,
            symbol2=symbol2)

def find_edge_from_split(root, split_to_find, mask=-1):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.find_edge_from_split()' function has moved to 'dendropy.calculate.treesplit.find_edge_from_split()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.find_edge_from_split(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.find_edge_from_split(...)")
    return treesplit.find_edge_from_split(
            root=root,
            split_to_find=split_to_find,
            mask=mask)

def encode_splits(tree, create_dict=True, delete_outdegree_one=True):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.encode_splits()' function has moved to 'dendropy.calculate.treesplit.encode_splits()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.encode_splits(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.encode_splits(...)")
    return treesplit.encode_splits(
            tree=tree,
            create_dict=create_dict,
            delete_outdegree_one=delete_outdegree_one)

def is_compatible(split1, split2, mask):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.is_compatible()' function has moved to 'dendropy.calculate.treesplit.is_compatible()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.is_compatible(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.is_compatible(...)")
    return treesplit.is_compatible(
            split1=split1,
            split2=split2,
            mask=mask)

def delete_outdegree_one(tree):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.delete_outdegree_one()' function has moved to 'dendropy.calculate.treesplit.delete_outdegree_one()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.delete_outdegree_one(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.delete_outdegree_one(...)")
    return treesplit.delete_outdegree_one(tree=tree)

def lowest_bit_only(s):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.lowest_bit_only()' function has moved to 'dendropy.calculate.treesplit.lowest_bit_only()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.lowest_bit_only(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.lowest_bit_only(...)")
    return treesplit.lowest_bit_only(split=s)

def count_bits(a):
    deprecate.dendropy_deprecation_warning(
            preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.count_bits()' function has moved to 'dendropy.calculate.treesplit.count_bits()'.",
            old_construct="from dendropy import treesplit\nd = treesplit.count_bits(...)",
            new_construct="from dendropy.calculate import treesplit\nd = treesplit.count_bits(...)")
    return treesplit.count_bits(split=a)

class SplitDistribution(treesplit.SplitDistribution):
    def __init__(self, taxon_set=None, split_set=None):
        deprecate.dendropy_deprecation_warning(
                preamble="Deprecated since DendroPy 4: The 'dendropy.treesplit.SplitDistribution' class has moved to 'dendropy.calculate.treesplit.SplitDistribution'.",
                old_construct="from dendropy import treesplit\nm = treesplit.SplitDistribution(...)",
                new_construct="from dendropy.calculate import treesplit\nm = treesplit.SplitDistribution(...)")
        treesplit.SplitDistribution.__init__(self,
                taxon_namespace=taxon_set)
        if split_set:
            for split in split_set:
                self.add_split_count(split, count=1)




