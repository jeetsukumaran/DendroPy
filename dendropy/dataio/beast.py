#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
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
Readers/writers for BEAST data.
"""

import re
from dendropy.dataio import ioclient
from dendropy.dataio.nexusreader_py import NexusReader

###############################################################################
## Constants

BEAST_NODE_INFO_PATTERN = re.compile(r'(.+?)=({.+?,.+?}|.+?)(,|$)')

BEAST_SUMMARY_TREE_NODE_FIELDS = [
        'height',
        'height_median',
        'height_95hpd',
        'height_range',
        'length',
        'length_median',
        'length_95hpd',
        'length_range',
        'posterior']

BEAST_SUMMARY_FIELDS_TO_ATTR_MAP = {
        'height_95%_HPD' : 'height_95hpd',
        'length_95%_HPD' : 'length_95hpd',
        }

###############################################################################
## tree_source_iter

def summary_tree_source_iter(stream, **kwargs):
    """
    Iterates over a NEXUS-formatted source of trees given by file-like object
    `stream`

    The following optional keyword arguments are recognized:

        - `taxon_set`: TaxonSet object to use when reading data
        - `as_rooted=True` (or `as_unrooted=False`): interprets trees as rooted
        - `as_unrooted=True` (or `as_rooted=False`): interprets trees as unrooted
        - `default_as_rooted=True` (or `default_as_unrooted=False`): interprets
           all trees as rooted if rooting not given by `[&R]` or `[&U]` comments
        - `default_as_unrooted=True` (or `default_as_rooted=False`): interprets
           all trees as rooted if rooting not given by `[&R]` or `[&U]` comments
        - `edge_len_type`: specifies the type of the edge lengths (int or float)
        - `encode_splits`: specifies whether or not split bitmasks will be
           calculated and attached to the edges.
        - `finish_node_func`: is a function that will be applied to each node
           after it has been constructed
        - `allow_duplicate_taxon_labels` : if True, allow duplicate labels
           on trees
        - `ignore_missing_node_info` : if True, then no errors will be thrown if
           tree nodes do not have the the required information.

    Only trees will be returned, and any and all character data will
    be skipped. The iterator will span over multiple tree blocks,
    but, because our NEXUS data model implementation currently does
    not recognize multiple taxon collection definnitions, taxa in
    those tree blocks will be aggregated into the same `TaxonSet` (a
    new one created, or the one passed to this method via the
    `taxon_set` argument). This behavior is similar to how multiple
    tree blocks are handled by a full NEXUS data file read.
    """
    reader = ioclient.get_reader('beast-summary-tree', **kwargs)
    for i, tree in enumerate(reader.tree_source_iter(stream)):
        yield tree

###############################################################################
## BeastNexusReader

class BeastSummaryTreeReader(NexusReader):
    """
    Encapsulates parsing of BEAST summary trees (NEXUS with supplementary
    information in node comments)
    """

    def parse_beast_tree_node_info(tree,
        set_node_attributes=True,
        value_type=float,
        create_field_if_missing=True,
        ignore_missing=False):
        """
        Parses the comment tokens associated with nodes in a BEAST summary tree,
        creating an attribute for each node on the tree, `beast_info`, which is a
        dictionary of key-value pairs. Multiple values (e.g., as given for HPD ranges, etc.)
        will be converted to a tuples.
        If `value_type` is not None, all values will be coerced to this type.
        If `set_node_attributes` is True, then all fields will be added as
        attributes of nodes (with some name-mapping to ensure legal Python names).
        If `create_field_if_missing` is True, then if any fields that expected but not
        found will be automatically created (and set to None).
        """
        for nd in tree.postorder_node_iter():
            beast_info = {}
            if nd.comments is None or len(nd.comments) == 0:
                if not ignore_missing:
                    raise ValueError("No comments found associated with node '%s'" % (str(nd)))
            else:
                # populate info dictionary
                node_comment = nd.comments[0][1:]
                for match_group in BEAST_NODE_INFO_PATTERN.findall(node_comment):
                    key, val = match_group[:2]
                    key = BEAST_SUMMARY_FIELDS_TO_ATTR_MAP.get(key, key)
                    if val.startswith('{'):
                        if value_type is not None:
                            val = [value_type(v) for v in val[1:-1].split(',')]
                        else:
                            val = val[1:-1].split(',')
                    else:
                        if value_type is not None:
                            val = value_type(val)
                    beast_info[key] = val
            # create missing fields
            if create_field_if_missing:
                for k in BEAST_SUMMARY_TREE_NODE_FIELDS:
                    if k not in beast_info:
                        beast_info[k] = None
            # assign to node
            nd.beast_info = beast_info
            # set attributes
            if set_node_attributes:
                for k,v in nd.beast_info.items():
                    setattr(nd, k, v)
        return tree
    parse_beast_tree_node_info = staticmethod(parse_beast_tree_node_info)

    def __init__(self, **kwargs):
        """
        __init__ recognizes the following keywords (in addition to those of `DataReader.__init__`):

            - `taxon_set`: TaxonSet object to use when reading data
            - `as_rooted=True` (or `as_unrooted=False`): interprets trees as rooted
            - `as_unrooted=True` (or `as_rooted=False`): interprets trees as unrooted
            - `default_as_rooted=True` (or `default_as_unrooted=False`): interprets
               all trees as rooted if rooting not given by `[&R]` or `[&U]` comments
            - `default_as_unrooted=True` (or `default_as_rooted=False`): interprets
               all trees as rooted if rooting not given by `[&R]` or `[&U]` comments
            - `edge_len_type`: specifies the type of the edge lengths (int or float)
            - `encode_splits`: specifies whether or not split bitmasks will be
               calculated and attached to the edges.
            - `finish_node_func`: is a function that will be applied to each node
               after it has been constructed
            - `allow_duplicate_taxon_labels` : if True, allow duplicate labels
               on trees
            - `ignore_missing_node_info` : if True, then no errors will be thrown if
               tree nodes do not have the the required information.
        """
        if 'edge_len_type' not in kwargs:
            kwargs['edge_len_type'] = float
        NexusReader.__init__(self, **kwargs)
        self.is_ignore_missing_node_info = kwargs.get('ignore_missing_node_info', False)

    def read(self, stream):
        """
        Instantiates and returns a DataSet object based on the
        NEXUS-formatted contents given in the file-like object `stream`.
        """
        if self.dataset is None:
            new_tree_list_idx = 0
        else:
            new_tree_list_idx = len(self.dataset.tree_lists)
        self.dataset = NexusReader.read(self, stream)
        for tree_list in self.dataset.tree_lists[new_tree_list_idx:]:
            for tree in tree_list:
                BeastSummaryTreeReader.parse_beast_tree_node_info(tree,
                        set_node_attributes=True,
                        value_type=float,
                        create_field_if_missing=True,
                        ignore_missing=self.is_ignore_missing_node_info)
        return self.dataset

    def tree_source_iter(self, stream):
        """
        Iterates over a NEXUS-formatted source of trees.
        Only trees will be returned, and any and all character data will
        be skipped. The iterator will span over multiple tree blocks,
        but, because our NEXUS data model implementation currently does
        not recognize multiple taxon collection definnitions, taxa in
        those tree blocks will be aggregated into the same `TaxonSet` (a
        new one created, or the one passed to this method via the
        `taxon_set` argument). This behavior is similar to how multiple
        tree blocks are handled by a full NEXUS data file read.
        """
        for tree in NexusReader.tree_source_iter(self, stream):
            BeastSummaryTreeReader.parse_beast_tree_node_info(tree,
                    set_node_attributes=True,
                    value_type=float,
                    create_field_if_missing=True,
                    ignore_missing=self.is_ignore_missing_node_info)
            yield tree
