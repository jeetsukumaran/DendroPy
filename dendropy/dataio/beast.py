#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

"""
Readers/writers for BEAST data.
"""

import re
from dendropy.dataio.nexusreader_py import NexusReader

BEAST_NODE_INFO_PATTERN = re.compile(r'(.+?)=({.+?,.+?}|.+?)(,|$)')

###############################################################################
## BeastNexusReader

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

class BeastSummaryTreeReader(NexusReader):
    """
    Encapsulates parsing of BEAST summary trees (NEXUS with supplementary
    information in node comments)
    """

    def parse_beast_tree_node_info(tree,
        set_node_attributes=True,
        value_type=float,
        create_field_if_missing=True):
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
            # grab comment
            if nd.comments is None or len(nd.comments) == 0:
                raise ValueError("No comments found associated with node '%s'" % (str(nd)))
            node_comment = nd.comments[0][1:]
            # populate info dictionary
            beast_info = {}
            for match_group in BEAST_NODE_INFO_PATTERN.findall(node_comment):
                key, val = match_group[:2]
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
                    k = k.replace('95%_HPD', '95hpd')
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
        """
        if 'edge_len_type' not in kwargs:
            kwargs['edge_len_type'] = float
        NexusReader.__init__(self, **kwargs)

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
                        create_field_if_missing=True)
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
        for tree in NexusReader.tree_source_iter(source):
            yield tree
