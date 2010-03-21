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
Specialized NEXUS tree iteration.
"""

from dendropy.dataio import newick
from dendropy.dataio import nexustokenizer
from dendropy.dataio import ioclient

###############################################################################
## tree_source_iter

def tree_source_iter(stream, **kwargs):
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

    Only trees will be returned, and any and all character data will
    be skipped. The iterator will span over multiple tree blocks,
    but, because our NEXUS data model implementation currently does
    not recognize multiple taxon collection definnitions, taxa in
    those tree blocks will be aggregated into the same `TaxonSet` (a
    new one created, or the one passed to this method via the
    `taxon_set` argument). This behavior is similar to how multiple
    tree blocks are handled by a full NEXUS data file read.
    """
    reader = ioclient.get_reader('nexus', **kwargs)
    for i, tree in enumerate(reader.tree_source_iter(stream)):
        yield tree

def generalized_tree_source_iter(stream, **kwargs):
    """
    Diagnoses and handles both NEXUS and NEWICK files.
    """
    stream_tokenizer = nexustokenizer.NexusTokenizer(stream)
    token = stream_tokenizer.read_next_token_ucase()
    schema = None
    if token.upper() == "#NEXUS":
        schema = "nexus"
    else:
        if token == "(":
            schema = "newick"
    try:
        stream_tokenizer.stream_handle.seek(0)
    except IOError:
        raise TypeError("File schema of non-random access source (such as stdin) must be specified in advance.")
    if schema == "nexus":
        return tree_source_iter(stream, **kwargs)
    elif schema == "newick":
        return ioclient.tree_source_iter(stream, 'newick', **kwargs)
    else:
        raise TypeError("Cannot diagnose file schema based on first token found: '%s' (looking for '#NEXUS' or '(')")
