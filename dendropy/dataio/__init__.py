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
Infrastructure for phylogenetic data object serialization and deserialization.
Provides support for reading/parsing and formatting/writing phylogenetic data
in various formats.
"""

import os
from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)
from dendropy.utility import error
from dendropy.dataio import ioclient
from dendropy.dataio import newick
from dendropy.dataio import nexus
from dendropy.dataio import fasta
from dendropy.dataio import phylip
from dendropy.dataio import nexml

# syntax is:
#   dataschema.register(<FORMAT NAME>, <READER TYPE>, <WRITER TYPE>, <TREE ITERATOR>)

DEFAULT_NEXUS_READER = None

if "DENDROPY_BYPASS_NCL" not in os.environ:
    from dendropy.dataio import ncl
    if ncl.DENDROPY_NCL_AVAILABILITY:
        _LOG.debug("Using NCL as the default NEXUS parser")
        DEFAULT_NEXUS_READER = ncl.NexusReader

if DEFAULT_NEXUS_READER is None:
    _LOG.debug("NCL wrapper is not available or disabled: using native Python NEXUS parser")
    DEFAULT_NEXUS_READER = nexus.NexusReader

ioclient.register("nexus", DEFAULT_NEXUS_READER, nexus.NexusWriter, nexus.tree_source_iter)
ioclient.register("nexus-native", nexus.NexusReader, nexus.NexusWriter, nexus.tree_source_iter)
ioclient.register("newick", newick.NewickReader, newick.NewickWriter, newick.tree_source_iter)
ioclient.register("nexus/newick", None, None, nexus.generalized_tree_source_iter)
ioclient.register("fasta", fasta.FastaReader, fasta.FastaWriter, None)
ioclient.register("dnafasta", fasta.DNAFastaReader, fasta.FastaWriter, None)
ioclient.register("rnafasta", fasta.RNAFastaReader, fasta.FastaWriter, None)
ioclient.register("proteinfasta", fasta.ProteinFastaReader, fasta.FastaWriter, None)
ioclient.register("phylip", phylip.PhylipReader, phylip.PhylipWriter, None)
ioclient.register("nexml", nexml.NexmlReader, nexml.NexmlWriter, None)

from dendropy.dataio.ioclient import get_reader, get_writer, tree_source_iter, multi_tree_source_iter

