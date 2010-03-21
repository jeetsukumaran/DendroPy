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
from dendropy.dataio import ioclient
from dendropy.dataio import newick
from dendropy.dataio import nexusreader_py
from dendropy.dataio import nexuswriter
from dendropy.dataio import fasta
from dendropy.dataio import phylip
from dendropy.dataio import nexml
from dendropy.dataio.ioclient import get_reader, get_writer, tree_source_iter, multi_tree_source_iter

_LOG = messaging.get_logger(__name__)

###############################################################################
## Data Schema Handlers
##
## Syntax is:
##   ioclient.register(<FORMAT NAME>, <READER TYPE>, <WRITER TYPE>, <TREE ITERATOR>)
##
ioclient.register("nexus", nexusreader_py.NexusReader, nexuswriter.NexusWriter, nexusreader_py.tree_source_iter)
ioclient.register("nexus-native", nexusreader_py.NexusReader, nexuswriter.NexusWriter, nexusreader_py.tree_source_iter)
ioclient.register("newick", newick.NewickReader, newick.NewickWriter, newick.tree_source_iter)
ioclient.register("nexus/newick", None, None, nexusreader_py.generalized_tree_source_iter)
ioclient.register("fasta", fasta.FastaReader, fasta.FastaWriter, None)
ioclient.register("dnafasta", fasta.DNAFastaReader, fasta.FastaWriter, None)
ioclient.register("rnafasta", fasta.RNAFastaReader, fasta.FastaWriter, None)
ioclient.register("proteinfasta", fasta.ProteinFastaReader, fasta.FastaWriter, None)
ioclient.register("phylip", phylip.PhylipReader, phylip.PhylipWriter, None)
ioclient.register("nexml", nexml.NexmlReader, nexml.NexmlWriter, None)

###############################################################################
## NEXUS Parser Implementation Selection
##

def disable_ncl():
    _LOG.debug('Disabling Nexus Class Library bindings: using native Python NEXUS parser')
    from dendropy.dataio import nexusreader_py
    ioclient.register("nexus", nexusreader_py.NexusReader, nexuswriter.NexusWriter, nexusreader_py.tree_source_iter)

def enable_ncl():
    from dendropy.dataio import nexusreader_ncl
    if nexusreader_ncl.DENDROPY_NCL_AVAILABILITY:
        _LOG.debug('Enabling Nexus Class Library bindings: using NCL NEXUS parser')
        ioclient.register("nexus", nexusreader_ncl.NexusReader, nexuswriter.NexusWriter, nexusreader_ncl.tree_source_iter)
    else:
        _LOG.debug('Nexus Class Library bindings are not available: using native Python NEXUS parser')

if "DENDROPY_ENABLE_NCL" in os.environ:
    enable_ncl()
else:
    disable_ncl()
