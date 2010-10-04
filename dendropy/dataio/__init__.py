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
Infrastructure for phylogenetic data object serialization and deserialization.
Provides support for reading/parsing and formatting/writing phylogenetic data
in various formats.
"""

import os
from dendropy.utility import messaging
from dendropy.dataio import ioclient
from dendropy.dataio import newick
from dendropy.dataio import nexusreader_py
from dendropy.dataio import nexustreeiter
from dendropy.dataio import nexuswriter
from dendropy.dataio import fasta
from dendropy.dataio import phylip
from dendropy.dataio import nexml
from dendropy.dataio import beast
from dendropy.dataio.ioclient import get_reader, get_writer, tree_source_iter, multi_tree_source_iter

_LOG = messaging.get_logger(__name__)

###############################################################################
## Data Schema Handlers
##
## Syntax is:
##   ioclient.register(<FORMAT NAME>, <READER TYPE>, <WRITER TYPE>, <TREE ITERATOR>)
##
ioclient.register("nexus", nexusreader_py.NexusReader, nexuswriter.NexusWriter, nexustreeiter.tree_source_iter)
ioclient.register("newick", newick.NewickReader, newick.NewickWriter, newick.tree_source_iter)
ioclient.register("nexus/newick", None, None, nexustreeiter.generalized_tree_source_iter)
ioclient.register("fasta", fasta.FastaReader, fasta.FastaWriter, None)
ioclient.register("dnafasta", fasta.DNAFastaReader, fasta.FastaWriter, None)
ioclient.register("rnafasta", fasta.RNAFastaReader, fasta.FastaWriter, None)
ioclient.register("proteinfasta", fasta.ProteinFastaReader, fasta.FastaWriter, None)
ioclient.register("phylip", phylip.PhylipReader, phylip.PhylipWriter, None)
ioclient.register("nexml", nexml.NexmlReader, nexml.NexmlWriter, None)
ioclient.register("beast-summary-tree", beast.BeastSummaryTreeReader, None, beast.summary_tree_source_iter)

###############################################################################
## NEXUS Parser Implementation Selection
##

def disable_ncl():
    _LOG.debug('Disabling Nexus Class Library bindings: using native Python NEXUS parser')
    from dendropy.dataio import nexusreader_py
    ioclient.register("nexus", nexusreader_py.NexusReader, nexuswriter.NexusWriter, nexustreeiter.tree_source_iter)

def enable_ncl():
    from dendropy.dataio import nexusreader_ncl
    if nexusreader_ncl.DENDROPY_NCL_AVAILABILITY:
        _LOG.debug('Enabling Nexus Class Library bindings: using NCL NEXUS parser')
        ioclient.register("nexus", nexusreader_ncl.NexusReader, nexuswriter.NexusWriter, nexustreeiter.tree_source_iter)
    else:
        _LOG.debug('Nexus Class Library bindings are not available: using native Python NEXUS parser')

## set default ##
if "DENDROPY_ENABLE_NCL" in os.environ:
    enable_ncl()

