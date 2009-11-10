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

from dendropy.utility import error
#from dendropy.dataobject.taxon import TaxonSet
from dendropy.dataio import dataformat
from dendropy.dataio import newick
from dendropy.dataio import nexus
from dendropy.dataio import fasta
from dendropy.dataio import phylip
from dendropy.dataio import nexml

# syntax is:
#   dataformat.register(<FORMAT NAME>, <READER TYPE>, <WRITER TYPE>, <TREE ITERATOR>, <TREE (LIST) WRITER>)
dataformat.register("newick", newick.NewickReader, newick.NewickWriter, newick.tree_source_iter, newick.write_tree_list)
dataformat.register("nexus", nexus.NexusReader, nexus.NexusWriter, nexus.tree_source_iter, nexus.write_tree_list)
dataformat.register("nexus/newick", None, None, nexus.generalized_tree_source_iter, None)
dataformat.register("fasta", None, fasta.FastaWriter, None, None)
dataformat.register("dnafasta", fasta.DNAFastaReader, fasta.FastaWriter, None, None)
dataformat.register("rnafasta", fasta.RNAFastaReader, fasta.FastaWriter, None, None)
dataformat.register("proteinfasta", fasta.ProteinFastaReader, fasta.FastaWriter, None, None)
dataformat.register("phylip", None, phylip.PhylipWriter, None, None)
dataformat.register("nexml", nexml.NexmlReader, nexml.NexmlWriter, None, None)

from dendropy.dataio.dataformat import get_reader, get_writer, tree_source_iter, write_tree_list

def multi_tree_source_iter(sources, format, **kwargs):
    """
    Diagnoses and handles multiple files in either NEXUS or NEWICK formats, as
    given by `sources` which should be a list of file-like objects or strings
    specifying file paths. Not that unless a TaxonSet object is explicitly
    passed using the 'taxon_set' keyword argument, the trees in each file will
    be associated with their own distinct, independent taxon sets: almost
    certainly not the desired behavior.
    """
#    if "taxon_set" not in kwargs:
#        kwargs["taxon_set"] = TaxonSet()
    for s in sources:
        if isinstance(s, str):
            src = open(s, "rU")
        else:
            src = s
        for t in tree_source_iter(src, format, **kwargs):
            yield t
