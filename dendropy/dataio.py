#! /usr/bin/env python

############################################################################
##  dataio.py
##
##  Part of the DendroPy phylogenetic computation library.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
##  with this programm. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Data i/o.
"""
import sys

import dendropy
from dendropy import datasets
from dendropy import trees
from dendropy import taxa
from dendropy import nexml
from dendropy import nexus
from dendropy import phylip
from dendropy import fasta

def get_dataset(file=None):
    """
    General-purpose multi-format file reader.

    I've tried:
        (a) an extension-based approach
            -- failed due to inconsistency in user naming schemes
        (b) scanning to the first non-blank line and checking for signature patterns ("#NEXUS" etc.)
            -- failed due to flexibility in breaking up some patterns across multiple lines
        (c) scanning entire files for signature patterns
            -- works, but expensive due to loading entire file into buffer
        (d) brute-force: try reading files through various methods, until one works
            -- yuck, and not sure if there is any performance advantage over (c), but will do
               for now ..
    """
    try:
        if dendropy.GLOBAL_DEBUG:
            sys.stderr.write("Trying NEXML format ..\n")
        return from_nexml(fileobj=file_handle)
    except Exception, e:
        if dendropy.GLOBAL_DEBUG:
            sys.stderr.write("NEXML parse failed: %s\n" % e)
        file_handle.seek(0)
        try:
            if dendropy.GLOBAL_DEBUG:
                sys.stderr.write("Trying NEXUS format ..\n")
            return from_nexus(fileobj=file_handle)
        except Exception, e:
            if dendropy.GLOBAL_DEBUG:
                sys.stderr.write("NEXUS parse failed: %s\n" % e)
            file_handle.seek(0)
            try:
                if dendropy.GLOBAL_DEBUG:
                    sys.stderr.write("Trying NEWICK format ..\n")
                return from_newick(fileobj=file_handle)
            except Exception, e:
                if dendropy.GLOBAL_DEBUG:
                    sys.stderr.write("NEWICK parse failed: %s\n" % e)
                file_handle.seek(0)
                raise Exception("Unrecognized file format")

def from_nexml(file=None):
    """
    Reads a nexml file and returns a corresponding Dataset object.
    """
    nexml_reader = nexml.NexmlReader()
    return nexml_reader.get_dataset(file=file)

def from_nexus(file=None):
    """
    Reads a NEXUS file and returns a correspoding Dataset object.
    """
    nexus_reader = nexus.NexusReader()
    return nexus_reader.get_dataset(file=file)

# def from_newick(file=None):
#     """
#     Reads a Newick file and returns a corresponding Dataset object.
#     """
#     newick_reader = newick.NewickTreeReader()
#     file_handle = datasets.Reader.get_file_handle(file=file)
#     trees_block = newick_reader.read_trees(fileobj=file_handle, trees_block=None)
#     dataset = datasets.Dataset()
#     dataset.add_trees_block(trees_block=trees_block)
#     return dataset

def to_nexml_file(dataset, destination):
    """
    Writes out a complete nexml document representation of the dataset.
    """
    nexml_writer = nexml.NexmlWriter()
    nexml_writer.store_dataset(dataset=dataset, destination=destination)

def to_nexml_string(dataset):
    """
    Returns nexml string representation of the dataset.
    """
    nexml_writer = nexml.NexmlWriter()
    return nexml_writer.compose_dataset(dataset=dataset)

def to_nexus_file(dataset, destination, simple=False):
    """
    Writes out a complete NEXUS document representation of the dataset.
    """
    nexus_writer = nexus.NexusWriter(simple=simple)
    nexus_writer.store_dataset(dataset=dataset, destination=destination)

def to_nexus_string(dataset):
    """
    Returns NEXUS string representation of the dataset.
    """
    nexus_writer = nexus.NexusWriter()
    return nexus_writer.compose_dataset(dataset=dataset)

# def to_newick_file(dataset, destination):
#     """
#     Writes a Newick file representation of all TreesBlocks in given dataset.
#     """
#     newick_writer = newick.NewickTreeWriter()
#     newick_writer.store_dataset(dataset, destination)
# 
# def to_newick_string(dataset):
#     """
#     Returns a Newick string representation of all TreesBlocks in given dataset.
#     """
#     newick_writer = newick.NewickTreeWriter()
#     return newick_writer.compose_dataset(dataset)

def to_phylip_file(dataset, destination):
    """
    Writes a phylip representation of the first char block in given dataset.
    """
    phylip_writer = phylip.PhylipWriter()
    phylip_writer.store_dataset(dataset, destination)

def to_phylip_string(dataset):
    """
    Returns a phylip representation of the first char block in given dataset.
    """
    phylip_writer = phylip.PhylipWriter()
    return phylip_writer.compose_dataset(dataset)

def to_fasta_file(dataset, destination):
    """
    Writes a fasta representation of the first char block in given dataset.
    """
    fasta_writer = fasta.FastaWriter()
    fasta_writer.store_dataset(dataset, destination)

def to_fasta_string(dataset):
    """
    Returns a fasta representation of the first char block in given dataset.
    """
    fasta_writer = fasta.FastaWriter()
    return fasta_writer.compose_dataset(dataset)

def from_nexus_to_phylip_file(nexus, phylip):
    """
    Reads a nexus file, saves as phylip.
    """
    dataset = from_nexus(filepath=nexus)
    to_phylip_file(dataset, destination=phylip)
    