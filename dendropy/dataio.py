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
from dendropy import newick
from dendropy import phylip
from dendropy import fasta

def get_dataset(filepath=None, fileobj=None, text=None):
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
               for now ...
    """
    file_handle = datasets.Reader.get_file_handle(filepath=filepath, fileobj=fileobj, text=text)
    try:
        if dendropy.GLOBAL_DEBUG:
            sys.stderr.write("Trying NEXML format ...\n")
        return from_nexml(fileobj=file_handle)        
    except Exception, e:
        if dendropy.GLOBAL_DEBUG:
            sys.stderr.write("NEXML parse failed: %s\n" % e)    
        file_handle.seek(0)
        try:
            if dendropy.GLOBAL_DEBUG:
                sys.stderr.write("Trying NEXUS format ...\n")        
            return from_nexus(fileobj=file_handle)
        except Exception, e:
            if dendropy.GLOBAL_DEBUG:
                sys.stderr.write("NEXUS parse failed: %s\n" % e)
            file_handle.seek(0)
            try:
                if dendropy.GLOBAL_DEBUG:
                    sys.stderr.write("Trying NEWICK format ...\n")
                return from_newick(fileobj=file_handle)
            except Exception, e:
                if dendropy.GLOBAL_DEBUG:
                    sys.stderr.write("NEWICK parse failed: %s\n" % e)
                file_handle.seek(0)
                raise Exception("Unrecognized file format")

def iterate_over_trees(filepath=None, fileobj=None, text=None):
    """
    Generator to iterate over trees in data file.
    Primary goal is to be memory efficient, storing no more than one tree 
    at a time. Speed might have to be sacrificed for this!
    """
    
    taxa_block = taxa.TaxaBlock()        
    nexus_reader = nexus.NexusReader()
    filehandle = datasets.Reader.get_file_handle(filepath=filepath, fileobj=fileobj, text=text)
    nexus_reader.filehandle = filehandle
    token = nexus_reader.read_next_token_ucase()
    
    if token == "#NEXUS":
        file_format = "NEXUS"
        while not nexus_reader.eof:        
            token = nexus_reader.read_next_token_ucase()
            while token != None and token != 'BEGIN' and not nexus_reader.eof:
                token = nexus_reader.read_next_token_ucase()
            token = nexus_reader.read_next_token_ucase()
            if token == 'TREES':
                trees_block = nexus_reader.trees_block_factory()
                trees_block.taxa_block = nexus_reader.get_default_taxa_block()
                nexus_reader.dataset.add_trees_block(trees_block=trees_block)
                nexus_reader.skip_to_semicolon() # move past BEGIN command
                while not (token == 'END' or token == 'ENDBLOCK') and not nexus_reader.eof and not token==None:
                    token = nexus_reader.read_next_token_ucase()
                    if token == 'TRANSLATE':
                        nexus_reader.parse_translate_statement()                         
                    if token == 'TREE':
                        tree = nexus_reader.parse_tree_statement(taxa_block)  
                        trees_block.pop()
                        yield tree
                nexus_reader.skip_to_semicolon() # move past END command    
            else:
                # unknown block
                while not (token == 'END' or token == 'ENDBLOCK') and not nexus_reader.eof and not token==None:
                    #print token
                    nexus_reader.skip_to_semicolon()
                    token = nexus_reader.read_next_token_ucase()         
        
    else:
        ### if not NEXUS, assume NEWICK ###      
        
        filehandle.seek(0)
        file_format = "NEWICK"

        # Load entire file at once and then parse ...
#         statement_block = filehandle.read()
#         statement_block = statement_block.replace('\n','').replace('\r','')
#         for statement in statement_block.split(';'):
#             # -- parse statement and yield tree --

        # Read stream byte-by-byte and parse as we go ...
        while True:
            statement = []
            ch = filehandle.read(1)
            while ch != '' and ch != ';':
                if ch not in ['\n', '\r']:
                    statement.append(ch)
                ch = filehandle.read(1)            
            if statement:                
                statement = ''.join(statement).replace('\n','').replace('\r','').strip()            
                newick_parser = newick.NewickTreeParser()
                tree = newick_parser.parse_tree_statement(statement, taxa_block)
                yield tree
            if ch == '':
                break
            
def tree_iter(filepath=None, fileobj=None, text=None):
    return nexus_tree_iter(filepath=filepath, fileobj=fileobj, text=text)
                
def nexus_tree_iter(filepath=None, fileobj=None, text=None):
    nexus_reader = nexus.NexusReader()
    return nexus_reader.tree_iter(filepath=filepath, fileobj=fileobj, text=text)
    
def newick_tree_iter(filepath=None, fileobj=None, text=None):
    newick_reader = newick.NewickTreeReader()
    return newick_reader.tree_iter(filepath=filepath, fileobj=fileobj, text=text)    
        
def from_nexml(filepath=None, fileobj=None, text=None):
    """
    Reads a nexml file and returns a corresponding Dataset object.
    """
    nexml_reader = nexml.NexmlReader()
    return nexml_reader.get_dataset(filepath=filepath, fileobj=fileobj, text=text)

def from_nexus(filepath=None, fileobj=None, text=None):
    """
    Reads a NEXUS file and returns a correspoding Dataset object.
    """
    nexus_reader = nexus.NexusReader()
    return nexus_reader.get_dataset(filepath=filepath, fileobj=fileobj, text=text)
        
def from_newick(filepath=None, fileobj=None, text=None):
    """
    Reads a Newick file and returns a corresponding Dataset object.
    """
    newick_reader = newick.NewickTreeReader()
    file_handle = datasets.Reader.get_file_handle(filepath=filepath, fileobj=fileobj, text=text)
    trees_block = newick_reader.read_trees(fileobj=file_handle, trees_block=None)
    dataset = datasets.Dataset()    
    dataset.add_trees_block(trees_block=trees_block)
    return dataset

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
    
def to_newick_file(dataset, destination):
    """
    Writes a Newick file representation of all TreesBlocks in given dataset.
    """
    newick_writer = newick.NewickTreeWriter()
    newick_writer.store_dataset(dataset, destination)

def to_newick_string(dataset):
    """
    Returns a Newick string representation of all TreesBlocks in given dataset.
    """
    newick_writer = newick.NewickTreeWriter()
    return newick_writer.compose_dataset(dataset)

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
    