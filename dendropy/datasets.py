#! /usr/bin/env python

############################################################################
##  datasets.py
##
##  Part of the DendroPy library for phylogenetic computing.
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
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

"""
Handles the full collection of phylogenetic data: character matrices,
trees, models etc.
"""

import os
import StringIO

from dendropy import taxa
from dendropy import characters
from dendropy import trees

class Dataset(object):
    "Top-level data structure."

    def __init__(self, taxa_blocks=None, char_blocks=None, trees_blocks=None):
        "Instantiates collections of taxa, blocks, trees, and models."
        if taxa_blocks is None:
            self.taxa_blocks = []
        else:
            self.taxa_blocks = taxa_blocks
        if char_blocks is None:
            self.char_blocks = []
        else:
            self.char_blocks = char_blocks
        if trees_blocks is None:
            self.trees_blocks = []
        else:
            self.trees_blocks = trees_blocks

    def normalize_taxa_blocks(self):
        """
        Builds up list of taxon blocks by collecting taxon blocks
        referenced in self's char_blocks and trees_blocks.
        """
        self.taxa_blocks = []
        for matrix in self.char_blocks:
            self.normalize_taxa_linked(matrix)
        for trees_block in self.trees_blocks:
            self.normalize_taxa_linked(trees_block)

    def normalize_taxa_linked(self, taxa_linked):
        """
        `taxa_linked` is a trees_block or char_block or some other
        object with a `taxa_block` attribute that points to a
        TaxaBlock object. This searches the current collection of
        taxon blocks to see if the referred taxon block exists, as
        determined by the oid. If it does, then the referer's
        taxa_block is set to point to it. If not, the taxa_block is
        added to the collection.
        """
        for taxa_block in self.taxa_blocks:
            if taxa_block.oid == taxa_linked.taxa_block.oid:
                taxa_linked.taxa_block = taxa_block
                return
        self.taxa_blocks.append(taxa_linked.taxa_block)

    def find_taxa_block(self, oid=None, label=None):
        """
        Returns taxon block based on element id or label, whichever is
        given and found first.
        """
        for taxa_block in self.taxa_blocks:
            if (oid and taxa_block.oid == oid) \
               or (label and taxa_block.label == label):
                return taxa_block
        return None

    def add_taxa_block(self, taxa_block=None, oid=None, label=None, taxa_block_factory=None):
        """
        Adds (and returns) new taxa block object, creating one using
        the default factory if not given.
        """
        if taxa_block is None:
            if taxa_block_factory is None:
                taxa_block = taxa.TaxaBlock(oid=oid, label=label)
            else:
                taxa_block = taxa_block_factory(oid=oid, label=label)
        self.taxa_blocks.append(taxa_block)
        return taxa_block

    def add_taxa_linked_block(self,
                              oid=None,
                              label=None,
                              taxa_block=None,
                              linked_block=None,
                              linked_block_factory=None,
                              normalize_taxa_blocks=True):
        """
        Adds (and returns) a tree block object, creating one using the
        default factory if not given.
        """
        if linked_block is None:
            if linked_block_factory is not None:
                linked_block = linked_block_factory(oid=oid, label=label)
                if taxa_block is not None:
                    linked_block.taxa_block = taxa_block
                else:
                    self.taxa_blocks.append(linked_block.taxa_block)
            else:
                raise Exception("Neither object nor method to create object given.")
        else:
            if oid is not None:
                linked_block.oid = oid
            if label is not None:
                linked_block.label = label
        if taxa_block is not None:
            linked_block.taxa_block = taxa_block
            ## normalize_taxa_linked takes care of this ...
#             if taxa_block not in self.taxa_blocks \
#                 and taxa_block.oid not in [tb.oid for tb in self.taxa_blocks]:
#                 self.taxa_blocks.append(taxa_block)
        if normalize_taxa_blocks:
            self.normalize_taxa_linked(linked_block)
        return linked_block

    def add_trees_block(self,
                       oid=None,
                       label=None,
                       taxa_block=None,
                       trees_block=None,
                       trees_block_factory=None,
                       normalize_taxa_blocks=True):
        """
        Adds (and returns) a tree block object, creating one using the
        default factory if not given.
        """
        if trees_block is None and trees_block_factory is None:
            trees_block_factory = trees.TreesBlock
        trees_block = self.add_taxa_linked_block(oid=oid,
                                                label=label,
                                                taxa_block=taxa_block,
                                                linked_block=trees_block,
                                                linked_block_factory=trees_block_factory,
                                                normalize_taxa_blocks=normalize_taxa_blocks)
        self.trees_blocks.append(trees_block)
        return trees_block

    def add_char_block(self,
                       oid=None,
                       label=None,
                       taxa_block=None,
                       char_block=None,
                       char_block_factory=None,
                       normalize_taxa_blocks=True):
        """
        Adds (and returns) a char block object, creating one using the
        default factory if not given if not given.
        """
        if char_block is None and char_block_factory is None:
            char_block_factory = characters.CharactersBlock        
        char_block = self.add_taxa_linked_block(oid=oid,
                                                label=label,
                                                taxa_block=taxa_block,
                                                linked_block=char_block,
                                                linked_block_factory=char_block_factory,
                                                normalize_taxa_blocks=normalize_taxa_blocks)
        self.char_blocks.append(char_block)
        return char_block
        
    ###########################################################################
    ## reader/writers
    
    def read(self, src, format):
        """
        Populates this dataset from `src`, given in `format`. `src`
        is a file descriptor object, `format` is one of the supported file
        format identifiers: 'NEXUS' (incl. 'NEWICK'), 'NEXML' etc.
        """
        from dendropy import dataio
        reader = dataio.get_reader(format)
        reader.read_dataset(src, self)
        return self
        
    def from_string(self, string, format):
        """
        Populates this dataset from `string`, given in `format`. `src`
        is a file descriptor object, `format` is one of the supported file
        format identifiers: 'NEXUS' (incl. 'NEWICK'), 'NEXML' etc.
        """
        src = StringIO.StringIO(string)
        return self.read(src, format)
        
    def read_trees(self, src, format):
        """
        Populates this dataset with trees from `src`, given in `format`.
        `src` is a file descriptor object, `format` is one of the
        supported file format identifiers: 'NEXUS' (incl. 'NEWICK'),
        'NEXML' etc. A (plain) list of all trees read (including those
        from multiple TreesBlocks will be returned). In
        single-taxon-block data formats (e.g., NEXUS, NEWICK), all trees
        will share the same existing TaxonBlock (which will be expanded
        to include new taxa in the trees, if any).
        """
        from dendropy import dataio
        reader = dataio.get_reader(format)
        reader.include_characters = False
        old_trees_block_len = len(self.trees_blocks)
        reader.read_dataset(src, self)
        new_trees_block_len = len(self.trees_blocks)
        if new_trees_block_len > old_trees_block_len:
            idxs = range(old_trees_block_len, new_trees_block_len)
            new_trees = []
            for idx in idxs:
                new_trees.extend(self.trees_blocks[idx])
            return new_trees
        else:
            return []

    def iterate_over_trees(self, src, format, taxa_block=None):
        from dendropy import dataio
        reader = dataio.get_reader(format)
        reader.include_characters = False
        if taxa_block is None:
            for tree in reader.iterate_over_trees(src, dataset=self):
                yield tree
        else:
            for tree in reader.iterate_over_trees(src, taxa_block=taxa_block, dataset=self):
                yield tree
 
    def trees_from_string(self, string, format):
        """
        Populates this dataset from `string`, given in `format`. `src`
        is a file descriptor object, `format` is one of the supported file
        format identifiers: 'NEXUS' (incl. 'NEWICK'), 'NEXML' etc.
        """
        src = StringIO.StringIO(string)
        return self.read_trees(src, format)            
            
    def write(self, dest, format):
        """
        Writes dataset to `dest`, a file descriptor object, in `format`.
        """
        from dendropy import dataio
        writer = dataio.get_writer(format)
        writer.write_dataset(self, dest)              

class Reader(object):
    """
    Interface for instantiation of Dataset objects from various
    formats, to be implemented by derived classes.
    """
    def __init__(self):
        "Initializes."
        # 0 = ignore all errors; 1 = print warning; 2 = raise exception
        self.error_level=0
        self.taxa_block_factory = taxa.TaxaBlock
        self.taxon_factory = taxa.Taxon
        #self.char_block_factory = characters.CharBlock
        self.trees_block_factory = trees.TreesBlock
        self.tree_factory = trees.Tree
        self.edge_factory = trees.Edge
        self.node_factory = trees.Node
        self.include_characters = True
        self.include_trees = True
        
    def read_dataset(self, file_obj, dataset=None):
        """
        Implementing classes should instantiate and return a Dataset
        object based on contents read from the file descriptor object
        `file_obj`. If `dataset` is given, then this will be populated,
        otherwise a new one will be created. If `taxa_block` is given, 
        then this will be used (and added to `dataset`.
        """
        raise NotImplementedError

    def read_characters(self, file_obj):
        """
        Instantiates and returns a list of CharacterBlock objects from a 
        file (descriptor).
        """
        dataset = Dataset()        
        dataset = self.read_dataset(file_obj=file_obj)
        return dataset.char_blocks

    def read_taxa(self, file_obj):
        """
        Instantiates and returns a list of TaxaBlock objects from a 
        file (descriptor).
        """
        dataset = Dataset()        
        dataset = self.read_dataset(file_obj=file_obj)
        return dataset.taxa_blocks

    def read_trees(self, file_obj):
        """
        Instantiates and returns a list of TreeBlock objects from a 
        file (descriptor).
        """
        dataset = self.read_dataset(file_obj=file_obj)
        return dataset.trees_blocks

class Writer(object):
    """
    Interface for composing and writing a representation of a DataSet
    object in various formats, to be implemented by derived classes.
    """

    def write_dataset(self, dataset, dest):
        """
        Writes a Dataset object to a full document-level
        representation of the format being implemented by the deriving
        class. `dest` is an output stream that support 'write'.
        """
        raise NotImplementedError

    def compose_dataset(self, dataset):
        """
        Returns a string representation of a DataSet as a fully-formed
        and formatted dataset document.
        """
        dataset_text = StringIO.StringIO()
        self.write_dataset(dataset, dataset_text)
        return dataset_text.getvalue()

