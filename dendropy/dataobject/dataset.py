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
Top-level phylogenetic data object: Dataset.
"""

from copy import deepcopy
from dendropy.utility import iosys
from dendropy.utility import containers
from dendropy.dataobject.base import DataObject
from dendropy.dataobject.taxon import TaxonSet
from dendropy.dataobject.tree import TreeList

###############################################################################
## Dataset

class Dataset(DataObject, iosys.Readable, iosys.Writeable):
    """
    The main data manager, consisting of a lists of taxa, trees, and character
    phylogenetic data objects, as well as methods to create, populate, access,
    and delete them.
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes a new `Dataset` object. Arguments accepted include
        another Dataset object (which will result in a deep copy of all
        taxa, trees and characters), or a file-like object opened for reading
        and a string specifying the format of the data in the file-like object,
        in which case the Dataset will be populated from data in the given
        file. Also accepts keywords `istream` and `format`.
        """
        DataObject.__init__(self)
        iosys.Writeable.__init__(self)
        iosys.Readable.__init__(self)
        self.taxon_sets = containers.OrderedSet()
        self.tree_lists = containers.OrderedSet()
        self.char_arrays = containers.OrderedSet()
        if len(args) > 0:
            if isinstance(args[0], Dataset):
                if len(args) > 1:
                    raise Exception("Only one argument accepted by Dataset() when instantiating from another Dataset")
                d = deepcopy(args[0])
                self.__dict__ = d.__dict__
            elif hasattr(args[0], "read"):
                if "istream" in kwargs:
                    raise Exception("Cannot specify both unnamed file object source ('%s') and named 'istream' source to Dataset" % (args[0]))
                istream = args[0]
                if len(args) > 2:
                    raise Exception("Dataset() accepts at most two non-keyword arguments")
                if len(args) < 2 and "format" not in kwargs:
                    raise Exception("Need to specify format if passing data source to Dataset()")
                elif len(args) == 2 and "format" in kwargs:
                    raise Exception("Cannot specify format both as non-keyworded argument ('%s') and keyworded argument" % args[1])
                elif len(args) == 2:
                    format = args[1]
                elif "format" in kwargs:
                    format = kwargs["format"]
                    del(kwargs["format"])
                self.read(istream, format, **kwargs)
            else:
                raise Exception("Invalid arguments passed to Dataset(): %s" % str(args))
        elif "istream" in kwargs:
            istream = kwargs["istream"]
            if "format" in kwargs:
                format = kwargs["format"]
                del(kwargs["istream"])
                del(kwargs["format"])
            else:
                raise Exception("Need to specify format if passing data source to Dataset()")
            self.read(istream, format, **kwargs)

    ###########################################################################
    ## CLONING

    def __deepcopy__(self, memo):
        o = self.__class__()
        memo[id(self)] = o

        for ts0 in self.taxon_sets:
            ts1 = o.new_taxon_set(label=ts0.label, oid=ts0.oid)
            memo[id(ts0)] = ts1
            for t in ts0:
                ts1.new_taxon(oid=t.oid, label=t.label)
                memo[id(t)] = ts1[-1]
        memo[id(self.taxon_sets)] = o.taxon_sets

        for tli, tl1 in enumerate(self.tree_lists):
            tl2 = o.new_tree_list(oid=tl1.oid, label=tl1.label, taxon_set=memo[id(tl1.taxon_set)])
            memo[id(tl1)] = tl2
            for ti, t1 in enumerate(tl1):
                t2 = deepcopy(t1, memo)
                t2.oid = t2.oid
                tl2.append(t2)
                memo[id(t1)] = t2

        for cai, ca1 in enumerate(self.char_arrays):
            ca2 = deepcopy(ca1, memo)
            ca2.oid = ca1.oid
            memo[id(ca1)] = ca2

        return o

    ###########################################################################
    ## I/O

    def read(self, istream, format, **kwargs):
        """
        Populates this `Dataset` object from a file-like object data
        source `istream`, formatted in `format`. `format` must be a
        recognized and supported phylogenetic data file format. If
        reading is not implemented for the format specified, then a
        `UnsupportedFormatError` is raised.

        The following optional keyword arguments are also recognized:
            - `exclude_trees` if True skips over tree data
            - `exclude_chars` if True skips over character data
            - `encode_splits` specifies whether or not split bitmasks will be
               calculated and attached to the edges.
            - `rooted` specifies the default rooting interpretation of the tree
               (see `dendropy.dataio.nexustokenizer` for details).
            - `finish_node_func` is a function that will be applied to each node
               after it has been constructed.
            - `edge_len_type` specifies the type of the edge lengths (int or float)

        Additional keyword arguments may be handled by various readers
        specialized to handle specific data formats.
        """
        from dendropy.utility import iosys
        from dendropy.dataio import get_reader
        kwargs["dataset"] = self
        reader = get_reader(format=format, **kwargs)
        return reader.read(istream, **kwargs)

    def write(self, ostream, format, **kwargs):
        """
        Writes this `Dataset` object to the file-like object `ostream`
        in `format`. `format` must be a recognized and supported
        phylogenetic data file format. If writing is not implemented for
        the format specified, then a `UnsupportedFormatError` is raised.

        The following optional keyword arguments are also recognized:
            - `exclude_trees` if True skips over tree data
            - `exclude_chars` if True skips over character data

        Additional keyword arguments may be handled by various writers
        specialized to handle specific data formats.
        """
        from dendropy.utility.iosys import require_format_from_kwargs
        from dendropy.dataio import get_writer
        kwargs["dataset"] = self
        writer = get_writer(format=format, **kwargs)
        writer.write(ostream, **kwargs)

    ###########################################################################
    ## DOMAIN DATA MANAGEMENT

    def get_taxon_set(self, **kwargs):
        """
        Returns an existing `TaxonSet` object in this `DataSet`, selected
        by keywords `oid` or `label`.
        """
        if 'oid' in kwargs:
            attr = 'oid'
            val = kwargs['oid']
        elif 'label' in kwargs:
            attr = 'label'
            val = kwargs['label']
        else:
            raise Exception("'get_taxon_set' only requires keywords 'oid' or 'label'")
        for t in self.taxon_sets:
            if getattr(t, attr) == val:
                return t
        return None

    def add_taxon_set(self, taxon_set):
        """
        Adds an existing `TaxonSet` object to this `Dataset`.
        """
        self.taxon_sets.add(taxon_set)
        return taxon_set

    def new_taxon_set(self, *args, **kwargs):
        """
        Creates a new `TaxonSet` object, according to the arguments given
        (passed to `TaxonSet()`), and adds it to this `DataSet`.
        """
        t = TaxonSet(*args, **kwargs)
        self.add_taxon_set(t)
        return t

    def add_tree_list(self, tree_list):
        "Accession of existing `TreeList` object into `tree_lists` of self."
        if tree_list.taxon_set not in self.taxon_sets:
            self.taxon_sets.add(tree_list.taxon_set)
        self.tree_lists.add(tree_list)
        return tree_list

    def new_tree_list(self, *args, **kwargs):
        "Creation and accession of new `TreeList` into `trees` of self."
        tree_list = TreeList(*args, **kwargs)
        return self.add_tree_list(tree_list)

    def add_char_array(self, char_array):
        "Accession of existing `CharacterArray` into `chars` of self."
        if char_array.taxon_set not in self.taxon_sets:
            self.taxon_sets.add(char_array.taxon_set)
        self.char_arrays.add(char_array)
        return char_array

    def new_char_array(self, char_array_type, *args, **kwargs):
        """
        Creation and accession of new `CharacterArray` (of class
        `char_array_type`) into `chars` of self."
        """
        char_array = char_array_type(*args, **kwargs)
        return self.add_char_array(char_array)
