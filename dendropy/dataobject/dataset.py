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
Top-level phylogenetic data object: DataSet.
"""

from copy import deepcopy
from cStringIO import StringIO
from dendropy.utility import iosys
from dendropy.utility import containers
from dendropy.utility import error
from dendropy.dataobject.base import DataObject
from dendropy.dataobject.taxon import TaxonSet
from dendropy.dataobject.tree import TreeList
from dendropy.dataobject.char import CharacterMatrix

###############################################################################
## DataSet

class DataSet(DataObject, iosys.Readable, iosys.Writeable):
    """
    The main data manager, consisting of a lists of taxa, trees, and character
    phylogenetic data objects, as well as methods to create, populate, access,
    and delete them.
    """

    def __init__(self, *args, **kwargs):
        """
        __init__ takes a new `DataSet` object from another DataSet object or by
        parsing an input stream.

        Can be invoked with:

            - a single unnamed argument that is an instance of a DataSet (which will result in a deep copy of all taxa, trees and characters),
            - 1 or more unnamed TaxonSet, TreeList or CharacterMatrix (which is like calling DataSet.add with each argument)
            - keyword arguments that supply a file-like object (`stream`) to be parsed and a string (`schema`) specifying the schema of the data in the file-like object.

        """
        DataObject.__init__(self)
        iosys.Writeable.__init__(self)
        iosys.Readable.__init__(self)
        self.taxon_sets = containers.OrderedSet()
        self.tree_lists = containers.OrderedSet()
        self.char_matrices = containers.OrderedSet()
        self.attached_taxon_set = None
        taxa = kwargs.get("taxon_set", None)
        attach_taxon_set = kwargs.get("attach_taxon_set", False)
        if attach_taxon_set or (taxa is not None):
            self.attach_taxon_set(taxa)
        else:
            self.attached_taxon_set = None
        stream = kwargs.get("stream")
        if len(args) > 0:
            if (stream is not None) or (kwargs.get("schema") is not None):
                raise error.MultipleInitializationSourceError(self.__class__.__name__, args[0])
            if len(args) == 1 and isinstance(args[0], DataSet):
                if attach_taxon_set or (taxa is not None):
                    raise error.MultipleInitializationSourceError("Cannot initialize DataSet from another DataSetobject taxon_set or attach_taxon_set are specified")
                d = deepcopy(args[0])
                self.__dict__ = d.__dict__
            else:
                for arg in args:
                    if isinstance(arg, DataSet):
                        raise error.MultipleInitializationSourceError("Cannot initialize DataSet from another DataSetobject when multiple other initialization objects are given")
                    else:
                        self.add(arg)
        elif stream is not None:
            if self.attached_taxon_set is not None:
                kwargs["taxon_set"] = self.attached_taxon_set
            self.process_source_kwargs(**kwargs)

    ###########################################################################
    ## CLONING

    def __deepcopy__(self, memo):
        o = self.__class__()
        for ts0 in self.taxon_sets:
            ts1 = o.new_taxon_set(label=ts0.label)
            memo[id(ts0)] = ts1
            for t in ts0:
                ts1.new_taxon(label=t.label)
                memo[id(t)] = ts1[-1]
        memo[id(self.taxon_sets)] = o.taxon_sets
        for tli, tl1 in enumerate(self.tree_lists):
            tl2 = o.new_tree_list(label=tl1.label, taxon_set=memo[id(tl1.taxon_set)])
            memo[id(tl1)] = tl2
            for ti, t1 in enumerate(tl1):
                t2 = deepcopy(t1, memo)
                tl2.append(t2)
                memo[id(t1)] = t2
        memo[id(self.tree_lists)] = o.tree_lists
        for cai, ca1 in enumerate(self.char_matrices):
            ca2 = deepcopy(ca1, memo)
            o.char_matrices.add(ca2)
            memo[id(ca1)] = ca2
        memo[id(self.char_matrices)] = o.char_matrices
        if self.attached_taxon_set is not None:
            o.attached_taxon_set = memo[id(self.attached_taxon_set)]
        else:
            o.attached_taxon_set = None
        return o

    ###########################################################################
    ## I/O and Representation

    def read(self, stream, schema, **kwargs):
        """
        Populates this `DataSet` object from a file-like object data
        source `stream`, formatted in `schema`. `schema` must be a
        recognized and supported phylogenetic data file schema. If
        reading is not implemented for the schema specified, then a
        `UnsupportedSchemaError` is raised.

        The following optional keyword arguments are also recognized:

            - `exclude_trees` if True skips over tree data
            - `exclude_chars` if True skips over character data
            - `encode_splits` specifies whether or not split bitmasks will be
               calculated and attached to the edges.
            - `finish_node_func` is a function that will be applied to each node
               after it has been constructed.

        The following keyword arguments are recognized when parsing NEXUS or
        NEWICK sources:

            - `taxon_set` TaxonSet object to use when reading data
            - `as_rooted=True` (or `as_unrooted=False`) interprets trees as rooted
            - `as_unrooted=True` (or `as_rooted=False`) interprets trees as rooted
            - `default_as_rooted` (or `default_as_unrooted=False`) interprets
               all trees as rooted if rooting not given by `[&R]` or `[&U]` comments
            - `default_as_unrooted` (or `default_as_rooted=False`) interprets
               all trees as rooted if rooting not given by `[&R]` or `[&U]` comments
            - `edge_len_type` specifies the type of the edge lengths (int or float)

        Additional keyword arguments may be handled by various readers
        specialized to handle specific data formats.
        """
        from dendropy.utility import iosys
        from dendropy.dataio import get_reader
        kwargs["dataset"] = self
        if self.attached_taxon_set is not None:
            if "taxon_set" not in kwargs:
                kwargs["taxon_set"] = self.attached_taxon_set
            elif kwargs["taxon_set"] is not self.attached_taxon_set:
                raise TypeError("DataSet object is already attached to a TaxonSet, but different TaxonSet passed to using 'taxon_set' keyword argument")
        reader = get_reader(schema=schema, **kwargs)
        try:
            reader.read(stream)
#        except error.DataParseError as x:
        except error.DataParseError, x:
            x.decorate_with_name(stream=stream)
            raise x

    def write(self, stream, schema, **kwargs):
        """
        Writes this `DataSet` object to the file-like object `stream`
        in `schema`. `schema` must be a recognized and supported
        phylogenetic data file schema. If writing is not implemented for
        the schema specified, then a `UnsupportedSchemaError` is raised.

        The following optional keyword arguments are also recognized:

            - `exclude_trees` if True skips over tree data
            - `exclude_chars` if True skips over character data

        Additional keyword arguments may be handled by various writers
        specialized to handle specific data formats.
        """
        from dendropy.utility.iosys import require_format_from_kwargs
        from dendropy.dataio import get_writer
        kwargs["dataset"] = self
#        if self.attached_taxon_set is not None:
#            if "taxon_set" not in kwargs:
#                kwargs["taxon_set"] = self.attached_taxon_set
#            elif kwargs["taxon_set"] is not self.attached_taxon_set:
#                raise TypeError("DataSet object is already attached to a TaxonSet, but different TaxonSet passed using 'taxon_set' keyword argument")
        writer = get_writer(schema=schema, **kwargs)
        writer.write(stream)

    def _subdescription(self, name, objs, depth, indent, itemize, output, **kwargs):
        if len(objs) == 0:
            return
        output.write('\n%s[%s]' % (indent*' ', name))
        for i, obj in enumerate(objs):
            output.write('\n')
            obj.description(depth=depth-1,
                       indent=indent+4,
                       itemize="[%d] " % (i),
                       output=output,
                       **kwargs)

    def description(self, depth=1, indent=0, itemize="", output=None):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return ""
        output_strio = StringIO()
        output_strio.write('DataSet object at %s' % hex(id(self)))
        if depth >= 1:
            output_strio.write(': %d Taxon Sets, %d Tree Lists, %d Character Matrices' %
                    (len(self.taxon_sets), len(self.tree_lists), len(self.char_matrices)))
        if depth >= 2:
            indent += 4
            self._subdescription('Taxon Sets', self.taxon_sets, depth, indent, itemize, output_strio)
            self._subdescription('Tree Lists', self.tree_lists, depth, indent, itemize, output_strio)
            self._subdescription('Character Matrices', self.char_matrices, depth, indent, itemize, output_strio)
        s =  output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

    ###########################################################################
    ## DOMAIN DATA MANAGEMENT

    def add(self, data_object, **kwargs):
        """
        Generic add for TaxonSet, TreeList or CharacterMatrix objects.
        """
        if isinstance(data_object, TaxonSet):
            self.add_taxon_set(data_object)
        elif isinstance(data_object, TreeList):
            self.add_tree_list(data_object)
        elif isinstance(data_object, CharacterMatrix):
            self.add_char_matrix(data_object)
        else:
            raise error.InvalidArgumentValueError("Cannot add object of type '%s' to DataSet" % type(data_object))

    def get_default_taxon_set(self, **kwargs):
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
            raise Exception("'get_default_taxon_set' only requires keywords 'oid' or 'label'")
        for t in self.taxon_sets:
            if getattr(t, attr) == val:
                return t
        return None

    def add_taxon_set(self, taxon_set):
        """
        Adds an existing `TaxonSet` object to this `DataSet`.
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

    def attach_taxon_set(self, taxon_set=None):
        """
        Forces all read() calls on this DataSet to use the same TaxonSet. If
        `taxon_set` If `taxon_set` is None, then a new TaxonSet will be
        created, added to self.taxa, and that is the TaxonSet that will be
        attached.
        """
        if taxon_set is None:
            taxon_set = self.new_taxon_set()
        elif taxon_set not in self.taxon_sets:
            self.add_taxon_set(taxon_set)
        self.attached_taxon_set = taxon_set
        return self.attached_taxon_set

    def detach_taxon_set(self):
        self.attached_taxon_set = None

    def unify_taxa(self, taxon_set=None, bind=True):
        """
        Reindices taxa across all subcomponents, mapping to single taxon set.
        """
        if len(self.taxon_sets) or len(self.tree_lists) or len(self.char_matrices):
            self.taxon_sets.clear()
            if taxon_set is None:
                taxon_set = self.new_taxon_set()
            for tree_list in self.tree_lists:
                tree_list.reindex_taxa(taxon_set=self.taxon_set, clear=False)
            for char_matrix in self.char_matrices:
                char_matrix.reindex_taxa(taxon_set=self.taxon_set, clear=False)
        if bind:
            self.attach_taxon_set(taxon_set)

    def add_tree_list(self, tree_list):
        "Accession of existing `TreeList` object into `tree_lists` of self."
        if self.attached_taxon_set is not None:
            tree_list.reindex_taxa(taxon_set=self.attached_taxon_set, clear=False)
        elif tree_list.taxon_set not in self.taxon_sets:
            self.taxon_sets.add(tree_list.taxon_set)
        self.tree_lists.add(tree_list)
        return tree_list

    def new_tree_list(self, *args, **kwargs):
        "Creation and accession of new `TreeList` into `trees` of self."
        if self.attached_taxon_set is not None:
            if "taxon_set" in kwargs and kwargs["taxon_set"] is not self.attached_taxon_set:
                raise TypeError("DataSet object is already attached to a TaxonSet, but different 'taxon_set' passed as argument")
            else:
                kwargs["taxon_set"] = self.attached_taxon_set
        tree_list = TreeList(*args, **kwargs)
        return self.add_tree_list(tree_list)

    def add_char_matrix(self, char_matrix):
        "Accession of existing `CharacterMatrix` into `chars` of self."
        if self.attached_taxon_set is not None:
            char_matrix.reindex_taxa(taxon_set=self.attached_taxon_set, clear=False)
        elif char_matrix.taxon_set not in self.taxon_sets:
            self.taxon_sets.add(char_matrix.taxon_set)
        self.char_matrices.add(char_matrix)
        return char_matrix

    def new_char_matrix(self, char_matrix_type, *args, **kwargs):
        """
        Creation and accession of new `CharacterMatrix` (of class
        `char_matrix_type`) into `chars` of self."
        """
        if self.attached_taxon_set is not None:
            if "taxon_set" in kwargs and kwargs["taxon_set"] is not self.attached_taxon_set:
                raise TypeError("DataSet object is already attached to a TaxonSet, but different 'taxon_set' passed as argument")
            else:
                kwargs["taxon_set"] = self.attached_taxon_set
        char_matrix = char_matrix_type(*args, **kwargs)
        return self.add_char_matrix(char_matrix)

    def get_tree_list(self, **kwargs):
        """
        Returns a TreeList object specified by one (and exactly one) of the
        following keyword arguments:

            - ``label``
            - ``oid``

        Raises ``KeyError`` if no matching ``TreeList`` is found, unless
        ``ignore_error`` is set to ``True``.
        """
        if "label" in kwargs and "oid" in kwargs:
            raise TypeError("Cannot specify both 'label' and 'oid'")
        elif "label" in kwargs:
            for t in self.tree_lists:
                if t.label == kwargs['label']:
                    return t
            if not kwargs.get("ignore_error", False):
                raise KeyError(kwargs['label'])
        elif "oid" in kwargs:
            for t in self.tree_lists:
                if t.oid == kwargs['oid']:
                    return t
            if not kwargs.get("ignore_error", False):
                raise KeyError(kwargs['oid'])
        else:
            raise TypeError("Must specify one of: 'label' or 'oid'")

