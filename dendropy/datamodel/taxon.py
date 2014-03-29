#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Taxon management.

Operation taxonomic unit concepts are essentially names for taxa in the "real
world". Operational taxonomic unit concepts are organized into taxonomic
namespaces. A taxonomic namespace is a self-contained and
functionally-complete collection of mutually-distinct operational taxonomic
unit concepts, and provide the semantic context in which operational taxonomic
units from across various data sources of different formats and provenances can
be related through correct interpretation of their taxon labels.

    * Operational taxonomic units are modeled by a :class:`Taxon` object.

    * Taxonomic namespaces, in which operational taxonomic units are organized,
      are modeled by a :class:`TaxonNamespace` object.

    * A :class:`TaxonNamespace` manages a collection of :class:`Taxon` objects, where each
      object represents a distinct operational taxonomic unit concept within
      the taxonomic namespace represented by that :class:`TaxonNamespace` object.

    * Each :class:`Taxon` object can belong to one and only one :class:`TaxonNamespace`:
      :class:`Taxon` objects are not shared across :class:`TaxonNamespace` objects.

    * Each :class:`Taxon` object has an attribute, `label`, whose (string) value
      is the name of the operational taxon unit concept that it represents.

    * Different :class:`Taxon` objects represent different operational taxonomic
      unit concepts, even if they have the same label value.

    * All client objects (:class:`TaxonNamespaceAssociated` objects) that reference
      the same :class:`TaxonNamespace` reference the same "universe" or domain of
      operational taxonomic unit concepts.

    * Operational taxonomic units from across different data sources are mapped
      to distinct :class:`Taxon` objects within a particular :class:`TaxonNamespace` based on
      matching the string values of labels of the :class:`Taxon` object.

    * A particular taxonomic unit concept in one data source will only be
      correctly related to the same taxonomic unit concept (i.e, the same
      :class:`Taxon` object) in another data source only if they have both
      been parsed with reference to the same taxonomic namespace (i.e., the
      same :class:`TaxonNamespace` has been used).

    * A :class:`TaxonNamespace` assigned an "accession index" to every :class:`Taxon` object
      added to it. This is a stable and unique number within the context of any
      given :class:`TaxonNamespace` object (though a :class:`Taxon` object may have different
      accession indexes in different :class:`TaxonNamespace` objects if it
      belongs to multiple namespaces). This number is will be used to
      calculate the "split bitmask" hash of the trivial split or external edge
      subtending the node to which this :class:`Taxon` object is assigned on a tree.
      The concept of a "split bitmask" hash is fundamental to DendroPy's tree
      operations. The split bitmask is a hash that uniquely identifies every
      split on a tree.  It is calculated by OR'ing the split bitmask of all the
      child splits of the given split. Terminal edges, of course, do not have
      child edges, and their split bitmask is given by the accession index of
      the :class:`Taxon` object at their head or target nodes.
"""


import warnings
import collections
import copy
try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3
from dendropy.datamodel import base
from dendropy.utility import text
from dendropy.utility import container
from dendropy.utility import error

##############################################################################
## Helper functions

def taxon_set_deprecation_warning():
    error.dump_stack()
    warnings.warn("`taxon_set` will no longer be supported in future releases; use `taxon_namespace` instead",
            FutureWarning, stacklevel=4)

def process_kwargs_for_taxon_namespace(kwargs_dict, default=None):
    if "taxon_set" in kwargs_dict:
        if "taxon_namespace" in kwargs_dict:
            raise TypeError("Cannot specify both 'taxon_namespace' and 'taxon_set' (legacy support) simultaneously")
        else:
            taxon_set_deprecation_warning()
            return kwargs_dict.pop("taxon_set", default)
    else:
        return kwargs_dict.pop("taxon_namespace", default)

##############################################################################
## TaxonAssociated

class TaxonAssociated(base.DataObject):

    def __init__(self, **kwargs):
        base.DataObject.__init__(self, label=kwargs.pop('label', None))
        if "taxon_namespace" not in kwargs or kwargs["taxon_namespace"] is None:
            self.taxon_namespace = TaxonNamespace()
        else:
            self.taxon_namespace = kwargs["taxon_namespace"]

##############################################################################
## TaxonNamespaceAssociated

class TaxonNamespaceAssociated(base.DataObject):
    """
    Provides infrastructure for the maintenance of references to taxa.
    """

    def __init__(self, **kwargs):
        base.DataObject.__init__(self, label=kwargs.pop('label', None))
        if "taxon_namespace" not in kwargs or kwargs["taxon_namespace"] is None:
            self.taxon_namespace = TaxonNamespace()
        else:
            self.taxon_namespace = kwargs["taxon_namespace"]

    ## for legacy
    def _get_taxon_set(self):
        # raise NotImplementedError("'taxon_set' is no longer supported: use 'taxon_namespace' instead")
        taxon_set_deprecation_warning()
        return self.taxon_namespace
    def _set_taxon_set(self, v):
        # raise NotImplementedError("'taxon_set' is no longer supported: use 'taxon_namespace' instead")
        taxon_set_deprecation_warning()
        self.taxon_namespace = v
    def _del_taxon_set(self):
        # raise NotImplementedError("'taxon_set' is no longer supported: use 'taxon_namespace' instead")
        taxon_set_deprecation_warning()
    taxon_set = property(_get_taxon_set, _set_taxon_set, _del_taxon_set)

    def reindex_taxa(self, taxon_namespace=None, clear=False):
        """
        Rebuilds `taxon_namespace` from scratch, or assigns :class:`Taxon` objects from
        given :class:`TaxonNamespace` object `taxon_namespace` based on label values. Calls
        on `self.reindex_member_taxa()` to synchronize taxa.
        """
        if taxon_namespace is not None:
            self.taxon_namespace = taxon_namespace
        if clear:
            self.taxon_namespace.clear()
        self.reindex_subcomponent_taxa()
        return self.taxon_namespace

    def reindex_subcomponent_taxa():
        """
        Derived classes should override this to ensure that their various
        components, attributes and members all refer to the same :class:`TaxonNamespace`
        object as `self.taxon_namespace`, and that `self.taxon_namespace` has all
        the :class:`Taxon` objects in the various members.
        """
        pass

##############################################################################
## TaxonNamespace

class TaxonNamespace(base.DataObject, base.Annotable):
    """
    A collection of :class:`Taxon` objects representing a self-contained and complete
    domain of distinct operational taxonomic unit definitions.
    Provides the common semantic context in which operational taxonomic units
    referenced by various phylogenetic data objects (e.g., trees or alignments)
    can be related.
    """

    ###########################################################################
    ## Life-cycle

    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------

        \*args : positional arguments, optional
            Accepts a single iterable as an optional positional argument.  If a
            :class:`TaxonNamespace` object is passed as the positional argument, then
            clones or deep-copies of its member :class:`Taxon` objects will be added
            to this one.  If any other iterable is passed as the positional
            argument, then each string in the iterable will result in a new
            :class:`Taxon` object being constructed and added to the namespace with
            the string as its label (name), while each Taxon object in the
            iterable will be added to the namespace directly.

        \*\*kwargs : keyword arguments
            label : string
                The label or name for this namespace.
            is_mutable : boolean, optional (default = `True`)
                If `True` (default), then :class:`Taxon` objects can be added to this
                namespace. If `False`, then adding :class:`Taxon` objects will result
                in an error.

        Notes
        -----
        An empty :class:`TaxonNamespace` can be created (with optional) label and :class:`Taxon`
        objects added later:

        >>> tns = dendropy.TaxonNamespace(label="taxa")
        >>> t1 = Taxon("a")
        >>> tns.add_taxon(t1)
        >>> t2 = Taxon("b")
        >>> tns.add_taxon(t2)
        >>> tns.add_taxon("c")
        >>> tns
        <TaxonNamespace 0x106509090 'taxa': [<Taxon 0x10661f050 'a'>, <Taxon 0x10651c590 'b'>, <Taxon 0x106642a90 'c'>]>

        Alternatively, an iterable can be passed in as an initializer, and all
        :class:`Taxon` objects will be added directly while, for each string, a new
        :class:`Taxon` object will be created and added. So, the below are all equivalent
        to the above:

        >>> tns = dendropy.TaxonNamespace(["a", "b", "c"], label="taxa")

        >>> taxa = [Taxon(n) for n in ["a", "b", "c"]]
        >>> tns = dendropy.taxonnamespace(taxa, label="taxa")

        >>> t1 = Taxon("a")
        >>> t2 = Taxon("b")
        >>> taxa = [t1, t2, "c"]
        >>> tns = dendropy.TaxonNamespace(taxa, label="taxa")

        If a :class:`TaxonNamespace` object is passed as the
        initializer argument, a *shallow* copy of the object is constructed:

        >>> tns1 = dendropy.TaxonNamespace(["a", "b", "c"], label="taxa1")
        >>> tns1
        <TaxonNamespace 0x1097275d0 'taxa1': [<Taxon 0x109727610 'a'>, <Taxon 0x109727e10 'b'>, <Taxon 0x109727e90 'c'>]>
        >>> tns2 = dendropy.TaxonNamespace(tns1, label="2")
        >>> tns2
        <TaxonNamespace 0x109727d50 'taxa1': [<Taxon 0x109727610 'a'>, <Taxon 0x109727e10 'b'>, <Taxon 0x109727e90 'c'>]>

        Thus, while "`tns1`" and "`tns2`" are independent collections, and
        addition/deletion of :class:`Taxon` instances to one will not effect
        the other, the label of a :class:`Taxon` instance that is an element in
        one will of course effect the same instance if it is in the other:

        >>> print(tns1[0].label)
        >>> a
        >>> print(tns2[0].label)
        >>> a
        >>> tns1[0].label = "Z"
        >>> print(tns1[0].label)
        >>> Z
        >>> print(tns2[0].label)
        >>> Z

        If what is needed is a **full** or **deep-copy** of a
        :class:`TaxonNamespace`, including copies of the member :class:`Taxon`
        instances, then  this can be achieved using :func:`copy.deepcopy()` or
        :meth:`TaxonNamespace.clone()`:

        >>> tns1 = dendropy.TaxonNamespace(["a", "b", "c"], label="taxa1")
        >>> tns2 = tns1.clone()
        >>> tns3 = copy.deepcopy(tns2)
        """
        base.DataObject.__init__(self, label=kwargs.get('label'))
        base.Annotable.__init__(self)
        self.is_mutable = kwargs.get('is_mutable', True)
        self._accession_index_taxon_map = {}
        self._taxa = []
        self._taxon_accession_index_map = {}
        self._taxon_bitmask_map = {}
        self._bitmask_taxon_map = {}
        self._current_accession_count = 0
        if len(args) > 1:
            raise TypeError("TaxonNamespace() takes at most 1 non-keyword argument ({} given)".format(len(args)))
        elif len(args) == 1:
            for i in args[0]:
                if isinstance(i, Taxon):
                    self.add_taxon(i)
                else:
                    self.new_taxon(label=i)

    ###########################################################################
    ## Identity and Comparison

    def __str__(self):
        return "[{}]".format(", ".join([str(i) for i in self._taxa]))

    def __repr__(self):
        return "<{} {} '{}': [{}]>".format(self.__class__.__name__, hex(id(self)), self.label, ", ".join(repr(i) for i in self._taxa))

    def __hash__(self):
        return id(self)

    def __lt__(self, other):
        return self._taxa < o._taxa

    def __eq__(self, other):
        # enforce non-equivalence of non-identical namespaces
        return self is other
        # if not isinstance(other, self.__class__):
        #     return False
        # return (self.label == other.label
        #         and self._taxa == other._taxa
        #         and base.Annotable.__eq__(self, other))

    ###########################################################################
    ## Collection Iteration

    def __iter__(self):
        return iter(self._taxa)

    def __reversed__(self):
        return reversed(self._taxa)

    ###########################################################################
    ## Collection Data

    def __len__(self):
        """
        Returns number of :class:`Taxon` objects in this :class:`TaxonNamespace`.
        """
        return len(self._taxa)

    ###########################################################################
    ## Collection Access and Management

    def __getitem__(self, key):
        """
        Returns :class:`Taxon` object with index or slice given by `key`.
        """
        if isinstance(key, int) or isinstance(key, slice):
            return self._taxa[key]
        raise ValueError("'TaxonNamespace[]' now only accepts indexes or slices. To access Taxon objects by label, use 'TaxonNamespace.get_taxon()' or 'TaxonNamespace.findall()'")

    def __setitem__(self, key, value):
        raise NotImplementedError("Item assignment not supported")

    def __delitem__(self, key):
        self.remove_taxon(self[key])

    def __contains__(self, taxon):
        """
        Returns `True` if Taxon object `taxon` is in self.
        """
        # look-up in dictionary for O(1) instead of O(n) in list
        return taxon in self._taxon_accession_index_map

    ###########################################################################
    ## Adding Taxa

    def add_taxon(self, taxon):
        """
        Adds a new :class:`Taxon` object to `self`.

        If `taxon` is not already in the collection of :class:`Taxon` objects in this
        namespace, and this namespace is mutable, it is added to the
        collection. If it is already in the collection, then nothing happens.
        If it is not already in the collection, but the namespace is not
        mutable, then `TypeError` is raised.

        Parameters
        ----------
        taxon : :class:`Taxon`
            The :class:`Taxon` object to be accessioned or registered in this
            collection.

        Raises
        ------
        TypeError
            If this namespace is immutable (i.e.
            :attr:`TaxonNamespace.is_mutable` is `False`).

        """
        ### NOTE
        ### Previously, this was:
        ###
        ###     if taxon in self._taxa:
        ###
        ### Changing the membership lookup to dictionaries resulted in 10x
        ### increase in speed!!!!
        if taxon in self._taxon_accession_index_map:
            return
        if not self.is_mutable:
            raise TypeError("Taxon '{}' cannot be added to an immutable TaxonNamespace".format((taxon.label)))
        self._current_accession_count += 1
        self._taxa.append(taxon)
        self._accession_index_taxon_map[self._current_accession_count] = taxon
        self._taxon_accession_index_map[taxon] = self._current_accession_count

    def add_taxa(self, taxa):
        """
        Adds multiple :class:`Taxon` objects to self.

        Each :class:`Taxon` object in `taxa` that is not already in the collection of
        :class:`Taxon` objects in this namespace is added to it. If any of the :class:`Taxon`
        objects are already in the collection, then nothing happens. If the
        namespace is immutable, then :class:`TypeError` is raised when trying
        to add :class:`Taxon` objects.

        Parameters
        ----------
        taxa : collections.Iterable [:class:`Taxon`]
            A list of :class:`Taxon` objects to be accessioned or registered in this
            collection.

        Raises
        ------
        TypeError
            If this namespace is immutable (i.e. `TaxonNamespace.is_mutable` is
            `False`).
        """
        for t in taxa:
            self.add_taxon(t)

    def new_taxon(self, label):
        """
        Creates, adds, and returns a new :class:`Taxon` object with corresponding
        label.

        Parameters
        ----------
        label : string or string-like
            The name or label of the new operational taxonomic unit concept.

        Returns
        -------
        taxon: :class:`Taxon`
            The new :class:`Taxon` object,

        """
        if not self.is_mutable:
            raise TypeError("Taxon '{}' cannot be added to an immutable TaxonNamespace".format(label))
        taxon = Taxon(label=label)
        self.add_taxon(taxon)
        return taxon

    def new_taxa(self, labels):
        """
        Creates and add a new :class:`Taxon` with corresponding label for each label
        in `labels`. Returns list of :class:`Taxon` objects created.

        Parameters
        ----------
        labels : :py:class:`collections.Iterable` [string]
            The values of the `label` attributes of the new :class:`Taxon` objects to
            be created, added to this namespace collection, and returned.

        Returns
        -------
        taxa : :py:class:`collections.Iterable` [:class:`Taxon`]
            A list of :class:`Taxon` objects created and added.

        Raises
        ------
        TypeError
            If this namespace is immutable (i.e. `TaxonNamespace.is_mutable` is
            `False`).

        """
        if not self.is_mutable:
            raise TypeError("Taxon objects cannot be added to an immutable TaxonNamespace")
        taxa = []
        for label in labels:
            taxa.append(self.new_taxon(label=label))
        return taxa

    ###########################################################################
    ## Removing Taxa

    def remove_taxon(self, taxon):
        """
        Removes specified :class:`Taxon` object from the collection in this namespace.

        Parameters
        ----------
        taxon : a :class:`Taxon` object
            The :class:`Taxon` object to be removed.

        Raises
        ------
        ValueError
            If `taxon` is not in the collection of this namespace.
        """
        if taxon not in self._taxa:
            raise ValueError(taxon)
        self._taxa.remove(taxon)
        # assert taxon not in self._taxa
        while taxon in self._taxa:
            self._taxa.remove(taxon)
        idx = self._taxon_accession_index_map.pop(taxon, None)
        if idx is not None:
            self._accession_index_taxon_map.pop(idx, None)
            self._taxon_accession_index_map.pop(taxon, None)
        bm = self._taxon_bitmask_map.pop(taxon, None)
        if bm is not None:
            self._bitmask_taxon_map.pop(bm, None)
            self._taxon_accession_index_map.pop(taxon, None)

    def remove_taxon_label(self, label, case_insensitive=False):
        """
        Removes *all* :class:`Taxon` objects with label matching `label` from the
        collection in this namespace.

        Parameters
        ----------
        label : string or string-like
            The value of the :class:`Taxon` object label to remove.
        case_insensitive : boolean, optional
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the :class:`Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Raises
        ------
        LookupError
            If no :class:`Taxon` objects are found with matching label(s).

        See Also
        --------
        :meth:`TaxonNamespace.discard_taxon_labels` : Similar, but does not raise an error if no
            matching :class:`Taxon` objects are found.
        """
        taxa = self._lookup_label(label,
                case_insensitive=case_insensitive,
                multiple=True,
                error_if_not_found=True)
        for taxon in taxa:
            self.remove_taxon(taxon)

    def discard_taxon_label(self, label, case_insensitive=False):
        """
        Removes *all* :class:`Taxon` objects with label matching `label` from the
        collection in this namespace.

        Parameters
        ----------
        label : string or string-like
            The value of the :class:`Taxon` object label to remove.
        case_insensitive : boolean, optional
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the :class:`Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        See Also
        --------
        :meth:`TaxonNamespace.discard_taxon_labels` : Similar, but does not
            raise an error if no matching :class:`Taxon` objects are found.
        """
        taxa = self._lookup_label(label,
                case_insensitive=case_insensitive,
                multiple=True,
                error_if_not_found=False)
        if taxa is None:
            return
        for taxon in taxa:
            self.remove_taxon(taxon)

    def clear(self):
        """
        Removes all :class:`Taxon` objects from this namespace.
        """
        # self._taxa.clear() # Python 2 `list` class does not have `clear()` method
        del self._taxa[:]
        self._accession_index_taxon_map.clear()
        self._taxon_accession_index_map.clear()
        self._taxon_bitmask_map.clear()
        self._bitmask_taxon_map.clear()

    ###########################################################################
    ## Look-up and Retrieval of Taxa

    def _lookup_label(self,
            label,
            multiple=True,
            case_insensitive=False,
            error_if_not_found=False):
        """
        Return :class:`Taxon` object(s) with label matching `label`.
        If `multiple` is `True`, then a list of :class:`Taxon` objects with labels
        that match `label` are returned, otherwise just the first one is
        returned. If `case_insensitive` is `True`, then the matching is done
        without regard for case. If no :class:`Taxon` object is in the current the
        namespace that matches the criteria, then `None` is returned unless
        `error_if_not_found` is `False`, in which case :class:`LookupError` is raised.
        """
        taxa = []
        if case_insensitive:
            label = str(label).lower()
            for taxon in self._taxa:
                if label == str(taxon.label).lower():
                    if not multiple:
                        return taxon
                    else:
                        taxa.append(taxon)
        else:
            for taxon in self._taxa:
                if label == taxon.label:
                    if not multiple:
                        return taxon
                    else:
                        taxa.append(taxon)
        if len(taxa) > 0:
            return taxa
        elif error_if_not_found:
            raise LookupError(label)
        else:
            return None

    def findall(self, label, case_insensitive=False):
        """
        Return list of :class:`Taxon` object(s) with label matching `label`.

        Parameters
        ----------
        label : string or string-like
            The value which the `label` attribute of the :class:`Taxon` object(s)
            to be returned must match.
        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the :class:`Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------
        taxa : :py:class:`list` [:class:`Taxon`]
            A list containing zero or more :class:`Taxon` objects with labels
            matching `label`.

        Raises
        ------
        TypeError
            If this namespace is immutable (i.e. `TaxonNamespace.is_mutable` is
            `False`).

        """
        taxa = self._lookup_label(
                label=label,
                multiple=True,
                case_insensitive=case_insensitive,
                error_if_not_found=False)
        if taxa is None:
            return []
        else:
            return taxa

    def has_taxon_label(self, label, case_insensitive=False):
        """
        Checks for presence of a :class:`Taxon` object with the given label.

        Parameters
        ----------
        label : string or string-like
            The value of the :class:`Taxon` object label to match.
        case_insensitive : boolean, optional
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the :class:`Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------
        b : boolean
            `True` if there is at least one :class:`Taxon` object in this namespace
            with a label matching the value of `label`. Otherwise, `False`.
        """
        t = self._lookup_label(
                label=label,
                multiple=False,
                case_insensitive=case_insensitive,
                error_if_not_found=False)
        return t is not None

    def has_taxa_labels(self, labels, case_insensitive=False):
        """
        Checks for presence of :class:`Taxon` objects with the given labels.

        Parameters
        ----------
        labels : :py:class:`collections.Iterable` [string]
            The values of the :class:`Taxon` object labels to match.
        case_insensitive : boolean, optional
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the :class:`Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------
        b : boolean
            Returns `True` if, for every element in the iterable `labels`,
            there is at least one :class:`Taxon` object that has a label attribute
            that matches this. `False` otherwise.
        """
        for label in labels:
            f = self._lookup_label(label=label,
                    multiple=False,
                    case_insensitive=case_insensitive,
                    error_if_not_found=False)
            if f is None:
                return False
        return True

    def get_taxon(self, label, case_insensitive=False):
        """
        Retrieves a :class:`Taxon` object with the given label.

        If `case_insensitive` is `True`, then the label matching is made
        without regard for case.  If multiple :class:`Taxon` objects exist with labels
        that match `label`, then only the first one is returned.  If no :class:`Taxon`
        object is found in this namespace with the specified critieria, `None`
        is returned.

        Parameters
        ----------
        label : string or string-like
            The value which the `label` attribute of the :class:`Taxon` object
            to be returned must match.
        case_insensitive : boolean,
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the :class:`Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------
        taxon : :class:`Taxon` object or `None`
            The first :class:`Taxon` object in this namespace collection with a label
            matching `label`, or `None` if no such :class:`Taxon` object exists.
        """
        return self._lookup_label(label=label,
                multiple=False,
                case_insensitive=case_insensitive,
                error_if_not_found=False)

    def get_taxa(self, labels, case_insensitive=False):
        """
        Retrieves list of :class:`Taxon` objects with given labels.

        If `case_insensitive` is `True`, then the label matching is made
        without regard for case.

        Parameters
        ----------
        labels : :py:class:`collections.Iterable` [string]
            Any :class:`Taxon` object in this namespace collection that has a label
            attribute that matches any value in `labels` will be included in
            the list returned.
        case_insensitive : boolean, optional
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the :class:`Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------
        taxa : :py:class:`list` [:class:`Taxon`]
            A list containing zero or more :class:`Taxon` objects with labels
            matching `label`.
        """
        taxa = []
        for label in labels:
            tt = self._lookup_label(label=label,
                    multiple=True,
                    case_insensitive=case_insensitive,
                    error_if_not_found=False)
            if tt is None:
                continue
            for t in tt:
                if t not in taxa:
                    taxa.append(t)
        return taxa

    def require_taxon(self, label, case_insensitive=False):
        """
        Retrieves a :class:`Taxon` object with the given label, creating it if
        necessary.

        Retrieves a Taxon object with the label, `label`.  If
        `case_insensitive` is `True`, then the label matching is made without
        regard for case.  If multiple :class:`Taxon` objects exist with labels that
        match `label`, then only the first one is returned.  If no such :class:`Taxon`
        object exists in the current namespace and the :class:`TaxonNamespace` is NOT
        mutable, an exception is raised.  If no such :class:`Taxon` object exists in
        the current namespace and :class:`TaxonNamespace` is mutable, then a new
        :class:`Taxon` is created, added, and returned.

        Parameters
        ----------
        label : string or string-like
            The value which the `label` attribute of the :class:`Taxon` object
            to be returned must match.
        case_insensitive : boolean, optional
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the :class:`Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------
        taxon : :class:`Taxon` object or `None`
            A :class:`Taxon` object in this namespace collection with a label
            matching `label`.

        Raises
        ------
        TypeError
            If no :class:`Taxon` object is currently in the collection with a label
            matching the input `label` and the `is_mutable` attribute of self
            is `False`.
        """
        taxon = self._lookup_label(label=label,
                multiple=False,
                case_insensitive=case_insensitive,
                error_if_not_found=False)
        if taxon is not None:
            return taxon
        if not self.is_mutable:
            raise TypeError("Taxon '{}' not in TaxonNamespace, and cannot be created because TaxonNamespace is immutable".format(label))
        taxon = self.new_taxon(label=label)
        return taxon

    ###########################################################################
    ## Taxon Ordering

    def sort(self, key=None, reverse=False):
        """
        Sorts :class:`Taxon` objects in collection. If `key` is not given, defaults
        to sorting by label (i.e., `key = lambda x: x.label`).

        Parameters
        ----------
        key : key function object, optional
            Function that takes a :class:`Taxon` object as an argument and
            returns the value that determines its sort order. Defaults to
            sorting by label.
        reverse : boolean, optional
            If `True`, sort will be in reverse order.
        """
        if key is None:
            key = lambda x: x.label
        self._taxa.sort(key=key, reverse=reverse)

    def reverse(self):
        """
        Reverses order of :class:`Taxon` objects in collection.
        """
        self._taxa.reverse()

    ###########################################################################
    ## Summarization of Collection

    def labels(self):
        """
        Returns list of labels of all :class:`Taxon` objects in `self`.

        Returns
        -------
        labels : :py:class:`list` [string]
            List of :attr:`Taxon.label` values of :class:`Taxon` objects in
            `self`.
        """
        return [t.label for t in self._taxa]

    def label_taxon_map(self, case_insensitive_keys=False):
        """
        Returns dictionary with taxon labels as keys and corresponding :class:`Taxon`
        objects as values.

        No attempt is made to handle collisions.

        Parameters
        ----------
        case_insensitive_keys : boolean, optional
            If `False` (default), then normal Python `dict` object will be
            returned, resulting in case-sensitive keys. If `True`, then a
            :class:`CaseInsensitiveDict` object will return, allowing for
            case-insensitive lookups.

        Returns
        -------
        d : :py:class:`dict` or :class:`CaseInsensitiveDict`
            Dictionary with :attr:`Taxon.label` values of :class:`Taxon` objects in
            `self` as keys and corresponding :class:`Taxon` objects as values.
        """
        if case_insensitive_keys:
            d = container.CaseInsensitiveDict()
        else:
            d = {}
        for t in self._taxa:
            d[t.label] = t
        return d

    ###########################################################################
    ## Split Management

    def complement_split_bitmask(self, split_bitmask):
        """
        Returns complement of the given split bitmask.

        Parameters
        ----------
        split_bitmask : integer
            Split bitmask hash to be complemented.

        Returns
        -------
        h : integer
            Complement of `split`.
        """
        return (~split) & self.all_taxa_bitmask()

    def all_taxa_bitmask(self):
        """
        Returns mask of all taxa.

        Returns
        -------
        h : integer
            Bitmask spanning all :class:`Taxon` objects in self.
        """
        #return pow(2, len(self)) - 1
        b = 1 << self._current_accession_count
        return b - 1

    def taxon_bitmask(self, taxon):
        """
        Returns bitmask value of split hash for split subtending node with
        `taxon`.

        Parameters
        ----------
        taxon : :class:`Taxon`
            :class:`Taxon` object for which to calculate split hash bitmask.

        Returns
        -------
        h : integer
            Split hash bitmask value for node associated with :class:`Taxon` object `taxon`.
        """
        try:
            return self._taxon_bitmask_map[taxon]
        except KeyError:
            i = self._taxon_accession_index_map[taxon]
            m = 1 << i
            self._taxon_bitmask_map[taxon] = m
            self._bitmask_taxon_map[m] = taxon
            return m

    def taxa_bitmask(self, **kwargs):
        """
        Retrieves the list of split hash bitmask values representing all taxa
        specified by keyword-specified list of taxon objects (`taxa=`) or
        labels (`labels=`).

        Parameters
        ----------
        \*\*kwargs : keyword arguments
            Requires one of:

                taxa : :py:class:`collections.Iterable` [:class:`Taxon`]
                    Iterable of :class:`Taxon` objects.
                labels : :py:class:`collections.Iterable` [string]
                    Iterable of :class:`Taxon` label values.

        Returns
        -------
        b : :py:class:`list` [integer]
            List of split hash bitmask values for specified :class:`Taxon`
            objects.
        """
        if "taxa" in kwargs:
            taxa = kwargs["taxa"]
        else:
            taxa = self.get_taxa(**kwargs)
        bitmask = 0
        for taxon in taxa:
            bitmask |= self.taxon_bitmask(taxon)
        return bitmask

    def split_bitmask_string(self, split_bitmask):
        """
        Returns bitstring representation of split_bitmask.

        Parameters
        ----------
        split_bitmask : integer
            Split hash bitmask value to be represented as a string.

        Returns
        -------
        s : string
            String representation of the split hash bitmask value passed as an
            argument.
        """
        return "{}".format(text.int_to_bitstring(split_bitmask).rjust(len(self._taxon_accession_index_map), "0"))

    def split_taxa_list(self, split_bitmask, index=0):
        """
        Returns list of :class:`Taxon` objects represented by split
        `split_bitmask`.

        Parameters
        ----------
        split_bitmask : integer
            Split hash bitmask value.
        index : integer, optional
            Start from this :class:`Taxon` object instead of the first
            :class:`Taxon` object in the collection.

        Returns
        -------
        taxa : :py:class:`list` [:class:`Taxon`]
            List of :class:`Taxon` objects specified or spanned by
            `split_bitmask`.
        """
        taxa = []
        while split_bitmask:
            if split_bitmask & 1:
                taxa.append(self._taxon_accession_index_map[index])
            split_bitmask = split_bitmask >> 1
            index += 1
        return taxa

    def split_as_newick_string(self,
            split_bitmask,
            preserve_spaces=False,
            quote_underscores=True):
        """
        Represents a split as a newick string.

        Parameters
        ----------
        split_bitmask : integer
            Split hash bitmask value.
        preserve_spaces : boolean, optional
            If `False` (default), then spaces in taxon labels will be replaced
            by underscores. If `True`, then taxon labels with spaces will be
            wrapped in quotes.
        quote_underscores : boolean, optional
            If `True` (default), then taxon labels with underscores will be
            wrapped in quotes. If `False`, then the labels will not be wrapped
            in quotes.

        Returns
        -------
        s : string
            NEWICK representation of split specified by `split_bitmask`.
        """
        return text.split_as_newick_string(split, self, preserve_spaces=preserve_spaces, quote_underscores=quote_underscores)

##############################################################################
## TaxonSet

class TaxonSet(TaxonNamespace):
    """
    This class is present for (temporary!) legacy support of code written under
    DendroPy 3.x.  It will be removed in future versions. All new code should
    be written using :class:`TaxonNamespace`. Old code needs to be updated to use
    :class:`TaxonNamespace`.
    """
    def __new__(cls):
        error.dump_stack()
        warnings.warn(":class:`TaxonSet` will no longer be supported in future releases; use :class:`TaxonNamespace` instead",
                FutureWarning, stacklevel=3)
        o = super(TaxonSet, cls).__new__(cls)
        return o

##############################################################################
## Taxon

class Taxon(base.DataObject, base.Annotable):
    """
    A taxon associated with a sequence or a node on a tree.
    """

    def __init__(self, label):
        """
        Parameters
        ----------
        label : string or :class:`Taxon` object
            Label or name of this operational taxonomic unit concept. If a
            string, then the `label` attribute of `self` is set to this value.
            If a :class:`Taxon` object, then the `label` attribute of `self` is
            set to the same value as the `label` attribute the other
            :class:`Taxon` object and all annotations/metadata are copied.
        """
        if isinstance(label, Taxon):
            other_taxon = label
            label = other_taxon.label
            memo={id(other_taxon):self}
            for k in other_taxon.__dict__:
                if k != "_annotations":
                    self.__dict__[k] = copy.deepcopy(other_taxon.__dict__[k], memo=memo)
            self.deep_copy_annotations_from(other_taxon, memo=memo)
        else:
            label = str(label)
            base.DataObject.__init__(self, label=label)

    def __copy__(self):
        raise TypeError("Shallow copies of Taxon objects are not allowed")

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        try:
            o = memo[id(self)]
        except KeyError:
            # o = type(self).__new__(self.__class__)
            o = self.__class__.__new__(self.__class__)
            memo[id(self)] = o
        for k in self.__dict__:
            if k != "_annotations":
                o.__dict__[k] = copy.deepcopy(self.__dict__[k], memo)
        o.deep_copy_annotations_from(self, memo)
        return o

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def __lt__(self, other):
        return self.label < other.label

    def __str__(self):
        "String representation of self = taxon name."
        return "'{}'".format(self._label)

    def __repr__(self):
        return "<{} {} '{}'>".format(self.__class__.__name__, hex(id(self)), self._label)

    def description(self, depth=1, indent=0, itemize="", output=None, **kwargs):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return ""
        output_strio = StringIO()
        if self._label is None:
            label = "<Unnamed Taxon>"
        else:
            label = "'{}'".format(self._label)
        output_strio.write('{}{} Taxon object at {}: {}' % (indent*' ', itemize, hex(id(self)), label))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s
