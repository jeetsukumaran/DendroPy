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

    * Operational taxonomic units are modeled by a `Taxon` object.

    * Taxonomic namespaces, in which operational taxonomic units are organized,
      are modeled by a `TaxonNamespace` object.

    * A `TaxonNamespace` manages a collection of `Taxon` objects, where each
      object represents a distinct operational taxonomic unit concept within
      the taxonomic namespace represented by that `TaxonNamespace` object.

    * Each `Taxon` object can belong to one and only one `TaxonNamespace`:
      `Taxon` objects are not shared across `TaxonNamespace` objects.

    * Each `Taxon` object has an attribute, `label`, whose (string) value
      is the name of the operational taxon unit concept that it represents.

    * Different `Taxon` objects represent different operational taxonomic
      unit concepts, even if they have the same label value.

    * All client objects (`TaxonNamespaceAssociated` objects) that reference
      the same `TaxonNamespace` reference the same "universe" or domain of
      operational taxonomic unit concepts.

    * Operational taxonomic units from across different data sources are mapped
      to distinct `Taxon` objects within a particular `TaxonNamespace` based on
      matching the string values of labels of the `Taxon` object.

    * A particular taxonomic unit concept in one data source will only be
      correctly related to the same taxonomic unit concept (i.e, the same
      `Taxon` object) in another data source only if they have both
      been parsed with reference to the same taxonomic namespace (i.e., the
      same `TaxonNamespace` has been used).

    * A `TaxonNamespace` assigned an "accession index" to every `Taxon` object
      added to it. This is a stable and unique number within the context of any
      given `TaxonNamespace` object (though a `Taxon` object may have different
      accession indexes in different `TaxonNamespace` objects if it
      belongs to multiple namespaces). This number is will be used to
      calculate the "split bitmask" hash of the trivial split or external edge
      subtending the node to which this `Taxon` object is assigned on a tree.
      The concept of a "split bitmask" hash is fundamental to DendroPy's tree
      operations. The split bitmask is a hash that uniquely identifies every
      split on a tree.  It is calculated by OR'ing the split bitmask of all the
      child splits of the given split. Terminal edges, of course, do not have
      child edges, and their split bitmask is given by the accession index of
      the `Taxon` object at their head or target nodes.
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
        Rebuilds `taxon_namespace` from scratch, or assigns `Taxon` objects from
        given `TaxonNamespace` object `taxon_namespace` based on label values. Calls
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
        components, attributes and members all refer to the same `TaxonNamespace`
        object as `self.taxon_namespace`, and that `self.taxon_namespace` has all
        the `Taxon` objects in the various members.
        """
        pass

##############################################################################
## TaxonNamespace

class TaxonNamespace(base.DataObject, base.Annotable):
    """
    A collection of `Taxon` objects representing a self-contained and complete
    domain of distinct operational taxonomic unit definitions.
    Provides the common semantic context in which operational taxonomic units
    referenced by various phylogenetic data objects (e.g., trees or alignments)
    can be related.
    """

    ## Construction, Copying and Life-cycle ##

    def __init__(self, *args, **kwargs):
        """
        Constructs a TaxonNamespace object.

        Parameters
        ----------

        *args : positional arguments
            Accepts a single iterable as an optional positional argument.  If a
            `TaxonNamespace` object is passed as the positional argument, then
            clones or deep-copies of its member `Taxon` objects will be added
            to this one.  If any other iterable is passed as the positional
            argument, then each string in the iterable will result in a new
            `Taxon` object being constructed and added to the namespace with
            the string as its label (name), while each Taxon object in the
            iterable will be added to the namespace directly.

        **kwargs : {label, is_mutable}
            `label` : string
                The label or name for this namespace.
            `full_taxon_object_model` : boolean (default = `True`)
                If `True`, then operational taxonomic unit concepts are
                represented by `Taxon` objects. This is rich model, supporting,
                e.g. metadata annotation, complex attributes, etc. If `False`,
                then operational taxonomic unit concepts are represented by
                simple strings. Using this results in *much* faster
                reading/writing operations.
            `is_mutable` : boolean, optional (default = `True`)
                If `True` (default), then `Taxon` objects can be added to this
                namespace. If `False`, then adding `Taxon` objects will result
                in an error.

        Notes
        -----

        An empty `TaxonNamespace` can be created (with optional) label and `Taxon`
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
        `Taxon` objects will be added directly while, for each string, a new
        `Taxon` object will be created and added. So, the below are all equivalent
        to the above:

        >>> tns = dendropy.TaxonNamespace(["a", "b", "c"], label="taxa")

        >>> taxa = [Taxon(n) for n in ["a", "b", "c"]]
        >>> tns = dendropy.taxonnamespace(taxa, label="taxa")

        >>> t1 = Taxon("a")
        >>> t2 = Taxon("b")
        >>> taxa = [t1, t2, "c"]
        >>> tns = dendropy.TaxonNamespace(taxa, label="taxa")

        If a `TaxonNamespace` object is passed as the initializer argument,
        then each `Taxon` object in the original `TaxonNamespace` will be
        *fully* cloned, i.e., will become completely separate and independent
        operational taxonomic concepts:

        >>> tns1 = dendropy.TaxonNamespace(["a", "b", "c"], label="taxa1")
        >>> tns2 = dendropy.TaxonNamespace(tns1, label="2")
        >>> tns1
        <TaxonNamespace 0x1097275d0 'taxa1': [<Taxon 0x109727610 'a'>, <Taxon 0x109727e10 'b'>, <Taxon 0x109727e90 'c'>]>
        >>> tns2
        <TaxonNamespace 0x109727d50 '2': [<Taxon 0x10972a0d0 'a'>, <Taxon 0x10972a090 'b'>, <Taxon 0x10972ac10 'c'>]>
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
            if isinstance(args[0], TaxonNamespace):
                self.label = kwargs.get('label', args[0].label)
                memo = {}
                for t in args[0]:
                    self.add_taxon(t.fullcopy(memo=memo))
                self.copy_annotations_from(args[0], memo=memo)
            else:
                for i in args[0]:
                    if isinstance(i, Taxon):
                        self.add_taxon(i)
                    else:
                        self.new_taxon(label=i)

    def __str__(self):
        return "[{}]".format(", ".join([str(i) for i in self._taxa]))

    def __repr__(self):
        return "<{} {} '{}': [{}]>".format(self.__class__.__name__, hex(id(self)), self.label, ", ".join(repr(i) for i in self._taxa))

    def __hash__(self):
        return hash( (t for t in self._taxa) )

    def __iter__(self):
        return iter(self._taxa)

    def __reversed__(self):
        return reversed(self._taxa)

    # def __deepcopy__(self, memo):
    #     """
    #     By default, deep copies of non-DataSet data objects do *not* deep-copy
    #     the corresponding TaxonNamespace, and id's of all TaxonNamespace
    #     objects are mapped to themselves. This can be overridden by
    #     pre-populating memo with appropriate clones.
    #     """
    #     memo[id(self)] = self
    #     return self

    # def fullcopy(self, memo=None):
    #     """
    #     "Truly" deep-copy or clone the TaxonNamespace: make copies of Taxon
    #     objects and return new TaxonNamespace.
    #     """
    #     raise NotImplementedError
    #     if memo is None:
    #         memo = {}
    #     o = base.AnnotatedDataObject.__deepcopy__(self, memo)
    #     for taxon in self._taxa:
    #         o.add_taxon(taxon.fullcopy(memo=memo))
    #     for k in self.__dict__:
    #         # o.__dict__[copy.deepcopy(k, memo)] = copy.deepcopy(self.__dict__[k], memo)
    #         o.__dict__[k] = copy.deepcopy(self.__dict__[k], memo)
    #     memo[id(self)] = o
    #     return o

    ## Container Interface ##

    def __len__(self):
        """
        Returns number of `Taxon` objects in this `TaxonNamespace`.
        """
        return len(self._taxa)

    def __getitem__(self, key):
        """
        Returns `Taxon` object with index or slice given by `key`.
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
        return taxon in self._taxa

    def _lookup_label(self,
            label,
            multiple=True,
            case_insensitive=False,
            error_if_not_found=False):
        """
        Return `Taxon` object(s) with label matching `label`.
        If `multiple` is `True`, then a list of `Taxon` objects with labels
        that match `label` are returned, otherwise just the first one is
        returned. If `case_insensitive` is `True`, then the matching is done
        without regard for case. If no `Taxon` object is in the current the
        namespace that matches the criteria, then `None` is returned unless
        `error_if_not_found` is `False`, in which case `LookupError` is raised.
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

    # def _resolve_taxon_request(self,
    #         key,
    #         multiple=True,
    #         case_insensitive=False,
    #         error_if_not_found=True):
    #     """
    #     If `key` is a `Taxon` object, return `key` as-is.  Otherwise, return the
    #     `Taxon` object with label value `key`.  If `Taxon` object cannot be
    #     resolved, then `None` is returned if `raise_error` is `True`,
    #     otherwise an IndexError or LookupError is raised.
    #     """
    #     if isinstance(key, Taxon):
    #         if key in self._taxa:
    #             if multiple:
    #                 return [key]
    #             else:
    #                 return key
    #         elif error_if_not_found:
    #             raise LookupError(key)
    #         else:
    #             return None
    #     else:
    #         return self.findall(
    #                 label=key,
    #                 multiple=multiple,
    #                 case_insensitive=case_insensitive,
    #                 error_if_not_found=error_if_not_found)

    def add_taxon(self, taxon):
        """
        Adds a new Taxon object to self.

        If `taxon` is not already in the collection of `Taxon` objects in this
        namespace, and this namespace is mutable, it is added to the
        collection. If it is already in the collection, then nothing happens.
        If it is not already in the collection, but the namespace is not
        mutable, then a `TypeError` is raised.

        Parameters
        ----------

        taxon : `Taxon` object (if full taxon model) or string (if not a full taxon model)
            The `Taxon` object to be accessioned or registered in this
            collection.

        Raises
        ------
        TypeError
            If this namespace is immutable (i.e. `TaxonNamespace.is_mutable` is
            `False`).

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
        Adds multiple `Taxon` objects to self.

        Each `Taxon` object in `taxa` that is not already in the collection of
        `Taxon` objects in this namespace is added to it. If any of the `Taxon`
        objects are already in the collection, then nothing happens.  If the
        namespace is immutable, then `TypeError` is raised when trying to add
        `Taxon` objects.

        Parameters
        ----------

        taxa : an iterable of `Taxon` objects or strings (if not a full taxon model)
            A list of `Taxon` objects to be accessioned or registered in this
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
        Creates, adds, and returns a new `Taxon` object with corresponding
        label (if a full taxon model) or simply adds the label as-is (if not a
        full taxon model).

        Parameters
        ----------

        label : string or string-like
            The name or label of the new operational taxonomic unit concept.

        Returns
        -------

        Taxon object or string
            A new `Taxon` object (if a full taxon model) or string (if not a
            full taxon model).

        """
        if not self.is_mutable:
            raise TypeError("Taxon '{}' cannot be added to an immutable TaxonNamespace".format(label))
        taxon = Taxon(label=label)
        self.add_taxon(taxon)
        return taxon

    def new_taxa(self, labels):
        """
        Creates and add a new `Taxon` with corresponding label for each label
        in `labels`. Returns list of `Taxon` objects created.

        Parameters
        ----------

        labels : an interable of string or string-like
            The values of the `label` attributes of the new `Taxon` objects to
            be created, added to this namespace collection, and returned.

        Returns
        -------

        List of `Taxon` objects.
            A list of the newly-created `Taxon` objects.

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


    def findall(self, label, case_insensitive=False):
        """
        Return list of `Taxon` object(s) with label matching `label`.

        Parameters
        ----------

        label : string or string-like
            The value which the `label` attribute of the `Taxon` object(s)
            to be returned must match.

        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the `Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------

        A list of zero or more `Taxon` objects with labels matching `label`.

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
        Checks for presence of a `Taxon` object with the given label.

        Parameters
        ----------

        label : string or string-like
            The value of the `Taxon` object label to match.

        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the `Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------
        boolean
            `True` if there is at least one `Taxon` object in this namespace
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
        Checks for presence of `Taxon` objects with the given labels.

        Parameters
        ----------

        labels : iterable
            The values of the `Taxon` object labels to match.

        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the `Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------
        boolean
            Returns `True` if, for every element in the iterable `labels`,
            therer is at least one `Taxon` object that has a label attribute
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
        Retrieves a `Taxon` object with the given label.

        If `case_insensitive` is `True`, then the label matching is made
        without regard for case.  If multiple `Taxon` objects exist with labels
        that match `label`, then only the first one is returned.  If no `Taxon`
        object is found in this namespace with the specified critieria, `None`
        is returned.

        Parameters
        ----------

        label : string or string-like
            The value which the `label` attribute of the `Taxon` object
            to be returned must match.

        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the `Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------

        Taxon object or None
            The first `Taxon` object in this namespace collection with a label
            matching `label`, or `None` if no such `Taxon` object exists.

        """
        return self._lookup_label(label=label,
                multiple=False,
                case_insensitive=case_insensitive,
                error_if_not_found=False)

    def get_taxa(self, labels, case_insensitive=False):
        """
        Retrieves list of `Taxon` objects with given labels.

        If `case_insensitive` is `True`, then the label matching is made
        without regard for case.

        Parameters
        ----------

        labels : iterable
            Any `Taxon` object in this namespace collection that has a label
            attribute that matches any value in `labels` will be included in
            the list returned.

        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the `Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------

        A list of zero or more `Taxon` objects, each of which has a `label`
        attribute that matches one of the values in `labels`.

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
        Retrieves a `Taxon` object with the given label, creating it if
        necessary.

        Retrieves a Taxon object with the label, `label`.  If
        `case_insensitive` is `True`, then the label matching is made without
        regard for case.  If multiple `Taxon` objects exist with labels that
        match `label`, then only the first one is returned.  If no such `Taxon`
        object exists in the current namespace and the `TaxonNamespace` is NOT
        mutable, an exception is raised.  If no such `Taxon` object exists in
        the current namespace and `TaxonNamespace` is mutable, then a new
        `Taxon` is created, added, and returned.

        Parameters
        ----------

        label : string or string-like
            The value which the `label` attribute of the `Taxon` object
            to be returned must match.

        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the `Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Returns
        -------

        Taxon object or None
            A `Taxon` object in this namespace collection with a label
            matching `label`.

        Raises
        ------

        TypeError
            If no `Taxon` object is currently in the collection with a label
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

    def remove_taxon(self, taxon):
        """
        Removes specified `Taxon` object from the collection in this namespace.

        Parameters
        ----------

        taxon : a `Taxon` object
            The `Taxon` object to be removed.

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
        Removes *all* `Taxon` objects with label matching `label` from the
        collection in this namespace.

        Parameters
        ----------

        label : string or string-like
            The value of the `Taxon` object label to remove.

        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the `Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        Raises
        ------

        LookupError
            If no `Taxon` objects are found with matching label(s).

        See Also
        --------

        discard_taxon_labels : Similar, but does not raise an error if no
                               matching `Taxon` objects are found.

        """
        taxa = self._lookup_label(label,
                case_insensitive=case_insensitive,
                multiple=True,
                error_if_not_found=True)
        for taxon in taxa:
            self.remove_taxon(taxon)

    def discard_taxon_label(self, label, case_insensitive=False):
        """
        Removes *all* `Taxon` objects with label matching `label` from the
        collection in this namespace.

        Parameters
        ----------

        label : string or string-like
            The value of the `Taxon` object label to remove.

        case_insensitive : boolean, optional (default = `False`)
            If `False` (default), then the label matching is done as-is. If
            `True`, then both the `label` argument as well as the `Taxon`
            object's `label` attribute are coerced into lower-case label
            strings before checking for a match.

        See Also
        --------

        remove_taxon_labels : Similar, but raises an error if no
                              matching `Taxon` objects are found.

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
        Removes all `Taxon` objects from this namespace.
        """
        # self._taxa.clear() # Python 2 `list` class does not have `clear()` method
        del self._taxa[:]
        self._accession_index_taxon_map.clear()
        self._taxon_accession_index_map.clear()
        self._taxon_bitmask_map.clear()
        self._bitmask_taxon_map.clear()

    def sort(self, key=None):
        """
        Sorts `Taxon` objects in collection. If `key` is not given, defaults
        to sorting by label (i.e., `key = lambda x: x.label`).
        """
        if key is None:
            key = lambda x: x.label
        self._taxa.sort(key=key)

    def reverse(self):
        """
        Reverses order of `Taxon` objects in collection.
        """
        self._taxa.reverse()

    def labels(self):
        "Convenience method to return all taxa labels."
        return [t.label for t in self._taxa]

    def label_taxon_map(self, case_insensitive_keys=False):
        """
        Returns dictionary with taxon labels as keys and corresponding `Taxon`
        objects as values.

        No attempt is made to handle collisions.

        Parameters
        ----------

        case_insensitive_keys : boolean
            If `False` (default), then normal Python `dict` object will be
            returned, resulting in case-sensitive keys. If `True`, then a
            `CaseInsensitiveDict` object will return, allowing for
            case-insensitive lookups.

        Returns
        -------
        dict or CaseInsensitiveDict
        """
        if case_insensitive_keys:
            d = container.CaseInsensitiveDict()
        else:
            d = {}
        for t in self._taxa:
            d[t.label] = t
        return d

    def complement_split_bitmask(self, split):
        "Returns complement of the split bitmask."
        return (~split) & self.all_taxa_bitmask()

    def all_taxa_bitmask(self):
        "Returns mask of all taxa."
        #return pow(2, len(self)) - 1
        b = 1 << self._current_accession_count
        return b - 1

    def taxon_bitmask(self, taxon):
        """
        Returns unique bitmask of given taxon. Will raise index error if
        taxon does not exist.
        """
        try:
            return self._taxon_bitmask_map[taxon]
        except KeyError:
            i = self._taxon_accession_index_map[taxon]
            m = 1 << i
            self._taxon_bitmask_map[taxon] = m
            self._bitmask_taxon_map[m] = taxon
            return m

    def get_taxa_bitmask(self, **kwargs):
        """
        Retrieves bitmask represent all taxa specified by keyword-specified
        list of taxon objects (`taxa=`) or labels (`labels=`).
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
        "Returns bitstring representation of split_bitmask."
        return "{}".format(text.int_to_bitstring(split_bitmask).rjust(len(self), "0"))

    def split_taxa_list(self, split_bitmask, index=0):
        "Returns list of taxa represented by split."
        taxa = []
        while split_bitmask:
            if split_bitmask & 1:
                taxa.append(self[index])
            split_bitmask = split_bitmask >> 1
            index += 1
        return taxa

    def split_as_newick_string(self, split, preserve_spaces=False, quote_underscores=True):
        """
        Represents a split as a newick string.
        """
        return text.split_as_newick_string(split, self, preserve_spaces=preserve_spaces, quote_underscores=quote_underscores)

##############################################################################
## TaxonSet
class TaxonSet(TaxonNamespace):
    """
    This class is present for (temporary!) legacy support of code written under
    DendroPy 3.x.  It will be removed in future versions. All new code should
    be written using `TaxonNamespace`. Old code needs to be updated to use
    `TaxonNamespace`.
    """
    def __new__(cls):
        error.dump_stack()
        warnings.warn("`TaxonSet` will no longer be supported in future releases; use `TaxonNamespace` instead",
                FutureWarning, stacklevel=3)
        o = super(TaxonSet, cls).__new__(cls)
        return o

##############################################################################
## SimpleTaxonNamespace

# class SimpleTaxonNamespace(TaxonNamespace):
#     """
#     A collection of `Taxon` objects representing a self-contained and complete
#     domain of distinct operational taxonomic unit definitions.
#     Provides the common semantic context in which operational taxonomic units
#     referenced by various phylogenetic data objects (e.g., trees or alignments)
#     can be related.
#     """

#     ## Construction, Copying and Life-cycle ##

#     def __init__(self, *args, **kwargs):
#         """
#         Constructs a TaxonNamespace object.

#         Parameters
#         ----------

#         *args : positional arguments
#             Accepts a single iterable as an optional positional argument.  If a
#             `TaxonNamespace` object is passed as the positional argument, then
#             clones or deep-copies of its member `Taxon` objects will be added
#             to this one.  If any other iterable is passed as the positional
#             argument, then each string in the iterable will result in a new
#             `Taxon` object being constructed and added to the namespace with
#             the string as its label (name), while each Taxon object in the
#             iterable will be added to the namespace directly.

#         **kwargs : {label, is_mutable}
#             `label` : string
#                 The label or name for this namespace.
#             `is_mutable` : boolean, optional (default = `True`)
#                 If `True` (default), then `Taxon` objects can be added to this
#                 namespace. If `False`, then adding `Taxon` objects will result
#                 in an error.

#         """
#         TaxonNamespace.__init__(self, **kwargs)

#     def _lookup_label(self,
#             label,
#             multiple=True,
#             case_insensitive=False,
#             error_if_not_found=False):
#         """
#         Return `Taxon` object(s) with label matching `label`.

#         """
#         if case_insensitive:
#             raise NotImplementedError
#         if label in self._taxa:
#             return label
#         elif error_if_not_found:
#             raise LookupError(label)
#         else:
#             return None

#     def new_taxon(self, label):
#         """
#         Creates, adds, and returns a new `Taxon` object with corresponding
#         label (if a full taxon model) or simply adds the label as-is (if not a
#         full taxon model).

#         Parameters
#         ----------

#         label : string or string-like
#             The name or label of the new operational taxonomic unit concept.

#         Returns
#         -------

#         Taxon object or string
#             A new `Taxon` object (if a full taxon model) or string (if not a
#             full taxon model).

#         """
#         self.add_taxon(label)
#         return label

#     def sort(self, key=None):
#         """
#         Sorts `Taxon` objects in collection. If `key` is not given, defaults
#         to sorting by label (i.e., `key = lambda x: x.label`).
#         """
#         if key is None:
#             key = lambda x: x
#         self._taxa.sort(key=key)

#     def reverse(self):
#         """
#         Reverses order of `Taxon` objects in collection.
#         """
#         self._taxa.reverse()

#     def labels(self):
#         "Convenience method to return all taxa labels."
#         return list(self._taxa)

#     def label_taxon_map(self, case_insensitive_keys=False):
#         """
#         Returns dictionary with taxon labels as keys and corresponding `Taxon`
#         objects as values.

#         No attempt is made to handle collisions.

#         Parameters
#         ----------

#         case_insensitive_keys : boolean
#             If `False` (default), then normal Python `dict` object will be
#             returned, resulting in case-sensitive keys. If `True`, then a
#             `CaseInsensitiveDict` object will return, allowing for
#             case-insensitive lookups.

#         Returns
#         -------
#         dict or CaseInsensitiveDict
#         """
#         return dict(zip(self._taxa, self._taxa))

##############################################################################
## Taxon

class Taxon(base.DataObject, base.Annotable):
    """
    A taxon associated with a sequence or a node on a tree.
    """

    def __init__(self, label):
        if isinstance(label, Taxon):
            label = label.label
        else:
            label = str(label)
        base.DataObject.__init__(self, label=label)
        base.Annotable.__init__(self)

    def __lt__(self, other):
        return self.label < other.label

    def __deepcopy__(self, memo):
        """
        By default, deep copies of non-DataSet data objects do *not* deep-copy
        the taxa, and id's of all taxon set objects are mapped to self. This
        can be overridden by pre-populating memo with appropriate clones.
        """
        memo[id(self)] = self
        return self

    def fullcopy(self, memo=None):
        if memo is None:
            memo = {}
        o = base.AnnotatedDataObject.__deepcopy__(self, memo)
        memo[id(self)] = o
        return o

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
