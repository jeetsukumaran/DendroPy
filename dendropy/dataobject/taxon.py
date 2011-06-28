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
This module provides classes and methods for managing taxa.
"""

import sys
import math
from cStringIO import StringIO
from dendropy.dataobject import base
from dendropy.utility import textutils
from dendropy.utility import containers
from dendropy.utility import error

def new_taxon_set(ntax=10, label_func=None):
    """
    Generates a new set of Taxon objects. `label_func` should be a function
    which can take an integer representing a new Taxon object's index in the
    new TaxonSet as an argument and return an appropriate (unique) label for
    the taxon. Note that an alternative way of instantiating a set of taxa
    would be to call `TaxonSet` with a list of labels as an argument.
    """
    taxon_set = TaxonSet()
    if label_func is None:
        label_idx_length = int(math.log(ntax, 10)) + 1
        label_template = "T%%0%dd" % (label_idx_length)
        label_func = lambda x: label_template % x
    for i in range(ntax):
        taxon_set.new_taxon(label=label_func(i+1))
    return taxon_set

class TaxonLinked(base.IdTagged):
    """
    Provides infrastructure for maintaining link/reference to a Taxon
    object.
    """

    def __init__(self, **kwargs):
        base.IdTagged.__init__(self,
                               label=kwargs.get("label", None),
                               oid=kwargs.get("oid", None))
        self.taxon = kwargs.get("taxon", None)

    def __deepcopy__(self, memo):
        """
        By default, deep copies of non-DataSet data objects do *not* deep-copy
        the taxa, and id's of all taxon set objects are mapped to self. This
        can be overridden by pre-populating memo with appropriate clones.
        """
        if id(self.taxon) in memo:
            o = self.__class__(taxon=memo[id(self.taxon)], label=self.label)
        else:
            o = self.__class__(taxon=self.taxon, label=self.label)
            memo[id(self.taxon)] = o.taxon
        memo[id(self)] = o
        memo[id(self._oid)] = o._oid
        return o

class TaxonSetLinked(base.IdTagged):
    """
    Provides infrastructure for the maintenance of references to taxa.
    """

    def __init__(self, **kwargs):
        base.IdTagged.__init__(self,
                               label=kwargs.get("label", None),
                               oid=kwargs.get("oid", None))
        if "taxon_set" not in kwargs or kwargs["taxon_set"] is None:
            self.taxon_set = TaxonSet()
        else:
            self.taxon_set = kwargs["taxon_set"]

    def reindex_taxa(self, taxon_set=None, clear=False):
        """
        Rebuilds `taxon_set` from scratch, or assigns `Taxon` objects from
        given `TaxonSet` object `taxon_set` based on label values. Calls
        on `self.reindex_member_taxa()` to synchronize taxa.
        """
        if taxon_set is not None:
            self.taxon_set = taxon_set
        if clear:
            self.taxon_set.clear()
        self.reindex_subcomponent_taxa()
        return self.taxon_set

    def reindex_subcomponent_taxa():
        """
        Derived classes should override this to ensure that their various
        components, attributes and members all refer to the same `TaxonSet`
        object as `self.taxon_set`, and that `self.taxon_set` has all
        the `Taxon` objects in the various members.
        """
        pass

    def __deepcopy__(self, memo):
        """
        By default, deep copies of non-DataSet data objects do *not* deep-copy
        the taxa, and id's of all taxon set objects are mapped to self. This
        can be overridden by pre-populating memo with appropriate clones.
        """
        if id(self.taxon_set) in memo:
            o = self.__class__(taxon_set=memo[id(self.taxon_set)], label=self.label)
        else:
            o = self.__class__(taxon_set=self.taxon_set, label=self.label)
            memo[id(self.taxon_set)] = o.taxon_set
        memo[id(self)] = o
        memo[id(self._oid)] = o._oid
        for i, t in enumerate(self.taxon_set):
            if id(t) not in memo:
                memo[id(t)] = o.taxon_set[i]
        return o

class TaxonSet(containers.OrderedSet, base.IdTagged):
    """
    Primary manager for collections of `Taxon` objects.
    """

    def _to_taxon(s):
        if isinstance(s, Taxon):
            return Taxon(label=s.label)
        if isinstance(s, str):
            return Taxon(label=s)
        raise ValueError("Cannot convert %s to Taxon" % str(s))
    _to_taxon = staticmethod(_to_taxon)

    def __init__(self, *args, **kwargs):
        """
        __init__ handles keyword arguments: `oid`, `label` or "is_mutable".
        If an iterable is passed as the first argument, then for every
        string in the iterable a Taxon object with the string is
        constructed and added to the set, while for every Taxon object
        in the iterable a new (distinct) Taxon object with the same
        label is constructed and added to the set.
        """
        containers.OrderedSet.__init__(self)
        base.IdTagged.__init__(self, oid=kwargs.get('oid'), label=kwargs.get('label'))
        if len(args) > 1:
            raise TypeError("TaxonSet() takes at most 1 non-keyword argument (%d given)" % len(args))
        elif len(args) == 1:
            if isinstance(args[0], TaxonSet):
                self.label = kwargs.get('label', args[0].label)
            for i in args[0]:
                if isinstance(i, Taxon):
                    self.add(i)
                else:
                    self.add(TaxonSet._to_taxon(i))
        self._is_mutable = kwargs.get('is_mutable', True) # immutable constraints not fully implemented -- only enforced at the add_taxon stage)

    def __deepcopy__(self, memo):
        o = self.__class__(list(self), label=self.label, is_mutable=self._is_mutable)
        memo[id(self)] = o
        memo[id(self._oid)] = o._oid
        return o

    def __getitem__(self, i):
        if isinstance(i, int):
            return containers.OrderedSet.__getitem__(self, i)
        else:
            for t in self:
                if t.label is i:
                    return t
        raise KeyError(i)

    def lock(self):
        self._is_mutable = False

    def unlock(self):
        self._is_mutable = True

    def get_is_locked(self):
        return self._is_mutable

    def set_is_locked(self, v):
        self._is_mutable = bool(v)
    is_locked = property(get_is_locked, set_is_locked)

    def has_taxon(self, **kwargs):
        """
        Returns True if taxon `taxon`, or with `oid` or `label`,
        exists (supplied by keywords; matches any)
        """
        if "taxon" not in kwargs and "oid" not in kwargs and "label" not in kwargs:
            raise TypeError("Need to specify oid or label.")
        req_taxon = kwargs.get("taxon", None)
        oid = kwargs.get("oid", None)
        label = kwargs.get("label", None)
        for self_taxon in self:
            if (req_taxon is not None and req_taxon is self_taxon) \
                or (oid is not None and self_taxon.oid == oid) \
                or (label is not None and self_taxon.label == label):
                return True
        return False

    def has_taxa(self, **kwargs):
        """
        Returns True if all taxon given by keyword argument `taxa` in self,
        or at least one Taxon object exists in self with oid or label
        for every oid given in list of oid's by keyword arg `oids`, or
        every label in list of `labels` given by keyword arg `labels`.
        """
        if "taxa" not in kwargs and "oids" not in kwargs and "labels" not in kwargs:
            raise TypeError("Need to specify `taxa`, `oids` or `labels` list.")
        taxa = set(kwargs.get("taxa",  []))
        oids = set(kwargs.get("oids", []))
        labels = set(kwargs.get("labels", []))
        taxon_oids = set([t.oid for t in self])
        taxon_labels = set([t.label for t in self])
        return taxa.issubset(self) \
            and oids.issubset(taxon_oids) \
            and labels.issubset(taxon_labels)

    def get_taxon(self, **kwargs):
        """
        Retrieves taxon object with given id OR label (if both are given, the
        first match found is returned). If taxon does not exist then
        None is returned. Also accepts `case_insensitive` as a keyword
        argument; if `label` is used as a selection criterion, and
        `case_insensitive` is True [default=False], then matching is case
        insensitive.
        """
        if "oid" not in kwargs and "label" not in kwargs:
            raise TypeError("Need to specify Taxon oid or label.")
        oid = kwargs.get("oid", None)
        label = kwargs.get("label", None)
        ci = kwargs.get("case_insensitive", False)
        if ci:
            label_lower = label.lower()
        for taxon in self:
            if (oid is not None and taxon.oid == oid) \
                    or (label is not None and taxon.label == label) \
                    or (ci and label is not None and label_lower == taxon.label.lower()):
                return taxon
        return None

    def require_taxon(self, **kwargs):
        """
        Retrieves taxon object with given id OR label (if both are given, the
        first match found is returned). If taxon does not exist and the
        `TaxonSet` is not mutable, an exception is raised.  If taxon does not
        exist and the `TaxonSet` is mutable, then a new taxon is created,
        added, and returned. Also accepts `case_insensitive` as a keyword
        argument; if `label` is used as a selection criterion, and
        `case_insensitive` is True [default=False], then matching is case
        insensitive.
        """
        taxon = self.get_taxon(**kwargs)
        if taxon is not None:
            return taxon
        label = kwargs.get("label")
        if self._is_mutable:
            taxon = Taxon(label=label, oid=kwargs.get("oid"))
            self.add(taxon)
            return taxon
        if label:
            s = '"%s" ' % label
        else:
            s = ''
        raise KeyError("Taxon %snot in TaxonSet, and cannot be created because TaxonSet is immutable" % s)

    def get_taxa(self, **kwargs):
        """
        Retrieves list of taxon objects with given id OR label (if both are
        given, the any match is included). If taxon does not
        exist then an empty list is returned.
        """
        if "oids" not in kwargs and "labels" not in kwargs:
            raise TypeError("Need to specify taxa oid's or labels")
        oids = kwargs.get("oids", [])
        labels = kwargs.get("labels", [])
        taxa = []
        for oid in oids:
            t = self.get_taxon(oid=oid)
            if t:
                taxa.append(t)
        for label in labels:
            t = self.get_taxon(label=label)
            if t:
                taxa.append(t)
        return taxa

    def add_taxon(self, taxon):
        """
        Adds taxon to self.
        """
        if not self._is_mutable:
            raise KeyError("Taxon %s:'%s' cannot be added to an immutable TaxonSet." % (taxon.oid, taxon.label))
        self.add(taxon)

    def new_taxon(self, label=None, oid=None, error_if_label_exists=False):
        "Creates and add a new `Taxon` if not already in the taxon index."
        if not self._is_mutable:
            raise KeyError("Taxon %s:'%s' cannot be added to an immutable TaxonSet" % (oid, label))
        if error_if_label_exists and self.get_taxon(label=label, taxon_required=False) is not None:
            raise KeyError("Taxon with label %s:'%s' already defined in the TaxonSet" % (oid, label))
        taxon = Taxon(label=label, oid=oid)
        self.add(taxon)
        return taxon

    def clear(self):
        "Removes all taxa from this list."
        for t in self:
            self.remove(t)

    def labels(self):
        "Convenience method to return all taxa labels."
        return [str(taxon.label) for taxon in self]

    def complement_split_bitmask(self, split):
        "Returns complement of the split bitmask."
        return (~split) & self.all_taxa_bitmask()

    def all_taxa_bitmask(self):
        "Returns mask of all taxa."
        #return pow(2, len(self)) - 1
        b = 1 << len(self)
        return b - 1

    def taxon_bitmask(self, taxon):
        """
        Returns unique bitmask of given taxon. Will raise index error if
        taxon does not exist.
        """
        try:
            return taxon.split_bitmask
        except AttributeError:
            pass
        try:
            i = self.index(taxon)
            m = 1 << i
            taxon.split_bitmask = m
            return m
        except ValueError:
            raise IndexError("Taxon with ID '%s' and label '%s' not found"
                             % (str(taxon.oid), str(taxon.label)))

    def get_taxa_bitmask(self, **kwargs):
        """
        Retrieves bitmask represent all taxa specified by keyword-specified list
        of taxon objects (`taxa=`), labels (`labels=`) or oids (`oids=`).
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
        return "%s" % textutils.int_to_bitstring(split_bitmask).rjust(len(self), "0")

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
        return textutils.split_as_newick_string(split, self, preserve_spaces=preserve_spaces, quote_underscores=quote_underscores)

    def partition(self, **kwargs):
        """
        Returns a TaxonSetPartition object, corresponding to this object
        partition according to criteria/values given in keyword arguments:

            ``membership_func``
                A function that takes a ``Taxon`` object as an argument and
                returns a a population membership identifier or flag
                (e.g., a string, an integer) .

            ``membership_attr_name``
                Name of an attribute of ``Taxon`` objects that serves as an
                identifier for subset membership.

            ``membership_dict``
                A dictionary with ``Taxon`` objects as keys and population
                membership identifier or flag as values (e.g., a string,
                an integer).

            ``membership_lists``
                A container of containers of ``Taxon`` objects, with every
                ``Taxon`` object in ``taxon_set`` represented once and only
                once in the sub-containers.
        """
        return TaxonSetPartition(self, **kwargs)

    def __str__(self):
        return "TaxonSet(%s)" % (", ".join([str(taxon) for taxon in self]))

    def __repr__(self):
        return "<%s object at %s>" % ("TaxonSet", hex(id(self)))

    def description(self, depth=1, indent=0, itemize="", output=None, **kwargs):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return ""
        output_strio = StringIO()
        if self.label is None:
            label = " (%s)" % self.oid
        else:
            label = " (%s: '%s')" % (self.oid, self.label)
        output_strio.write('%s%sTaxonSet object at %s%s'
                % (indent*' ',
                   itemize,
                   hex(id(self)),
                   label))
        if depth >= 1:
            output_strio.write(': %d Taxa' % len(self))
            if depth >= 2 and len(self) > 0:
                for i, t in enumerate(self):
                    output_strio.write('\n')
                    t.description(depth=depth-1, indent=indent+4, itemize="[%d]" % (i), output=output_strio, **kwargs)
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

class Taxon(base.IdTagged):
    """
    A taxon associated with a sequence or a node on a tree.
    """

    def __deepcopy__(self, memo):
        "Should not be copied"
        memo[id(self)] = self
        return self

    def cmp(taxon1, taxon2):
        "Compares taxon1 and taxon2 based on label."
        return cmp(str(taxon1.label), str(taxon2.label))

    cmp = staticmethod(cmp)

    def __init__(self, *args, **kwargs):
        """
        __init__ can take the kwargs needed by base.IdTagged, or the label keyword
        can be inferred from the label of an unnamed argument
        """
        if len(args) > 1:
            raise TypeError("Taxon() takes at most 1 non-keyword argument (%d given)" % len(args))
        elif len(args) == 1:
            kwargs['label'] = args[0].label
        base.IdTagged.__init__(self, **kwargs)

    def __str__(self):
        "String representation of self = taxon name."
        return "%s" % str(self.label)

    def __repr__(self):
        return "<Taxon object at %s: %s>" % (hex(id(self)), str(self))

    def description(self, depth=1, indent=0, itemize="", output=None, **kwargs):
        """
        Returns description of object, up to level `depth`.
        """
        if depth is None or depth < 0:
            return ""
        output_strio = StringIO()
        if self.label is None:
            label = "<Unnamed Taxon>"
        else:
            label = "'%s'" % self.label
        output_strio.write('%s%s Taxon object at %s (%s): %s' % (indent*' ', itemize, hex(id(self)), self.oid, label))
        s = output_strio.getvalue()
        if output is not None:
            output.write(s)
        return s

class TaxonSetPartition(TaxonSetLinked):
    """
    Manages a partition of a TaxonSet (i.e., a set of mutually-exclusive
    and exhaustive subsets of a TaxonSet).
    """

    def __init__(self, taxon_set, **kwargs):
        """
        __init__ uses one of the following keyword arguments:

            - `membership_func`
                A function that takes a ``Taxon`` object as an argument and
                returns a a population membership identifier or flag
                (e.g., a string, an integer) .
            - `membership_attr_name`
                Name of an attribute of ``Taxon`` objects that serves as an
                identifier for subset membership.
            - `membership_dict`
                A dictionary with ``Taxon`` objects as keys and population
                membership identifier or flag as values (e.g., a string,
                an integer).
            - `membership_lists`
                A container of containers of ``Taxon`` objects, with every
                ``Taxon`` object in ``taxon_set`` represented once and only
                once in the sub-containers.

        If none of these are specified, defaults to a partition consisting of
        a single subset with all the objects in ``taxon_set``.
        """
        TaxonSetLinked.__init__(self, taxon_set=taxon_set, **kwargs)
        self.subset_map = {}
        if taxon_set is not None:
            if len(kwargs) > 0:
                self.apply(**kwargs)
            else:
                ss = TaxonSet(self.taxon_set)
                self.subset_map = { self.taxon_set.label : ss}

    def subsets(self):
        """
        Return subsets of partition.
        """
        return set(self.subset_map.values())

    def __len__(self):
        """
        Number of subsets.
        """
        return len(self.subset_map)

    def __iter__(self):
        """
        Iterate over subsets.
        """
        for k, v in self.subset_map.items():
            yield v

    def __getitem__(self, label):
        """
        Get subset with specified label.
        """
        return self.subset_map[label]

    def apply(self, **kwargs):
        """
        Builds the subsets of the linked TaxonSet resulting from the
        partitioning scheme specified by one of the following keyword arguments:

            ``membership_func``
                A function that takes a ``Taxon`` object as an argument and
                returns a a population membership identifier or flag
                (e.g., a string, an integer).

            ``membership_attr_name``
                Name of an attribute of ``Taxon`` objects that serves as an
                identifier for subset membership.

            ``membership_dict``
                A dictionary with ``Taxon`` objects as keys and population
                membership identifier or flag as values (e.g., a string,
                an integer).

            ``membership_lists``
                A container of containers of ``Taxon`` objects, with every
                ``Taxon`` object in ``taxon_set`` represented once and only
                once in the sub-containers.
        """
        if "membership_func" in kwargs:
            self.apply_membership_func(kwargs["membership_func"])
        elif  "membership_attr_name" in kwargs:
            self.apply_membership_attr_name(kwargs["membership_attr_name"])
        elif  "membership_dict" in kwargs:
            self.apply_membership_dict(kwargs["membership_dict"])
        elif "membership_lists" in kwargs:
            self.apply_membership_lists(kwargs["membership_lists"])
        else:
            raise TypeError("Must specify partitioning scheme using one of: " \
                + "'membership_func', 'membership_dict', or 'membership_lists'")

    def apply_membership_func(self, mfunc):
        """
        Constructs subsets based on function ``mfunc``, which should take a
        ``Taxon`` object as an argument and return a population membership
        identifier or flag (e.g., a string, an integer).
        """
        self.subset_map = {}
        for t in self.taxon_set:
            subset_id = mfunc(t)
            if subset_id not in self.subset_map:
                self.subset_map[subset_id] = TaxonSet(label=subset_id)
            self.subset_map[subset_id].add(t)
        return self.subsets()

    def apply_membership_attr_name(self, attr_name):
        """
        Constructs subsets based on attribute ``attr_name`` of each
        ``Taxon`` object.
        """
        return self.apply_membership_func(lambda x: getattr(x, attr_name))

    def apply_membership_dict(self, mdict):
        """
        Constructs subsets based on dictionary ``mdict``, which should be
        dictionary with ``Taxon`` objects as keys and population membership
        identifier or flag as values (e.g., a string, an integer).
        """
        return self.apply_membership_func(lambda x: mdict[x])

    def apply_membership_lists(self, mlists, subset_labels=None):
        """
        Constructs subsets based on list ``mlists``, which should be an interable
        of iterables of ``Taxon`` objects, with every ``Taxon`` object in
        ``taxon_set`` represented once and only once in the sub-containers.
        """
        if subset_labels is not None:
            if len(subset_labels) != len(mlists):
                raise ValueError('Length of subset label list must equal to number of subsets')
        else:
            subset_labels = range(len(mlists))
        self.subset_map = {}
        for lidx, mlist in enumerate(mlists):
            subset_id = subset_labels[lidx]
            self.subset_map[subset_id] = TaxonSet(label=subset_id)
            for i, t in enumerate(mlist):
                self.subset_map[subset_id].add(t)
        return self.subsets()

class TaxonSetMapping(base.IdTagged):
    """
    A many-to-one mapping of ``Taxon`` objects (e.g., gene taxa to population/species taxa).
    """

    @staticmethod
    def create_contained_taxon_mapping(containing_taxon_set,
            num_contained,
            contained_taxon_label_prefix=None,
            contained_taxon_label_separator=' ',
            contained_taxon_label_func=None):
        """
        Creates and returns a TaxonSetMapping object that maps multiple
        "contained" Taxon objects (e.g., genes) to Taxon objects in
        `containing_taxon_set` (e.g., populations or species).

            `containing_taxon_set`
                A TaxonSet object that defines a Taxon for each population or
                species.

            `num_contained`
                The number of genes per population of species. The value of
                this attribute can be a scalar integer, in which case each
                species or population taxon will get the same fixed number
                of genes. Or it can be a list, in which case the list has
                to have as many elements as there are members in
                `containing_taxon_set`, and each element will specify the
                number of genes that the corresponding species or population
                Taxon will get.

            `contained_taxon_label_prefix`
                If specified, then each gene Taxon label will begin with this.
                Otherwise, each gene Taxon label will begin with the same label
                as its corresponding species/population taxon label.

            `contained_taxon_label_separator`
                String used to separate gene Taxon label prefix from its index.

            `contained_taxon_label_func`
                If specified, should be a function that takes two arguments: a
                Taxon object from `containing_taxon_set` and an integer
                specifying the contained gene index. It should return a string
                which will be used as the label for the corresponding gene
                taxon. If not None, this will bypass the
                `contained_taxon_label_prefix` and
                `contained_taxon_label_separator` arguments.
        """
        if isinstance(num_contained, int):
            _num_contained = [num_contained] * len(containing_taxon_set)
        else:
            _num_contained = num_contained
        contained_to_containing = {}
        contained_taxa = TaxonSet()
        for cidx, containing_taxon in enumerate(containing_taxon_set):
            num_new = _num_contained[cidx]
            for new_idx in range(num_new):

                if contained_taxon_label_func is not None:
                    label = contained_taxon_label_func(containing_taxon,
                            new_idx)
                else:
                    label = "%s%s%d" % (containing_taxon.label,
                            contained_taxon_label_separator,
                            new_idx+1)
                contained_taxon = Taxon(label=label)
                contained_to_containing[contained_taxon] = containing_taxon
                contained_taxa.append(contained_taxon)
        contained_to_containing_map = TaxonSetMapping(domain_taxon_set=contained_taxa,
                range_taxon_set=containing_taxon_set,
                mapping_dict=contained_to_containing)
        return contained_to_containing_map


    def __init__(self, **kwargs):
        """
        __init__ uses one of the following keyword arguments:

            - `mapping_func`
                A function that takes a ``Taxon`` object from the domain taxa
                as an argument and returns the corresponding ``Taxon`` object
                from the range taxa. If this argument is given, then a
                ``TaxonSet`` or some other container of ``Taxon`` objects needs
                to be passed using the ``taxon_set`` argument.
            - `mapping_attr_name`
                Name of an attribute of ``Taxon`` object of the domain taxa
                that references the corresponding ``Taxon`` object from the
                range taxa. If this argument is given, then a ``TaxonSet`` or
                some other container of ``Taxon`` objects needs to be passed
                using the ``taxon_set`` argument.
            - `mapping_dict`
                A dictionary with ``Taxon`` objects from the domain taxa as
                keys, and the corresponding ``Taxon`` object from the range
                taxa as values.
        """
        base.IdTagged.__init__(self, **kwargs)
        self.forward = {}
        self.reverse = {}
        if "mapping_func" in kwargs:
            if "domain_taxon_set" not in kwargs:
                raise TypeError("Must specify 'domain_taxon_set'")
            self.apply_mapping_func(kwargs["mapping_func"],
                    domain_taxon_set=kwargs["domain_taxon_set"],
                    range_taxon_set=kwargs.get("range_taxon_set", None))
        elif "mapping_attr_name" in kwargs:
            if "domain_taxon_set" not in kwargs:
                raise TypeError("Must specify 'domain_taxon_set'")
            self.apply_mapping_attr_name(kwargs["mapping_attr_name"],
                    domain_taxon_set=kwargs["domain_taxon_set"],
                    range_taxon_set=kwargs.get("range_taxon_set", None))
        elif "mapping_dict" in kwargs:
            self.apply_mapping_dict(kwargs["mapping_dict"],
                    domain_taxon_set=kwargs.get("domain_taxon_set", None),
                    range_taxon_set=kwargs.get("range_taxon_set", None))
        else:
            raise TypeError("Must specify at least one of: 'mapping_func', 'mapping_attr_name', or 'mapping_dict'")

    def __len__(self):
        """
        Number of subsets.
        """
        return len(self.forward)

    def __iter__(self):
        """
        Iterate over subsets.
        """
        for k in self.forward:
            yield k

    def items(self):
        return self.forward.items()

    def keys(self):
        return self.forward.keys()

    def __getitem__(self, taxon):
        """
        Get mapping for specified taxon.
        """
        return self.forward[taxon]

    def _get_domain_taxon_set(self):
        return self._domain_taxon_set

    def _set_domain_taxon_set(self, taxa):
        if taxa and not isinstance(taxa, TaxonSet):
            self._domain_taxon_set = TaxonSet(taxa)
        else:
            self._domain_taxon_set = taxa

    domain_taxon_set = property(_get_domain_taxon_set, _set_domain_taxon_set)

    def _get_range_taxon_set(self):
        return self._range_taxon_set

    def _set_range_taxon_set(self, taxa):
        if taxa and not isinstance(taxa, TaxonSet):
            self._range_taxon_set = TaxonSet(taxa)
        else:
            self._range_taxon_set = taxa

    range_taxon_set = property(_get_range_taxon_set, _set_range_taxon_set)

    def apply_mapping_func(self, mfunc, domain_taxon_set, range_taxon_set=None):
        """
        Constructs forward and reverse mapping dictionaries based on ``mfunc``,
        which should take a ``Taxon`` object in ``domain_taxon_set`` as an argument
        and return another ``Taxon`` object.
        """
        self.forward = {}
        self.reverse = {}
        self.domain_taxon_set = domain_taxon_set
        if range_taxon_set is None:
            self.range_taxon_set = TaxonSet()
        else:
            self.range_taxon_set = range_taxon_set
        for dt in self.domain_taxon_set:
            rt = mfunc(dt)
            if rt not in self.range_taxon_set:
                self.range_taxon_set.add(rt)
            self.forward[dt] = rt
            try:
                self.reverse[rt].add(dt)
            except KeyError:
                self.reverse[rt] = set([dt])

    def apply_mapping_attr_name(self, attr_name, domain_taxon_set, range_taxon_set=None):
        """
        Constructs mapping based on attribute ``attr_name`` of each
        ``Taxon`` object in ``domain_taxon_set``.
        """
        return self.apply_mapping_func(lambda x: getattr(x, attr_name), domain_taxon_set=domain_taxon_set, range_taxon_set=range_taxon_set)

    def apply_mapping_dict(self, mdict, domain_taxon_set=None, range_taxon_set=None):
        """
        Constructs mapping based on dictionary ``mdict``, which should have
        domain taxa as keys and range taxa as values.
        """
        if domain_taxon_set is None:
            domain_taxon_set = TaxonSet(mdict.keys())
        return self.apply_mapping_func(lambda x: mdict[x], domain_taxon_set=domain_taxon_set, range_taxon_set=range_taxon_set)

    def mesquite_association_rows(self):
        rows = []
        for rt in self.reverse:
            x1 = textutils.escape_nexus_token(rt.label)
            dt_labels = [dt.label for dt in self.reverse[rt]]
            dt_labels.sort()
            x2 = " ".join([textutils.escape_nexus_token(d) for d in dt_labels])
            rows.append("        %s / %s" % (x1, x2))
        return ",\n".join(rows)

    def write_mesquite_association_block(self, out, domain_taxon_set_title=None, range_taxon_set_title=None):
        """
        For debugging purposes ...
        """
        out.write("BEGIN TaxaAssociation;\n")
        if self.label:
            title = self.label
        else:
            title = self.oid
        out.write("    TITLE %s;\n"  % textutils.escape_nexus_token(title))
        if domain_taxon_set_title is None:
            if self.domain_taxon_set.label:
                domain_taxon_set_title = self.domain_taxon_set.label
            else:
                domain_taxon_set_title = self.domain_taxon_set.oid
        if range_taxon_set_title is None:
            if self.range_taxon_set.label:
                range_taxon_set_title = self.range_taxon_set.label
            else:
                range_taxon_set_title = self.range_taxon_set.oid
        out.write("    TAXA %s, %s;\n" % (
            textutils.escape_nexus_token(range_taxon_set_title),
            textutils.escape_nexus_token(domain_taxon_set_title)
            ))
        out.write("    ASSOCIATES\n")
        out.write(self.mesquite_association_rows() + "\n")
        out.write("    ;\n")
        out.write("END;\n")
