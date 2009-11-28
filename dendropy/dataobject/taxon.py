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
This module provides classes and methods for managing taxa.
"""

import sys
import math
from cStringIO import StringIO
from dendropy.dataobject import base
from dendropy.utility import texttools
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
        Handles keyword arguments: `oid`, `label` or "is_mutable".
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
            for i in args[0]:
                self.add(TaxonSet._to_taxon(i))
        self._is_mutable = kwargs.get('is_mutable', True) # immutable constraints not fully implemented -- only enforced at the add_taxon stage)

    def __deepcopy__(self, memo):
        o = self.__class__(list(self), label=self.label, is_mutable=self._is_mutable)
        memo[id(self)] = o
        memo[id(self._oid)] = o._oid
        return o

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
            raise TypeError("Need to specify oid or Label.")
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
        Retrieves taxon object with given id OR label (if both are
        given, the first match found is returned). If taxon does not
        exist then None is returned.
        """
        if "oid" not in kwargs and "label" not in kwargs:
            raise TypeError("Need to specify Taxon oid or Label.")
        oid = kwargs.get("oid", None)
        label = kwargs.get("label", None)
        for taxon in self:
            if (oid is not None and taxon.oid == oid) \
                or (label is not None and taxon.label == label):
                return taxon
        return None

    def get_taxa(self, **kwargs):
        """
        Retrieves list of taxon objects with given id OR label (if both are
        given, the any match is included). If taxon does not
        exist then an empty list is returned.
        """
        if "oid" not in kwargs and "label" not in kwargs:
            raise TypeError("Need to specify Taxon oid or Label.")
        oid = kwargs.get("oid", None)
        label = kwargs.get("label", None)
        taxa = []
        for taxon in self:
            if (oid is not None and taxon.oid == oid) \
                or (label is not None and taxon.label == label):
                taxa.append(taxon)
        return taxa

    def require_taxon(self, **kwargs):
        """
        Retrieves taxon object with given id OR label (if both are
        given, the first match found is returned). If taxon does not
        exist and the `TaxonSet` is not mutable, an exception is raised.
        If taxon does not exist and the `TaxonSet` is mutable, then a
        new taxon is created, added, and returned.
        """
        taxon = self.get_taxon(**kwargs)
        if taxon is not None:
            return taxon
        elif self._is_mutable:
            taxon = Taxon(label=kwargs.get("label", None), oid=kwargs.get("oid", None))
            self.add(taxon)
            return taxon
        else:
            raise KeyError("Taxon not in TaxonSet, and cannot be created because TaxonSet is immutable.")

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
            raise KeyError("Taxon %s:'%s' cannot be added to an immutable TaxonSet." % (oid, label))
        if error_if_label_exists and self.get_taxon(label=label, taxon_required=False) is not None:
            raise KeyError("Taxon with label %s:'%s' already defined in the TaxonSet." % (oid, label))
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
        return "%s" % texttools.int_to_bitstring(split_bitmask).rjust(len(self), "0")

    def __str__(self):
        return "TaxonSet(%s)" % (", ".join([str(taxon) for taxon in self]))

    def __repr__(self):
        return "<TaxonSet object at %s>" % (hex(id(self)))

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

    def __init__(self, label=None, oid=None):
        "Initializes by calling base class."
        base.IdTagged.__init__(self, label=label, oid=oid)

    def __str__(self):
        "String representation of self = taxon name."
        return "'%s'" % str(self.label)

    def __repr__(self):
        return "<Taxon object at %s>" % (hex(id(self)))

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
