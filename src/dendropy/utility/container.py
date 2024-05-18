#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Various data structures.
"""

import collections
import copy
import csv
from dendropy.utility import deprecate
from collections.abc import MutableMapping

###############################################################################
## OrderedSet

class OrderedSet(object):

    """
    Ordered collection of unique objects with O(1) look-ups, addition, and
    deletion.
    """

    def __init__(self, *args):
        self._item_list = []
        self._item_set = set()
        if len(args) > 1:
            raise TypeError("OrderedSet expected at most 1 arguments, got 2")
        elif len(args) == 1:
            for a in args[0]:
                if a not in self._item_set:
                    self._item_set.add(a)
                    self._item_list.append(a)

    def __copy__(self, memo=None):
        o = OrderedSet(self._item_list)
        return o

    def __deepcopy__(self, memo=None):
        other = self.__class__()
        memo[id(self)] = other
        memo[id(self._item_set)] = other._item_set
        memo[id(self._item_list)] = other._item_list
        for item in self._item_list:
            c = copy.deepcopy(item, memo)
            memo[id(item)] = c
            other._item_set.add(c)
            other._item_list.append(c)
        for k in self.__dict__:
            if k in other.__dict__:
                continue
            other.__dict__[k] = copy.deepcopy(self.__dict__[k], memo)
            memo[id(self.__dict__[k])] = other.__dict__[k]
        return other

    def __len__(self):
        return len(self._item_list)

    def __getitem__(self, index):
        """
        Returns value at ``index``.
        Note takes *index* of than value as key.
        """
        return self._item_list[index]

    def __setitem__(self, index, value):
        """
        Sets value at ``index``.
        Note takes *index* of than value as key.
        """
        item = self._item_list[index]
        self._item_set.remove(item)
        self._item_set.add(value)
        self._item_list[index] = value

    def __delitem__(self, index):
        """
        Deletes value at ``index``.
        Note takes *index* of than value as key.
        """
        self._item_set.remove(self._item_list[index])
        del self._item_list[index]

    def discard(self, key):
        """
        Deletes value of ``key`` from ``self``.
        No error if no value of ``key`` is not in ``self``.
        """
        if key in self._item_set:
            self._item_set.remove(key)
            self._item_list.remove(key)

    def remove(self, key):
        """
        Deletes value of ``key`` from ``self``.
        KeyErrorif no value of ``key`` is not in ``self``.
        """
        self._item_set.remove(key)
        self._item_list.remove(key)

    def __iter__(self):
        """
        Returns iterator over values in ``self``.
        """
        return iter(self._item_list)

    def next(self):
        """
        Returns iterator over values in ``self``.
        """
        return self.__iter__()

    def __reversed__(self):
        """
        Returns ``OrderedSet`` with values in reversed order.
        """
        return OrderedSet(reversed(self._item_list))

    def __add__(self, other):
        """
        Returns ``OrderedSet`` consisting of union of values in ``self``
        and ``other``.
        """
        v = self._item_list + other._item_list
        return OrderedSet(v)

    def index(self, value):
        """
        Returns index of element with value of ``value``.
        """
        return self._item_list.index(value)

    def __contains__(self, value):
        """
        Returns |True| if ``value`` is in ``self`` or |False| otherwise.
        """
        return value in self._item_set

    def add(self, value):
        """
        Adds a new element, ``value``, to ``self`` if ``value`` is not already in
        ``self``.
        """
        if value not in self._item_set:
            self._item_set.add(value)
            self._item_list.append(value)
            return value
        else:
            return None

    def update(self, other):
        """
        Updates ``self`` with values in ``other`` for each value in ``other`` that is
        not already in ``self``.
        """
        for i in other:
            if i not in self._item_set:
                self._item_set.add(i)
                self._item_list.append(i)

    def __str__(self):
        return "[{}]".format(", ".join([str(i) for i in self._item_list]))

    def __repr__(self):
        return "{}([{}])".format(self.__class__.__name__,
            ", ".join([str(i) for i in self._item_list]))

    def __hash__(self):
        return id(self)
    #     return hash( (t for t in self._item_list) )

    def __eq__(self, o):
        return self._item_list == o._item_list

    def __lt__(self, o):
        return self._item_list < o._item_list

    def pop(self, last=True):
        """
        Removes and return value in ``self``.
        By default, removes last value.
        """
        if not self._item_set:
            raise KeyError('OrderedSet is empty')
        if last:
            key = self._item_list[-1]
        else:
            key = self._item_list[0]
        self.discard(key)
        return key

    def clear(self):
        self._item_set = set()
        self._item_list = []

###############################################################################
## NormalizedBitmaskDict

class NormalizedBitmaskDict(collections.OrderedDict):
    """
    Keys, {K_i}, are longs. ``fill_bitmask`` must be provided before elements can be
    added or removed from dictionary. All keys are normalized such that the
    least- significant bit is '0'. That is, if the key's least-significant bit
    is '0', it is added as-is, otherwise it is complemented by XOR'ing it with
    'fill_bitmask'.
    """

    @staticmethod
    def least_significant_set_bit(s):
        """
        Returns least-significant bit in integer 's' that is set.
        """
        m = s & (s - 1)
        return m ^ s

    # this is for the least-significant-bit-is-1 normalization convention
    # def normalize(key, fill_bitmask, lowest_relevant_bit):
    #     if key & lowest_relevant_bit:
    #         return key & fill_bitmask                # keep least-significant bit to 1
    #     else:
    #         return (~key) & fill_bitmask             # force least-significant bit as 1
    # normalize = staticmethod(normalize)

    # this is for the least-significant-bit-is-0 normalization convention
    @staticmethod
    def normalize(key, fill_bitmask, lowest_relevant_bit):
        if key & lowest_relevant_bit:
            return (~key) & fill_bitmask             # force least-significant bit to 0
        else:
            return key & fill_bitmask                # keep least-significant bit as 0

    def __init__(self, other=None, fill_bitmask=None):
        """
        Parameters
        ----------
        fill_bitmask : integer
            A bitmask where all possible bits that can be set to 1 are set to 1.
            When representing a taxon namespaces, with 8 taxa, for example,
            this would be 0b11111111. Incomplete leaf-sets on trees need to
            having the missing taxa bits set to 0. For example, for a tree
            missing taxa 2, 3, and 5, ``fill_bitmask`` would be 0b11101001.

        """
        if not other and not fill_bitmask:
            deprecate.dendropy_deprecation_warning(message="Deprecated since DendroPy 5: either other or fill_bitmask must be specified.")

        collections.OrderedDict.__init__(self)
        self.lowest_relevant_bit = NormalizedBitmaskDict.least_significant_set_bit(fill_bitmask)
        self.fill_bitmask = fill_bitmask
        if other is not None:
            if isinstance(other, NormalizedBitmaskDict):
                self.fill_bitmask = other.fill_bitmask
                self.lowest_relevant_bit = other.lowest_relevant_bit
            if isinstance(other, dict):
                for key, val in other.items():
                    self[key] = val

    def __deepcopy__(self, memo):
        o = NormalizedBitmaskDict(fill_bitmask=self.fill_bitmask)
        memo[id(self)] = o
        o.fill_bitmask = self.fill_bitmask
        o.lowest_relevant_bit = self.lowest_relevant_bit
        for key, val in self.items():
            o[key] = copy.deepcopy(val, memo)
        return o

    def normalize_key_and_assign_value(self, key, value):
        "*Almost* like __setitem__(), but returns value of normalized key to calling code."
        normalized_key = self.normalize_key(key)
        dict.__setitem__(self, normalized_key, value)
        return normalized_key

    def normalize_key(self, key):
        return NormalizedBitmaskDict.normalize(key, self.fill_bitmask, self.lowest_relevant_bit)

    def __setitem__(self, key, value):
        "Sets item with normalized key."
        self.normalize_key_and_assign_value(key, value)

    def __getitem__(self, key):
        "Gets an item by its key."
        key = self.normalize_key(key)
        return dict.__getitem__(self, key)

    def __delitem__(self, key):
        "Remove item with normalized key."
        key = self.normalize_key(key)
        dict.__delitem__(self, key)

    def __contains__(self, key):
        "Returns true if has normalized key."
        key = self.normalize_key(key)
        return dict.__contains__(self, key)

    def pop(self, key, alt_val=None):
        "a.pop(k[, x]):  a[k] if k in a, else x (and remove k)"
        key = self.normalize_key(key)
        return dict.pop(self, key)

    def get(self, key, def_val=None):
        "Gets an item by its key, returning default if key not present."
        key = self.normalize_key(key)
        return dict.get(self, key, def_val)

    def setdefault(self, key, def_val=None):
        "Sets the default value to return if key not present."
        dict.setdefault(self, self.normalize_key(key), def_val)

    def update(self, other):
        """
        updates (and overwrites) key/value pairs:
        k = { 'a':'A', 'b':'B', 'c':'C'}
        q = { 'c':'C', 'd':'D', 'f':'F'}
        k.update(q)
        {'a': 'A', 'c': 'C', 'b': 'B', 'd': 'D', 'f': 'F'}
        """
        for key, val in other.items():
            dict.__setitem__(self, self.normalize_key(key), val)

    def fromkeys(self, iterable, value=None):
        "Creates a new dictionary with keys from seq and values set to value."
        raise NotImplementedError

###############################################################################
# CaseInsensitiveDict
#
# From:
#        https://github.com/kennethreitz/requests
#
# Copyright 2014 Kenneth Reitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

class CaseInsensitiveDict(MutableMapping):
    """
    A case-insensitive ``dict``-like object.

    Implements all methods and operations of
    ``collections.MutableMapping`` as well as dict's ``copy``. Also
    provides ``lower_items``.

    All keys are expected to be strings. The structure remembers the
    case of the last key to be set, and ``iter(instance)``,
    ``keys()``, ``items()``, ``iterkeys()``, and ``iteritems()``
    will contain case-sensitive keys. However, querying and contains
    testing is case insensitive:

        cid = CaseInsensitiveDict()
        cid['Accept'] = 'application/json'
        cid['aCCEPT'] == 'application/json'  # True
        list(cid) == ['Accept']  # True

    For example, ``headers['content-encoding']`` will return the
    value of a ``'Content-Encoding'`` response header, regardless
    of how the header name was originally stored.

    If the constructor, ``.update``, or equality comparison
    operations are given keys that have equal ``.lower()``s, the
    behavior is undefined.

    """
    def __init__(self, data=None, **kwargs):
        self._store = dict()
        if data is None:
            data = {}
        self.update(data, **kwargs)

    def __setitem__(self, key, value):
        # Use the lowercased key for lookups, but store the actual
        # key alongside the value.
        self._store[key.lower()] = (key, value)

    def __getitem__(self, key):
        return self._store[key.lower()][1]

    def __delitem__(self, key):
        del self._store[key.lower()]

    def __iter__(self):
        return (casedkey for casedkey, mappedvalue in self._store.values())

    def __len__(self):
        return len(self._store)

    def lower_items(self):
        """Like iteritems(), but with all lowercase keys."""
        return (
            (lowerkey, keyval[1])
            for (lowerkey, keyval)
            in self._store.items()
        )

    def __eq__(self, other):
        if isinstance(other, collections.Mapping):
            other = CaseInsensitiveDict(other)
        else:
            return NotImplementedError
        # Compare insensitively
        return dict(self.lower_items()) == dict(other.lower_items())

    # Copy is required
    def copy(self):
        return CaseInsensitiveDict(self._store.values())

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, dict(self.items()))

# CaseInsensitiveDict
###############################################################################

###############################################################################
## OrderedCaselessDict

class OrderedCaselessDict(dict):
    """
    Inherits from dict. Maintains two sets of keys: the first the keys
    belonging to dict, which actually accesses the container
    items. This is always cast to lower() whenever it is called, thus
    ensuring that keys are always of the same case. The second set of
    keys it maintains locally in an list, thus maintaining the order
    in which they were added. The second set of keys is not cast to
    lower(), which means that client code can always recover the
    original 'canonical' casing of the keys.

    ONLY TAKES STRING KEYS!
    """

    def __init__(self, other=None):
        """
        __init__ creates the local set of keys, and then initializes self with
        arguments, if any, by using the superclass methods, keeping
        the ordered keys in sync.
        """
        super(OrderedCaselessDict, self).__init__()
        self._ordered_keys = []
        if other is not None:
            if isinstance(other, dict):
                for key, val in other.items():
                    if key.lower() not in self:
                        self._ordered_keys.append(str(key))
                    super(OrderedCaselessDict, \
                          self).__setitem__(key.lower(), val)
            else:
                for key, val in other:
                    if key.lower() not in self:
                        self._ordered_keys.append(str(key))
                    super(OrderedCaselessDict, \
                          self).__setitem__(key.lower(), val)

    def __deepcopy__(self, memo):
        o = self.__class__()
        memo[id(self)] = o
        for key, val in self.items():
            o[key] = copy.deepcopy(val, memo)
        return o

    def copy(self):
        "Returns a shallow copy of self."
        return self.__class__(self)

    def iterkeys(self):
        "Returns an iterator over self's ordered keys."
        return iter(self._ordered_keys)

    def itervalues(self):
        "Returns an iterator over self's key, value pairs."
        for key in self.iterkeys():
            yield self[key.lower()]

    def iteritems(self):
        "Returns an iterator over self's values."
        for key in self.iterkeys():
            yield (key, self[key.lower()])

    def items(self):
        "Returns key, value pairs in key-order."
        return [(key, self[key]) for key in self.iterkeys()]

    def values(self):
        "Returns list of key, value pairs."
        return [v for v in self.itervalues()]

    def __iter__(self):
        "Returns an iterator over self's ordered keys."
        return self.iterkeys()

    def __getitem__(self, key):
        "Gets an item using a case-insensitive key."
        return super(OrderedCaselessDict, self).__getitem__(key.lower())

    def __setitem__(self, key, value):
        "Sets an item using a case-insensitive key,"
        if key.lower() not in self:
            self._ordered_keys.append(str(key))
        super(OrderedCaselessDict, self).__setitem__(key.lower(), value)

    def __delitem__(self, key):
        "Remove item with specified key."
        del(self._ordered_keys[self.index(key)])
        super(OrderedCaselessDict, \
              self).__delitem__(key.lower())

    def __contains__(self, key):
        "Returns true if has key, regardless of case."
        return super(OrderedCaselessDict, self).__contains__(key.lower())

    def pop(self, key, alt_val=None):
        "a.pop(k[, x]):  a[k] if k in a, else x (and remove k)"
        if key.lower() in self:
            val = self[key]
            self.__delitem__(key.lower())
            return val
        else:
            return alt_val

    def popitem(self):
        "a.popitem()  remove and last (key, value) pair"
        key = self._ordered_keys[-1]
        item = (key, self[key.lower()])
        self.__delitem__(key)
        return item

    def caseless_keys(self):
        "Returns a copy of the ordered list of keys."
        return [k.lower() for k in self._ordered_keys]

    def index(self, key):
        """
        Return the index of (caseless) key.
        Raise KeyError if not found.
        """
        count = 0
        for k in self._ordered_keys:
            if k.lower() == key.lower():
                return count
            count = count + 1
        raise KeyError(key)

    def keys(self):
        "Returns a copy of the ordered list of keys."
        return list(self._ordered_keys)

    def clear(self):
        "Deletes all items from the dictionary."
        self._ordered_keys = []
        super(OrderedCaselessDict, self).clear()

    def has_key(self, key):
        "Returns true if has key, regardless of case."
        return key.lower() in self

    def get(self, key, def_val=None):
        "Gets an item by its key, returning default if key not present."
        return super(OrderedCaselessDict, self).get(key.lower(), def_val)

    def setdefault(self, key, def_val=None):
        "Sets the default value to return if key not present."
        return super(OrderedCaselessDict, self).setdefault(key.lower(), def_val)

    def update(self, other):
        """
        updates (and overwrites) key/value pairs:
        k = { 'a':'A', 'b':'B', 'c':'C'}
        q = { 'c':'C', 'd':'D', 'f':'F'}
        k.update(q)
        {'a': 'A', 'c': 'C', 'b': 'B', 'd': 'D', 'f': 'F'}
        """
        for key, val in other.items():
            if key.lower() not in self:
                self._ordered_keys.append(str(key))
            super(OrderedCaselessDict, self).__setitem__(key.lower(), val)

    def fromkeys(self, iterable, value=None):
        "Creates a new dictionary with keys from seq and values set to value."
        ocd = OrderedCaselessDict()
        for key in iterable:
            if key.lower() not in self:
                self[key] = value
        return ocd

###############################################################################
## FrozenOrderedDict

class FrozenOrderedDict(collections.OrderedDict):

    class ImmutableTypeError(TypeError):
        def __init__(self, *args, **kwargs):
            super(FrozenOrderedDict.ImmutableTypeError, self).__init__(*args, **kwargs)

    def __init__(self, *args, **kwargs):
        self._is_frozen = False
        super(FrozenOrderedDict, self).__init__(*args, **kwargs)
        self._is_frozen = True

    def __setitem__(self, key, value):
        if self._is_frozen:
            raise FrozenOrderedDict.ImmutableTypeError("{} is immutable".format(self.__class__.__name__))
        else:
            super(FrozenOrderedDict, self).__setitem__(key, value)

    def __delitem__(self, key):
        if self._is_frozen:
            raise FrozenOrderedDict.ImmutableTypeError("{} is immutable".format(self.__class__.__name__))
        else:
            super(FrozenOrderedDict, self).__delitem__(key)

    def pop(self, key, alt_val=None):
        if self._is_frozen:
            raise FrozenOrderedDict.ImmutableTypeError("{} is immutable".format(self.__class__.__name__))
        else:
            super(FrozenOrderedDict, self).pop(key, alt_val)

    def popitem(self):
        if self._is_frozen:
            raise FrozenOrderedDict.ImmutableTypeError("{} is immutable".format(self.__class__.__name__))
        else:
            super(FrozenOrderedDict, self).popitem()

    def clear(self):
        if self._is_frozen:
            raise FrozenOrderedDict.ImmutableTypeError("{} is immutable".format(self.__class__.__name__))
        else:
            super(FrozenOrderedDict, self).clear()

    def update(self, other):
        if self._is_frozen:
            raise FrozenOrderedDict.ImmutableTypeError("{} is immutable".format(self.__class__.__name__))
        else:
            super(FrozenOrderedDict, self).update(other)

    def fromkeys(self, iterable, value=None):
        if self._is_frozen:
            raise FrozenOrderedDict.ImmutableTypeError("{} is immutable".format(self.__class__.__name__))
        else:
            super(FrozenOrderedDict, self).fromkeys(iterable, value)

    def __deepcopy__(self, memo):
        temp = FrozenOrderedDict()
        temp._is_frozen = False
        for k in self:
            k2 = copy.deepcopy(k, memo)
            memo[id(k)] = k2
            v = self[k]
            v2 = copy.deepcopy(v, memo)
            memo[id(v)] = v2
            temp[k2] = v2
        temp._is_frozen = True
        return temp

    def __copy__(self):
        temp = FrozenOrderedDict()
        temp._is_frozen = False
        for k in self:
            temp[k] = self[k]
        temp._is_frozen = True
        return temp

##############################################################################
## DataTable

class DataTable(object):

    @classmethod
    def from_csv(cls,
            src,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            is_trim_space=True,
            default_data_type=None,
            column_data_types=None,
            label_transform_fn=None,
            **csv_reader_kwargs
            ):
        """
        Returns table from a full-configured csv.Reader instance.

        Parameters
        -----------
        csv_reader : a csv.Reader instance
            The source of the data.
        is_first_row_column_names : bool
            If True, then the first row is interpreted as column names;
            otherwise, treated as data row.
        is_first_column_row_names : bool
            If True, then the first column is interpreted as row names;
            otherwise, treated as data column.
        is_trim_space : bool
            If True, then will strip space from both sides of all tokens before
            processing. Note that only spaces will be trimmed, not tabs.
        default_data_type : type or function object
            Any callable that, when passed a value, returns the coerced-to-type
            equivalent for the data.
        column_data_types : dict
            A dictionary where the key are the column names and the
            corresponding value the type (see ``default_data_type`` for
            description).

        Returns
        -------
        t: a DataTable instance
            Returns table from a full-configured csv.Reader instance.

        """
        if isinstance(src, str):
            with open(src, "r") as fsrc:
                return cls._from_csv_file(
                    src=fsrc,
                    is_first_row_column_names=is_first_row_column_names,
                    is_first_column_row_names=is_first_column_row_names,
                    is_trim_space=is_trim_space,
                    default_data_type=default_data_type,
                    column_data_types=column_data_types,
                    label_transform_fn=label_transform_fn,
                    **csv_reader_kwargs)
        else:
            return cls._from_csv_file(
                src=src,
                is_first_row_column_names=is_first_row_column_names,
                is_first_column_row_names=is_first_column_row_names,
                is_trim_space=is_trim_space,
                default_data_type=default_data_type,
                column_data_types=column_data_types,
                label_transform_fn=label_transform_fn,
                **csv_reader_kwargs)

    @classmethod
    def _from_csv_file(cls,
            src,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            is_trim_space=True,
            default_data_type=None,
            column_data_types=None,
            label_transform_fn=None,
            **csv_reader_kwargs
            ):
        ncols = None
        data_table = cls()
        if is_first_row_column_names:
            first_data_row_offset = 1
        else:
            first_data_row_offset = 0
        if is_first_column_row_names:
            first_data_column_offset = 1
        else:
            first_data_column_offset = 0
        if column_data_types is None:
            column_data_types = {}
        if label_transform_fn is None:
            label_transform_fn = lambda x: x
        csv_reader = csv.reader(src, **csv_reader_kwargs)
        for row_idx, row in enumerate(csv_reader):
            if ncols is None:
                ncols = len(row)
            else:
                if len(row) == 1 and row[0].strip() == "" and is_trim_space: # blank row
                    continue
                assert ncols == len(row), "{} != {}".format(ncols, len(row))
            for cell_idx, cell in enumerate(row):
                if is_trim_space:
                    cell = cell.strip(" ")
                if row_idx == 0 and is_first_row_column_names:
                    if cell_idx == 0 and is_first_column_row_names:
                        continue
                    data_table.add_column(
                            column_name=label_transform_fn(cell),
                            data_type=column_data_types.get(cell, default_data_type))
                elif cell_idx == 0 and is_first_column_row_names:
                    data_table.add_row(label_transform_fn(cell))
                else:
                    if row_idx == 0 and not is_first_row_column_names:
                        data_table.add_column(data_type=column_data_types.get(cell, default_data_type))
                    if cell_idx == 0 and not is_first_column_row_names:
                        data_table.add_row()
                    effective_row_idx = row_idx - first_data_row_offset
                    assert effective_row_idx < len(data_table._row_names), "{}: {}".format(effective_row_idx, data_table._row_names)
                    effective_column_idx = cell_idx - first_data_column_offset
                    assert effective_column_idx < len(data_table._column_names), "{}: {}".format(effective_column_idx, data_table._column_names)
                    data_table[effective_row_idx, effective_column_idx] = cell
        return data_table

    def __init__(self):
        self._row_names = []
        self._row_name_set = set()
        self._column_names = []
        self._column_data_types = {}
        self._column_name_set = set()
        self._data = {}

    def add_column(self, column_name=None, pos=None, data_type=None):
        column_name = self._validate_new_column_name(column_name)
        assert column_name not in self._column_name_set
        if pos is None:
            pos = len(self._column_names)
        self._column_names.insert(pos, column_name)
        self._column_name_set.add(column_name)
        self._column_data_types[column_name] = data_type

    def add_row(self, row_name=None, pos=None):
        row_name = self._validate_new_row_name(row_name)
        assert row_name not in self._row_name_set
        if pos is None:
            pos = len(self._row_names)
        self._row_names.insert(pos, row_name)
        self._row_name_set.add(row_name)

    def __getitem__(self, key):
        row_name = self._dereference_key(
                key=key[0],
                name_list=self._row_names,
                name_set=self._row_name_set)
        column_name = self._dereference_key(
                key=key[1],
                name_list=self._column_names,
                name_set=self._column_name_set)
        if row_name not in self._data:
            return None
        if column_name not in self._data[row_name]:
            return None
        return self._data[row_name][column_name]

    def __setitem__(self, key, value):
        row_name = self._dereference_key(
                key=key[0],
                name_list=self._row_names,
                name_set=self._row_name_set)
        column_name = self._dereference_key(
                key=key[1],
                name_list=self._column_names,
                name_set=self._column_name_set)
        if row_name not in self._data:
            self._data[row_name] = {}
        if self._column_data_types[column_name] is not None:
            value = self._column_data_types[column_name](value)
        self._data[row_name][column_name] = value

    def row_name_iter(self):
        for row_name in self._row_names:
            yield row_name

    def column_name_iter(self):
        for column_name in self._column_names:
            yield column_name

    def row_value_iter(self, column_name):
        column_name = self._dereference_key(
                key=column_name,
                name_list=self._column_names,
                name_set=self._column_name_set)
        for row_name in self._row_names:
            yield self[row_name, column_name]

    def column_value_iter(self, row_name):
        row_name = self._dereference_key(
                key=row_name,
                name_list=self._row_names,
                name_set=self._row_name_set)
        for column_name in self._column_names:
            yield self[row_name, column_name]

    def write_csv(self,
            out,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            missing_data_value="NA",
            **csv_writer_kwargs
            ):
        if isinstance(out, str):
            dest = open(out, "w")
        else:
            dest = out
        if "delimiter" not in csv_writer_kwargs:
            csv_writer_kwargs["delimiter"] = ","
        csv_writer = csv.writer(dest, csv_writer_kwargs)
        if is_first_row_column_names:
            header = []
            if is_first_column_row_names:
                header.append("")
            header.extend(self._column_names)
            csv_writer.writerow(header)
        for row_name in self._row_names:
            row = []
            if is_first_column_row_names:
                row.append(row_name)
            for column_name in self._column_names:
                v = self[row_name, column_name]
                if v is None:
                    v = missing_data_value
                row.append(v)
            csv_writer.writerow(row)

    def _validate_new_column_name(self, column_name=None):
        if column_name is None:
            column_name = "V{}".format(len(self._column_names))
        else:
            column_name = str(column_name)
        return column_name

    def _validate_new_row_name(self, row_name=None):
        if row_name is None:
            row_name = "V{}".format(len(self._row_names))
        else:
            row_name = str(row_name)
        return row_name

    def _dereference_key(self, key, name_list, name_set):
        if isinstance(key, int):
            return name_list[key]
        else:
            if key not in name_set:
                raise KeyError(key)
            else:
                return key

    def num_rows(self):
        return len(self._row_names)

    def num_columns(self):
        return len(self._column_names)

    def __len__(self):
        return self.num_rows()

###############################################################################
## Generic Container Interace (for reference)

class _ContainerInterface(list):

    """
    The following methods can be defined to implement container objects. Containers usually are
    sequences (such as lists or tuples) or mappings (like dictionaries), but can represent other
    containers as well. The first set of methods is used either to emulate a sequence or to
    emulate a mapping; the difference is that for a sequence, the allowable keys should be the
    integers k for which 0 <= k < N where N is the length of the sequence, or slice objects,
    which define a range of items. (For backwards compatibility, the method __getslice__() (see
    below) can also be defined to handle simple, but not extended slices.) It is also recommended
    that mappings provide the methods keys(), values(), items(), has_key(), get(), clear(),
    setdefault(), iterkeys(), itervalues(), iteritems(), pop(), popitem(), copy(), and update()
    behaving similar to those for Python's standard dictionary objects. The UserDict module
    provides a DictMixin class to help create those methods from a base set of __getitem__(),
    __setitem__(), __delitem__(), and keys(). Mutable sequences should provide methods append(),
    count(), index(), extend(), insert(), pop(), remove(), reverse() and sort(), like Python
    standard list objects. Finally, sequence types should implement addition (meaning
    concatenation) and multiplication (meaning repetition) by defining the methods __add__(),
    __radd__(), __iadd__(), __mul__(), __rmul__() and __imul__() described below; they should not
    define __coerce__() or other numerical operators. It is recommended that both mappings and
    sequences implement the __contains__() method to allow efficient use of the in operator; for
    mappings, in should be equivalent of has_key(); for sequences, it should search through the
    values. It is further recommended that both mappings and sequences implement the __iter__()
    method to allow efficient iteration through the container; for mappings, __iter__() should be
    the same as iterkeys(); for sequences, it should iterate through the values.
    """

    def __len__(self):
        """
        Called to implement the built-in function len(). Should return the length of the
        object, an integer >= 0. Also, an object that doesn't define a __nonzero__() method
        and whose __len__() method returns zero is considered to be false in a Boolean
        context.
        """
        pass

    def __getitem__(self, key):
        """
        Called to implement evaluation of self[key]. For sequence types, the accepted keys
        should be integers and slice objects. Note that the special interpretation of
        negative indices (if the class wishes to emulate a sequence type) is up to the
        __getitem__() method. If key is of an inappropriate type, TypeError may be raised; if
        of a value outside the set of indices for the sequence (after any special
        interpretation of negative values), IndexError should be raised. For mapping types,
        if key is missing (not in the container), KeyError should be raised.

        Note

        for loops expect that an IndexError will be raised for illegal indices to allow proper
        detection of the end of the sequence.
        """
        pass

    def __setitem__(self, key, value):
        """
        Called to implement assignment to self[key]. Same note as for __getitem__(). This should only be implemented for mappings if the objects support changes to the values for keys, or if new keys can be added, or for sequences if elements can be replaced. The same exceptions should be raised for improper key values as for the __getitem__() method.
        """
        pass

    def __delitem__(self, key):
        """
        Called to implement deletion of self[key]. Same note as for __getitem__(). This should only be implemented for mappings if the objects support removal of keys, or for sequences if elements can be removed from the sequence. The same exceptions should be raised for improper key values as for the __getitem__() method.
        """
        pass

    def __iter__(self):
        """
        This method is called when an iterator is required for a container. This method
        should return a new iterator object that can iterate over all the objects in the
        container. For mappings, it should iterate over the keys of the container, and should
        also be made available as the method iterkeys().

        Iterator objects also need to implement this method; they are required to return
        themselves. For more information on iterator objects, see Iterator Types.
        """
        pass

    def __reversed__(self):
        """
        Called (if present) by the reversed() builtin to implement reverse iteration. It
        should return a new iterator object that iterates over all the objects in the
        container in reverse order.

        If the __reversed__() method is not provided, the reversed() builtin will fall back
        to using the sequence protocol (__len__() and __getitem__()). Objects that support
        the sequence protocol should only provide __reversed__() if they can provide an
        implementation that is more efficient than the one provided by reversed().

        New in version 2.6.
        """
        pass

    def __contains__(self, item):
        """
        The membership test operators (in and not in) are normally implemented as an
        iteration through a sequence. However, container objects can supply the following
        special method with a more efficient implementation, which also does not require the
        object be a sequence.

        Called to implement membership test operators. Should return true if item is in self,
        false otherwise. For mapping objects, this should consider the keys of the mapping
        rather than the values or the key-item pairs.
        """
        pass

###############################################################################
## RecastingIterator

class RecastingIterator(object):
    """
    Given an iterator I_X that returns objects of type X {x1, x2, x3,
    ... etc.}, and a function F(X), that takes objects of type X as an
    argument and returns objects of type Y, F(X) = Y, this class will
    act as an iterator that returns objects of type Y, I_Y, given an
    iterator on X. The 'function' given can be a class if the class's
    constructor takes a single argument of type X.
    """
    def __init__(self, source_iter, casting_fn=None, filter_fn=None):
        """
        ``source_iter`` is an iterator. ``casting_fn`` is a function
        that takes objects returned by ``source_iter`` and returns other
        objects. ``filter_fn`` is what will be applied to the SOURCE object
        to decide if it will be returned.
        """
        self.source_iter = iter(source_iter)
        self.casting_fn = casting_fn
        self.filter_fn = filter_fn
    def __iter__(self):
        "Returns self."
        return self
    def next(self):
        """
        Gets next item from the underlying iterator, and if
        filter_fn returns True on it, applies casting_fn to it and
        returns it.
        """
        while True:
            source_next = self.source_iter.next()
            if self.filter_fn is None or self.filter_fn(source_next):
                if self.casting_fn is not None:
                    return self.casting_fn(source_next)
                else:
                    return source_next

###############################################################################
## ItemAttributeProviderList

class ItemAttributeProxyList(list):
    """
    This list will return the attribute of its elements instead of the
    elements themselves.
    """

    def __init__(self, attr_name, *args):
        """
        __init__ calls the list.__init__ with all unnamed args.

        ``attr_name`` is the name of the attribute or property that should be
        returned.
        """
        self.bound_attr_name = attr_name
        list.__init__(self, *args)

    def __getitem__(self, *args):
        item = list.__getitem__(self, *args)
#         assert hasattr(item, self.bound_attr_name)
        return getattr(item, self.bound_attr_name)

    def __iter__(self):
        """
        Iterates over all elements in self, returning their `<bound_attr_name>`
        attribute.
        """
        for item in list.__iter__(self):
            yield getattr(item, self.bound_attr_name)

    def aggregate(self):
        """
        Returns a shallow-copy list of the `<bound_attr_name>` attribute of self's
        members.
        """
        return [getattr(item, self.bound_attr_name) for item in self]

###############################################################################
## ItemSublistProxyList

class ItemSublistProxyList(list):
    """
    This list will return a the elements of the bound list attribute of its own
    elements instead of the element themselves.
    """

    def __init__(self, attr_name, *args):
        """
       __init__ calls the list.__init__ with all unnamed args.

        ``attr_name`` is the name of the attribute or property that should be
        returned.
        """
        self.bound_attr_name = attr_name
        list.__init__(self, *args)

    def __str__(self):
        return str([str(i) for i in self])

    def __getitem__(self, *args):
        req_idx = int(args[0])
        for idx, i in enumerate(self):
            if idx == req_idx:
                return i
        raise IndexError("list index out of range: %d (max=%d)" % (req_idx, idx))

    def __len__(self):
        return len([i for i in self])

    def __iter__(self):
        """
        Iterates over all elements in self, returning their `<bound_attr_name>`
        attribute.
        """
        for item in list.__iter__(self):
            for sublist_item in getattr(item, self.bound_attr_name):
                yield sublist_item

    def aggregate(self):
        """
        Returns a shallow-copy list of the `<bound_attr_name>` attribute of self's
        members.
        """
        return [i for i in self]

