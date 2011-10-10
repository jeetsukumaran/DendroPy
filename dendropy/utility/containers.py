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
Various data structures.
"""

import copy
import sys

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

        `attr_name` is the name of the attribute or property that should be
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

        `attr_name` is the name of the attribute or property that should be
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

###############################################################################
## OrderedSet

class OrderedSet(list):

    """
    Ordered collection of unique objects.
    """

    def __init__(self, *args):
        list.__init__(self, *args)

    def __setitem__(self, *args):
        raise Exception("'OrderedSet' object does not support item assignment\n")

    def __str__(self):
        return "[%s]" % ", ".join([str(i) for i in self])

    def __repr__(self):
        return "OrderedSet([%s])" % ", ".join([str(i) for i in self])

    def __hash__(self):
        return hash( (t for t in self) )

    def __cmp__(self, o):
        return list.__cmp__(self, o)

    def add(self, x):
        if x not in self:
            list.append(self, x)
            return x
        else:
            return None

    def extend(self, t):
        for x in t:
            self.append(x)

###############################################################################
## NormalizedBitmaskDict

class NormalizedBitmaskDict(dict):
    """
    Keys, {K_i}, are longs. `mask` must be provided before elements can be
    added removed from dictionary. All keys are normalized such that the right-
    most bit is '1'. That is, if the key's right-most bit is '1', it is added
    as-is, otherwise it is complemented by XOR'ing it with 'mask'.
    """

    def normalize(key, mask):
        if key & 1 == 0:
            return (~key) & mask
        return key
    normalize = staticmethod(normalize)

    def __init__(self, other=None, mask=None):
        "__init__ assigns `mask`, and then populates from `other`, if given."
        dict.__init__(self)
        self.mask = mask
        if other is not None:
            if isinstance(other, NormalizedBitmaskDict):
                self.mask = other.mask
            if isinstance(other, dict):
                for key, val in other.items():
                    self[key] = val

    def __deepcopy__(self, memo):
        o = NormalizedBitmaskDict()
        memo[id(self)] = o
        o.mask = self.mask
        for key, val in self.items():
            o[key] = copy.deepcopy(val, memo)
        return o

    def normalize_key(self, key):
        return NormalizedBitmaskDict.normalize(key, self.mask)

    def __setitem__(self, key, value):
        "Sets item with normalized key."
        dict.__setitem__(self, self.normalize_key(key), value)

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
## OrderedDict

try:
    from collections import OrderedDict
    #from collections import DUMMY_FOR_TESTING
except ImportError:

    class OrderedDict(dict):

        def __init__(self, other=None):
            """
            __init__ creates the local set of keys, and then initializes self with
            arguments, if any, by using the superclass methods, keeping
            the ordered keys in sync.
            """
            super(OrderedDict, self).__init__()
            self._ordered_keys = []
            if other is not None:
                if isinstance(other, dict):
                    for key, val in other.items():
                        if key.lower() not in self:
                            self._ordered_keys.append(key)
                        super(OrderedDict, self).__setitem__(key, val)
                else:
                    for key, val in other:
                        if key.lower() not in self:
                            self._ordered_keys.append(key)
                        super(OrderedDict, self).__setitem__(key, val)

        def copy(self):
            "Returns a shallow copy of self."
            return self.__class__(self)

        def iterkeys(self):
            "Returns an iterator over self's ordered keys."
            return iter(self._ordered_keys)

        def itervalues(self):
            "Returns an iterator over self's key, value pairs."
            for key in self.iterkeys():
                yield self[key]

        def iteritems(self):
            "Returns an iterator over self's values."
            for key in self.iterkeys():
                yield (key, self[key])

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
            return super(OrderedDict, self).__getitem__(key)

        def __setitem__(self, key, value):
            "Sets an item using a case-insensitive key,"
            if key not in self:
                self._ordered_keys.append(key)
            super(OrderedDict, self).__setitem__(key, value)

        def __delitem__(self, key):
            "Remove item with specified key."
            del(self._ordered_keys[self.index(key)])
            super(OrderedDict, \
                self).__delitem__(key)

        def __contains__(self, key):
            "Returns true if has key, regardless of case."
            return super(OrderedDict, self).__contains__(key)

        def pop(self, key, alt_val=None):
            "a.pop(k[, x]):  a[k] if k in a, else x (and remove k)"
            if key in self:
                val = self[key]
                self.__delitem__(key)
                return val
            else:
                return alt_val

        def popitem(self):
            "a.popitem()  remove and last (key, value) pair"
            key = self._ordered_keys[-1]
            item = (key, self[key])
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
                if k.lower() == key:
                    return count
                count = count + 1
            raise KeyError(key)

        def keys(self):
            "Returns a copy of the ordered list of keys."
            return list(self._ordered_keys)

        def clear(self):
            "Deletes all items from the dictionary."
            self._ordered_keys = []
            super(OrderedDict, self).clear()

        def has_key(self, key):
            "Returns true if has key, regardless of case."
            return key in self

        def get(self, key, def_val=None):
            "Gets an item by its key, returning default if key not present."
            return super(OrderedDict, self).get(key, def_val)

        def setdefault(self, key, def_val=None):
            "Sets the default value to return if key not present."
            return super(OrderedDict, self).setdefault(key, def_val)

        def update(self, other):
            """
            updates (and overwrites) key/value pairs:
            k = { 'a':'A', 'b':'B', 'c':'C'}
            q = { 'c':'C', 'd':'D', 'f':'F'}
            k.update(q)
            {'a': 'A', 'c': 'C', 'b': 'B', 'd': 'D', 'f': 'F'}
            """
            for key, val in other.items():
                if key not in self:
                    self._ordered_keys.append(key)
                super(OrderedDict, self).__setitem__(key, val)

        def fromkeys(self, iterable, value=None):
            "Creates a new dictionary with keys from seq and values set to value."
            ocd = OrderedDict()
            for key in iterable:
                if key not in self:
                    self[key] = value
            return ocd

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
    def __init__(self, source_iter, casting_func=None, filter_func=None):
        """
        `source_iter` is an iterator. `casting_func` is a function
        that takes objects returned by `source_iter` and returns other
        objects. `filter_func` is what will be applied to the SOURCE object
        to decide if it will be returned.
        """
        self.source_iter = iter(source_iter)
        self.casting_func = casting_func
        self.filter_func = filter_func
    def __iter__(self):
        "Returns self."
        return self
    def next(self):
        """
        Gets next item from the underlying iterator, and if
        filter_func returns True on it, applies casting_func to it and
        returns it.
        """
        while True:
            source_next = self.source_iter.next()
            if self.filter_func is None or self.filter_func(source_next):
                if self.casting_func is not None:
                    return self.casting_func(source_next)
                else:
                    return source_next

