#! /usr/bin/env python

############################################################################
##  util.py
##
##  Part of the DendroPy phylogenetic library.
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
This module contains various utility functions and methods
"""

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
        """
        Returns self.
        """
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
        Creates the local set of keys, and then initializes self with
        arguments, if any, by using the superclass methods, keeping
        the ordered keys in sync.
        """
        super(OrderedCaselessDict, self).__init__()
        self.__ordered_keys = []
        if other:
            if isinstance(other, dict):
                for key, val in other.items():
                    if key.lower() not in self:
                        self.__ordered_keys.append(str(key))
                    super(OrderedCaselessDict, \
                          self).__setitem__(key.lower(), val)
            else:
                for key, val in other:
                    if key.lower() not in self:
                        self.__ordered_keys.append(str(key))
                    super(OrderedCaselessDict, \
                          self).__setitem__(key.lower(), val)
    
    def copy(self):
        """
        Returns a shallow copy of self.
        """
        return self.__class__(self)
        
    def iterkeys(self):
        """
        Returns an iterator over self's ordered keys.
        """
        return iter(self.__ordered_keys)
    
    def itervalues(self):
        """
        Returns an iterator over self's key, value pairs.
        """
        for key in self.iterkeys():
            yield self[key.lower()]
    
    def iteritems(self):
        """
        Returns an iterator over self's values.
        """
        for key in self.iterkeys():
            yield (key, self[key.lower()])

    def items(self):
        """
        Returns key, value pairs in key-order.
        """
        return [(key, self[key]) for key in self.iterkeys()]

    def values(self):
        """
        Returns list of key, value pairs.
        """
        return [v for v in self.itervalues()]
    
    def __iter__(self):
        """
        Returns an iterator over self's ordered keys.
        """
        return self.iterkeys()
    
    def __repr__(self):
        """
        Returns a representation of self's ordered keys.
        """
        return "%s([%s])" \
               % (self.__class__.__name__, ', \
               '.join(["(%r, %r)" % item for item in self.iteritems()]))

    def __str__(self):
        """
        Returns a string representation of self.
        """
        return "{%s}" \
               % (', '.join(["(%r, %r)" % item for item in self.iteritems()]),)

    def __getitem__(self, key):
        """
        Gets an item using a case-insensitive key.
        """
        return super(OrderedCaselessDict, self).__getitem__(key.lower())

    def __setitem__(self, key, value):
        """
        Sets an item using a case-insensitive key,
        """
        if key.lower() not in self:
            self.__ordered_keys.append(str(key))
        super(OrderedCaselessDict, self).__setitem__(key.lower(), value)

    def __delitem__(self, key):
        """
        Remove item with specified key.
        """
        del(self.__ordered_keys[self.index(key)])
        super(OrderedCaselessDict, \
              self).__delitem__(key.lower())                

    def __contains__(self, key):
        """
        Returns true if has key, regardless of case.
        """
        return super(OrderedCaselessDict, self).__contains__(key.lower())

    def pop(self, key, alt_val=None):
        """
        a.pop(k[, x]):  a[k] if k in a, else x (and remove k)
        """
        if key.lower() in self:
            val = self[key]
            self.__delitem__(key.lower())
            return val
        else:
            return alt_val
        
    def popitem(self):
        """
        a.popitem()  remove and last (key, value) pair
        """
        key = self.__ordered_keys[-1]
        item = (key, self[key.lower()])
        self.__delitem__(key)
        return item

    def caseless_keys(self):
        """
        Returns a copy of the ordered list of keys.
        """
        return [k.lower() for k in self.__ordered_keys]
    
    def index(self, key):
        """
        Return the index of (caseless) key.
        Raise KeyError if not found.
        """
        count = 0
        for k in self.__ordered_keys:
            if k.lower() == key.lower():
                return count
            count = count + 1
        raise KeyError(key)

    def keys(self):
        """
        Returns a copy of the ordered list of keys.
        """
        return list(self.__ordered_keys)

    def clear(self):
        """
        Deletes all items from the dictionary.
        """
        self.__ordered_keys = []
        super(OrderedCaselessDict, self).clear()

    def has_key(self, key):
        """
        Returns true if has key, regardless of case.
        """
        return key.lower() in self

    def get(self, key, def_val=None):
        """
        Gets an item by its key, returning default if key not present.
        """
        return super(OrderedCaselessDict, self).get(key.lower(), def_val)

    def setdefault(self, key, def_val=None):
        """
        Sets the default value to return if key not present.
        """
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
                self.__ordered_keys.append(str(key))
            super(OrderedCaselessDict, self).__setitem__(key.lower(), val)

    def fromkeys(self, iterable, value=None):
        """
        Creates a new dictionary with keys from seq and values set to value.
        """
        ocd = OrderedCaselessDict()
        for key in iterable:
            if key.lower() not in self:
                self[key] = value
        return ocd
    
