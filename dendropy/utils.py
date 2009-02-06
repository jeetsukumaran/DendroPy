#! /usr/bin/env python

############################################################################
##  util.py
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
This module contains various utility functions and methods.
"""

import os
import fnmatch

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

class ComplementingDict(dict):
    """
    Keys, {K_i}, are longs. Elements are accessible through both K and K ^ mask
    if mask is defined and is_complementing is True.
    """

    def __init__(self, other=None, mask=None):
        """
        Creates the local set of keys, and then initializes self with
        arguments, if any, by using the superclass methods, keeping
        the ordered keys in sync.
        """
        dict.__init__(self)
        self.mask = mask
        self.is_complementing = False
        if other:
            if isinstance(other, ComplementingDict):
                self.mask = other.mask
            if isinstance(other, dict):                
                for key, val in other.items():
                    self[key] = val
                    
    def _check_complement(self):
        return self.mask is None or not self.is_complementing
    check_complement = property(_check_complement)        

    def __getitem__(self, key):
        """
        Gets an item using a case-insensitive key.
        """
        item = self.get(key)
        if item is None:
            if not self.check_complement:
                raise KeyError(key)
            ckey = key ^ mask
            item = self.get(ckey)
            if item is None:
                raise KeyError(key)
        return item

    def __delitem__(self, key):
        """
        Remove item with specified key.
        """
        if dict.__contains__(self, key):
            dict.__delitem__(self, key)
        else:
            if not self.check_complement:            
                raise KeyError(key)        
            ckey = key ^ self.mask
            if dict.__contains__(self, ckey):
                dict.__delitem__(self, ckey)
            else:
                raise KeyError(key)

    def __contains__(self, key):
        """
        Returns true if has key.
        """
        if dict.__contains__(self, key):
            return True
        if self.check_complement:
            return False
        if dict.__contains__(self, (key ^ mask)):
            return True
        return false            

    def pop(self, key, alt_val=None):
        """
        a.pop(k[, x]):  a[k] if k in a, else x (and remove k)
        """
        if self.__contains__(key):
            val = self[key]
            del(self[key])
        else:
            return alt_val

    def complemented_keys(self):
        """
        Returns a copy of complemented keys.        
        """
        if self.mask is not None:
            return [ (k ^ self.mask) for k in self]
        else:
            raise Exception("mask not set")
    
    def index(self, key):
        """
        Return the index of key.
        Raise KeyError if not found.
        """
        count = 0
        for idx, k in enumerate(self):
            if (k == key) or (self.check_complement and (k ^ self.mask) == key):
                return count
        raise KeyError(key)

    def get(self, key, def_val=None):
        """
        Gets an item by its key, returning default if key not present.
        """
        if dict.__contains__(self, key) or not self.check_complement:
            return dict.get(self, key, def_val)
        elif self.check_complement:
            return dict.get(self, key ^ self.mask, def_val)
        return def_val  
         
    def setdefault(self, key, def_val=None):
        """
        Sets the default value to return if key not present.
        """
        dict.setdefault(self, def_val)

    def update(self, other):
        """
        updates (and overwrites) key/value pairs:
        k = { 'a':'A', 'b':'B', 'c':'C'}
        q = { 'c':'C', 'd':'D', 'f':'F'}
        k.update(q)
        {'a': 'A', 'c': 'C', 'b': 'B', 'd': 'D', 'f': 'F'}
        """
        for key, val in other.items():
            if not dict.__contains__(self, key) and not self.check_complement:
                self[key] = val
            elif self.check_complement and dict.__contains__(self, key ^ self.mask):
                self[key ^ self.mask] = val
            else:
                self[key] = val

    def fromkeys(self, iterable, value=None):
        """
        Creates a new dictionary with keys from seq and values set to value.
        """
        raise NotImplementedError
#         cd = ComplementingDict(mask=self.mask)
#         for key in iterable:
#             if not dict.__contains__(self, key) and not self.check_complement:
#                 self[key] = value
#             elif self.check_complement and not dict.__contains__(self, key ^ self.mask):
#                 self[key] = value
#         return cd
    
def pretty_print_timedelta(timedelta):
    hours, mins, secs = str(timedelta).split(":")
    return("%s hour(s), %s minute(s), %s second(s)" % (hours, mins, secs))
            
def glob_match(pathname, pattern, respect_case=False, complement=False):
    if respect_case:
        if fnmatch.fnmatchcase(pathname, pattern):
            if complement:
                return False
            else:
                return True
        else:
            if complement:
                return True
            else:
                return False
    else:
        pathname = pathname.lower()
        pattern = pattern.lower()
        if fnmatch.fnmatch(pathname, pattern):
            if complement:
                return False
            else:
                return True
        else:
            if complement:
                return True
            else:
                return False

def find_files(top, recursive=True,
               filename_filter=None,
               dirname_filter=None,
               excludes=None,
               complement=False,
               respect_case=False,
               expand_vars=True,
               include_hidden=True):
    if expand_vars:
        top = os.path.abspath(os.path.expandvars(os.path.expanduser(top)))
    if excludes == None:
        excludes = []
    filepaths = []
    if os.path.exists(top):
        for fpath in os.listdir(top):
            abspath = os.path.abspath(os.path.join(top, fpath))
            if os.path.isfile(abspath):
                if (include_hidden or not fpath.startswith('.')) \
                    and (not filename_filter or glob_match(fpath, filename_filter, respect_case, complement)):
                    to_exclude = False
                    for e in excludes:
                        if glob_match(fpath, e, respect_case):
                            to_exclude = True
                    if not to_exclude:
                        filepaths.append(abspath)
            elif os.path.isdir(abspath):
                if recursive:
                    if (include_hidden or not fpath.startswith('.')) \
                        and (not dirname_filter or (glob_match(fpath, dirname_filter, respect_case, complement))): 
                        filepaths.extend(find_files(abspath,
                                                     recursive=recursive,
                                                     filename_filter=filename_filter,
                                                     dirname_filter=dirname_filter,
                                                     excludes=excludes,
                                                     complement=complement,
                                                     respect_case=respect_case,
                                                     expand_vars=False))
    filepaths.sort()
    return filepaths

def sample_mean_var_ml(x):
    """Returns the sample mean and variance of x using the ML estimator of the 
    sample variance.
    """
    n = len(x)
    assert(n > 0)
    if n == 1:
        return x[0], 0
    s = 0.0
    ss = 0.0
    for i in x:
            s += i
            ss += i*i
    mu = s/n
    var =  (ss/n) - mu*mu
    return mu, var

def sample_mean_var_unbiased(x):
    """Returns the sample mean and variance of x using the common unbiased
    estimator of the sample variance.
    """
    n = len(x)
    assert(n > 0)
    if n == 1:
        return x[0], float('Inf')
    mean, v = sample_mean_var_ml(x)
    var = v*n/(n-1)
    return mean, var
