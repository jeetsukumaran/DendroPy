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
##  Parts of this file (the threading code) were written by Mark T. Holder
##      as a part of CIPRES (event_consumer.py in the PIPRes python lib).
##
############################################################################

"""
This module contains various utility functions and methods.
"""
import sys
import time
import os
import copy
import subprocess
import fnmatch

from threading import Event, Thread, Lock
from dendropy import get_logger
_LOG = get_logger('dendropy.utils')

DEFAULT_SLEEP_INTERVAL=0.1
class LineReadingThread(Thread):
    """A thread that will read the input stream - designed to work with a file 
    thas is being written. Note that if the file does not end with a newline
    and the keep_going() method does not return False, then the thread will not
    terminate
    
    self.keep_going()
    is called with each line. (sub classes should override).
    
    LineReadingThread.__init__ must be called by subclasses.
    """
    def __init__(self, 
                lineCallback=None, 
                stream=None, 
                filename="", 
                stop_event=None, 
                sleep_interval=DEFAULT_SLEEP_INTERVAL, 
                store_lines=False, 
                is_file=True, 
                subproc=None, 
                *args,
                **kwargs):
        """`lineCallback` is the callable that takes a string that is each line, and 
        returns False to stop reading.  This is a way of using the class without
        sub-classing and overriding keep_going
        
        `stream` is in input file-like object
        `filename` can be sent instead of `stream`, it should be the path to the file to read.
        `stop_event` is an Event, that will kill the thread if it is triggered.
        `sleep_interval` is the interval to sleep while waiting for a new tree to appear.
        other arguments are passed to the Thread.__init__()
        """
        self.stream = stream
        self.filename = filename 
        self.lineCallback = lineCallback
        self.unfinished_line = None
        self.stop_event = stop_event
        self.sleep_interval = sleep_interval
        self.store_lines = store_lines
        if store_lines:
            self.line_list_lock = Lock()
        self.is_file = is_file
        self.lines = []
        self.subproc = subproc
        self.stop_on_subproc_exit = kwargs.get('stop_on_subproc_exit', False)
        Thread.__init__(self,group=None, target=None, name=None, 
                        args=tuple(*args), kwargs=dict(**kwargs))

    def wait_for_file_to_appear(self, filename):
        """Blocks until the file `filename` appears or stop_event is triggered.
        
        Returns True if `filename` exists.
        
        Checks for the stop_event *before* checking for the file existence. 
        (paup_wrap and raxml_wrap threads depend on this behavior).
        """
        while True:
            if (self.stop_event is not None) and self.stop_event.isSet():
                return False
            if os.path.exists(filename):
                return True
            #_LOG.debug("Waiting for %s" %filename)
            time.sleep(self.sleep_interval)

    def open_file_when_exists(self, filename):
        """Blocks until the file `filename` appears and then returns a file 
        object opened in rU mode.
        
        Returns None if the stop event is triggered.
        """
        if self.wait_for_file_to_appear(filename):
            return open(filename, "rU")
        return None


    def run(self):
        if self.stream is None:
            if not self.filename:
                _LOG.debug('"stream" and "filename" both None when LineReadingThread.run called')
                return
            self.stream = self.open_file_when_exists(self.filename)
            if self.stream is None:
                return
        self._read_stream()

    def keep_going(self, line):
        _LOG.debug("In keep_going: " + line)
        if self.store_lines:
            self.line_list_lock.acquire()
            try:
                self.lines.append(line)
                _LOG.debug("self.lines = %s" % str(self.lines))
            finally:
                self.line_list_lock.release()
        if self.lineCallback is None:
            r = True
            if self.subproc:
                _LOG.debug("subproc is not None")
                if self.subproc.returncode is None:
                    _LOG.debug("subproc.returncode is None")
                else:
                    _LOG.debug("subproc.returncode is %d" % self.subproc.returncode)
                    if self.store_lines:
                        _LOG.debug("about to call readlines")
                        line = line + "\n".join(self.stream.readlines())
                    r = False
            else:
                _LOG.debug("subproc is None")
            return r
        return self.lineCallback(line)

    def _read_stream(self):
        self.unfinished_line = ""
        while True:
            if (self.stop_event is not None) and self.stop_event.isSet():
                # when we terminate because of an event setting,
                # we pass any unfinished_line line that we have to 
                if not self.unfinished_line is None:
                    self.keep_going(self.unfinished_line)
                break
            _LOG.debug("about to readline")
            line = self.stream.readline()
            if not line:
                _LOG.debug("line is empty")
                if self.stop_on_subproc_exit:
                    self.subproc.poll()
                    if self.subproc is not None:
                        _LOG.debug("subproc is not None")
                        if self.subproc.returncode is not None:
                            _LOG.debug("subproc.returncode is %d" % self.subproc.returncode)
                            _LOG.debug("%s" % repr(self.stream))
                            l = "".join(self.stream.readlines())
                            if l:
                                self.keep_going(l)
                            break
                        else:
                            _LOG.debug("subproc.returncode is None")
                    else:
                        _LOG.debug("subproc is None")
            else:
                _LOG.debug('line is "%s"' % line)
            if not line.endswith("\n"):
                if self.unfinished_line:
                    self.unfinished_line = self.unfinished_line + line
                else:
                    self.unfinished_line = line
                time.sleep(self.sleep_interval)
            else:
                if self.unfinished_line:
                    line = self.unfinished_line + line
                self.unfinished_line = ""
                if not self.keep_going(line):
                    break
        _LOG.debug("LineReadingThread exiting")

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
        self._ordered_keys = []
        if other is not None:
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
        "Returns a shallow copy of self."
        return self.__class__(self)
        
    def iterkeys(self):
        "Returns an iterator over self's ordered keys."
        return iter(self.__ordered_keys)
    
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
    
    def __repr__(self):
        "Returns a representation of self's ordered keys."
        return "%s([%s])" \
               % (self.__class__.__name__, ', \
               '.join(["(%r, %r)" % item for item in self.iteritems()]))

    def __str__(self):
        "Returns a string representation of self."
        return "{%s}" \
               % (', '.join(["(%r, %r)" % item for item in self.iteritems()]),)

    def __getitem__(self, key):
        "Gets an item using a case-insensitive key."
        return super(OrderedCaselessDict, self).__getitem__(key.lower())

    def __setitem__(self, key, value):
        "Sets an item using a case-insensitive key,"
        if key.lower() not in self:
            self.__ordered_keys.append(str(key))
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
        return [k.lower() for k in self.__ordered_keys]
    
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
        return list(self.__ordered_keys)

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
                self.__ordered_keys.append(str(key))
            super(OrderedCaselessDict, self).__setitem__(key.lower(), val)

    def fromkeys(self, iterable, value=None):
        "Creates a new dictionary with keys from seq and values set to value."
        ocd = OrderedCaselessDict()
        for key in iterable:
            if key.lower() not in self:
                self[key] = value
        return ocd

class NormalizedBitmaskDict(dict):
    """
    Keys, {K_i}, are longs. `mask` must be provided before elements can be
    added removed from dictionary. All keys are normalized such that the right-
    most bit is '1'. That is, if the key's right-most bit is '1', it is added
    as-is, otherwise it is complemented by XOR'ing it with 'mask'.
    """
    
    def normalize(key, mask):
        if key & 1:
            return (~key) & mask
        return key
    normalize = staticmethod(normalize)        
        
    def __init__(self, other=None, mask=None):
        "Assigns mask, and then populates from `other`, if given."
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
    
def format_dict_table(rows, column_names=None, max_column_width=None, border_style=2):
    """
    Returns a string representation of a tuple of dictionaries in a
    table format. This method can read the column names directly off the
    dictionary keys, but if a tuple of these keys is provided in the
    'column_names' variable, then the order of column_names will follow
    the order of the fields/keys in that variable.
    """
    if column_names or len(rows) > 0:
        lengths = {}
        rules = {}
        if column_names:
            column_list = column_names
        else:
            try:
                column_list = rows[0].keys()
            except:
                column_list = None
        if column_list:
            # characters that make up the table rules
            border_style = int(border_style)
            #border_style = 0
            if border_style >= 1:
                vertical_rule = ' | '
                horizontal_rule = '-'
                rule_junction = '-+-'
            else:
                vertical_rule = '  '
                horizontal_rule = ''
                rule_junction = ''                
            if border_style >= 2:
                left_table_edge_rule = '| '
                right_table_edge_rule = ' |'
                left_table_edge_rule_junction = '+-'
                right_table_edge_rule_junction = '-+'
            else:
                left_table_edge_rule = ''
                right_table_edge_rule = ''
                left_table_edge_rule_junction = ''
                right_table_edge_rule_junction = ''
            
            if max_column_width:
                column_list = [c[:max_column_width] for c in column_list]
                trunc_rows = []
                for row in rows:
                    new_row = {}
                    for k in row.keys():
                        new_row[k[:max_column_width]] = str(row[k])[:max_column_width]
                    trunc_rows.append(new_row)
                rows = trunc_rows
                
            for col in column_list:
                rls = [len(str(row[col])) for row in rows]
                lengths[col] = max(rls+[len(col)])
                rules[col] = horizontal_rule*lengths[col]
                
            template_elements = ["%%(%s)-%ss" % (col, lengths[col]) for col in column_list]
            row_template = vertical_rule.join(template_elements)
            border_template = rule_junction.join(template_elements)
            full_line = left_table_edge_rule_junction + (border_template % rules) + right_table_edge_rule_junction 
            display = []
            if border_style > 0:
                display.append(full_line)
            display.append(left_table_edge_rule + (row_template % dict(zip(column_list, column_list))) + right_table_edge_rule)
            if border_style > 0:
                display.append(full_line)
            for row in rows:       
                display.append(left_table_edge_rule + (row_template % row) + right_table_edge_rule)
            if border_style > 0:        
                display.append(full_line)
            return "\n".join(display)
        else:
            return ''
    else:
        return ''

def get_git_tag(dirpath=os.path.curdir):
    p = subprocess.Popen(['git describe --tag'], 
            shell=True,
        cwd=os.path.abspath(os.path.expandvars(dirpath)),            
            stdin=subprocess.PIPE, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE)
    if p.stderr.read() == 'fatal: Not a git repository\n':
        return '[NOT A GIT REPOSITORY]'
    tags = [t for t in p.stdout.read().split("\n") if t]
    if len(tags) == 0:
        return 'UNSPECIFIED'
    else:        
        return tags[-1]
        
def get_git_commit(dirpath=os.path.curdir):
    # git rev-parse HEAD
    # git show --quiet --pretty=format:%H
    # git log -1 --pretty=format:%H
    # git rev-list -1 HEAD
    p = subprocess.Popen(['git show --quiet --pretty=format:%H'], 
        shell=True, 
        cwd=os.path.abspath(os.path.expandvars(dirpath)),
        stdin=subprocess.PIPE, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE)
    stdout = p.stdout.read()
    if len(stdout) == 0:
        return "[NOT A GIT REPOSITORY]"
    else:
        return stdout.replace("\n","")
        
def get_git_branch(dirpath=os.path.curdir):
    p = subprocess.Popen(['git branch'], 
        shell=True, 
        cwd=os.path.abspath(os.path.expandvars(dirpath)),
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE)
    branches = [t for t in p.stdout.read().split("\n") if t]
    if len(branches) == 0:
        return "[NOT A GIT REPOSITORY]"
    else:
        for b in branches:
            if b.startswith('* '):
                return b[2:]
        return "[UNIDENTIFIABLE]"

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
