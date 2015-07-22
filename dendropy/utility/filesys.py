#! /usr/bin/env python

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
Various utilities in support of filesystem interaction.
"""

import fnmatch
import time
import os
import sys
import re
from threading import Event, Thread, Lock

from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)

###############################################################################
## Universal Opening

def pre_py34_open(file,
        mode='r',
        buffering=-1,
        encoding=None,
        errors=None,
        newline=None,
        closefd=True,
        opener=None):
    if encoding is not None:
        raise NotImplementedError
    if errors is not None:
        raise NotImplementedError
    if newline is None:
        if mode.startswith("r"):
            mode = mode + "U"
    else:
        raise NotImplementedError
    if closefd is not True:
        raise NotImplementedError
    if opener is not None:
        raise NotImplementedError
    return open(
            file,
            mode=mode,
            buffering=buffering)

###############################################################################
## LineReadingThread

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
        """
        __init__ processes the following arguments:

            - ``lineCallback`` is the callable that takes a string that is each line, and returns False to stop reading.  This is a way of using the class without sub-classing and overriding keep_going
            - ``stream`` is in input file-like object
            - ``filename`` can be sent instead of ``stream``, it should be the path to the file to read.
            - ``stop_event`` is an Event, that will kill the thread if it is triggered.
            - ``sleep_interval`` is the interval to sleep while waiting for a new tree to appear.

        All other arguments are passed to the Thread.__init__()
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
        """Blocks until the file ``filename`` appears or stop_event is triggered.

        Returns True if ``filename`` exists.

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
        """Blocks until the file ``filename`` appears and then returns a file
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

###############################################################################
## File Finding

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

def find_files(top,
               recursive=True,
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


# from http://snippets.dzone.com/posts/show/6313
def find_executable(executable, path=None):
    """Try to find 'executable' in the directories listed in 'path' (a
    string listing directories separated by 'os.pathsep'; defaults to
    os.environ['PATH']).  Returns the complete filename or None if not
    found
    """
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if os.name == 'os2':
        (base, ext) = os.path.splitext(executable)
        # executable files on OS/2 can have an arbitrary extension, but
        # .exe is automatically appended if no dot is present in the name
        if not ext:
            executable = executable + ".exe"
    elif sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (base, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    else:
        return None

###############################################################################
## Dealing with streams that may not have been opened with the universal
## newline option

def get_lines(stream):
    """
    Parse stream into lines, dealing with all line break conventions.
    """
    s = stream.read()
    return re.split(r'\r\n|\n|\r', s)


