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
Wraps up source version control system system information.
"""

import os
import subprocess
import datetime
from dendropy.utility import textprocessing
from dendropy.utility import processio

class Revision(object):
    """
    Provides (Git) version control information
    about a project.
    """

    class VcsUnavailable(Exception):
        def __init__(self, *args, **kwargs):
            Exception.__init__(self, *args, **kwargs)

    class NonRepositoryError(Exception):
        def __init__(self, *args, **kwargs):
            Exception.__init__(self, *args, **kwargs)

    class NonBranchException(Exception):
        def __init__(self, *args, **kwargs):
            Exception.__init__(self, *args, **kwargs)

    class UntaggedException(Exception):
        def __init__(self, *args, **kwargs):
            Exception.__init__(self, *args, **kwargs)

    def __init__(self, repo_path, release=None, vcs_app_path='git'):
        self.vcs_app_path = vcs_app_path
        self.repo_path = repo_path
        self.release = release
        self._commit_id = None
        self._commit_date = None
        self._commit_tag = None
        self._branch_name = None
        self._description = None
        self._long_description = None
        self._is_available = None

    def __str__(self):
        parts = []
        if self.branch_name:
            parts.append("%s-" % self.branch_name)
        if self.commit_id:
            parts.append(self.commit_id[:10])
        if self.commit_date:
            parts.append(", %s" % str(self.commit_date))
        if parts:
            return "".join(parts)
        else:
            return ""

    def __repr__(self):
        return "<%s: '%s'>" % (self.__class__.__name__, self.__str__())

    @property
    def commit_id(self):
        if self._commit_id is None:
            self.update()
        return self._commit_id

    @property
    def commit_date(self):
        if self._commit_date is None:
            self.update()
        return self._commit_date

    @property
    def commit_tag(self):
        if self._commit_tag is None:
            self.update()
        return self._commit_tag

    @property
    def branch_name(self):
        if self._branch_name is None:
            self.update()
        return self._branch_name

    @property
    def description(self):
        if self._description is None:
            self.update()
        return self._description

    @property
    def long_description(self):
        if self._long_description is None:
            self.update()
        return self._long_description

    @property
    def is_available(self):
        if self._is_available is None:
            self.update()
        return self._is_available

    def update(self, repo_path=None):
        if repo_path is not None:
            self.repo_path = repo_path
        if not self.repo_path or not self._vcs_available():
            self._commit_id = None
            self._commit_date = None
            self._commit_tag = None
            self._branch_name = None
            self._description = None
            self._long_description = None
            self._is_available = False
            return
        self._commit_id = self.get_commit_id()
        self._commit_date = self.get_datetime()
        self._commit_tag = self.get_commit_tag()
        self._branch_name = self.get_branch_name()
        self._description = self.get_description()
        self._long_description = self._build_long_description()
        self._is_available = True

    def _run_vcs(self, cmd):
        if textprocessing.is_str_type(cmd):
            cmd = self.vcs_app_path + " " + cmd
        else:
            cmd.insert(0, self.vcs_app_path)
        try:
            p = subprocess.Popen(cmd,
                shell=True,
                cwd=os.path.abspath(self.repo_path),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            stdout, stderr = processio.communicate(p)
            retcode = p.returncode
        except OSError as e:
            return -999, "", str(e)
        return retcode, stdout, stderr

    def _vcs_available(self):
        retcode, stdout, stderr = self._run_vcs("--version")
        if retcode != 0:
            return False
        retcode, stdout, stderr = self._run_vcs("status")
        if "fatal: Not a git repository" in stderr:
            return False
        return True

    def get_commit_id(self):
        # cmd = "show --quiet --pretty=format:'%H' HEAD"
        cmd = "rev-parse --short HEAD"
        retcode, stdout, stderr = self._run_vcs(cmd)
        return stdout.replace('\n', '')

    def get_datetime(self):
        cmd = "show --quiet --pretty=format:'%at' HEAD"
        retcode, stdout, stderr = self._run_vcs(cmd)
        if stdout:
            try:
                return datetime.datetime.fromtimestamp(float(stdout.replace('\n', '').replace("'", "").replace('"','')))
            except ValueError:
                return None
        else:
            return None

    def get_commit_tag(self):
        cmd = "name-rev --name-only --tags --no-undefined HEAD"
        retcode, stdout, stderr = self._run_vcs(cmd)
        if "fatal: cannot describe" in stderr:
            return None
        else:
            return stdout.strip('\n')

    def get_branch_name(self):
        # git name-rev --name-only HEAD
        cmd = "symbolic-ref HEAD"
        retcode, stdout, stderr = self._run_vcs(cmd)
        if retcode != 0:
            if "fatal: ref HEAD is not a symbolic ref" in stderr:
                return "(no branch)"
            else:
                return None
        else:
            return stdout.replace('\n', '').split('/')[-1]

    def get_description(self):
        cmd = "describe --tags --long --always --abbrev=12"
        retcode, stdout, stderr = self._run_vcs(cmd)
        if retcode != 0:
            if "fatal: No names found, cannot describe anything." in stderr:
                return "(unnamed)"
            else:
                return None
        else:
            return stdout.replace('\n', '')

    def get_long_description(self):
        parts = []
        if self.commit_id:
            parts.append(self.commit_id)
        if self.branch_name:
            parts.append("on branch '%s'" % self.branch_name)
        if self.commit_date:
            parts.append("committed on %s" % str(self.commit_date))
        if parts:
            return ", ".join(parts)
        else:
            return None

    def _build_long_description(self):
        parts = []
        if self._commit_id:
            parts.append(self._commit_id)
        if self._branch_name:
            parts.append("on branch '%s'" % self._branch_name)
        if self._commit_date:
            parts.append("committed on %s" % str(self._commit_date))
        if parts:
            return ", ".join(parts)
        else:
            return None


