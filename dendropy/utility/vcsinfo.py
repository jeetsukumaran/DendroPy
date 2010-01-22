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
Wraps up source version control system system information.
"""

import os
import sys
from cStringIO import StringIO
import subprocess
import datetime

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
            parts.append("%s/" % self.branch_name)
        if self.commit_id:
            parts.append(self.commit_id[:10])
        if self.commit_date:
            parts.append(", %s" % str(self.commit_date))
        if parts:
            return "".join(parts)
        else:
            return None

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
        self._long_description = self.get_long_description()
        self._is_available = True

    def _run_vcs(self, cmd):
        if isinstance(cmd, str):
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
            retcode = p.wait()
        except OSError, e:
            sys.stderr.write(p.stderr.read())
            raise VersionControlInfo.VcsUnavailable("Unable to execute '%s': '%s' not found" % (cmd, self.vcs_app_path))
        stdout = p.stdout.read()
        stderr = p.stderr.read()
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
        cmd = "show --quiet --pretty=format:'%H' HEAD"
        retcode, stdout, stderr = self._run_vcs(cmd)
        return stdout.replace('\n', '')

    def get_datetime(self):
        cmd = "show --quiet --pretty=format:'%at' HEAD"
        retcode, stdout, stderr = self._run_vcs(cmd)
        if stdout:
            return datetime.datetime.fromtimestamp(float(stdout.replace('\n', '')))
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
#            if "fatal: ref HEAD is not a symbolic ref" in stderr:
#                raise VersionControlInfo.NonBranchException("Current branch is invalid: detached HEAD")
            return None
        else:
            return stdout.replace('\n', '').split('/')[-1]

    def get_description(self):
        cmd = "describe --tags --long --always --abbrev=12"
        retcode, stdout, stderr = self._run_vcs(cmd)
        if retcode != 0:
#            if "fatal: No names found, cannot describe anything." in stderr:
#                raise VersionControlInfo.UntaggedException("Revision is untagged")
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
