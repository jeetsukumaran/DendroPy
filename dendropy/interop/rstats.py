#! /usr/bin/env python

import os
import sys
import subprocess
from dendropy.utility import metavar
from dendropy.utility import libexec
from dendropy.utility import processio

RSCRIPT_EXECUTABLE = os.environ.get(metavar.DENDROPY_RSCRIPT_PATH_ENVAR, "Rscript")
if RSCRIPT_EXECUTABLE == "NONE":
    DENDROPY_RSCRIPT_INTEROPERABILITY = False
else:
    DENDROPY_RSCRIPT_INTEROPERABILITY = True
rsubprocess_pipe_path = libexec.filepath("rsubprocess.R")

class RService(object):

    @staticmethod
    def call(r_commands,
            cwd=None,
            env=None,
            rscript_path=RSCRIPT_EXECUTABLE,
            ):
        """
        Executes a sequence of commans in R and returns the results.

        Parameters
        ----------
        r_commands : iterable of strings
            A list or some other iterable of strings of R commands.
        cwd : string
            Set the working directory of the PAUP* process to this directory.
        env : dictionary
            Environmental variables to set for the PAUP* process.
        rscript_path : string
            Path to the Rscript executable.

        Returns
        -------
        returncode : exit value of the R process
        stdout : string
            Contents of the R process standard output.
        stderr : string
            Contents of the R process standard error.
        """
        if not isinstance(r_commands, str):
            r_commands = "\n".join(r_commands)
        r_commands += "\n"
        p = subprocess.Popen(
                [RSCRIPT_EXECUTABLE, rsubprocess_pipe_path],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=cwd,
                env=env,
                )
        stdout, stderr = processio.communicate(p, r_commands)
        return p.returncode, stdout, stderr

def call(*args, **kwargs):
    return RService.call(*args, **kwargs)

