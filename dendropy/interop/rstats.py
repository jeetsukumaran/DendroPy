#! /usr/bin/env python

import os
import sys
import subprocess
from dendropy.utility import error
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
            ignore_error_returncode=False,
            cwd=None,
            env=None,
            rscript_path=RSCRIPT_EXECUTABLE,
            ):
        """
        Executes a sequence of commans in R and returns the results.

        Note that newlines ('\n') and other special characters will be
        converted before being passed to the R interpreter, so need to
        be escaped or entered as raw string expressions.

        That is, instead of, e.g.:

            returncode, stdout, stderr = RService.call([
                "cat('hello, world\n')",
            ])

        use this:

            returncode, stdout, stderr = RService.call([
                "cat('hello, world\\n')",
            ])

        or:

            returncode, stdout, stderr = RService.call([
                r"cat('hello, world\n')",
            ])

        Parameters
        ----------
        r_commands : iterable of strings
            A list or some other iterable of strings of R commands.
        ignore_error_returncode : bool
            If `True`, then a non-0 return code from the PAUP process will not
            result in an exception being raised.
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
        invocation_command = [RSCRIPT_EXECUTABLE, rsubprocess_pipe_path]
        p = subprocess.Popen(
                invocation_command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=cwd,
                env=env,
                )
        stdout, stderr = processio.communicate(p, r_commands)
        if (p.returncode != 0 and not ignore_error_returncode):
            raise error.ExternalServiceError(
                    service_name="Rscript",
                    invocation_command=invocation_command,
                    service_input=r_commands,
                    returncode = p.returncode,
                    stdout=stdout,
                    stderr=stderr)
        return p.returncode, stdout, stderr

def call(*args, **kwargs):
    return RService.call(*args, **kwargs)

