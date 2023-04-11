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

import os
import sys
import subprocess
from dendropy.utility import error
from dendropy.utility import metavar
from dendropy.utility import libexec
from dendropy.utility import processio
from dendropy.utility import textprocessing

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
        r"""
        Executes a sequence of commands in R and returns the results. All the
        noise is sunk into the stderr return variable, and just the output
        comes out cleanly in the stdout return variable.

        Parameters
        ----------
        r_commands : iterable of strings
            A list or some other iterable of strings of R commands.
        ignore_error_returncode : bool
            If |True|, then a non-0 return code from the R process will not
            result in an exception being raised.
        cwd : string
            Set the working directory of the R process to this directory.
        env : dictionary
            Environmental variables to set for the R process.
        rscript_path : string
            Path to the Rscript executable.

        Returns
        -------
        returncode : exit value of the R process
        stdout : string
            Contents of the R process standard output.
        stderr : string
            Contents of the R process standard error.

        Examples
        --------

        Build up a script (``s``) to calculate a range of values, print them
        to the standard output, and then post-process this to extract the
        values::

            import itertools
            from dendropy.interop import rstats

            bb = [0.01, 0.05, 0.10, 0.50, 1.0]
            cc = [0.01, 0.05, 0.10, 0.50, 1.0]
            ee = [0.0, 0.1, 0.2]

            # store commands of script as a list
            # to be passed to the ``call()``
            s = []

            # set options, load required libraries, etc.
            s.append("options(digits=22)")
            s.append("library(PBD)")

            # build up list of commands in script
            params = []
            for b, c, e in itertools.product(bb, cc, ee):
                s.append("print(pbd_durspec_mean(pars=c({},{},{})))".format(b, c, e))

            # execute script
            returncode, stdout, stderr  = rstats.call(s)

            # peek at the results
            print(stdout)

            # [1] 69.31472
            # [1] 9.853723
            # [1] 4.981369
            # [1] 0.9950331
            # ...

            # post-process the stdout to extract values
            results = [float(x.split(" ")[1]) for x in stdout.split("\n") if x]

        Notes
        -----

        Note that newlines ('\n') and other special characters will be
        converted before being passed to the R interpreter, so need to
        be escaped or entered as raw string expressions.

        That is, instead of, e.g.::

            returncode, stdout, stderr = RService.call([
                "cat('hello, world\n')",
            ])

        use this::

            returncode, stdout, stderr = RService.call([
                "cat('hello, world\n')",
            ])

        or::

            returncode, stdout, stderr = RService.call([
                r"cat('hello, world\n')",
            ])

        """
        if not textprocessing.is_str_type(r_commands):
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

