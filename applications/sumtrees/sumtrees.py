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
Summarizations collections of trees, e.g., MCMC samples from a posterior
distribution, non-parametric bootstrap replicates, mapping posterior
probability, support, or frequency that splits/clades are found in the source
set of trees onto a target tree.
"""

import os
import sys
import re
import argparse

if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
try:
    # Python 3
    import queue
except ImportError:
    # Python 2.7
    import Queue as queue
import multiprocessing

from dendropy.utility.cli import confirm_overwrite, show_splash

_program_name = "SumTrees"
_program_subtitle = "Phylogenetic Tree Summarization"
_program_date = "Jan 31 2015"
_program_version = "Version 4.0.0 (%s)" % _program_date
_program_author = "Jeet Sukumaran and Mark T. Holder"
_program_contact = "jeetsukumaran@gmail.com"
_program_copyright = """\
Copyright (C) 2008-2014 Jeet Sukumaran and Mark T. Holder.
License GPLv3+: GNU GPL version 3 or later.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law."""

_program_citation = """\
Sukumaran, J and MT Holder. {prog_name}: {prog_subtitle}. {prog_version}. Available at https://github.com/jeetsukumaran/DendroPy.
""".format(prog_name=_program_name, prog_subtitle=_program_subtitle, prog_version=_program_version)

def main():
    parser = argparse.ArgumentParser(description=__doc__)

    run_optgroup = parser.add_argument_group("Program Run Options")
    run_optgroup.add_argument("-m", "--multiprocessing",
            dest="multiprocess",
            metavar="NUM-PROCESSES",
            default=None,
            help="Run in parallel mode with up to a maximum of NUM-PROCESSES processes " \
                 "(specify '*' to run in as many processes as there are cores on the "\
                 "local machine).")
    run_optgroup.add_argument("-g", "--log-frequency",
            type=int,
            metavar="LOG-FREQUENCY",
            default=500,
            help="Tree processing progress logging frequency (default: %(default)s; set to 0 to suppress).")
    run_optgroup.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Suppress ALL logging, progress and feedback messages.")
    run_optgroup.add_argument("--ignore-missing-support",
            action="store_true",
            default=False,
            help="Ignore missing support tree files (at least one must exist!).")

    args = parser.parse_args()

    if not args.quiet:
        show_splash(prog_name=_program_name,
                prog_subtitle=_program_subtitle,
                prog_version=_program_version,
                prog_author=_program_author,
                prog_copyright=_program_copyright,
                dest=sys.stderr,
                include_citation=True,
                include_copyright=False,
                additional_citations=[_program_citation],
                )

if __name__ == '__main__':
    main()
