#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2026 Jeet Sukumaran, Mark T. Holder, and Matthew Andres Moreno.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Moreno, M. A., Holder, M. T., & Sukumaran, J. (2024). DendroPy 5: a
##     mature Python library for phylogenetic computing. Journal of Open
##     Source Software, 9(101), 6943, https://doi.org/10.21105/joss.06943
##
##############################################################################

"""
Wrapper around calls to MUSCLE
"""

import subprocess
from dendropy.utility import processio

def muscle_align(char_matrix, muscle_args=None, muscle_path='muscle'):
    cmd = [muscle_path]
    if muscle_args:
        cmd = cmd + muscle_args
    p = subprocess.Popen(cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    stdout, stderr = processio.communicate(p, char_matrix.as_string("fasta"))
    if p.returncode:
        raise Exception(stderr)
    d = char_matrix.__class__.get_from_string(stdout,
            "fasta",
            taxon_namespace=char_matrix.taxon_namespace)
    return d


