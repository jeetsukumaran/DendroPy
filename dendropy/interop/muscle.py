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
Wrapper around calls to MUSCLE
"""

import dendropy
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


