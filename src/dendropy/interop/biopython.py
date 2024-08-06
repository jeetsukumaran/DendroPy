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
Wrappers for interacting with the Biopython library.
"""

from dendropy.utility import messaging
_LOG = messaging.get_logger(__name__)
import dendropy

DENDROPY_BIOPYTHON_INTEROPERABILITY = False
try:
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna, generic_rna, generic_nucleotide, generic_protein
    DENDROPY_BIOPYTHON_INTEROPERABILITY = True
except ImportError:
    _LOG.warn("BioPython could not be imported: BioPython interoperability not available")
else:

    def as_biopython_object(o):
        """
        Returns ``o`` as an biopython object.
        """
        if isinstance(o, dendropy.CharacterMatrix):
            if isinstance(o, dendropy.DnaCharacterMatrix):
                bpa = generic_dna
            elif isinstance(o, dendropy.RnaCharacterMatrix):
                bpa = generic_rna
            elif isinstance(o, dendropy.ProteinCharacterMatrix):
                bpa = generic_protein
            else:
                raise ValueError("Character data type not supported in Biopython: '%s'" % type(o))
        else:
            raise ValueError("Invalid object type for conversion to Biopython: '%s'" % type(o))

