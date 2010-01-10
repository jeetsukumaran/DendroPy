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
Wrappers for interacting with the Biopython library.
"""

import tempfile
import re
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
        Returns `o` as an biopython object.
        """
        if isinstance(o, dendropy.CharacterMatrix):
            if len(o.state_alphabets) > 1:
                raise ValueError('Multi-alphabet matrices not supported in Biopython')
            sa = o.state_alphabets[0]
            if isinstance(sa, dendropy.DnaCharacterMatrix):
                bpa = generic_dna
            elif isinstance(sa, dendropy.RnaCharacterMatrix):
                bpa = generic_rna
            elif isinstance(sa, dendropy.ProteinCharacterMatrix):
                bpa = generic_protein
            else:
                raise ValueError("Character data type not supported in Biopython: '%s'" % type(sa))
        else:
            raise ValueError("Invalid object type for conversion to Biopython: '%s'" % type(o))
