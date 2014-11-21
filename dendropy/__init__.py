#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Imports into the `dendropy` namespace all fundamental
classes and methods for instantiating objects in the
`dendropy.dataobject` subpackage to for usage by client code.
"""

###############################################################################
## Populate the 'dendropy' namespace

from dendropy.datamodel.taxonmodel import Taxon
from dendropy.datamodel.taxonmodel import TaxonNamespace
from dendropy.datamodel.taxonmodel import TaxonNamespacePartition
from dendropy.datamodel.taxonmodel import TaxonNamespaceMapping
from dendropy.datamodel.taxonmodel import TaxonSet # Legacy
from dendropy.datamodel.treemodel import Bipartition
from dendropy.datamodel.treemodel import Edge
from dendropy.datamodel.treemodel import Node
from dendropy.datamodel.treemodel import Tree
from dendropy.datamodel.treemodel import TreeList
from dendropy.datamodel.treemodel import SplitDistribution
from dendropy.datamodel.treemodel import TreeArray
from dendropy.datamodel.charstatemodel import StateAlphabet
from dendropy.datamodel.charstatemodel import DNA_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import RNA_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import NUCLEOTIDE_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import PROTEIN_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import BINARY_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import RESTRICTION_SITES_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import INFINITE_SITES_STATE_ALPHABET
from dendropy.datamodel.charmatrixmodel import CharacterMatrix
from dendropy.datamodel.charmatrixmodel import CharacterDataSequence
from dendropy.datamodel.charmatrixmodel import DnaCharacterMatrix
from dendropy.datamodel.charmatrixmodel import RnaCharacterMatrix
from dendropy.datamodel.charmatrixmodel import NucleotideCharacterMatrix
from dendropy.datamodel.charmatrixmodel import ProteinCharacterMatrix
from dendropy.datamodel.charmatrixmodel import RestrictionSitesCharacterMatrix
from dendropy.datamodel.charmatrixmodel import InfiniteSitesCharacterMatrix
from dendropy.datamodel.charmatrixmodel import StandardCharacterMatrix
from dendropy.datamodel.charmatrixmodel import ContinuousCharacterMatrix
from dendropy.datamodel.datasetmodel import DataSet

###############################################################################
## Legacy Support

from dendropy.legacy import coalescent
from dendropy.legacy import continuous
from dendropy.legacy import treecalc
from dendropy.legacy import popgensim
from dendropy.legacy import popgenstat
from dendropy.legacy import reconcile
from dendropy.legacy import seqmodel
from dendropy.legacy import seqsim
from dendropy.legacy import treecalc
from dendropy.legacy import treemanip
from dendropy.legacy import treesim
from dendropy.legacy import treesplit
from dendropy.legacy import treesum

###############################################################################
## PACKAGE METADATA
import collections
version_info = collections.namedtuple("dendropy_version_info",
        ["major", "minor", "micro", "releaselevel"])(
                major=4,
                minor=0,
                micro=0,
                releaselevel="development"
                )
__project__ = "DendroPy"
__version__ = ".".join(str(s) for s in version_info[:3])
__author__ = "Jeet Sukumaran and Mark T. Holder"
__copyright__ = "Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder."
__citation__ = "Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library for phylogenetic computing. Bioinformatics 26: 1569-1571."
PACKAGE_VERSION = __version__ # for backwards compatibility (with sate)

def revision():
    import os
    from dendropy.utility import vcsinfo
    try:
        try:
            __homedir__ = __path__[0]
        except AttributeError:
            __homedir__ = os.path.dirname(os.path.abspath(__file__))
        except IndexError:
            __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except OSError:
        __homedir__ = None
    except:
        __homedir__ = None
    __revision__ = vcsinfo.Revision(repo_path=__homedir__)
    return __revision__

def description():
    __revision__ = revision()
    if __revision__.is_available:
        revision_text = " ({})".format(__revision__)
    else:
        revision_text = ""
    return "{} {}{}".format(__project__, __version__, revision_text)

if __name__ == "__main__":
    import sys
    sys.stdout.write("{}\n".format(description()))


