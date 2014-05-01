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

import sys
import os

###############################################################################
## Populate the 'dendropy' namespace

from dendropy.datamodel.taxon import Taxon
from dendropy.datamodel.taxon import TaxonNamespace
from dendropy.datamodel.tree import Edge
from dendropy.datamodel.tree import Node
from dendropy.datamodel.tree import Tree
from dendropy.datamodel.tree import TreeList
from dendropy.datamodel.char import DnaCharacterMatrix
from dendropy.datamodel.char import RnaCharacterMatrix
from dendropy.datamodel.char import ProteinCharacterMatrix
from dendropy.datamodel.char import StandardCharacterMatrix
from dendropy.datamodel.char import RestrictionSitesCharacterMatrix
from dendropy.datamodel.char import InfiniteSitesCharacterMatrix

###############################################################################
## PACKAGE METADATA

__project__ = "DendroPy"
__version__ = "4.0.0"
__author__ = "Jeet Sukumaran and Mark T. Holder"
__copyright__ = "Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder."
__citation__ = "Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library for phylogenetic computing. Bioinformatics 26: 1569-1571."
PACKAGE_VERSION = __version__ # for backwards compatibility (with sate)

def revision():
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
    sys.stdout.write("{}\n".format(description()))


