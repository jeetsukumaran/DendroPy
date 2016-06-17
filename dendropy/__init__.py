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

import sys

###############################################################################
## Populate the 'dendropy' namespace

from dendropy.dataio.nexusprocessing import get_rooting_argument
from dendropy.datamodel.basemodel import Annotation
from dendropy.datamodel.basemodel import AnnotationSet
from dendropy.datamodel.taxonmodel import Taxon
from dendropy.datamodel.taxonmodel import TaxonNamespace
from dendropy.datamodel.taxonmodel import TaxonNamespacePartition
from dendropy.datamodel.taxonmodel import TaxonNamespaceMapping
from dendropy.datamodel.taxonmodel import TaxonSet # Legacy
from dendropy.datamodel.treemodel import Bipartition
from dendropy.datamodel.treemodel import Edge
from dendropy.datamodel.treemodel import Node
from dendropy.datamodel.treemodel import Tree
from dendropy.datamodel.treecollectionmodel import TreeList
from dendropy.datamodel.treecollectionmodel import SplitDistribution
from dendropy.datamodel.treecollectionmodel import TreeArray
from dendropy.datamodel.charstatemodel import StateAlphabet
from dendropy.datamodel.charstatemodel import DNA_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import RNA_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import NUCLEOTIDE_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import PROTEIN_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import BINARY_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import RESTRICTION_SITES_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import INFINITE_SITES_STATE_ALPHABET
from dendropy.datamodel.charstatemodel import new_standard_state_alphabet
from dendropy.datamodel.charmatrixmodel import CharacterDataSequence
from dendropy.datamodel.charmatrixmodel import CharacterMatrix
from dendropy.datamodel.charmatrixmodel import DnaCharacterDataSequence
from dendropy.datamodel.charmatrixmodel import DnaCharacterMatrix
from dendropy.datamodel.charmatrixmodel import NucleotideCharacterDataSequence
from dendropy.datamodel.charmatrixmodel import NucleotideCharacterMatrix
from dendropy.datamodel.charmatrixmodel import RnaCharacterDataSequence
from dendropy.datamodel.charmatrixmodel import RnaCharacterMatrix
from dendropy.datamodel.charmatrixmodel import ProteinCharacterDataSequence
from dendropy.datamodel.charmatrixmodel import ProteinCharacterMatrix
from dendropy.datamodel.charmatrixmodel import RestrictionSitesCharacterDataSequence
from dendropy.datamodel.charmatrixmodel import RestrictionSitesCharacterMatrix
from dendropy.datamodel.charmatrixmodel import InfiniteSitesCharacterDataSequence
from dendropy.datamodel.charmatrixmodel import InfiniteSitesCharacterMatrix
from dendropy.datamodel.charmatrixmodel import StandardCharacterDataSequence
from dendropy.datamodel.charmatrixmodel import StandardCharacterMatrix
from dendropy.datamodel.charmatrixmodel import ContinuousCharacterDataSequence
from dendropy.datamodel.charmatrixmodel import ContinuousCharacterMatrix
from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix
from dendropy.datamodel.datasetmodel import DataSet
from dendropy.utility.error import ImmutableTaxonNamespaceError
from dendropy.utility.error import DataParseError
from dendropy.utility.error import UnsupportedSchemaError
from dendropy.utility.error import UnspecifiedSchemaError
from dendropy.utility.error import UnspecifiedSourceError
from dendropy.utility.error import TooManyArgumentsError
from dendropy.utility.error import InvalidArgumentValueError
from dendropy.utility.error import MultipleInitializationSourceError
from dendropy.utility.error import TaxonNamespaceIdentityError
from dendropy.utility.error import TaxonNamespaceReconstructionError
from dendropy.utility.error import UltrametricityError
from dendropy.utility.error import TreeSimTotalExtinctionException
from dendropy.utility.error import SeedNodeDeletionException
from dendropy.utility import deprecate


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
                minor=2,
                micro=0,
                releaselevel=""
                )
__project__ = "DendroPy"
__version__ = ".".join(str(s) for s in version_info[:4] if s != "")
__author__ = "Jeet Sukumaran and Mark T. Holder"
__copyright__ = "Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder."
__citation__ = "Sukumaran, J and MT Holder. 2010. DendroPy: a Python library for phylogenetic computing. Bioinformatics 26: 1569-1571."
PACKAGE_VERSION = __version__ # for backwards compatibility (with sate)

def _get_revision_object():
    from dendropy.utility import vcsinfo
    __revision__ = vcsinfo.Revision(repo_path=homedir())
    return __revision__

def revision_description():
    __revision__ = _get_revision_object()
    if __revision__.is_available:
        revision_text = " ({})".format(__revision__)
    else:
        revision_text = ""
    return revision_text

def name():
    return "{} {}{}".format(__project__, __version__, revision_description())

def homedir():
    import os
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
    return __homedir__

def description(dest=None):
    import sys
    import site
    if dest is None:
        dest = sys.stdout
    fields = collections.OrderedDict()
    fields["DendroPy version"] = name()
    fields["DendroPy location"] = homedir()
    fields["Python version"] = sys.version.replace("\n", "")
    fields["Python executable"] = sys.executable
    try:
        fields["Python site packages"] = site.getsitepackages()
    except:
        pass
    max_fieldname_len = max(len(fieldname) for fieldname in fields)
    for fieldname, fieldvalue in fields.items():
        dest.write("{fieldname:{fieldnamewidth}}: {fieldvalue}\n".format(
            fieldname=fieldname,
            fieldnamewidth=max_fieldname_len + 2,
            fieldvalue=fieldvalue))

def description_text():
    from dendropy.utility.textprocessing import StringIO
    s = StringIO()
    description(s)
    return s.getvalue()

def citation_info(include_preamble=True, width=76):
    import textwrap
    citation_lines = []
    if include_preamble:
        citation_preamble =(
                            "If any stage of your work or analyses relies"
                            " on code or programs from this library, either"
                            " directly or indirectly (e.g., through usage of"
                            " your own or third-party programs, pipelines, or"
                            " toolkits which use, rely on, incorporate, or are"
                            " otherwise primarily derivative of code/programs"
                            " in this library), please cite:"
                        )
        citation_lines.extend(textwrap.wrap(citation_preamble, width=width))
        citation_lines.append("")
    citation = textwrap.wrap(
            __citation__,
            width=width,
            initial_indent="  ",
            subsequent_indent="    ",
            )
    citation_lines.extend(citation)
    return citation_lines

def tree_source_iter(*args, **kwargs):
    s = "No longer supported in DendroPy 4: Instead of 'tree_source_iter()', use 'Tree.yield_from_files()' instead"
    raise NotImplementedError(s)

def multi_tree_source_iter(*args, **kwargs):
    s = "No longer supported in DendroPy 4: Instead of 'multi_tree_source_iter()', use 'Tree.yield_from_files()' instead"
    raise NotImplementedError(s)

if __name__ == "__main__":
    description(sys.stdout)



