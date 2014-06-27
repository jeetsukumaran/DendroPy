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


import collections
from dendropy.dataio import newickreader
from dendropy.dataio import newickwriter
from dendropy.dataio import fastareader
from dendropy.dataio import nexusreader
from dendropy.dataio import nexuswriter
from dendropy.utility import container

_IOServices = collections.namedtuple(
        "_IOServices",
        ["reader", "writer", "tree_iterator"]
        )

_IO_SERVICE_REGISTRY = container.CaseInsensitiveDict()
_IO_SERVICE_REGISTRY["newick"] = _IOServices(newickreader.NewickReader, newickwriter.NewickWriter, None)
_IO_SERVICE_REGISTRY["fasta"] = _IOServices(fastareader.FastaReader, None, None)
_IO_SERVICE_REGISTRY["nexus"] = _IOServices(nexusreader.NexusReader, nexuswriter.NexusWriter, None)
_IO_SERVICE_REGISTRY["dnafasta"] = _IOServices(fastareader.DnaFastaReader, None, None)
_IO_SERVICE_REGISTRY["rnafasta"] = _IOServices(fastareader.RnaFastaReader, None, None)
_IO_SERVICE_REGISTRY["proteinfasta"] = _IOServices(fastareader.ProteinFastaReader, None, None)
# ioclient.register("nexus", nexusreader_py.NexusReader, nexuswriter.NexusWriter, nexustreeiter.tree_source_iter)
# ioclient.register("newick", newick.NewickReader, newick.NewickWriter, newick.tree_source_iter)
# ioclient.register("nexus/newick", None, None, nexustreeiter.generalized_tree_source_iter)
# ioclient.register("fasta", fasta.FastaReader, fasta.FastaWriter, None)
# ioclient.register("dnafasta", fasta.DNAFastaReader, fasta.FastaWriter, None)
# ioclient.register("rnafasta", fasta.RNAFastaReader, fasta.FastaWriter, None)
# ioclient.register("proteinfasta", fasta.ProteinFastaReader, fasta.FastaWriter, None)
# ioclient.register("phylip", phylip.PhylipReader, phylip.PhylipWriter, None)
# ioclient.register("nexml", nexml.NexmlReader, nexml.NexmlWriter, None)
# ioclient.register("beast-summary-tree", beast.BeastSummaryTreeReader, None, beast.summary_tree_source_iter)

def get_reader(schema, **kwargs):
    try:
        reader_type =_IO_SERVICE_REGISTRY[schema].reader
        if reader_type is None:
            raise KeyError
        reader = reader_type(**kwargs)
        return reader
    except KeyError:
        raise NotImplementedError("'{}' is not a supported data reading schema".format(schema))

def get_writer(
        schema,
        **kwargs):
    try:
        writer_type =_IO_SERVICE_REGISTRY[schema].writer
        if writer_type is None:
            raise KeyError
        writer = writer_type(**kwargs)
        return writer
    except KeyError:
        raise NotImplementedError("'{}' is not a supported data writing schema".format(schema))

