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
Implementation of NEWICK-schema tree iterator.
"""

from dendropy.dataio import ioservice
from dendropy.dataio import newickreader
from dendropy.dataio import nexusprocessing

class NewickTreeDataYielder(ioservice.TreeDataYielder):

    def __init__(self,
            files=None,
            taxon_namespace=None,
            tree_type=None,
            **kwargs):
        r"""

        Parameters
        ----------
        files : iterable of sources
            Iterable of sources, which can either be strings specifying file
            paths or file-like objects open for reading. If a source element is
            a stringm then it is assumed to be a path to a file. Otherwise, the
            source is assumed to be a file-like object.
        taxon_namespace : |TaxonNamespace| instance
            The operational taxonomic unit concept namespace to use to manage
            taxon definitions.
        \*\*kwargs : keyword arguments
            These will be passed directly to the base `newickreader.NexusReader`
            class. See `newickreader.NexusReader` for details.
        """
        ioservice.TreeDataYielder.__init__(self,
                files=files,
                taxon_namespace=taxon_namespace,
                tree_type=tree_type)
        self.newick_reader = newickreader.NewickReader(**kwargs)

    ###########################################################################
    ## Implementation of DataYielder interface

    def _yield_items_from_stream(self, stream):
        nexus_tokenizer = nexusprocessing.NexusTokenizer(stream,
                preserve_unquoted_underscores=self.newick_reader.preserve_unquoted_underscores)
        taxon_symbol_mapper = nexusprocessing.NexusTaxonSymbolMapper(
                taxon_namespace=self.attached_taxon_namespace,
                enable_lookup_by_taxon_number=False,
                case_sensitive=self.newick_reader.case_sensitive_taxon_labels)
        while True:
            tree = self.newick_reader._parse_tree_statement(
                    nexus_tokenizer=nexus_tokenizer,
                    tree_factory=self.tree_factory,
                    taxon_symbol_map_fn=taxon_symbol_mapper.require_taxon_for_symbol)
            if tree is None:
                break
            yield tree
