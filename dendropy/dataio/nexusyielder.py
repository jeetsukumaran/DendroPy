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
Implementation of NEXUS-schema tree iterator.
"""

import sys
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
from dendropy.utility import textprocessing
from dendropy.dataio import ioservice
from dendropy.dataio import nexusreader
from dendropy.dataio import nexusprocessing

class NexusTreeDataYielder(
        ioservice.TreeDataYielder,
        nexusreader.NexusReader):

    def __init__(self,
            files=None,
            taxon_namespace=None,
            tree_type=None,
            **kwargs):
        """

        Parameters
        ----------
        files : iterable of sources
            Iterable of sources, which can either be strings specifying file
            paths or file-like objects open for reading. If a source element is
            a string (``isinstance(i,str) == True``), then it is assumed to be
            a path to a file. Otherwise, the source is assumed to be a file-like
            object.
        taxon_namespace : |TaxonNamespace| instance
            The operational taxonomic unit concept namespace to use to manage
            taxon definitions.
        \*\*kwargs : keyword arguments
            These will be passed directly to the base `nexusreader.NexusReader`
            class. See `nexusreader.NexusReader` for details.
        """
        ioservice.TreeDataYielder.__init__(self,
                files=files,
                taxon_namespace=taxon_namespace,
                tree_type=tree_type)
        self.assume_newick_if_not_nexus = kwargs.pop("assume_newick_if_not_nexus", False)
        kwargs["attached_taxon_namespace"] = self.attached_taxon_namespace
        nexusreader.NexusReader.__init__(self, **kwargs)
        self.exclude_chars = True
        self.exclude_trees = False

    ###########################################################################
    ## Implementation of DataYielder interface

    def _yield_items_from_stream(self, stream):
        if self._nexus_tokenizer is None:
            self.create_tokenizer(stream,
                preserve_unquoted_underscores=self.preserve_underscores)
        else:
            self._nexus_tokenizer.set_stream(stream)
        token = self._nexus_tokenizer.next_token()
        if token.upper() != "#NEXUS":
            if self.assume_newick_if_not_nexus:
                taxon_symbol_mapper = self._get_taxon_symbol_mapper(
                        taxon_namespace=self.attached_taxon_namespace,
                        enable_lookup_by_taxon_number=False,
                        )
                while True:
                    tree = self._build_tree_from_newick_tree_string(
                            tree_factory=self.tree_factory,
                            taxon_symbol_mapper=taxon_symbol_mapper)
                    if tree is None:
                        break
                    yield tree
            else:
                raise self._nexus_error("Expecting '#NEXUS', but found '{}'".format(token),
                        nexusreader.NexusReader.NotNexusFileError)
        while not self._nexus_tokenizer.is_eof():
            token = self._nexus_tokenizer.next_token_ucase()
            while token != None and token != 'BEGIN' and not self._nexus_tokenizer.is_eof():
                token = self._nexus_tokenizer.next_token_ucase()
            self._nexus_tokenizer.process_and_clear_comments_for_item(
                    self._global_annotations_target,
                    self.extract_comment_metadata)
            token = self._nexus_tokenizer.next_token_ucase()
            if token == 'TAXA':
                self._parse_taxa_block()
            elif token == 'TREES':
                for tree in self._yield_from_trees_block():
                    yield tree
            elif token == 'BEGIN':
                raise self._nexus_error("'BEGIN' found without completion of previous block",
                        nexusreader.NexusReader.IncompleteBlockError)
            else:
                # unknown block
                token = self._consume_to_end_of_block(token)

    ###########################################################################
    ## Supporting Functions

    def _yield_from_trees_block(self):
        """
        Expectations:
            - current token: "TREES" [part of "BEGIN TREES"]
        """
        token = self._nexus_tokenizer.cast_current_token_to_ucase()
        if token != "TREES":
            raise self._nexus_error("Expecting 'TREES' token, but instead found '{}'".format(token))
        if self.exclude_trees:
            self._consume_to_end_of_block(self._nexus_tokenizer.current_token)
            return
        self._nexus_tokenizer.skip_to_semicolon() # move past "BEGIN TREES" command
        link_title = None
        taxon_namespace = None
        taxon_symbol_mapper = None
        trees_block = None
        block_title = None
        while ((not self._nexus_tokenizer.is_eof())
                and token is not None
                and token != 'END'
                and token != 'ENDBLOCK'):
            token = self._nexus_tokenizer.next_token_ucase()
            if token == 'LINK':
                link_title = self._parse_link_statement().get("taxa")
            elif token == 'TITLE':
                block_title = self._parse_title_statement()
                token = "" # clear; repopulate at start of loop
            elif token == 'TRANSLATE':
                if taxon_namespace is None:
                    taxon_namespace = self._get_taxon_namespace(link_title)
                taxon_symbol_mapper = self._parse_translate_statement(taxon_namespace)
                token = "" # clear; repopulate at start of loop
            elif token == 'TREE':
                if taxon_namespace is None:
                    taxon_namespace = self._get_taxon_namespace(link_title)
                if taxon_symbol_mapper is None:
                    taxon_symbol_mapper = self._get_taxon_symbol_mapper(taxon_namespace=taxon_namespace)
                pre_tree_comments = self._nexus_tokenizer.pull_captured_comments()
                tree_factory = self.tree_factory
                while True:
                    ## After the following, the current token
                    ## will be the token immediately following
                    ## the terminating semi-colon of a tree
                    ## statement. Typically, this will be
                    ## 'TREE' if there is another tree, or
                    ## 'END'/'ENDBLOCK'.
                    tree = self._parse_tree_statement(
                            tree_factory=tree_factory,
                            taxon_symbol_mapper=taxon_symbol_mapper)
                    yield tree
                    if self._nexus_tokenizer.is_eof() or not self._nexus_tokenizer.current_token:
                        break
                    if self._nexus_tokenizer.cast_current_token_to_ucase() != "TREE":
                        token = self._nexus_tokenizer.current_token
                        break
            elif token == 'BEGIN':
                raise self._nexus_error("'BEGIN' found without completion of previous block",
                        nexusreader.NexusReader.IncompleteBlockError)
        self._nexus_tokenizer.skip_to_semicolon() # move past END command
        raise StopIteration

class NexusNewickTreeDataYielder(NexusTreeDataYielder):

    def __init__(self,
            files=None,
            taxon_namespace=None,
            tree_type=None,
            **kwargs):
        kwargs["assume_newick_if_not_nexus"] = kwargs.get("assume_newick_if_not_nexus", True)
        NexusTreeDataYielder.__init__(self,
                files=files,
                taxon_namespace=taxon_namespace,
                tree_type=tree_type,
                **kwargs)
