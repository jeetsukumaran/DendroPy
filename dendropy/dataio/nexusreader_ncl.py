#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
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
Facultative use of NCL for NEXUS parsing.
"""

import os
from dendropy.utility.messaging import get_logger
_LOG = get_logger("dataio.ncl")

DENDROPY_NCL_AVAILABILITY = False
try:
    from nexusclasslib import nclwrapper
    DENDROPY_NCL_AVAILABILITY = True
except ImportError:
    DENDROPY_NCL_AVAILABILITY = False
else:

    import os
    from threading import Thread, Event
    from dendropy import dataobject
    from dendropy.dataio import nexusreader_py
    from dendropy.dataio import nexustokenizer
    from dendropy.utility import iosys

    if "DENDROPY_ENABLE_NCL_WARNINGS" in os.environ:
        DENDROPY_NCL_WARNING_LEVEL = nclwrapper.NxsReader.SKIPPING_CONTENT_WARNING
    else:
        DENDROPY_NCL_WARNING_LEVEL = nclwrapper.NxsReader.SUPPRESS_WARNINGS_LEVEL

    def _ncl_datatype_enum_to_dendropy(d):
        e = nclwrapper.NxsCharactersBlock
        if d == e.dna:
            return dataobject.DnaCharacterMatrix
        if d == e.nucleotide:
            return dataobject.NucleotideCharacterMatrix
        if d == e.rna:
            return dataobject.RnaCharacterMatrixk
        if d == e.protein:
            return dataobject.ProteinCharacterMatrix
        if (d == e.continuous):
            return dataobject.ContinuousCharacterMatrix
        if d == e.standard:
            return dataobject.StandardCharacterMatrix
        if (d == e.mixed) or (d == e.codon):
            s = d == e.continuous and "continuous" or (d == e.mixed and "mixed" or "codon")
            raise NotImplementedError("%s datatype not supported" % s)
        sys.exit(0)

    class NCLTreeStream(nclwrapper.NxsTreeStream):
        """Simple thread-safe class that waits for `need_tree_event', and signals the
        presence of a new tree by `ready_event`"""
        def __init__(self, need_tree_event, ready_event, die_event):
            self.need_tree_event = need_tree_event
            self.ready_event = ready_event
            self.tree_tokens = None
            self.ncl_taxa_block =  None
            self.exception = None
            self.die_event = die_event
            nclwrapper.NxsTreeStream.__init__(self)

        def handleTree(self, ftd, tb):
            t = ftd.GetTreeTokens()
            rooted_flag = ftd.IsRooted()
            if self.die_event.isSet():
                raise RuntimeError("exception in calling thread")
            self.need_tree_event.wait()
            self.need_tree_event.clear()
            if self.die_event.isSet():
                raise RuntimeError("exception in calling thread")
            try:
                self.ncl_taxa_block =  tb.GetTaxaBlockPtr()
                tb_iid = self.ncl_taxa_block.GetInstanceIdentifierString()
                self.tree_tokens =  t
                self.rooted_flag = rooted_flag
            except Exception, v:
                _LOG.debug("NCLTreeStream Exception: %s" % str(v))
                self.exception = v
                self.ncl_taxa_block = None
                self.tree_tokens = None
                self.rooted_flag = None
                self.ready_event.set()
                raise v

            self.ready_event.set()
            return False

    class NCLTreeStreamThread(Thread):
        """Subclass of thread ,that uses a NCLTreeStream to get trees
        from NCL one at a time"""
        def __init__(self, file_path, need_tree_event, ready_event, die_event, format="NEXUS", **kwargs):
            self.nts = NCLTreeStream(need_tree_event, ready_event, die_event)
            self.file_path = file_path
            self.format = format
            self.exception = None
            self.done = False
            self.reader = nclwrapper.MultiFormatReader()
            self.reader.SetWarningOutputLevel(DENDROPY_NCL_WARNING_LEVEL)
            self.reader.cullIdenticalTaxaBlocks(True)
            self.die_event = die_event
            Thread.__init__(self,
                            group=None,
                            target=None,
                            name=None,
                            args=tuple(),
                            kwargs=dict(**kwargs))
        def run(self):
            self.done = False
            self.exception = None
            try:
                if not self.die_event.isSet():
                    self.nts.ReadFilepath(self.file_path, self.format, self.reader)
            except Exception, v:
                if self.nts.exception:
                    self.exception = self.nts.exception
                else:
                    self.exception = v
                _LOG.debug("NCLTreeStreamThread Exception: %s" % str(self.exception))
            else:
                self.nts.need_tree_event.wait()
                self.nts.need_tree_event.clear()
            self.done = True
            self.nts.tree_tokens = None
            self.nts.taxa_block = None
            self.rooted_flag = None
            self.nts.ready_event.set()

    class ListOfTokenIterator(object):
        def __init__(self, tokens):
            self.tokens_iter = iter(tokens)
            self.eof = False
            self.queued = None
            self.tree_rooted = None
            self.comments = None

        def clear_comments(self):
            self.tree_rooted = None
            self.comments = None

        def tree_rooted_comment(self):
            "This is a hack, and only works if you have just one tree in the token stream"
            return self.tree_rooted

        def __iter__(self):
            return self
        def read_next_token(self, ignore_punctuation=None):
            if not self.eof:
                try:
                    return self.tokens_iter.next()
                except StopIteration:
                    self.eof = True
            return ""

        def read_next_token_ucase(self):
            t = self.read_next_token()
            if t:
                return t.upper()
        def syntax_exception(self, msg):
            return SyntaxException(message=msg)

    class NCLBasedReader(iosys.DataReader):
        "Encapsulates loading and parsing of a NEXUS format file."

        def __init__(self, schema="NEXUS", **kwargs):
            iosys.DataReader.__init__(self, **kwargs)
            self.encode_splits = kwargs.get("encode_splits", False)
            self.rooting_interpreter = kwargs.get("rooting_interpreter", nexustokenizer.RootingInterpreter(**kwargs))
            self.finish_node_func = kwargs.get("finish_node_func", None)
            self.allow_duplicate_taxon_labels = kwargs.get("allow_duplicate_taxon_labels", False)
            self.preserve_underscores = kwargs.get('preserve_underscores', False)
            self.suppress_internal_node_taxa = kwargs.get("suppress_internal_node_taxa", False)
            self.finish_node_func = None
            self.format = schema
            self._prev_taxa_block = None
            self.ncl_taxa_to_native = {}
            self._taxa_to_fill = None
            self.tree_translate_dicts = {}

        def _get_fp(self, stream):
            "Returns filepath and True if the file that `stream` refers to exists on the filesystem"
            try:
                n = stream.name
                use_ncl = os.path.exists(n)
                return n, use_ncl
            except AttributeError:
                return "", False

        def read(self, stream):
            """
            Instantiates and returns a DataSet object based on the
            NEXUS-formatted contents read from the file descriptor object
            `stream`.
            """
            n, use_ncl = self._get_fp(stream)
            if not use_ncl:
                pure_python_reader = nexusreader_py.NexusReader(
                    encode_splits = self.encode_splits,
                    rooting_interpreter = self.rooting_interpreter,
                    finish_node_func = self.finish_node_func,
                    allow_duplicate_taxon_labels = self.allow_duplicate_taxon_labels,
                    preserve_underscores = self.preserve_underscores,
                    suppress_internal_node_taxa = self.suppress_internal_node_taxa,
                    taxon_set = self.attached_taxon_set,
                    dataset = self.dataset
                )
                return pure_python_reader.read(stream)
            return self.read_filepath_into_dataset(n)

        def read_filepath_into_dataset(self, file_path):

            _LOG.debug("Creating MultiFormatReader")
            ncl_nxs_reader_handle = nclwrapper.MultiFormatReader()
            _LOG.debug("Setting MultiFormatReader's WarningOutput Level")
            ncl_nxs_reader_handle.SetWarningOutputLevel(DENDROPY_NCL_WARNING_LEVEL)
            _LOG.debug("Calling MultiFormatReader.cullIdenticalTaxaBlocks(True)")
            ncl_nxs_reader_handle.cullIdenticalTaxaBlocks(True)

            if self.dataset is None:
                self.dataset = dataobject.DataSet()

            if self.attached_taxon_set is not None and len(self.attached_taxon_set) == 0:
                self._taxa_to_fill = self.attached_taxon_set
            else:
                self._taxa_to_fill = None
            if self.attached_taxon_set is not None:
                self._register_taxa_context(ncl_nxs_reader_handle, [self.attached_taxon_set])

            _LOG.debug("Calling MultiFormatReader.ReadFilepath(%s, %s)" % (file_path, self.format))
            ncl_nxs_reader_handle.ReadFilepath(file_path, self.format)

            _LOG.debug("Calling MultiFormatReader.GetNumTaxaBlocks()")
            num_taxa_blocks = ncl_nxs_reader_handle.GetNumTaxaBlocks()
            for i in xrange(num_taxa_blocks):
                _LOG.debug("Calling MultiFormatReader.GetTaxaBlock(%d)" % i)
                ncl_tb = ncl_nxs_reader_handle.GetTaxaBlock(i)
                taxa_block = self._ncl_taxa_block_to_native(ncl_tb)
                self.dataset.add(taxa_block)

                #nab = ncl_nxs_reader_handle.GetNumAssumptionsBlocks(ncl_tb)
                #for k in xrange(nab):
                #    a = ncl_nxs_reader_handle.GetAssumptionsBlock(ncl_tb, k)
                #    cs = a.GetTaxSetNames()
                #    print "TaxSets have the names " , str(cs)

                _LOG.debug("Calling MultiFormatReader.GetNumCharactersBlocks()")
                num_char_blocks = ncl_nxs_reader_handle.GetNumCharactersBlocks(ncl_tb)
                for j in xrange(num_char_blocks):
                    _LOG.debug("Calling MultiFormatReader.GetCharactersBlock(taxablock, %d)" % j)
                    ncl_cb = ncl_nxs_reader_handle.GetCharactersBlock(ncl_tb, j)
                    char_block = self._ncl_characters_block_to_native(taxa_block, ncl_cb, ncl_nxs_reader_handle)
                    if char_block:
                        self.dataset.add(char_block)
                _LOG.debug("Calling MultiFormatReader.GetNumTreesBlocks()")
                ntrb = ncl_nxs_reader_handle.GetNumTreesBlocks(ncl_tb)
                for j in xrange(ntrb):
                    trees_block = dataobject.TreeList()
                    trees_block.taxon_set = taxa_block
                    _LOG.debug("Calling MultiFormatReader.GetTreesBlock(%d)" % j)
                    ncl_trb = ncl_nxs_reader_handle.GetTreesBlock(ncl_tb, j)
                    for k in xrange(ncl_trb.GetNumTrees()):
                        ftd = ncl_trb.GetFullTreeDescription(k)
                        tokens = ftd.GetTreeTokens()
                        rooted_flag = ftd.IsRooted()
                        t = self._ncl_tree_tokens_to_native_tree(ncl_tb, taxa_block, tokens, rooted_flag=rooted_flag)
                        if t:
                            trees_block.append(t)
                    self.dataset.add(trees_block)
            return self.dataset

        def tree_source_iter(self, stream):
            """
            Generator to iterate over trees in data file.
            Primary goal is to be memory efficient, storing no more than one tree
            at a time. Speed might have to be sacrificed for this!
            """

            n, use_ncl = self._get_fp(stream)
            if not use_ncl:
                pure_python_reader = nexusreader_py.NexusReader(
                    encode_splits = self.encode_splits,
                    rooting_interpreter = self.rooting_interpreter,
                    finish_node_func = self.finish_node_func,
                    allow_duplicate_taxon_labels = self.allow_duplicate_taxon_labels,
                    preserve_underscores = self.preserve_underscores,
                    suppress_internal_node_taxa = self.suppress_internal_node_taxa,
                    taxon_set = self.attached_taxon_set,
                    dataset = self.dataset
                )
                for tree in pure_python_reader.tree_source_iter(stream):
                    yield tree
                return

            need_tree_event = Event()
            tree_ready_event = Event()
            die_event = Event()
            ntst = NCLTreeStreamThread(n, need_tree_event=need_tree_event, ready_event=tree_ready_event, die_event=die_event, format=self.format)

            if self.dataset is None:
                self.dataset = dataobject.DataSet()
#            if self.attached_taxon_set is not None and len(self.attached_taxon_set) == 0:
#                self._taxa_to_fill = self.attached_taxon_set
#            else:
#                self._taxa_to_fill = None
#            if self.attached_taxon_set is not None:
#                self._register_taxa_context(ntst.reader, [self.attached_taxon_set])

            ncl_streamer = ntst.nts
            ntst.start()
            try:
                need_tree_event.set()
                self.curr_tree_tokens = None
                self.curr_tree = None
                while True:
                    if ntst.done:
                        break
                    tree_ready_event.wait()
                    tree_ready_event.clear()
                    ncl_taxa_block = ncl_streamer.ncl_taxa_block

                    self.curr_tree_tokens = ncl_streamer.tree_tokens
                    if self.curr_tree_tokens is None:
                        break
                    rooted_flag = ncl_streamer.rooted_flag
                    ncl_streamer.tree_tokens = None
                    need_tree_event.set()
                    self.curr_tree = self._ncl_tree_tokens_to_native_tree(ncl_taxa_block, self.attached_taxon_set, self.curr_tree_tokens, rooted_flag=rooted_flag)
                    if self.curr_tree:
                        yield self.curr_tree
                del self.curr_tree_tokens
                del self.curr_tree
            except Exception, v:
                _LOG.debug("%s" % str(v))
                die_event.set()
                need_tree_event.set()
                raise
            if ntst.exception:
                raise ntst.exception
            die_event.set()

        def _register_taxa_context(self, ncl_reader, incoming_taxa_blocks):
            if not incoming_taxa_blocks:
                return

            num_taxa_blocks = ncl_reader.GetNumTaxaBlocks()
            _LOG.debug("Registering previously read taxa blocks.  Currently %d.\nIncoming = %s" % (num_taxa_blocks, str(incoming_taxa_blocks)))
            existing_taxa_blocks = []
            for i in xrange(num_taxa_blocks):
                ncl_tb = ncl_reader.GetTaxaBlock(i)
                labels = list(ncl_tb.GetAllLabels())
                existing_taxa_blocks.append(ncl_tb, labels)

            to_add = []
            for tb in incoming_taxa_blocks:
                if tb or True:
                    found = False
                    l = [i.label for i in tb]
                    for k, v in existing_taxa_blocks:
                        if l == v:
                            found = True
                            self.ncl_taxa_to_native[k.GetInstanceIdentifierString()] = tb
                            break
                    if not found:
                        to_add.append(tb)

            for tb in to_add:
                tn = tuple([i.label for i in tb])
                _LOG.debug("RegisterTaxa(%s)" % str(tn))
                ncl_tb = ncl_reader.RegisterTaxa(tn)
                if ncl_tb:
                    iid = ncl_tb.GetInstanceIdentifierString()
                    self.ncl_taxa_to_native[iid] = tb

        def _ncl_taxa_block_to_native(self, ncl_tb):
            _LOG.debug("Converting NCL taxa block to native")
            _LOG.debug("calling NxsTaxaBlock.GetInstanceIdentifierString()")
            tbiid = ncl_tb.GetInstanceIdentifierString()
            _LOG.debug("got %s" % tbiid)
            taxa_block = self.ncl_taxa_to_native.get(tbiid)
            if taxa_block is not None:
                return taxa_block

            _LOG.debug("calling NxsTaxaBlock.GetAllLabels()")
            labels = ncl_tb.GetAllLabels()
            _LOG.debug("labels = %s" % ' '.join(labels))
            if self._taxa_to_fill is None:
                taxa_block =  dataobject.TaxonSet(labels)
            else:
                taxa_block = self._taxa_to_fill
                self._taxa_to_fill = None
                taxa_block.extend([dataobject.Taxon(label=i) for i in labels])
            self.ncl_taxa_to_native[tbiid] = taxa_block
            return taxa_block

        def _ncl_tree_tokens_to_native_tree(self, ncl_tb, taxa_block, tree_tokens, rooted_flag=None):
            if not tree_tokens:
                return None
            iid = ncl_tb.GetInstanceIdentifierString()
            if taxa_block is None:
                taxa_block = self._ncl_taxa_block_to_native(ncl_tb)
#            self.taxa_block = taxa_block
            lti = ListOfTokenIterator(tree_tokens)
            lti.tree_rooted = rooted_flag

            if iid not in self.tree_translate_dicts:
                self.tree_translate_dicts[ncl_tb] = {}
                for n, t in enumerate(taxa_block):
                    self.tree_translate_dicts[ncl_tb][str(n + 1)] = t
                    if self.encode_splits:
                        t.clade_mask = (1 << n)
            return nexustokenizer.tree_from_token_stream(lti,
                                            taxon_set=taxa_block,
                                            translate_dict=self.tree_translate_dicts[ncl_tb],
                                            encode_splits=self.encode_splits,
                                            rooting_interpreter=self.rooting_interpreter,
                                            finish_node_func=self.finish_node_func)

        def _ncl_characters_block_to_native(self, taxa_block, ncl_cb, ncl_nxs_reader_handle):
            """
            Processes a FORMAT command. Assumes that the file reader is
            positioned right after the "FORMAT" token in a FORMAT command.
            """
            raw_matrix = ncl_cb.GetRawDiscreteMatrixRef()
            if ncl_cb.IsMixedType():
                _LOG.warn("Mixed datatype character blocks are not supported in Dendropy.  Skipping...")
                return None
            char_block_type = _ncl_datatype_enum_to_dendropy(ncl_cb.GetDataType())
            mapper = ncl_cb.GetDatatypeMapperForCharRef(0)
            symbols = mapper.GetSymbols()
            state_codes_mapping = mapper.GetPythonicStateVectors()

            char_block = char_block_type()
            char_block.taxon_set = taxa_block
            if isinstance(char_block, dataobject.StandardCharacterMatrix):
                sa = dataobject.get_state_alphabet_from_symbols(
                        symbols=symbols,
                        gap_symbol='-',
                        missing_symbol='?'
                )
                char_block.state_alphabets = [sa]
                char_block.default_state_alphabet = char_block.state_alphabets[0]
            symbol_state_map = char_block.default_state_alphabet.symbol_state_map()

            ncl_numeric_code_to_state = []
            for s in symbols:
                ncl_numeric_code_to_state.append(symbol_state_map[s])
            for sc in state_codes_mapping[len(symbols):-2]:
                search = set()
                for fundamental_state in sc:
                    search.add(ncl_numeric_code_to_state[fundamental_state])
                found = False
                for sym, state in symbol_state_map.iteritems():
                    ms = state.member_states
                    if ms:
                        possible = set(ms)
                        if possible == search:
                            found = True
                            ncl_numeric_code_to_state.append(state)
                            break
                if not found:
                    raise ValueError("NCL datatype cannot be coerced into datatype because ambiguity code for %s is missing " % str(search))
            ncl_numeric_code_to_state.append(symbol_state_map['-'])
            ncl_numeric_code_to_state.append(symbol_state_map['?'])

            assert (len(raw_matrix) == len(taxa_block))
            for row_ind, taxon in enumerate(taxa_block):
                v = dataobject.CharacterDataVector(taxon=taxon)
                raw_row = raw_matrix[row_ind]
                char_block[taxon] = v
                if not self.exclude_chars:
                    for c in raw_row:
                        state = ncl_numeric_code_to_state[c]
                        v.append(dataobject.CharacterDataCell(value=state))

            #dataset.characters_blocks.append(char_block)
            supporting_exsets = False
            supporting_charset_exsets = False

            if supporting_exsets:
                s = ncl_cb.GetExcludedIndexSet()
                print "Excluded chars =", str(nclwrapper.NxsSetReader.GetSetAsVector(s))
            if supporting_charset_exsets:
                _LOG.debug("Calling MultiFormatReader.GetNumTaxaBlocks()")
                nab = ncl_nxs_reader_handle.GetNumAssumptionsBlocks(ncl_cb)
                for k in xrange(nab):
                    _LOG.debug("Calling MultiFormatReader.GetNumTaxaBlocks()")
                    a = ncl_nxs_reader_handle.GetAssumptionsBlock(ncl_cb, k)
                    cs = a.GetCharSetNames()
                    print "CharSets have the names " , str(cs)
            return char_block

    NexusReader = NCLBasedReader
