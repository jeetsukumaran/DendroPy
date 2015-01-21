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
Summarizations collections of trees, e.g., MCMC samples from a posterior
distribution, non-parametric bootstrap replicates, mapping posterior
probability, support, or frequency that splits/clades are found in the source
set of trees onto a target tree.
"""

import os
import sys
import re
import argparse
import collections
import datetime

if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
try:
    # Python 3
    import queue
except ImportError:
    # Python 2.7
    import Queue as queue
import multiprocessing

import dendropy
from dendropy.utility import cli
from dendropy.utility import constants
from dendropy.utility import messaging

##############################################################################
## Preamble

_program_name = "SumTrees"
_program_subtitle = "Phylogenetic Tree Summarization"
_program_date = "Jan 31 2015"
_program_version = "4.0.0 (%s)" % _program_date
_program_author = "Jeet Sukumaran and Mark T. Holder"
_program_contact = "jeetsukumaran@gmail.com"
_program_copyright = """\
Copyright (C) 2008-2014 Jeet Sukumaran and Mark T. Holder.
License GPLv3+: GNU GPL version 3 or later.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law."""

_program_citation = """\
Sukumaran, J and MT Holder. {prog_name}: {prog_subtitle}. {prog_version}. Available at https://github.com/jeetsukumaran/DendroPy.
""".format(prog_name=_program_name, prog_subtitle=_program_subtitle, prog_version=_program_version)

##############################################################################
## Primary Processing

class TreeProcessingWorker(multiprocessing.Process):

    def __init__(self,
            work_queue,
            results_queue,
            source_schema,
            taxon_labels,
            tree_offset,
            process_idx,
            is_source_trees_rooted,
            rooting_interpretation,
            preserve_underscores,
            ignore_edge_lengths,
            ignore_node_ages,
            use_tree_weights,
            ultrametricity_precision,
            num_processes,
            log_frequency,
            messenger,
            messenger_lock,
            ):
        multiprocessing.Process.__init__(self)
        self.work_queue = work_queue
        self.results_queue = results_queue
        self.source_schema = source_schema
        self.taxon_labels = taxon_labels
        self.taxon_namespace = dendropy.TaxonNamespace(self.taxon_labels)
        self.taxon_namespace.is_mutable = False
        self.tree_offset = tree_offset
        self.process_idx = process_idx
        self.is_source_trees_rooted = is_source_trees_rooted
        self.rooting_interpretation =rooting_interpretation
        self.preserve_underscores = preserve_underscores
        self.ignore_edge_lengths = ignore_edge_lengths
        self.ignore_node_ages = ignore_node_ages
        self.use_tree_weights = use_tree_weights
        self.ultrametricity_precision = ultrametricity_precision
        self.num_processes = num_processes
        self.log_frequency = log_frequency
        self.messenger = messenger
        self.messenger_lock = messenger_lock
        self.kill_received = False
        self.tree_array = dendropy.TreeArray(
                taxon_namespace=self.taxon_namespace,
                is_rooted_trees=self.is_source_trees_rooted,
                ignore_edge_lengths=self.ignore_edge_lengths,
                ignore_node_ages=self.ignore_node_ages,
                use_tree_weights=self.use_tree_weights,
                ultrametricity_precision=self.ultrametricity_precision,
                )

    def send_message(self, msg, level, wrap=True):
        if self.messenger is None:
            return
        if self.messenger.messaging_level > level or self.messenger.silent:
            return
        msg = "Thread %d: %s" % (self.process_idx+1, msg)
        self.messenger_lock.acquire()
        try:
            self.messenger.log(msg, level=level, wrap=wrap)
        finally:
            self.messenger_lock.release()

    def send_info(self, msg, wrap=True):
        self.send_message(msg, messaging.ConsoleMessenger.INFO_MESSAGING_LEVEL, wrap=wrap)

    def send_warning(self, msg, wrap=True):
        self.send_message(msg, messaging.ConsoleMessenger.WARNING_MESSAGING_LEVEL, wrap=wrap)

    def send_error(self, msg, wrap=True):
        self.send_message(msg, messaging.ConsoleMessenger.ERROR_MESSAGING_LEVEL, wrap=wrap)

    def run(self):
        while not self.kill_received:
            try:
                tree_source = self.work_queue.get_nowait()
            except queue.Empty:
                break
            self.send_info("Received task: '{}'".format(tree_source), wrap=False)

            self.tree_array.read_from_files(
                files=[tree_source],
                schema=self.source_schema,
                rooting=self.rooting_interpretation,
                tree_offset=self.tree_offset,
                preserve_underscores=self.preserve_underscores,
                store_tree_weights=self.use_tree_weights,
                ignore_unrecognized_keyword_arguments=True,
                )

            # with open(source, "rU") as fsrc:
            #     pass
            # for tidx, tree in enumerate(dendropy.Tree.yield_from_files(
            #         [fsrc],
            #         schema=self.schema,
            #         taxon_namespace=self.taxon_namespace,
            #         rooting=self.rooting_interpretation,
            #         store_tree_weights=self.weighted_trees)):
            #     assert tree.taxon_namespace is self.taxon_namespace
            #     if tidx >= self.tree_offset:
            #         if (self.log_frequency == 1) or (tidx > 0 and self.log_frequency > 0 and tidx % self.log_frequency == 0):
            #             self.send_info("(processing) '%s': tree at offset %d" % (source, tidx), wrap=False)
            #         self.split_distribution.count_splits_on_tree(tree, is_splits_encoded=False)
            #         if self.calc_tree_probs:
            #             self.topology_counter.count(tree,
            #                     is_splits_encoded=True)
            #     else:
            #         if (self.log_frequency == 1) or (tidx > 0 and self.log_frequency > 0 and tidx % self.log_frequency == 0):
            #             self.send_info("(processing) '%s': tree at offset %d (skipping)" % (source, tidx), wrap=False)
            #     if self.kill_received:
            #         break

            if self.kill_received:
                break
            self.send_info("Completed task: '{}'".format(tree_source), wrap=False)
        if self.kill_received:
            self.send_warning("Terminating in response to kill request")
        else:
            self.results_queue.put(self.tree_array)
            # self.result_topology_hash_map_queue.put(self.topology_counter.topology_hash_map)

class SumTrees(object):

    def __init__(self,
            is_source_trees_rooted,
            ignore_edge_lengths,
            ignore_node_ages,
            use_tree_weights,
            ultrametricity_precision,
            num_processes,
            log_frequency,
            messenger,
            ):
        self._is_source_trees_rooted = None
        self._rooting_interpretation = None
        self.is_source_trees_rooted = is_source_trees_rooted
        self.ignore_edge_lengths = ignore_edge_lengths
        self.ignore_node_ages = ignore_node_ages
        self.use_tree_weights = use_tree_weights
        self.ultrametricity_precision = ultrametricity_precision
        self.num_processes = num_processes
        self.log_frequency = log_frequency
        self.messenger = messenger

    def _get_is_source_trees_rooted(self):
        return self._is_source_trees_rooted
    def _set_is_source_trees_rooted(self, is_source_trees_rooted):
        self._is_source_trees_rooted = is_source_trees_rooted
        if self._is_source_trees_rooted is True:
            self._rooting_interpretation = "force-rooted"
        elif self._is_source_trees_rooted is False:
            self._rooting_interpretation = "force-unrooted"
        else:
            self._rooting_interpretation = "default-unrooted"
    is_source_trees_rooted = property(_get_is_source_trees_rooted, _set_is_source_trees_rooted)

    def info_message(self, msg, wrap=True):
        if self.messenger:
            self.messenger.info(msg, wrap=wrap)

    def warning_message(self, msg, wrap=True):
        if self.messenger:
            self.messenger.warning(msg, wrap=wrap)

    def error_message(self, msg, wrap=True):
        if self.messenger:
            self.messenger.error(msg, wrap=wrap)

    def process_trees(self,
            tree_sources,
            schema,
            taxon_namespace=None,
            tree_offset=0,
            preserve_underscores=False,
            ):
        if self.num_processes is None or self.num_processes <= 1:
            tree_array = self.serial_process_trees(
                    tree_sources=tree_sources,
                    schema=schema,
                    taxon_namespace=taxon_namespace,
                    tree_offset=tree_offset,
                    preserve_underscores=preserve_underscores,
                    )
        else:
            tree_array = self.parallel_process_trees(
                    tree_sources=tree_sources,
                    schema=schema,
                    taxon_namespace=taxon_namespace,
                    tree_offset=tree_offset,
                    preserve_underscores=preserve_underscores,
                    )
        return tree_array

    def serial_process_trees(self,
            tree_sources,
            schema,
            taxon_namespace=None,
            tree_offset=0,
            preserve_underscores=False,
            ):
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        self.info_message("Running in serial mode")
        tree_array = dendropy.TreeArray(
                taxon_namespace=taxon_namespace,
                is_rooted_trees=self.is_source_trees_rooted,
                ignore_edge_lengths=self.ignore_edge_lengths,
                ignore_node_ages=self.ignore_node_ages,
                use_tree_weights=self.use_tree_weights,
                ultrametricity_precision=self.ultrametricity_precision,
                )
        if not self.log_frequency:
            tree_array.read_from_files(
                files=tree_sources,
                schema=schema,
                rooting=self._rooting_interpretation,
                tree_offset=tree_offset,
                store_tree_weights=self.use_tree_weights,
                preserve_underscores=preserve_underscores,
                ignore_unrecognized_keyword_arguments=True,
                )
        else:
            def _log_progress(current_tree_offset, aggregate_tree_idx):
                if (
                        self.messenger is not None
                        and self.log_frequency == 1
                        or current_tree_offset == tree_offset
                        or (aggregate_tree_idx > 0 and self.log_frequency > 0 and aggregate_tree_idx % self.log_frequency == 0)
                        ):
                    if current_tree_offset >= tree_offset:
                        # coda = " (processing: {trees_stored} trees processed out of {trees_read} trees read]".format(
                        #         trees_read=aggregate_tree_idx,
                        #         trees_stored=len(tree_array),
                        #         )
                        coda = " (processing)"
                    else:
                        coda = " (skipping)"
                    self.info_message("'{source_name}': tree at offset {current_tree_offset}{coda}".format(
                        source_name=source_name,
                        current_tree_offset=current_tree_offset,
                        coda=coda,
                        ), wrap=False)
            tree_yielder = dendropy.Tree.yield_from_files(
                    tree_sources,
                    schema=schema,
                    taxon_namespace=taxon_namespace,
                    store_tree_weights=self.use_tree_weights,
                    preserve_underscores=preserve_underscores,
                    rooting=self._rooting_interpretation,
                    ignore_unrecognized_keyword_arguments=True,
                    )
            current_source_index = None
            current_tree_offset = None
            for aggregate_tree_idx, tree in enumerate(tree_yielder):
                current_yielder_index = tree_yielder.current_file_index
                if current_yielder_index != current_source_index:
                    current_source_index = current_yielder_index
                    current_tree_offset = 0
                    source_name = tree_yielder.current_file_name
                    if source_name is None:
                        source_name = "<stdin>"
                    self.info_message("Processing {} of {}: '{}'".format(current_source_index+1, len(tree_sources), source_name), wrap=False)
                if current_tree_offset >= tree_offset:
                    tree_array.add_tree(tree=tree, is_bipartitions_updated=False)
                    _log_progress(current_tree_offset, aggregate_tree_idx)
                    # if len(tree_array._split_distribution.tree_rooting_types_counted) > 1:
                    #     mixed_tree_rootings_in_source_error(messenger)
                else:
                    _log_progress(current_tree_offset, aggregate_tree_idx)
                current_tree_offset += 1
        return tree_array

    def parallel_process_trees(self,
            tree_sources,
            schema,
            tree_offset=0,
            preserve_underscores=False,
            taxon_namespace=None,
            ):
        # describe
        self.info_message("Running in multiprocessing mode (up to {} processes)".format(self.num_processes))
        # taxon definition
        if taxon_namespace is not None:
            self.info_message("Using taxon names provided by user")
        else:
            tdfpath = tree_sources[0]
            self.info_message("Pre-loading taxon names based on first tree in source '{}'".format(tdfpath))
            taxon_namespace = self.discover_taxa(tdfpath, schema, preserve_underscores=preserve_underscores)
        taxon_labels = [t.label for t in taxon_namespace]
        self.info_message("{} taxa defined: [{}]".format(
                len(taxon_labels),
                ', '.join(["'{}'".format(t) for t in taxon_labels]),
                ))
        # load up queue
        self.info_message("Creating work queue")
        work_queue = multiprocessing.Queue()
        for f in tree_sources:
            work_queue.put(f)

        # launch processes
        self.info_message("Launching worker processes")
        tree_array_queue = multiprocessing.Queue()
        messenger_lock = multiprocessing.Lock()
        for idx in range(self.num_processes):
            tree_processing_worker = TreeProcessingWorker(
                    work_queue=work_queue,
                    results_queue=tree_array_queue,
                    source_schema=schema,
                    taxon_labels=taxon_labels,
                    tree_offset=tree_offset,
                    process_idx=idx,
                    is_source_trees_rooted=self.is_source_trees_rooted,
                    rooting_interpretation=self._rooting_interpretation,
                    preserve_underscores=preserve_underscores,
                    ignore_edge_lengths=self.ignore_edge_lengths,
                    ignore_node_ages=self.ignore_node_ages,
                    use_tree_weights=self.use_tree_weights,
                    ultrametricity_precision=self.ultrametricity_precision,
                    num_processes=self.num_processes,
                    messenger=self.messenger,
                    messenger_lock=messenger_lock,
                    log_frequency=self.log_frequency)
            tree_processing_worker.start()

        # collate results
        result_count = 0
        master_tree_array = dendropy.TreeArray(
                taxon_namespace=taxon_namespace,
                is_rooted_trees=self.is_source_trees_rooted,
                ignore_edge_lengths=self.ignore_edge_lengths,
                ignore_node_ages=self.ignore_node_ages,
                use_tree_weights=self.use_tree_weights,
                ultrametricity_precision=self.ultrametricity_precision,
                )
        while result_count < self.num_processes:
            worker_tree_array = tree_array_queue.get()
            master_tree_array.update(worker_tree_array)
            result_count += 1
        self.info_message("Recovered results from all worker processes")
        return master_tree_array

    def discover_taxa(self,
            treefile,
            schema,
            preserve_underscores):
        """
        Reads first tree in treefile, and assumes that is sufficient to populate a
        taxon set object fully, which it then returns.
        """
        for tree in dendropy.Tree.yield_from_files([treefile],
                schema=schema,
                preserve_underscores=preserve_underscores,
                ignore_unrecognized_keyword_arguments=True,
                ):
            return tree.taxon_namespace

##############################################################################
## Preprocessing

def preprocess_tree_sources(args, messenger):
    tree_sources = []
    for fpath in args.tree_sources:
        if fpath == "-":
            if args.source_format is None:
                messenger.error("Format of source trees must be specified using '-i' or '--source-format' flags when reading trees from standard input")
                sys.exit(1)
            elif args.source_format.lower() == "nexus/newick":
                messenger.error("The 'nexus/newick' format is not supported when reading trees from standard input")
                sys.exit(1)
            if len(args.tree_sources) > 1:
                messenger.error("Cannot specify multiple sources when reading from standard input")
            return []
        else:
            fpath = os.path.expanduser(os.path.expandvars(fpath))
            if not os.path.exists(fpath):
                missing_msg = "Support file not found: '{}'".format(fpath)
                if args.ignore_missing_support:
                    messenger.warning(missing_msg)
                else:
                    messenger.error(missing_msg )
                    messenger.error("Terminating due to missing support files. "
                            "Use the '--ignore-missing-support' option to continue even "
                            "if some files are missing.")
                    sys.exit(1)
            else:
                tree_sources.append(fpath)
    if len(tree_sources) == 0:
        messenger.error("No valid sources of input trees specified. "
                + "Please provide the path to at least one (valid and existing) file "
                + "containing tree samples to summarize.")
        sys.exit(1)
    return tree_sources

##############################################################################
## Front-End

def citation(args):
    show_splash(dest=sys.stdout)
    sys.exit(0)

def show_splash(dest=None):
    if dest is None:
        dest = sys.stderr
    cli.show_splash(
            prog_name=_program_name,
            prog_subtitle=_program_subtitle,
            prog_version=_program_version,
            prog_author=_program_author,
            prog_copyright=_program_copyright,
            include_citation=True,
            include_copyright=False,
            additional_citations=[_program_citation],
            dest=dest,
            )
    dest.write("\n")

def print_usage_examples(dest=None):
    if dest is None:
        dest = sys.stdout
    examples = ("""\
Summarize a set of tree files using a 95% rule consensus tree, with support for
clades expressed as proportions (posterior probabilities) on internal node
labels and branch lengths the mean across all trees, dropping the first 200
trees in each file as a burn-in, and saving the result to "``result.tre``"::

    $ sumtrees.py \\
            --summary-tree consensus \\
            --min-clade-freq=0.95 \\
            --edges mean-length \\
            --burnin=200 \\
            --support-as-labels \\
            --output=result.tre \\
            treefile1.tre treefile2.tre treefile3.tre

To use a different type of summary tree, e.g., the tree that maximizes the
product of posterior probabilities, you can specify 'mct' for the
'--summary-tree' option:

    $ sumtrees.py \\
            --summary-tree mct \\
            --min-clade-freq=0.95 \\
            --edges mean-length \\
            --burnin=200 \\
            --support-as-labels \\
            --output=result.tre \\
            treefile1.tre treefile2.tre treefile3.tre

If the input trees are ultrametric and you want to set the node ages to the
median node age, specify the '--ultrametic' option and set the '--edges'
argument to 'median-age':

    $ sumtrees.py \\
            --ultrametric \\
            --summary-tree mct \\
            --edges median-age \\
            --burnin=200 \\
            --output=result.tre \\
            treefile1.tre treefile2.tre treefile3.tre

Calculate support for nodes on a specific tree, "``best.tre``" as given by a
set of tree files, with support reported as percentages rounded to integers,
and saving the result to "``result.tre``"::

    $ sumtrees.py \\
            --target=best.tre \\
            --decimals=0 \\
            --percentages \\
            --output=result.tre \\
            treefile1.tre treefile2.tre treefile3.tre

""")
    dest.write(examples + "\n")

def print_description(dest=None):
    import site
    if dest is None:
        dest = sys.stdout
    fields = collections.OrderedDict()
    fields["DendroPy version"] = dendropy.description()
    fields["DendroPy location"] = dendropy.homedir()
    fields["Python version"] = sys.version.replace("\n", "")
    fields["Python executable"] = sys.executable
    fields["Python site packages"] = site.getsitepackages()
    max_fieldname_len = max(len(fieldname) for fieldname in fields)
    for fieldname, fieldvalue in fields.items():
        dest.write("{fieldname:{fieldnamewidth}}: {fieldvalue}\n".format(
            fieldname=fieldname,
            fieldnamewidth=max_fieldname_len + 2,
            fieldvalue=fieldvalue))

def main():
    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=cli.CustomFormatter,
            add_help=False,
            )
    source_options = parser.add_argument_group("Source Options")
    source_options.add_argument("-i","--source-format",
            metavar="FORMAT",
            default=None,
            choices=["nexus/newick", "nexus", "newick", "phylip", "nexml"],
            help="Format of the source trees (defaults to handling either NEXUS or NEWICK through inspection; it is more efficient to explicitly specify the format if it is known).")
    source_options.add_argument("-b", "--burnin",
            type=int,
            default=0,
            help="Number of trees to skip from the beginning of *each* tree file when counting support (default: %(default)s).")
    source_options.add_argument("--rooted",
            dest="is_source_trees_rooted",
            action="store_true",
            default=None,
            help="Treat source trees as rooted.")
    source_options.add_argument("--unrooted",
            dest="is_source_trees_rooted",
            action="store_false",
            default=None,
            help="Treat source trees as unrooted.")
    source_options.add_argument("-u", "--ultrametric",
            dest="is_source_trees_ultrametric",
            action="store_true",
            default=None,
            help="Assume source trees are ultrametric (implies '--rooted'; will result in node ages being summarized; will result in error if trees are not ultrametric).")
    source_options.add_argument("-y", "--ultrametricity-precision",
            action="store_true",
            default=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            help="Precision to use when validating ultrametricity (default: %(default)s; specify '0' to disable validation).")
    source_options.add_argument("--weighted-trees",
            action="store_true",
            default=False,
            help="Use weights of trees (as indicated by '[&W m/n]' comment token) to weight contribution of splits found on each tree to overall split frequencies.")
    source_options.add_argument("--preserve-underscores",
            action="store_true",
            default=False,
            help="Do not convert unquoted/unprotected underscores to spaces when reading NEXUS/NEWICK format trees.")
    source_options.add_argument("--taxon-name-file",
            metavar="FILEPATH",
            default=None,
            help=(
                "Path to file listing all the taxon names or labels that"
                " will be found across the entire set of source trees."
                " This file should be a plain text file with a single"
                " name list on each line. This file is only read when the"
                " '-m' or '--multiprocessing' operation is specific. When"
                " parallel processing using the '-m'/'--multiprocessing'"
                " option, all taxon names need to be defined in advance"
                " of any actual tree processing. By default this is done"
                " by reading the first tree in the first tree source,"
                " extracting the taxon names. At best, this is wasteful,"
                " as it involves an extraneous reading of the tree."
                " At worst, this can be wasteful AND errorneous, if the"
                " first tree does not contain all the taxa. Explicitly"
                " specifying the taxon names can avoid these issues."
                ))

    summary_tree_options = parser.add_argument_group("Target Tree Topology Options")
    summary_tree_options.add_argument(
            "-s", "--summary-tree-target",
            default=None,
            metavar="{consensus,mct,msct}",
            help="\n".join((
                "R}Map support and other information from the ",
                "source trees to a topology summarized from ",
                "the source trees under different criteria:  ",
                "- 'consensus' : consensus tree [DEFAULT];",
                "                The minimum frequency threshold can",
                "                be specified using the '-f' or",
                "                '--min-clade-freq' flags.",
                "- 'mct'       : maximum clade credibility tree;",
                "                Tree from the source set that ",
                "                maximizes the *product* of clade",
                "                posterior probabilities.",
                "- 'msct'      : maximum *sum* clade credibility tree;",
                "                Tree from the source set that ",
                "                maximizes the *sum* of clade ",
                "                posterior probabilities.",
                )))
    summary_tree_options.add_argument(
            "-t", "--target-tree-filepath",
            default=None,
            metavar="FILE",
            help=(
                  "Map support and other information from the "
                  "source trees to a topology or topologies "
                  "given by the tree(s) described in FILE."
                 ))

    target_tree_supplemental_options = parser.add_argument_group("Target Tree Supplemental Options")
    target_tree_supplemental_options.add_argument("-f", "--min-consensus-freq",
            type=float,
            default=constants.GREATER_THAN_HALF,
            metavar="#.##",
            help=(
                "If using a consensus tree summarization strategy, then "
                "this is the minimum frequency or probability for a clade "
                "or a split to be included in the resulting tree "
                "(default: > 0.5)."
                ))

    target_tree_rooting_options = parser.add_argument_group("Target Tree Rooting Options")
    target_tree_rooting_options.add_argument("--root-target-at-outgroup",
            metavar="TAXON-LABEL",
            default=None,
            help="Root target tree(s) using specified taxon as outgroup.")
    target_tree_rooting_options.add_argument("--root-target-at-midpoint",
            action="store_true",
            default=None,
            help="Root target tree(s) at midpoint.")

    edge_summarization_options = parser.add_argument_group("Target Tree Edge Options")
    edge_summarization_choices = ["mean-length", "median-length", "mean-age", "median-age", "keep", "clear"]
    edge_summarization_options.add_argument("-e", "--edges",
            dest="edge_summarization",
            # metavar="<%s>" % ("".join(edge_summarization_choices)),
            choices=edge_summarization_choices,
            default=None,
            help="\n".join((
                "R} Set the edge lengths of the summary target tree(s):",
                "- 'mean-age'       : Edge lengths will be adjusted so that",
                "                     the age of subtended nodes will be equal",
                "                     to the mean age of the corresponding",
                "                     split or clade in the source trees. This",
                "                     option requires that the source trees",
                "                     are ultrametric (i.e., '-u' or",
                "                     '--ultrametric' must be specified).",
                "                     If no external summary tree targets",
                "                     are specified and the input source",
                "                     trees are specified as ultrametric,",
                "                     then this is the default option.",
                "- 'mean-length'    : Edge lengths will be set to the mean",
                "                     of the lengths of the corresponding",
                "                     split or clade in the source trees.",
                "                     If no external summary tree targets",
                "                     are specified and the input source",
                "                     trees are not ultrametric, then this",
                "                     is the default option.",
                "- 'keep'           : Do not change the existing edge lengths.",
                "                     This is the default if target tree(s) are",
                "                     sourced from an external file using the",
                "                     '-t' or '--target' option",
                "- 'median-age'     : Edge lengths will be adjusted so that",
                "                     the age of subtended nodes will be equal",
                "                     to the median age of the corresponding",
                "                     split or clade in the source trees. This",
                "                     option requires that the source trees",
                "                     are ultrametric (i.e., '-u' or",
                "                     '--ultrametric' must be specified).",
                "- 'median-length'  : Edge lengths will be set to the median",
                "                     of the lengths of the corresponding",
                "                     split or clade in the source trees.",
                "                     If no external summary tree targets",
                "                     are specified and the input source",
                "                     trees are not ultrametric, then this",
                "                     is the default option.",
                "- 'clear'          : Edge lengths will be cleared from the",
                "                     target trees if they are present.",
                )))
    edge_summarization_options.add_argument("--collapse-negative-edges",
            action="store_true",
            default=False,
            help="(If setting edge lengths) force parent node ages to be at least as old as its oldest child when summarizing node ages.")


    support_summarization_options = parser.add_argument_group("Support Summarization Options")
    support_summarization_options.add_argument("-l","--support-as-labels",
            action="store_const",
            dest="support_annotation_target",
            default=1,
            const=1,
            help="Indicate split/clade/branch support as internal node labels.")
    support_summarization_options.add_argument("-v","--support-as-lengths",
            action="store_const",
            dest="support_annotation_target",
            default=1,
            const=2,
            help="Indicate split/clade/branch support as branch lengths.")
    support_summarization_options.add_argument("-x","--no-support-labels",
            action="store_const",
            dest="support_annotation_target",
            default=1,
            const=0,
            help=("Do not indicate support with internal node labels or edge"
                  " lengths. Support will still be indicated as node metadata "
                  " annotations unless '--no-annotations' is specified "
                 ))
    support_summarization_options.add_argument("-p", "--percentages",
            action="store_true",
            dest="support_as_percentages",
            default=False,
            help="Indicate branch support as percentages (otherwise, will report as proportions by default).")
    support_summarization_options.add_argument("-d", "--decimals",
            dest="support_label_decimals",
            type=int,
            metavar="#",
            default=8,
            help="Number of decimal places in indication of support values (default: %(default)s).")

    other_summarization_options = parser.add_argument_group("Other Summarization Options")
    #other_summarization_options.add_argument("--with-node-ages",
    #        action="store_true",
    #        dest="summarize_node_ages",
    #        default=None,
    #        help="summarize node ages as well as edge lengths (implies '--rooted' and '--ultrametric'; automatically enabled if '--ultrametric' is specified; will result in error if trees are not ultrametric)")
    other_summarization_options.add_argument("--trprobs", "--calc-tree-probabilities",
            dest="trprobs_filepath",
            default=None,
            metavar="FILEPATH",
            help=("If specified, a file listing tree (topologies) and the "
                 "frequencies of their occurrences will be saved to FILEPATH."))
    other_summarization_options.add_argument("--extract-edges",
            dest="split_edge_map_filepath",
            default=None,
            metavar="FILEPATH",
            help=("if specified, a tab-delimited file of splits and their edge "
                  "lengths across input trees will be saved to FILEPATH"))
    other_summarization_options.add_argument("--no-summarize-node-ages",
            action="store_false",
            dest="summarize_node_ages",
            default=None,
            help="Do not calculate/summarize node ages, even if '--ultrametric' is specified.")
    other_summarization_options.add_argument("--no-summarize-edge-lengths",
            action="store_false",
            dest="summarize_edge_lengths",
            default=None,
            help="Do not summarize edge lengths.")
    other_summarization_options.add_argument("--no-annotations",
            action="store_true",
            default=False,
            help="Do NOT annotate nodes and edges with any summarization information metadata.")

    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument("-o","--output-tree-filepath",
            metavar="FILEPATH",
            default=None,
            help="Path to output file (if not given, will print to standard output).")
    output_options.add_argument("-F","--output-format",
            metavar="FORMAT",
            default="nexus",
            choices=["nexus", "newick", "phylip", "nexml"],
            help="Format of the output tree file (default: '%(default)s').")
    output_options.add_argument("--no-taxa-block",
            action="store_false",
            default=True,
            help="When writing NEXUS format output, do not include a taxa block in the output treefile (otherwise will create taxa block by default).")
    output_options.add_argument("--no-meta-comments",
            action="store_false",
            dest="include_meta_comments",
            default=True,
            help="Do not include initial file comment annotating details of scoring operation")
    output_options.add_argument("-c", "--additional-comments",
            action="store",
            dest="additional_comments",
            default=None,
            help="Additional comments to be added to the summary file.")
    output_options.add_argument("-r", "--replace",
            action="store_true",
            dest="replace",
            default=False,
            help="Replace/overwrite output file without asking if it already exists.")

    run_options = parser.add_argument_group("Program Run Options")
    run_options.add_argument("-m", "--multiprocessing",
            dest="multiprocess",
            metavar="NUM-PROCESSES",
            default=None,
            help=(
                 "Run in parallel mode with up to a maximum of NUM-PROCESSES processes "
                 "(specify 'max' or '#' to run in as many processes as there are cores on the "
                 "local machine)."
                 ))
    run_options.add_argument("-g", "--log-frequency",
            type=int,
            metavar="LOG-FREQUENCY",
            default=500,
            help="Tree processing progress logging frequency (default: %(default)s; set to 0 to suppress).")
    run_options.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Suppress ALL logging, progress and feedback messages.")
    run_options.add_argument("--ignore-missing-support",
            action="store_true",
            default=False,
            help="Ignore missing support tree files (note that at least one must exist).")

    information_options = parser.add_argument_group("Program Information Options")
    information_options.add_argument("-h", "--help",
            action="store_true",
            default=False,
            help="Show help information for program and exit.")
    information_options.add_argument("--citation",
            action="store_true",
            default=False,
            help="Show citation information for program and exit.")
    information_options.add_argument("--usage-examples",
            action="store_true",
            default=False,
            help="Show usage examples of program and exit.")
    information_options.add_argument("--describe",
            action="store_true",
            default=False,
            help="Show information regarding your DendroPy and Python installations and exit.")

    parser.add_argument("tree_sources",
            nargs="*",
            metavar="TREE-FILEPATH",
            help= (
                "Source(s) of trees to summarize. At least one valid"
                " source of trees must be provided. Use '-' to specify"
                " reading from standard input (note that this requires"
                " the input file format to be explicitly set using"
                " the '-i' or '--source-format' option)."
            ))

    args = parser.parse_args()

    if args.citation:
        citation(args)

    if not args.quiet:
        show_splash()

    if args.help:
        parser.print_help(sys.stdout)
        sys.exit(0)

    if args.usage_examples:
        print_usage_examples(sys.stdout)
        sys.exit(0)

    if args.describe:
        print_description(sys.stdout)
        sys.exit(0)

    ######################################################################
    ## Set up messenger

    if args.quiet:
        messaging_level = messaging.ConsoleMessenger.ERROR_MESSAGING_LEVEL
    else:
        messaging_level = messaging.ConsoleMessenger.INFO_MESSAGING_LEVEL
    messenger = messaging.ConsoleMessenger(name="SumTrees", messaging_level=messaging_level)

    ######################################################################
    ## Support Files

    if len(args.tree_sources) == 0:
        parser.print_usage()
        sys.stdout.write("\n")
        sys.stdout.write("Type 'sumtrees.py --help' for details on usage.\n")
        sys.stdout.write("Type 'sumtrees.py --usage-examples' for examples of usage.\n")
        sys.exit(0)

    tree_sources = preprocess_tree_sources(args, messenger)
    if tree_sources is None or len(tree_sources) == 0:
        tree_sources = [sys.stdin]
        messenger.info("Reading trees from standard input")
    else:
        messenger.info("{} source(s) to be processed".format(len(tree_sources)))

    ######################################################################
    ## Target Validation

    if args.summary_tree_target is not None and args.target_tree_filepath is not None:
        messenger.error("Cannot specify both '-s'/'--summary-tree-target' and '-t'/'--target-tree-filepath' simultaneously")
    elif args.target_tree_filepath is not None:
        target_tree_filepath = os.path.expanduser(os.path.expandvars(args.target_tree_filepath))
        if not os.path.exists(target_tree_filepath):
            messenger.error("Target tree file not found: '{}'".format(target_tree_filepath))
            sys.exit(1)
    elif args.summary_tree_target is not None:
        target_tree_filepath = None

    ######################################################################
    ## Tree Ultrametricity and Rooting State

    if args.edge_summarization == "mean-age" or args.edge_summarization == "median-age":
        if args.is_source_trees_ultrametric is None:
            messenger.info("Edge summarization strategy '{}' requires ultrametric source trees: assuming source trees are ultrametric".format(args.edge_summarization))
            args.is_source_trees_ultrametric = True
        elif args.is_source_trees_ultrametric is False:
            messenger.error("Edge summarization strategy '{}' requires ultrametric source trees, but source trees are specified as non-ultrametric".format(args.edge_summarization))
            sys.exit(1)
        if args.summarize_node_ages is False:
            messenger.error("Edge summarization strategy '{}' requires node ages to be summarized, but '--no-summarize-node-ages' specified".format(args.edge_summarization))
            sys.exit(1)
        args.summarize_node_ages = True
    if args.edge_summarization == "mean-lengths" or args.edge_summarization == "median-lengths":
        if args.summarize_edge_lengths is False:
            messenger.error("Edge summarization strategy '{}' requires edge lengths to be summarized, but '--no-summarize-edge-lengths' specified".format(args.edge_summarization))
            sys.exit(1)

    if args.is_source_trees_ultrametric:
        if args.is_source_trees_rooted is False:
            messenger.error("Ultrametric source trees imply rooted trees, but source trees are explicitly specified as unrooted".format(args.edge_summarization))
            sys.exit(1)
        args.is_source_trees_rooted = True
        if args.summarize_node_ages is None:
            args.summarize_node_ages = True
    else:
        if args.summarize_node_ages is True:
            if args.is_source_trees_ultrametric is None:
                messenger.info("Summarization of node ages ('--summarize-node-ages') require ultrametric trees: assuming source trees are ultrametric")
                args.is_source_trees_ultrametric = True
                if args.is_source_trees_rooted is False:
                    messenger.error("Ultrametric source trees imply rooted trees, but source trees are explicitly specified as unrooted".format(args.edge_summarization))
                    sys.exit(1)
                args.is_source_trees_rooted = True
            elif args.is_source_trees_ultrametric is False:
                messenger.error("Summarization of node ages ('--summarize-node-ages') require ultrametric trees, but source trees are specified as non-ultrametric")
                sys.exit(1)
        else:
            args.sumamrize_node_ages = False

    ######################################################################
    ## Output File Setup

    if args.output_tree_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(args.output_tree_filepath))
        if cli.confirm_overwrite(filepath=output_fpath, replace_without_asking=args.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)

    if args.trprobs_filepath:
        trprobs_filepath = os.path.expanduser(os.path.expandvars(args.trprobs_filepath))
        if cli.confirm_overwrite(filepath=trprobs_filepath, replace_without_asking=args.replace):
            trprobs_dest = open(trprobs_filepath, "w")
        else:
            sys.exit(1)
        args.calc_tree_probs = True
    else:
        trprobs_dest = None
        args.calc_tree_probs = False

    if args.split_edge_map_filepath:
        split_edge_map_filepath = os.path.expanduser(os.path.expandvars(args.split_edge_map_filepath))
        if confirm_overwrite(filepath=split_edge_map_filepath, replace_without_asking=args.replace):
            cli.split_edge_map_dest = open(split_edge_map_filepath, "w")
        else:
            sys.exit(1)
    else:
        split_edge_map_dest = None

    ######################################################################
    ## Multiprocessing Setup

    if len(tree_sources) > 1 and args.multiprocess is not None:
        num_cpus = multiprocessing.cpu_count()
        if (
                args.multiprocess.lower() == "max"
                or args.multiprocess == "#"
                or args.multiprocess == "*"
            ):
            num_processes = num_cpus
        elif args.multiprocess == "@":
            num_processes = len(tree_sources)
        else:
            try:
                num_processes = int(args.multiprocess)
            except ValueError:
                messenger.error("'{}' is not a valid number of processes (must be a positive integer)".format(args.multiprocess))
                sys.exit(1)
            if num_processes > num_cpus:
                messenger.warning("Number of requested processes ({}) exceeds number of CPU's ({})".format(num_processes, num_cpus))
        if num_processes <= 0:
            messenger.error("Maximum number of processes set to {}: cannot run SumTrees with less than 1 process".format(num_processes))
            sys.exit(1)
    else:
        if args.multiprocess is not None and args.multiprocess > 1:
            messenger.info("Number of valid sources is less than 2: forcing serial processing")
        num_processes = 1

    ######################################################################
    ## Format

    if args.source_format is None:
        schema = "nexus/newick"
    else:
        schema = args.source_format.lower()

    ######################################################################
    ## Taxon Discovery

    if args.taxon_name_file is not None:
        with open(os.path.expanduser(os.path.expandvars(args.taxon_name_file)), "r") as tnf:
            taxon_labels = [name.strip() for name in tnf.read().split("\n") if name]
            taxon_labels = [name for name in taxon_labels if name]
        taxon_namespace = dendropy.TaxonNamespace(taxon_labels)
    else:
        taxon_namespace = None

    ######################################################################
    ## Main Work

    start_time = datetime.datetime.now()
    sumtrees = SumTrees(
            is_source_trees_rooted=args.is_source_trees_rooted,
            ignore_edge_lengths=not args.summarize_edge_lengths,
            ignore_node_ages=not args.summarize_node_ages,
            use_tree_weights=args.weighted_trees,
            ultrametricity_precision=args.ultrametricity_precision,
            num_processes=num_processes,
            log_frequency=args.log_frequency if not args.quiet else 0,
            messenger=messenger,
            )
    tree_array = sumtrees.process_trees(
            tree_sources=tree_sources,
            schema=schema,
            taxon_namespace=taxon_namespace,
            tree_offset=args.burnin,
            preserve_underscores=args.preserve_underscores,
            )

    ######################################################################
    ## Post-Processing

    t = tree_array.consensus_tree()
    # print(t.as_string("nexus"))
    print(len(tree_array))

if __name__ == '__main__':
    main()
