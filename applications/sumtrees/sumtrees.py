#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
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
Summarizations collections of trees, e.g., MCMC samples from a posterior
distribution, non-parametric bootstrap replicates, mapping posterior
probability, support, or frequency that splits/clades are found in the source
set of trees onto a target tree.
"""

import os
import sys
import re
import getpass
import argparse
import collections
import datetime
import platform
import socket
import math
import csv
import json

try:
    # Python 3
    import queue
except ImportError:
    # Python 2.7
    import Queue as queue
import multiprocessing

import dendropy
if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
from dendropy.utility import cli
from dendropy.utility import constants
from dendropy.utility import error
from dendropy.utility import messaging
from dendropy.utility import timeprocessing
from dendropy.utility import bitprocessing
from dendropy.utility import textprocessing

##############################################################################
## Preamble

_program_name = "SumTrees"
_program_subtitle = "Phylogenetic Tree Summarization"
_program_version = dendropy.__version__
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
## Primary Analyzing

def _read_into_tree_array(
        tree_array,
        tree_sources,
        schema,
        taxon_namespace,
        rooting,
        tree_offset,
        use_tree_weights,
        preserve_underscores,
        info_message_func,
        error_message_func,
        log_frequency,
        debug_mode,
        ):
    if not log_frequency:
        tree_array.read_from_files(
            files=tree_sources,
            schema=schema,
            rooting=rooting,
            tree_offset=tree_offset,
            store_tree_weights=use_tree_weights,
            preserve_underscores=preserve_underscores,
            ignore_unrecognized_keyword_arguments=True,
            )
    else:
        def _log_progress(source_name, current_tree_offset):
            if (
                    info_message_func is not None
                    and (
                        (log_frequency == 1)
                        or (tree_offset > 0 and current_tree_offset == tree_offset)
                        or (current_tree_offset >= 0 and log_frequency > 0 and (current_tree_offset % log_frequency) == 0)
                        )
                    ):
                if current_tree_offset >= tree_offset:
                    coda = " (analyzing)"
                else:
                    coda = " (burning-in)"
                info_message_func("'{source_name}': tree at offset {current_tree_offset}{coda}".format(
                    source_name=source_name,
                    current_tree_offset=current_tree_offset,
                    coda=coda,
                    ), wrap=False)
        tree_yielder = dendropy.Tree.yield_from_files(
                tree_sources,
                schema=schema,
                taxon_namespace=taxon_namespace,
                store_tree_weights=use_tree_weights,
                preserve_underscores=preserve_underscores,
                rooting=rooting,
                ignore_unrecognized_keyword_arguments=True,
                )
        current_source_index = None
        current_tree_offset = None
        try:
            for aggregate_tree_idx, tree in enumerate(tree_yielder):
                current_yielder_index = tree_yielder.current_file_index
                if current_yielder_index != current_source_index:
                    current_source_index = current_yielder_index
                    current_tree_offset = 0
                    source_name = tree_yielder.current_file_name
                    if source_name is None:
                        source_name = "<stdin>"
                    if len(tree_sources) > 1:
                        info_message_func("Analyzing {} of {}: '{}'".format(current_source_index+1, len(tree_sources), source_name), wrap=False)
                    else:
                        info_message_func("Analyzing: '{}'".format(source_name), wrap=False)
                if current_tree_offset >= tree_offset:
                    tree_array.add_tree(tree=tree, is_bipartitions_updated=False)
                    _log_progress(source_name, current_tree_offset)
                else:
                    _log_progress(source_name, current_tree_offset)
                current_tree_offset += 1
        except (Exception, KeyboardInterrupt) as e:
            if debug_mode and not isinstance(e, KeyboardInterrupt):
                raise
            e.exception_tree_source_name = tree_yielder.current_file_name
            e.exception_tree_offset = current_tree_offset
            raise e

class TreeAnalysisWorker(multiprocessing.Process):

    def __init__(self,
            name,
            work_queue,
            results_queue,
            source_schema,
            taxon_labels,
            tree_offset,
            is_source_trees_rooted,
            preserve_underscores,
            ignore_edge_lengths,
            ignore_node_ages,
            use_tree_weights,
            ultrametricity_precision,
            taxon_label_age_map,
            log_frequency,
            messenger,
            messenger_lock,
            debug_mode,
            ):
        multiprocessing.Process.__init__(self, name=name)
        self.work_queue = work_queue
        self.results_queue = results_queue
        self.source_schema = source_schema
        self.taxon_labels = taxon_labels
        self.taxon_namespace = dendropy.TaxonNamespace(self.taxon_labels)
        self.taxon_namespace.is_mutable = False
        self.tree_offset = tree_offset
        self.is_source_trees_rooted = is_source_trees_rooted
        self.rooting_interpretation = dendropy.get_rooting_argument(is_rooted=self.is_source_trees_rooted)
        self.preserve_underscores = preserve_underscores
        self.ignore_edge_lengths = ignore_edge_lengths
        self.ignore_node_ages = ignore_node_ages
        self.use_tree_weights = use_tree_weights
        self.ultrametricity_precision = ultrametricity_precision
        self.taxon_label_age_map = taxon_label_age_map
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
                taxon_label_age_map=self.taxon_label_age_map,
                )
        self.tree_array.worker_name = self.name
        self.num_tasks_received = 0
        self.num_tasks_completed = 0
        self.debug_mode = debug_mode

    def send_message(self, msg, level, wrap=True):
        if self.messenger is None:
            return
        if self.messenger.messaging_level > level or self.messenger.silent:
            return
        msg = "{}: {}".format(self.name, msg)
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
            self.num_tasks_received += 1
            # self.send_info("Received task {task_count}: '{task_name}'".format(
            self.send_info("Received task: '{task_name}'".format(
                task_count=self.num_tasks_received,
                task_name=tree_source), wrap=False)
            # self.tree_array.read_from_files(
            #     files=[tree_source],
            #     schema=self.source_schema,
            #     rooting=self.rooting_interpretation,
            #     tree_offset=self.tree_offset,
            #     preserve_underscores=self.preserve_underscores,
            #     store_tree_weights=self.use_tree_weights,
            #     ignore_unrecognized_keyword_arguments=True,
            #     )
            try:
                _read_into_tree_array(
                        tree_array=self.tree_array,
                        tree_sources=[tree_source],
                        schema=self.source_schema,
                        taxon_namespace=self.taxon_namespace,
                        rooting=self.rooting_interpretation,
                        tree_offset=self.tree_offset,
                        use_tree_weights=self.use_tree_weights,
                        preserve_underscores=self.preserve_underscores,
                        info_message_func=self.send_info,
                        error_message_func=self.send_error,
                        log_frequency=self.log_frequency,
                        debug_mode=self.debug_mode,
                        )
            except (KeyboardInterrupt, Exception) as e:
                e.worker_name = self.name
                self.results_queue.put(e)
                break
            if self.kill_received:
                break
            self.num_tasks_completed += 1
            # self.send_info("Completed task {task_count}: '{task_name}'".format(
            self.send_info("Completed task: '{task_name}'".format(
                task_count=self.num_tasks_received,
                task_name=tree_source), wrap=False)
        if self.kill_received:
            self.send_warning("Terminating in response to kill request")
        else:
            self.results_queue.put(self.tree_array)

class TreeProcessor(object):

    def __init__(self,
            is_source_trees_rooted,
            ignore_edge_lengths,
            ignore_node_ages,
            use_tree_weights,
            ultrametricity_precision,
            taxon_label_age_map,
            num_processes,
            log_frequency,
            messenger,
            debug_mode,
            ):
        self.is_source_trees_rooted = is_source_trees_rooted
        self.rooting_interpretation = dendropy.get_rooting_argument(is_rooted=self.is_source_trees_rooted)
        self.ignore_edge_lengths = ignore_edge_lengths
        self.ignore_node_ages = ignore_node_ages
        self.use_tree_weights = use_tree_weights
        self.ultrametricity_precision = ultrametricity_precision
        self.taxon_label_age_map = taxon_label_age_map
        self.num_processes = num_processes
        self.log_frequency = log_frequency
        self.messenger = messenger
        self.debug_mode = debug_mode

    def info_message(self, msg, wrap=True, prefix=""):
        if self.messenger:
            self.messenger.info(msg, wrap=wrap, prefix=prefix)

    def warning_message(self, msg, wrap=True, prefix=""):
        if self.messenger:
            self.messenger.warning(msg, wrap=wrap, prefix=prefix)

    def error_message(self, msg, wrap=True, prefix=""):
        if self.messenger:
            self.messenger.error(msg, wrap=wrap, prefix=prefix)

    def analyze_trees(self,
            tree_sources,
            schema,
            taxon_namespace=None,
            tree_offset=0,
            preserve_underscores=False,
            ):
        if self.num_processes is None or self.num_processes <= 1:
            tree_array = self.serial_analyze_trees(
                    tree_sources=tree_sources,
                    schema=schema,
                    taxon_namespace=taxon_namespace,
                    tree_offset=tree_offset,
                    preserve_underscores=preserve_underscores,
                    )
        else:
            tree_array = self.parallel_analyze_trees(
                    tree_sources=tree_sources,
                    schema=schema,
                    taxon_namespace=taxon_namespace,
                    tree_offset=tree_offset,
                    preserve_underscores=preserve_underscores,
                    )
        return tree_array

    def serial_analyze_trees(self,
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
                taxon_label_age_map=self.taxon_label_age_map,
                )
        _read_into_tree_array(
                tree_array=tree_array,
                tree_sources=tree_sources,
                schema=schema,
                taxon_namespace=taxon_namespace,
                rooting=self.rooting_interpretation,
                tree_offset=tree_offset,
                use_tree_weights=self.use_tree_weights,
                preserve_underscores=preserve_underscores,
                info_message_func=self.info_message,
                error_message_func=self.error_message,
                log_frequency=self.log_frequency,
                debug_mode=self.debug_mode,
                )
        return tree_array

    def parallel_analyze_trees(self,
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
        self.info_message("{} taxa defined: {}".format( len(taxon_labels), taxon_labels))
        # max_idx_width = int(math.floor(math.log(len(taxon_labels), 10))) + 1
        # idx_col_width = (2 * max_idx_width) + 6
        # for tidx, taxon_label in enumerate(taxon_labels):
        #     index_col = "{idx:>{idx_col_width}}".format(
        #         idx=" ({}/{}): ".format(tidx+1, len(taxon_labels)),
        #         idx_col_width=idx_col_width,
        #         )
        #     self.info_message(taxon_label, prefix=index_col)


        # load up queue
        self.info_message("Creating work queue")
        work_queue = multiprocessing.Queue()
        for f in tree_sources:
            work_queue.put(f)

        # launch processes
        self.info_message("Launching {} worker processes".format(self.num_processes))
        results_queue = multiprocessing.Queue()
        messenger_lock = multiprocessing.Lock()
        workers = []
        for idx in range(self.num_processes):
            # self.info_message("Launching {} of {} worker processes".format(idx+1, self.num_processes))
            tree_analysis_worker = TreeAnalysisWorker(
                    name="Process-{}".format(idx+1),
                    work_queue=work_queue,
                    results_queue=results_queue,
                    source_schema=schema,
                    taxon_labels=taxon_labels,
                    tree_offset=tree_offset,
                    is_source_trees_rooted=self.is_source_trees_rooted,
                    preserve_underscores=preserve_underscores,
                    ignore_edge_lengths=self.ignore_edge_lengths,
                    ignore_node_ages=self.ignore_node_ages,
                    use_tree_weights=self.use_tree_weights,
                    ultrametricity_precision=self.ultrametricity_precision,
                    taxon_label_age_map=self.taxon_label_age_map,
                    messenger=self.messenger,
                    messenger_lock=messenger_lock,
                    log_frequency=self.log_frequency,
                    debug_mode=self.debug_mode)
            tree_analysis_worker.start()
            workers.append(tree_analysis_worker)

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
        try:
            while result_count < self.num_processes:
                result = results_queue.get()
                if isinstance(result, Exception) or isinstance(result, KeyboardInterrupt):
                    self.info_message("Exception raised in worker process '{}'".format(result.worker_name))
                    raise result
                master_tree_array.update(result)
                self.info_message("Recovered results from worker process '{}'".format(result.worker_name))
                result_count += 1
                # self.info_message("Recovered results from {} of {} worker processes".format(result_count, self.num_processes))
        except (Exception, KeyboardInterrupt) as e:
            for worker in workers:
                worker.terminate()
            raise
        self.info_message("All {} worker processes terminated".format(self.num_processes))
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
## Output

def _write_trees(trees,
        output_dest,
        args,
        file_comments):
    if args.output_tree_format == "newick" or args.output_tree_format == "phylip":
        trees.write_to_stream(
                output_dest,
                "newick",
                suppress_rooting=False,
                suppress_edge_lengths=True if args.edge_length_summarization == "clear" else False,
                unquoted_underscores=True if args.preserve_underscores else False,
                preserve_spaces=True if args.preserve_underscores else False,
                store_tree_weights=True,
                suppress_annotations=args.suppress_annotations,
                suppress_item_comments=args.clear_item_comments,
                )
    elif args.output_tree_format == "nexus":
        trees.write_to_stream(
                output_dest,
                "nexus",
                suppress_rooting=False,
                suppress_edge_lengths=True if args.edge_length_summarization == "clear" else False,
                unquoted_underscores=True if args.preserve_underscores else False,
                preserve_spaces=True if args.preserve_underscores else False,
                store_tree_weights=True,
                suppress_annotations=args.suppress_annotations,
                suppress_item_comments=args.clear_item_comments,
                simple=args.no_taxa_block,
                file_comments=file_comments,
                )
    elif args.output_tree_format == "nexml":
        if file_comments:
            trees.comments = file_comments
        trees.write_to_stream(
                output_dest,
                "nexml",
                )
    else:
        raise ValueError(args.output_tree_format)

##############################################################################
## Front-End

def print_citation(args):
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
            --summary-target consensus \\
            --min-clade-freq=0.95 \\
            --edges mean-length \\
            --burnin=200 \\
            --support-as-labels \\
            --output=result.tre \\
            treefile1.tre treefile2.tre treefile3.tre

To use a different type of summary tree, e.g., the tree that maximizes the
product of posterior probabilities, you can specify 'mcct' for the
'--summary-tree' option:

    $ sumtrees.py \\
            --summary-target mcct \\
            --min-clade-freq=0.95 \\
            --edges mean-length \\
            --burnin=200 \\
            --support-as-labels \\
            --output=result.tre \\
            treefile1.tre treefile2.tre treefile3.tre

If the input trees are ultrametric and you want to set the node ages to the
median node age, set the '--edges' argument to 'median-age':

    $ sumtrees.py \\
            --summary-target mcct \\
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
    return dendropy.description(dest=dest)

def main():

    ######################################################################
    ## Start Recording Total Job Time

    main_time_start = datetime.datetime.now()

    ######################################################################
    ## CLI

    parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=cli.CustomFormatter,
            add_help=False,
            )
    source_options = parser.add_argument_group("Source Options")
    source_options.add_argument("tree_sources",
            nargs="*",
            metavar="TREE-FILEPATH",
            help= (
                "Source(s) of trees to summarize. At least one valid"
                " source of trees must be provided. Use '-' to specify"
                " reading from standard input (note that this requires"
                " the input file format to be explicitly set using"
                " the '--source-format' option)."
            ))
    source_options.add_argument("-i", "--input-format", "--source-format",
            metavar="FORMAT",
            default=None,
            choices=["nexus/newick", "nexus", "newick", "phylip", "nexml"],
            help=(
                 "Format of all input trees (defaults to handling either NEXUS"
                 " or NEWICK through inspection; it is more efficient to"
                 " explicitly specify the format if it is known)."
                 ))
    source_options.add_argument("-b", "--burnin",
            type=int,
            default=0,
            help=(
                 "Number of trees to skip from the beginning of *each* tree "
                 "file when counting support (default: %(default)s)."
                 ))
    source_options.add_argument("--force-rooted", "--rooted",
            dest="is_source_trees_rooted",
            action="store_true",
            default=None,
            help="Treat source trees as rooted.")
    source_options.add_argument("--force-unrooted", "--unrooted",
            dest="is_source_trees_rooted",
            action="store_false",
            default=None,
            help="Treat source trees as unrooted.")
    source_options.add_argument("--weighted-trees",
            action="store_true",
            default=False,
            help=(
                "Use weights of trees (as indicated by '[&W m/n]' comment token) "
                "to weight contribution of splits found on each tree to overall "
                "split frequencies."
                ))
    source_options.add_argument("--preserve-underscores",
            action="store_true",
            default=False,
            help=(
                "Do not convert unprotected (unquoted) underscores to spaces"
                " when reading NEXUS/NEWICK format trees."
                ))
    source_options.add_argument("-v", "--ultrametricity-precision", "--edge-weight-epsilon", "--branch-length-epsilon",
            type=float,
            default=constants.DEFAULT_ULTRAMETRICITY_PRECISION,
            help="Precision to use when validating ultrametricity (default: %(default)s; specify '0' to disable validation).")
    source_options.add_argument(
            "--taxon-name-filepath", "--taxon-names-filepath",
            metavar="FILEPATH",
            default=None,
            help=(
                "Path to file listing all the taxon names or labels that"
                " will be found across the entire set of source trees."
                " This file should be a plain text file with a single"
                " name list on each line. This file is only read when"
                " multiprocessing ('-M' or '-m') is requested. When"
                " multiprocessing using the '-M' or '-m' options,"
                " all taxon names need to be defined in advance"
                " of any actual tree analysis. By default this is done"
                " by reading the first tree in the first tree source"
                " and extracting the taxon names. At best, this is,"
                " inefficient, as it involves an extraneous reading of"
                " the tree. At worst, this can be errorneous, if the"
                " first tree does not contain all the taxa. Explicitly"
                " providing the taxon names via this option can avoid"
                " these issues."
                ))
    source_options.add_argument("--tip-ages", "--tip-ages-filepath",
            dest="tip_ages_filepath",
            metavar="FILEPATH",
            default=None,
            help=(
                "Path to file providing ages (i.e., time from present) of"
                " tips. For format of this file, see '--tip-age-format'."
                " If not specified, or for any taxon omitted from the data,"
                " an age of 0.0 will be assumed."
                ))
    source_options.add_argument("--tip-ages-format",
            dest="tip_ages_format",
            default="tsv",
            choices=["tsv", "csv", "json",],
            metavar="FORMAT",
            help=cli.CustomFormatter.format_definition_list_help(
                    preamble=
                        (
                        "Format of the tip date data (default: '%(default)s'):"
                        ),
                    definitions=
                        (
                            ("'tsv'",
                                "A tab-delimited file. This should consist of two columns"
                                " separated by tabs. The first column lists the taxon labels"
                                " (matching the taxon label of the input trees EXACTLY)"
                                " and the second column lists the ages of the"
                                " corresponding tips."
                            ),
                            ("'csv'",
                                "A comma-delimited file. This should consist of two columns"
                                " separated by commas. The first column lists the taxon labels"
                                " (matching the taxon label of the input trees EXACTLY)"
                                " and the second column lists the ages of the"
                                " corresponding tips."
                            ),
                            ("'json'",
                                "A JSON file. This should specify a single"
                                " dictionary at the top-level with keys being taxon"
                                " labels (matching the taxon labels of the input "
                                " trees EXACTLY) and values being the ages of the "
                                " corresponding tips."
                            ),
                        )
                ))
    source_options.add_argument("--no-trim-tip-age-labels",
            action="store_false",
            default=True,
            help="By default, whitespace will be trimmed from the labels"
                 " found in the tip ages data source. Specifing this option"
                 " suppresses this.")

    target_tree_options = parser.add_argument_group("Target Tree Topology Options")
    target_tree_options.add_argument(
            "-t", "--target-tree-filepath",
            default=None,
            metavar="FILE",
            help=(
                    "Summarize support and other information from the source"
                    " trees to topology or topologies given by the tree(s)"
                    " described in FILE. If no use-specified target topologies"
                    " are given, then a summary topology will be used as the"
                    " target. Use the '-s' or '--summary-target' to specify the"
                    " type of summary tree to use."
                 ))
    target_tree_options.add_argument(
            "-s", "--summary-target",
            default=None,
            choices=["consensus", "mcct", "msct"],
            metavar="SUMMARY-TYPE",
            help=cli.CustomFormatter.format_definition_list_help(
                    preamble=
                        (
                        "Construct and summarize support and other information "
                        "from the source trees to one of the following summary "
                        "topologies: "
                        ),
                    definitions=
                        (
                            ("'consensus'",
                                "A consensus tree. The minimum frequency       "
                                "threshold of clades to be included            "
                                "can be specified using the '-f' or            "
                                "'--min-clade-freq' flags. This is the DEFAULT "
                                "if a user- specified target tree is not given "
                                "through the '-t' or '--target-tree-filepath'       "
                                "options.                                      "
                            ),
                            ("'mcct'",
                                "The maximum clade credibility tree.   "
                                "The tree from the source set that     "
                                "maximizes the *product* of clade      "
                                "posterior probabilities.              "
                            ),
                            ("'msct'",
                                "The maximum sum of clade credibilities tree.   "
                                "The tree from the source set that     "
                                "maximizes the *sum* of clade      "
                                "posterior probabilities.              "
                            ),
                        )
                ))
    target_tree_supplemental_options = parser.add_argument_group("Target Tree Supplemental Options")
    target_tree_supplemental_options.add_argument("-f", "--min-clade-freq", "--min-freq", "--min-split-freq", "--min-consensus-freq",
            dest="min_clade_freq",
            type=float,
            # default=constants.GREATER_THAN_HALF,
            default=None,
            metavar="#.##",
            help=(
                "Minimum frequency or "
                "probability for a clade or a split to be included in "
                "the summary target trees. "
                "If user-defined or non-consensus trees are specified "
                "as summary targets "
                "and a explicit value is provided for this argument, "
                "then clades with support values below this threshold "
                "will be collapsed. "
                "If a consensus tree summary target is specified, "
                "then clades with support values below this threshold "
                "will not be included, and this threshold takes on a "
                "default value of greater than 0.5 if not explicitly "
                "specified."))
    target_tree_supplemental_options.add_argument("--allow-unknown-target-tree-taxa",
            action="store_true",
            default=False,
            help=(
                "Do not fail with error if target tree(s) have taxa not"
                " previously encountered in source trees or defined in"
                " the taxon discovery file."
                ))
    # target_tree_supplemental_options.add_argument(
    #         "--collapse-edges-with-less-than-minimum-support",
    #         dest="collapse_edges_with_less_than_minimum_support",
    #         action="store_true",
    #         default=False,
    #         help=(
    #             "Collapse edges that have support values less than that specifed "
    #             "by '-f', '--min-freq'. By default, the '--min-freq' option "
    #             "only applies for consensus tree summary targets. Specifying "
    #             "this option will apply it to all other types of summary targets "
    #             "including user-specified target trees."
    #             ))

    target_tree_rooting_options = parser.add_argument_group("Target Tree Rooting Options")
    target_tree_rooting_options.add_argument("--root-target-at-outgroup",
            dest="root_target_at_outgroup",
            metavar="TAXON-LABEL",
            default=None,
            help="Root target tree(s) using specified taxon as outgroup.")
    target_tree_rooting_options.add_argument("--root-target-at-midpoint",
            action="store_true",
            default=None,
            help="Root target tree(s) at midpoint.")
    target_tree_rooting_options.add_argument("--set-outgroup",
            dest="set_outgroup",
            metavar="TAXON-LABEL",
            default=None,
            help="Rotate the target trees such the specified taxon is in the outgroup position, but do not explicitly change the target tree rooting.")

    edge_length_summarization_options = parser.add_argument_group("Target Tree Edge Options")
    edge_length_summarization_choices = ["mean-length", "median-length", "mean-age", "median-age", "support", "keep", "clear",]
    edge_length_summarization_options.add_argument(
            "-e",
            "--set-edges",
            "--edges",
            dest="edge_length_summarization",
            metavar="STRATEGY",
            choices=edge_length_summarization_choices,
            default=None,
            help=cli.CustomFormatter.format_definition_list_help(
                    preamble=
                        (
                        "Set the edge lengths of the target or        "
                        "summary trees based on the specified         "
                        "summarization STRATEGY:                      "
                        ),
                    definitions=
                        (
                            ("'mean-length'",
                                "Edge lengths will be set to the mean  "
                                "of the lengths of the corresponding   "
                                "split or clade in the source trees.   "
                            ),
                            ("'median-length'",
                                " Edge lengths will be set to the      "
                                " median of the lengths of the         "
                                " corresponding split or clade in the  "
                                " source trees.                        "
                            ),
                            ("'mean-age'",
                                "Edge lengths will be adjusted so      "
                                "that the age of subtended nodes will  "
                                "be equal to the mean age of the       "
                                "corresponding split or clade in the   "
                                "source trees. Source trees will need  "
                                "to to be ultrametric for this option. "
                            ),
                            ("'median-age'",
                                "Edge lengths will be adjusted so     "
                                "that the age of subtended nodes will "
                                "be equal to the median age of the    "
                                "corresponding split or clade in the  "
                                "source trees. Source trees will need "
                                "to to be ultrametric for this option."
                                            ),
                            ("support",
                                "Edge lengths will be set to the       "
                                "support value for the split           "
                                "represented by the edge.              "
                            ),
                            ("'keep'",
                                "Do not change the existing edge       "
                                "lengths. This is the DEFAULT if       "
                                "target tree(s) are sourced from an    "
                                "external file using the '-t' or       "
                                "'--target-tree-filepath' option            "
                            ),
                            ("'clear'",
                                "Edge lengths will be cleared from the "
                                "target trees if they are present.     "
                            ),
                        ),
                    coda="\n".join((
                            "<pre>Note the default settings varies according to the ",
                            "following, in order of preference:                  ",
                            "(1) If target trees are specified using the '-t' or ",
                            "    '--target-tree-filepath' option, then the default edge ",
                            "    summarization strategy is: 'keep'. ",
                            "(2) If target trees are not specified, but the ",
                            "    '--summarize-node-ages' option is specified, ",
                            "    then the default edge summarization strategy is: ",
                            "    'mean-age'. ",
                            "(3) If no target trees are specified and the ",
                            "    node ages are NOT specified to be summarized, ",
                            "    then the default edge summarization strategy is: ",
                            "    'mean-length'. ",
                        ))
                ))
    edge_length_summarization_options.add_argument("--force-minimum-edge-length",
            default=None,
            type=float,
            help="(If setting edge lengths) force all edges to be at least this length.")
    edge_length_summarization_options.add_argument("--collapse-negative-edges",
            action="store_true",
            default=False,
            help="(If setting edge lengths) force parent node ages to be at least as old as its oldest child when summarizing node ages.")

    node_summarization_options = parser.add_argument_group("Target Tree Annotation Options")
    node_summarization_options.add_argument(
            "--summarize-node-ages", "--ultrametric", "--node-ages",
            action="store_true",
            dest="summarize_node_ages",
            default=None,
            help="Assume that source trees are ultrametic and summarize node ages (distances from tips).")
    node_summarization_options.add_argument("-l","--labels",
            dest="node_labels",
            default="support",
            choices=["support", "keep", "clear",],
            help=cli.CustomFormatter.format_definition_list_help(
                preamble="Set the node labels of the summary or target tree(s):",
                definitions=(
                    ("'support'",
                        "Node labels will be set to the support value  "
                        "for the clade represented by the node. This is "
                        "the DEFAULT.                                  "
                    ),
                    ("'keep'",
                        "Do not change the existing node labels."
                    ),
                    ("'clear'",
                        "Node labels will be cleared from the target   "
                        "trees if they are present.                    "
                    )
                )
                ))
    node_summarization_options.add_argument("--suppress-annotations", "--no-annotations",
            action="store_true",
            default=False,
            help=(
                "Do NOT annotate nodes and edges with any summarization information metadata such as."
                "support values, edge length and/or node age summary statistcs, etc."
                ))

    support_expression_options = parser.add_argument_group("Support Expression Options")
    support_expression_options.add_argument("-p", "--percentages",
            action="store_true",
            dest="support_as_percentages",
            default=False,
            help="Indicate branch support as percentages (otherwise, will report as proportions by default).")
    support_expression_options.add_argument("-d", "--decimals",
            dest="support_label_decimals",
            type=int,
            metavar="#",
            default=8,
            help="Number of decimal places in indication of support values (default: %(default)s).")
    # other_summarization_options.add_argument("--no-summarize-edge-lengths",
    #         action="store_false",
    #         dest="summarize_edge_lengths",
    #         default=None,
    #         help="Do not summarize edge lengths.")

    output_options = parser.add_argument_group("Output Options")
    output_options.add_argument("-o","--output-tree-filepath", "--output",
            metavar="FILEPATH",
            default=None,
            help="Path to output file (if not specified, will print to standard output).")
    output_options.add_argument("-F","--output-tree-format",
            default=None,
            choices=["nexus", "newick", "phylip", "nexml"],
            help="Format of the output tree file (if not specifed, defaults to input format, if this has been explicitly specified, or 'nexus' otherwise).")
    output_options.add_argument("-x", "--extended-output",
            dest="extended_output_prefix",
            default=None,
            metavar="PREFIX",
            help=cli.CustomFormatter.format_definition_list_help(
                    preamble=
                        (
                        "If specified, extended summarization information "
                        "will be generated, consisting of the following "
                        "files:"
                        ),
                    definitions=
                        (
                            # ("'<PREFIX>.summary.trees'",
                            #     "Summary or target trees onto which the "
                            #     "summarization information from the source set "
                            #     "has been mapped."
                            # ),
                            ("'<PREFIX>.topologies.trees'",
                                "A collection of topologies found in the sources "
                                "reported with their associated posterior "
                                "probabilities as metadata annotations."
                            ),
                            ("'<PREFIX>.bipartitions.trees'",
                                "A collection of bipartitions, each represented as "
                                "a tree, with associated information as metadata"
                                "annotations."
                            ),
                            ("'<PREFIX>.bipartitions.tsv'",
                                "Table listing bipartitions as a group pattern as "
                                "the key column, and information regarding each "
                                "the bipartitions as the remaining columns."
                            ),
                            ("'<PREFIX>.edge-lengths.tsv'",
                                "List of bipartitions and "
                                "corresponding edge lengths. Only "
                                "generated if edge lengths are "
                                "summarized. "
                            ),
                            ("'<PREFIX>.node-ages.tsv'",
                                "List of bipartitions and corresponding "
                                "ages. Only generated if node ages are "
                                "summarized. "
                            ),
                        )
                        ))
    output_options.add_argument("--no-taxa-block",
            action="store_true",
            default=False,
            help="When writing NEXUS format output, do not include a taxa block in the output treefile (otherwise will create taxa block by default).")
    output_options.add_argument("--no-analysis-metainformation", "--no-meta-comments",
            dest="suppress_analysis_metainformation",
            action="store_true",
            default=False,
            help="Do not include meta-information describing the summarization parameters and execution details.")
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

    deprecated_output_options = parser.add_argument_group("Deprecated Output Options")
    deprecated_output_options.add_argument("--trprobs", "--calc-tree-probabilities",
            dest="trprobs_filepath",
            default=None,
            metavar="FILEPATH",
            help=argparse.SUPPRESS,
            )
    deprecated_output_options.add_argument("--support-as-labels",
            action="store_true",
            default=None,
            help=argparse.SUPPRESS,
            )
    deprecated_output_options.add_argument("--extract-edges",
            dest="split_edge_map_filepath",
            default=None,
            metavar="FILEPATH",
            help=argparse.SUPPRESS,
            )

    multiprocessing_options = parser.add_argument_group("Parallel Processing Options")
    multiprocessing_options.add_argument("-M", "--maximum-multiprocessing",
            action="store_const",
            const="max",
            dest="multiprocess",
            help=(
                 "Run in parallel mode using as many processors as available, up to the number of sources."
                 ))
    multiprocessing_options.add_argument("-m", "--multiprocessing",
            dest="multiprocess",
            metavar="NUM-PROCESSES",
            help=(
                 "Run in parallel mode with up to a maximum of NUM-PROCESSES processes "
                 "('max' or '#' means to run in as many processes as there are cores on the "
                 "local machine; i.e., same as specifying '-M' or '--maximum-multiprocessing')."
                 ))

    logging_options = parser.add_argument_group("Program Logging Options")
    logging_options.add_argument("-g", "--log-frequency",
            type=int,
            metavar="LOG-FREQUENCY",
            default=500,
            help="Tree processing progress logging frequency (default: %(default)s; set to 0 to suppress).")
    logging_options.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Suppress ALL logging, progress and feedback messages.")

    error_options = parser.add_argument_group("Program Error Options")
    error_options.add_argument("--ignore-missing-support",
            action="store_true",
            default=False,
            help="Ignore missing support tree files (note that at least one must exist).")

    information_options = parser.add_argument_group("Program Information Options")
    information_options.add_argument("-h", "--help",
            action="store_true",
            default=False,
            help="Show help information for program and exit.")
    information_options.add_argument("--version","--citation",
            dest="citation",
            action="store_true",
            default=False,
            help="Show version and citation information for program and exit.")
    information_options.add_argument("--usage-examples",
            action="store_true",
            default=False,
            help="Show usage examples of program and exit.")
    information_options.add_argument("--describe",
            action="store_true",
            default=False,
            help="Show information regarding your DendroPy and Python installations and exit.")

    parser.add_argument('--debug-mode', action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args()

    ######################################################################
    ## add stuff here: incorporate into CLI later
    args.clear_item_comments = False

    ######################################################################
    ## Information (Only) Operations

    if args.citation:
        print_citation(args)
        sys.exit(0)

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

    processing_report_lines = []
    def _message_and_log(msg, wrap=True, prefix=""):
        messenger.info(msg, wrap=wrap, prefix=prefix)
        processing_report_lines.append(msg)
    def _bulleted_message_and_log(msg, prefix="- "):
        messenger.info(msg, wrap=True, prefix=prefix)
        processing_report_lines.append(prefix + msg)

    ######################################################################
    ## Set up some common messages

    mixed_rooting_solution = (
                            "Re-run SumTrees using the"
                            " '--force-rooted' or the"
                            " '--force-unrooted' option to force a"
                            " consistent rooting state for all"
                            " trees."
                            )

    ######################################################################
    ## Support Files

    if len(args.tree_sources) == 0:
        parser.print_usage()
        sys.stdout.write("\n")
        sys.stdout.write("Type 'sumtrees.py --help' for details on usage.\n")
        sys.stdout.write("Type 'sumtrees.py --usage-examples' for examples of usage.\n")
        sys.exit(0)

    tree_sources = []
    ignored_sources = []
    for fpath in args.tree_sources:
        if fpath == "-":
            if args.input_format is None:
                messenger.error("Format of source trees must be specified using '--source-format' flag when reading trees from standard input")
                sys.exit(1)
            elif args.input_format.lower() == "nexus/newick":
                messenger.error("The 'nexus/newick' format is not supported when reading trees from standard input")
                sys.exit(1)
            if len(args.tree_sources) > 1:
                messenger.error("Cannot specify multiple sources when reading from standard input")
            tree_sources = None
            break
        else:
            if args.input_format is None:
                args.input_format = "nexus/newick"
            else:
                args.input_format = args.input_format.lower()
            fpath = os.path.expanduser(os.path.expandvars(fpath))
            if not os.path.exists(fpath):
                missing_msg = "Ignoring missing source file: '{}'".format(fpath)
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
    if tree_sources is None:
        tree_sources = [sys.stdin]
        messenger.info("Reading trees from standard input")
        processing_report_lines.append("Trees read from standard input source")
    elif len(tree_sources) == 0:
            messenger.error("No valid sources of input trees specified. "
                    + "Please provide the path to at least one (valid and existing) file "
                    + "containing tree samples to summarize.")
            sys.exit(1)
    else:
        # messenger.info("{} source(s) to be analyzed and summarized: {}".format(
        #     len(tree_sources),
        #     tree_sources))
        messenger.info("Trees to be read from {} source(s):".format(len(tree_sources)))
        processing_report_lines.append("Trees read from {} source(s):".format(len(tree_sources)))
        max_idx_width = int(math.floor(math.log(len(tree_sources), 10))) + 1
        idx_col_width = (2 * max_idx_width) + 6
        for tidx, tree_source in enumerate(tree_sources):
            # index_col = "{idx:>{idx_col_width}}".format(
            #     idx=" ({}/{}): ".format(tidx+1, len(tree_sources)),
            #     idx_col_width=idx_col_width,
            #     )
            # if tree_source is sys.stdin:
            #     _bulleted_message_and_log("<standard input>", prefix=index_col)
            # else:
            #     _bulleted_message_and_log("'" + tree_source + "'", prefix=index_col)
            if tree_source is sys.stdin:
                _bulleted_message_and_log("<standard input>")
            else:
                _bulleted_message_and_log("'" + tree_source + "'")
    if args.burnin:
        messenger.info("{} initial trees to be discarded/ignored as burn-in from *each* source".format(args.burnin))
        processing_report_lines.append("{} initial trees discarded/ignored as burn-in from *each* source".format(args.burnin))

    ######################################################################
    ## Target Validation

    if args.summary_target is not None and args.target_tree_filepath is not None:
        messenger.error("Cannot specify both '-s'/'--summary-tree-target' and '-t'/'--target-tree-filepath' simultaneously")
    elif args.target_tree_filepath is not None:
        target_tree_filepath = os.path.expanduser(os.path.expandvars(args.target_tree_filepath))
        if not os.path.exists(target_tree_filepath):
            messenger.error("Target tree file not found: '{}'".format(target_tree_filepath))
            sys.exit(1)
    else:
        target_tree_filepath = None
        if args.summary_target is None:
            args.summary_target = "consensus"

    ######################################################################
    ## Tree Rooting

    if args.root_target_at_outgroup is not None or args.root_target_at_midpoint:
        if not args.is_source_trees_rooted:
            messenger.info("Rooting directive specified for target tree(s): source trees will also be treated as rooted")
            args.is_source_trees_rooted = True

    num_target_rooting_directives = 0
    if args.root_target_at_outgroup is not None:
        num_target_rooting_directives += 1
    if args.set_outgroup is not None:
        num_target_rooting_directives += 1
    if args.root_target_at_midpoint:
        num_target_rooting_directives += 1
    if num_target_rooting_directives > 1:
        messenger.error("Only one target tree rooting directive can be specified")
        sys.exit(1)

    ######################################################################
    ## Node Age Summarization

    if args.edge_length_summarization in ("mean-age", "median-age"):
        args.summarize_node_ages = True

    ######################################################################
    ## Tip Ages

    if args.tip_ages_filepath is None:
        taxon_label_age_map = None
    else:
        if not args.summarize_node_ages:
            messenger.error("Tip age data file specified, but node ages are not going be analyzed."
                            " Specify '--summarize-node-ages' or '--edges=mean-age' or --edges=median-age'"
                            " to analyze node ages.")
            sys.exit(1)
        tip_ages_filepath = os.path.expanduser(os.path.expandvars(args.tip_ages_filepath))
        messenger.info("Tip ages will be read from: '{}'".format(tip_ages_filepath))
        with open(tip_ages_filepath, "r") as src:
            if args.tip_ages_format == "csv" or args.tip_ages_format == "tsv":
                raw_tip_data_map = collections.OrderedDict()
                if args.tip_ages_format == "csv":
                    delimiter=","
                    delimiter_desc = "',' (comma)"
                else:
                    delimiter="\t"
                    delimiter_desc = "'\\t' (the tab character)"
                csv_reader = csv.reader(src, delimiter=delimiter)
                for row_idx, row in enumerate(csv_reader):
                    try:
                        taxon_label, age = row
                    except ValueError as e:
                        messenger.error("Tip age data file '{}', line {}: error reading columns. Perhaps wrong delimiter was specified? Expecting {} as a delimiter. Use '--tip-ages-format' to change this.".format(
                            tip_ages_filepath, row_idx, delimiter_desc))
                        sys.exit(1)
                    try:
                        raw_tip_data_map[taxon_label] = float(age)
                    except ValueError as e:
                        messenger.error("Tip age data file '{}', line {} (label = '{}'): invalid age value: {}".format(
                            tip_ages_filepath, row_idx, taxon_label, age))
                        sys.exit(1)
            elif args.tip_ages_format == "json":
                raw_tip_data_map = collections.OrderedDict(json.load(src))
            else:
                messenger.error("Unrecognized or unsupported format for tip age data: '{}'".format(args.tip_ages_format))
                sys.exit(1)
        taxon_label_age_map = {}
        for taxon_label in raw_tip_data_map:
            try:
                age = float(raw_tip_data_map[taxon_label])
            except ValueError:
                messenger.error("Tip age data file '{}': taxon with label '{}': invalid value for age: {}".format(
                    tip_ages_filepath, taxon_label, raw_tip_data_map[taxon_label]))
                sys.exit(1)
            if not args.no_trim_tip_age_labels:
                taxon_label = taxon_label.strip()
            if args.input_format in ("nexus/newick", "nexus", "newick"):
                if not args.preserve_underscores:
                    if not (
                            (taxon_label[0] == taxon_label[-1] == "'")
                            or (taxon_label[0] == taxon_label[-1] == '"')
                            ):
                        taxon_label = taxon_label.replace("_", " ")
            taxon_label_age_map[taxon_label] = age

    ######################################################################
    ## Ultrametricity Precision

    if args.ultrametricity_precision == 0:
        # underlying function expects negative value to disable, but SumTrees
        # API uses 0
        args.ultrametricity_precision = -1

    ######################################################################
    ## Output File Setup

    # legacy
    if args.trprobs_filepath:
        messenger.error(
                "The '--trprobs' or '--calc-tree-probabilities' "
                "option is no longer supported directly. Use '-x' or "
                "'--extended-output' to specify an extended suite of "
                "output, which includes the topology probabilities. "
                )
        sys.exit(1)
    if args.split_edge_map_filepath:
        messenger.error(
                "The '--extract-edges' option is no longer supported. "
                "Use '-x' or '--extended-output' to specify an "
                "extended suite of output which includes this "
                "information. "
                )
        sys.exit(1)

    # output format
    if args.output_tree_format is None:
        if args.input_format is None or args.input_format == "nexus/newick":
            args.output_tree_format = "nexus"
        else:
            args.output_tree_format = args.input_format

    # primary output
    if args.output_tree_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(args.output_tree_filepath))
        if cli.confirm_overwrite(filepath=output_fpath, replace_without_asking=args.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)

    # extended output
    extended_output_paths = {}
    if args.extended_output_prefix is not None:
        if not args.extended_output_prefix.endswith("."):
            args.extended_output_prefix += "."
        for results_key, suffix in (
                    ("summary-trees", "summary.trees"),
                    ("topologies", "topologies.trees"),
                    ("bipartition-trees", "bipartitions.trees"),
                    ("bipartition-table", "bipartitions.tsv"),
                    ("edge-lengths", "edge-lengths.tsv"),
                    ("node-ages", "node-ages.tsv"),
                ):
            full_path = args.extended_output_prefix + suffix
            # if full_path.endswith("trees") and args.output_tree_format == "nexml":
            #     full_path += ".nexml"
            if cli.confirm_overwrite(
                    filepath=full_path,
                    replace_without_asking=args.replace):
                extended_output_paths[results_key] = full_path
            else:
                sys.exit(1)

    ######################################################################
    ## Multiprocessing Setup

    num_cpus = multiprocessing.cpu_count()
    if len(tree_sources) > 1 and args.multiprocess is not None:
        if (
                args.multiprocess.lower() == "max"
                or args.multiprocess == "#"
                or args.multiprocess == "*"
            ):
            num_processes = min(num_cpus, len(tree_sources))
        # elif args.multiprocess == "@":
        #     num_processes = len(tree_sources)
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
        if len(tree_sources) > 1 and num_cpus > 1:
            messenger.info(
                    ("Multiple processors ({num_cpus}) available:"
                    " consider using the '-M' or '-m' options to"
                    " parallelize processing of trees"
                    ).format(num_cpus=num_cpus))
        num_processes = 1

    ######################################################################
    ## Taxon Discovery

    if args.taxon_name_filepath is not None:
        with open(os.path.expanduser(os.path.expandvars(args.taxon_name_filepath)), "r") as tnf:
            taxon_labels = [name.strip() for name in tnf.read().split("\n") if name]
            taxon_labels = [name for name in taxon_labels if name]
        taxon_namespace = dendropy.TaxonNamespace(taxon_labels)
    else:
        taxon_namespace = None

    ######################################################################
    ## Main Work

    tree_processor = TreeProcessor(
            is_source_trees_rooted=args.is_source_trees_rooted,
            ignore_edge_lengths=False,
            ignore_node_ages=not args.summarize_node_ages,
            use_tree_weights=args.weighted_trees,
            ultrametricity_precision=args.ultrametricity_precision,
            taxon_label_age_map=taxon_label_age_map,
            num_processes=num_processes,
            log_frequency=args.log_frequency if not args.quiet else 0,
            messenger=messenger,
            debug_mode=args.debug_mode,
            )
    analysis_time_start = datetime.datetime.now()
    # messenger.info("Processing of source trees starting at {}".format(
    #     analysis_time_start,
    #     ))
    try:
        tree_array = tree_processor.analyze_trees(
                tree_sources=tree_sources,
                schema=args.input_format,
                taxon_namespace=taxon_namespace,
                tree_offset=args.burnin,
                preserve_underscores=args.preserve_underscores,
                )
        if tree_array.split_distribution.is_mixed_rootings_counted():
            raise TreeArray.IncompatibleRootingTreeArrayUpdate("Mixed rooting states detected in source trees")
    except KeyboardInterrupt as e:
        raise e
    except Exception as exception_object:
        if isinstance(exception_object, error.MixedRootingError) or isinstance(exception_object, dendropy.TreeArray.IncompatibleRootingTreeArrayUpdate):
            error_message_epilog = mixed_rooting_solution
        else:
            error_message_epilog = ""
        message = []
        if hasattr(exception_object, "exception_tree_source_name"):
            source_name = exception_object.exception_tree_source_name
            source_offset = exception_object.exception_tree_offset
            subparts = []
            if source_name is not None:
                subparts.append("'{}'".format(source_name))
                if source_offset is not None:
                    subparts.append(", tree at offset {}".format(source_offset))
            else:
                if source_offset is not None:
                    subparts.append("Tree at offset {}".format(source_offset))
            message.append("".join(subparts) + ":")
        message.append(str(exception_object))
        if error_message_epilog:
            if not message[-1].endswith("."):
                message[-1] = message[-1] + "."
            message.append(error_message_epilog)
        message = " ".join(message)
        messenger.error(message)
        if args.debug_mode:
            raise
        sys.exit(1)
    analysis_time_end = datetime.datetime.now()
    analysis_time_delta =  analysis_time_end-analysis_time_start

    ### Validate/match taxa with tip ages
    ### We do this here so as to avoid reporting a "job complete"
    ### if there are errors
    if taxon_label_age_map:
        taxon_age_map = {}
        for taxon_label, age in taxon_label_age_map.items():
            taxon = tree_array.taxon_namespace.get_taxon(taxon_label)
            if taxon is None:
                messenger.error("Unable to find taxon with label matching tip data label '{}' in the taxon namespace: {}".format(
                    taxon_label,
                    [t.label for t in tree_array.taxon_namespace]))
                sys.exit(1)
            taxon_age_map[taxon] = taxon_label_age_map[taxon_label]

    messenger.info("Analysis of source trees completed in: {}".format(timeprocessing.pretty_timedelta(analysis_time_delta),
        wrap=False,
        ))

    ######################################################################
    ## Post-Processing

    ### post-analysis reports

    if len(tree_array) == 0:
        messenger.error("No trees retained for processing (is the burn-in too high?)")
        sys.exit(1)

    _message_and_log("Total of {} trees analyzed for summarization:".format(len(tree_array)))
    if args.weighted_trees:
        _bulleted_message_and_log("All trees were treated as weighted (default weight = 1.0).")
    else:
        _bulleted_message_and_log("All trees were treated as unweighted")
    if args.is_source_trees_rooted is None:
        if tree_array.split_distribution.is_all_counted_trees_rooted():
            _bulleted_message_and_log("All trees were rooted")
        elif tree_array.split_distribution.is_all_counted_trees_strictly_unrooted():
            _bulleted_message_and_log("All trees were unrooted")
        elif tree_array.split_distribution.is_all_counted_trees_treated_as_unrooted():
            _bulleted_message_and_log("All trees were assumed to be unrooted")
    elif args.is_source_trees_rooted is True:
        _bulleted_message_and_log("All trees were treated as rooted")
    else:
        _bulleted_message_and_log("All trees were treated as unrooted")
    # if args.is_source_trees_ultrametric and args.ultrametricity_precision:
    #     _bulleted_message_and_log("Trees were ultrametric within an error of {}".format(args.ultrametricity_precision))
    # elif args.is_source_trees_ultrametric:
    #     _bulleted_message_and_log("Trees were expected to be ultrametric (not verified)")
    _bulleted_message_and_log("{} unique taxa across all trees".format(len(tree_array.taxon_namespace)))
    num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits = tree_array.split_distribution.splits_considered()
    if args.weighted_trees:
        _bulleted_message_and_log("{} unique splits with a total weight of {}".format(num_unique_splits, num_splits))
        _bulleted_message_and_log("{} unique non-trivial splits with a total weight of {}".format(num_nt_unique_splits, num_nt_splits))
    else:
        _bulleted_message_and_log("{} unique splits out of a total of {} splits".format(num_unique_splits, int(num_splits)))
        _bulleted_message_and_log("{} unique non-trivial splits counted out of a total of non-trivial {} splits".format(num_nt_unique_splits, int(num_nt_splits)))

    ###  tip ages
    if not taxon_label_age_map:
        pass
        # _message_and_log("All tips assumed to be contemporaneous with age of 0.0")
    else:
        _message_and_log("Tips assigned the following ages for node age analysis:")
        taxon_labels = [taxon.label for taxon in tree_array.taxon_namespace]
        max_taxon_label_length = max([len(x) for x in taxon_labels])
        taxon_age_template = "{{:>{}}} : {{}}".format(max_taxon_label_length)
        for taxon in tree_array.taxon_namespace:
            _bulleted_message_and_log(taxon_age_template.format(taxon.label, taxon_label_age_map.get(taxon.label, 0.0)))

    ### build target tree(s)
    target_trees = dendropy.TreeList(taxon_namespace=tree_array.taxon_namespace)
    if target_tree_filepath is None:
        args.include_external_splits_when_scoring_clade_credibility_tree = False
        if args.include_external_splits_when_scoring_clade_credibility_tree:
            coda = ", including tip clades"
        else:
            coda = ""
        if args.summary_target is None:
            args.summary_target = "consensus"
        if args.summary_target == "consensus":
            if args.min_clade_freq is None:
                min_freq = constants.GREATER_THAN_HALF
            else:
                min_freq = args.min_clade_freq
            tree = tree_array.consensus_tree(min_freq=min_freq, summarize_splits=False)
            msg = "Summarized onto consensus tree with minimum clade frequency threshold of {}:".format(min_freq)
        elif args.summary_target == "mcct" or args.summary_target == "mcc":
            tree = tree_array.maximum_product_of_split_support_tree(
                    include_external_splits=args.include_external_splits_when_scoring_clade_credibility_tree,
                    summarize_splits=False)
            msg = "Summarized onto Maximum Credibility Tree (i.e., tree given in sources that maximizes the product of clade credibilities{}):".format(coda)
        elif args.summary_target == "msct":
            tree = tree_array.maximum_sum_of_split_support_tree(
                    include_external_splits=args.include_external_splits_when_scoring_clade_credibility_tree,
                    summarize_splits=False)
            msg = "Summarized onto Maximum Sum of Credibilities Tree (i.e., tree given in sources that maximizes the sum of clade credibilities{}):".format(coda)
        else:
            raise ValueError(args.summary_target)
        _message_and_log(msg, wrap=True)
        target_trees.append(tree)
    else:
        try:
            if not args.allow_unknown_target_tree_taxa:
                tree_array.taxon_namespace.is_mutable = False
            # we go through the yielder because it can handle the 'nexus/newick'
            # schema; TreeList.get_from_*() etc. does not (yet)
            is_target_trees_rooted = None
            for tree_idx, tree in enumerate(dendropy.Tree.yield_from_files(
                    files=[target_tree_filepath],
                    schema=args.input_format,
                    rooting=dendropy.get_rooting_argument(is_rooted=args.is_source_trees_rooted),
                    preserve_underscores=args.preserve_underscores,
                    taxon_namespace=target_trees.taxon_namespace,
                    )):
                if args.root_target_at_outgroup is not None or args.root_target_at_midpoint:
                    tree.is_rooted = True
                if tree.is_rooted is not tree_array.is_rooted_trees:
                    messenger.error("Target trees rooting state do not match source trees rooting state. " + mixed_rooting_solution)
                    sys.exit(1)
                if tree_idx > 0:
                    if tree.is_rooted is not is_target_trees_rooted:
                        messenger.error("Mixed rooting states detected in target trees. " + mixed_rooting_solution)
                        sys.exit(1)
                is_target_trees_rooted = tree.is_rooted
                tree.encode_bipartitions()
                target_trees.append(tree)
        except (Exception, KeyboardInterrupt) as e:
            if isinstance(e, dendropy.utility.error.ImmutableTaxonNamespaceError):
                message = "Target trees have one or more taxon names not seen in sources: {}".format(e)
            else:
                message = str(e)
            messenger.error(message)
            if args.debug_mode:
                raise
            sys.exit(1)
        if len(target_trees) > 1:
            msg = "Summarizing onto {} target trees".format(len(target_trees))
        else:
            msg = "Summarizing onto target tree".format(len(target_trees))
        msg += " defined in '{}':".format(target_tree_filepath)
        _message_and_log(msg, wrap=False)

    ###  set up summarization regime

    split_summarization_kwargs = {}
    if not args.support_as_percentages:
        _bulleted_message_and_log("Support values expressed as proportions or probabilities")
        if args.support_label_decimals < 2:
            messenger.warning("Reporting support by proportions require that support will be reported to at least 2 decimal places")
            args.support_label_decimals = 2
    else:
        _bulleted_message_and_log("Support values expressed as percentages")
    split_summarization_kwargs["support_as_percentages"] = args.support_as_percentages
    split_summarization_kwargs["support_label_decimals"] = args.support_label_decimals
    if args.support_as_labels:
        split_summarization_kwargs["set_support_as_node_label"] = True
    if args.node_labels == "support":
        split_summarization_kwargs["set_support_as_node_label"] = True
    else:
        split_summarization_kwargs["set_support_as_node_label"] = False

    if args.edge_length_summarization is None:
        if target_tree_filepath:
            args.edge_length_summarization = "keep"
        elif args.summarize_node_ages:
            args.edge_length_summarization = "mean-age"
        else:
            args.edge_length_summarization = "mean-length"
    if args.edge_length_summarization == "mean-length":
        _bulleted_message_and_log("Edge lengths on target trees set to mean of edge lengths in sources")
    elif args.edge_length_summarization == "median-length":
        _bulleted_message_and_log("Edge lengths on target trees set to median of edge lengths in sources")
    elif args.edge_length_summarization == "mean-age":
        _bulleted_message_and_log("Node ages on target trees set to mean of node ages in sources")
    elif args.edge_length_summarization == "median-age":
        _bulleted_message_and_log("Node ages on target trees set to median of node ages in sources")
    elif args.edge_length_summarization == "support":
        _bulleted_message_and_log("Edge lengths on target trees set to support values of corresponding split")
    elif args.edge_length_summarization == "keep":
        _bulleted_message_and_log("Edge lengths as given on target trees")
    elif args.edge_length_summarization == "clear":
        _bulleted_message_and_log("Edge lengths cleared from target trees")
    else:
        raise ValueError(args.edge_length_summarization)
    split_summarization_kwargs["set_edge_lengths"] = args.edge_length_summarization

    split_summarization_kwargs["error_on_negative_edge_lengths"] = False
    if args.collapse_negative_edges:
        split_summarization_kwargs["minimum_edge_length"] = 0.0
        _bulleted_message_and_log("Negative edge lengths collapsed to 0.0 (may result in non-ultrametric trees)")
    elif args.force_minimum_edge_length is not None:
        split_summarization_kwargs["minimum_edge_length"] = args.force_minimum_edge_length
        _bulleted_message_and_log("Edge lengths less than {val} set to {val} (may result in non-ultrametric trees)".format(val=args.force_minimum_edge_length))

    if args.suppress_annotations:
        split_summarization_kwargs["add_support_as_node_attribute"] = False
        split_summarization_kwargs["add_support_as_node_annotation"] = False
        split_summarization_kwargs["add_node_age_summaries_as_node_attributes"] = False
        split_summarization_kwargs["add_node_age_summaries_as_node_annotations"] = False
        split_summarization_kwargs["add_edge_length_summaries_as_edge_attributes"] = False
        split_summarization_kwargs["add_edge_length_summaries_as_edge_annotations"] = False
        _bulleted_message_and_log("Metadata annotations NOT added to target trees as metadata".format())
    else:
        split_summarization_kwargs["add_support_as_node_attribute"] = True
        split_summarization_kwargs["add_support_as_node_annotation"] = True
        split_summarization_kwargs["add_node_age_summaries_as_node_attributes"] = True
        split_summarization_kwargs["add_node_age_summaries_as_node_annotations"] = True
        split_summarization_kwargs["add_edge_length_summaries_as_edge_attributes"] = True
        split_summarization_kwargs["add_edge_length_summaries_as_edge_annotations"] = True
        _bulleted_message_and_log("Support and other summarization annotations added to target trees as metadata".format())

    for tree in target_trees:
        tree_array.summarize_splits_on_tree(
                tree=tree,
                is_bipartitions_updated=True,
                **split_summarization_kwargs)
        if args.node_labels == "clear":
            for nd in tree:
                nd.label = None

    # collapse if below minimum threshold
    if args.min_clade_freq is not None and args.summary_target != "consensus":
        msg = "Collapsing clades or splits with support frequency less than {}".format(args.min_clade_freq)
        for tree in target_trees:
            tree_array.collapse_edges_with_less_than_minimum_support(
                    tree=tree,
                    min_freq=args.min_clade_freq,)
        _message_and_log(msg, wrap=False)

    ###  rooting

    if args.root_target_at_outgroup is not None or args.set_outgroup is not None:
        if args.root_target_at_outgroup is not None:
            outgroup_label = args.root_target_at_outgroup
        elif args.set_outgroup is not None:
            outgroup_label = args.set_outgroup
        if args.input_format in ("nexus/newick", "nexus", "newick"):
            if not args.preserve_underscores:
                outgroup_label = outgroup_label.replace("_", " ")
        for tree in target_trees:
            outgroup_node = tree.find_node_with_taxon_label(outgroup_label)
            if outgroup_node is None:
                messenger.error("Cannot locate node with outgroup taxon '{}' on target tree".format(outgroup_label))
                sys.exit(1)
            tree.to_outgroup_position(
                    outgroup_node=outgroup_node,
                    update_bipartitions=True,
                    suppress_unifurcations=True)
            if args.root_target_at_outgroup is not None:
                tree.is_rooted = True
        if args.root_target_at_outgroup is not None:
            _bulleted_message_and_log("Target tree(s) rerooted using outgroup: '{}'".format(outgroup_label))
        elif args.set_outgroup is not None:
            _bulleted_message_and_log("Target tree(s) rotated to set outgroup: '{}'".format(outgroup_label))
    elif args.root_target_at_midpoint:
        for tree in target_trees:
            tree.reroot_at_midpoint(
                    update_bipartitions=True,
                    suppress_unifurcations=True)
        _bulleted_message_and_log("Target tree(s) rerooted at midpoint")

    main_time_end = datetime.datetime.now()

    ###################################################
    #  Primary Output

    ## set up file-level annotations
    final_run_report = []
    final_run_report.append("Started at: {}".format(main_time_start.isoformat(' ')))
    final_run_report.append("Ended at: {}".format(main_time_end.isoformat(' ')))
    final_run_report.append("Total elapsed time: {}".format(
        timeprocessing.pretty_timedelta(main_time_end-main_time_start),
        ))
    final_run_report.append("Actual analysis time: {}".format(
        timeprocessing.pretty_timedelta(analysis_time_delta),
        ))

    if not args.suppress_analysis_metainformation:
        summarization_metainfo = []

        summarization_metainfo.append("")
        summarization_metainfo.append("Summarization Information")
        summarization_metainfo.append("-------------------------")
        summarization_metainfo.extend(processing_report_lines)

        summarization_metainfo.append("")
        summarization_metainfo.append("Program Information")
        summarization_metainfo.append("-------------------")
        summarization_metainfo.append("{} {} by {}".format(_program_name, _program_version, _program_author))
        summarization_metainfo.append(dendropy.description_text())

        summarization_metainfo.append("Execution Information")
        summarization_metainfo.append("---------------------")
        try:
            username = getpass.getuser()
        except:
            username = "<user>"
        summarization_metainfo.append("Executed on {} by {}@{}".format(platform.node(), username, socket.gethostname()))
        summarization_metainfo.append("Working directory: '{}'".format(os.getcwd()))
        summarization_metainfo.extend(final_run_report)

        summarization_metainfo.append("")
        summarization_metainfo.append("Citation Information")
        summarization_metainfo.append("--------------------")
        summarization_metainfo.append("")
        citation = cli.compose_citation_for_program(
                prog_name=_program_name,
                prog_version=_program_version,
                additional_citations=[_program_citation],
                include_preamble=False,
                include_epilog=False,
                )
        summarization_metainfo.extend(citation)
        summarization_metainfo.append("")

        if args.additional_comments:
            summarization_metainfo.append("")
            summarization_metainfo.append("Additional Remarks")
            summarization_metainfo.append("------------------")
            summarization_metainfo.append(args.additional_comments)
    else:
        summarization_metainfo = []

    ### PRIMARY OUTPUT
    if not args.suppress_analysis_metainformation:
        primary_output_metainfo = []
        primary_output_metainfo.append("=============")
        primary_output_metainfo.append("Summary Trees")
        primary_output_metainfo.append("=============")
        primary_output_metainfo.append("")
        primary_output_metainfo.append("Summary trees generated by SumTrees.")
        primary_output_metainfo.extend(summarization_metainfo)
    else:
        primary_output_metainfo = []
    if hasattr(output_dest, "name"):
        messenger.info("Writing primary results to: '{}'".format(output_dest.name))
    else:
        messenger.info("Writing primary results to standard output".format(output_dest.name))
    _write_trees(trees=target_trees,
            output_dest=output_dest,
            args=args,
            file_comments=primary_output_metainfo)

    ### EXTENDED OUTPUT
    if extended_output_paths:

        messenger.info("Calculating extended summarization results")

        #### get data: topologies
        topologies = tree_array.topologies(
                sort_descending=True,
                frequency_attr_name="frequency",
                frequency_annotation_name="frequency",
                )

        #### get data: bipartitions
        all_taxa_bitmask = tree_array.taxon_namespace.all_taxa_bitmask()
        seen_split_bitmasks = set()
        all_bipartitions = collections.OrderedDict()
        bipartition_table = []
        bipartitions_as_trees = dendropy.TreeList(taxon_namespace=tree_array.taxon_namespace)
        # bipartition_stats_fieldname_map
        # biparitition_table_fieldnames = [
        #         "bipartitionId",
        #         "bipartitionGroup",
        #         "frequency",
        # ]
        # for stat_fieldname in SplitDistribution.SUMMARY_STATS_FIELDNAMES:
        #     f = textprocessing.camel_case("{}_{}".format(summary_stat_prefix, stat_fieldname))
        #     bipartition_table_fieldnames.append(f)
        _inf = float("inf")
        def _add_split_bitmask_data(split_bitmask):

            # do not add if already accessioned
            if split_bitmask in seen_split_bitmasks:
                return
            seen_split_bitmasks.add(split_bitmask)

            # create bipartition from split
            bipartition = dendropy.Bipartition(
                    leafset_bitmask=split_bitmask,
                    tree_leafset_bitmask=all_taxa_bitmask,
                    is_rooted=tree_array.is_rooted_trees,
                    is_mutable=False,
                    compile_bipartition=True)
            bipartition_newick_str = bipartition.leafset_as_newick_string(
                    tree_array.taxon_namespace,
                    preserve_spaces=True if args.preserve_underscores else False,
                    quote_underscores=False if args.preserve_underscores else True,
                    )

            # bipartition table
            bipartition_data = collections.OrderedDict()
            bipartition_data["bipartitionGroup"] = bipartition.leafset_as_bitstring(
                    symbol0=".",
                    symbol1="*",
                    reverse=True,
                    )
            bipartition_data["bipartitionId"] = bipartition.split_bitmask
            bipartition_data["bipartitionBitmask"] = bipartition.split_as_bitstring(
                    symbol0="0",
                    symbol1="1",
                    reverse=False,
                    )
            bipartition_data["bipartitionLeafset"] = bipartition.leafset_as_bitstring(
                    symbol0="0",
                    symbol1="1",
                    reverse=False,
                    )
            bipartition_data["count"] = tree_array.split_distribution.split_counts[split_bitmask]
            bipartition_data["frequency"] = tree_array.split_distribution[split_bitmask]
            for summary_stat_prefix, summary_source in (
                    ("edge_length", tree_array.split_distribution.split_edge_length_summaries),
                    ("node_age", tree_array.split_distribution.split_node_age_summaries),
                    ):
                if not summary_source:
                    continue
                for stat_fieldname in dendropy.SplitDistribution.SUMMARY_STATS_FIELDNAMES:
                    f = textprocessing.camel_case("{}_{}".format(summary_stat_prefix, stat_fieldname))
                    if split_bitmask in summary_source:
                        value = summary_source[split_bitmask].get(stat_fieldname, 0.0)
                    else:
                        value = None
                    if value is None:
                        if stat_fieldname in ("hpd95", "quant_5_95", "range"):
                            value = (0.0, 0.0)
                        else:
                            value = 0.0
                    elif value == _inf:
                        value = 0.0
                    if isinstance(value, list) or isinstance(value, tuple):
                        for sub_f, sub_value in zip(("Min", "Max"), sorted(value)):
                            bipartition_data[f+sub_f] = sub_value
                    else:
                        bipartition_data[f] = value
            bipartition_data["newick"] = '"{}"'.format(bipartition_newick_str)
            bipartition_table.append(bipartition_data)

            # bipartition as tree
            tree = dendropy.Tree.get_from_string(
                    bipartition_newick_str,
                    "newick",
                    taxon_namespace=tree_array.taxon_namespace,
                    rooting=dendropy.get_rooting_argument(is_rooted=tree_array.is_rooted_trees),
                    extract_comment_metadata=False,
                    )
            tree.label = "Bipartition{}".format(bipartition_data["bipartitionId"])
            # tree.label = "Bipartition{}".format(bipartition.split_as_bitstring())
            tree.weight = bipartition_data["frequency"]
            # tree.seed_node.annotations.add_new("bipartitionId",
            #         '"{}"'.format(bipartition_data["bipartitionId"]))
            # tree_array.summarize_splits_on_tree(
            #         tree=tree,
            #         is_bipartitions_updated=False,
            #         **split_summarization_kwargs)
            for key in bipartition_data:
                if key in ("newick", ):
                    continue
                value = bipartition_data[key]
                if key in ("bipartitionId", "bipartitionBitmask", "bitpartitionLeafset"):
                    # FigTree cannot cast bigger integers values to float
                    value = '"{}"'.format(value)
                tree.seed_node.annotations.add_new(
                        textprocessing.snake_case(key),
                        value)
            bipartitions_as_trees.append(tree)

            all_bipartitions[bipartition] = bipartition_data
            return bipartition

        # this is to preserve order seen in Mr. Bayes
        _add_split_bitmask_data(all_taxa_bitmask)
        for taxon in tree_array.taxon_namespace:
            split_bitmask = tree_array.taxon_namespace.taxon_bitmask(taxon)
            _add_split_bitmask_data(split_bitmask)

        # add the rest in order
        sd_split_bitmasks = list(tree_array.split_distribution.split_counts.keys())
        sd_split_bitmasks.sort(key=lambda x: tree_array.split_distribution.split_counts[x], reverse=True)
        for split_bitmask in sd_split_bitmasks:
            _add_split_bitmask_data(split_bitmask)

        #### EXTENDED OUTPUT: topologies / trprobs
        if not args.suppress_analysis_metainformation:
            metainfo = []
            metainfo.append("======================")
            metainfo.append("Topology Probabilities")
            metainfo.append("======================")
            metainfo.append("")
            metainfo.append("\n".join((
                    "Topologies in the source set of trees, listing in",
                    "descending order of frequency with an indication ",
                    "of their individual frequencies ('frequency') and",
                    "cumulative frequencies ('cumulative_frequency'). ",
                    )))
            metainfo.extend(summarization_metainfo)
        else:
            metainfo = []
        output_path = extended_output_paths["topologies"]
        messenger.info("Writing topologies to: '{}'".format(output_path))
        cumulative_frequency = 0.0
        for tree in topologies:
            tree.weight = tree.frequency
            cumulative_frequency += tree.frequency
            tree.cumulative_frequency = cumulative_frequency
            tree.annotations.add_bound_attribute("cumulative_frequency")
            # tree_array.summarize_splits_on_tree(
            #         tree=tree,
            #         is_bipartitions_updated=True,
            #         support_as_percentages=args.support_as_percentages,
            #         support_label_decimals=args.support_as_percentages,
            #         add_support_as_node_annotation=not args.suppress_annotations,
            #         add_node_age_summaries_as_node_attributes=False,
            #         add_node_age_summaries_as_node_annotations=False,
            #         add_edge_length_summaries_as_edge_attributes=False,
            #         add_edge_length_summaries_as_edge_annotations=False,
            #         )
        with open(output_path, "w") as out:
            _write_trees(trees=topologies,
                    output_dest=out,
                    args=args,
                    file_comments=metainfo)

        #### EXTENDED OUTPUT: bipartition trees
        if not args.suppress_analysis_metainformation:
            metainfo = []
            metainfo.append("============")
            metainfo.append("Bipartitions")
            metainfo.append("============")
            metainfo.append("")
            metainfo.append("\n".join((
                    "Bipartitions in the source set of trees, represented ",
                    "as trees, with information summarized from the source ",
                    "set of trees annotated as metadata.",
                    )))
            metainfo.extend(summarization_metainfo)
        else:
            metainfo = []
        output_path = extended_output_paths["bipartition-trees"]
        messenger.info("Writing bipartition trees to: '{}'".format(output_path))
        with open(output_path, "w") as out:
            _write_trees(trees=bipartitions_as_trees,
                    output_dest=out,
                    args=args,
                    file_comments=metainfo)

        #### EXTENDED OUTPUT: bipartition table
        output_path = extended_output_paths["bipartition-table"]
        messenger.info("Writing bipartition table to: '{}'".format(output_path))
        sample_row = list(bipartition_table[0].keys())
        with open(output_path, "w") as out:
            writer = csv.DictWriter(
                    out,
                    fieldnames=sample_row,
                    lineterminator=os.linesep,
                    delimiter="\t",
                    )
            writer.writeheader()
            writer.writerows(bipartition_table)

        #### EXTENDED OUTPUT: edge lengths and node ages
        bipartition_table_keys_to_import = [
                "bipartitionGroup",
                "bipartitionId",
                "bipartitionBitmask",
                "bipartitionLeafset",
                "frequency",
                ]

        #### EXTENDED OUTPUT: edge lengths
        if tree_array.split_distribution.split_edge_lengths:
            rows = []
            for b in all_bipartitions:
                entry = collections.OrderedDict()
                for key in bipartition_table_keys_to_import:
                    entry[key] = all_bipartitions[b][key]
                # entry["bipartitionId"] = all_bipartitions[b]["bipartitionId"]
                try:
                    value_list = tree_array.split_distribution.split_edge_lengths[b.split_bitmask]
                except KeyError:
                    value_list = []
                entry["edgeCount"] = len(value_list)
                entry["edgeLengths"] = ",".join(str(v) for v in value_list)
                rows.append(entry)
            output_path = extended_output_paths["edge-lengths"]
            messenger.info("Writing edge set to: '{}'".format(output_path))
            with open(output_path, "w") as out:
                writer = csv.DictWriter(
                        out,
                        fieldnames=bipartition_table_keys_to_import + ["edgeCount", "edgeLengths"],
                        lineterminator=os.linesep,
                        delimiter="\t",
                        )
                writer.writeheader()
                writer.writerows(rows)

        #### EXTENDED OUTPUT: node ages
        if tree_array.split_distribution.split_node_ages:
            rows = []
            for b in all_bipartitions:
                entry = collections.OrderedDict()
                for key in bipartition_table_keys_to_import:
                    entry[key] = all_bipartitions[b][key]
                try:
                    value_list = tree_array.split_distribution.split_node_ages[b.split_bitmask]
                except KeyError:
                    value_list = []
                entry["nodeCount"] = len(value_list)
                entry["nodeAges"] = ",".join(str(v) for v in value_list)
                rows.append(entry)
            output_path = extended_output_paths["node-ages"]
            messenger.info("Writing edge set to: '{}'".format(output_path))
            with open(output_path, "w") as out:
                writer = csv.DictWriter(
                        out,
                        fieldnames=bipartition_table_keys_to_import + ["nodeCount", "nodeAges"],
                        lineterminator=os.linesep,
                        delimiter="\t",
                        )
                writer.writeheader()
                writer.writerows(rows)

    ###################################################
    #  WRAP UP

    messenger.info("Summarization completed")
    messenger.info_lines(final_run_report)
    messenger.silent = True

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\n(Terminating due to user interrupt signal)")
