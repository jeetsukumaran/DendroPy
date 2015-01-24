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
from dendropy.utility import error
from dendropy.utility import messaging
from dendropy.utility import timeprocessing

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
        def _log_progress(source_name, current_tree_offset, aggregate_tree_idx):
            if (
                    info_message_func is not None
                    and (
                        (log_frequency == 1)
                        or (tree_offset > 0 and current_tree_offset == tree_offset)
                        or (aggregate_tree_idx > 0 and log_frequency > 0 and (aggregate_tree_idx % log_frequency) == 0)
                        )
                    ):
                if current_tree_offset >= tree_offset:
                    coda = " (processing)"
                else:
                    coda = " (skipping)"
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
                        info_message_func("Processing {} of {}: '{}'".format(current_source_index+1, len(tree_sources), source_name), wrap=False)
                    else:
                        info_message_func("Processing: '{}'".format(source_name), wrap=False)
                if current_tree_offset >= tree_offset:
                    tree_array.add_tree(tree=tree, is_bipartitions_updated=False)
                    _log_progress(source_name, current_tree_offset, aggregate_tree_idx)
                else:
                    _log_progress(source_name, current_tree_offset, aggregate_tree_idx)
                current_tree_offset += 1
        except (Exception, KeyboardInterrupt) as e:
            e.exception_tree_source_name = tree_yielder.current_file_name
            e.exception_tree_offset = current_tree_offset
            raise e

class TreeProcessingWorker(multiprocessing.Process):

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
            num_processes,
            log_frequency,
            messenger,
            messenger_lock,
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
        self.tree_array.worker_name = self.name
        self.num_tasks_received = 0
        self.num_tasks_completed = 0

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
            num_processes,
            log_frequency,
            messenger,
            ):
        self.is_source_trees_rooted = is_source_trees_rooted
        self.rooting_interpretation = dendropy.get_rooting_argument(is_rooted=self.is_source_trees_rooted)
        self.ignore_edge_lengths = ignore_edge_lengths
        self.ignore_node_ages = ignore_node_ages
        self.use_tree_weights = use_tree_weights
        self.ultrametricity_precision = ultrametricity_precision
        self.num_processes = num_processes
        self.log_frequency = log_frequency
        self.messenger = messenger

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
                )
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
        self.info_message("Launching {} worker processes".format(self.num_processes))
        results_queue = multiprocessing.Queue()
        messenger_lock = multiprocessing.Lock()
        workers = []
        for idx in range(self.num_processes):
            # self.info_message("Launching {} of {} worker processes".format(idx+1, self.num_processes))
            tree_processing_worker = TreeProcessingWorker(
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
                    num_processes=self.num_processes,
                    messenger=self.messenger,
                    messenger_lock=messenger_lock,
                    log_frequency=self.log_frequency)
            tree_processing_worker.start()
            workers.append(tree_processing_worker)

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
## Preprocessing

def preprocess_tree_sources(args, messenger):
    tree_sources = []
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
            return []
        else:
            if args.input_format is None:
                args.input_format = "nexus/newick"
            else:
                args.input_format = args.input_format.lower()
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

    ######################################################################
    ## Start Recording Total Job Time

    main_start = datetime.datetime.now()

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
            help="Format of all input trees (defaults to handling either NEXUS or NEWICK through inspection; it is more efficient to explicitly specify the format if it is known).")
    source_options.add_argument("-b", "--burnin",
            type=int,
            default=0,
            help="Number of trees to skip from the beginning of *each* tree file when counting support (default: %(default)s).")
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
    source_options.add_argument("-u", "--ultrametric",
            dest="is_source_trees_ultrametric",
            action="store_true",
            default=None,
            help="Assume source trees are ultrametric (implies '--force-rooted'; will result in node ages being summarized; will result in error if trees are not ultrametric).")
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
                " multiprocessing ('-M' or '-m') is requested. When"
                " parallel processing using the '-M' or '-m'"
                " options, all taxon names need to be defined in advance"
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
                "(default: > 0.5)."))
    target_tree_supplemental_options.add_argument("--allow-unknown-target-tree-taxa",
            action="store_true",
            default=False,
            help=(
                "Do not fail with error if target tree(s) have taxa not"
                " previously encountered in source trees or defined in"
                " the taxon discovery file."
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
    edge_summarization_choices = ["mean-length", "median-length", "mean-age", "median-age", "support", "keep", "clear",]
    edge_summarization_options.add_argument("-e", "--edges",
            dest="edge_summarization",
            # metavar="<%s>" % ("".join(edge_summarization_choices)),
            choices=edge_summarization_choices,
            default=None,
            help="\n".join((
                "R}Set the edge lengths of the summary/target tree(s):",
                "- 'mean-length'    : Edge lengths will be set to the mean",
                "                     of the lengths of the corresponding",
                "                     split or clade in the source trees.",
                "                     If no external summary tree targets",
                "                     are specified and the input source",
                "                     trees are not ultrametric, then this",
                "                     is the DEFAULT option.",
                "- 'median-length'  : Edge lengths will be set to the median",
                "                     of the lengths of the corresponding",
                "                     split or clade in the source trees.",
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
                "                     then this is the DEFAULT option.",
                "- 'median-age'     : Edge lengths will be adjusted so that",
                "                     the age of subtended nodes will be equal",
                "                     to the median age of the corresponding",
                "                     split or clade in the source trees. This",
                "                     option requires that the source trees",
                "                     are ultrametric (i.e., '-u' or",
                "                     '--ultrametric' must be specified).",
                "- 'support'        : Edge lengths will be set to the support",
                "                     value for the split represented by the ",
                "                     edge.",
                "- 'keep'           : Do not change the existing edge lengths.",
                "                     This is the DEFAULT if target tree(s) are",
                "                     sourced from an external file using the",
                "                     '-t' or '--target' option",
                "- 'clear'          : Edge lengths will be cleared from the",
                "                     target trees if they are present.",
                "Note the default settings varies depending according to the",
                "following, in order of preference:",
                "(1) If target trees are specified using the '-t' or '--target'",
                "    option, then the default is: 'keep'.",
                "(2) If no target trees are specified, but the source trees ",
                "    are specified to be ultrametric using the '-u' or ",
                "    '--ultrametric' option, then the default is: 'mean-age'.",
                "(3) If no target trees are specified and the source trees ",
                "    are NOT specified to be ultrametric then the default ",
                "    is: 'mean-length'.",
                )))
    edge_summarization_options.add_argument("--collapse-negative-edges",
            action="store_true",
            default=False,
            help="(If setting edge lengths) force parent node ages to be at least as old as its oldest child when summarizing node ages.")

    node_summarization_options = parser.add_argument_group("Target Tree Node Options")
    node_summarization_options.add_argument("-l","--labels",
            dest="node_labels",
            default="support",
            choices=["support", "keep", "clear",],
            help="\n".join((
                "R}Set the node labels of the summary target tree(s):",
                "- 'support'        : Node labels will be set to the support",
                "                     value for the clade represented by the ",
                "                     node. This is the DEFAULT.",
                "- 'keep'           : Do not change the existing node labels.",
                "- 'clear'          : Node labels will be cleared from the",
                "                     target trees if they are present.",
                )))
    node_summarization_options.add_argument("--no-node-annotations",
            action="store_true",
            default=False,
            help=(
                "Do NOT annotate nodes and edges with any summarization information metadata such as."
                "support values, edge length and/or node age summary statistcs, etc."
                ))

    support_summarization_options = parser.add_argument_group("Support Summarization Options")
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

    parser.add_argument('--debug-mode', action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args()

    ######################################################################
    ## Information (Only) Operations

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
    else:
        target_tree_filepath = None
        if args.summary_tree_target is None:
            args.summary_tree_target = "consensus"

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

    if args.taxon_name_file is not None:
        with open(os.path.expanduser(os.path.expandvars(args.taxon_name_file)), "r") as tnf:
            taxon_labels = [name.strip() for name in tnf.read().split("\n") if name]
            taxon_labels = [name for name in taxon_labels if name]
        taxon_namespace = dendropy.TaxonNamespace(taxon_labels)
    else:
        taxon_namespace = None

    ######################################################################
    ## Main Work

    tree_processor = TreeProcessor(
            is_source_trees_rooted=args.is_source_trees_rooted,
            ignore_edge_lengths=not args.summarize_edge_lengths,
            ignore_node_ages=not args.summarize_node_ages,
            use_tree_weights=args.weighted_trees,
            ultrametricity_precision=args.ultrametricity_precision,
            num_processes=num_processes,
            log_frequency=args.log_frequency if not args.quiet else 0,
            messenger=messenger,
            )
    processing_time_start = datetime.datetime.now()
    # messenger.info("Processing of source trees starting at {}".format(
    #     processing_time_start,
    #     ))
    try:
        tree_array = tree_processor.process_trees(
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
    processing_time_end = datetime.datetime.now()
    messenger.info("Processing of source trees completed in: {}".format(
        timeprocessing.pretty_timedelta(processing_time_end-processing_time_start),
        wrap=False,
        ))

    ######################################################################
    ## Post-Processing

    ### post-analysis reports
    report_lines = []
    def _report(msg):
        messenger.info(msg)
        report_lines.append(msg)

    _report("{} trees considered in total for summarization".format(len(tree_array)))
    if args.weighted_trees:
        _report("Trees were treated as weighted (default weight = 1.0).")
    else:
        _report("Trees were treated as unweighted")
    _report("{} unique taxa across all trees".format(len(tree_array.taxon_namespace)))
    if args.is_source_trees_rooted is None:
        if tree_array.split_distribution.is_all_counted_trees_rooted():
            _report("All trees were rooted")
        elif tree_array.split_distribution.is_all_counted_trees_strictly_unrooted():
            _report("All trees were unrooted")
        elif tree_array.split_distribution.is_all_counted_trees_treated_as_unrooted():
            _report("All trees were assumed to be unrooted")
    elif args.is_source_trees_rooted is True:
        _report("All trees were treated as rooted")
    else:
        _report("All trees were treated as unrooted")
    if args.is_source_trees_ultrametric and args.ultrametricity_precision:
        _report("Trees were ultrametric within an error of {}".format(args.ultrametricity_precision))
    elif args.is_source_trees_ultrametric:
        _report("Trees were expected to be ultrametric (not verified)")
    num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits = tree_array.split_distribution.splits_considered()
    _report("{} unique splits counted".format(num_unique_splits))
    _report("{} unique non-trivial splits counted".format(num_nt_unique_splits))

    ### build target tree
    if target_tree_filepath is None:
        pass
    else:
        try:
            if not args.allow_unknown_target_tree_taxa:
                tree_array.taxon_namespace.is_mutable = False
            # we go through the yielder because it can handle the 'nexus/newick'
            # schema; TreeList.get_from_*() etc. does not (yet)
            target_trees = dendropy.TreeList(taxon_namespace=tree_array.taxon_namespace)
            is_target_trees_rooted = None
            for tree_idx, tree in enumerate(dendropy.Tree.yield_from_files(
                    files=[target_tree_filepath],
                    schema=args.input_format,
                    rooting=dendropy.get_rooting_argument(is_rooted=args.is_source_trees_rooted),
                    preserve_underscores=args.preserve_underscores,
                    taxon_namespace=target_trees.taxon_namespace,
                    )):
                if tree.is_rooted is not tree_array.is_rooted_trees:
                    messenger.error("Target trees rooting state do not match source trees rooting state. " + mixed_rooting_solution)
                    sys.exit(1)
                if tree_idx > 0:
                    if tree.is_rooted is not is_target_trees_rooted:
                        messenger.error("Mixed rooting states detected in target trees. " + mixed_rooting_solution)
                        sys.exit(1)
                is_target_trees_rooted = tree.is_rooted
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
        msg += " defined in: '{}'".format(target_tree_filepath)
        messenger.info(msg, wrap=False)

    ###  set up summarization regime

    split_summarization_kwargs = {}
    if not args.support_as_percentages:
        _report("Support values expressed as percentages")
        if args.support_label_decimals < 2:
            messenger.warning("Reporting support by proportions require that support will be reported to at least 2 decimal places")
            args.support_label_decimals = 2
    else:
        _report("Support values expressed as proportions or probabilities")
    split_summarization_kwargs["support_as_percentages"] = args.support_as_percentages
    split_summarization_kwargs["support_label_decimals"] = args.support_label_decimals

    if args.edge_summarization is None:
        pass
    if args.edge_summarization == "mean_length":
        _report("Edge lengths on target trees set to mean of edge lengths in sources")
    elif args.edge_summarization == "median_length":
        _report("Edge lengths on target trees set to median of edge lengths in sources")
    elif args.edge_summarization == "mean_age":
        _report("Node ages on target trees set to mean of node ages in sources")
    elif args.edge_summarization == "median_age":
        _report("Node ages on target trees set to median of node ages in sources")
    elif args.edge_summarization == "support":
        _report("Edge lengths on target trees set to support values of corresponding split")
    elif args.edge_summarization == "keep":
        _report("Edge lengths on target trees are not modified")
    elif args.edge_summarization == "clear":
        _report("Edge lengths on target trees are cleared")
    else:
        raise ValueError(args.edge_summarization)
    split_summarization_kwargs["set_edge_lengths"] = args.edge_summarization

    # split_summarization_kwargs["set_edge_lengths"] = None
    # split_summarization_kwargs["add_support_as_node_attribute"] = None
    # split_summarization_kwargs["add_support_as_node_annotation"] = None
    # split_summarization_kwargs["set_support_as_node_label"] = None
    # split_summarization_kwargs["add_node_age_summaries_as_node_attributes"] = None
    # split_summarization_kwargs["add_node_age_summaries_as_node_annotations"] = None
    # split_summarization_kwargs["add_edge_length_summaries_as_edge_attributes"] = None
    # split_summarization_kwargs["add_edge_length_summaries_as_edge_annotations"] = None
    # split_summarization_kwargs["support_label_decimals"] = None
    # split_summarization_kwargs["support_as_percentages"] = None
    # split_summarization_kwargs["support_label_compose_func"] = None
    # split_summarization_kwargs["primary_fieldnames"] = None
    # split_summarization_kwargs["summary_stats_fieldnames"] = None
    # split_summarization_kwargs["node_age_summaries_fieldnames"] = None
    # split_summarization_kwargs["edge_length_summaries_fieldnames"] = None
    # split_summarization_kwargs["fieldnames"] = None


    # root target tree(s)
    # decorate target tree(s)
    # other stuff



    # comments = []
    # comments.extend(report)
    # messenger.info("Split counting completed:")
    # messenger.info_lines(report, prefix="  ~ ")





    real_value_format_specifier = None

    t = tree_array.consensus_tree()
    # print(t.as_string("nexus"))
    print(len(tree_array))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\n(Terminating due to user interrupt signal)")
