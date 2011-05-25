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
CLI wrapper for tree summarization.
"""

import os
import sys
import textwrap
from optparse import OptionParser
from optparse import OptionGroup

import datetime
import time
import socket
try:
    import getpass
except:
    pass
import platform
try:
    import Queue
    import multiprocessing
    _MP = True
except ImportError:
    _MP = False

import dendropy
from dendropy import treesplit
from dendropy import treesum
from dendropy.dataio import tree_source_iter
from dendropy.dataio import multi_tree_source_iter
from dendropy.dataio import newick
from dendropy.utility.messaging import ConsoleMessenger
from dendropy.utility.cli import confirm_overwrite, show_splash
from dendropy.utility import statistics

_program_name = "SumTrees"
_program_subtitle = "Phylogenetic Tree Split Support Summarization"
_program_date = "May 05 2011"
_program_version = "Version 3.3.0 (%s)" % _program_date
_program_author = "Jeet Sukumaran and Mark T. Holder"
_program_contact = "jeetsukumaran@gmail.com"
_program_copyright = "Copyright (C) 2008 Jeet Sukumaran.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \
                 "There is NO WARRANTY,\nto the extent permitted by law."
if _MP:
    class SplitCountingWorker(multiprocessing.Process):

        def __init__(self,
                work_queue,
                result_split_dist_queue,
                result_topology_hash_map_queue,
                schema,
                taxon_labels,
                is_rooted,
                ignore_node_ages,
                calc_tree_probs,
                weighted_trees,
                tree_offset,
                process_idx,
                messenger,
                messenger_lock,
                log_frequency=1000):
            multiprocessing.Process.__init__(self)
            self.work_queue = work_queue
            self.result_split_dist_queue = result_split_dist_queue
            self.result_topology_hash_map_queue = result_topology_hash_map_queue
            self.schema = schema
            self.taxon_labels = list(taxon_labels)
            self.taxon_set = dendropy.TaxonSet(self.taxon_labels)
            self.split_distribution = treesplit.SplitDistribution(taxon_set=self.taxon_set)
            self.split_distribution.is_rooted = is_rooted
            self.split_distribution.ignore_node_ages = ignore_node_ages
            self.is_rooted = is_rooted
            self.calc_tree_probs = calc_tree_probs
            # we do not store the tree directly in the topology hash map, because it cannot be
            # pickled (yet); instead, we store the normalized newick string; if memory usage
            # is a concern (and not speed), we might be try to use the normalized newick topology
            # hash as the primary topology hash function as well, instead of the default split sets
            self.topology_counter = treesum.TopologyCounter(tree_store_func=\
                    treesum.TopologyCounter.normalized_newick_topology_hash)
            self.weighted_trees = weighted_trees
            self.tree_offset = tree_offset
            self.process_idx = process_idx
            self.messenger = messenger
            self.messenger_lock = messenger_lock
            self.log_frequency = log_frequency
            self.kill_received = False

        def send_message(self, msg, level, wrap=True):
            if self.messenger is None:
                return
            if self.messenger.messaging_level > level or self.messenger.silent:
                return
            msg = "Thread %d: %s" % (self.process_idx+1, msg)
            self.messenger_lock.acquire()
            try:
                self.messenger.send(msg, level=level, wrap=wrap)
            finally:
                self.messenger_lock.release()

        def send_info(self, msg, wrap=True):
            self.send_message(msg, ConsoleMessenger.INFO_MESSAGING_LEVEL, wrap=wrap)

        def send_warning(self, msg, wrap=True):
            self.send_message(msg, ConsoleMessenger.WARNING_MESSAGING_LEVEL, wrap=wrap)

        def send_error(self, msg, wrap=True):
            self.send_message(msg, ConsoleMessenger.ERROR_MESSAGING_LEVEL, wrap=wrap)

        def run(self):
            while not self.kill_received:
                try:
                    source = self.work_queue.get_nowait()
                except Queue.Empty:
                    break
                self.send_info("Received task: '%s'." % source, wrap=False)
                fsrc = open(source, "rU")
                for tidx, tree in enumerate(tree_source_iter(fsrc,
                        schema=self.schema,
                        taxon_set=self.taxon_set,
                        as_rooted=self.is_rooted,
                        store_tree_weights=self.weighted_trees)):
                    if tidx >= self.tree_offset:
                        if (self.log_frequency == 1) or (tidx > 0 and self.log_frequency > 0 and tidx % self.log_frequency == 0):
                            self.send_info("(processing) '%s': tree at offset %d" % (source, tidx), wrap=False)
                        treesplit.encode_splits(tree)
                        self.split_distribution.count_splits_on_tree(tree)
                        if self.calc_tree_probs:
                            self.topology_counter.count(tree,
                                    tree_splits_encoded=True)
                    else:
                        if (self.log_frequency == 1) or (tidx > 0 and self.log_frequency > 0 and tidx % self.log_frequency == 0):
                            self.send_info("(processing) '%s': tree at offset %d (skipping)" % (source, tidx), wrap=False)
                    if self.kill_received:
                        break
                if self.kill_received:
                    break
                self.send_info("Completed task: '%s'." % (source), wrap=False)
            if self.kill_received:
                self.send_warning("Terminating in response to kill request.")
            else:
                self.result_split_dist_queue.put(self.split_distribution)
                self.result_topology_hash_map_queue.put(self.topology_counter.topology_hash_map)

def discover_taxa(treefile, schema):
    """
    Reads first tree in treefile, and assumes that is sufficient to populate a
    taxon set object fully, which it then returns.
    """
    if isinstance(treefile, str):
        tdf = open(treefile, "rU")
    else:
        tdf = treefile
    tt = None
    for tree in tree_source_iter(tdf, schema=schema):
        tt = tree
        break
    taxon_set = tt.taxon_set
    return taxon_set

def process_sources_parallel(
        num_processes,
        support_filepaths,
        schema,
        is_rooted,
        ignore_node_ages,
        calc_tree_probs,
        weighted_trees,
        tree_offset,
        log_frequency,
        messenger):
    """
    Returns a SplitDistribution object summarizing all trees found in
    `support_filepaths`.
    """

    # describe
    messenger.send_info("Running in multiprocessing mode (up to %d processes)." % num_processes)
    messenger.send_info("%d sources to be processed." % (len(support_filepaths)))

    # pre-discover taxa
    tdfpath = support_filepaths[0]
    messenger.send_info("Pre-loading taxa based on '%s' ..." % tdfpath)
    taxon_set = discover_taxa(tdfpath, schema)
    taxon_labels = [str(t) for t in taxon_set]
    messenger.send_info("Found %d taxa: [%s]" % (len(taxon_labels), (', '.join(["'%s'" % t for t in taxon_labels]))))

    # load up queue
    messenger.send_info("Creating work queue ...")
    work_queue = multiprocessing.Queue()
    for f in support_filepaths:
        work_queue.put(f)

    # launch processes
    messenger.send_info("Launching worker processes ...")
    result_split_dist_queue = multiprocessing.Queue()
    result_topology_hash_map_queue = multiprocessing.Queue()
    messenger_lock = multiprocessing.Lock()
    for idx in range(num_processes):
        sct = SplitCountingWorker(work_queue,
                result_split_dist_queue=result_split_dist_queue,
                result_topology_hash_map_queue=result_topology_hash_map_queue,
                schema=schema,
                taxon_labels=taxon_labels,
                is_rooted=is_rooted,
                ignore_node_ages=ignore_node_ages,
                calc_tree_probs=calc_tree_probs,
                weighted_trees=weighted_trees,
                tree_offset=tree_offset,
                process_idx=idx,
                messenger=messenger,
                messenger_lock=messenger_lock,
                log_frequency=log_frequency)
        sct.start()

    # collate results
    result_count = 0
    split_distribution = treesplit.SplitDistribution(taxon_set=taxon_set)
    split_distribution.is_rooted = is_rooted
    split_distribution.ignore_node_ages = ignore_node_ages
    topology_counter = treesum.TopologyCounter()
    while result_count < num_processes:
        result_split_dist = result_split_dist_queue.get()
        split_distribution.update(result_split_dist)
        result_topology_hash_map = result_topology_hash_map_queue.get()
        topology_counter.update_topology_hash_map(result_topology_hash_map)
        result_count += 1
    messenger.send_info("Recovered results from all worker processes.")
    return split_distribution, topology_counter

def process_sources_serial(
        support_filepaths,
        schema,
        is_rooted,
        ignore_node_ages,
        calc_tree_probs,
        weighted_trees,
        tree_offset,
        log_frequency,
        messenger):
    """
    Returns a SplitDistribution object summarizing all trees found in
    `support_filepaths`.
    """
    messenger.send_info("Running in serial mode.")
    taxon_set = dendropy.TaxonSet()
    split_distribution = treesplit.SplitDistribution(taxon_set=taxon_set)
    split_distribution.ignore_node_ages = ignore_node_ages
    split_distribution.is_rooted = is_rooted
    topology_counter = treesum.TopologyCounter()

    if support_filepaths is None or len(support_filepaths) == 0:
        messenger.send_info("Reading trees from standard input.")
        srcs = [sys.stdin]
    else:
        messenger.send_info("%d source(s) to be processed." % len(support_filepaths))

        # do not want to have all files open at the same time
        #srcs = [open(f, "rU") for f in support_filepaths]

        # store filepaths, to open individually in loop
        srcs = support_filepaths

    for sidx, src in enumerate(srcs):

        # hack needed because we do not want to open all input files at the
        # same time; if not a file object, assume it is a file path and create
        # corresponding file object
        if not isinstance(src, file):
            src = open(src, "rU")

        name = getattr(src, "name", "<stdin>")
        messenger.send_info("Processing %d of %d: '%s'" % (sidx+1, len(srcs), name), wrap=False)
        for tidx, tree in enumerate(tree_source_iter(src,
                schema=schema,
                taxon_set=taxon_set,
                store_tree_weights=weighted_trees,
                as_rooted=is_rooted)):
            if tidx >= tree_offset:
                if (log_frequency == 1) or (tidx > 0 and log_frequency > 0 and tidx % log_frequency == 0):
                    messenger.send_info("(processing) '%s': tree at offset %d" % (name, tidx), wrap=False)
                treesplit.encode_splits(tree)
                split_distribution.count_splits_on_tree(tree)
                topology_counter.count(tree, tree_splits_encoded=True)
            else:
                if (log_frequency == 1) or (tidx > 0 and log_frequency > 0 and tidx % log_frequency == 0):
                    messenger.send_info("(processing) '%s': tree at offset %d (skipping)" % (name, tidx), wrap=False)
        try:
            src.close()
        except ValueError:
            # "I/O operation on closed file" if we try to close sys.stdin
            pass

    messenger.send_info("Serial processing of %d source(s) completed." % len(srcs))
    return split_distribution, topology_counter

def main_cli():

    description =  "%s %s %s" % (_program_name, _program_version, _program_subtitle)
    usage = "%prog [options] TREES-FILE [TREES-FILE [TREES-FILE [...]]"

    parser = OptionParser(usage=usage, add_help_option=True, version = _program_version, description=description)

    sum_tree_optgroup = OptionGroup(parser, "Source Treatment Options")
    parser.add_option_group(sum_tree_optgroup)
    sum_tree_optgroup.add_option("-b", "--burnin",
            action="store",
            dest="burnin",
            type="int",
            default=0,
            help='number of trees to skip from the beginning of *each tree file* when counting support [default=%default]')

    source_tree_optgroup = OptionGroup(parser, "Source Tree Options")
    parser.add_option_group(source_tree_optgroup)
    source_tree_optgroup.add_option("--rooted",
            action="store_true",
            dest="rooted_trees",
            default=None,
            help="treat trees as rooted")
    source_tree_optgroup.add_option("--unrooted",
            action="store_false",
            dest="rooted_trees",
            default=None,
            help="treat trees as unrooted")
    source_tree_optgroup.add_option("--weighted-trees",
            action="store_true",
            dest="weighted_trees",
            default=False,
            help="use weights of trees as indicated by '[&W m/n]' comment to weight contribution of splits found on each tree to overall split frequencies")
    source_tree_optgroup.add_option("--from-newick-stream",
            action="store_true",
            dest="from_newick_stream",
            default=False,
            help="support trees will be streamed in newick format")
    source_tree_optgroup.add_option("--from-nexus-stream",
            action="store_true",
            dest="from_nexus_stream",
            default=False,
            help="support trees will be streamed in NEXUS format")

    support_summarization_optgroup = OptionGroup(parser, "Support Summarization Options")
    parser.add_option_group(support_summarization_optgroup)
    support_summarization_optgroup.add_option("-l","--support-as-labels",
            action="store_true",
            dest="support_as_labels",
            default=True,
            help="indicate branch support as internal node labels [default=%default]")
    support_summarization_optgroup.add_option("-v","--support-as-lengths",
            action="store_false",
            dest="support_as_labels",
            default=True,
            help="indicate branch support as branch lengths (otherwise support will be indicated by internal node labels)")
    support_summarization_optgroup.add_option("-p", "--percentages",
            action="store_true",
            dest="support_as_percentages",
            default=False,
            help="indicate branch support as percentages (otherwise, will report as proportions by default)")
    support_summarization_optgroup.add_option("-d", "--decimals",
            dest="support_label_decimals",
            type="int",
            metavar="#",
            default=2,
            help="number of decimal places in indication of support values [default=%default]")

    edge_summarization_optgroup = OptionGroup(parser, "Edge Summarization Options")
    parser.add_option_group(edge_summarization_optgroup)
    edge_summarization_optgroup.add_option("--with-node-ages",
            action="store_true",
            dest="calc_node_ages",
            default=False,
            help="summarize node ages as well as edge lengths (implies rooted trees and requires all trees to be ultrametric)")
    edge_summarization_choices = ["mean-length", "median-length", "mean-age", "median-age", "keep"]
    edge_summarization_optgroup.add_option("-e", "--edges",
            type="choice",
            dest="edge_summarization",
            metavar="<%s>" % ("|".join(edge_summarization_choices)),
            choices=edge_summarization_choices,
            default=None,
            help="""\
one of: %s; set edge lengths of target tree(s) to mean/median lengths/ages of
corresponding splits or edges of input trees (note that using 'mean-age' or
'median-age' require rooted ultrametric input trees")""" % (str(edge_summarization_choices)))
    edge_summarization_optgroup.add_option("--collapse-negative-edges",
            action="store_true",
            dest="collapse_negative_edges",
            default=False,
            help="(if setting edge lengths) force parent node ages to be at least as old as its oldest child when summarizing node ages")
    edge_summarization_optgroup.add_option("--no-extended-summary",
            action="store_true",
            dest="suppress_extended_summary",
            default=False,
            help="do not calculate ranges, 5%/95 quartiles, 95% HPD's etc. of edge lengths and node ages")

    target_tree_optgroup = OptionGroup(parser, 'Target Tree Options')
    parser.add_option_group(target_tree_optgroup)
    target_tree_optgroup.add_option("-t","--target",
            dest="target_tree_filepath",
            default=None,
            help="path to optional target, model or best topology tree file (Newick or NEXUS format) "
            + "to which support will be mapped; "
            + "if not given, then a majority-rule clade consensus tree will be constructed based on the "
            + "all the trees given in the support tree files (except for those discarded as burn-ins), "
            + "and this will be used as the target tree")
    target_tree_optgroup.add_option("-f", "--min-clade-freq",
            dest="min_clade_freq",
            type="float",
            default=0.50,
            metavar="#.##",
            help="minimum frequency or probability for a clade or a split to be "\
                    + "included in the consensus tree, if used [default=%default]")
    target_tree_optgroup.add_option("--no-branch-lengths",
            action="store_true",
            dest="no_branch_lengths",
            default=False,
            help="by default, if using a consensus tree as the target tree, branch lengths will be the mean of the lengths " \
                    + "of the given branch across all trees considered; this option forces branch " \
                    + "lengths to be unspecified (obviously, this is only applicable if you do not ask the support to be mapped as "  \
                    + "branch lengths)")

    output_filepath_optgroup = OptionGroup(parser, "Output File Options")
    parser.add_option_group(output_filepath_optgroup)
    output_filepath_optgroup.add_option("-o","--output",
            dest="output_filepath",
            default=None,
            help="path to output file (if not given, will print to standard output)")
    output_filepath_optgroup.add_option("--no-taxa-block",
            action="store_false",
            dest="include_taxa_block",
            default=True,
            help="do not include a taxa block in the output treefile (otherwise will create taxa block by default)")
    output_filepath_optgroup.add_option("--no-meta-comments",
            action="store_false",
            dest="include_meta_comments",
            default=True,
            help="do not include initial file comment annotating details of scoring operation")
    output_filepath_optgroup.add_option("-c", "--additional-comments",
            action="store",
            dest="additional_comments",
            default=None,
            help="additional comments to be added to the summary file")
    output_filepath_optgroup.add_option("--to-newick",
            action="store_true",
            dest="to_newick_format",
            default=False,
            help="save results in NEWICK (PHYLIP) format (default is to save in NEXUS format)")
    output_filepath_optgroup.add_option("--to-phylip",
            action="store_true",
            dest="to_newick_format",
            default=False,
            help="same as --newick")
    output_filepath_optgroup.add_option("-r", "--replace",
            action="store_true",
            dest="replace",
            default=False,
            help="replace/overwrite output file without asking if it already exists ")

    other_optgroup = OptionGroup(parser, "Other Options")
    parser.add_option_group(other_optgroup)

    other_optgroup.add_option("--trprobs", "--calc-tree-probabilities",
            dest="trprobs_filepath",
            default=None,
            metavar="FILEPATH",
            help="if specified, a file listing tree (topologies) and the " \
                    + "frequencies of their occurrences will be saved to FILEPATH")
    other_optgroup.add_option("--extract-edges",
            dest="split_edges_filepath",
            default=None,
            metavar="FILEPATH",
            help="if specified, a tab-delimited file of splits and their edge " \
                    + "lengths across input trees will be saved to FILEPATH")

    run_optgroup = OptionGroup(parser, "Program Run Options")
    parser.add_option_group(run_optgroup)
    if _MP:
        run_optgroup.add_option("-m", "--multiprocessing",
                action="store",
                dest="multiprocess",
                metavar="NUM-PROCESSES",
                default=None,
                help="run in parallel mode with up to a maximum of NUM-PROCESSES processes " \
                        + "(specify '*' to run in as many processes as there are cores on the "\
                        + "local machine)")

    run_optgroup.add_option("-g", "--log-frequency",
            type="int",
            metavar="LOG-FREQUENCY",
            dest="log_frequency",
            default=500,
            help="tree processing progress logging frequency (default=%default; set to 0 to suppress)")
    run_optgroup.add_option("-q", "--quiet",
            action="store_true",
            dest="quiet",
            default=False,
            help="suppress ALL logging, progress and feedback messages")
    run_optgroup.add_option("--ignore-missing-support",
            action="store_true",
            dest="ignore_missing_support",
            default=False,
            help="ignore missing support tree files (at least one must exist!)")
    run_optgroup.add_option("--ignore-missing-target",
            action="store_true",
            dest="ignore_missing_target",
            default=False,
            help="ignore missing target tree file (will construct majority rule consensus tree if missing)")

    (opts, args) = parser.parse_args()
    if opts.quiet:
        messaging_level = ConsoleMessenger.ERROR_MESSAGING_LEVEL
    else:
        messaging_level = ConsoleMessenger.INFO_MESSAGING_LEVEL
    messenger = ConsoleMessenger(name="SumTrees", messaging_level=messaging_level)

    # splash
    if not opts.quiet:
        show_splash(prog_name=_program_name,
                prog_subtitle=_program_subtitle,
                prog_version=_program_version,
                prog_author=_program_author,
                prog_copyright=_program_copyright,
                dest=sys.stderr,
                extended=False)

    ###################################################
    # Support file idiot checking

    support_filepaths = []
    if len(args) > 0:
        for fpath in args:
            fpath = os.path.expanduser(os.path.expandvars(fpath))
            if not os.path.exists(fpath):
                if opts.ignore_missing_support:
                    messenger.send_warning("Support file not found: '%s'" % fpath)
                else:
                    messenger.send_error("Terminating due to missing support files. "
                           + "Use the '--ignore-missing-support' option to continue even "
                           + "if some files are missing.")
                    sys.exit(1)
            else:
                support_filepaths.append(fpath)
        if len(support_filepaths) == 0:
            messenger.send_error("No valid sources of input trees specified. "
                    + "Please provide the path to at least one (valid and existing) file "
                    + "containing tree samples to summarize.")
            sys.exit(1)
    else:
        if not opts.from_newick_stream and not opts.from_nexus_stream:
            messenger.send_info("No sources of input trees specified. "
                    + "Please provide the path to at least one (valid and existing) file "
                    + "containing tree samples to summarize. See '--help' for other options.")
            sys.exit(1)

    ###################################################
    # Lots of other idiot-checking ...

    # target tree
    if opts.target_tree_filepath is not None:
        target_tree_filepath = os.path.expanduser(os.path.expandvars(opts.target_tree_filepath))
        if not os.path.exists(target_tree_filepath):
            if opts.ignore_missing_target:
                if not opts.quiet:
                    messenger.send_warning("Target tree file not found: '%s': using majority-rule consensus tree instead." % target_tree_filepath)
                target_tree_filepath = None
            else:
                messenger.send_error("Target tree file not found: '%s'" % target_tree_filepath)
                sys.exit(1)
    else:
        target_tree_filepath = None

    ### TODO: these will be command-line options in the future
    ### here we just set it
    assert not hasattr(opts, 'outgroup')
    opts.outgroup = None
    assert not hasattr(opts, 'root_target')
    opts.root_target = None

    ### TODO: idiot-check edge length summarization
    # edge lengths
    if opts.edge_summarization:
        opts.edge_summarization = opts.edge_summarization.lower()
        if opts.edge_summarization not in edge_summarization_choices:
            messenger.send_error("'%s' is not a valid edge summarization choice; must be one of: %s" % (opts.edge_summarization, edge_summarization_choices))
            sys.exit(1)
    if opts.edge_summarization == "mean-age" or opts.edge_summarization == "median-age":
        opts.calc_node_ages = True
        opts.rooted_trees = True
    else:
        if opts.calc_node_ages:
            opts.rooted_trees = True

    # output
    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if confirm_overwrite(filepath=output_fpath, replace_without_asking=opts.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)

    if opts.trprobs_filepath:
        trprobs_filepath = os.path.expanduser(os.path.expandvars(opts.trprobs_filepath))
        if confirm_overwrite(filepath=trprobs_filepath, replace_without_asking=opts.replace):
            trprobs_dest = open(trprobs_filepath, "w")
        else:
            sys.exit(1)
        opts.calc_tree_probs = True
    else:
        trprobs_dest = None
        opts.calc_tree_probs = False

    if opts.split_edges_filepath:
        split_edges_filepath = os.path.expanduser(os.path.expandvars(opts.split_edges_filepath))
        if confirm_overwrite(filepath=split_edges_filepath, replace_without_asking=opts.replace):
            split_edges_dest = open(split_edges_filepath, "w")
        else:
            sys.exit(1)
    else:
        split_edges_dest = None

    if opts.from_newick_stream:
        schema = "newick"
    elif opts.from_nexus_stream:
        schema = "nexus"
    else:
        schema = 'nexus/newick'

    ###################################################
    # Main work begins here: Count the splits

    start_time = datetime.datetime.now()
    master_split_distribution = None
    if (support_filepaths is not None and len(support_filepaths) > 1) \
            and _MP \
            and opts.multiprocess:
        if opts.multiprocess is not None:
            if opts.multiprocess == "*":
                num_processes = multiprocessing.cpu_count()
            elif  opts.multiprocess == "@":
                num_processes = len(support_filepaths)
            else:
                try:
                    num_processes = int(opts.multiprocess)
                except ValueError:
                    messenger.send_error("'%s' is not a valid number of processes (must be a positive integer)." % opts.multiprocess)
                    sys.exit(1)
            if num_processes <= 0:
                messenger.send_error("Maximum number of processes set to %d: cannot run SumTrees with less than 1 process" % num_processes)
                sys.exit(1)
            if num_processes == 1:
                messenger.send_warning("Running in parallel processing mode but limited to only 1 process: probably more efficient to run in serial mode!")

        master_split_distribution, master_topology_counter = process_sources_parallel(
                num_processes=num_processes,
                support_filepaths=support_filepaths,
                schema=schema,
                is_rooted=opts.rooted_trees,
                ignore_node_ages=not opts.calc_node_ages,
                calc_tree_probs=opts.calc_tree_probs,
                weighted_trees=opts.weighted_trees,
                tree_offset=opts.burnin,
                log_frequency=opts.log_frequency,
                messenger=messenger)
    else:
        if (_MP and opts.multiprocess is not None and len(support_filepaths) == 1):
            messenger.send_warning("Parallel processing mode requested but only one source specified: defaulting to serial mode.")
        if opts.from_newick_stream or opts.from_nexus_stream:
            support_filepaths = None
        master_split_distribution, master_topology_counter = process_sources_serial(
                support_filepaths=support_filepaths,
                schema=schema,
                is_rooted=opts.rooted_trees,
                ignore_node_ages=not opts.calc_node_ages,
                calc_tree_probs=opts.calc_tree_probs,
                weighted_trees=opts.weighted_trees,
                tree_offset=opts.burnin,
                log_frequency=opts.log_frequency,
                messenger=messenger)

    ###################################################
    # Compose post-counting report

    # if not splits counted or the taxon set was not populated for any reason,
    # we just produce an empty block so we don't crash as we report nothing of interest
    if master_split_distribution.taxon_set is None:
        assert(master_split_distribution.total_trees_counted == 0)
        master_split_distribution.taxon_set = dendropy.TaxonSet()

    # taxon set to handle target trees
    master_taxon_set = master_split_distribution.taxon_set

    report = []
    report.append("%d trees considered in total for split support assessment." % (master_split_distribution.total_trees_counted))
    if opts.rooted_trees is None:
        report.append("Tree rooting as given by tree statement (defaults to unrooted).")
    elif opts.rooted_trees:
        report.append("Trees treated as rooted.")
    else:
        report.append("Trees treated as unrooted.")
    if opts.weighted_trees:
        report.append("Trees treated as weighted (default weight = 1.0).")
    else:
        report.append("Trees treated as unweighted.")
    n_taxa = len(master_taxon_set)
    report.append("%d unique taxa across all trees." % n_taxa)
    num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits = master_split_distribution.splits_considered()
    report.append("%d unique splits out of %d total splits counted." % (num_unique_splits, num_splits))
    report.append("%d unique non-trivial splits out of %d total non-trivial splits counted." % (num_nt_unique_splits, num_nt_splits))

    comments = []
    comments.extend(report)
    messenger.send_info("Split counting completed:")
    messenger.send_info_lines(report, prefix=" - ")

    ###################################################
    #  Target tree and mapping

    if not opts.support_as_percentages and opts.support_label_decimals < 2:
        messenger.send_warning("Reporting support by proportions require that support will be reported to at least 2 decimal places")
        opts.support_label_decimals = 2

    tsum = treesum.TreeSummarizer()
    tsum.support_as_labels = opts.support_as_labels
    tsum.support_as_percentages = opts.support_as_percentages
    tsum.support_label_decimals = opts.support_label_decimals
    tsum.weighted_splits = opts.weighted_trees

    if opts.support_as_percentages:
        support_units = "Percentage"
    else:
        support_units = "Proportion (frequency or probability)"
    if opts.support_as_labels:
        support_show = "node labels"
    else:
        support_show = "branch lengths"
    support_summarization = "%s of support for each split indicated by %s" % (support_units, support_show)

    tt_trees = []
    if target_tree_filepath is not None:
        messenger.send_info("Mapping support to target tree ...")
        for tree in tree_source_iter(stream=open(target_tree_filepath, 'r'),
                schema="nexus/newick",
                taxon_set=master_taxon_set,
                as_rooted=opts.rooted_trees):
            if opts.root_target:
                if opts.outgroup:
                    pass
                else:
                    tree.root_at_midpoint(splits=True)
            if opts.rooted_trees and not tree.is_rooted:
                messenger.send_error("Support trees are treated as rooted, but target tree is unrooted. Root target tree(s) and re-run, or run using the '--root-target' flag.")
                sys.exit(1)
            stree = tsum.map_split_support_to_tree(tree,
                    master_split_distribution)
            tt_trees.append(stree)
        messenger.send_info("Parsed '%s': %d tree(s) in file" % (target_tree_filepath, len(tt_trees)))
        comments.append("Split support mapped to trees in:")
        comments.append("  - '%s' (%d trees)" % (os.path.abspath(target_tree_filepath), len(tt_trees)))
        if opts.root_target:
            if opts.outgroup:
                comments.append("Target tree(s) rooted using outgroup: %s." % opts.outgroup)
            else:
                comments.append("Target tree(s) rooted at midpoint.")
        comments.append(support_summarization + '.')
    else:
        messenger.send_info("Constructing clade consensus tree ...")
        if opts.min_clade_freq > 1.0:
            messenger.send_warning("Maximum frequency threshold for clade inclusion is 1.0: reset to 1.0.")
            min_freq = 1.0
        else:
            min_freq = opts.min_clade_freq
        stree = tsum.tree_from_splits(master_split_distribution,
                min_freq=min_freq,
                include_edge_lengths=not opts.no_branch_lengths)
        if opts.root_target:
            stree.reroot_at_midpoint(splits=True)
        report = []
        report.append("Consensus tree (%f clade frequency threshold) constructed from splits." % min_freq)
        if not opts.edge_summarization:
            if opts.calc_node_ages:
                tsum.summarize_node_ages_on_tree(tree=stree,
                        split_distribution=master_split_distribution,
                        set_edge_lengths=True,
                        collapse_negative_edges=opts.collapse_negative_edges,
                        allow_negative_edges=not opts.collapse_negative_edges,
                        summarization_func=statistics.median)
                report.append("Consensus tree node ages set to median ages of corresponding nodes of input trees.")
            else:
                report.append("Consensus tree edge lengths set to mean of corresponding edges of input trees.")
        tt_trees.append(stree)
        if opts.root_target:
            if opts.outgroup:
                report.append("Consensus tree rooted using outgroup: %s." % opts.outgroup)
            else:
                report.append("Consensus tree rooted at midpoint.")
        report.append(support_summarization + ".")
        messenger.send_info_lines(report)
        comments.extend(report)

    if not opts.suppress_extended_summary:
        messenger.send_info("Summarizing node ages and lengths ...")
        for stree in tt_trees:
            tsum.annotate_nodes_and_edges(tree=stree, split_distribution=master_split_distribution)

    if opts.edge_summarization is not None:
        if opts.edge_summarization.startswith('mean'):
            summary_func_desc = "mean"
            summarization_func = lambda x: statistics.mean_and_sample_variance(x)[0]
        else:
            summary_func_desc = "desc"
            summarization_func = statistics.median
        if opts.edge_summarization.endswith("age"):
            messenger.send_info("Mapping node ages ...")
            comments.append("Setting node ages of output tree(s) to %s ages of corresponding nodes of input trees." % summary_func_desc)
            if opts.collapse_negative_edges:
                comments.append("Parent node ages coerced to be at least as old as oldest daughter node age.")
                collapse_negative_edges = True
                allow_negative_edges = False
            else:
                comments.append("Parent node ages not adjusted: negative edge lengths allowed.")
                collapse_negative_edges = False
                allow_negative_edges = True
            for stree in tt_trees:
                tsum.summarize_node_ages_on_tree(tree=stree,
                        split_distribution=master_split_distribution,
                        set_edge_lengths=True,
                        collapse_negative_edges=collapse_negative_edges,
                        allow_negative_edges=allow_negative_edges,
                        summarization_func=summarization_func)
        elif opts.edge_summarization.endswith("length"):
            messenger.send_info("Mapping edge lengths ...")
            comments.append("Setting edge lengths of output tree(s) to %s length of corresponding edges of input trees." % summary_func_desc)
            for stree in tt_trees:
                tsum.summarize_edge_lengths_on_tree(tree=stree,
                        split_distribution=master_split_distribution,
                        summarization_func=summarization_func)

    end_time = datetime.datetime.now()

    ###################################################
    #  RESULTS

    messenger.send_info("Writing results ...")

    final_run_report = []
    final_run_report.append("Began at: %s." % (start_time.isoformat(' ')))
    final_run_report.append("Ended at: %s." % (end_time.isoformat(' ')))
    hours, mins, secs = str(end_time-start_time).split(":")
    run_time = "Run time: %s hour(s), %s minute(s), %s second(s)." % (hours, mins, secs)
    final_run_report.append(run_time)

    output_dataset = dendropy.DataSet(dendropy.TreeList(tt_trees, taxon_set=master_taxon_set))
    if opts.to_newick_format:
        output_dataset.write(output_dest, "newick")
    else:
        if opts.include_taxa_block:
            simple = False
        else:
            simple = True
        if opts.include_meta_comments:
            comment = []
            try:
                username = getpass.getuser()
            except:
                username = "a user"
            comment.append("%s %s by %s." % (_program_name, _program_version, _program_author))
            comment.append("Using DendroPy Version %s by Jeet Sukumaran and Mark T. Holder."
                % dendropy.__version__)
            python_version = sys.version.replace("\n", "").replace("[", "(").replace("]",")")
            comment.append("Running under Python %s on %s." % (python_version, sys.platform))
            comment.append("Executed on %s by %s@%s." % (platform.node(),  username, socket.gethostname()))
            if support_filepaths is not None and len(support_filepaths) > 0:
                comment.append("Basis of split support:")
                for support_file in support_filepaths:
                    comment.append("  - '%s'" % os.path.abspath(support_file))
            else:
                comment.append("Basis of split support: trees read from standard input.")
            comment.extend(final_run_report)
            comment.extend(comments)
        if opts.additional_comments:
            comment.append("\n")
            comment.append(opts.additional_comments)
        output_dataset.write(output_dest, "nexus", simple=simple, comment=comment)

    messenger.send_info("Writing tree probabilities ...")
    if trprobs_dest:
        tree_list = dendropy.TreeList(taxon_set=master_split_distribution.taxon_set)
        tree_freqs = master_topology_counter.calc_freqs()
        cumulative_prob = 0.0
        for idx, (tree, freq) in enumerate(tree_freqs.items()):
            # in parallel mode, the topology is stored as a normalized newick string;
            # in serial mode, we store the tree directly
            # here, we deal with either possibility
            ### TODO: check to see if there are any performance gains in using
            ### the normalized newick string approach; we are currently forced to use
            ### the normalized newick string tree store function in parallel mode
            ### due to the fact that Tree objects cannot be pickled. If this is actually
            ### advantageous in other ways (e.g. memory usage), and the negatives (e.g.,
            ### loss of speed) are not too bad, we should use this for serial mode as well
            if isinstance(tree, str):
                tree_list.read_from_string(tree, 'newick')
                tree = tree_list[-1]
            else:
                tree_list.append(tree)
            count = freq[0]
            prop = freq[1]
            cumulative_prob += prop
            tree.probability = prop
            tree.count = count
            tree.cumulative_probability = cumulative_prob
            tree.annotate('count')
            tree.annotate('probability')
            tree.annotate('cumulative_probability')
            tree.label = "Tree%d" % (idx+1)
        tree_list.write_to_stream(trprobs_dest,
                'nexus',
                edge_lengths=False,
                write_rooting=False,
                intenral_labels=False,
                annotations_as_comments=True,
                annotations_as_nhx=False,
                write_item_comments=False
                )

    messenger.send_info("Writing split edge lengths ...")
    if split_edges_dest:
        for split in master_split_distribution.splits:
            row = []
            row.append(master_split_distribution.taxon_set.split_as_newick_string(split))
            for edge_length in master_split_distribution.split_edge_lengths[split]:
                row.append("%s" % edge_length)
            split_edges_dest.write("%s\n" % ("\t".join(row)))

    if not opts.output_filepath:
        pass
    else:
        messenger.send_info("Results written to: '%s'." % (output_fpath))

    ###################################################
    #  WRAP UP
    messenger.send_info("Summarization completed.")
    messenger.send_info_lines(final_run_report)
    messenger.silent = True

if __name__ == '__main__':
    main_cli()
