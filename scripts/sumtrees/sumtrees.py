#! /usr/bin/env python

############################################################################
##  sumtrees.py
##
##  Copyright 2008 Jeet Sukumaran.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

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

_program_name = "SumTrees"
_program_subtitle = "Phylogenetic Tree Split Support Summarization"
_program_date = "Aug 22 2010"
_program_version = "Version 3.0.0 (%s)" % _program_date
_program_author = "Jeet Sukumaran and Mark T. Holder"
_program_contact = "jeetsukumaran@gmail.com"
_program_copyright = "Copyright (C) 2008 Jeet Sukumaran.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \
                 "There is NO WARRANTY,\nto the extent permitted by law."

class SplitCountingThread(multiprocessing.Process):

    def pickle_result(split_distribution):
        result = (split_distribution.total_trees_counted,
                split_distribution.splits,
                split_distribution.split_counts,
                split_distribution.split_edge_lengths,
                split_distribution.split_node_ages)
        return result
    pickle_result = staticmethod(pickle_result)

    def unpickle_result(pickled_result, split_distribution=None):
        if split_distribution is None:
            split_distribution = treesplit.SplitDistribution()
        split_distribution.total_trees_counted += pickled_result[0]
        for split in pickled_result[1]:
            if split not in split_distribution.splits:
                split_distribution.splits.append(split)
                if split in pickled_result[2]:
                    split_distribution.split_counts[split] = 0
                if split in pickled_result[3]:
                    split_distribution.split_edge_lengths[split] = []
                if split in pickled_result[4]:
                   split_distribution.split_node_ages[split] = []
        for s, c in pickled_result[2].items():
            split_distribution.split_counts[s] += c
        for s, c in pickled_result[3].items():
            split_distribution.split_edge_lengths[s].extend(c)
        for s, c in pickled_result[4].items():
            split_distribution.node_ages[s].extend(c)
        return split_distribution
    unpickle_result = staticmethod(unpickle_result)

    def __init__(self,
            work_queue,
            result_queue,
            schema,
            taxon_labels,
            is_rooted,
            tree_offset,
            thread_idx,
            messenger,
            messenger_lock,
            log_frequency=1000):
        multiprocessing.Process.__init__(self)
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.schema = schema
        self.taxon_labels = list(taxon_labels)
        self.taxon_set = dendropy.TaxonSet(self.taxon_labels)
        self.split_distribution = treesplit.SplitDistribution(taxon_set=self.taxon_set)
        self.split_distribution.is_rooted = is_rooted
        self.tree_offset = tree_offset
        self.thread_idx = thread_idx
        self.messenger = messenger
        self.messenger_lock = messenger_lock
        self.log_frequency = log_frequency
        self.kill_received = False

    def send_message(self, msg, level, wrap=True):
        if self.messenger is None:
            return
        if self.messenger.messaging_level > level or self.messenger.silent:
            return
        msg = "Thread %d: %s" % (self.thread_idx+1, msg)
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
            for tidx, tree in enumerate(tree_source_iter(fsrc, schema=self.schema, taxon_set=self.taxon_set)):
                if tidx >= self.tree_offset:
                    if (self.log_frequency == 1) or (tidx > 0 and self.log_frequency > 0 and tidx % self.log_frequency == 0):
                        self.send_info("(processing) '%s': tree at offset %d" % (source, tidx), wrap=False)
                    treesplit.encode_splits(tree)
                    self.split_distribution.count_splits_on_tree(tree)
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
            result = SplitCountingThread.pickle_result(self.split_distribution)
            self.result_queue.put(result)

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
        num_threads,
        support_filepaths,
        schema,
        is_rooted,
        tree_offset,
        log_frequency,
        messenger):
    """
    Returns a SplitDistribution object summarizing all trees found in
    `support_filepaths`.
    """

    # describe
    messenger.send_info("Running in multithreaded mode (up to %d threads)." % num_threads)
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

    # launch threads
    messenger.send_info("Launching worker threads ...")
    result_queue = multiprocessing.Queue()
    messenger_lock = multiprocessing.Lock()
    for idx in range(num_threads):
        sct = SplitCountingThread(work_queue,
                result_queue,
                schema=schema,
                taxon_labels=taxon_labels,
                is_rooted=is_rooted,
                tree_offset=tree_offset,
                thread_idx=idx,
                messenger=messenger,
                messenger_lock=messenger_lock,
                log_frequency=log_frequency)
        sct.start()

    # collate results
    thread_result_count = 0
    split_distribution = treesplit.SplitDistribution(taxon_set=taxon_set)
    split_distribution.is_rooted = is_rooted
    while thread_result_count < num_threads:
        result = result_queue.get()
        SplitCountingThread.unpickle_result(result, split_distribution)
        thread_result_count += 1
    messenger.send_info("Recovered results from all worker threads.")
    return split_distribution

def process_sources_serial(
        support_filepaths,
        schema,
        is_rooted,
        tree_offset,
        log_frequency,
        messenger):
    """
    Returns a SplitDistribution object summarizing all trees found in
    `support_filepaths`.
    """
    messenger.send_info("Running in single-threaded mode.")
    taxon_set = dendropy.TaxonSet()
    split_distribution = treesplit.SplitDistribution(taxon_set=taxon_set)
    if support_filepaths is None or len(support_filepaths) == 0:
        messenger.send_info("Reading trees from standard input.")
        srcs = [sys.stdin]
    else:
        messenger.send_info("%d source(s) to be processed." % len(support_filepaths))
        srcs = [open(f, "rU") for f in support_filepaths]
    for sidx, src in enumerate(srcs):
        name = getattr(src, "name", "<stdin>")
        messenger.send_info("Processing %d of %d: '%s'" % (sidx+1, len(srcs), name), wrap=False)
        for tidx, tree in enumerate(tree_source_iter(src, schema=schema, taxon_set=taxon_set, is_rooted=is_rooted)):
            if tidx >= tree_offset:
                if (log_frequency == 1) or (tidx > 0 and log_frequency > 0 and tidx % log_frequency == 0):
                    messenger.send_info("(processing) '%s': tree at offset %d" % (name, tidx), wrap=False)
                treesplit.encode_splits(tree)
                split_distribution.count_splits_on_tree(tree)
            else:
                if (log_frequency == 1) or (tidx > 0 and log_frequency > 0 and tidx % log_frequency == 0):
                    messenger.send_info("(processing) '%s': tree at offset %d (skipping)" % (name, tidx), wrap=False)
    messenger.send_info("Serial processing of %d source(s) completed." % len(srcs))
    return split_distribution

def main_cli():

    description =  "%s %s %s" % (_program_name, _program_version, _program_subtitle)
    usage = "%prog [options] TREES-FILE [TREES-FILE [TREES-FILE [...]]"

    parser = OptionParser(usage=usage, add_help_option=True, version = _program_version, description=description)

    sum_tree_optgroup = OptionGroup(parser, "Summarization Options")
    parser.add_option_group(sum_tree_optgroup)
    sum_tree_optgroup.add_option("-b", "--burnin",
            action="store",
            dest="burnin",
            type="int",
            default=0,
            help='number of trees to skip from the beginning of *each tree file* when counting support [default=%default]')

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

    source_tree_optgroup = OptionGroup(parser, "Source Tree Options")
    parser.add_option_group(source_tree_optgroup)
    source_tree_optgroup.add_option("--rooted",
            action="store_true",
            dest="rooted_trees",
            default=False,
            help="treat trees as rooted")
    source_tree_optgroup.add_option("--unrooted",
            action="store_false",
            dest="rooted_trees",
            default=False,
            help="treat trees as unrooted")
    source_tree_optgroup.add_option("--from-newick-stream",
            action="store_true",
            dest="from_newick_stream",
            default=False,
            help="support trees will be streamed in Newick format")
    source_tree_optgroup.add_option("--from-nexus-stream",
            action="store_true",
            dest="from_nexus_stream",
            default=False,
            help="support trees will be streamed in NEXUS format")

    output_tree_optgroup = OptionGroup(parser, "Output Tree Options")
    parser.add_option_group(output_tree_optgroup)
    output_tree_optgroup.add_option("-l","--support-as-labels",
            action="store_true",
            dest="support_as_labels",
            default=True,
            help="indicate branch support as internal node labels [default=%default]")
    output_tree_optgroup.add_option("-v","--support-as-lengths",
            action="store_false",
            dest="support_as_labels",
            default=True,
            help="indicate branch support as branch lengths (otherwise support will be indicated by internal node labels)")
    output_tree_optgroup.add_option("-p", "--percentages",
            action="store_true",
            dest="support_as_percentages",
            default=False,
            help="indicate branch support as percentages (otherwise, will report as proportions by default)")
    output_tree_optgroup.add_option("-d", "--decimals",
            dest="support_label_decimals",
            type="int",
            metavar="#",
            default=2,
            help="number of decimal places in indication of support values [default=%default]")

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

    other_optgroup.add_option("-e","--split-edges",
            dest="split_edges_filepath",
            default=None,
            metavar="FILEPATH",
            help="if specified, a tab-delimited file of splits and their edge " \
                    + "lengths across runs will be saved to FILEPATH")

    run_optgroup = OptionGroup(parser, "Program Run Options")
    parser.add_option_group(run_optgroup)
    if _MP:
        run_optgroup.add_option("-m", "--multithreaded",
                action="store",
                dest="multithreaded",
                metavar="NUM-THREADS",
                default=None,
                help="run in multithreaded (parallel) mode: process tree sources in separate " \
                        + "threads, with up to a maximum of NUM-THREADS parallel threads " \
                        + "(specify '*' to run in as many threads as there are cores on the "\
                        + "local machine)")
        #run_optgroup.add_option("-x", "--max-threads",
        #        dest="max_threads",
        #        metavar="MAX-THREADS",
        #        type="int",
        #        default=None,
        #        help="limit number of threads launched to MAX-THREADS (implies '-m'/'--multithreaded')")
    run_optgroup.add_option("-q", "--quiet",
            action="store_true",
            dest="quiet",
            default=False,
            help="suppress progress messages")
    run_optgroup.add_option("-g", "--log-frequency",
            type="int",
            metavar="LOG-FREQUENCY",
            dest="log_frequency",
            default=500,
            help="tree processing progress logging frequency (default=%default; set to 0 to suppress)")
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

    # output
    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if confirm_overwrite(filepath=output_fpath, replace_without_asking=opts.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)

    if opts.split_edges_filepath:
        split_edges_filepath = os.path.expanduser(os.path.expandvars(opts.split_edges_filepath))
        if confirm_overwrite(filepath=output_fpath, replace_without_asking=opts.replace):
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
            and opts.multithreaded:
        if opts.multithreaded is not None:
            if opts.multithreaded == "*":
                num_threads = multiprocessing.cpu_count()
            elif  opts.multithreaded == "@":
                num_threads = len(support_filepaths)
            else:
                try:
                    num_threads = int(opts.multithreaded)
                except ValueError:
                    messenger.send_error("'%s' is not a valid number of threads (must be a positive integer)." % opts.multithreaded)
                    sys.exit(1)
            if num_threads <= 0:
                messenger.send_error("Maximum number of threads set to %d: cannot run SumTrees with less than 1 thread" % num_threads)
                sys.exit(1)
            if num_threads == 1:
                messenger.send_warning("Running in multithreaded mode but limited to only 1 thread: probably more efficient to run in serial mode!")

        master_split_distribution = process_sources_parallel(
                num_threads=num_threads,
                support_filepaths=support_filepaths,
                schema=schema,
                is_rooted=opts.rooted_trees,
                tree_offset=opts.burnin,
                log_frequency=opts.log_frequency,
                messenger=messenger)
    else:
        if (_MP and opts.multithreaded is not None and len(support_filepaths) == 1):
            messenger.send_warning("Multithreaded mode requested but only one source specified: defaulting to single-threaded mode.")
        if opts.from_newick_stream or opts.from_nexus_stream:
            support_filepaths = None
        master_split_distribution = process_sources_serial(
                support_filepaths=support_filepaths,
                schema=schema,
                is_rooted=opts.rooted_trees,
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
    if opts.rooted_trees:
        report.append("Trees treated as rooted.")
    else:
        report.append("Trees treated as unrooted.")
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
    tsum.ignore_node_ages = True # until a more efficient implementation is developed

    if opts.support_as_percentages:
        support_units = "Percentage"
    else:
        support_units = "Proportion (frequency or probability)"
    if opts.support_as_labels:
        support_show = "node labels"
    else:
        support_show = "branch lengths"
    support_indication = "%s of support for each split indicated by %s" % (support_units, support_show)

    tt_trees = []
    if target_tree_filepath is not None:
        messenger.send_info("Mapping support to target tree ...")
        for tree in tree_source_iter(stream=open(target_tree_filepath, 'r'),
                schema="nexus/newick",
                taxon_set=master_taxon_set):
            tt_trees.append(tsum.map_split_support_to_tree(tree, master_split_distribution))
        messenger.send_info("Parsed '%s': %d tree(s) in file" % (target_tree_filepath, len(tt_trees)))
        comments.append("Split support mapped to trees in:")
        comments.append("  - '%s' (%d trees)" % (os.path.abspath(target_tree_filepath), len(tt_trees)))
        comments.append(support_indication + '.')
    else:
        messenger.send_info("Constructing clade consensus tree ...")
        if opts.min_clade_freq > 1.0:
            messenger.send_warning("Maximum frequency threshold for clade inclusion is 1.0: reset to 1.0.")
            min_freq = 1.0
        else:
            min_freq = opts.min_clade_freq
        tt_trees.append(tsum.tree_from_splits(master_split_distribution,
                                              min_freq=min_freq,
                                              include_edge_lengths=not opts.no_branch_lengths))
        report = []
        report.append("Consensus tree (%f clade frequency threshold) constructed from splits." % min_freq)
        report.append(support_indication + ".")
        messenger.send_info_lines(report)
        comments.extend(report)

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

    if split_edges_dest:
        for split in master_split_distribution.splits:
            row = []
            row.append(master_taxa_block.split_as_newick_string(split))
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
    try:
        main_cli()
    except (KeyboardInterrupt, EOFError), e:
        sys.stderr.write("SumTrees: Terminating (user-abort).\n")
        sys.exit(1)
    except Exception, e:
        sys.stderr.write("SumTrees: Error encountered: %s : %s.\n" % (str(type(e)), str(e)))
        raise # reraise exception, with correct traceback
