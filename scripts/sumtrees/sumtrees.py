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

import dendropy
from dendropy import treesplit
from dendropy import treesum
from dendropy.dataio import tree_source_iter
from dendropy.dataio import multi_tree_source_iter
from dendropy.dataio import newick
from dendropy.utility.messaging import ConsoleMessenger
from dendropy.utility.cli import confirm_overwrite, show_splash

_program_name = 'SumTrees'
_program_subtitle = 'Phylogenetic Tree Split Support Summarization'
_program_date = 'Feb 1 2010'
_program_version = 'Version 2.0.2 (%s)' % _program_date
_program_author = 'Jeet Sukumaran and Mark T. Holder'
_program_contact = 'jeetsukumaran@gmail.com'
_program_copyright = "Copyright (C) 2008 Jeet Sukumaran.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \
                 "There is NO WARRANTY,\nto the extent permitted by law."

def main_cli():

    description =  '%s %s %s' % (_program_name, _program_version, _program_subtitle)
    usage = "%prog [options] <TREES FILE> [<TREES FILE> [<TREES FILE> [...]]"

    parser = OptionParser(usage=usage, add_help_option=True, version = _program_version, description=description)

    sum_tree_optgroup = OptionGroup(parser, 'Summarization Options')
    parser.add_option_group(sum_tree_optgroup)
    sum_tree_optgroup.add_option('-b', '--burnin',
                        action='store',
                        dest='burnin',
                        type='int', # also 'float', 'string' etc.
                        default=0,
                        help='number of trees to skip from the beginning of *each tree file* when counting support [default=%default]')

    target_tree_optgroup = OptionGroup(parser, 'Target Tree Options')
    parser.add_option_group(target_tree_optgroup)
    target_tree_optgroup.add_option('-t','--target',
                  dest='target_tree_filepath',
                  default=None,
                  help="path to optional target, model or best topology tree file (Newick or NEXUS format) "
                       + "to which support will be mapped; "
                       + "if not given, then a majority-rule clade consensus tree will be constructed based on the "
                       + "all the trees given in the support tree files (except for those discarded as burn-ins), "
                       + "and this will be used as the target tree")
    target_tree_optgroup.add_option('-f', '--min-clade-freq',
                      dest='min_clade_freq',
                      type='float',
                      default=0.50,
                      metavar='#.##',
                      help="minimum frequency or probability for a clade or a split to be included in the consensus tree, if used [default=%default]")
    target_tree_optgroup.add_option('--no-branch-lengths',
                      action='store_true',
                      dest='no_branch_lengths',
                      default=False,
                      help="by default, if using a consensus tree as the target tree, branch lengths will be the mean of the lengths " \
                          + "of the given branch across all trees considered; this option forces branch " \
                          + "lengths to be unspecified (obviously, this is only applicable if you do not ask the support to be mapped as "  \
                          + "branch lengths)")

    source_tree_optgroup = OptionGroup(parser, 'Source Tree Options')
    parser.add_option_group(source_tree_optgroup)
    source_tree_optgroup.add_option('--rooted',
                      action='store_true',
                      dest='rooted_trees',
                      default=False,
                      help="treat trees as rooted")
    source_tree_optgroup.add_option('--unrooted',
                      action='store_false',
                      dest='rooted_trees',
                      default=False,
                      help="treat trees as unrooted")
    source_tree_optgroup.add_option('--from-newick-stream',
                      action='store_true',
                      dest='from_newick_stream',
                      default=False,
                      help="support trees will be streamed in Newick format")
    source_tree_optgroup.add_option('--from-nexus-stream',
                      action='store_true',
                      dest='from_nexus_stream',
                      default=False,
                      help="support trees will be streamed in NEXUS format")

    output_tree_optgroup = OptionGroup(parser, 'Output Tree Options')
    parser.add_option_group(output_tree_optgroup)
    output_tree_optgroup.add_option('-l','--support-as-labels',
                      action='store_true',
                      dest='support_as_labels',
                      default=True,
                      help="indicate branch support as internal node labels [default=%default]")
    output_tree_optgroup.add_option('-v','--support-as-lengths',
                      action='store_false',
                      dest='support_as_labels',
                      default=True,
                      help="indicate branch support as branch lengths (otherwise support will be indicated by internal node labels)")
    output_tree_optgroup.add_option('-p', '--percentages',
                      action='store_true',
                      dest='support_as_percentages',
                      default=False,
                      help="indicate branch support as percentages (otherwise, will report as proportions by default)")
    output_tree_optgroup.add_option('-d', '--decimals',
                      dest='support_label_decimals',
                      type='int',
                      metavar='#',
                      default=2,
                      help="number of decimal places in indication of support values [default=%default]")

    output_filepath_optgroup = OptionGroup(parser, 'Output File Options')
    parser.add_option_group(output_filepath_optgroup)
    output_filepath_optgroup.add_option('-o','--output',
                  dest='output_filepath',
                  default=None,
                  help="path to output file (if not given, will print to standard output)")
    output_filepath_optgroup.add_option('--no-taxa-block',
                      action='store_false',
                      dest='include_taxa_block',
                      default=True,
                      help="do not include a taxa block in the output treefile (otherwise will create taxa block by default)")
    output_filepath_optgroup.add_option('--no-meta-comments',
                      action='store_false',
                      dest='include_meta_comments',
                      default=True,
                      help="include initial file comment annotating details of scoring operation")
    output_filepath_optgroup.add_option('-m', '--additional_comments',
                      action='store',
                      dest='additional_comments',
                      default=None,
                      help="additional comments to be added to the summary file")
    output_filepath_optgroup.add_option('--to-newick',
                      action='store_true',
                      dest='to_newick_format',
                      default=False,
                      help="save results in NEWICK (PHYLIP) format (default is to save in NEXUS format)")
    output_filepath_optgroup.add_option('--to-phylip',
                      action='store_true',
                      dest='to_newick_format',
                      default=False,
                      help="same as --newick")
    output_filepath_optgroup.add_option('-r', '--replace',
                      action='store_true',
                      dest='replace',
                      default=False,
                      help="replace/overwrite output file without asking if it already exists ")

    other_optgroup = OptionGroup(parser, 'Other Options')
    parser.add_option_group(other_optgroup)

    other_optgroup.add_option('-e','--split-edges',
                  dest='split_edges_filepath',
                  default=None,
                  metavar='FILEPATH',
                  help="if specified, a tab-delimited file of splits and their edge " \
                    + "lengths across runs will be saved to FILEPATH")

    run_optgroup = OptionGroup(parser, 'Program Run Options')
    parser.add_option_group(run_optgroup)
    run_optgroup.add_option('-q', '--quiet',
                      action='store_true',
                      dest='quiet',
                      default=False,
                      help="suppress progress messages")
    run_optgroup.add_option('--ignore-missing-support',
                      action='store_true',
                      dest='ignore_missing_support',
                      default=False,
                      help="ignore missing support tree files (at least one must exist!)")
    run_optgroup.add_option('--ignore-missing-target',
                      action='store_true',
                      dest='ignore_missing_target',
                      default=False,
                      help="ignore missing target tree file (will construct majority rule consensus tree if missing)")

    (opts, args) = parser.parse_args()
    messenger = ConsoleMessenger(quiet=opts.quiet)

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
    if len(args) == 0 and (opts.from_newick_stream or opts.from_nexus_stream):
        if not opts.quiet:
            sys.stderr.write("(reading trees from standard input)")
        support_file_objs = [sys.stdin]
    else:
        missing = False
        for fpath in args:
            fpath = os.path.expanduser(os.path.expandvars(fpath))
            if not os.path.exists(fpath):
                messenger.send_error('Support file not found: "%s"' % fpath)
                missing = True
            else:
                support_filepaths.append(fpath)
        if missing:
            messenger.send("")
            if opts.ignore_missing_support:
                pass
            else:
                messenger.send_formatted('Terminating due to missing support files. '
                       + 'Use the "--ignore-missing-support" option to continue even '
                       + 'if some files are missing.', force=True)
                sys.exit(1)
        if len(support_filepaths) == 0:
            messenger.send_formatted("No sources of support specified or could be found. "
            + "Please provide the path to at least one (valid and existing) file "
            + "containing tree samples "
            + "to summarize.", force=True)
            sys.exit(1)

        support_file_objs = [open(f, "r") for f in support_filepaths]

    ###################################################
    # Lots of other idiot-checking ...

    # target tree
    if opts.target_tree_filepath is not None:
        target_tree_filepath = os.path.expanduser(os.path.expandvars(opts.target_tree_filepath))
        if not os.path.exists(target_tree_filepath):
            messenger.send_error('Target tree file not found: "%s"\n' % target_tree_filepath)
            if opts.ignore_missing_target:
                if not opts.quiet:
                    messenger.send('Will construct and use majority-rule consensus tree instead.\n')
                target_tree_filepath = None
            else:
                sys.exit(1)
    else:
        target_tree_filepath = None

    # output
    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if confirm_overwrite(output_fpath, messenger, opts.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)

    if opts.split_edges_filepath:
        split_edges_filepath = os.path.expanduser(os.path.expandvars(opts.split_edges_filepath))
        if confirm_overwrite(split_edges_filepath, messenger, opts.replace):
            split_edges_dest = open(split_edges_filepath, "w")
        else:
            sys.exit(1)
    else:
        split_edges_dest = None


    ###################################################
    # Main work begins here: Count the splits

    start_time = datetime.datetime.now()

    comments = []
    tsum = treesum.TreeSummarizer()
    tsum.support_as_labels = opts.support_as_labels
    tsum.support_as_percentages = opts.support_as_percentages
    if not opts.support_as_percentages and opts.support_label_decimals < 2:
        messenger.send_error("(WARNING: proportions require that support will be reported to at least 2 decimal places)")
        opts.support_label_decimals = 2
    tsum.support_label_decimals = opts.support_label_decimals
    tsum.ignore_node_ages = True # until a more efficient implementation is developed

    messenger.send("### COUNTING SPLITS ###\n")
    if opts.from_newick_stream:
        file_format = "newick"
    elif opts.from_nexus_stream:
        file_format = "nexus"
    else:
        file_format = 'nexus/newick'

    taxon_set = dendropy.TaxonSet()
    if file_format == "nexus":
        encode_splits = True
    else:
        encode_splits = False # cannot encode while parsing newick
    tree_source = multi_tree_source_iter(sources=support_file_objs,
                                         schema=file_format,
                                         tree_offset=opts.burnin,
                                         write_progress=messenger.write,
                                         taxon_set=taxon_set,
                                         encode_splits=encode_splits)
    split_distribution = treesplit.SplitDistribution()
    split_distribution.is_rooted = opts.rooted_trees
    tsum.count_splits_on_trees(tree_source,
        split_distribution=split_distribution,
        trees_splits_encoded=encode_splits)
    if split_distribution.taxon_set is None:
        assert(tsum.total_trees_counted == 0)
        split_distribution.taxon_set = dendropy.TaxonSet() # we just produce an empty block so we don't crash as we report nothing of interest
    report = []
#    report.append("%d trees read from %d files." % (tree_source.total_trees_read, len(support_filepaths)))
#    report.append("%d trees from each file requested to be ignored for burn-in." % (opts.burnin))
#    report.append("%d trees ignored in total." % (tree_source.total_trees_ignored))
    report.append("%d trees considered in total for split support assessment." % (tsum.total_trees_counted))
    if opts.rooted_trees:
        report.append("Trees treated as rooted.")
    else:
        report.append("Trees treated as unrooted.")
    n_taxa = len(split_distribution.taxon_set)
    report.append("%d unique taxa across all trees." % n_taxa)
    num_splits, num_unique_splits, num_nt_splits, num_nt_unique_splits = split_distribution.splits_considered()
    report.append("%d unique splits out of %d total splits counted." % (num_unique_splits, num_splits))
    report.append("%d unique non-trivial splits out of %d total non-trivial splits counted." % (num_nt_unique_splits, num_nt_splits))

    comments.extend(report)
    messenger.send("---")
    messenger.send_multi(report)
    messenger.send("")

    ###################################################
    #  Target tree and mapping

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
        messenger.send("### MAPPING SUPPORT TO TARGET TREE(S) ###\n")
        for tree in tree_source_iter(stream=open(target_tree_filepath, 'r'), schema="nexus/newick", taxon_set=taxon_set):
            tt_trees.append(tsum.map_split_support_to_tree(tree, split_distribution))
        messenger.send('Parsed "%s": %d tree(s) in file' % (target_tree_filepath, len(tt_trees)))
        comments.append('Split support mapped to trees in:')
        comments.append('  - "%s" (%d trees)' % (os.path.abspath(target_tree_filepath), len(tt_trees)))
        comments.append(support_indication + ".")
    else:
        messenger.send("### CONSTRUCTING CLADE CONSENSUS TREE ###\n")
        if opts.min_clade_freq > 1.0:
            messenger.send("Maximum frequency threshold for clade inclusion is 1.0: reset to 1.0.", force=True)
            min_freq = 1.0
        else:
            min_freq = opts.min_clade_freq
        tt_trees.append(tsum.tree_from_splits(split_distribution,
                                              min_freq=min_freq,
                                              include_edge_lengths=not opts.no_branch_lengths))
        report = []
        report.append('Consensus tree (%f clade frequency threshold) constructed from splits.' % min_freq)
        report.append(support_indication + ".")
        messenger.send_multi(report)
        comments.extend(report)
    messenger.send("")

    end_time = datetime.datetime.now()

    ###################################################
    #  RESULTS

    messenger.send("### RESULTS ###\n")

    final_run_report = []
    final_run_report.append("Began at: %s." % (start_time.isoformat(' ')))
    final_run_report.append("Ended at: %s." % (end_time.isoformat(' ')))
    hours, mins, secs = str(end_time-start_time).split(":")
    run_time = "Run time: %s hour(s), %s minute(s), %s second(s)." % (hours, mins, secs)
    final_run_report.append(run_time)

    output_dataset = dendropy.DataSet(dendropy.TreeList(tt_trees, taxon_set=taxon_set))

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
            comment.append("Basis of split support:")
            for support_file in support_filepaths:
                comment.append('  - "%s"' % os.path.abspath(support_file))
            comment.extend(final_run_report)
            comment.extend(comments)
        if opts.additional_comments:
            comment.append("\n")
            comment.append(opts.additional_comments)
        output_dataset.write(output_dest, "nexus", simple=simple, comment=comment)

    if split_edges_dest:
        for split in split_distribution.splits:
            row = []
            row.append(split_distribution.taxa_block.split_as_newick_string(split))
            for edge_length in split_distribution.split_edge_lengths[split]:
                row.append("%s" % edge_length)
            split_edges_dest.write("%s\n" % ("\t".join(row)))

    if not opts.output_filepath:
        pass
    else:
        messenger.send('Results written to: "%s".' % (output_fpath))
    messenger.send("")

    ###################################################
    #  WRAP UP
    messenger.send("### DONE ###\n")
    messenger.send_multi(final_run_report)

if __name__ == '__main__':
    try:
        main_cli()
    except (KeyboardInterrupt, EOFError), e:
        sys.stderr.write("Terminating (user-abort).\n")
        sys.exit(1)
    except Exception, e:
        sys.stderr.write("Error encountered: %s : %s.\n" % (str(type(e)), str(e)))
        raise # reraise exception, with correct traceback
