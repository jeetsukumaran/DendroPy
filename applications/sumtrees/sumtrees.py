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

if not (sys.version_info.major >= 3 and sys.version_info.minor >= 4):
    from dendropy.utility.filesys import pre_py34_open as open
try:
    # Python 3
    import queue
except ImportError:
    # Python 2.7
    import Queue as queue
import multiprocessing

from dendropy.utility import cli
from dendropy.utility import constants

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

def main():
    parser = argparse.ArgumentParser(
            description=__doc__,
            # formatter_class=argparse.HelpFormatter,
            formatter_class=cli.CustomFormatter,
            )
    source_options = parser.add_argument_group("Source Options")
    source_options.add_argument("-i","--source-format",
            metavar="FORMAT",
            default="nexus/newick",
            choices=["nexus", "newick", "phylip", "nexml"],
            help="Format of the source trees (defaults to handling either NEXUS or NEWICK through inspection; it is more efficient to explicitly specify the format if it is known).")
    source_options.add_argument("-b", "--burnin",
            type=int,
            default=0,
            help="Number of trees to skip from the beginning of *each* tree file when counting support (default: %(default)s).")
    source_options.add_argument("--rooted",
            action="store_true",
            default=None,
            help="Treat source trees as rooted.")
    source_options.add_argument("--unrooted",
            action="store_false",
            default=None,
            help="Treat source trees as unrooted.")
    source_options.add_argument("-u", "--ultrametric",
            action="store_true",
            default=False,
            help="Assume source trees are ultrametric (implies '--rooted'; will result in node ages being summarized; will result in error if trees are not ultrametric).")
    source_options.add_argument("-y", "--ultrametricity-precision",
            action="store_true",
            default=constants.DEFAULT_ULTRAMETRICITY_CHECK_PRECISION,
            help="Precision to use when validating ultrametricity (default: %(default)s; specify '0' to disable validation).")
    source_options.add_argument("--weighted-trees",
            action="store_true",
            default=False,
            help="Use weights of trees (as indicated by '[&W m/n]' comment token) to weight contribution of splits found on each tree to overall split frequencies.")

    summary_tree_options = parser.add_argument_group("Summary Target Tree Options")
    summary_tree_options.add_argument(
            "-s", "--summary-tree-type",
            default=None,
            metavar="{consensus,mct,msct}",
            help="\n".join((
                "R}The summary tree to construct or select as a target ",
                "for mapping support and other information from the ",
                "source trees:",
                "- 'consensus' : consensus tree;",
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
    summary_tree_options.add_argument("-t", "--target-trees-filepath",
            default=None,
            metavar="FILE",
            help=(
                 "Instead of constructing a summary tree as a target, use the "
                 "tree or trees in the file specified by FILE as target(s) "
                 "to which to map support and other information summarized "
                 "by trees in the source set."
                 ))
    summary_tree_options.add_argument("-f", "--min-consensus-freq",
            type=float,
            default=constants.GREATER_THAN_HALF,
            metavar="#.##",
            help=(
                "If using a consensus tree summarization strategy, then "
                "this is the minimum frequency or probability for a clade "
                "or a split to be included in the resulting tree "
                "(default: > 0.05)."
                ))
    summary_tree_options.add_argument("--root-target-at-midpoint",
            action="store_true",
            default=None,
            help="Root target tree(s) at midpoint.")
    summary_tree_options.add_argument("--root-target-outgroup",
            metavar="TAXON-LABEL",
            default=None,
            help="Root target tree(s) using specified taxon as outgroup.")

    support_summarization_optgroup = parser.add_argument_group("Support Summarization Options")
    support_summarization_optgroup.add_argument("-l","--support-as-labels",
            action="store_const",
            dest="support_annotation_target",
            default=1,
            const=1,
            help="Indicate split/clade/branch support as internal node labels.")
    support_summarization_optgroup.add_argument("-v","--support-as-lengths",
            action="store_const",
            dest="support_annotation_target",
            default=1,
            const=2,
            help="Indicate split/clade/branch support as branch lengths.")
    support_summarization_optgroup.add_argument("-x","--no-support",
            action="store_const",
            dest="support_annotation_target",
            default=1,
            const=0,
            help=("Do not indicate support with internal node labels or edge"
                  " lengths. Support will still be indicated as node metadata "
                  " annotations unless '--no-summary-metadata' is specified "
                 ))
    support_summarization_optgroup.add_argument("-p", "--percentages",
            action="store_true",
            dest="support_as_percentages",
            default=False,
            help="Indicate branch support as percentages (otherwise, will report as proportions by default).")
    support_summarization_optgroup.add_argument("-d", "--decimals",
            dest="support_label_decimals",
            type=int,
            metavar="#",
            default=8,
            help="Number of decimal places in indication of support values (default: %(default)s).")

    edge_summarization_optgroup = parser.add_argument_group("Edge Length Summarization Options")
    edge_summarization_choices = ["mean-length", "median-length", "mean-age", "median-age", "keep", "clear"]
    edge_summarization_optgroup.add_argument("-e", "--edges",
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
                "                     sourced from an external file using the '-t'",
                "                     or '--target' option",
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
    edge_summarization_optgroup.add_argument("--collapse-negative-edges",
            action="store_true",
            default=False,
            help="(If setting edge lengths) force parent node ages to be at least as old as its oldest child when summarizing node ages.")

    other_summarization_optgroup = parser.add_argument_group("Other Summarization Options")
    #other_summarization_optgroup.add_argument("--with-node-ages",
    #        action="store_true",
    #        dest="calc_node_ages",
    #        default=None,
    #        help="summarize node ages as well as edge lengths (implies '--rooted' and '--ultrametric'; automatically enabled if '--ultrametric' is specified; will result in error if trees are not ultrametric)")
    other_summarization_optgroup.add_argument("--trprobs", "--calc-tree-probabilities",
            dest="trprobs_filepath",
            default=None,
            metavar="FILEPATH",
            help=("If specified, a file listing tree (topologies) and the "
                 "frequencies of their occurrences will be saved to FILEPATH."))
    other_summarization_optgroup.add_argument("--extract-edges",
            dest="split_edge_map_filepath",
            default=None,
            metavar="FILEPATH",
            help=("if specified, a tab-delimited file of splits and their edge "
                  "lengths across input trees will be saved to FILEPATH"))
    other_summarization_optgroup.add_argument("--no-node-ages",
            action="store_false",
            default=None,
            help="Do not calculate/summarize node ages, even if '--ultrametric' is specified")
    other_summarization_optgroup.add_argument("--no-summary-metadata",
            action="store_true",
            default=False,
            help="Do not annotate nodes and edges with any summarization information.")

    output_options = parser.add_argument_group("Output File Options")
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
            help="replace/overwrite output file without asking if it already exists ")

    run_optgroup = parser.add_argument_group("Program Run Options")
    run_optgroup.add_argument("-m", "--multiprocessing",
            dest="multiprocess",
            metavar="NUM-PROCESSES",
            default=None,
            help="Run in parallel mode with up to a maximum of NUM-PROCESSES processes " \
                 "(specify '*' to run in as many processes as there are cores on the "\
                 "local machine).")
    run_optgroup.add_argument("-g", "--log-frequency",
            type=int,
            metavar="LOG-FREQUENCY",
            default=500,
            help="Tree processing progress logging frequency (default: %(default)s; set to 0 to suppress).")
    run_optgroup.add_argument("-q", "--quiet",
            action="store_true",
            default=False,
            help="Suppress ALL logging, progress and feedback messages.")
    run_optgroup.add_argument("--ignore-missing-support",
            action="store_true",
            default=False,
            help="Ignore missing support tree files (at least one must exist!).")

    args = parser.parse_args()

    if not args.quiet:
        cli.show_splash(prog_name=_program_name,
                prog_subtitle=_program_subtitle,
                prog_version=_program_version,
                prog_author=_program_author,
                prog_copyright=_program_copyright,
                dest=sys.stderr,
                include_citation=True,
                include_copyright=False,
                additional_citations=[_program_citation],
                )

if __name__ == '__main__':
    main()
