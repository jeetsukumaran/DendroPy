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
Marge labels of splits/branches from different input trees onto corresponding
splits/branches of a single tree.
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
from dendropy.mathlib import statistics

_program_name = "SumLabels"
_program_subtitle = "Phylogenetic Tree Label Concatenation"
_program_date = "Mar 20 2012"
_program_version = "Version 1.0.0 (%s)" % _program_date
_program_author = "Jeet Sukumaran"
_program_contact = "jeetsukumaran@gmail.com"
_program_copyright = "Copyright (C) 2012 Jeet Sukumaran.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \
                 "There is NO WARRANTY,\nto the extent permitted by law."

def main_cli():

    description =  "%s %s %s" % (_program_name, _program_version, _program_subtitle)
    usage = "%prog [options] -t <TARGET-TREE-FILE> <TREES-FILE> [TREES-FILE [TREES-FILE [...]]"

    parser = OptionParser(usage=usage, add_help_option=True, version = _program_version, description=description)

    parser.add_option("-t","--target",
            dest="target_tree_filepath",
            default=None,
            help="path to file with tree (Newick or NEXUS format) "
            + "to which labels will be written")
    parser.add_option("--preserve-target-labels",
            action="store_true",
            dest="preserve_target_labels",
            default=False,
            help="keep any existing labels on target tree (by default, these will be cleared before writing the new labels)")
    parser.add_option("--rooted",
            action="store_true",
            dest="rooted_trees",
            default=None,
            help="treat trees as rooted")
    parser.add_option("--unrooted",
            action="store_false",
            dest="rooted_trees",
            default=None,
            help="treat trees as unrooted")
    parser.add_option("--ignore-missing-source",
            action="store_true",
            dest="ignore_missing_source",
            default=False,
            help="ignore missing source tree files (at least one must exist!)")
    parser.add_option("-o","--output",
            dest="output_filepath",
            default=None,
            help="path to output file (if not given, will print to standard output)")
    parser.add_option("-s","--separator",
            dest="separator",
            default="/",
            help="string to use to separate labels from different source trees (default='%default')")
    parser.add_option("--no-taxa-block",
            action="store_false",
            dest="include_taxa_block",
            default=True,
            help="do not include a taxa block in the output treefile (otherwise will create taxa block by default)")
    parser.add_option("-c", "--additional-comments",
            action="store",
            dest="additional_comments",
            default=None,
            help="additional comments to be added to the summary file")
    parser.add_option("--to-newick",
            action="store_true",
            dest="to_newick_format",
            default=False,
            help="save results in NEWICK (PHYLIP) format (default is to save in NEXUS format)")
    parser.add_option("--to-phylip",
            action="store_true",
            dest="to_newick_format",
            default=False,
            help="same as --newick")
    parser.add_option("-r", "--replace",
            action="store_true",
            dest="replace",
            default=False,
            help="replace/overwrite output file without asking if it already exists ")
    parser.add_option("-q", "--quiet",
            action="store_true",
            dest="quiet",
            default=False,
            help="suppress ALL logging, progress and feedback messages")

    (opts, args) = parser.parse_args()
    if opts.quiet:
        messaging_level = ConsoleMessenger.ERROR_MESSAGING_LEVEL
    else:
        messaging_level = ConsoleMessenger.INFO_MESSAGING_LEVEL
    messenger = ConsoleMessenger(name="SumLabels", messaging_level=messaging_level)

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
    # Source file idiot checking

    source_filepaths = []
    if len(args) > 0:
        for fpath in args:
            fpath = os.path.expanduser(os.path.expandvars(fpath))
            if not os.path.exists(fpath):
                if opts.ignore_missing_source:
                    messenger.send_warning("Source file not found: '%s'" % fpath)
                else:
                    messenger.send_error("Terminating due to missing source files. "
                           + "Use the '--ignore-missing-source' option to continue even "
                           + "if some files are missing.")
                    sys.exit(1)
            else:
                source_filepaths.append(fpath)
        if len(source_filepaths) == 0:
            messenger.send_error("No valid sources of input trees specified. "
                    + "Please provide the path to at least one (valid and existing) file "
                    + "containing trees")
            sys.exit(1)
    else:
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
            messenger.send_error("Target tree file not found: '%s'" % target_tree_filepath)
            sys.exit(1)
    else:
        messenger.send_error("Target tree file not specified: use the '-t' or '--target' option to provide path to target tree")
        sys.exit(1)

    # output
    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if confirm_overwrite(filepath=output_fpath, replace_without_asking=opts.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)

    # taxon set to handle target trees
    master_taxon_set = dendropy.TaxonSet()
    is_rooted = opts.rooted_trees
    messenger.send_info("Reading target tree: '%s'" % target_tree_filepath)
    target_tree = None
    for tree in tree_source_iter(
            open(target_tree_filepath, "rU"),
            schema='nexus/newick',
            taxon_set=master_taxon_set,
            as_rooted=is_rooted):
        target_tree = tree
        break
    split_labels = {}
    for src_fpath in source_filepaths:
        messenger.send_info("Reading source tree(s) from: '%s'" % src_fpath)
        for tree in tree_source_iter(
                open(src_fpath, "rU"),
                schema='nexus/newick',
                taxon_set=master_taxon_set,
                as_rooted=is_rooted):
            tree.update_splits()
            for split, edge in tree.split_edges.items():
                label = edge.head_node.label
                print label
                if not label:
                    continue
                try:
                    split_labels[split].append(label)
                except KeyError:
                    split_labels[split] = [label]
    messenger.send_info("Mapping labels")
    target_tree.update_splits()
    for split, edge in target_tree.split_edges.items():
        label = []
        if opts.preserve_target_labels and edge.head_node.label:
            label.append(edge.head_node.label)
        elif not opts.preserve_target_labels:
            edge.head_node.label = None
        if split in split_labels:
            label.extend(split_labels[split])
        else:
            pass
            # messenger.send_warning("Split on target tree not found in source trees: ignoring")
        if label:
            edge.head_node.label = opts.separator.join(label)
    output_dataset = dendropy.DataSet(dendropy.TreeList([target_tree], taxon_set=master_taxon_set))
    if opts.to_newick_format:
        output_dataset.write(output_dest,
                "newick",
                suppress_rooting=False,
                suppress_edge_lengths=False,
                unquoted_underscores=False,
                preserve_spaces=False,
                store_tree_weights=False,
                suppress_annotations=False,
                annotations_as_nhx=False,
                suppress_item_comments=False,
                suppress_leaf_taxon_labels=False,
                suppress_leaf_node_labels=True,
                suppress_internal_taxon_labels=False,
                suppress_internal_node_labels=False,
                node_label_element_separator=' ',
                node_label_compose_func=None)
    else:
        if opts.include_taxa_block:
            simple = False
        else:
            simple = True
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
        if opts.additional_comments:
            comment.append("\n")
            comment.append(opts.additional_comments)
        output_dataset.write(output_dest,
                "nexus",
                simple=simple,
                file_comments=comment,
                suppress_rooting=False,
                unquoted_underscores=False,
                preserve_spaces=False,
                store_tree_weights=False,
                suppress_annotations=False,
                annotations_as_nhx=False,
                suppress_item_comments=False,
                suppress_leaf_taxon_labels=False,
                suppress_leaf_node_labels=True,
                suppress_internal_taxon_labels=False,
                suppress_internal_node_labels=False,
                node_label_element_separator=' ',
                node_label_compose_func=None)
    if not opts.output_filepath:
        pass
    else:
        messenger.send_info("Results written to: '%s'." % (output_fpath))

if __name__ == '__main__':
    main_cli()
