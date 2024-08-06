#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
Merge labels of splits/branches from different input trees onto corresponding
splits/branches of a single tree.
"""

import os
import sys
import argparse

import socket
try:
    import getpass
except:
    pass
import platform

import dendropy
from dendropy.utility.messaging import ConsoleMessenger
from dendropy.utility.cli import confirm_overwrite, show_splash
from dendropy.utility import deprecate

_program_name = "SumLabels"
_program_subtitle = "Phylogenetic Tree Label Concatenation"
_program_date = "Jan 20 2017"
_program_version = "Version 2.0.0 (%s)" % _program_date
_program_author = "Jeet Sukumaran"
_program_contact = "jeetsukumaran@gmail.com"
_program_copyright = "Copyright (C) 2017 Jeet Sukumaran.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \
                 "There is NO WARRANTY,\nto the extent permitted by law."

def main_cli():

    if "bin/sumlabels.py" in sys.argv[0]:
        deprecate.dendropy_deprecation_warning(
            message="`sumlabels.py` entrypoint is deprecated since DendroPy 5. "
            "Use `sumlabels` entrypoint (i.e., without `.py`) instead.",
        )

    description =  "%s %s %s" % (_program_name, _program_version, _program_subtitle)

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
            "sources",
            metavar="TREEFILE",
            nargs="+")
    parser.add_argument("-t","--target",
            dest="target_tree_filepath",
            default=None,
            help="path to file with tree (Newick or NEXUS format) "
            + "to which labels will be written")
    parser.add_argument("--preserve-target-labels",
            action="store_true",
            dest="preserve_target_labels",
            default=False,
            help="keep any existing labels on target tree (by default, these will be cleared before writing the new labels)")
    parser.add_argument("--rooted",
            action="store_true",
            dest="rooted_trees",
            default=None,
            help="treat trees as rooted")
    parser.add_argument("--unrooted",
            action="store_false",
            dest="rooted_trees",
            default=None,
            help="treat trees as unrooted")
    parser.add_argument("--ignore-missing-source",
            action="store_true",
            dest="ignore_missing_source",
            default=False,
            help="ignore missing source tree files (at least one must exist!)")
    parser.add_argument("-o","--output",
            dest="output_filepath",
            default=None,
            help="path to output file (if not given, will print to standard output)")
    parser.add_argument("-s","--separator",
            dest="separator",
            default="/",
            help="string to use to separate labels from different source trees (default='%(default)s')")
    parser.add_argument("--no-taxa-block",
            action="store_false",
            dest="include_taxa_block",
            default=True,
            help="do not include a taxa block in the output treefile (otherwise will create taxa block by default)")
    parser.add_argument("-c", "--additional-comments",
            action="store",
            dest="additional_comments",
            default=None,
            help="additional comments to be added to the summary file")
    parser.add_argument("--to-newick",
            action="store_true",
            dest="to_newick_format",
            default=False,
            help="save results in NEWICK (PHYLIP) format (default is to save in NEXUS format)")
    parser.add_argument("--to-phylip",
            action="store_true",
            dest="to_newick_format",
            default=False,
            help="same as --to-newick")
    parser.add_argument("-r", "--replace",
            action="store_true",
            dest="replace",
            default=False,
            help="replace/overwrite output file without asking if it already exists ")
    parser.add_argument("-q", "--quiet",
            action="store_true",
            dest="quiet",
            default=False,
            help="suppress ALL logging, progress and feedback messages")

    args = parser.parse_args()
    if args.quiet:
        messaging_level = ConsoleMessenger.ERROR_MESSAGING_LEVEL
    else:
        messaging_level = ConsoleMessenger.INFO_MESSAGING_LEVEL
    messenger = ConsoleMessenger(name="SumLabels", messaging_level=messaging_level)

    # splash
    if not args.quiet:
        show_splash(prog_name=_program_name,
                prog_subtitle=_program_subtitle,
                prog_version=_program_version,
                prog_author=_program_author,
                prog_copyright=_program_copyright,
                dest=sys.stderr,
                )

    ###################################################
    # Source file idiot checking

    source_filepaths = []
    if len(args.sources) > 0:
        for fpath in args.sources:
            fpath = os.path.expanduser(os.path.expandvars(fpath))
            if not os.path.exists(fpath):
                if args.ignore_missing_source:
                    messenger.send_warning("Source file not found: '%s'" % fpath)
                else:
                    messenger.error("Terminating due to missing source files. "
                           + "Use the '--ignore-missing-source' option to continue even "
                           + "if some files are missing.")
                    sys.exit(1)
            else:
                source_filepaths.append(fpath)
        if len(source_filepaths) == 0:
            messenger.error("No valid sources of input trees specified. "
                    + "Please provide the path to at least one (valid and existing) file "
                    + "containing trees")
            sys.exit(1)
    else:
        messenger.info("No sources of input trees specified. "
                + "Please provide the path to at least one (valid and existing) file "
                + "containing tree samples to summarize. See '--help' for other options.")
        sys.exit(1)

    ###################################################
    # Lots of other idiot-checking ...

    # target tree
    if args.target_tree_filepath is not None:
        target_tree_filepath = os.path.expanduser(os.path.expandvars(args.target_tree_filepath))
        if not os.path.exists(target_tree_filepath):
            messenger.error("Target tree file not found: '%s'" % target_tree_filepath)
            sys.exit(1)
    else:
        messenger.error("Target tree file not specified: use the '-t' or '--target' option to provide path to target tree")
        sys.exit(1)

    # output
    if args.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(args.output_filepath))
        if confirm_overwrite(filepath=output_fpath, replace_without_asking=args.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)

    # taxon set to handle target trees
    master_taxon_namespace = dendropy.TaxonNamespace()
    is_rooted = args.rooted_trees
    messenger.info("Reading target tree: '%s'" % target_tree_filepath)
    target_tree = None
    if is_rooted:
        rooting = "force-rooted"
    else:
        rooting = None
    for tree in dendropy.Tree.yield_from_files(
            [target_tree_filepath, "rU",],
            schema='nexus/newick',
            taxon_namespace=master_taxon_namespace,
            rooting=rooting):
        target_tree = tree
        break
    bipartition_labels = {}
    for src_fpath in source_filepaths:
        messenger.info("Reading source tree(s) from: '%s'" % src_fpath)
        for tree in dendropy.Tree.yield_from_files(
                [src_fpath,],
                schema='nexus/newick',
                taxon_namespace=master_taxon_namespace,
                rooting=rooting):
            tree.encode_bipartitions()
            for bipartition, edge in tree.bipartition_edge_map.items():
                label = edge.head_node.label
                if not label:
                    continue
                try:
                    bipartition_labels[bipartition].append(label)
                except KeyError:
                    bipartition_labels[bipartition] = [label]
    messenger.info("Mapping labels")
    target_tree.encode_bipartitions()
    for bipartition, edge in target_tree.bipartition_edge_map.items():
        label = []
        if args.preserve_target_labels and edge.head_node.label:
            label.append(edge.head_node.label)
        elif not args.preserve_target_labels:
            edge.head_node.label = None
        if bipartition in bipartition_labels:
            label.extend(bipartition_labels[bipartition])
        else:
            pass
            # messenger.send_warning("Split on target tree not found in source trees: ignoring")
        if label:
            edge.head_node.label = args.separator.join(label)
    output_dataset = dendropy.DataSet()
    tree_list = output_dataset.new_tree_list(taxon_namespace=master_taxon_namespace)
    tree_list.append(target_tree)
    if args.to_newick_format:
        output_dataset.write(
                file=output_dest,
                schema="newick",
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
                )
    else:
        if args.include_taxa_block:
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
        if args.additional_comments:
            comment.append("\n")
            comment.append(args.additional_comments)
        output_dataset.write(
                file=output_dest,
                schema="nexus",
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
                )
    if not args.output_filepath:
        pass
    else:
        messenger.info("Results written to: '%s'." % (output_fpath))

if __name__ == '__main__':
    main_cli()
