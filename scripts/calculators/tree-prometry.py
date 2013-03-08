#! /usr/bin/env python


import os
import sys
import optparse
import dendropy
import Queue
import multiprocessing
from cStringIO import StringIO
from dendropy.dataio import tree_source_iter
from dendropy.dataio import multi_tree_source_iter
from dendropy.utility.messaging import ConsoleMessenger
from dendropy.utility.cli import confirm_overwrite, show_splash
from dendropy import prometry

_program_name = "Tree-Prometry"
_program_subtitle = "Tree profile distance calculator"
_program_date = "Mar 02 2013"
_program_version = "Version 1.0.0 (%s)" % _program_date
_program_author = "Jeet Sukumaran"
_program_contact = "jeetsukumaran@gmail.com"
_program_copyright = "Copyright (C) 2013 Jeet Sukumaran.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \

def read_trees(
        tree_filepaths,
        schema,
        measure_node_ages,
        ultrametricity_precision,
        log_frequency,
        messenger):
    messenger.send_info("Running in serial mode.")

    tree_offset = 0
    taxon_set = dendropy.TaxonSet()

    if tree_filepaths is None or len(tree_filepaths) == 0:
        messenger.send_info("Reading trees from standard input.")
        srcs = [sys.stdin]
    else:
        messenger.send_info("%d source(s) to be processed." % len(tree_filepaths))

        # do not want to have all files open at the same time
        #srcs = [open(f, "rU") for f in tree_filepaths]

        # store filepaths, to open individually in loop
        srcs = tree_filepaths

    profile_matrices = []
    edge_len_profiles = prometry.ProfileMatrix("Edge.Lengths")
    profile_matrices.append(edge_len_profiles)

    master_idx = 0

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
                taxon_set=taxon_set)):
            master_idx += 1
            if tidx >= tree_offset:
                if (log_frequency == 1) or (tidx > 0 and log_frequency > 0 and tidx % log_frequency == 0):
                    messenger.send_info("(reading) '%s': tree at offset %d" % (name, tidx), wrap=False)
                # label = "[{:03d}] File {}, Tree {}".format(master_idx, sidx+1, tidx+1, tree.label)
                label = "File {}, Tree {}".format(sidx+1, tidx+1, tree.label)
                edge_len_profiles.add(
                        index=master_idx,
                        label=label,
                        profile_data=[e.length for e in tree.postorder_edge_iter()])
            else:
                if (log_frequency == 1) or (tidx > 0 and log_frequency > 0 and tidx % log_frequency == 0):
                    messenger.send_info("(reading) '%s': tree at offset %d (skipping)" % (name, tidx), wrap=False)
        try:
            src.close()
        except ValueError:
            # "I/O operation on closed file" if we try to close sys.stdin
            pass

    messenger.send_info("Serial processing of %d source(s) completed." % len(srcs))
    return profile_matrices


def main_cli():

    description =  "%s %s %s" % (_program_name, _program_version, _program_subtitle)
    usage = "%prog [options] TREES-FILE [TREES-FILE [TREES-FILE [...]]"

    parser = optparse.OptionParser(usage=usage,
            add_help_option=True,
            version = _program_version,
            description=description)

    source_trees_optgroup = optparse.OptionGroup(parser, "Metric Options")
    parser.add_option_group(source_trees_optgroup)
    source_trees_optgroup.add_option("--from-newick-stream",
            action="store_true",
            dest="from_newick_stream",
            default=False,
            help="trees will be streamed in newick format")
    source_trees_optgroup.add_option("--from-nexus-stream",
            action="store_true",
            dest="from_nexus_stream",
            default=False,
            help="trees will be streamed in NEXUS format")

    metrics_optgroup = optparse.OptionGroup(parser, "Metric Options")
    parser.add_option_group(metrics_optgroup)
    # metrics_optgroup.add_option("--measure-edge-lengths",
    #         action="store_true",
    #         default=True,
    #         help="use edge lengths")
    # metrics_optgroup.add_option("--measure-patristic-distance",
    #         action="store_true",
    #         default=True,
    #         help="use pairwise tip-to-tip path distances")
    metrics_optgroup.add_option("--measure-ages",
            action="store_true",
            default=False,
            help="use node ages (distance from tips); assumes ultrametric trees")
    metrics_optgroup.add_option("--normalize-global",
            action="store_true",
            default=False,
            help="normalize all measures to maximum sum of measures")

    output_filepath_optgroup = optparse.OptionGroup(parser, "Output File Options")
    parser.add_option_group(output_filepath_optgroup)
    output_filepath_optgroup.add_option("-o","--output",
            dest="output_filepath",
            default=None,
            help="path to output file (if not given, will print to standard output)")
    output_filepath_optgroup.add_option("-r", "--replace",
            action="store_true",
            dest="replace",
            default=False,
            help="replace/overwrite output file without asking if it already exists ")

    run_optgroup = optparse.OptionGroup(parser, "Program Run Options")
    parser.add_option_group(run_optgroup)
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


    (opts, args) = parser.parse_args()

    (opts, args) = parser.parse_args()
    if opts.quiet:
        messaging_level = ConsoleMessenger.ERROR_MESSAGING_LEVEL
    else:
        messaging_level = ConsoleMessenger.INFO_MESSAGING_LEVEL
    messenger = ConsoleMessenger(name="Tree-Prometry", messaging_level=messaging_level)

    # splash
    if not opts.quiet:
        show_splash(prog_name=_program_name,
                prog_subtitle=_program_subtitle,
                prog_version=_program_version,
                prog_author=_program_author,
                prog_copyright=_program_copyright,
                dest=sys.stderr,
                extended=False)

    tree_filepaths = []
    if len(args) > 0:
        for fpath in args:
            fpath = os.path.expanduser(os.path.expandvars(fpath))
            if not os.path.exists(fpath):
                messenger.send_error("Terminating due to missing tree file: '{}'".format(fpath))
                sys.exit(1)
            else:
                tree_filepaths.append(fpath)
        if len(tree_filepaths) == 0:
            messenger.send_error("No valid sources of input trees specified. "
                    + "Please provide the path to at least one (valid and existing) file "
                    + "containing trees to profile.")
            sys.exit(1)
    else:
        if not opts.from_newick_stream and not opts.from_nexus_stream:
            messenger.send_info("No sources of input trees specified. "
                    + "Please provide the path to at least one (valid and existing) file "
                    + "containing trees to profile. See '--help' for other options.")
            sys.exit(1)

    profiles = read_trees(
            tree_filepaths=tree_filepaths,
            schema="nexus/newick",
            measure_node_ages=opts.measure_ages,
            ultrametricity_precision=1e-6,
            log_frequency=opts.log_frequency,
            messenger=messenger,
            )

    prometry.summarize_profile_matrices(profiles, sys.stdout)

if __name__ == "__main__":
    main_cli()
