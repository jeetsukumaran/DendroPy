#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
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
Benchmarking (raw) NEXUS tokenizer performance.
"""

import sys
import os
import timeit
import argparse

from dendropy.utility import messaging
from dendropy.test.support import pathmap

from dendropy.dataio import nexusprocessing

TREE_FILENAMES = [
    "APG_Angiosperms.nexus",
    "APG_Angiosperms.newick",
    "GEBA.tree.nexus",
    "GEBA.tree.newick",
    "feb032009.trees.nexus",
    "feb032009.trees.newick",
    "Bininda-emonds_2007_mammals.nexus",
    "Bininda-emonds_2007_mammals.newick",
    "Jetz_et_al_2012_Aves.sample.tree.nexus",
    "Jetz_et_al_2012_Aves.sample.tree.newick",
    "Smith_2001_angiosperms.nexus",
    "Smith_2001_angiosperms.newick",
        ]
CHAR_FILENAMES = [
    "actinopterygii.chars.nexus",
    "angiosperms.chars.nexus",
        ]

def tokenizing_fn_factory(src_paths, verbose=False):
    def f():
        for src_path in src_paths:
            if verbose:
                sys.stderr.write("  .. {}\n".format(src_path))
            src = open(src_path, "rU")
            nt = nexusprocessing.NexusTokenizer(src)
            for token in nt:
                pass
    return f

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--target-file",
            type=str,
            dest="target_files",
            default=[],
            action="append",
            help="Path to file to be tokenized; option may be specified multiple times for multiple files.")
    parser.add_argument("-t", "--target-type",
            type=str,
            dest="target_types",
            default=[],
            choices=["trees", "chars", "all"],
            action="append",
            help="Input data file types (default='all' if '-f'/'--file' argument not given); option may be specified multiple times.")
    parser.add_argument("-r", "--repeat",
            type=int,
            default=10,
            help="Repeat each tokenization this number of times (default=%(default)s).")
    parser.add_argument("--delimited-output",
            action="store_true",
            default=False,
            help="Output in tab-delimited instead of aligned format")
    args = parser.parse_args()

    messenger = messaging.ConsoleMessenger(name="-benchmark")

    src_descs = []
    src_paths = []
    results = []

    if args.target_files:
        for f in args.target_files:
            ff = os.path.expanduser(os.path.expandvars(f))
            src_paths.append(ff)
            src_descs.append( ("User", f) )

    if not args.target_types and not args.target_files:
        messenger.info("No sources specified: adding default benchmark target set")
        args.target_types = ["all"]

    if "all" in args.target_types or "trees" in args.target_types:
        for f in TREE_FILENAMES:
            ff = pathmap.tree_source_path(f)
            src_paths.append(ff)
            src_descs.append( ("Trees", f) )

    if "all" in args.target_types or "chars" in args.target_types:
        for f in CHAR_FILENAMES:
            ff = pathmap.char_source_path(f)
            src_paths.append(ff)
            src_descs.append( ("Alignment", f) )

    for src_path, src_desc in zip(src_paths, src_descs):
        messenger.info("Processing: '{}'".format(src_desc[1]))
        t = timeit.Timer(tokenizing_fn_factory([src_path]))
        result = min(t.repeat(args.repeat, 1))
        messenger.info("Best time (of {} repetions): {:.10f} seconds".format(args.repeat, result))
        results.append(result)

    messenger.info("Benchmarking complete: all files processed")

    if args.delimited_output:
        result_template = "{}\t{}\t{:.10f}\n"
        header_template = "{}\t{}\t{}\n"
    else:
        max_len1 = max(len(r[0]) for r in src_descs)
        max_len2 = max(len(r[1]) for r in src_descs)
        col1 = "{{:{}}}".format(max_len1)
        col2 = "{{:{}}}".format(max_len2)
        result_template = "[" + col1 + "]  " + col2 + "  {:.10f}\n"
        header_template = col1 + "    " + col2 + "  {}\n"
    sys.stdout.write(header_template.format("Type", "File", "Seconds"))
    for result, src_desc in zip(results, src_descs):
        sys.stdout.write(result_template.format(src_desc[0], src_desc[1], result))

if __name__ == "__main__":
    main()




