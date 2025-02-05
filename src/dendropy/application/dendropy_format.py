#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import dendropy

usage = """\
dendropy-format --from [FORMAT] --to [FORMAT] [OPTIONS] <SOURCE-FILE> > DEST-FILE
"""

def convert(args):
    if args.input_format is None:
        sys.exit("Please specify source format")
    if args.source_file is None or args.source_file == "-":
        sys.stderr.write("(Reading from standard input)\n")
        src = sys.stdin
    else:
        src = open(os.path.expanduser(os.path.expandvars(args.source_file)))
    read_kwargs = {}
    if args.input_format == "phylip-strict":
        args.input_format = "phylip"
        read_kwargs["strict"] = True
        read_kwargs["multispace_delimiter"] = False
        read_kwargs["interleaved"] = False
    elif args.input_format == "phylip-relaxed-singlespace":
        args.input_format = "phylip"
        read_kwargs["strict"] = False
        read_kwargs["multispace_delimiter"] = False
        read_kwargs["interleaved"] = False
    elif args.input_format == "phylip-relaxed-multispace":
        args.input_format = "phylip"
        read_kwargs["strict"] = False
        read_kwargs["multispace_delimiter"] = True
        read_kwargs["interleaved"] = False
    if args.input_format == "phylip-strict-interleaved":
        args.input_format = "phylip"
        read_kwargs["strict"] = True
        read_kwargs["multispace_delimiter"] = False
        read_kwargs["interleaved"] = True
    elif args.input_format == "phylip-relaxed-singlespace-interleaved":
        args.input_format = "phylip"
        read_kwargs["strict"] = False
        read_kwargs["multispace_delimiter"] = False
        read_kwargs["interleaved"] = True
    elif args.input_format == "phylip-relaxed-multispace-interleaved":
        args.input_format = "phylip"
        read_kwargs["strict"] = False
        read_kwargs["multispace_delimiter"] = True
        read_kwargs["interleaved"] = True
    if args.input_format in ("fasta", "phylip"):
        read_kwargs["data_type"] = args.data_type
    with src:
        ds = dendropy.DataSet.get(
                file=src,
                schema=args.input_format,
                **read_kwargs
                )
    write_kwargs = {}
    if args.output_format is None:
        args.output_format = args.input_format
    if args.output_format == "phylip-strict":
        args.output_format = "phylip"
        write_kwargs["strict"] = True
    if args.output_format == "nexus" or args.output_format == "newick":
        if args.unquoted_underscores:
            write_kwargs["unquoted_underscores"] = True
    if args.recode_uncertain is not None:
        operational_state_alphabet = dendropy.DNA_STATE_ALPHABET
        if args.recode_uncertain == "gap":
            convert_to = operational_state_alphabet.gap
        elif args.recode_uncertain == "missing":
            convert_to = operational_state_alphabet.missing
        else:
            raise ValueError(args.recode_uncertain)
        convert_from = set()
        for s in operational_state_alphabet._polymorphic_states:
            convert_from.add(s)
        for s in operational_state_alphabet._ambiguous_states:
            convert_from.add(s)
        for char_matrix in ds.char_matrices:
            for taxon in char_matrix:
                seq = char_matrix[taxon]
                for idx, c in enumerate(seq):
                    if c in convert_from:
                        seq[idx] = convert_to
    dest = sys.stdout
    with dest:
        ds.write(
                file=dest,
                schema=args.output_format,
                **write_kwargs
                )

def to_nexus(args):
    args.output_format = "nexus"
    convert(args)

def to_newick(args):
    args.output_format = "newick"
    convert(args)

def to_nexus(args):
    args.output_format = "nexus"
    convert(args)

def to_nexml(args):
    args.output_format = "nexml"
    convert(args)

def to_fasta(args):
    args.output_format = "fasta"
    convert(args)

def to_phylip(args):
    args.output_format = "phylip"
    convert(args)

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(usage=usage)
    source_options = parser.add_argument_group("Source")
    source_options.add_argument(
            "source_file",
            default=None,
            help="Path to source data.")
    source_options.add_argument(
            "-f", "--from",
            dest="input_format",
            metavar="FORMAT",
            default=None,
            choices=[
                    "fasta",
                    "newick",
                    "nexml",
                    "nexus",
                    "phylip-strict",
                    "phylip-relaxed-singlespace",
                    "phylip-relaxed-multispace",
                    "phylip-strict-interleaved",
                    "phylip-relaxed-singlespace-interleaved",
                    "phylip-relaxed-multispace-interleaved",
                    ],
            help="Format of data source.")
    source_options.add_argument(
            "-d", "--data-type",
            dest="data_type",
            default=None,
            choices=[
                    "dna",
                    "rna",
                    "standard",
                    ],
            help="Type of data")
    parents = [source_options]
    destination_options = parser.add_argument_group("Destination")
    destination_options.add_argument(
            "-t", "--to",
            dest="output_format",
            metavar="FORMAT",
            default=None,
            choices=[
                    "fasta",
                    "newick",
                    "nexml",
                    "nexus",
                    "phylip",
                    "phylip-strict",
                    ],
            help="Format of data source.")
    destination_options.add_argument(
            "-u", "--unquoted-underscores",
            action="store_true",
            default=None,
            help="[NEXUS/Newick:] Do not quote labels with undescores.",
            )
    destination_options.add_argument(
            "--recode-uncertain",
            choices=[
                "missing",
                "gap",
                ],
            default=None,
            help="Recode ambiguous or uncertain characters as missing ('?') or gaps ('-')",
            )
    args = parser.parse_args()
    convert(args)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\n(Terminating due to user interrupt signal)")
