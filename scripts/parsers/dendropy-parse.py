#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import dendropy
from dendropy.utility import texttools

class ParseTarget(object):

    def __init__(self, filepath):
        self.filepath = filepath
        self.fullpath = os.path.expanduser(os.path.expandvars(ipath))
        self.parse_start = None
        self.parse_end = None
        self.parse_success = None

    def register_parse_start(self):
        self.parse_start = datetime.datetime.now()

    def register_parse_end(self, success):
        self.parse_end = datetime.datetime.now()
        self.parse_success = success

    def report(self, out):
        out.write("%s: " % self.filepath)
        if self.parse_success:
            out.write("PASS")
        elif self.parse_success is not None:
            out.write("FAIL")
        else:
            out.write("NONE")
        out.write(": %s\n" % texttools.pretty_timedelta(self.parse_end-self.parse_start))

if __name__ == "__main__":

    parser = optparse.OptionParser(
        usage="%prog [options] file [file [file ...]]",
        description="""\
Reads one or more phylogenetic data files of various formats, optionally
rewriting to standard output. Reports parse times to standard error.
Primary purpose is for testing/profiling parsing operations in DendroPy."""
        )

    parser.add_option('-s', '--schema',
            action="store",
            choices=['nexus',
                     'newick',
                     'nexus/newick',
                     'phylip',
                     'fasta',
                     'nexml',],
            default=None,
            help="format of input data files (mandatory)")

    parser.add_option('-d', '--datatype',
            action="store",
            choices=['dna',
                     'rna',
                     'protein',
                     'standard'],
            default=None,
            help="type of data (required for FASTA  and PHYLIP formats)")

    parser.add_option('--interleaved-phylip',
        action="store_true",
        default=False,
        help="if parsing a PHYLIP format, treat data as interleaved")

    parser.add_option('--strict-phylip',
        action="store_true",
        default=False,
        help="if parsing a PHYLIP format, treat data as strict")

    parser.add_option('-w', '--write',
            action="store",
            choices=['nexus',
                     'newick',
                     'phylip',
                     'phylip-strict',
                     'fasta',
                     'fasta-wrapped',
                     'nexml',],
            default=None,
            metavar='FORMAT',
            help="if given, data will be written to standard output in this format")

    parser.add_option('-e', '--exit-on-error',
        action="store_true",
        default=False,
        help="quit on first parse error encountered; otherwise report to standard error and continue")

    parser.add_option('-n', '--filenames-only',
        action="store_true",
        default=False,
        help="only report paths of files being parsed to standard error")

    parser.add_option('-p', '--passes-only',
        action="store_true",
        default=False,
        help="only report successful parses to standard error")

    parser.add_option('-f', '--fails-only',
        action="store_true",
        default=False,
        help="only report failed parses to standard error")

    parser.add_option('-H', '--hide-error-messages',
        action="store_true",
        default=False,
        help="do not report parse error messages to standard error")

    opts, args = parser.parse_args()

    if len(args) == 0:
        sys.exit("No input files specified")

    if opts.schema is None:
        sys.exit("Schema not specified: '-s' or '--schema' argument required")

    if opts.passes_only and opts.fails_only:
        sys.exit("'--passes-only' and '--fails--only' are mutually-exclusive options")

    total_args = len(args)
    results = []
    global_start_time = datetime.datetime.now()
    for idx, ipath in enumerate(args):
        parse_target = ParseTarget(ipath)
        parse_target.register_parse_start()
        try:
            d = dendropy.DataSet.get_from_path(parse_target.fullpath,
                    schema=opts.schema,
                    datatype=opts.datatype,
                    interleaved=opts.interleaved_phylip,
                    strict=opts.strict_phylip)
            parse_target.register_parse_end(True)
        except Exception, e:
            parse_target.register_parse_end(False)
            if opts.exit_on_error:
                raise
            else:
                if not opts.passes_only:
                    if opts.filenames_only:
                        sys.stderr.write("%s\n" % ipath)
                    else:
                        parse_target.report(sys.stderr)
                        if not opts.hide_error_messages:
                            sys.stderr.write("    %s\n"  % str(e))
        else:
            if not opts.fails_only:
                if opts.filenames_only:
                    sys.stderr.write("%s\n" % ipath)
                else:
                    parse_target.report(sys.stderr)
