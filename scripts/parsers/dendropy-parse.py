#! /usr/bin/env python

import os
import sys
import optparse
import datetime
import dendropy
from dendropy.utility import textutils

__version__ = "1.0.0"
__title__ = "DENDROPY-PARSE v%s" % __version__

class ParseTarget(object):

    def __init__(self, filepath):
        self.filepath = filepath
        self.fullpath = os.path.expanduser(os.path.expandvars(ipath))
        self.parse_start_time = None
        self.parse_end_time = None
        self.parse_success = None

    def register_parse_start_time(self):
        self.parse_start_time = datetime.datetime.now()

    def register_parse_end_time(self, success):
        self.parse_end_time = datetime.datetime.now()
        self.parse_success = success

    def report(self, out):
        if not out:
            return
        out.write("%s\n    " % self.filepath)
        if self.parse_success:
            out.write("PASS")
        elif self.parse_success is not None:
            out.write("FAIL")
        else:
            out.write("NONE")
        out.write(" (%s)\n" % (self.parse_end_time-self.parse_start_time))

    @property
    def parse_time(self):
        return self.parse_end_time - self.parse_start_time

if __name__ == "__main__":

    parser = optparse.OptionParser(
        version=__title__,
        usage="%prog [options] file [file [file ...]]",
        description="""\
Reads one or more phylogenetic data files of various formats, optionally
rewriting to standard output. Reports parse times to standard error.
Primary purpose is for testing/profiling parsing operations in DendroPy."""
        )

    parser.add_option('--help-schemas',
        action="store_true",
        default=False,
        help="show available input and output schemas")

    parser.add_option('-e', '--exit-on-error',
        action="store_true",
        default=False,
        help="quit on first parse error encountered")

    input_schemas = ['nexus',
                     'newick',
                     'nexus/newick',
                     'phylip',
                     'fasta',
                     'nexml',]

    output_schemas = ['nexus',
                     'newick',
                     'phylip',
                     'phylip-strict',
                     'fasta',
                     'fasta-wrapped',
                     'nexml',]

    datatypes = ['dna',
                 'rna',
                 'protein',
                 'standard',
                 'restriction',
                 'infinite']

    input_opts = optparse.OptionGroup(parser, 'Input')
    parser.add_option_group(input_opts)

    input_opts.add_option('-s', '--schema',
        action="store",
        choices=input_schemas + [s.upper() for s in input_schemas],
        default=None,
        help="format of input data files (mandatory)")

    input_opts.add_option('-d', '--datatype',
        action="store",
        choices=datatypes + [d.upper() for d in datatypes],
        default=None,
        help="type of data (required for FASTA  and PHYLIP formats)")

    input_opts.add_option('--interleaved-phylip',
        action="store_true",
        default=False,
        help="if parsing a PHYLIP format, treat data as interleaved")

    input_opts.add_option('--strict-phylip',
        action="store_true",
        default=False,
        help="if parsing a PHYLIP format, treat data as strict")

    report_opts = optparse.OptionGroup(parser, 'Reporting')
    parser.add_option_group(report_opts)

    report_opts.add_option('-o', '--stdout',
        dest="report_to_stdout",
        action="store_true",
        default=False,
        help="report to standard output instead of standard error")

    report_opts.add_option('-n', '--names-only',
        action="store_true",
        default=False,
        help="only report paths of files being parsed")

    report_opts.add_option('-p', '--passes-only',
        action="store_true",
        default=False,
        help="only report successful parses")

    report_opts.add_option('-f', '--fails-only',
        action="store_true",
        default=False,
        help="only report failed parses")

    report_opts.add_option('-H', '--hide-error-messages',
        action="store_true",
        default=False,
        help="do not report parse error messages")

    report_opts.add_option('-q', '--quiet',
        dest="no_report",
        action="store_true",
        default=False,
        help="do not report anything")

    output_opts = optparse.OptionGroup(parser, 'Output [optional]')
    parser.add_option_group(output_opts)

    output_opts.add_option('-w', '--write',
        dest="rewrite",
        default=None,
        metavar="FILEPATH",
        help="rewrite data to FILEPATH")

    output_opts.add_option('--wo', '--write-stdout',
        dest="rewrite",
        action="store_const",
        const="&1",
        help="rewrite data to standard output")

    output_opts.add_option('--we', '--write-stderr',
        dest="rewrite",
        action="store_const",
        const="&2",
        help="rewrite data to standard error")

    output_opts.add_option('-S', '--rewrite-schema',
        action="store",
        choices=output_schemas + [s.upper() for s in output_schemas],
        dest="rewrite_schema",
        metavar="SCHEMA",
        help="schema to use if rewriting data (defaults to input schema)")

    opts, args = parser.parse_args()

    if opts.help_schemas:
        sys.stdout.write("%s\n" % __title__)
        sys.stdout.write("\nInput Schemas:\n    %s\n" % ("\n    ".join(input_schemas)))
        sys.stdout.write("\nOutput Schemas:\n    %s\n" % ("\n    ".join(output_schemas)))
        sys.exit(0)

    if len(args) == 0:
        sys.exit("Input file(s) not specified")

    if opts.schema is None:
        sys.exit("Schema not specified: '-s' or '--schema' argument required")

    if opts.report_to_stdout and opts.no_report:
        sys.exit("'--stdout' and '--quiet' are mutually-exclusive options")

    if opts.passes_only and opts.fails_only:
        sys.exit("'--passes-only' and '--fails--only' are mutually-exclusive options")

    if opts.rewrite:
        if opts.rewrite == "&1":
            rewrite_dest = sys.stdout
        elif opts.rewrite == "&2":
            rewrite_dest = sys.stderr
        else:
            rewrite_dest = open(os.path.expanduser(os.path.expandvars(opts.rewrite)), "r")
    else:
        rewrite_dest = None

    if opts.report_to_stdout:
        report_stream = sys.stdout
    elif opts.no_report:
        report_stream = open(os.devnull, 'w')
    else:
        report_stream = sys.stderr

    total_args = len(args)
    fails = 0
    passes = 0
    parse_targets = []
    global_start_time = datetime.datetime.now()
    for idx, ipath in enumerate(args):
        parse_target = ParseTarget(ipath)
        parse_targets.append(parse_target)
        input_stream = open(parse_target.fullpath, 'rU')
        parse_target.register_parse_start_time()
        try:
            d = dendropy.DataSet(
                    stream=input_stream,
                    schema=opts.schema,
                    data_type=opts.datatype,
                    interleaved=opts.interleaved_phylip,
                    strict=opts.strict_phylip)
            parse_target.register_parse_end_time(True)
        except Exception, e:
            parse_target.register_parse_end_time(False)
            fails += 1
            if opts.exit_on_error:
                raise
            else:
                if not opts.passes_only and report_stream:
                    if opts.names_only:
                        report_stream.write("%s\n" % ipath)
                    else:
                        report_stream.write("[%d/%d]: FAIL (#%d): %s\n" % (idx+1, total_args, fails, parse_target.fullpath))
                        if not opts.hide_error_messages:
                            report_stream.write("%s\n"  % str(e))
                        report_stream.write("\n")
        else:
            passes += 1
            if not opts.fails_only and report_stream:
                if opts.names_only:
                    report_stream.write("%s\n" % ipath)
                else:
                    report_stream.write("[%d/%d]: PASS (#%d): %s\n" % (idx+1, total_args, passes, parse_target.fullpath))
                    report_stream.write("Parse time: %s\n" % (parse_target.parse_time))
                    report_stream.write("\n")
            if rewrite_dest is not None:
                kwargs = {}
                if opts.rewrite_schema is None:
                    rewrite_schema = opts.schema
                elif opts.rewrite_schema.lower() == "phylip-strict":
                    rewrite_schema = "phylip"
                    kwargs['strict'] = True
                elif opts.rewrite_schema.lower() == "fasta-wrapped":
                    rewrite_schema = "fasta"
                    kwargs['wrap'] = True
                else:
                    rewrite_schema = opts.rewrite_schema
                d.write_to_stream(rewrite_dest, schema=rewrite_schema, **kwargs)
    global_end_time = datetime.datetime.now()
    if not opts.names_only:
        report_stream.write('--\n')
        report_stream.write('%d files processed: %s passes, %s fails\n' % (total_args, passes, fails))
        report_stream.write('Total run time: %s\n' % (global_end_time-global_start_time))
        total_parse_time = datetime.timedelta()
        for p in parse_targets:
            total_parse_time += p.parse_time
        report_stream.write('Total parse time: %s\n' % total_parse_time)
