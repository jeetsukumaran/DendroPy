#! /usr/bin/env python

############################################################################
##  cattrees.py
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
Concatenate multiple tree files.
"""

import os
import sys
from optparse import OptionParser
from optparse import OptionGroup
import textwrap

from dendropy.dataio import tree_source_iter
from dendropy.utility.cli import confirm_overwrite, show_splash
from dendropy.utility.messaging import ConsoleMessenger

_program_name = 'CatTrees'
_program_subtitle = 'Phylogenetic Tree File Concatenation'
_program_date = 'Nov 15 2009'
_program_version = 'Version 2.0.0 (%s)' % _program_date
_program_author = 'Jeet Sukumaran'
_program_contact = 'jeetsukumaran@gmail.com'
_program_copyright = "Copyright (C) 2008 Jeet Sukumaran.\n" \
                 "License GPLv3+: GNU GPL version 3 or later.\n" \
                 "This is free software: you are free to change\nand redistribute it. " \
                 "There is NO WARRANTY,\nto the extent permitted by law."

def main_cli():

    description =  '%s %s %s' % (_program_name, _program_version, _program_subtitle)
    usage = "%prog [options] <TREES FILE> [<TREES FILE> [<TREES FILE> [...]]"

    parser = OptionParser(usage=usage, add_help_option=True, version = _program_version, description=description)

    input_optgroup = OptionGroup(parser, 'Input File Options')
    parser.add_option_group(input_optgroup)
    input_optgroup.add_option('-b', '--burnin',
                        action='store',
                        dest='burnin',
                        type='int', # also 'float', 'string' etc.
                        default=0,
                        help='number of trees to skip from the beginning of *each tree file* when counting support [default=%default]')
    input_optgroup.add_option('-s', '--stride',
                        action='store',
                        dest='stride',
                        metavar="STRIDE",
                        type='int', # also 'float', 'string' etc.
                        default=1,
                        help='resample rate: only include one out of every STRIDE trees')

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
    output_filepath_optgroup.add_option('--newick',
                      action='store_true',
                      dest='phylip_format',
                      default=False,
                      help="save results in NEWICK (PHYLIP) format (default is to save in NEXUS format)")
    output_filepath_optgroup.add_option('--phylip',
                      action='store_true',
                      dest='phylip_format',
                      default=False,
                      help="same as --newick")
    output_filepath_optgroup.add_option('-r', '--replace',
                      action='store_true',
                      dest='replace',
                      default=False,
                      help="replace/overwrite output file without asking if it already exists ")

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
    # Tree file idiot checking

    tree_filepaths = []
    missing = False
    for fpath in args:
        fpath = os.path.expanduser(os.path.expandvars(fpath))
        if not os.path.exists(fpath):
            messenger.send_error('File not found: "%s"' % fpath)
            missing = True
        else:
            tree_filepaths.append(fpath)
    if missing:
        messenger.send("")
        if opts.ignore_missing_support:
            pass
        else:
            messenger.send_formatted('Terminating due to missing tree files. '
                   + 'Use the "--ignore-missing-support" option to continue even '
                   + 'if some files are missing.', force=True)
            sys.exit(1)
    if len(tree_filepaths) == 0:
        messenger.send_formatted("No sources of trees specified or could be found. "
        + "Please specify path(s) to tree files to concatenate.", force=True)
        sys.exit(1)

    tree_file_objs = [open(f, "r") for f in tree_filepaths]

    ###################################################
    # Other prepping...

    # output
    if opts.output_filepath is None:
        output_dest = sys.stdout
    else:
        output_fpath = os.path.expanduser(os.path.expandvars(opts.output_filepath))
        if confirm_overwrite(output_fpath, messenger, opts.replace):
            output_dest = open(output_fpath, "w")
        else:
            sys.exit(1)

    ###################################################
    # write nexus header if neccessary

    if opts.phylip_format:
        pass
    else:
        output_dest.write("#NEXUS\n\n")
        output_dest.write("Begin Trees;\n")

    ###################################################
    # Main work begins here

    report = []
    total_trees_added = 0
    for tree_filepath_idx, tree_filepath in enumerate(tree_filepaths):
        messenger.send("-- Reading tree source %d of %d: %s" \
            % (tree_filepath_idx+1, len(tree_filepaths), tree_filepath))
        trees_added = 0
        for tree_count, tree in enumerate(tree_source_iter(stream=open(tree_filepath, "rU"), format='nexus/newick')):
            if tree_count >= opts.burnin and not (tree_count % opts.stride):
                trees_added += 1
                if opts.phylip_format:
                    output_dest.write(tree.as_string(format="newick"))
                else:
                    output_dest.write("tree %d = %s" % (trees_added, tree.as_string(format="newick")))
        total_trees_added += trees_added
        message = ("%s: %d trees in file, sampling 1 tree of every %d trees after %d tree burn-in: %d trees added (current total = %d trees)" \
            % (tree_filepath, tree_count+1, opts.stride, opts.burnin, trees_added, total_trees_added))
        report.append(message)
        messenger.send("   " + message)

    if opts.phylip_format:
        pass
    else:
        output_dest.write("End;\n")
        if opts.include_meta_comments:
            output_dest.write("\n")
            output_dest.write("[Total of %d trees sourced from:]\n" % total_trees_added)
            maxlen = max([len(tf) for tf in report])
            for tf in report:
                output_dest.write("[ %s ]\n" % tf.ljust(maxlen))
        if opts.additional_comments:
            nexus_writer.comment.append("\n")
            nexus_writer.comment.append(opts.additional_comments)

if __name__ == '__main__':
    try:
        main_cli()
    except (KeyboardInterrupt, EOFError), e:
        sys.stderr.write("Terminating (user-abort).\n")
        sys.exit(1)
    except Exception, e:
        sys.stderr.write("Error encountered: %s : %s.\n" % (str(type(e)), str(e)))
        raise # reraise exception, with correct traceback
