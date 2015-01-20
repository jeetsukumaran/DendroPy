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
Support for CLI operations.
"""

import re
import argparse
import os
import sys
if sys.hexversion < 0x03000000:
    input_str = raw_input
else:
    input_str = input
import textwrap

import dendropy

def confirm_overwrite(filepath,
                      replace_without_asking=False,
                      file_desc="Output",
                      out=sys.stdout):
    if os.path.exists(filepath):
        if replace_without_asking:
            overwrite = 'y'
        else:
            out.write('%s file already exists: "%s"\n' % (file_desc, filepath))
            overwrite = input_str("Overwrite (y/N)? ")
        if not overwrite.lower().startswith("y"):
            return False
        else:
            return True
    else:
        return True

def show_splash(
        prog_name,
        prog_subtitle,
        prog_version,
        prog_author,
        prog_copyright,
        dest=sys.stderr,
        include_citation=True,
        include_copyright=False,
        additional_citations=None,
        width=70,
        ):
    wrap_width = width - 2
    dendropy_description = dendropy.description()
    lines = []
    lines.append("^"+prog_name)
    lines.append("^"+prog_subtitle)
    lines.append("^Version {}".format(prog_version))
    lines.append("^By {}".format(prog_author))
    lines.append("^Using: {}".format(dendropy_description))
    if include_copyright:
        copyright_lines = []
        copyright_lines.append("{sub_bar1}")
        copyright_text = textwrap.wrap(
                prog_copyright,
                width=wrap_width,
                )
        copyright_lines.extend(copyright_text)
        lines.extend(copyright_lines)
    if include_citation:
        citation_lines = []
        citation_lines.append("{sub_bar1}")
        citation_lines.append("^Citation")
        citation_lines.append("^~~~~~~~~")
        citation_lines.extend(dendropy.citation_info(width=width))
        lines.extend(citation_lines)
        if additional_citations:
            lines.append("")
            lines.append("You should also cite the following")
            for additional_citation in additional_citations:
                lines.append("")
                c = textwrap.wrap(
                        additional_citation,
                        width=wrap_width,
                        initial_indent="  ",
                        subsequent_indent="    ",
                        )
                lines.extend(c)
        lines.append("")
        extra = (
                "Note that in the interests of scientific reproducibility, you "
                "should note in the text of your publications not only the "
                "specific version of the {prog_name} program, but also the "
                "DendroPy library used in your analysis. "
                "For your information, you are running {dendropy_desc}."
                ).format( prog_name=prog_name,
                        prog_version=prog_version,
                        dendropy_desc=dendropy_description,
                        python_version=sys.version)
        lines.extend(textwrap.wrap(extra))
    max_splash_text_width = max(len(i) for i in lines)
    top_bar = "/{}\\".format("=" * (max_splash_text_width + 2))
    bottom_bar = "\\{}/".format("=" * (max_splash_text_width + 2))
    sub_bar1 = "-" * (max_splash_text_width + 2)
    dest.write(top_bar + "\n")
    for line in lines:
        if line == "{sub_bar1}":
            dest.write("+{}+\n".format(sub_bar1))
        else:
            if line.startswith("^"):
                line = line[1:]
                align_char = "^"
            else:
                align_char = "<"
            dest.write("| {:{align_char}{width}} |\n".format(line, align_char=align_char, width=max_splash_text_width))
    dest.write(bottom_bar + "\n")


# from http://stackoverflow.com/a/22157136/268330
class CustomFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
       if text.startswith('R}'):
           return text[2:].splitlines()
       return argparse.HelpFormatter._split_lines(self, text, width)
