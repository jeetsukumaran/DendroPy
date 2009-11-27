#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
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
###############################################################################

"""
Various text-manipulating and formatting utilities.
"""

import re
import sys
import time

###############################################################################
## Properly protects a NEXUS token. Placed here instead of `nexustokenizer` so that
## it is available to the entire library within needing to import `nexustokenizer`.

def escape_nexus_token(label, spaces_to_underscores=False):
    if label is None:
        return ""
    if spaces_to_underscores and not re.search('[\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>\0]', label):
        label = label.replace(' ', '_').replace('\t', '_')
    elif re.search('[\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>\0\t\n\r ]', label):
        s = label.split("'")
        if len(s) == 1:
            return "'" + label + "'"
        return "'%s'" % "''".join(s)
    return label

###############################################################################
## Allows string objects to be annotated/decorated with attributes.

class RichString(str):

    def __new__(cls, *args):
        return str.__new__(cls, *args)

###############################################################################
## Various formatters and pretty-printer

if sys.version_info[1] > 5 or sys.version_info[0] > 2:
    def int_to_bitstring(n):
        assert n >= 0
        return bin(n)[2:]
else:
    def int_to_bitstring(n):
        m = 1
        sl = []
        while m <= n:
            if m & n:
                sl.insert(0, '1')
            else:
                sl.insert(0, '0')
            m <<= 1
        return "".join(sl)

def pretty_timestamp(t=None, style=0):
    if t is None:
        t = time.localtime()
    if style == 0:
        return time.strftime("%Y-%m-%d", t)
    else:
        return time.strftime("%Y%m%d%H%M%S", t)

def pretty_timedelta(timedelta):
    hours, mins, secs = str(timedelta).split(":")
    return("%s hour(s), %s minute(s), %s second(s)" % (hours, mins, secs))

def format_dict_table(rows, column_names=None, max_column_width=None, border_style=2):
    """
    Returns a string representation of a tuple of dictionaries in a
    table format. This method can read the column names directly off the
    dictionary keys, but if a tuple of these keys is provided in the
    'column_names' variable, then the order of column_names will follow
    the order of the fields/keys in that variable.
    """
    if column_names or len(rows) > 0:
        lengths = {}
        rules = {}
        if column_names:
            column_list = column_names
        else:
            try:
                column_list = rows[0].keys()
            except:
                column_list = None
        if column_list:
            # characters that make up the table rules
            border_style = int(border_style)
            #border_style = 0
            if border_style >= 1:
                vertical_rule = ' | '
                horizontal_rule = '-'
                rule_junction = '-+-'
            else:
                vertical_rule = '  '
                horizontal_rule = ''
                rule_junction = ''
            if border_style >= 2:
                left_table_edge_rule = '| '
                right_table_edge_rule = ' |'
                left_table_edge_rule_junction = '+-'
                right_table_edge_rule_junction = '-+'
            else:
                left_table_edge_rule = ''
                right_table_edge_rule = ''
                left_table_edge_rule_junction = ''
                right_table_edge_rule_junction = ''

            if max_column_width:
                column_list = [c[:max_column_width] for c in column_list]
                trunc_rows = []
                for row in rows:
                    new_row = {}
                    for k in row.keys():
                        new_row[k[:max_column_width]] = str(row[k])[:max_column_width]
                    trunc_rows.append(new_row)
                rows = trunc_rows

            for col in column_list:
                rls = [len(str(row[col])) for row in rows]
                lengths[col] = max(rls+[len(col)])
                rules[col] = horizontal_rule*lengths[col]

            template_elements = ["%%(%s)-%ss" % (col, lengths[col]) for col in column_list]
            row_template = vertical_rule.join(template_elements)
            border_template = rule_junction.join(template_elements)
            full_line = left_table_edge_rule_junction + (border_template % rules) + right_table_edge_rule_junction
            display = []
            if border_style > 0:
                display.append(full_line)
            display.append(left_table_edge_rule + (row_template % dict(zip(column_list, column_list))) + right_table_edge_rule)
            if border_style > 0:
                display.append(full_line)
            for row in rows:
                display.append(left_table_edge_rule + (row_template % row) + right_table_edge_rule)
            if border_style > 0:
                display.append(full_line)
            return "\n".join(display)
        else:
            return ''
    else:
        return ''
