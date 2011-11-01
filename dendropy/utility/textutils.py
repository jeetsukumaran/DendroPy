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
Various text-manipulating and formatting utilities.
"""

import re
import sys
import time

###############################################################################
## NEWICK/NEXUS format support. Placed here instead of `nexustokenizer` so that
## it is available to the entire library without needing to import `nexustokenizer`.

def escape_nexus_token(label, preserve_spaces=False, quote_underscores=True):
    """
    Properly protects a NEXUS token.
    """
    if label is None:
        return ""
    if not preserve_spaces \
            and "_" not in label \
            and not re.search('[\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>\0\t\n]', label):
        label = label.replace(' ', '_').replace('\t', '_')
    elif re.search('[\(\)\[\]\{\}\\\/\,\;\:\=\*\'\"\`\+\-\<\>\0\t\n\r ]', label) \
        or quote_underscores and "_" in label:
        s = label.split("'")
        if len(s) == 1:
            return "'" + label + "'"
        return "'%s'" % "''".join(s)
    return label

def split_as_newick_string(split, taxon_set, preserve_spaces=False, quote_underscores=True):
    """
    Represents a split as a newick string.
    """
    taxlabels = [escape_nexus_token(label, preserve_spaces=preserve_spaces, quote_underscores=quote_underscores) for label in taxon_set.labels()]

    # do not do the root
    if split == 0 or (split == taxon_set.all_taxa_bitmask()):
        return "(%s)" % (",".join(taxlabels))

    idx = 0
    left = []
    right = []
    while split >= 0 and idx < len(taxlabels):
        if split & 1:
            left.append(taxlabels[idx])
        else:
            right.append(taxlabels[idx])
        idx += 1
        split = split >> 1
    assert ( len(left) + len(right) ) == len(taxlabels)
    return "((%s), (%s))" % (", ".join(left), ", ".join(right))

###############################################################################
## Allows string objects to be annotated/decorated with attributes.

class RichString(str):

    def __new__(cls, *args):
        return str.__new__(cls, *args)

###############################################################################
## Various formatters and pretty-printer

def unique_taxon_label_map(taxa, taxon_label_map=None, max_label_len=0, logger=None):
    """
    Given a list of taxa, returns a dictionary with the Taxon objects as
    keys and string labels as values, where the labels are guaranteed to
    be unique. If `taxon_label_map` is pre-populated (as <Taxon> : 'label'),
    then those labels will be used as the basis for the label composition,
    otherwise the original taxon object label will be used. `max_label_len`
    can be used to restrict the maximum length of the labels.
    """
    if taxon_label_map is None:
        taxon_label_map = {}
        for t in taxa:
            taxon_label_map[t] = t.label
    labels = []
    for t in taxon_label_map:
        label = taxon_label_map[t]
        idx = 1
        if label in labels:
            candidate_label = label
            while candidate_label in labels:
                idx += 1
                if max_label_len > 0:
                    k = max_label_len - len(str(idx))
                    if k < 1:
                        raise ValueError("Unable to make labels unique with maximum label length of %d" % max_label_len)
                    candidate_label = label[:k] + str(idx)
                else:
                    candidate_label = label + str(idx)
            label = candidate_label
        labels.append(label)
        taxon_label_map[t] = label
    return taxon_label_map

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
            if border_style == 0:
                vertical_rule = '  '
                horizontal_rule = ''
                rule_junction = ''
            elif border_style == 1:
                vertical_rule = ' '
                horizontal_rule = '-'
                rule_junction = '-'
            else:
                vertical_rule = ' | '
                horizontal_rule = '-'
                rule_junction = '-+-'
            if border_style >= 3:
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
