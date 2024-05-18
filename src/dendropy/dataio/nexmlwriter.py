#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
Serialization of NeXML-formatted data.
"""

from io import StringIO
import json
import collections
from dendropy.dataio import ioservice

############################################################################
## Local Module Methods

def _safe_unicode(obj, *args):
    """ return the unicode representation of obj """
    try:
        return str(obj, *args)
    except UnicodeDecodeError:
        # obj is byte string
        ascii_text = str(obj).encode('string_escape')
        return str(ascii_text)

def _safe_str(obj):
    """ return the byte string representation of obj """
    try:
        return str(obj)
    except UnicodeEncodeError:
        # obj is unicode
        return str(obj).encode('unicode_escape')

def _protect_attr(x):
#     return cgi.escape(x)
    return json.dumps(_safe_str(x))

def _to_nexml_indent_items(items, indent="", indent_level=0):
    """
    Renders list of items into a string of lines in which each line is
    indented appropriately.
    """
    return '\n'.join(["%s%s" % (indent * indent_level, str(item)) \
                     for item in items])

def _to_nexml_chartype(chartype):
    """
    Returns nexml characters element attribute corresponding to given
    chartype.
    """
#     if chartype == dendropy.DNA_CHARTYPE:
#         return "nex:DnaSeqs"
#     if chartype == dendropy.RNA_CHARTYPE:
#         return "nex:RnaSeqs"
    return None

def _to_nexml_tree_length_type(length_type):
    """
    Returns attribute string for nexml tree type depending on whether
    ``length_type`` is an int or a float.
    """
    if length_type == int:
        return "nex:IntTree"
    elif length_type == float:
        return "nex:FloatTree"
    else:
        raise Exception('Unrecognized value class %s' % length_type)

############################################################################
## NexmlWriter

class NexmlWriter(ioservice.DataWriter):
    "Implements the DataWriter interface for handling NEXML files."

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------

        markup_as_sequences : boolean
            If |True|, then character data will be marked up as sequences
            instead of individual cells. Defaults to |False|.
        suppress_unreferenced_taxon_namespaces: boolean, default: |False|
            If |True|, then when writing |DataSet| objects, any
            |TaxonNamespace| object in the DataSet's ``taxon_namespaces``
            collection will *not* be written as a "TAXA" block if it is not
            referenced by any character matrix (``char_matrices``) or tree list
            (``tree_lists``).
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.

        """

        # base
        ioservice.DataWriter.__init__(self)

        # customization
        self.markup_as_sequences = kwargs.pop("markup_as_sequences", False)
        self.suppress_unreferenced_taxon_namespaces = kwargs.pop("suppress_unreferenced_taxon_namespaces", False)
        self.check_for_unused_keyword_arguments(kwargs)

        # book-keeping
        self.indent = "    "
        self._prefix_uri_tuples = set()
        self._taxon_namespaces_to_write = []
        self._taxon_namespace_id_map = {}
        self._object_xml_id = {}
        self._taxon_id_map = {}
        self._node_id_map = {}
        self._state_alphabet_id_map = {}
        self._state_id_map = {}

    def _write(self,
            stream,
            taxon_namespaces=None,
            tree_lists=None,
            char_matrices=None,
            global_annotations_target=None):

        # reset book-keeping
        self._taxon_namespaces_to_write = []
        self._taxon_namespace_id_map = {}
        self._taxon_id_map = {}
        self._node_id_map = {}
        self._state_alphabet_id_map = {}
        self._state_id_map = {}

        # Destination:
        # Writing to buffer instead of directly to output
        # stream so that all namespaces referenced in metadata
        # can be written
        body = StringIO()

        # comments and metadata
        self._write_annotations_and_comments(global_annotations_target, body, 1)

        # Taxon namespace discovery
        candidate_taxon_namespaces = collections.OrderedDict()
        if self.attached_taxon_namespace is not None:
            candidate_taxon_namespaces[self.attached_taxon_namespace] = True
        elif taxon_namespaces is not None:
            if self.suppress_unreferenced_taxon_namespaces:
                # preload to preserve order
                for tns in taxon_namespaces:
                    candidate_taxon_namespaces[tns] = False
            else:
                for tns in taxon_namespaces:
                    candidate_taxon_namespaces[tns] = True
        for data_collection in (tree_lists, char_matrices):
            if data_collection is not None:
                for i in data_collection:
                    if self.attached_taxon_namespace is None or i.taxon_namespace is self.attached_taxon_namespace:
                        candidate_taxon_namespaces[i.taxon_namespace] = True
        self._taxon_namespaces_to_write = [tns for tns in candidate_taxon_namespaces if candidate_taxon_namespaces[tns]]

        for tns in self._taxon_namespaces_to_write:
            self._write_taxon_namespace(tns, body)

        if char_matrices:
            for char_matrix in char_matrices:
                self._write_char_matrix(char_matrix=char_matrix, dest=body)

        if tree_lists:
            for tree_list in tree_lists:
                self._write_tree_list(tree_list=tree_list, dest=body)

        self._write_to_nexml_open(stream, indent_level=0)
        stream.write(body.getvalue())
        self._write_to_nexml_close(stream, indent_level=0)

    def _write_taxon_namespace(self, taxon_namespace, dest, indent_level=1):
        self._taxon_namespace_id_map[taxon_namespace] = self._get_nexml_id(taxon_namespace)
        dest.write(self.indent * indent_level)
        parts = []
        parts.append('otus')
        parts.append('id="%s"' % self._taxon_namespace_id_map[taxon_namespace])
        if taxon_namespace.label:
            parts.append('label=%s' % _protect_attr(taxon_namespace.label))
        dest.write("<%s>\n" % ' '.join(parts))
        self._write_annotations_and_comments(taxon_namespace, dest, indent_level=indent_level+1)
        for taxon in taxon_namespace:
            dest.write(self.indent * (indent_level+1))
            parts = []
            parts.append('otu')
            self._taxon_id_map[taxon] = self._get_nexml_id(taxon)
            parts.append('id="%s"' % self._taxon_id_map[taxon])
            if taxon.label:
                parts.append('label=%s' % _protect_attr(taxon.label))
            if taxon.has_annotations or (hasattr(taxon, "comments") and taxon.comments):
                dest.write("<%s>\n" % ' '.join(parts))
                # self.write_extensions(taxon, dest, indent_level=indent_level+2)
                self._write_annotations_and_comments(taxon, dest, indent_level=indent_level+2)
                dest.write(self.indent * (indent_level+1))
                dest.write("</otu>\n")
            else:
                dest.write("<%s />\n" % ' '.join(parts))
        dest.write(self.indent * indent_level)
        dest.write('</otus>\n')

    def _write_tree_list(self, tree_list, dest, indent_level=1):
        dest.write(self.indent * indent_level)
        parts = []
        parts.append('trees')
        parts.append('id="%s"' % self._get_nexml_id(tree_list))
        if tree_list.label:
            parts.append('label=%s' % _protect_attr(tree_list.label))
        parts.append('otus="%s"' % self._taxon_namespace_id_map[tree_list.taxon_namespace])
        dest.write("<%s>\n" % ' '.join(parts))
        if tree_list.has_annotations or (hasattr(tree_list, "comments") and tree_list.comments):
            self._write_annotations_and_comments(tree_list, dest,
                    indent_level=indent_level+1)
        for tree in tree_list:
            self._write_tree(tree=tree, dest=dest, indent_level=2)
        dest.write(self.indent * indent_level)
        dest.write('</trees>\n')

    def _compose_state_definition(self, state, state_alphabet, indent_level, member_state=False):
        "Writes out state definition."
        parts = []
        if state not in self._state_id_map:
            self._state_id_map[state] = self._get_nexml_id(state)
        if member_state:
            parts.append('%s<member state="%s"/>'
                                % (self.indent * indent_level, self._state_id_map[state]))
        elif state.state_denomination == state_alphabet.FUNDAMENTAL_STATE:
            parts.append('%s<state id="%s" symbol="%s" />'
                                % (self.indent * indent_level, self._state_id_map[state], state.symbol))
        else:
            if state.state_denomination == state_alphabet.AMBIGUOUS_STATE:
                tag = "uncertain_state_set"
            else:
                tag = "polymorphic_state_set"

            parts.append('%s<%s id="%s" symbol="%s">'
                            % (self.indent * indent_level, tag, self._state_id_map[state], state.symbol))
            for member in state.member_states:
                parts.extend(self._compose_state_definition(member, state_alphabet, indent_level+1, member_state=True))
            parts.append("%s</%s>" % ((self.indent * indent_level), tag))
        return parts

    def _write_char_matrix(self, char_matrix, dest, indent_level=1):
        dest.write(self.indent * indent_level)
        parts = []
        parts.append('characters')
        parts.append('id="%s"' % self._get_nexml_id(char_matrix))
        if char_matrix.label:
            parts.append('label=%s' % _protect_attr(char_matrix.label))
        parts.append('otus="%s"' % self._taxon_namespace_id_map[char_matrix.taxon_namespace])
        if char_matrix.data_type == "dna":
            xsi_datatype = 'nex:Dna'
        elif char_matrix.data_type == "rna":
            xsi_datatype = 'nex:Rna'
        elif char_matrix.data_type == "protein":
            xsi_datatype = 'nex:Protein'
        elif char_matrix.data_type == "restriction":
            xsi_datatype = 'nex:Restriction'
        elif char_matrix.data_type == "standard":
            xsi_datatype = 'nex:Standard'
        elif char_matrix.data_type == "continuous":
            xsi_datatype = 'nex:Continuous'
        else:
            raise Exception("Unrecognized character block data type.")
        if self.markup_as_sequences:
            xsi_markup = 'Seqs'
        else:
            xsi_markup = 'Cells'
        xsi_type = xsi_datatype + xsi_markup
        parts.append('xsi:type="%s"' % xsi_type)
        dest.write("<%s>\n" % ' '.join(parts))

        if char_matrix.has_annotations or (hasattr(char_matrix, "comments") and char_matrix.comments):
            self._write_annotations_and_comments(char_matrix, dest, indent_level=indent_level+1)

        cell_char_type_id_map = self._write_format_section(char_matrix, dest, indent_level=indent_level+1)

        dest.write("%s<matrix>\n" % (self.indent * (indent_level+1)))

        # with new data model,  char_matrix == taxon_seq_map!
        # if char_matrix.taxon_seq_map.has_annotations:
        #     self._write_annotations_and_comments(char_matrix.taxon_seq_map, dest, indent_level=indent_level+1)

        for taxon in char_matrix:
            char_vector = char_matrix[taxon]
            # for col_idx, (char_value, cell_char_type, cell_annotations) in enumerate(char_vector):
            dest.write(self.indent*(indent_level+2))
            parts = []
            parts.append('row')
            parts.append('id="%s"' % self._get_nexml_id(char_vector))
            if taxon is not None:
                parts.append('otu="%s"' % self._taxon_id_map[taxon])
            dest.write("<%s>\n" % ' '.join(parts))
            if char_vector.has_annotations or (hasattr(char_vector, "comments") and char_vector.comments):
                self._write_annotations_and_comments(char_vector, dest, indent_level=indent_level+3)
            if self.markup_as_sequences:
                if char_matrix.data_type in ("dna", "rna", "protein", "restriction", "aa", "amino-acid"):
                    separator = ''
                else:
                    # standard or continuous
                    separator = ' '
                print_count = 1
                dest.write("{}<seq>".format(self.indent * (indent_level+3)))
                for cidx, c in enumerate(char_vector):
                    s = str(c)
                    if not s:
                        raise TypeError("Character %d in char_vector '%s' does not have a symbol defined for its character state:" % (cidx, char_vector.default_oid) \
                                    + " this matrix cannot be written in sequence format (set 'markup_as_sequences' to False)'")
                    if print_count == 1:
                        dest.write("\n{}".format(self.indent * (indent_level+4)))
                    else:
                        dest.write(separator)
                    dest.write(s)
                    if print_count == 58:
                        print_count = 1
                    else:
                        print_count += 1
                dest.write("\n{}</seq>\n".format(self.indent * (indent_level+3)))
            else:
                for col_idx, (char_value, cell_char_type, cell_annotations) in enumerate(char_vector.cell_iter()):
                    parts = []
                    parts.append('%s<cell' % (self.indent*(indent_level+3)))
                    parts.append('char="%s"' % cell_char_type_id_map[ (taxon, col_idx) ])
                    if char_matrix.data_type == "continuous":
                        v = str(char_value)
                    else:
                        v = self._state_id_map[char_value]
                    parts.append('state="%s"' % v)
                    dest.write(' '.join(parts))
                    if cell_annotations is not None:
                        dest.write('>\n')
                        self._write_annotation_set(cell_annotations, dest, indent_level=indent_level+4)
                        dest.write('%s</cell>' % (self.indent*(indent_level+3)))
                    else:
                        dest.write('/>\n')
            dest.write(self.indent * (indent_level+2))
            dest.write('</row>\n')
        dest.write("%s</matrix>\n" % (self.indent * (indent_level+1)))
        dest.write(self.indent * indent_level)
        dest.write('</characters>\n')

    def _write_tree(self, tree, dest, indent_level=0):
        """
        Writes a single DendroPy Tree object as a NEXML nex:tree
        element.
        """
        parts = []
        parts.append('tree')
        parts.append('id="%s"' % self._get_nexml_id(tree))
        if hasattr(tree, 'label') and tree.label:
            parts.append('label=%s' % _protect_attr(tree.label))
        if hasattr(tree, 'length_type') and tree.length_type:
            parts.append('xsi:type="%s"' % _to_nexml_tree_length_type(tree.length_type))
        else:
            parts.append('xsi:type="nex:FloatTree"')
        parts = ' '.join(parts)
        dest.write('%s<%s>\n'
                   % (self.indent * indent_level, parts))
        if tree.has_annotations or (hasattr(tree, "comments") and tree.comments):
            self._write_annotations_and_comments(tree, dest,
                    indent_level=indent_level+1)
        for node in tree.preorder_node_iter():
            self._write_node(
                    node=node,
                    dest=dest,
                    is_root=tree.is_rooted and node is tree.seed_node,
                    indent_level=indent_level+1)
        for edge in tree.preorder_edge_iter():
            self._write_edge(
                    edge=edge,
                    dest=dest,
                    is_root=tree.is_rooted and node is tree.seed_node,
                    indent_level=indent_level+1)
        dest.write('%s</tree>\n' % (self.indent * indent_level))

    def _write_to_nexml_open(self, dest, indent_level=0):
        "Writes the opening tag for a nexml element."
        parts = []
        parts.append('<?xml version="1.0" encoding="ISO-8859-1"?>')
        parts.append('<nex:nexml')
        parts.append('%sversion="0.9"' % (self.indent * (indent_level+1)))
        ensured_namespaces = [
            ["", "http://www.nexml.org/2009"],
            ["xsi", "http://www.w3.org/2001/XMLSchema-instance"],
            ["xml", "http://www.w3.org/XML/1998/namespace"],
            ["nex", "http://www.nexml.org/2009"],
            ["xsd", "http://www.w3.org/2001/XMLSchema#"],
            # ["dendropy", "http://pypi.org/project/DendroPy/"],
                ]
        # parts.append('%sxmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"' \
        #              % (self.indent * (indent_level+1)))
        # parts.append('%sxmlns:xml="http://www.w3.org/XML/1998/namespace"' \
        #              % (self.indent * (indent_level+1)))
        parts.append('%sxsi:schemaLocation="http://www.nexml.org/2009 ../xsd/nexml.xsd"'
                     % (self.indent * (indent_level+1)))
        # parts.append('%sxmlns="http://www.nexml.org/2009"'
        #              % (self.indent * (indent_level+1)))
        # parts.append('%sxmlns:nex="http://www.nexml.org/2009"'
        #              % (self.indent * (indent_level+1)))
        seen_prefixes = {}
        for prefix, uri in self._prefix_uri_tuples:
            if prefix in seen_prefixes:
                if seen_prefixes[prefix] != uri:
                    raise ValueError("Prefix '%s' mapped to multiple namespaces: '%s', '%s'" % (
                        prefix,
                        uri,
                        seen_prefixes[prefix]))
            seen_prefixes[prefix] = uri
            if prefix:
                prefix = ":" + prefix
            parts.append('%sxmlns%s="%s"'
                        % (self.indent * (indent_level+1), prefix, uri))
        for prefix, uri in ensured_namespaces:
            if prefix in seen_prefixes:
                continue
            if prefix:
                prefix = ":" + prefix
            parts.append('%sxmlns%s="%s"'
                        % (self.indent * (indent_level+1), prefix, uri))
        parts.append('>\n')
        dest.write('\n'.join(parts))

    def _write_to_nexml_close(self, dest, indent_level=0):
        "Closing tag for a nexml element."
        dest.write('%s</nex:nexml>\n' % (self.indent*indent_level))

    def _write_node(self, node, dest, is_root, indent_level=0):
        "Writes out a NEXML node element."
        parts = []
        parts.append('<node')
        self._node_id_map[node] = self._get_nexml_id(node)
        parts.append('id="%s"' % self._node_id_map[node])
        if hasattr(node, 'label') and node.label:
            parts.append('label=%s' % _protect_attr(node.label))
        if hasattr(node, 'taxon') and node.taxon:
            parts.append('otu="%s"' % self._taxon_id_map[node.taxon])
        if is_root:
            parts.append('root="true"')
        parts = ' '.join(parts)
        dest.write('%s%s' % ((self.indent * indent_level), parts))
        if node.has_annotations or (hasattr(node, "comments") and node.comments):
            dest.write('>\n')
            self._write_annotations_and_comments(node, dest, indent_level=indent_level+1)
            dest.write('%s</node>\n' % (self.indent * indent_level))
        else:
            dest.write(' />\n')

    def _write_edge(self, edge, dest, is_root, indent_level=0):
        "Writes out a NEXML edge element."
        if edge and edge.head_node:
            parts = []
            if edge.tail_node is not None:
                tag = "edge"
                parts.append('<%s' % tag)
            else:
                # EDGE-ON-ROOT:
                tag = "rootedge"
                parts.append('<%s' % tag)
            parts.append('id="%s"' % self._get_nexml_id(edge))
            # programmatically more efficent to do this in above
            # block, but want to maintain this tag order ...
            if edge.tail_node is not None:
                parts.append('source="%s"' % self._node_id_map[edge.tail_node])
            if edge.head_node is not None:
                parts.append('target="%s"' % self._node_id_map[edge.head_node])
            if hasattr(edge, 'length') and edge.length is not None:
                parts.append('length="%s"' % edge.length)
            if hasattr(edge, 'label') and edge.label:
                parts.append('label=%s' % _protect_attr(edge.label))

            # only write if we have more than just the 'edge' and '/' bit
            if len(parts) > 2:
                parts = ' '.join(parts)
                dest.write('%s%s' % ((self.indent * indent_level), parts))
                if edge.has_annotations or (hasattr(edge, "comments") and edge.comments):
                    dest.write('>\n')
                    self._write_annotations_and_comments(edge, dest, indent_level=indent_level+1)
                    dest.write('%s</%s>\n' % ((self.indent * indent_level), tag))
                else:
                    dest.write(' />\n')

    def _write_annotations_and_comments(self, item, dest, indent_level=0):
        if item is not None:
            self._write_annotations(item, dest, indent_level=indent_level)
            self._write_comments(item, dest, indent_level=indent_level, newline=True)

    def _write_comments(self, commented, dest, indent_level=0, newline=False):
        if hasattr(commented, "comments") and commented.comments:
            if newline:
                post = "\n"
            else:
                post = ""
            for comment in commented.comments:
                dest.write('%s<!-- %s -->%s' % ((self.indent * indent_level),
                    comment, post))

    def _write_annotations(self, annotated, dest, indent_level=0):
        "Writes out annotations for an Annotable object."
        # import sys
        if hasattr(annotated, "annotations"):
            self._write_annotation_set(annotated.annotations, dest, indent_level)

    def _write_annotation_set(self, annotation_set, dest, indent_level=0):
        for annote in annotation_set:
            if annote.is_hidden:
                continue
            dest.write(self._compose_annotation_xml(annote,
                    indent=self.indent,
                    indent_level=indent_level,
                    prefix_uri_tuples=self._prefix_uri_tuples))
            dest.write("\n")

    def _compose_char_type_xml_for_continuous_type(self, indent_level, char_type_id=None):
        if char_type_id is None:
            char_type_id = self._get_nexml_id(object())
        s = ('%s<char id="%s" />'
            % ((self.indent*(indent_level)), char_type_id))
        return char_type_id, s

    def _compose_char_type_xml_for_state_alphabet(self, state_alphabet, indent_level, char_type_id=None):
        if state_alphabet:
            char_type_state = ' states="%s" ' % self._state_alphabet_id_map[state_alphabet]
        else:
            char_type_state = ' '
        if char_type_id is None:
            char_type_id = self._get_nexml_id(object())
        s = ('%s<char id="%s"%s/>'
            % ((self.indent*(indent_level)), char_type_id, char_type_state))
        return char_type_id, s

    def _compose_char_type_xml_for_character_type(self, character_type, indent_level):
        state_alphabet = character_type.state_alphabet
        return self._compose_char_type_xml_for_state_alphabet(
                state_alphabet,
                indent_level=indent_level,
                char_type_id=self._get_nexml_id(character_type))

    def _get_state_alphabet_for_char_matrix(self, char_matrix):
        sa = None
        if char_matrix.default_state_alphabet is not None:
            sa = char_matrix.default_state_alphabet
        elif len(char_matrix.state_alphabets) == 1:
            sa = char_matrix.state_alphabets[0]
        elif len(char_matrix.state_alphabets) > 1:
            raise TypeError("Multiple state alphabets defined for this matrix with no default specified")
        elif len(char_matrix.state_alphabets) == 0:
            raise TypeError("No state alphabets defined for this matrix")
        return sa

    def _write_format_section(self, char_matrix, dest, indent_level):
        format_section_parts = []
        if hasattr(char_matrix, "state_alphabets"): #isinstance(char_matrix, dendropy.StandardCharacterMatrix):
            for state_alphabet in char_matrix.state_alphabets:
                self._state_alphabet_id_map[state_alphabet] = self._get_nexml_id(state_alphabet)
                format_section_parts.append('%s<states id="%s">'
                    % (self.indent * (indent_level+1), self._state_alphabet_id_map[state_alphabet]))
                for state in state_alphabet:
                    if state.state_denomination == state_alphabet.FUNDAMENTAL_STATE:
                        format_section_parts.extend(self._compose_state_definition(state, state_alphabet, indent_level+3))
                for state in state_alphabet:
                    if state.state_denomination == state_alphabet.POLYMORPHIC_STATE:
                        format_section_parts.extend(self._compose_state_definition(state, state_alphabet, indent_level+3))
                for state in state_alphabet:
                    if state.state_denomination == state_alphabet.AMBIGUOUS_STATE:
                        format_section_parts.extend(self._compose_state_definition(state, state_alphabet, indent_level+3))
                format_section_parts.append('%s</states>' % (self.indent * (indent_level+1)))
        cell_char_type_id_map = {}
        char_type_ids_written = set()
        for taxon in char_matrix:
            char_vector = char_matrix[taxon]
            for col_idx, (char_value, cell_char_type, cell_annotations) in enumerate(char_vector.cell_iter()):
                if cell_char_type is None:
                    if char_matrix.data_type == "continuous":
                        char_type_id, char_type_xml = self._compose_char_type_xml_for_continuous_type(indent_level=indent_level+1)
                    else:
                        sa = self._get_state_alphabet_for_char_matrix(char_matrix)
                        assert sa is not None
                        char_type_id, char_type_xml = self._compose_char_type_xml_for_state_alphabet(sa, indent_level=indent_level+1)
                else:
                    char_type_id, char_type_xml = self._compose_char_type_xml_for_character_type(cell_char_type, indent_level=indent_level+1)
                if char_type_id not in char_type_ids_written:
                    format_section_parts.append(char_type_xml)
                    char_type_ids_written.add(char_type_id)
                cell_char_type_id_map[ (taxon, col_idx) ] = char_type_id
        if format_section_parts:
            dest.write("%s<format>\n" % (self.indent*(indent_level)))
            dest.write(('\n'.join(format_section_parts)) + '\n')
            dest.write("%s</format>\n" % (self.indent*(indent_level)))
        return cell_char_type_id_map

    def _get_nexml_id(self, o):
        try:
            return self._object_xml_id[o]
        except KeyError:
            oid = "d{}".format(len(self._object_xml_id))
            self._object_xml_id[o] = oid
            return oid

    def _compose_annotation_xml(self,
            annote,
            indent="",
            indent_level=0,
            prefix_uri_tuples=None):
        parts = ["%s<meta" % (indent * indent_level)]
        value = annote.value
        # if value is not None:
        #     value = _protect_attr(value)
        # else:
        #     value = None
        if isinstance(value, list) or isinstance(value, tuple):
            value = _protect_attr(" ".join(str(v) for v in value))
        elif value is not None:
            value = _protect_attr(value)
        key = annote.prefixed_name
        # assert ":" in key
        if annote.annotate_as_reference:
            parts.append('xsi:type="nex:ResourceMeta"')
            parts.append('rel="%s"' % key)
            if value is not None:
                parts.append('href=%s' % value)
        else:
            parts.append('xsi:type="nex:LiteralMeta"')
            parts.append('property="%s"' % key)
            if value is not None:
                parts.append('content=%s' % value)
            else:
                parts.append('content=""')
        if annote.datatype_hint:
            parts.append('datatype="%s"'% annote.datatype_hint)
        parts.append('id="%s"' % self._get_nexml_id(annote))
        if prefix_uri_tuples is not None:
            prefix_uri_tuples.add((annote.name_prefix, annote.namespace))
        if len(annote.annotations) > 0:
            parts.append(">")
            for a in annote.annotations:
                parts.append("\n" + self._compose_annotation_xml(a, indent=indent, indent_level=indent_level+1, prefix_uri_tuples=prefix_uri_tuples))
            parts.append("\n%s</meta>" % (indent * indent_level))
        else:
            parts.append("/>")
        return " ".join(parts)

