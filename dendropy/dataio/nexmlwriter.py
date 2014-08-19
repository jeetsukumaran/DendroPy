#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2014 Jeet Sukumaran and Mark T. Holder.
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
Serialization of NeXML-formatted data.
"""

import json
import collections
from dendropy.dataio import ioservice

try:
    from StringIO import StringIO # Python 2 legacy support: StringIO in this module is the one needed (not io)
except ImportError:
    from io import StringIO # Python 3

############################################################################
## Local Module Methods

def _get_nexml_id(obj):
    return "d{}".format(id(obj))

def _safe_unicode(obj, *args):
    """ return the unicode representation of obj """
    try:
        return unicode(obj, *args)
    except UnicodeDecodeError:
        # obj is byte string
        ascii_text = str(obj).encode('string_escape')
        return unicode(ascii_text)

def _safe_str(obj):
    """ return the byte string representation of obj """
    try:
        return str(obj)
    except UnicodeEncodeError:
        # obj is unicode
        return unicode(obj).encode('unicode_escape')

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
    `length_type` is an int or a float.
    """
    if length_type == int:
        return "nex:IntTree"
    elif length_type == float:
        return "nex:FloatTree"
    else:
        raise Exception('Unrecognized value class %s' % length_type)

def _compose_annotation_xml(annote, indent="", indent_level=0, prefix_uri_tuples=None):
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
    parts.append('id="%s"' % _get_nexml_id(annote))
    if prefix_uri_tuples is not None:
        prefix_uri_tuples.add((annote.name_prefix, annote.namespace))
    if len(annote.annotations) > 0:
        parts.append(">")
        for a in annote.annotations:
            parts.append("\n" + _compose_annotation_xml(a, indent=indent, indent_level=indent_level+1, prefix_uri_tuples=prefix_uri_tuples))
        parts.append("\n%s</meta>" % (indent * indent_level))
    else:
        parts.append("/>")
    return " ".join(parts)

############################################################################
## NexmlWriter

class NexmlWriter(ioservice.DataWriter):
    "Implements the DataWriter interface for handling NEXML files."

    def __init__(self, **kwargs):
        """

        Keyword Arguments
        -----------------

        markup_as_sequences : boolean
            If `True`, then character data will be marked up as sequences
            instead of individual cells. Defaults to `False`.

        """

        # base
        ioservice.DataWriter.__init__(self)

        # customization
        self.markup_as_sequences = kwargs.pop("markup_as_sequences", False)
        self.check_for_unused_keyword_arguments(kwargs)

        # book-keeping
        self.indent = "    "
        self._prefix_uri_tuples = set()
        self._taxon_namespaces_to_write = []
        self._taxon_namespace_id_map = {}
        self._taxon_id_map = {}
        self._node_id_map = {}

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

        # self.write_char_matrices(char_matrices=self.dataset.char_matrices, dest=body)
        if tree_lists:
            for tree_list in tree_lists:
                self._write_tree_list(tree_list=tree_list, dest=body)

        self.write_to_nexml_open(stream, indent_level=0)
        stream.write(body.getvalue())
        self.write_to_nexml_close(stream, indent_level=0)

    def _write_taxon_namespace(self, taxon_namespace, dest, indent_level=1):
        self._taxon_namespace_id_map[taxon_namespace] = _get_nexml_id(taxon_namespace)
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
            self._taxon_id_map[taxon] = _get_nexml_id(taxon)
            parts.append('id="%s"' % self._taxon_id_map[taxon])
            if taxon.label:
                parts.append('label=%s' % _protect_attr(taxon.label))
            if taxon.has_annotations:
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
        parts.append('id="%s"' % _get_nexml_id(tree_list))
        if tree_list.label:
            parts.append('label=%s' % _protect_attr(tree_list.label))
        parts.append('otus="%s"' % self._taxon_namespace_id_map[tree_list.taxon_namespace])
        dest.write("<%s>\n" % ' '.join(parts))
        if tree_list.has_annotations:
            self._write_annotations_and_comments(tree_list, dest,
                    indent_level=indent_level+1)
        for tree in tree_list:
            self._write_tree(tree=tree, dest=dest, indent_level=2)
        dest.write(self.indent * indent_level)
        dest.write('</trees>\n')

    def compose_state_definition(self, state, indent_level, member_state=False):
        "Writes out state definition."
        parts = []
        if member_state:
            parts.append('%s<member state="%s"/>'
                                % (self.indent * indent_level, state.default_oid))
        elif state.multistate == dendropy.StateAlphabetElement.SINGLE_STATE:
            parts.append('%s<state id="%s" symbol="%s" />'
                                % (self.indent * indent_level, state.default_oid, state.symbol))
        else:
            if state.multistate == dendropy.StateAlphabetElement.AMBIGUOUS_STATE:
                tag = "uncertain_state_set"
            else:
                tag = "polymorphic_state_set"

            parts.append('%s<%s id="%s" symbol="%s">'
                            % (self.indent * indent_level, tag, state.default_oid, state.symbol))
            for member in state.member_states:
                parts.extend(self.compose_state_definition(member, indent_level+1, member_state=True))
            parts.append("%s</%s>" % ((self.indent * indent_level), tag))
        return parts

    def write_char_matrices(self, char_matrices, dest, indent_level=1):
        "Writes out character matrices."
        for idx, char_matrix in enumerate(char_matrices):
            dest.write(self.indent * indent_level)
            parts = []
            parts.append('characters')
            parts.append('id="%s"' % char_matrix.default_oid)
            if char_matrix.label:
                parts.append('label=%s' % _protect_attr(char_matrix.label))
            parts.append('otus="%s"' % char_matrix.taxon_namespace.default_oid)
            if isinstance(char_matrix, dendropy.DnaCharacterMatrix):
                xsi_datatype = 'nex:Dna'
            elif isinstance(char_matrix, dendropy.RnaCharacterMatrix):
                xsi_datatype = 'nex:Rna'
            elif isinstance(char_matrix, dendropy.ProteinCharacterMatrix):
                xsi_datatype = 'nex:Protein'
            elif isinstance(char_matrix, dendropy.RestrictionSitesCharacterMatrix):
                xsi_datatype = 'nex:Restriction'
            elif isinstance(char_matrix, dendropy.StandardCharacterMatrix):
                xsi_datatype = 'nex:Standard'
            elif isinstance(char_matrix, dendropy.ContinuousCharacterMatrix):
                xsi_datatype = 'nex:Continuous'
            else:
                raise Exception("Unrecognized character block data type.")
            if char_matrix.markup_as_sequences:
                xsi_markup = 'Seqs'
            else:
                xsi_markup = 'Cells'
            xsi_type = xsi_datatype + xsi_markup
            parts.append('xsi:type="%s"' % xsi_type)
            dest.write("<%s>\n" % ' '.join(parts))

            # annotate
#             self.write_extensions(char_matrix, dest, indent_level=indent_level+1)
            self.write_annotations(char_matrix, dest, indent_level=indent_level+1)
            self.write_comments(char_matrix, dest, indent_level=indent_level+1, newline=True)
            state_alphabet_parts = []
            if hasattr(char_matrix, "state_alphabets"): #isinstance(char_matrix, dendropy.StandardCharacterMatrix):
                for state_alphabet in char_matrix.state_alphabets:
                    state_alphabet_parts.append('%s<states id="%s">'
                        % (self.indent * (indent_level+2), state_alphabet.default_oid))
                    for state in state_alphabet:
                        if state.multistate == dendropy.StateAlphabetElement.SINGLE_STATE:
                            state_alphabet_parts.extend(self.compose_state_definition(state, indent_level+3))
                    for state in state_alphabet:
                        if state.multistate == dendropy.StateAlphabetElement.POLYMORPHIC_STATE:
                            state_alphabet_parts.extend(self.compose_state_definition(state, indent_level+3))
                    for state in state_alphabet:
                        if state.multistate == dendropy.StateAlphabetElement.AMBIGUOUS_STATE:
                            state_alphabet_parts.extend(self.compose_state_definition(state, indent_level+3))
                    state_alphabet_parts.append('%s</states>' % (self.indent * (indent_level+2)))

            chartypes_to_add = []
            column_chartype_map = {}
            cell_chartype_map = {}
            for taxon, row in char_matrix.taxon_seq_map.items():
                for col_idx, cell in enumerate(row):
                    chartype = None

                    if hasattr(char_matrix, "state_alphabets"):
                        if col_idx in column_chartype_map:
                            chartype = column_chartype_map[col_idx]
                        elif not hasattr(cell, 'character_type') or cell.character_type is None:
                            chartype_oid = "c%d" % col_idx
                            if char_matrix.default_state_alphabet is not None:
                                chartype = dendropy.CharacterType(state_alphabet=char_matrix.default_state_alphabet, oid=chartype_oid)
                            elif len(char_matrix.state_alphabets) == 1:
                                chartype = dendropy.CharacterType(state_alphabet=char_matrix.state_alphabets[0], oid=chartype_oid)
                            elif len(char_matrix.state_alphabets) > 1:
                                raise TypeError("Character cell %d for taxon %s ('%s') does not have a state alphabet mapping given by the" % (col_idx, taxon.default_oid, taxon.label)\
                                        + " 'character_type' property, and multiple state alphabets are defined for the containing" \
                                        + " character matrix ('%s') with no default specified" % char_matrix.default_oid)
                            elif len(char_matrix.state_alphabets) == 0:
                                raise TypeError("Character cell %d for taxon %s ('%s') does not have a state alphabet mapping given by the" % (col_idx, taxon.default_oid, taxon.label)\
                                        + " 'character_type' property, and no state alphabets are defined for the containing" \
                                        + " character matrix" % char_matrix.default_oid)
                        else:
                            chartype = dendropy.CharacterType(state_alphabet=cell.character_type.state_alphabet, oid="c%d" % col_idx)
                    else:
                        if col_idx not in column_chartype_map:
                            chartype = dendropy.CharacterType(oid="c%d" % col_idx)
                        else:
                            chartype = column_chartype_map[col_idx]

                    if chartype is not None:
                        if chartype not in chartypes_to_add:
                            chartypes_to_add.append(chartype)
                        cell_chartype_map[cell] = chartype
                        column_chartype_map[col_idx] = chartype
                    else:
                        raise TypeError("Cannot create character type mapping: character cell %d for taxon %s ('%s') in character matrix '%s'" %  (col_idx, taxon.default_oid, taxon.label, char_matrix))

            character_types_parts = []
            for column in chartypes_to_add:
                if column.state_alphabet:
                    chartype_state = ' states="%s" ' % column.state_alphabet.default_oid
                else:
                    chartype_state = ' '
                character_types_parts.append('%s<char id="%s"%s/>'
                    % ((self.indent*(indent_level+1)), column.default_oid, chartype_state))

            if state_alphabet_parts or character_types_parts:
                dest.write("%s<format>\n" % (self.indent*(indent_level+1)))
                if state_alphabet_parts:
                    dest.write(('\n'.join(state_alphabet_parts)) + '\n')
                if character_types_parts:
                    dest.write(('\n'.join(character_types_parts)) + '\n')
                    pass
                dest.write("%s</format>\n" % (self.indent*(indent_level+1)))


            dest.write("%s<matrix>\n" % (self.indent * (indent_level+1)))

#             self.write_extensions(char_matrix.taxon_seq_map, dest, indent_level=indent_level+1)
            self.write_annotations(char_matrix.taxon_seq_map, dest, indent_level=indent_level+1)
            self.write_comments(char_matrix.taxon_seq_map, dest, indent_level=indent_level+1, newline=True)

            for taxon, row in char_matrix.taxon_seq_map.items():
                dest.write(self.indent*(indent_level+2))
                parts = []
                parts.append('row')
                parts.append('id="%s"' % row.default_oid)
                if taxon:
                    parts.append('otu="%s"' % taxon.default_oid)
                dest.write("<%s>\n" % ' '.join(parts))

#                 self.write_extensions(row, dest, indent_level=indent_level+3)
                self.write_annotations(row, dest, indent_level=indent_level+3)
                self.write_comments(row, dest, indent_level=indent_level+3, newline=True)

                if ( (self.markup_as_sequences is not None and self.markup_as_sequences is False)
                        or (hasattr(char_matrix, 'markup_as_sequences') and not char_matrix.markup_as_sequences)
                        ):
                    for idx, cell in enumerate(row):
                        parts = []
                        parts.append('%s<cell' % (self.indent*(indent_level+3)))
                        parts.append('char="%s"' % cell_chartype_map[cell].default_oid)
                        if hasattr(cell, "value") and hasattr(cell.value, "oid"):
                            v = cell.value.default_oid
                        else:
                            v = str(cell.value)
                        parts.append('state="%s"' % v)
                        dest.write(' '.join(parts))
                        if isinstance(cell, dendropy.AnnotatedDataObject) and len(cell.annotations) > 0:
                            dest.write('>\n')
                            # self.write_extensions(cell, dest, indent_level=indent_level+4)
                            self.write_annotations(cell, dest, indent_level=indent_level+4)
                            self.write_comments(cell, dest, indent_level=indent_level+4, newline=True)
                            dest.write('%s</cell>' % (self.indent*(indent_level+3)))
                        else:
                            dest.write('/>\n')
                else:
                    ### actual sequences get written here ###
                    if isinstance(char_matrix, dendropy.DnaCharacterMatrix) \
                        or isinstance(char_matrix, dendropy.RnaCharacterMatrix) \
                        or isinstance(char_matrix, dendropy.ProteinCharacterMatrix) \
                        or isinstance(char_matrix, dendropy.RestrictionSitesCharacterMatrix):
                        separator = ''
                        break_long_words = True
                    else:
                        # Standard or Continuous
                        separator = ' '
                        break_long_words = False

                    if isinstance(char_matrix, dendropy.DiscreteCharacterMatrix):
                        seq_symbols = []
                        for cidx, c in enumerate(row):
                            if c.value.symbol is None:
                                raise TypeError("Character %d in row '%s' does not have a symbol defined for its character state:" % (cidx, row.default_oid) \
                                            + " this matrix cannot be written in sequence format (set 'markup_as_sequences' to False)'")
                            seq_symbols.append(c.value.symbol)
                    else:
                        seq_symbols = [str(c) for c in row]
                    seqlines = textwrap.fill(separator.join(seq_symbols),
                                           width=70,
                                           initial_indent=self.indent*(indent_level+3) + "<seq>",
                                           subsequent_indent=self.indent*(indent_level+4),
                                           break_long_words=break_long_words)
                    seqlines = seqlines + "</seq>\n"
                    dest.write(seqlines)
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
        parts.append('id="%s"' % _get_nexml_id(tree))
        if hasattr(tree, 'label') and tree.label:
            parts.append('label=%s' % _protect_attr(tree.label))
        if hasattr(tree, 'length_type') and tree.length_type:
            parts.append('xsi:type="%s"' % _to_nexml_tree_length_type(tree.length_type))
        else:
            parts.append('xsi:type="nex:FloatTree"')
        parts = ' '.join(parts)
        dest.write('%s<%s>\n'
                   % (self.indent * indent_level, parts))
        if tree.has_annotations:
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

    def write_to_nexml_open(self, dest, indent_level=0):
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
            # ["dendropy", "http://packages.python.org/DendroPy/"],
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

    def write_to_nexml_close(self, dest, indent_level=0):
        "Closing tag for a nexml element."
        dest.write('%s</nex:nexml>\n' % (self.indent*indent_level))

    def _write_node(self, node, dest, is_root, indent_level=0):
        "Writes out a NEXML node element."
        parts = []
        parts.append('<node')
        self._node_id_map[node] = _get_nexml_id(node)
        parts.append('id="%s"' % self._node_id_map[node])
        if hasattr(node, 'label') and node.label:
            parts.append('label=%s' % _protect_attr(node.label))
        if hasattr(node, 'taxon') and node.taxon:
            parts.append('otu="%s"' % self._taxon_id_map[node.taxon])
        if is_root:
            parts.append('root="true"')
        parts = ' '.join(parts)
        dest.write('%s%s' % ((self.indent * indent_level), parts))
        if node.has_annotations:
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
            parts.append('id="%s"' % _get_nexml_id(edge))
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
                if edge.has_annotations:
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
            for annote in annotated.annotations:
                # sys.stderr.write("{}\t\t{}\n".format(annote.name_prefix, annote.namespace))
                if annote.is_hidden:
                    continue
                # self._prefix_uri_tuples.add((annote.name_prefix, annote.namespace))
                dest.write(_compose_annotation_xml(annote, indent=self.indent, indent_level=indent_level, prefix_uri_tuples=self._prefix_uri_tuples))
                dest.write("\n")
