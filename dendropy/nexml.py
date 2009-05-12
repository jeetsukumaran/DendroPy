#! /usr/bin/env python

############################################################################
##  nexml.py
##
##  Part of the DendroPy library for phylogenetic computing.
##
##  Copyright 2008 Jeet Sukumaran and Mark T. Holder.
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
This module wraps routines needed for reading and writing trees in
NEXML format.
"""

import time
import textwrap
# import cgi
from xml.sax.saxutils import quoteattr
from dendropy import base
from dendropy import datasets
from dendropy import taxa
from dendropy import characters
from dendropy import trees
from dendropy import xmlparser

def protect_attr(x):
#     return cgi.escape(x)
    return quoteattr(x)

############################################################################
## Local Module Methods

def _to_nexml_indent_items(items, indent="", indent_level=0):
    """
    Renders list of items into a string of lines in which each line is
    indented appropriately.
    """
    return '\n'.join(["%s%s" % (indent * indent_level, str(item)) \
                     for item in items])

def _to_nexml_dict(annotes_dict, indent="", indent_level=0):
    "Composes a nexml dict entry, given a python dictionary."
    main_indent = indent * indent_level    
    parts = []
    parts.append('%s<dict id="%s">' % (main_indent,  annotes_dict.oid))
    keyvals = _to_nexml_dict_keyvalues(annotes_dict=annotes_dict,
                                    indent=indent,
                                    indent_level=indent_level+1)
    parts.append(_to_nexml_indent_items(keyvals, indent=indent, indent_level=0))
    parts.append('%s</dict>' % (main_indent))        
    return parts

def _to_nexml_dict_keyvalues(annotes_dict, indent="", indent_level=0):
    """
    Returns a list of lines corresponding to a nexml rendering of a
    dictionary.
    """
    parts = []
    subindent = indent * (indent_level + 0)
    for key, value in annotes_dict.items():
#         parts.append('%s<key>%s</key>' % (subindent, key))
#         anvalue = _to_nexml_dict_value(value=value[0],
#                                        type_hint=value[1],
#                                        indent=indent,
#                                        indent_level=indent_level)

        annote_value = value[0]
        type_hint = value[1]
        if type_hint is None:
            value_type = _to_nexml_dict_value_type(annote_value)
        else:
            value_type = type_hint
        if value_type == 'boolean':
            annote_value = str(annote_value==True).lower()
        if isinstance(annote_value, list):
            value_str = '%s<%s id="%s">%s</%s>' % (subindent,
                                           value_type,
                                           key,
                                           ' '.join([str(item) for item in annote_value]),
                                           value_type)
            parts.append(value_str)
        elif isinstance(annote_value, dict):
            parts.append(_to_nexml_indent_items(_to_nexml_dict(annote_value, indent=indent, indent_level=indent_level), indent, indent_level=0))
        else:
            parts.append('%s<%s id="%s">%s</%s>' % (subindent, value_type, key, str(annote_value), value_type))
    return parts    

def _to_nexml_dict_value_type(value):
    """
    Figures out the value type, and returns and appropriate nexml
    string corresponding to it.
    """
    value_type = 'any'
    if isinstance(value, list):
        # assumes rest of the vector is same type as the first element
        value_type = _to_nexml_dict_value_type(value[0]) + 'vector'
    elif isinstance(value, dict):
        value_type = 'dict'
    else:
        if type(value) == bool:
            value_type = 'boolean'
        elif type(value) == float:
            value_type = 'float'
        elif type(value) == int:
            value_type = 'integer'
        elif type(value) == str:
            value_type = 'string'
        else:
            value_type = 'any'
    return value_type

def _to_nexml_chartype(chartype):
    """
    Returns nexml characters element attribute corresponding to given
    chartype.
    """
#     if chartype == characters.DNA_CHARTYPE:
#         return "nex:DnaSeqs"
#     if chartype == characters.RNA_CHARTYPE:
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

def _from_nexml_tree_length_type(type_attr):
    """
    Given an attribute string read from a nexml tree element, returns
    the python class of the edge length attribute.
    """
    if type_attr == "nex:IntTree":
        return int
    else:
        return float

def _from_nexml_dict_value(value, value_type):
    """
    A text representation of a value of type `type`, where `type`
    is specified in terms of an nexml element, returns the Python
    representation of the value.
    """
    parsed_value = None
    value = value.strip()
    if value_type == "integer":
        try:
            parsed_value = int(value)
        except ValueError:
            raise Exception("Could not parse integer value")
    elif value_type == "float":
        try:
            parsed_value = float(value)
        except ValueError:
            raise Exception("Could not parse float value")
    elif value_type == "boolean":
        try:
            parsed_value = bool(value)
        except ValueError:
            raise Exception("Could not parse boolean value")
    elif value_type == "string":
        try:
            parsed_value = str(value)
        except ValueError:
            raise Exception("Could not parse string value")
    else:
        # what else to do?
        parsed_value = value
    return parsed_value


############################################################################
## NexmlReader

class NexmlReader(datasets.Reader):
    "Implements thinterface for handling NEXML files."
    
    def __init__(self):
        """
        `tree_factory` is a DendroPy TreeFactory class or derived
        object.
        """
        datasets.Reader.__init__(self)
        self.load_time = None
        self.parse_time = None

    ## Implementation of the datasets.Reader interface ##

    def read_dataset(self, file_obj, dataset=None):
        """
        Instantiates and returns a DataSet object based on the
        NEXML-formatted contents read from the file descriptor object
        `fileobj`. If `dataset` is given, its factory methods will be
        used to instantiate objects.
        """
        start = time.clock()
        xml_doc = xmlparser.xml_document(file_obj=file_obj)
        self.load_time = time.clock() - start
        start = time.clock()
        dataset = self.parse_dataset(xml_doc, dataset)
        self.parse_time = time.clock() - start
        return dataset

    ## Following methods are class-specific ###

    def parse_dataset(self, xml_doc, dataset):
        """
        Given an xml_document, parses the XmlElement representation of
        taxon sets, character matrices, and trees into a DataSet object.
        """
        if dataset is None:
            dataset = datasets.Dataset()
        self.parse_taxa_blocks(xml_doc, dataset)
        if self.include_characters:
            self.parse_char_blocks(xml_doc, dataset)
        if self.include_trees:            
            self.parse_trees_blocks(xml_doc, dataset)
        return dataset
        
    def parse_taxa_blocks(self, xml_doc, dataset):
        """
        Given an xml_document, parses the XmlElement representation of
        taxon sets into a TaxaBlocks objects.
        """
        nxt = _NexmlTaxaParser(self.taxa_block_factory, self.taxon_factory)
        for taxa_element in xml_doc.getiterator('otus'):
            taxa_block = nxt.parse_taxa(taxa_element, dataset) 
        
    def parse_char_blocks(self, xml_doc, dataset):
        """
        Given an xml_document, parses the XmlElement representation of
        character sequences into a list of CharacterMatrix objects.
        """
        nxc = _NexmlCharBlockParser()
        for char_block_element in xml_doc.getiterator('characters'):
            nxc.parse_char_block(char_block_element, dataset)

    def parse_trees_blocks(self, xml_doc, dataset):
        """
        Given an xml_document object, parses the XmlElement structural
        representations of a set of NEXML treeblocks (`nex:trees`) and
        returns a TreesBlocks object corresponding to the NEXML.
        """
        nx_tree_parser = _NexmlTreesParser(self.trees_block_factory, self.tree_factory, self.node_factory, self.edge_factory)
        for trees_idx, trees_element in enumerate(xml_doc.getiterator('trees')):
            for tree in nx_tree_parser.parse_trees(trees_element, dataset, trees_idx, add_to_trees_block=True):
                pass

    def iterate_over_trees(file_obj, taxa_block=None, dataset=None):
        "Generator to iterate over trees in file without retaining any in memory."
        xml_doc = xmlparser.xml_document(file_obj=file_obj)
        if dataset is None:
            dataset = datasets.Dataset() or dataset
        if taxa_block is None:
            taxa_block = taxa.TaxaBlock()
        if not (taxa_block in dataset.taxa_blocks):
            dataset.taxa_blocks.append(taxa_block)
        self.parse_taxa_blocks(xml_doc, dataset)
        nx_tree_parser = _NexmlTreesParser()
        for trees_idx, trees_element in enumerate(xml_doc.getiterator('trees')):
            for tree in nx_tree_parser.parse_trees(trees_element, dataset, trees_idx, add_to_trees_block=False):
                yield tree    

class _NexmlElementParser(object):
    "Base parser class: wraps around annotations/dictionary element handling."
    
    def __init__(self):
        "Right now, does nothing ..."
        pass

    def parse_annotations(self, annotated, nxelement):
        """
        Given an nexml element, this looks for a 'dict' child element
        and passes it to the dictionary parse if found. Results are
        placed as attributes of `annotated`.
        """
        pass
#         xml_dict = nxelement.find('dict')
#         if xml_dict:
#             return self.parse_dict(annotated=annotated, xml_dict=xml_dict)

    def parse_dict(self, annotated, xml_dict):
        """
        This parses an xml_dict and sets the attributes of annotable
        correspondingly.
        """
        xml_keys = []
        xml_values = []            
        for child in xml_dict.getchildren():
            if child.tag == 'key':
                xml_keys.append(child)
            else:
                xml_values.append(child)
        if len(xml_keys) > 0 or len(xml_values) > 0:
            if len(xml_keys) == len(xml_values):
                xml_keyvals = dict(zip(xml_keys, xml_values))
                self.parse_keyvals(annotated, xml_keyvals)
            else:
                raise Exception("Unequal numbers of keys and values in annotations")                    

    def parse_keyvals(self, annotated, xml_keyvals):
        """
        Given a dictionary where the keys are nexml dict key
        XmlElements and the values are nexl dict value XmlElements
        corresponding to those keys, this will parse the elements into
        the attributes of an Annotable object.
        """
        for xml_key, xml_value in xml_keyvals.items():
            an_key = xml_key.text
            an_value = None
            if xml_value.tag == 'dict':
                subannotable = base.Annotated()
                self.parse_dict(subannotable, xml_value)
                an_value = subannotable
            elif xml_value.tag.count('vector'):
                an_value = []
                vector_text = xml_value.text
                vector_text = vector_text.strip('\n').strip('\r').strip()
                vector_type = xml_value.tag.replace('vector', '')
                if vector_type == 'dict':
                    ## must handle it here:
                    ## loop through child elements of xml_value,
                    ## parsing the dicts and building up a list of
                    ## Annotable objects
                    raise NotImplementedError
                else:
                    vector_items = vector_text.split()
                    for item in vector_items:
                        an_value.append(_from_nexml_dict_value(item, vector_type))
            else:
                an_value = _from_nexml_dict_value(xml_value.text, xml_value.tag)
            if an_key is not None and an_value is not None:
                setattr(annotated, an_key, an_value)
                annotated.annotate(an_key)

class _NexmlTreesParser(_NexmlElementParser):
    "Parses an XmlElement representation of NEXML format tree blocks."

    def __init__(self, trees_block_factory=None, tree_factory=None, node_factory=None, edge_factory=None):
        "Must be given tree factory to create trees."
        super(_NexmlTreesParser, self).__init__()
        if trees_block_factory is None:
            self.trees_block_factory = trees.TreesBlock
        else:
            self.trees_block_factory = trees_block_factory
        if tree_factory is None:
            self.tree_factory = trees.Tree
        else:
            self.tree_factory = tree_factory
        if node_factory is None:
            self.node_factory = trees.Node
        else:
            self.node_factory = node_factory
        if edge_factory is None:
            self.edge_factory = trees.TreesBlock
        else:
            self.edge_factory = edge_factory

    def parse_trees(self, nxtrees, dataset, trees_idx=None, add_to_trees_block=True):
        """
        Given an XmlElement object representing a NEXML treeblock,
        self.nxtrees (corresponding to a `nex:trees` element), this
        will construct and return a TreesBlock object defined by the
        underlying NEXML. If `add_to_trees_block` is False, then each tree,
        *IS NOT ADDED TO THE DATASET*.
        """
        oid = nxtrees.get('id', "Trees" + str(trees_idx))
        label = nxtrees.get('label', None)
        taxa_id = nxtrees.get('otus', None)
        if taxa_id is None:
            raise Exception("Taxa block not specified for trees block \"%s\"" % oid)
        taxa_block = dataset.find_taxa_block(oid = taxa_id)
        if not taxa_block:
            raise Exception("Taxa block \"%s\" not found" % taxa_id)
        trees_block = dataset.add_trees_block(taxa_block=taxa_block, 
                                            trees_block=self.trees_block_factory(oid=oid, label=label))
        self.parse_annotations(annotated=trees_block, nxelement=nxtrees)                                            
        tree_counter = 0
        for tree_element in nxtrees.getiterator('tree'):
            tree_counter = tree_counter + 1
            oid = tree_element.get('id', tree_counter)
            label = tree_element.get('label', '')
            treeobj = self.tree_factory(oid=oid, label=label)
            treeobj.taxa_block = taxa_block
            tree_type_attr = tree_element.get('{http://www.w3.org/2001/XMLSchema-instance}type')
            treeobj.length_type = _from_nexml_tree_length_type(tree_type_attr)
            self.parse_annotations(annotated=treeobj, nxelement=tree_element)
            nodes = self.parse_nodes(tree_element, taxa_block=treeobj.taxa_block, node_factory=self.node_factory)
            edges = self.parse_edges(tree_element, length_type=treeobj.length_type, edge_factory=self.edge_factory)
            for edge in edges.values():
                # EDGE-ON-ROOT:
                # allow "blank" tail nodes: so we only enforce
                # this check if tail node id is specified
                if edge.tail_node_id and edge.tail_node_id not in nodes:
                    msg = 'Edge "%s" specifies a non-defined ' \
                          'source node ("%s")\nCurrent nodes: %s' % (edge.oid,
                                                                     edge.tail_node_id,
                                                                     (','.join([n for n in nodes])))
                    raise Exception(msg)
                if edge.head_node_id not in nodes:
                    msg = 'Edge "%s" specifies a non-defined ' \
                          'target node ("%s")\nCurrent nodes: %s' % (edge.oid,
                                                                     edge.head_node_id,
                                                                     (','.join([n.oid for n in nodes])))
                    raise Exception(msg)

                if edge.head_node_id and edge.tail_node_id:
                    head_node = nodes[edge.head_node_id]
                    head_node.edge = edge
                    tail_node = nodes[edge.tail_node_id]                    
                    tail_node.add_child(head_node)
                elif edge.head_node_id and not edge.tail_node_id:
                    head_node = nodes[edge.head_node_id]
                    head_node.edge = edge

            # find node(s) without parent
            parentless = []
            for node in nodes.values():
                if node.parent_node == None:
                    parentless.append(node)
                    
            # If one parentless node found, this is the root: we use
            # it as the tree head node. If multiple parentless nodes
            # are found, then we add them all as child_nodes of the
            # existing head node. If none, then we have some sort of
            # cyclicity, and we are not dealing with a tree.
            if len(parentless) == 1:
                treeobj.seed_node = parentless[0]
            elif len(parentless) > 1:
                for node in parentless:
                    treeobj.seed_node.add_child(node)
            else:
                raise Exception("Structural error: tree must be acyclic.")
                
            rootedge = self.parse_root_edge(tree_element, length_type=treeobj.length_type, edge_factory=self.edge_factory)
            if rootedge:
                if rootedge.head_node_id not in nodes:
                    msg = 'Edge "%s" specifies a non-defined ' \
                          'target node ("%s")\nCurrent nodes: %s' % (edge.oid,
                                                                     edge.head_node_id,
                                                                     (','.join([n.oid for n in nodes])))
                    raise Exception(msg)
                else:
                    nodes[rootedge.head_node_id].edge = rootedge
                    ### should we make this node the seed node by rerooting the tree here? ###
            else:
                treeobj.seed_node.edge = None
            if add_to_trees_block:
               trees_block.append(treeobj)
            yield treeobj

    def parse_nodes(self, tree_element, taxa_block, node_factory):
        """
        Given an XmlElement representation of a NEXML tree element,
        (`nex:tree`) this will return a dictionary of DendroPy Node
        objects with the node_id as the key.
        """
        nodes = {}
        for nxnode in tree_element.getiterator('node'):
            node_id = nxnode.get('id', None)
            nodes[node_id] = node_factory()
            nodes[node_id].oid = node_id
            nodes[node_id].label = nxnode.get('label', None)
            taxon_id = nxnode.get('otu', None)
            if taxon_id is not None:
                taxon = taxa_block.get_taxon(oid=taxon_id)
                if not taxon:
                    raise Exception('Taxon with id "%s" not defined in taxa block "%s"' % (taxon_id, taxa_block.oid))
                nodes[node_id].taxon = taxon
            self.parse_annotations(annotated=nodes[node_id], nxelement=nxnode)
        return nodes
        
    def parse_root_edge(self, tree_element, length_type, edge_factory):
        "Returns the edge subtending the root node, or None if not defined."
        rootedge = tree_element.find('rootedge')
        if rootedge:
            edge = edge_factory()
            edge.head_node_id = rootedge.get('target', None)
            edge.oid = rootedge.get('id', 'e' + str(id(edge)))
            edge_length_str = length_type(rootedge.get('length', '0.0'))
            edge.rootedge = True
            edge_length = None
            try:
                edge_length = length_type(edge_length_str)
            except:
                msg = 'Edge "%s" `length` attribute is not a %s' \
                      % (edge.oid, str(length_type))
                raise Exception(msg)
            edge.length = edge_length
            self.parse_annotations(annotated=edge, nxelement=rootedge)            
            return edge
        else:
            return None

    def parse_edges(self, tree_element, length_type, edge_factory):
        """
        Given an XmlElement representation of a NEXML tree element
        this will return a dictionary of DendroPy Edge objects created with
        the edge factory method, self.new_edge, with the oid as
        key. As at this stage, this method knows nothing about defined
        nodes, the Edge tail_node and head_node properties of the
        Edge are not set, but the tail_node_id and head_node_id are.
        """
        edges = {}
        edge_counter = 0        
        for nxedge in tree_element.getiterator('edge'):
            edge = edge_factory()
            edge_counter = edge_counter + 1
            edge.tail_node_id = nxedge.get('source', None)
            edge.head_node_id = nxedge.get('target', None)
            edge.oid = nxedge.get('id', 'e' + str(edge_counter))
            edge_length_str = length_type(nxedge.get('length', 0.0))

            if not edge.tail_node_id:
                msg = 'Edge %d ("%s") does not have a source' \
                      % (edge_counter, edge.oid)
                raise Exception(msg)

            if not edge.head_node_id:
                msg = 'Edge %d ("%s") does not have a target' \
                      % (edge_counter, edge.oid)
                raise Exception(msg)
            edge_length = None
            try:
                edge_length = length_type(edge_length_str)
            except:
                msg = 'Edge %d ("%s") `length` attribute is not a %s' \
                      % (edge_counter, edge.oid, str(length_type))
                raise Exception(msg)
            edge.length = edge_length
            self.parse_annotations(annotated=edge, nxelement=nxedge)            
            edges[edge.oid] = edge
        return edges

class _NexmlTaxaParser(_NexmlElementParser):
    "Parses an XmlElement representation of NEXML taxa blocks."

    def __init__(self, taxa_block_factory=None, taxon_factory=None):
        "Does nothing too useful right now."
        super(_NexmlTaxaParser, self).__init__()
        if taxa_block_factory is None:
            self.taxa_block_factory = taxa.TaxaBlock
        else:
            self.taxa_block_factory = taxa_block_factory   
        if taxon_factory is None:
            self.taxon_factory = taxa.Taxon
        else:
            self.taxon_factory = taxon_factory             

    def parse_taxa(self, nxtaxa, dataset):
        """
        Given an XmlElement representing a nexml taxa block, this
        instantiates and returns a corresponding DendroPy Taxa object.
        """
        oid = nxtaxa.get('id', None)
        label = nxtaxa.get('label', None)
        taxa_block = self.taxa_block_factory(oid=oid, label=label)
        self.parse_annotations(annotated=taxa_block, nxelement=nxtaxa) 
        for idx, nxtaxon in enumerate(nxtaxa.getiterator('otu')):
            taxon = self.taxon_factory(nxtaxon.get('id', "s" + str(idx) ), nxtaxon.get('label', "Taxon" + str(idx)))
            self.parse_annotations(annotated=taxon, nxelement=nxtaxon)
            taxa_block.append(taxon)
        dataset.taxa_blocks.append(taxa_block)
        
class _NexmlCharBlockParser(_NexmlElementParser):
    "Parses an XmlElement representation of NEXML taxa blocks."

    def __init__(self):
        "Does nothing too useful right now."
        super(_NexmlCharBlockParser, self).__init__()
#         if char_block_factory is None:
#             self.char_block_factory = characters.CharBlock()
#         else:
#             self.char_block_factory = char_block_factory
            
    def parse_ambiguous_state(self, nxambiguous, state_alphabet):
        """
        Parses an XmlElement represent an ambiguous discrete character state,
        ("uncertain_state_set")
        and returns a corresponding StateAlphabetElement object.
        """
        state = characters.StateAlphabetElement(oid=nxambiguous.get('id', None),
                                                label=nxambiguous.get('label', None),
                                                symbol=nxambiguous.get('symbol', None),
                                                token=nxambiguous.get('token', None))
        state.member_states = []                                                    
        for nxmember in nxambiguous.getiterator('member'):
            member_state_id = nxmember.get('state', None)
            member_state = state_alphabet.get_state('oid', member_state_id)
            state.member_states.append(member_state)   
        state.multistate = characters.StateAlphabetElement.AMBIGUOUS_STATE    
        return state
            
    def parse_polymorphic_state(self, nxpolymorphic, state_alphabet):
        """
        Parses an XmlElement represent a polymorphic discrete character state, 
        ("polymorphic_state_set")
        and returns a corresponding StateAlphabetElement object.
        """
        state = characters.StateAlphabetElement(oid=nxpolymorphic.get('id', None),
                                                label=nxpolymorphic.get('label', None),
                                                symbol=nxpolymorphic.get('symbol', None),
                                                token=nxpolymorphic.get('token', None))
        state.member_states = []                                                    
        for nxmember in nxpolymorphic.getiterator('member'):
            member_state_id = nxmember.get('state', None)
            member_state = state_alphabet.get_state('oid', member_state_id)
            state.member_states.append(member_state)
        for nxambiguous in nxpolymorphic.getiterator('uncertain_state_set'):
            state.member_states.append(self.parse_ambiguous_state(nxambiguous, state_alphabet))
        state.multistate = characters.StateAlphabetElement.POLYMORPHIC_STATE
        return state
                  
    def parse_state_alphabet(self, nxstates):
        """
        Given an XmlElement representing a nexml definition of (discrete or standard) states 
        ("states"), this returns a corresponding StateAlphabet object.
        """
        
        state_alphabet = characters.StateAlphabet(oid=nxstates.get('id', None),
                                                         label=nxstates.get('label', None))
        for nxstate in nxstates.getiterator('state'):
            state = characters.StateAlphabetElement(oid=nxstate.get('id', None),
                                                    label=nxstate.get('label', None),
                                                    symbol=nxstate.get('symbol', None),
                                                    token=nxstate.get('token', None))        
            state_alphabet.append(state)
        for nxstate in nxstates.getiterator('uncertain_state_set'):
            state_alphabet.append(self.parse_ambiguous_state(nxstate, state_alphabet))
        for nxstate in nxstates.getiterator('polymorphic_state_set'):
            state_alphabet.append(self.parse_polymorphic_state(nxstate, state_alphabet))        
        return state_alphabet
            
    def parse_characters_format(self, nxformat, char_block):
        """
        Given an XmlElement format element ("format"), this parses the 
        state definitions (if any) and characters (column definitions, if any),
        and populates the given char_block accordingly.
        """
        if nxformat is not None:
            for nxstates in nxformat.getiterator('states'):
                char_block.state_alphabets.append(self.parse_state_alphabet(nxstates))
            for nxchars in nxformat.getiterator('char'):
                col = characters.ColumnType(oid=nxchars.get('id', None))
                char_state_set_id = nxchars.get('states')
                if char_state_set_id is not None:
                    state_alphabet = None
                    for state_sets in char_block.state_alphabets:
                        if state_sets.oid == char_state_set_id:
                            state_alphabet = state_sets
                            break
                    if state_alphabet is None:
                        raise Exception("State set '%s' not defined" % char_state_set_id)
                    col.state_alphabet = state_alphabet
                elif hasattr(char_block, "default_state_alphabet") \
                    and char_block.default_state_alphabet is not None:
                    col.state_alphabet = char_block.default_state_alphabet
                char_block.column_types.append(col)
                
    def create_standard_character_alphabet(self, char_block, symbol_list=None):
        """
        Returns a standard character state alphabet based on symbol_list. 
        Defaults to '0' - '9' if not specified.
        """
        if symbol_list is None:
            symbol_list = [str(i) for i in xrange(10)]
        state_alphabet = characters.StateAlphabet()
        for s in symbol_list:
            state_alphabet.append(characters.StateAlphabetElement(symbol=s))
        char_block.state_alphabets.append(state_alphabet)
        char_block.default_state_alphabet = state_alphabet

    def parse_char_block(self, nxchars, dataset):
        """
        Given an XmlElement representing a nexml characters block, this
        instantiates and returns a corresponding DendroPy CharacterMatrix object.
        """
        nxchartype = nxchars.get('{http://www.w3.org/2001/XMLSchema-instance}type', None)
        if nxchartype.startswith('nex:Dna'):
            char_block = characters.DnaCharactersBlock()
        elif nxchartype.startswith('nex:Rna'):
            char_block = characters.RnaCharactersBlock()
        elif nxchartype.startswith('nex:Protein'):
            char_block = characters.ProteinCharactersBlock()            
        elif nxchartype.startswith('nex:Restriction'):
            char_block = characters.RestrictionSitesCharactersBlock()
        elif nxchartype.startswith('nex:Standard'):
            char_block = characters.StandardCharactersBlock()
        elif nxchartype.startswith('nex:Continuous'):
            char_block = characters.ContinuousCharactersBlock()
        else:
            raise NotImplementedError('Character Block %s (\"%s\"): Character type "%s" not supported.' 
                % (char_block.oid, char_block.label, nxchartype))
            
        oid = nxchars.get('id', None)
        label = nxchars.get('label', None)
        char_block.oid = oid
        char_block.label = label   
          
        taxa_id = nxchars.get('otus', None)
        if taxa_id is None:
            raise Exception("Character Block %s (\"%s\"): Taxa block not specified for trees block \"%s\"" % (char_block.oid, char_block.label, char_block.oid))
        taxa_block = dataset.find_taxa_block(oid = taxa_id)
        if not taxa_block:
            raise Exception("Character Block %s (\"%s\"): Taxa block \"%s\" not found" % (char_block.oid, char_block.label, taxa_id))
        char_block.taxa_block = taxa_block
        self.parse_annotations(annotated=char_block, nxelement=nxchars)                
        
        nxformat = nxchars.find('format')        
        if nxformat is not None:
            self.parse_characters_format(nxformat, char_block)
        elif isinstance(char_block, characters.StandardCharactersBlock):
            # default to all integers < 10 as symbols            
            self.create_standard_character_alphabet(char_block)            
            
        matrix = nxchars.find('matrix')
        self.parse_annotations(annotated=char_block.matrix, nxelement=matrix)
        
        if char_block.column_types:
            id_column_map = char_block.id_column_map()
            column_ids = [char.oid for char in char_block.column_types]
        else:
            id_column_map = {}
            column_ids = [] 
            
        for nxrow in matrix.getiterator('row'):
            row_id = nxrow.get('id', None)
            label = nxrow.get('label', None)
            taxon_id = nxrow.get('otu', None)
            taxon = taxa_block.get_taxon(oid=taxon_id)
            if not taxon:
                raise Exception('Character Block %s (\"%s\"): Taxon with id "%s" not defined in taxa block "%s"' % (char_block.oid, char_block.label, taxon_id, taxa_block.oid))                   
                
            character_vector = characters.CharacterDataVector(oid=row_id, label=label, taxon=taxon)
            self.parse_annotations(annotated=character_vector, nxelement=nxrow)
            
            if isinstance(char_block, characters.ContinuousCharactersBlock):
                if nxchartype.endswith('Seqs'):
                    char_block.markup_as_sequences = True
                    seq = nxrow.findtext('seq')
                    if seq is not None:
                        seq = seq.replace('\n\r', ' ').replace('\r\n', ' ').replace('\n', ' ').replace('\r',' ')
                        for char in seq.split(' '):
                            char = char.strip()
                            if char:
                                character_vector.append(characters.CharacterDataCell(value=float(char)))
                else:
                    char_block.markup_as_sequences = False                
                    for nxcell in nxrow.getiterator('cell'):
                        column_id = nxcell.get('char', None)
                        pos_idx = column_ids.index(column_id)
#                         column = id_column_map[column_id]
#                         state = column.state_id_map[cell.get('state', None)]
                        cell = characters.CharacterDataCell(value=float(nxcell.get('state')), column_type=id_column_map[column_id])
                        self.parse_annotations(annotated=cell, nxelement=nxcell)                        
                        character_vector.set_cell_by_index(pos_idx, cell)
            else:
                if nxchartype.endswith('Seqs'):
                    char_block.markup_as_sequences = True            
#                     symbol_state_map = char_block.default_state_alphabet.symbol_state_map()
                    seq = nxrow.findtext('seq')
                    if seq is not None:
                        seq = seq.replace(' ', '').replace('\n', '').replace('\r', '')
                        col_idx = 0
                        for char in seq:
                            symbol_state_map = char_block.column_types[col_idx].state_alphabet.symbol_state_map()                                                
                            if char in symbol_state_map:
                                col_idx += 1
                                state = symbol_state_map[char]
                            else:
                                raise NameError('Character Block %s (\"%s\"): State with symbol "%s" in sequence "%s" not defined' % (char_block.oid, char_block.label, char, seq))
                            character_vector.append(characters.CharacterDataCell(value=state))
                else:
                    char_block.markup_as_sequences = False                
                    id_state_maps = {}
                    for nxcell in nxrow.getiterator('cell'):
                        column_id = nxcell.get('char', None)
                        column = id_column_map[column_id]
                        pos_idx = column_ids.index(column_id)
                        if column_id not in id_state_maps:
                            id_state_maps[column_id] = column.state_alphabet.id_state_map()
                        state = id_state_maps[column_id][nxcell.get('state')]
                        cell = characters.CharacterDataCell(value=state, column_type=column)
                        self.parse_annotations(annotated=cell, nxelement=nxcell)                        
                        character_vector.set_cell_by_index(pos_idx, cell)

            char_block[taxon] = character_vector 

        dataset.char_blocks.append(char_block)
                
class NexmlWriter(datasets.Writer):
    "Implements the DataWriter interface for handling NEXML files."

    def __init__(self):
        "Calls the base class constructor."
        datasets.Writer.__init__(self)
        self.indent = "    "

    ### datasets.Writer interface  ###

    def write_dataset(self, dataset, dest):
        """
        Writes a list of DendroPy Tree objects to a full NEXML
        document.
        """
        self.write_to_nexml_open(dest, indent_level=0)
        self.write_extensions(dataset, dest)
        self.write_taxa_blocks(taxa_blocks=dataset.taxa_blocks, dest=dest)
        self.write_char_blocks(char_blocks=dataset.char_blocks, dest=dest)
        self.write_trees_blocks(trees_blocks=dataset.trees_blocks, dest=dest)
        self.write_to_nexml_close(dest, indent_level=0)

    ### class-specific  ###
            
    def write_taxa_blocks(self, taxa_blocks, dest, indent_level=1):
        "Writes out TaxaBlocks."
        for idx, taxa_block in enumerate(taxa_blocks):
            dest.write(self.indent * indent_level)
            parts = []
            parts.append('otus')
            if taxa_block.oid is not None:
                parts.append('id="%s"' % taxa_block.oid)
            else:
                raise Exception("Taxa block given without ID")
            if taxa_block.label:
                parts.append('label=%s' % protect_attr(taxa_block.label))
            dest.write("<%s>\n" % ' '.join(parts))
            
            # annotate
            self.write_extensions(taxa_block, dest, indent_level=indent_level+1)
            if isinstance(taxa_block, base.Annotated) and taxa_block.has_annotations():
                self.write_annotations(taxa_block, dest, indent_level=indent_level+1)
                
            for taxon in taxa_block:
                dest.write(self.indent * (indent_level+1))
                parts = []
                parts.append('otu')
                if taxon.oid is not None:
                    parts.append('id="%s"' % taxon.oid)
                else:
                    raise Exception("Taxon without ID")
                if taxon.label:
                    parts.append('label=%s' % protect_attr(taxon.label))
                if isinstance(taxon, base.Annotated) and taxon.has_annotations():
                    dest.write("<%s>\n" % ' '.join(parts))
                    self.write_extensions(taxon, dest, indent_level=indent_level+2)
                    self.write_annotations(taxon, dest, indent_level=indent_level+2)
                    dest.write(self.indent * (indent_level+1))
                    dest.write("</otu>\n")
                else:   
                    dest.write("<%s />\n" % ' '.join(parts))
            dest.write(self.indent * indent_level)                
            dest.write('</otus>\n')

    def write_trees_blocks(self, trees_blocks, dest, indent_level=1):
        "Writes out TreesBlocks."
        for idx, trees_block in enumerate(trees_blocks):
            dest.write(self.indent * indent_level)
            parts = []
            parts.append('trees')
            if trees_block.oid is not None:
                parts.append('id="%s"' % trees_block.oid)
            else:
                raise Exception("Tree block given without ID")
            if trees_block.label:
                parts.append('label=%s' % protect_attr(trees_block.label))
            parts.append('otus="%s"' % trees_block.taxa_block.oid)
            dest.write("<%s>\n" % ' '.join(parts))
            
            # annotate
            self.write_extensions(trees_block, dest, indent_level=indent_level+1)
            if isinstance(trees_block, base.Annotated) and trees_block.has_annotations():
                self.write_annotations(trees_block, dest, indent_level=indent_level+1)            
            
            for tree in trees_block:
                self.write_tree(tree=tree, dest=dest, indent_level=2)
            dest.write(self.indent * indent_level)                
            dest.write('</trees>\n')

    def compose_state_definition(self, state, indent_level, member_state=False):
        "Writes out state definition."
        parts = []
        if member_state:
            parts.append('%s<member state="%s"/>' 
                                % (self.indent * indent_level, state.oid))
        elif state.multistate == characters.StateAlphabetElement.SINGLE_STATE:
            parts.append('%s<state id="%s" symbol="%s" />' 
                                % (self.indent * indent_level, state.oid, state.symbol))
        else:
            if state.multistate == characters.StateAlphabetElement.AMBIGUOUS_STATE:
                tag = "uncertain_state_set"
            else:
                tag = "polymorphic_state_set"
                
            parts.append('%s<%s id="%s" symbol="%s">' 
                            % (self.indent * indent_level, tag, state.oid, state.symbol))
            for member in state.member_states:
                parts.extend(self.compose_state_definition(member, indent_level+1, member_state=True))
            parts.append("%s</%s>" % ((self.indent * indent_level), tag))
        return parts        
                                    
    def write_char_blocks(self, char_blocks, dest, indent_level=1):
        "Writes out character matrices."
        for idx, char_block in enumerate(char_blocks):
            dest.write(self.indent * indent_level)
            parts = []
            parts.append('characters')
            if char_block.oid is not None:
                parts.append('id="%s"' % char_block.oid)
            else:
                raise Exception("Character block without ID")
            if char_block.label:
                parts.append('label=%s' % protect_attr(char_block.label))
            parts.append('otus="%s"' % char_block.taxa_block.oid)                    
            if isinstance(char_block, characters.DnaCharactersBlock):
                xsi_datatype = 'nex:Dna'
            elif isinstance(char_block, characters.RnaCharactersBlock):
                xsi_datatype = 'nex:Rna'            
            elif isinstance(char_block, characters.ProteinCharactersBlock):
                xsi_datatype = 'nex:Protein'       
            elif isinstance(char_block, characters.RestrictionSitesCharactersBlock):
                xsi_datatype = 'nex:Restriction'      
            elif isinstance(char_block, characters.StandardCharactersBlock):
                xsi_datatype = 'nex:Standard'  
            elif isinstance(char_block, characters.ContinuousCharactersBlock):
                xsi_datatype = 'nex:Continuous'                  
            else:
                raise Exception("Unrecognized character block data type.")                
            if char_block.markup_as_sequences:
                xsi_markup = 'Seqs'
            else:
                xsi_markup = 'Cells'                
            xsi_type = xsi_datatype + xsi_markup            
            parts.append('xsi:type="%s"' % xsi_type)                        
            dest.write("<%s>\n" % ' '.join(parts))
            
            # annotate
            self.write_extensions(char_block, dest, indent_level=indent_level+1)
            if isinstance(char_block, base.Annotated) and char_block.has_annotations():
                self.write_annotations(char_block, dest, indent_level=indent_level+1)            
            state_alphabet_parts = []
            if hasattr(char_block, "state_alphabets"): #isinstance(char_block, characters.StandardCharactersBlock):
                for state_alphabet in char_block.state_alphabets:
                    state_alphabet_parts.append('%s<states id="%s">' 
                        % (self.indent * (indent_level+2), state_alphabet.oid))
                    for state in state_alphabet:
                        if state.multistate == characters.StateAlphabetElement.SINGLE_STATE:
                            state_alphabet_parts.extend(self.compose_state_definition(state, indent_level+3))
                    for state in state_alphabet:
                        if state.multistate == characters.StateAlphabetElement.POLYMORPHIC_STATE:
                            state_alphabet_parts.extend(self.compose_state_definition(state, indent_level+3))
                    for state in state_alphabet:
                        if state.multistate == characters.StateAlphabetElement.AMBIGUOUS_STATE:
                            state_alphabet_parts.extend(self.compose_state_definition(state, indent_level+3))                        
                    state_alphabet_parts.append('%s</states>' % (self.indent * (indent_level+2)))
            
            column_types_parts = []
            if char_block.column_types:
                for column in char_block.column_types:
                    if column.state_alphabet:
                        column_state = ' states="%s" ' % column.state_alphabet.oid
                    else:
                        column_state = ' '
                    column_types_parts.append('%s<char id="%s"%s/>' 
                        % ((self.indent*(indent_level+1)), column.oid, column_state))
                
            if state_alphabet_parts or column_types_parts:
                dest.write("%s<format>\n" % (self.indent*(indent_level+1)))
                if state_alphabet_parts:
                    dest.write(('\n'.join(state_alphabet_parts)) + '\n')
                if column_types_parts:
                    dest.write(('\n'.join(column_types_parts)) + '\n')
                    pass
                dest.write("%s</format>\n" % (self.indent*(indent_level+1)))
            
           
            dest.write("%s<matrix>\n" % (self.indent * (indent_level+1)))
            
            self.write_extensions(char_block.matrix, dest, indent_level=indent_level+1)
            if isinstance(char_block.matrix, base.Annotated) and char_block.matrix.has_annotations():
                self.write_annotations(char_block.matrix, dest, indent_level=indent_level+1)            
            
            for taxon, row in char_block.matrix.items():
                dest.write(self.indent*(indent_level+2))
                parts = []
                parts.append('row')
                if row.oid is not None:
                    parts.append('id="%s"' % row.oid)
                else:
                    raise Exception("Row without ID")
                if taxon:
                    parts.append('otu="%s"' % taxon.oid)
                dest.write("<%s>\n" % ' '.join(parts))
                
                self.write_extensions(row, dest, indent_level=indent_level+3)
                if isinstance(row, base.Annotated) and row.has_annotations():
                    self.write_annotations(row, dest, indent_level=indent_level+3)            
                

                if char_block.markup_as_sequences:
                    ### actual sequences get written here ###
                    if isinstance(char_block, characters.DnaCharactersBlock) \
                        or isinstance(char_block, characters.RnaCharactersBlock) \
                        or isinstance(char_block, characters.ProteinCharactersBlock) \
                        or isinstance(char_block, characters.RestrictionSitesCharactersBlock):
                        separator = ''
                        break_long_words = True
                    else:
                        # Standard or Continuous
                        separator = ' '
                        break_long_words = False
                    
                    seqlines = textwrap.fill(separator.join([str(c) for c in row]),
                                           width=70,
                                           initial_indent=self.indent*(indent_level+3) + "<seq>",
                                           subsequent_indent=self.indent*(indent_level+4),
                                           break_long_words=break_long_words)
                    seqlines = seqlines + "</seq>\n"                
                    dest.write(seqlines)
                else:                    
                    for cell in row:
                        parts = []
                        parts.append('%s<cell' % (self.indent*(indent_level+3)))
                        if cell.column_type is not None:
                            parts.append('char="%s"' % cell.column_type.oid)
                        if hasattr(cell, "value") and hasattr(cell.value, "oid"):
                            v = cell.value.oid
                        else:
                            v = str(cell.value)
                        parts.append('state="%s"' % v)
                        dest.write(' '.join(parts))
                        if isinstance(cell, base.Annotated) and cell.has_annotations():
                            dest.write('>\n')
                            self.write_extensions(cell, dest, indent_level=indent_level+4)
                            self.write_annotations(cell, dest, indent_level=indent_level+4)            
                            dest.write('%s</cell>' % (self.indent*(indent_level+3)))
                        else:
                            dest.write('/>\n')
                dest.write(self.indent * (indent_level+2))
                dest.write('</row>\n')
            dest.write("%s</matrix>\n" % (self.indent * (indent_level+1)))
            dest.write(self.indent * indent_level)                
            dest.write('</characters>\n')
        
    def write_tree(self, tree, dest, indent_level=0):
        """
        Writes a single DendroPy Tree object as a NEXML nex:tree
        element.
        """
        parts = []
        parts.append('tree')
        if hasattr(tree, 'oid') and tree.oid is not None:
            parts.append('id="%s"' % tree.oid)
        else:
            parts.append('id="%s"' % ("Tree" + str(id(tree))))
        if hasattr(tree, 'label') and tree.label:
            parts.append('label=%s' % protect_attr(tree.label))
        if hasattr(tree, 'length_type') and tree.length_type:
            parts.append('xsi:type="%s"' % _to_nexml_tree_length_type(tree.length_type))
        else:
            parts.append('xsi:type="nex:FloatTree"')
        parts = ' '.join(parts)
        dest.write('%s<%s>\n'
                   % (self.indent * indent_level, parts))   
        # annotate
        self.write_extensions(tree, dest, indent_level=indent_level+1)
        if isinstance(tree, base.Annotated) and tree.has_annotations():
            self.write_annotations(tree, dest, indent_level=indent_level+1)            
            
        for node in tree.preorder_node_iter():
            self.write_node(node=node, dest=dest, indent_level=indent_level+1)
        for edge in tree.preorder_edge_iter():
            self.write_edge(edge=edge, dest=dest, indent_level=indent_level+1)
        dest.write('%s</tree>\n' % (self.indent * indent_level))            
        
    def write_to_nexml_open(self, dest, indent_level=0):
        "Writes the opening tag for a nexml element."
        parts = []
        parts.append('<?xml version="1.0" encoding="ISO-8859-1"?>')
        parts.append('<nex:nexml')
        parts.append('%sversion="0.8"' % (self.indent * (indent_level+1)))
        parts.append('%sxmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"' \
                     % (self.indent * (indent_level+1)))
        parts.append('%sxmlns:xml="http://www.w3.org/XML/1998/namespace"' \
                     % (self.indent * (indent_level+1)))
        parts.append('%sxsi:schemaLocation="http://www.nexml.org/1.0 ../xsd/nexml.xsd"'
                     % (self.indent * (indent_level+1)))
        parts.append('%sxmlns="http://www.nexml.org/1.0"'
                     % (self.indent * (indent_level+1))) 
        parts.append('%sxmlns:nex="http://www.nexml.org/1.0">\n'
                     % (self.indent * (indent_level+1)))                     
        dest.write('\n'.join(parts))

    def write_to_nexml_close(self, dest, indent_level=0):
        "Closing tag for a nexml element."
        dest.write('%s</nex:nexml>\n' % (self.indent*indent_level))

    def write_node(self, node, dest, indent_level=0):
        "Writes out a NEXML node element."
        parts = []
        parts.append('<node')
        parts.append('id="%s"' % node.oid)
        if hasattr(node, 'label') and node.label:
            parts.append('label=%s' % protect_attr(node.label))
        if hasattr(node, 'taxon') and node.taxon:
            parts.append('otu="%s"' % node.taxon.oid)
        parts = ' '.join(parts)
        dest.write('%s%s' % ((self.indent * indent_level), parts))
        if node.has_annotations():
            dest.write('>\n')
            self.write_extensions(node, dest, indent_level=indent_level+1)
            self.write_annotations(node, dest, indent_level=indent_level+1)
            dest.write('%s</node>\n' % (self.indent * indent_level))
        else:
            dest.write(' />\n')
        
    def write_edge(self, edge, dest, indent_level=0):
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
            if hasattr(edge, 'oid') and edge.oid:
                parts.append('id="%s"' % edge.oid)
            # programmatically more efficent to do this in above
            # block, but want to maintain this tag order ...
            if edge.tail_node is not None:
                parts.append('source="%s"' % edge.tail_node.oid)            
            if edge.head_node is not None:
                parts.append('target="%s"' % edge.head_node.oid)
            if hasattr(edge, 'length') and edge.length is not None:
                parts.append('length="%s"' % edge.length)

            # only write if we have more than just the 'edge' and '/' bit
            if len(parts) > 2:
                parts = ' '.join(parts)
                dest.write('%s%s' % ((self.indent * indent_level), parts))
                if edge.has_annotations():
                    dest.write('>\n')
                    self.write_extensions(edge, dest, indent_level=indent_level+1)
                    self.write_annotations(edge, dest,
                                           indent_level=indent_level+1)
                    dest.write('%s</%s>\n' % ((self.indent * indent_level), tag))
                else:
                    dest.write(' />\n')

    def write_annotations(self, annotated, dest, indent_level=0):
        "Writes out annotations for an Annotable object."
        pass
#         if hasattr(annotated, "annotations"):
#             annotes_dict = annotated.annotations()
#             if len(annotes_dict) > 0:
#                 parts = _to_nexml_dict(annotes_dict, self.indent, indent_level)
#                 parts = '\n'.join(parts)
#                 dest.write(parts + '\n')
            
    def write_extensions(self, element, dest, indent_level=0):           
        ### HACK TO SUPPORT RICH STRUCTURED METADATA ###  
        from xml.etree import ElementTree
        for e in element.extensions:
            dest.write(ElementTree.tostring(e))
            dest.write("\n")
 