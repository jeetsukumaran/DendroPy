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
Implementation of NeXML-schema data reader and writer.
"""

from xml.sax import saxutils
from cStringIO import StringIO
import time
import textwrap
import re

from dendropy.utility import containers
from dendropy.utility import iosys
from dendropy.utility import error
from dendropy.dataio import xmlparser
import dendropy

SUPPORTED_NEXML_NAMESPACES = ('http://www.nexml.org/1.0', 'http://www.nexml.org/2009')

############################################################################
## Local Module Methods

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
    return saxutils.quoteattr(_safe_str(x))

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

def _from_nexml_tree_length_type(type_attr):
    """
    Given an attribute string read from a nexml tree element, returns
    the python class of the edge length attribute.
    """
    if type_attr == "nex:IntTree":
        return int
    else:
        return float

def _compose_annotation_xml(annote, indent="", indent_level=0, prefix_uri_tuples=None):
    parts = ["%s<meta" % (indent * indent_level)]
    value = annote.value
    if value is not None:
        value = _protect_attr(value)
    else:
        value = None
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
    parts.append('id="%s"' % annote.default_oid)
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

class _AnnotationParser(object):

    def __init__(self, namespace_registry=None):
        self.namespace_registry = namespace_registry

    def parse_annotations(self, annotated, nxelement):
        attrib = nxelement.attrib
        xml_type = nxelement.parse_type()
        if xml_type == 'LiteralMeta':
            value = attrib.get("content", None)
            key = attrib.get("property", None)
            annotate_as_reference = False
        else:
            value = attrib.get("href", None)
            key = attrib.get("rel", None)
            annotate_as_reference = True
        datatype_hint = attrib.get("datatype", None)
        if key is None:
            raise ValueError("Could not determine property/rel for meta element: %s\n%s" % (nxelement, attrib))
        name_prefix, name = dendropy.Annotation.parse_prefixed_name(key)
        try:
            namespace = self.namespace_registry.prefix_namespace_map[name_prefix]
        except KeyError:
            raise ValueError("CURIE-standard prefix '%s' not defined in document: %s" % (name_prefix, self.namespace_registry))
        a = annotated.annotations.add_new(
                name=name,
                value=value,
                datatype_hint=datatype_hint,
                name_prefix=name_prefix,
                namespace=namespace,
                name_is_prefixed=False,
                annotate_as_reference=annotate_as_reference)
        top_annotations = [i for i in nxelement.findall_annotations()]
        for annotation in top_annotations:
            self.parse_annotations(a, annotation)
    # print _compose_annotation_xml(a)

############################################################################
## Abstraction Layer Around XML


class NexmlElement(xmlparser.XmlElement):

    def __init__(self, element, default_namespace=None):
        # if default_namespace is None:
        #     default_namespace = NexmlReader.DEFAULT_NEXML_NAMESPACE
        # else:
        #     default_namespace = default_namespace
        xmlparser.XmlElement.__init__(self,
                element=element,
                default_namespace=default_namespace)
        self.type_parse_pattern = re.compile(r"([A-Za-z0-9]+?):(.+)")

    ## Annotations ##

    def findall_annotations(self):
        return self.namespaced_findall("meta")

    ## TaxonSet Elements ##

    def iter_otus(self):
        return self.namespaced_getiterator("otus")

    def findall_otu(self):
        return self.namespaced_findall("otu")

    ## Characters ##

    def iter_characters(self):
        return self.namespaced_getiterator("characters")

    def find_char_format(self):
        return self.namespaced_find("format")

    def find_char_matrix(self):
        return self.namespaced_find("matrix")

    def findall_multistate_member(self):
        return self.namespaced_findall("member")

    def findall_polymorphic_state_set(self):
        return self.namespaced_findall("polymorphic_state_set")

    def findall_uncertain_state_set(self):
        return self.namespaced_findall("uncertain_state_set")

    def findall_char_state(self):
        return self.namespaced_findall("state")

    def findall_char_states(self):
        return self.namespaced_findall("states")

    def findall_char(self):
        return self.namespaced_findall("char")

    def findall_char_row(self):
        return self.namespaced_findall("row")

    def findall_char_cell(self):
        return self.namespaced_findall("cell")

    def find_char_seq(self):
        return self.namespaced_findtext("seq")

    ## Trees ##

    def iter_trees(self):
        return self.namespaced_getiterator("trees")

    def findall_tree(self):
        return self.namespaced_findall("tree")

    def findall_node(self):
        return self.namespaced_findall("node")

    def findall_edge(self):
        return self.namespaced_findall("edge")

    def find_rootedge(self):
        return self.namespaced_find("rootedge")

    def parse_type(self):
        type_value = self._element.get('{http://www.w3.org/2001/XMLSchema-instance}type', None)
        if type_value is None:
            raise ValueError("Type not specified for element '%s'" % self._element.get("id", None))
        m = self.type_parse_pattern.match(type_value)
        if m:
            return m.groups()[1]
        else:
            return type_value

############################################################################
## NexmlReader and Supporting Classes

class NexmlReader(iosys.DataReader, _AnnotationParser):
    "Implements thinterface for handling NEXML files."

    DEFAULT_NEXML_NAMESPACE = "http://www.nexml.org/2009"

    def __init__(self, **kwargs):
        """
        See `iosys.IOService.__init__` and `iosys.DataReader.__init__` for
        kwargs.
        ``default_namespace`` allows for parsing of NeXML files with
        different base or default namespaces.
        """
        iosys.DataReader.__init__(self, **kwargs)
        self.load_time = None
        self.parse_time = None
        self.id_taxon_set_map = {}
        self.id_taxon_map = {}
        self.default_namespace = kwargs.get("default_namespace", NexmlReader.DEFAULT_NEXML_NAMESPACE)
        self.case_sensitive_taxon_labels = kwargs.get('case_sensitive_taxon_labels', False)
        _AnnotationParser.__init__(self)

    ## Implementation of the datasets.Reader interface ##

    def subelement_factory(self, element):
        return NexmlElement(element, default_namespace=self.default_namespace)

    def read(self, stream):
        """
        Instantiates and returns a DataSet object based on the
        NEXML-formatted contents read from the file descriptor object
        `stream`.
        """
        start = time.clock()
        xml_doc = xmlparser.XmlDocument(file_obj=stream,
                subelement_factory=self.subelement_factory)
        # import xml.etree.ElementTree as xx
        self.namespace_registry = xml_doc.namespace_registry
        self.load_time = time.clock() - start
        start = time.clock()
        dataset = self.parse_dataset(xml_doc)
        self.parse_time = time.clock() - start
        return self.dataset

    ## Following methods are class-specific ###

    def parse_dataset(self, xml_doc):
        """
        Given an xml_document, parses the XmlElement representation of
        taxon sets, character matrices, and trees into a DataSet object.
        """
        if self.dataset is None:
            self.dataset = dendropy.DataSet()
        xml_root = xml_doc.root
        self.parse_taxon_sets(xml_root, self.dataset)
        if not self.exclude_chars:
            self.parse_char_matrices(xml_root, self.dataset)
        if not self.exclude_trees:
            self.parse_tree_lists(xml_root, self.dataset)
        top_annotations = [i for i in xml_root.findall_annotations()]
        for annotation in top_annotations:
            self.parse_annotations(self.dataset, annotation)
        return self.dataset

    def parse_taxon_sets(self, xml_root, dataset):
        """
        Given an xml_document, parses the XmlElement representation of
        taxon sets into a TaxonSets objects.
        """
        # nxt = _NexmlTaxaParser(self.namespace_registry)
        # if len(taxon_set_elements) > 1:
        #     if self.dataset is None:
        #         self.dataset = dendropy.DataSet()
        #     else:
        #         if self.dataset.attached_taxon_set is not None:
        #             raise TypeError('Multiple taxon sets in data source, but DataSet object is in attached (single) taxon set mode')
        #     for nxtaxa in taxon_set_elements:
        #         taxon_set = nxt.set_taxon_set_from_xml(nxtaxa, taxon_set=None)
        #         self.id_taxon_set_map[taxon_set.oid] = taxon_set
        # elif len(taxon_set_elements) == 1:
        #     if self.dataset is None:
        #         self.dataset = dendropy.DataSet()
        #     nxtaxa = taxon_set_elements[0]
        #     taxon_set = self.get_default_taxon_set(oid=nxtaxa.get('id', None), label=nxtaxa.get('label', None))
        #     taxon_set = nxt.set_taxon_set_from_xml(nxtaxa, taxon_set=taxon_set)
        #     self.id_taxon_set_map[taxon_set.oid] = taxon_set
        # else:
        #     raise error.DataParseError(message="No taxon definitions found in data source")
        taxon_set_elements = [i for i in xml_root.iter_otus()]
        if len(taxon_set_elements) == 0:
            raise error.DataParseError(message="No taxon definitions found in data source")
        if dataset.attached_taxon_set is not None:
            if len(taxon_set_elements) > 1:
                raise TypeError('Multiple taxon sets in data source, but DataSet object is in attached (single) taxon set mode')
            nxtaxa = taxon_set_elements[0]
            # taxon_set = self.get_default_taxon_set(oid=nxtaxa.get('id', None), label=nxtaxa.get('label', None))
            taxon_set = dataset.attached_taxon_set
            taxon_set = self.parse_taxon_set(nxtaxa, taxon_set=taxon_set, single_taxon_set_mode=True)
            self.id_taxon_set_map[taxon_set.oid] = taxon_set
            for t in taxon_set:
                self.id_taxon_map[(taxon_set.oid, t.oid)] = t
        else:
            for nxtaxa in taxon_set_elements:
                taxon_set = self.parse_taxon_set(nxtaxa, taxon_set=None)
                dataset.taxon_sets.append(taxon_set)
                self.id_taxon_set_map[taxon_set.oid] = taxon_set
                for t in taxon_set:
                    self.id_taxon_map[(taxon_set.oid, t.oid)] = t

    def parse_taxon_set(self, nxtaxa, taxon_set=None, single_taxon_set_mode=False):
        oid = nxtaxa.get('id', None)
        label = nxtaxa.get('label', None)
        if taxon_set is None:
            taxon_set = dendropy.TaxonSet(oid=oid, label=label)
        else:
            if oid is not None:
                taxon_set.oid = oid
            if label is not None:
                taxon_set.label = label
        annotations = [i for i in nxtaxa.findall_annotations()]
        for annotation in annotations:
            self.parse_annotations(taxon_set, annotation)
        if self.case_sensitive_taxon_labels:
            label_taxon_map = {}
        else:
            label_taxon_map = containers.OrderedCaselessDict()
        if single_taxon_set_mode:
            for t in taxon_set:
                label_taxon_map[t.label] = t
        for idx, nxtaxon in enumerate(nxtaxa.findall_otu()):
            taxon = None
            label = nxtaxon.get('label', "Taxon" + str(idx))
            oid = nxtaxon.get('id', "s" + str(idx) )
            if single_taxon_set_mode:
                taxon = label_taxon_map.get(label, None)
                if taxon is not None:
                    taxon.oid = oid
            if taxon is None:
                taxon = dendropy.Taxon(label=label, oid=oid)
                label_taxon_map[label] = taxon
                taxon_set.append(taxon)
            annotations = [i for i in nxtaxon.findall_annotations()]
            for annotation in annotations:
                self.parse_annotations(taxon, annotation)
        return taxon_set


    def parse_char_matrices(self, xml_root, dataset):
        """
        Given an xml_document, parses the XmlElement representation of
        character sequences into a list of CharacterMatrix objects.
        """
        nxc = _NexmlCharBlockParser(self.namespace_registry,
                self.id_taxon_set_map,
                self.id_taxon_map)
        for char_matrix_element in xml_root.iter_characters():
            nxc.parse_char_matrix(char_matrix_element, dataset)

    def parse_tree_lists(self, xml_root, dataset):
        """
        Given an xml_document object, parses the XmlElement structural
        representations of a set of NEXML treeblocks (`nex:trees`) and
        returns a TreeLists object corresponding to the NEXML.
        """
        nx_tree_parser = _NexmlTreesParser(self.namespace_registry,
                self.id_taxon_set_map,
                self.id_taxon_map)
        for trees_idx, trees_element in enumerate(xml_root.iter_trees()):
            for tree in nx_tree_parser.parse_trees(trees_element,
                    dataset,
                    trees_idx,
                    add_to_tree_list=True):
                pass

    def iterate_over_trees(file_obj, taxon_set=None, dataset=None):
        "Generator to iterate over trees in file without retaining any in memory."
        raise NotImplementedError
        # xml_doc = xmlparser.xml_document(file_obj=file_obj, namespace_list=SUPPORTED_NEXML_NAMESPACES)
        # if dataset is None:
        #     dataset = datasets.DataSet() or dataset
        # if taxon_set is None:
        #     taxon_set = taxa.TaxonSet()
        # if not (taxon_set in dataset.taxon_sets):
        #     dataset.taxon_sets.append(taxon_set)
        # self.parse_taxon_sets(xml_doc, dataset)
        # nx_tree_parser = _NexmlTreesParser(self.namespace_registry, self.id_taxon_set_map)
        # for trees_idx, trees_element in enumerate(xml_doc.getiterator('trees')):
        #     for tree in nx_tree_parser.parse_trees(trees_element, dataset, trees_idx, add_to_tree_list=False):
        #         yield tree

class _NexmlElementParser(_AnnotationParser):
    "Base parser class: wraps around annotations/dictionary element handling."

    def __init__(self, namespace_registry):
        _AnnotationParser.__init__(self, namespace_registry)

class _NexmlTreesParser(_NexmlElementParser):
    "Parses an XmlElement representation of NEXML schema tree blocks."

    def __init__(self, namespace_registry, id_taxon_set_map, id_taxon_map):
        _NexmlElementParser.__init__(self, namespace_registry)
        self.id_taxon_set_map = id_taxon_set_map
        self.id_taxon_map = id_taxon_map

    def parse_trees(self, nxtrees, dataset, trees_idx=None, add_to_tree_list=True):
        """
        Given an XmlElement object representing a NEXML treeblock,
        self.nxtrees (corresponding to a `nex:trees` element), this
        will construct and return a TreeList object defined by the
        underlying NEXML. If `add_to_tree_list` is False, then each tree,
        *IS NOT ADDED TO THE DATASET*.
        """
        oid = nxtrees.get('id', "Trees" + str(trees_idx))
        label = nxtrees.get('label', None)
        taxa_id = nxtrees.get('otus', None)
        if taxa_id is None:
            raise Exception("Taxa block not specified for trees block \"%s\"" % oid)
        # taxon_set = dataset.get_default_taxon_set(oid=taxa_id)
        taxon_set = self.id_taxon_set_map.get(taxa_id, None)
        if not taxon_set:
            raise Exception("Tree block '%s': Taxa block '%s' not found" % (oid, taxa_id))
        tree_list = dendropy.TreeList(oid=oid, label=label, taxon_set=taxon_set)
        dataset.add_tree_list(tree_list)
        annotations = [i for i in nxtrees.findall_annotations()]
        for annotation in annotations:
            self.parse_annotations(tree_list, annotation)
        tree_counter = 0
        for tree_element in nxtrees.findall_tree():
            tree_counter = tree_counter + 1
            oid = tree_element.get('id', tree_counter)
            label = tree_element.get('label', '')
            treeobj = dendropy.Tree (oid=oid, label=label)
            treeobj.taxon_set = taxon_set
            tree_type_attr = tree_element.get('{http://www.w3.org/2001/XMLSchema-instance}type')
            treeobj.length_type = _from_nexml_tree_length_type(tree_type_attr)
            annotations = [i for i in tree_element.findall_annotations()]
            for annotation in annotations:
                self.parse_annotations(treeobj, annotation)
            nodes = self.parse_nodes(tree_element, taxon_set=treeobj.taxon_set)
            edges = self.parse_edges(tree_element, length_type=treeobj.length_type)
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

            rootedge = self.parse_root_edge(tree_element, length_type=treeobj.length_type)
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

            if add_to_tree_list:
               tree_list.append(treeobj)

            yield treeobj

    def parse_nodes(self, tree_element, taxon_set):
        """
        Given an XmlElement representation of a NEXML tree element,
        (`nex:tree`) this will return a dictionary of DendroPy Node
        objects with the node_id as the key.
        """
        nodes = {}
        for nxnode in tree_element.findall_node():
            node_id = nxnode.get('id', None)
            nodes[node_id] = dendropy.Node()
            nodes[node_id].oid = node_id
            nodes[node_id].label = nxnode.get('label', None)
            taxon_id = nxnode.get('otu', None)
            if taxon_id is not None:
                try:
                    taxon = self.id_taxon_map[(taxon_set.oid, taxon_id)]
                except KeyError:
                    raise Exception('Taxon with id "%s" not defined in taxa block "%s"' % (taxon_id, taxon_set.oid))
                nodes[node_id].taxon = taxon
            annotations = [i for i in nxnode.findall_annotations()]
            for annotation in annotations:
                self.parse_annotations(nodes[node_id], annotation)
        return nodes

    def parse_root_edge(self, tree_element, length_type):
        "Returns the edge subtending the root node, or None if not defined."
        rootedge = tree_element.find_rootedge()
        if rootedge is not None:
            edge = dendropy.Edge()
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
            annotations = [i for i in rootedge.findall_annotations()]
            for annotation in annotations:
                self.parse_annotations(edge, annotation)
            return edge
        else:
            return None

    def parse_edges(self, tree_element, length_type):
        """
        Given an XmlElement representation of a NEXML tree element
        this will return a dictionary of DendroPy Edge objects.
        As at this stage, this method knows nothing about defined
        nodes, the Edge tail_node and head_node properties of the
        Edge are not set, but the tail_node_id and head_node_id are.
        """
        edges = {}
        edge_counter = 0
        for nxedge in tree_element.findall_edge():
            edge = dendropy.Edge()
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
            annotations = [i for i in nxedge.findall_annotations()]
            for annotation in annotations:
                self.parse_annotations(edge, annotation)
            edges[edge.oid] = edge
        return edges

class _NexmlCharBlockParser(_NexmlElementParser):
    "Parses an XmlElement representation of NEXML taxa blocks."

    def __init__(self, namespace_registry, id_taxon_set_map, id_taxon_map):
        _NexmlElementParser.__init__(self, namespace_registry)
        self.id_taxon_set_map = id_taxon_set_map
        self.id_taxon_map = id_taxon_map

    def parse_ambiguous_state(self, nxambiguous, state_alphabet):
        """
        Parses an XmlElement represent an ambiguous discrete character state,
        ("uncertain_state_set")
        and returns a corresponding StateAlphabetElement object.
        """
        state = dendropy.StateAlphabetElement(oid=nxambiguous.get('id', None),
                                                label=nxambiguous.get('label', None),
                                                symbol=nxambiguous.get('symbol', None),
                                                token=nxambiguous.get('token', None))
        state.member_states = []
        for nxmember in nxambiguous.findall_multistate_member():
            member_state_id = nxmember.get('state', None)
            member_state = state_alphabet.get_state('oid', member_state_id)
            state.member_states.append(member_state)
        state.multistate = dendropy.StateAlphabetElement.AMBIGUOUS_STATE
        return state

    def parse_polymorphic_state(self, nxpolymorphic, state_alphabet):
        """
        Parses an XmlElement represent a polymorphic discrete character state,
        ("polymorphic_state_set")
        and returns a corresponding StateAlphabetElement object.
        """
        state = dendropy.StateAlphabetElement(oid=nxpolymorphic.get('id', None),
                                                label=nxpolymorphic.get('label', None),
                                                symbol=nxpolymorphic.get('symbol', None),
                                                token=nxpolymorphic.get('token', None))
        state.member_states = []
        for nxmember in nxpolymorphic.findall_multistate_member():
            member_state_id = nxmember.get('state', None)
            member_state = state_alphabet.get_state('oid', member_state_id)
            state.member_states.append(member_state)
        for nxambiguous in nxpolymorphic.findall_uncertain_state_set():
            state.member_states.append(self.parse_ambiguous_state(nxambiguous, state_alphabet))
        state.multistate = dendropy.StateAlphabetElement.POLYMORPHIC_STATE
        return state

    def parse_state_alphabet(self, nxstates):
        """
        Given an XmlElement representing a nexml definition of (discrete or standard) states
        ("states"), this returns a corresponding StateAlphabet object.
        """

        state_alphabet = dendropy.StateAlphabet(oid=nxstates.get('id', None),
                                                         label=nxstates.get('label', None))
        for nxstate in nxstates.findall_char_state():
            state = dendropy.StateAlphabetElement(oid=nxstate.get('id', None),
                                                    label=nxstate.get('label', None),
                                                    symbol=nxstate.get('symbol', None),
                                                    token=nxstate.get('token', None))
            state_alphabet.append(state)
        for nxstate in nxstates.findall_uncertain_state_set():
            state_alphabet.append(self.parse_ambiguous_state(nxstate, state_alphabet))
        for nxstate in nxstates.findall_polymorphic_state_set():
            state_alphabet.append(self.parse_polymorphic_state(nxstate, state_alphabet))
        return state_alphabet

    def parse_characters_format(self, nxformat, char_matrix):
        """
        Given an XmlElement schema element ("format"), this parses the
        state definitions (if any) and characters (column definitions, if any),
        and populates the given char_matrix accordingly.
        """
        if nxformat is not None:
            for nxstates in nxformat.findall_char_states():
                char_matrix.state_alphabets.append(self.parse_state_alphabet(nxstates))
            for nxchars in nxformat.findall_char():
                col = dendropy.CharacterType(oid=nxchars.get('id', None))
                char_state_set_id = nxchars.get('states')
                if char_state_set_id is not None:
                    state_alphabet = None
                    for state_sets in char_matrix.state_alphabets:
                        if state_sets.oid == char_state_set_id:
                            state_alphabet = state_sets
                            break
                    if state_alphabet is None:
                        raise Exception("State set '%s' not defined" % char_state_set_id)
                    col.state_alphabet = state_alphabet
                elif hasattr(char_matrix, "default_state_alphabet") \
                    and char_matrix.default_state_alphabet is not None:
                    col.state_alphabet = char_matrix.default_state_alphabet
                char_matrix.character_types.append(col)

    def create_standard_character_alphabet(self, char_matrix, symbol_list=None):
        """
        Returns a standard character state alphabet based on symbol_list.
        Defaults to '0' - '9' if not specified.
        """
        if symbol_list is None:
            symbol_list = [str(i) for i in xrange(10)]
        state_alphabet = dendropy.StateAlphabet()
        for s in symbol_list:
            state_alphabet.append(dendropy.StateAlphabetElement(symbol=s))
        char_matrix.state_alphabets.append(state_alphabet)
        char_matrix.default_state_alphabet = state_alphabet

    def parse_char_matrix(self, nxchars, dataset):
        """
        Given an XmlElement representing a nexml characters block, this
        instantiates and returns a corresponding DendroPy CharacterMatrix object.
        """
        oid = nxchars.get('id', None)
        label = nxchars.get('label', None)
        nxchartype = nxchars.parse_type()
        if nxchartype.startswith('Dna'):
            char_matrix = dendropy.DnaCharacterMatrix()
            fixed_state_alphabet = True
        elif nxchartype.startswith('Rna'):
            char_matrix = dendropy.RnaCharacterMatrix()
            fixed_state_alphabet = True
        elif nxchartype.startswith('Protein'):
            char_matrix = dendropy.ProteinCharacterMatrix()
            fixed_state_alphabet = True
        elif nxchartype.startswith('Restriction'):
            char_matrix = dendropy.RestrictionSitesCharacterMatrix()
            fixed_state_alphabet = True
        elif nxchartype.startswith('Standard'):
            char_matrix = dendropy.StandardCharacterMatrix()
            fixed_state_alphabet = False
        elif nxchartype.startswith('Continuous'):
            char_matrix = dendropy.ContinuousCharacterMatrix()
            fixed_state_alphabet = False
        else:
            raise NotImplementedError('Character Block %s (\"%s\"): Character type "%s" not supported.'
                % (oid, label, nxchartype))

        char_matrix.oid = oid
        char_matrix.label = label

        taxa_id = nxchars.get('otus', None)
        if taxa_id is None:
            raise Exception("Character Block %s (\"%s\"): Taxa block not specified for trees block \"%s\"" % (char_matrix.oid, char_matrix.label, char_matrix.oid))
        # taxon_set = dataset.get_default_taxon_set(oid = taxa_id)
        taxon_set = self.id_taxon_set_map.get(taxa_id, None)
        if not taxon_set:
            raise Exception("Character Block %s (\"%s\"): Taxa block \"%s\" not found" % (char_matrix.oid, char_matrix.label, taxa_id))
        char_matrix.taxon_set = taxon_set
        annotations = [i for i in nxchars.findall_annotations()]
        for annotation in annotations:
            self.parse_annotations(char_matrix, annotation)

        nxformat = nxchars.find_char_format()
        if nxformat is not None:
            self.parse_characters_format(nxformat, char_matrix)
        elif isinstance(char_matrix, dendropy.StandardCharacterMatrix):
            # default to all integers < 10 as symbols
            self.create_standard_character_alphabet(char_matrix)

        matrix = nxchars.find_char_matrix()
        annotations = [i for i in matrix.findall_annotations()]
        for annotation in annotations:
            self.parse_annotations(char_matrix.taxon_seq_map, annotation)

        if char_matrix.character_types:
            id_chartype_map = char_matrix.id_chartype_map()
            chartypes = [char for char in char_matrix.character_types]
        else:
            id_chartype_map = {}
            chartypes = []
        for nxrow in matrix.findall_char_row():
            row_id = nxrow.get('id', None)

            label = nxrow.get('label', None)
            taxon_id = nxrow.get('otu', None)
            try:
                taxon = self.id_taxon_map[(taxon_set.oid, taxon_id)]
            except KeyError:
                raise error.DataParseError(message='Character Block %s (\"%s\"): Taxon with id "%s" not defined in taxa block "%s"' % (char_matrix.oid, char_matrix.label, taxon_id, taxon_set.oid))

            character_vector = dendropy.CharacterDataVector(oid=row_id, label=label, taxon=taxon)
            annotations = [i for i in nxrow.findall_annotations()]
            for annotation in annotations:
                self.parse_annotations(character_vector, annotation)

            if isinstance(char_matrix, dendropy.ContinuousCharacterMatrix):
                if nxchartype.endswith('Seqs'):
                    char_matrix.markup_as_sequences = True
                    seq = nxrow.find_char_seq()
                    if seq is not None:
                        seq = seq.replace('\n\r', ' ').replace('\r\n', ' ').replace('\n', ' ').replace('\r',' ')
                        col_idx = -1
                        for char in seq.split(' '):
                            char = char.strip()
                            if char:
                                col_idx += 1
                                if len(chartypes) <= col_idx:
                                    raise error.DataParseError(message="Character column/type ('<char>') not defined for character in position"\
                                        + " %d (matrix = '%s' row='%s', taxon='%s')" % (col_idx+1, char_matrix.oid, row_id, taxon.label))
                                cell = dendropy.CharacterDataCell(value=float(char), character_type=chartypes[col_idx])
                                character_vector.append(cell)
                else:
                    char_matrix.markup_as_sequences = False
                    for nxcell in nxrow.findall_char_cell():
                        chartype_id = nxcell.get('char', None)
                        if chartype_id is None:
                            raise error.DataParseError(message="'char' attribute missing for cell: cell markup must indicate character column type for character"\
                                        + " (matrix = '%s' row='%s', taxon='%s')" % (char_matrix.oid, row_id, taxon.label))
                        if chartype_id not in id_chartype_map:
                            raise error.DataParseError(message="Character type ('<char>') with id '%s' referenced but not found for character" % chartype_id \
                                        + " (matrix = '%s' row='%s', taxon='%s')" % (char_matrix.oid, row_id, taxon.label))
                        chartype = id_chartype_map[chartype_id]
                        pos_idx = chartypes.index(chartype)
#                         column = id_chartype_map[chartype_id]
#                         state = column.state_id_map[cell.get('state', None)]
                        cell = dendropy.CharacterDataCell(value=float(nxcell.get('state')), character_type=chartype)
                        annotations = [i for i in nxtaxon.findall_annotations]
                        for annotation in annotations:
                            self.parse_annotations(cell, annotation)
                        character_vector.set_cell_by_index(pos_idx, cell)
            else:
                if nxchartype.endswith('Seqs'):
                    char_matrix.markup_as_sequences = True
#                     symbol_state_map = char_matrix.default_state_alphabet.symbol_state_map()
                    seq = nxrow.find_char_seq()
                    if seq is not None:
                        seq = seq.replace(' ', '').replace('\n', '').replace('\r', '')
                        col_idx = -1
                        for char in seq:
                            col_idx += 1
                            symbol_state_map = char_matrix.character_types[col_idx].state_alphabet.symbol_state_map()
                            if char in symbol_state_map:
                                if len(chartypes) <= col_idx:
                                    raise error.DataParseError(message="Character column/type ('<char>') not defined for character in position"\
                                        + " %d (matrix = '%s' row='%s', taxon='%s')" % (col_idx+1, char_matrix.oid, row_id, taxon.label))
                                state = symbol_state_map[char]
                                character_type = chartypes[col_idx]
                                character_vector.append(dendropy.CharacterDataCell(value=state, character_type=character_type))
                            else:
                                raise error.DataParseError(message="Character Block '%s', row '%s', character position %s: State with symbol '%s' in sequence '%s' not defined" \
                                        % (char_matrix.oid, row_id, col_idx, char, seq))
                else:
                    char_matrix.markup_as_sequences = False
                    id_state_maps = {}
                    for nxcell in nxrow.findall_char_cell():
                        chartype_id = nxcell.get('char', None)
                        if chartype_id is None:
                            raise error.DataParseError(message="'char' attribute missing for cell: cell markup must indicate character column type for character"\
                                        + " (matrix = '%s' row='%s', taxon='%s')" % (char_matrix.oid, row_id, taxon.label))
                        if chartype_id not in id_chartype_map:
                            raise error.DataParseError(message="Character type ('<char>') with id '%s' referenced but not found for character" % chartype_id \
                                        + " (matrix = '%s' row='%s', taxon='%s')" % (char_matrix.oid, row_id, taxon.label))
                        chartype = id_chartype_map[chartype_id]
                        pos_idx = chartypes.index(chartype)
                        if chartype_id not in id_state_maps:
                            id_state_maps[chartype_id] = chartype.state_alphabet.id_state_map()
                        state = id_state_maps[chartype_id][nxcell.get('state')]
                        cell = dendropy.CharacterDataCell(value=state, character_type=chartype)
                        character_vector.set_cell_by_index(pos_idx, cell)

            char_matrix[taxon] = character_vector

        if fixed_state_alphabet:
            char_matrix.remap_to_default_state_alphabet_by_symbol(purge_other_state_alphabets=True)
        dataset.char_matrices.append(char_matrix)

class NexmlWriter(iosys.DataWriter):
    "Implements the DataWriter interface for handling NEXML files."

    def __init__(self, **kwargs):
        "Calls the base class constructor."
        iosys.DataWriter.__init__(self, **kwargs)
        self.indent = "    "
        self._prefix_uri_tuples = set()
        self.markup_as_sequences = kwargs.get("markup_as_sequences", None)

    def write(self, stream):
        """
        Writes a list of DendroPy Tree objects to a full NEXML
        document.
        """
        body = StringIO()
        self.write_annotations(self.dataset, body, indent_level=1)
        self.write_taxon_sets(taxon_sets=self.dataset.taxon_sets, dest=body)
        if not self.exclude_chars:
            self.write_char_matrices(char_matrices=self.dataset.char_matrices, dest=body)
        if not self.exclude_trees:
            self.write_tree_lists(tree_lists=self.dataset.tree_lists, dest=body)

        self.write_to_nexml_open(stream, indent_level=0)
        stream.write(body.getvalue())
        self.write_to_nexml_close(stream, indent_level=0)

    ### class-specific  ###

    def write_taxon_sets(self, taxon_sets, dest, indent_level=1):
        "Writes out TaxonSets."
        for idx, taxon_set in enumerate(taxon_sets):
            dest.write(self.indent * indent_level)
            parts = []
            parts.append('otus')
            parts.append('id="%s"' % taxon_set.default_oid)
            if taxon_set.label:
                parts.append('label=%s' % _protect_attr(taxon_set.label))
            dest.write("<%s>\n" % ' '.join(parts))

            # annotate
#             self.write_extensions(taxon_set, dest, indent_level=indent_level+1)
            self.write_annotations(taxon_set, dest, indent_level=indent_level+1)
            for taxon in taxon_set:
                dest.write(self.indent * (indent_level+1))
                parts = []
                parts.append('otu')
                parts.append('id="%s"' % taxon.default_oid)
                if taxon.label:
                    parts.append('label=%s' % _protect_attr(taxon.label))
                if isinstance(taxon, dendropy.AnnotatedDataObject) and len(taxon.annotations) > 0:
                    dest.write("<%s>\n" % ' '.join(parts))
                    # self.write_extensions(taxon, dest, indent_level=indent_level+2)
                    self.write_annotations(taxon, dest, indent_level=indent_level+2)
                    dest.write(self.indent * (indent_level+1))
                    dest.write("</otu>\n")
                else:
                    dest.write("<%s />\n" % ' '.join(parts))
            dest.write(self.indent * indent_level)
            dest.write('</otus>\n')

    def write_tree_lists(self, tree_lists, dest, indent_level=1):
        "Writes out TreeLists."
        for idx, tree_list in enumerate(tree_lists):
            dest.write(self.indent * indent_level)
            parts = []
            parts.append('trees')
            parts.append('id="%s"' % tree_list.default_oid)
            if tree_list.label:
                parts.append('label=%s' % _protect_attr(tree_list.label))
            parts.append('otus="%s"' % tree_list.taxon_set.default_oid)
            dest.write("<%s>\n" % ' '.join(parts))

            # annotate
#             self.write_extensions(tree_list, dest, indent_level=indent_level+1)
            self.write_annotations(tree_list, dest, indent_level=indent_level+1)

            for tree in tree_list:
                self.write_tree(tree=tree, dest=dest, indent_level=2)
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
            parts.append('otus="%s"' % char_matrix.taxon_set.default_oid)
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


                if ( (self.markup_as_sequences is not None and self.markup_as_sequences is False)
                        or (hasattr(char_matrix, 'markup_as_sequences') and not char_matrix.markup_as_sequences)
                        ):
                    for cell in row:
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

    def write_tree(self, tree, dest, indent_level=0):
        """
        Writes a single DendroPy Tree object as a NEXML nex:tree
        element.
        """
        parts = []
        parts.append('tree')
        if hasattr(tree, 'oid') and tree.default_oid is not None:
            parts.append('id="%s"' % tree.default_oid)
        else:
            parts.append('id="%s"' % ("Tree" + str(id(tree))))
        if hasattr(tree, 'label') and tree.label:
            parts.append('label=%s' % _protect_attr(tree.label))
        if hasattr(tree, 'length_type') and tree.length_type:
            parts.append('xsi:type="%s"' % _to_nexml_tree_length_type(tree.length_type))
        else:
            parts.append('xsi:type="nex:FloatTree"')
        parts = ' '.join(parts)
        dest.write('%s<%s>\n'
                   % (self.indent * indent_level, parts))
        # annotate
#         self.write_extensions(tree, dest, indent_level=indent_level+1)
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

    def write_node(self, node, dest, indent_level=0):
        "Writes out a NEXML node element."
        parts = []
        parts.append('<node')
        parts.append('id="%s"' % node.default_oid)
        if hasattr(node, 'label') and node.label:
            parts.append('label=%s' % _protect_attr(node.label))
        if hasattr(node, 'taxon') and node.taxon:
            parts.append('otu="%s"' % node.taxon.default_oid)
        parts = ' '.join(parts)
        dest.write('%s%s' % ((self.indent * indent_level), parts))
        if len(node.annotations) > 0:
            dest.write('>\n')
#             self.write_extensions(node, dest, indent_level=indent_level+1)
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
            if hasattr(edge, 'oid') and edge.default_oid:
                parts.append('id="%s"' % edge.default_oid)
            # programmatically more efficent to do this in above
            # block, but want to maintain this tag order ...
            if edge.tail_node is not None:
                parts.append('source="%s"' % edge.tail_node.default_oid)
            if edge.head_node is not None:
                parts.append('target="%s"' % edge.head_node.default_oid)
            if hasattr(edge, 'length') and edge.length is not None:
                parts.append('length="%s"' % edge.length)
            if hasattr(edge, 'label') and edge.label:
                parts.append('label=%s' % _protect_attr(edge.label))

            # only write if we have more than just the 'edge' and '/' bit
            if len(parts) > 2:
                parts = ' '.join(parts)
                dest.write('%s%s' % ((self.indent * indent_level), parts))
                if len(edge.annotations) > 0:
                    dest.write('>\n')
                    # self.write_extensions(edge, dest, indent_level=indent_level+1)
                    self.write_annotations(edge, dest,
                                           indent_level=indent_level+1)
                    dest.write('%s</%s>\n' % ((self.indent * indent_level), tag))
                else:
                    dest.write(' />\n')

    def write_annotations(self, annotated, dest, indent_level=0):
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
