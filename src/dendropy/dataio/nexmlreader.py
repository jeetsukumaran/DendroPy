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
Deserialization of NeXML-formatted data.
"""

import re
import collections
from dendropy.dataio import ioservice
from dendropy.utility import container
from dendropy.utility import textprocessing
from dendropy.utility import error
from dendropy.dataio import xmlprocessing
from dendropy.datamodel.charstatemodel import StateAlphabet, StateIdentity

SUPPORTED_NEXML_NAMESPACES = ('http://www.nexml.org/1.0', 'http://www.nexml.org/2009')

def _from_nexml_tree_length_type(type_attr):
    """
    Given an attribute string read from a nexml tree element, returns
    the python class of the edge length attribute.
    """
    if type_attr == "nex:IntTree":
        return int
    else:
        return float

class _AnnotationParser(object):

    def __init__(self, namespace_registry=None):
        self._namespace_registry = namespace_registry

    def _parse_annotations(self, annotated, nxelement):
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
        name_prefix, name = textprocessing.parse_curie_standard_qualified_name(key)
        try:
            namespace = self._namespace_registry.prefix_namespace_map[name_prefix]
        except KeyError:
            raise ValueError("CURIE-standard prefix '%s' not defined in document: %s" % (name_prefix, self._namespace_registry))
        if datatype_hint is not None:
            dt_prefix, dt = textprocessing.parse_curie_standard_qualified_name(datatype_hint)
            if dt_prefix is not None:
                try:
                    dt_namespace = self._namespace_registry.prefix_namespace_map[dt_prefix]
                except KeyError:
                    raise ValueError("CURIE-standard prefix '%s' not defined in document: %s" % (dt_prefix, self._namespace_registry))
                if dt_namespace.startswith("http://www.w3.org/2001/XMLSchema"):
                    value = self._coerce_to_xml_schema_type(value, dt)
                elif dt_namespace.startswith("http://www.nexml.org/1.0") or dt_namespace.startswith("http://www.nexml.org/2009"):
                    value = self._coerce_to_nexml_type(value, dt)
                elif dt_namespace.startswith("http://dendropy.org") or dt_namespace.startswith("http://packages.python.org/DendroPy") or dt_namespace.startswith("http://pypi.org/project/DendroPy/"):
                    value = self._coerce_to_dendropy_type(value, dt)
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
            self._parse_annotations(a, annotation)

    def _coerce_to_xml_schema_type(self, value, type_name):
        if type_name in ("boolean"):
            if value.lower() in ("1", "t", "true", "y", "yes",):
                return True
            else:
                return False
        elif type_name in ("decimal", "float", "double"):
            coerce_type = float
        elif type_name in ("byte", "int", "integer", "long", "negativeInteger", "nonNegativeInteger", "nonPositiveInteger", "short", "unsignedInt", "unsignedLong", "unsignedShort"):
            coerce_type = int
        else:
            return value
        return self._coerce_type(value, coerce_type)

    def _coerce_to_dendropy_type(self, value, type_name):
        if type_name in ("decimalRange"):
            try:
                value = [float(v) for v in value.split()]
            except KeyError:
                pass
        return value

    def _coerce_to_nexml_type(self, value, type_name):
        if type_name in ("ContinuousSeq"):
            try:
                value = [float(v) for v in value.split()]
            except KeyError:
                pass
        return value

    def _coerce_type(self, value, to_type):
        try:
            value = to_type(value)
        except ValueError:
            pass
        return value


class NexmlElement(xmlprocessing.XmlElement):

    def __init__(self, element, default_namespace=None):
        # if default_namespace is None:
        #     default_namespace = NexmlReader.DEFAULT_NEXML_NAMESPACE
        # else:
        #     default_namespace = default_namespace
        xmlprocessing.XmlElement.__init__(self,
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

class NexmlReader(ioservice.DataReader, _AnnotationParser):
    "Implements thinterface for handling NEXML files."

    ###########################################################################
    ## Life-cycle and Setup

    DEFAULT_NEXML_NAMESPACE = "http://www.nexml.org/2009"

    def __init__(self, **kwargs):
        """
        Keyword Arguments
        -----------------

        default_namespace : str
            Default namespace to use for elements.
        case_sensitive_taxon_labels: boolean, default: |False|
            If |True|, then case is respected when matching taxon names.
            Default is |False|: case is ignored.
        ignore_unrecognized_keyword_arguments : boolean, default: |False|
            If |True|, then unsupported or unrecognized keyword arguments will
            not result in an error. Default is |False|: unsupported keyword
            arguments will result in an error.

        """
        _AnnotationParser.__init__(self)
        ioservice.DataReader.__init__(self)
        self.default_namespace = kwargs.pop("default_namespace", NexmlReader.DEFAULT_NEXML_NAMESPACE)
        self.case_sensitive_taxon_labels = kwargs.pop('case_sensitive_taxon_labels', False)

        ### NOTE: following are not actually supported.
        ### They are here because some tests automatically add include them on calls.
        ### TODO: remove these keywords from generic tests
        self.suppress_internal_node_taxa = kwargs.pop("suppress_internal_node_taxa", True)
        self.suppress_leaf_node_taxa = kwargs.pop("suppress_external_node_taxa", False) # legacy (will be deprecated)
        self.suppress_leaf_node_taxa = kwargs.pop("suppress_leaf_node_taxa", self.suppress_leaf_node_taxa)
        self.data_type = kwargs.pop("data_type", None)

        self.check_for_unused_keyword_arguments(kwargs)

        # Set up parsing meta-variables
        self._id_taxon_namespace_map = {}
        self._id_taxon_map = {}
        self._taxon_namespace_factory = None
        self._tree_list_factory = None
        self._char_matrix_factory = None
        self._state_alphabet_factory = None
        self._global_annotations_target = None
        self._taxon_namespaces = []
        self._char_matrices = []
        self._tree_lists = []
        self._product = None

    ###########################################################################
    ## Reader Interface

    def _read(self,
            stream,
            taxon_namespace_factory=None,
            tree_list_factory=None,
            char_matrix_factory=None,
            state_alphabet_factory=None,
            global_annotations_target=None):
        xml_doc = xmlprocessing.XmlDocument(file_obj=stream,
                subelement_factory=self._subelement_factory)
        self._namespace_registry = xml_doc.namespace_registry
        self._taxon_namespace_factory = taxon_namespace_factory
        self._tree_list_factory = tree_list_factory
        self._char_matrix_factory = char_matrix_factory
        self._state_alphabet_factory = state_alphabet_factory
        self._global_annotations_target = global_annotations_target
        self._parse_document(xml_doc)
        self._product = self.Product(
                taxon_namespaces=self._taxon_namespaces,
                tree_lists=self._tree_lists,
                char_matrices=self._char_matrices)
        return self._product

    ###########################################################################
    ## Support

    def _subelement_factory(self, element):
        return NexmlElement(element, default_namespace=self.default_namespace)

    ###########################################################################
    ## Data Management

    def _new_taxon_namespace(self, label=None):
        if self.attached_taxon_namespace is not None:
            return self.attached_taxon_namespace
        taxon_namespace = self._taxon_namespace_factory(label=label)
        self._taxon_namespaces.append(taxon_namespace)
        return taxon_namespace

    def _new_char_matrix(self, data_type, taxon_namespace, label=None, **kwargs):
        char_matrix = self._char_matrix_factory(
                data_type,
                taxon_namespace=taxon_namespace,
                label=label,
                **kwargs)
        self._char_matrices.append(char_matrix)
        return char_matrix

    def _new_state_alphabet(self, *args, **kwargs):
        return self._state_alphabet_factory(*args, **kwargs)

    def _new_tree_list(self, taxon_namespace, label=None):
        tree_list = self._tree_list_factory(
                taxon_namespace=taxon_namespace,
                label=label)
        self._tree_lists.append(tree_list)
        return tree_list

    ###########################################################################
    ## Parsers

    ## Following methods are class-specific ###

    def _parse_document(self, xml_doc):
        xml_root = xml_doc.root
        self._parse_taxon_namespaces(xml_root)
        if self._char_matrix_factory is not None:
            self._parse_char_matrices(xml_root)
        if self._tree_list_factory is not None:
            self._parse_tree_lists(xml_root)
        if self._global_annotations_target is not None:
            top_annotations = [i for i in xml_root.findall_annotations()]
            for annotation in top_annotations:
                self._parse_annotations(self._global_annotations_target, annotation)

    def _parse_taxon_namespaces(self, xml_root):
        for nxtaxa in xml_root.iter_otus():
            taxon_namespace_label = nxtaxa.get('label', None)
            taxon_namespace = self._new_taxon_namespace(label=taxon_namespace_label)
            taxon_namespace_id = nxtaxa.get('id', id(taxon_namespace))
            self._id_taxon_namespace_map[taxon_namespace_id] = taxon_namespace
            annotations = [i for i in nxtaxa.findall_annotations()]
            for annotation in annotations:
                self._parse_annotations(taxon_namespace, annotation)
            if self.case_sensitive_taxon_labels:
                label_taxon_map = {}
            else:
                label_taxon_map = container.OrderedCaselessDict()
            if self.attached_taxon_namespace is not None:
                for t in taxon_namespace:
                    label_taxon_map[t.label] = t
            for idx, nxtaxon in enumerate(nxtaxa.findall_otu()):
                taxon = None
                taxon_label = nxtaxon.get('label', None)
                taxon_oid = nxtaxon.get('id', id(nxtaxon))
                if taxon_label is not None and self.attached_taxon_namespace is not None:
                    # taxon = label_taxon_map.get_taxon(
                    #         label=taxon_label,
                    #         case_sensitive=self.case_sensitive_taxon_labels)
                    try:
                        taxon = label_taxon_map[taxon_label]
                    except KeyError:
                        taxon = None
                if taxon is None:
                    taxon = taxon_namespace.new_taxon(label=taxon_label)
                annotations = [i for i in nxtaxon.findall_annotations()]
                for annotation in annotations:
                    self._parse_annotations(taxon, annotation)
                self._id_taxon_map[(taxon_namespace_id, taxon_oid)] = taxon

    def _parse_char_matrices(self, xml_root):
        nxc = _NexmlCharBlockParser(self._namespace_registry,
                self._id_taxon_namespace_map,
                self._id_taxon_map,
                self._new_char_matrix,
                self._state_alphabet_factory)
        for char_matrix_element in xml_root.iter_characters():
            nxc.parse_char_matrix(char_matrix_element)

    def _parse_tree_lists(self, xml_root):
        for trees_idx, trees_element in enumerate(xml_root.iter_trees()):
            self._parse_tree_list(trees_element, trees_idx, True)

    def _parse_tree_list(self, nxtrees, trees_idx=None, add_to_tree_list=True):
        trees_id = nxtrees.get('id', "Trees" + str(trees_idx))
        trees_label = nxtrees.get('label', None)
        otus_id = nxtrees.get('otus', None)
        if otus_id is None:
            raise Exception("Taxa block not specified for trees block '{}'".format(otus_id))
        taxon_namespace = self._id_taxon_namespace_map.get(otus_id, None)
        if not taxon_namespace:
            raise Exception("Tree block '{}': Taxa block '{}' not found".format(trees_id, otus_id))
        tree_list = self._new_tree_list(
                label=trees_label,
                taxon_namespace=taxon_namespace)
        annotations = [i for i in nxtrees.findall_annotations()]
        for annotation in annotations:
            self._parse_annotations(tree_list, annotation)
        tree_parser = _NexmlTreeParser(
                id_taxon_map=self._id_taxon_map,
                annotations_processor_fn=self._parse_annotations,
                )
        for tree_element in nxtrees.findall_tree():
            tree_obj = tree_list.new_tree()
            tree_parser.build_tree(tree_obj, tree_element, otus_id)

class _NexmlTreeParser(object):

    EdgeInfo = collections.namedtuple(
            "EdgeInfo",
            ["edge_id", "label", "tail_node_id", "head_node_id", "length", "annotations"]
            )

    def __init__(self,
            id_taxon_map,
            annotations_processor_fn,
            ):
        self._id_taxon_map = id_taxon_map
        self._process_annotations = annotations_processor_fn

    def build_tree(self, tree_obj, tree_element, otus_id):
        tree_obj.label = tree_element.get('label', '')
        tree_type_attr = tree_element.get('{http://www.w3.org/2001/XMLSchema-instance}type')
        tree_obj.length_type = _from_nexml_tree_length_type(tree_type_attr)
        tree_obj.is_rooted = False # default unless explicitly set
        annotations = [i for i in tree_element.findall_annotations()]
        for annotation in annotations:
            self._process_annotations(tree_obj, annotation)
        unparented_node_set, node_id_map, root_node_list = self._parse_nodes(
                tree_element,
                tree_obj.node_factory,
                otus_id)
        if len(root_node_list) > 1:
            raise Exception("Multiple root nodes defined but the DendroPy tree model currently only supported single-root trees")
        elif len(root_node_list) == 1:
            tree_obj.seed_node = root_node_list[0]
            tree_obj.is_rooted = True
        edges_info = self._parse_edges_info(tree_element, length_type=tree_obj.length_type)
        for edge_info in edges_info:
            if edge_info.head_node_id is None:
                raise Exception("Edge '{}' does not specify a head (target) node".format(edge_info.edge_id, edge_info.head_node_id))
            try:
                head_node = node_id_map[edge_info.head_node_id]
            except KeyError:
                raise Exception("Edge '{}' specifies a non-defined head (target) node '{}'".format(edge_info.edge_id, edge_info.head_node_id))
            if edge_info.tail_node_id is not None:
                try:
                    tail_node = node_id_map[edge_info.tail_node_id]
                except KeyError:
                    raise Exception("Edge '{}' specifies a non-defined tail (source) node '{}'".format(edge_info.edge_id, edge_info.tail_node_id))
                head_node.parent_node = tail_node
                assert head_node.edge.tail_node is tail_node
                unparented_node_set.remove(head_node)
            else:
                assert head_node.parent_node is None
            self._set_edge_details(head_node.edge, edge_info)

        # find node(s) without parent
        # unparented_node_sets = []
        # for node in nodes.values():
        #     if node.parent_node == None:
        #         unparented_node_sets.append(node)

        # If one unparented_node_sets node found, this is the root: we use
        # it as the tree head node. If multiple unparented_node_sets nodes
        # are found, then we add them all as child_nodes of the
        # existing head node. If none, then we have some sort of
        # cyclicity, and we are not dealing with a tree.
        if len(unparented_node_set) == 1:
            unparented_node = unparented_node_set.pop()
            if root_node_list:
                if root_node_list[0] is not unparented_node:
                    raise Exception("Tree already has an explictly defined root node, but node without parent found: {}".format(unparented_node))
            else:
                tree_obj.seed_node = unparented_node
        elif len(unparented_node_set) > 1:
            for node in unparented_node_set:
                tree_obj.seed_node.add_child(node)
        else:
            raise Exception("Structural error: tree must be acyclic.")

        root_edge_element = tree_element.find_rootedge()
        if root_edge_element is not None:
            root_edge_info = self._parse_edge_info(root_edge_element, tree_obj.length_type)
            root_node = node_id_map[root_edge_info.head_node_id]
            if root_node is not tree_obj.seed_node:
                raise Exception("Root edge does not subtend root node")
            root_edge = root_node.edge
            self._set_edge_details(root_edge, root_edge_info)
            # tree_obj.is_rooted = True

        return tree_obj

    def _parse_nodes(self, tree_element, node_factory, otus_id):
        """
        Given an XmlElement representation of a NEXML tree element,
        (`nex:tree`) this will return a dictionary of DendroPy Node
        objects with the node_id as the key.
        """
        node_set = set()
        node_id_map = {}
        root_node_list = []
        for nxnode in tree_element.findall_node():
            node_id = nxnode.get('id', None)
            node = node_factory()
            node_id_map[node_id] = node
            node_set.add(node)
            node_id_map[node_id].label = nxnode.get('label', None)
            taxon_id = nxnode.get('otu', None)
            if taxon_id is not None:
                try:
                    taxon = self._id_taxon_map[(otus_id, taxon_id)]
                except KeyError:
                    raise Exception("Taxon with id '{}' not defined '{}'".format(taxon_id, otus_id))
                node_id_map[node_id].taxon = taxon
            annotations = [i for i in nxnode.findall_annotations()]
            for annotation in annotations:
                self._process_annotations(node_id_map[node_id], annotation)
            rooting_state = nxnode.get('root', None)
            if rooting_state is not None and rooting_state.lower() in ("1", "t", "true"):
                root_node_list.append(node)
        return node_set, node_id_map, root_node_list

    def _parse_edges_info(self, tree_element, length_type):
        edges_info = []
        for nxedge in tree_element.findall_edge():
            edges_info.append(self._parse_edge_info(nxedge, length_type))
        return edges_info

    def _parse_edge_info(self, nxedge, length_type):
        edge_id = nxedge.get('id', None)
        edge_label = nxedge.get('label', None)
        tail_node_id = nxedge.get('source', None)
        # assert tail_node_id is not None
        head_node_id = nxedge.get('target', None)
        # assert head_node_id is not None
        edge_length_str = nxedge.get('length', 0.0)
        edge_length = None
        try:
            edge_length = length_type(edge_length_str)
        except ValueError:
            msg = "Edge {} 'length' attribute is not of type {}: '{}'".format(
                    edge_id, str(length_type), edge_length_str)
            raise Exception(msg)
        annotations = [i for i in nxedge.findall_annotations()]
        e = _NexmlTreeParser.EdgeInfo(
                edge_id=edge_id,
                label=edge_label,
                tail_node_id=tail_node_id,
                head_node_id=head_node_id,
                length=edge_length,
                annotations=annotations)
        return e

    def _set_edge_details(self, edge, edge_info):
        edge.length = edge_info.length
        edge.label = edge_info.label
        for annotation in edge_info.annotations:
            self._process_annotations(edge, annotation)

class _NexmlCharBlockParser(_AnnotationParser):
    "Parses an XmlElement representation of NEXML taxa blocks."

    def __init__(self,
            namespace_registry,
            id_taxon_namespace_map,
            id_taxon_map,
            char_matrix_factory,
            state_alphabet_factory,
            ):
        _AnnotationParser.__init__(self, namespace_registry)
        self._id_taxon_namespace_map = id_taxon_namespace_map
        self._id_taxon_map = id_taxon_map
        self._char_matrix_factory = char_matrix_factory
        self._state_alphabet_factory = state_alphabet_factory

        self._id_state_alphabet_map = {}
        self._id_state_map = {}
        self._id_chartype_map = {}
        self._char_types = []
        self._chartype_id_to_pos_map = {}

    def parse_char_matrix(self, nxchars):
        """
        Given an XmlElement representing a nexml characters block, this
        instantiates and returns a corresponding DendroPy CharacterMatrix object.
        """

        # clear
        self._id_state_alphabet_map = {}
        self._id_state_map = {}
        self._id_chartype_map = {}
        self._char_types = []
        self._chartype_id_to_pos_map = {}

        # initiaiize
        label = nxchars.get('label', None)
        char_matrix_oid = nxchars.get('oid', '')

        # set up taxa
        otus_id = nxchars.get('otus', None)
        if otus_id is None:
            raise Exception("Character Block %s (\"%s\"): Taxon namespace not specified" % (char_matrix_oid, label))
        taxon_namespace = self._id_taxon_namespace_map.get(otus_id, None)
        if not taxon_namespace:
            raise Exception("Character Block %s (\"%s\"): Specified taxon namespace not found" % (char_matrix_oid, label))

        # character matrix instantiation
        nxchartype = nxchars.parse_type()
        extra_kwargs = {}
        if nxchartype.startswith('Dna'):
            data_type = "dna"
        elif nxchartype.startswith('Rna'):
            data_type = "rna"
        elif nxchartype.startswith('Protein'):
            data_type = "protein"
        elif nxchartype.startswith('Restriction'):
            data_type = "restriction"
        elif nxchartype.startswith('Standard'):
            data_type = "standard"
            extra_kwargs["default_state_alphabet"] = None
        elif nxchartype.startswith('Continuous'):
            data_type = "continuous"
        else:
            raise Exception("Character Block %s (\"%s\"): Character type '%s' not supported" % (char_matrix_oid, label, nxchartype))
        char_matrix = self._char_matrix_factory(
                data_type,
                taxon_namespace=taxon_namespace,
                label=label,
                **extra_kwargs)

        # annotation processing
        annotations = [i for i in nxchars.findall_annotations()]
        for annotation in annotations:
            self._parse_annotations(char_matrix, annotation)

        # get state mappings
        nxformat = nxchars.find_char_format()
        if nxformat is not None:
            self.parse_characters_format(nxformat, data_type, char_matrix)
        elif data_type == "standard":
            self.create_standard_character_alphabet(char_matrix)

        nxmatrix = nxchars.find_char_matrix()
        annotations = [i for i in nxmatrix.findall_annotations()]
        for annotation in annotations:
            self._parse_annotations(char_matrix.taxon_seq_map, annotation)
        for nxrow in nxmatrix.findall_char_row():
            row_id = nxrow.get('id', None)
            label = nxrow.get('label', None)
            taxon_id = nxrow.get('otu', None)
            try:
                taxon = self._id_taxon_map[(otus_id, taxon_id)]
            except KeyError:
                raise error.DataParseError(message='Character Block %s (\"%s\"): Taxon with id "%s" not defined in taxa block "%s"' % (char_matrix.oid, char_matrix.label, taxon_id, otus_id))

            character_vector = char_matrix.new_sequence(taxon=taxon)
            annotations = [i for i in nxrow.findall_annotations()]
            for annotation in annotations:
                self._parse_annotations(character_vector, annotation)

            if data_type == "continuous":
                if nxchartype.endswith('Seqs'):
                    seq = nxrow.find_char_seq()
                    if seq is not None:
                        seq = seq.replace('\n\r', ' ').replace('\r\n', ' ').replace('\n', ' ').replace('\r',' ')
                        col_idx = -1
                        for char in seq.split(' '):
                            char = char.strip()
                            if char:
                                col_idx += 1
                                if len(self._char_types) <= col_idx:
                                    raise error.DataParseError(message="Character column/type ('<char>') not defined for character in position"\
                                        + " %d (matrix = '%s' row='%s', taxon='%s')" % (col_idx+1, char_matrix.oid, row_id, taxon.label))
                                character_vector.append(character_value=float(char), character_type=self._char_types[col_idx])
                else:
                    for nxcell in nxrow.findall_char_cell():
                        chartype_id = nxcell.get('char', None)
                        if chartype_id is None:
                            raise error.DataParseError(message="'char' attribute missing for cell: cell markup must indicate character column type for character"\
                                        + " (matrix = '%s' row='%s', taxon='%s')" % (char_matrix.oid, row_id, taxon.label))
                        if chartype_id not in self._id_chartype_map:
                            raise error.DataParseError(message="Character type ('<char>') with id '%s' referenced but not found for character" % chartype_id \
                                        + " (matrix = '%s' row='%s', taxon='%s')" % (char_matrix.oid, row_id, taxon.label))
                        chartype = self._id_chartype_map[chartype_id]
                        pos_idx = self._char_types.index(chartype)
#                         column = id_chartype_map[chartype_id]
#                         state = column.state_id_map[cell.get('state', None)]
                        # annotations = [i for i in nxcell.findall_annotations]
                        # for annotation in annotations:
                        #     self._parse_annotations(cell, annotation)
                        character_vector.append(character_value=float(nxcell.get('state')),
                                character_type=chartype)
            else:
                if nxchartype.endswith('Seqs'):
                    seq = nxrow.find_char_seq()
                    if seq is not None:
                        seq = seq.replace(' ', '').replace('\n', '').replace('\r', '')
                        col_idx = -1
                        for char in seq:
                            col_idx += 1
                            state_alphabet = char_matrix.character_types[col_idx].state_alphabet
                            try:
                                state = state_alphabet[char]
                            except KeyError:
                                raise error.DataParseError(message="Character Block row '%s', character position %s: State with symbol '%s' in sequence '%s' not defined" \
                                        % (row_id, col_idx, char, seq))
                            if len(self._char_types) <= col_idx:
                                raise error.DataParseError(message="Character column/type ('<char>') not defined for character in position"\
                                    + " %d (row='%s', taxon='%s')" % (col_idx+1, row_id, taxon.label))
                            character_type = self._char_types[col_idx]
                            character_vector.append(character_value=state,
                                    character_type=character_type)
                else:
                    for nxcell in nxrow.findall_char_cell():
                        chartype_id = nxcell.get('char', None)
                        if chartype_id is None:
                            raise error.DataParseError(message="'char' attribute missing for cell: cell markup must indicate character column type for character"\
                                        + " (matrix = '%s' row='%s', taxon='%s')" % (char_matrix_oid, row_id, taxon.label))
                        if chartype_id not in self._id_chartype_map:
                            raise error.DataParseError(message="Character type ('<char>') with id '%s' referenced but not found for character" % chartype_id \
                                        + " (matrix = '%s' row='%s', taxon='%s')" % (char_matrix_oid, row_id, taxon.label))
                        chartype = self._id_chartype_map[chartype_id]
                        state_alphabet = self._id_chartype_map[chartype_id].state_alphabet
                        pos_idx = self._chartype_id_to_pos_map[chartype_id]
                        state = self._id_state_map[ (state_alphabet, nxcell.get('state', None)) ]
                        character_vector.set_at(pos_idx,
                                character_value=state,
                                character_type=chartype)
                        # self._id_state_alphabet_map = {}
                        # self._id_state_map = {}
                        # self._id_chartype_map = {}

            char_matrix[taxon] = character_vector

        # if fixed_state_alphabet:
        #     char_matrix.remap_to_default_state_alphabet_by_symbol(purge_other_state_alphabets=True)

    def parse_ambiguous_state(self, nxstate, state_alphabet):
        """
        Parses an XmlElement represent an ambiguous discrete character state,
        ("uncertain_state_set")
        and returns a corresponding StateAlphabetElement object.
        """
        state_oid = nxstate.get('id', None)
        state_symbol = nxstate.get('symbol', None)
        token = nxstate.get('token', None)
        member_states = []
        for nxmember in nxstate.findall_multistate_member():
            member_state_id = nxmember.get('state', None)
            member_state = self._id_state_map[ (state_alphabet, member_state_id) ]
            member_states.append(member_state)
        state = state_alphabet.new_multistate(symbol=state_symbol,
                state_denomination=state_alphabet.AMBIGUOUS_STATE,
                member_states=member_states)
        assert (state_alphabet, state_oid) not in self._id_state_map
        self._id_state_map[ (state_alphabet, state_oid) ] = state
        if token is not None:
            state_alphabet.new_symbol_synonym(token, state_symbol)
        return state

    def parse_polymorphic_state(self, nxstate, state_alphabet):
        """
        Parses an XmlElement represent a polymorphic discrete character state,
        ("polymorphic_state_set")
        and returns a corresponding StateAlphabetElement object.
        """
        state_oid = nxstate.get('id', None)
        state_symbol = nxstate.get('symbol', None)
        token = nxstate.get('token', None)
        member_states = []
        for nxmember in nxstate.findall_multistate_member():
            member_state_id = nxmember.get('state', None)
            member_state = self._id_state_map[ (state_alphabet, member_state_id) ]
            member_states.append(member_state)
        for nxambiguous in nxstate.findall_uncertain_state_set():
            member_states.append(self.parse_ambiguous_state(nxstate, state_alphabet))
        state = state_alphabet.new_multistate(symbol=state_symbol,
                state_denomination=state_alphabet.POLYMORPHIC_STATE,
                member_states=member_states)
        assert (state_alphabet, state_oid) not in self._id_state_map
        self._id_state_map[ (state_alphabet, state_oid) ] = state
        if token is not None:
            state_alphabet.new_symbol_synonym(token, state_symbol)
        return state

    def parse_state_alphabet(self, nxstates):
        """
        Given an XmlElement representing a nexml definition of (discrete or standard) states
        ("states"), this returns a corresponding StateAlphabet object.
        """
        state_alphabet_oid = nxstates.get("id", None)
        state_alphabet = self._state_alphabet_factory()
        state_alphabet.autocompile_lookup_tables = False
        self._id_state_alphabet_map[state_alphabet_oid] = state_alphabet
        for nxstate in nxstates.findall_char_state():
            state_oid = nxstate.get('id', None)
            state_symbol = nxstate.get('symbol', None)
            state = state_alphabet.new_fundamental_state(symbol=state_symbol)
            token = nxstate.get('token', None)
            assert (state_alphabet, state_oid) not in self._id_state_map
            self._id_state_map[ (state_alphabet, state_oid) ] = state
            if token is not None:
                state_alphabet.new_symbol_synonym(token, state_symbol)
        for nxstate in nxstates.findall_uncertain_state_set():
            self.parse_ambiguous_state(nxstate, state_alphabet)
        for nxstate in nxstates.findall_polymorphic_state_set():
            self.parse_polymorphic_state(nxstate, state_alphabet)
        state_alphabet.autocompile_lookup_tables = True
        state_alphabet.compile_lookup_mappings()
        return state_alphabet

    def parse_characters_format(self, nxformat, data_type, char_matrix):
        """
        Given an XmlElement schema element ("format"), this parses the
        state definitions (if any) and characters (column definitions, if any),
        and populates the given char_matrix accordingly.
        """
        # if data_type == "standard":
        #     for nxstates in nxformat.findall_char_states():
        #         char_matrix.state_alphabets.append(self.parse_state_alphabet(nxstates))
        # else:
        #     pass
        if data_type in ("dna", "rna", "protein", "restriction"):
            # fixed alphabet: map to existing states
            for nxstates in nxformat.findall_char_states():
                state_alphabet_oid = nxstates.get("id", None)
                if state_alphabet_oid is not None:
                    self._id_state_alphabet_map[state_alphabet_oid] = char_matrix.default_state_alphabet
                for nxstate in nxstates.findall_char_state():
                    state_oid = nxstate.get('id', None)
                    label = nxstate.get('label', None)
                    symbol = nxstate.get('symbol', None)
                    token = nxstate.get('token', None)
                    try:
                        state = char_matrix.default_state_alphabet[symbol]
                    except KeyError:
                        raise Exception("'{}' is not a recognized symbol for the state alphabet for the '{}' data type".format(symbol, data_type))
                    assert (char_matrix.default_state_alphabet, state_oid) not in self._id_state_map
                    self._id_state_map[ (char_matrix.default_state_alphabet, state_oid) ] = state
                for nxstate in nxstates.findall_polymorphic_state_set():
                    state_oid = nxstate.get('id', None)
                    symbol = nxstate.get('symbol', None)
                    try:
                        state = char_matrix.default_state_alphabet[symbol]
                    except KeyError:
                        raise Exception("'{}' is not a recognized symbol for the state alphabet for the '{}' data type".format(symbol, data_type))
                    assert (char_matrix.default_state_alphabet, state_oid) not in self._id_state_map
                    self._id_state_map[ (char_matrix.default_state_alphabet, state_oid) ] = state
                for nxstate in nxstates.findall_uncertain_state_set():
                    state_oid = nxstate.get('id', None)
                    symbol = nxstate.get('symbol', None)
                    try:
                        state = char_matrix.default_state_alphabet[symbol]
                    except KeyError:
                        raise Exception("'{}' is not a recognized symbol for the state alphabet for the '{}' data type".format(symbol, data_type))
                    assert (char_matrix.default_state_alphabet, state_oid) not in self._id_state_map
                    self._id_state_map[ (char_matrix.default_state_alphabet, state_oid) ] = state
        else:
            for nxstates in nxformat.findall_char_states():
                char_matrix.state_alphabets.append(self.parse_state_alphabet(nxstates))
        for nxchars in nxformat.findall_char():
            col = char_matrix.new_character_type()
            char_state_set_id = nxchars.get('states')
            if char_state_set_id is not None:
                state_alphabet = self._id_state_alphabet_map.get(char_state_set_id, None)
                if state_alphabet is None:
                    raise Exception("State set '%s' not defined" % char_state_set_id)
                col.state_alphabet = state_alphabet
            elif hasattr(char_matrix, "default_state_alphabet") \
                and char_matrix.default_state_alphabet is not None:
                col.state_alphabet = char_matrix.default_state_alphabet
            char_matrix.character_types.append(col)
            chartype_id = nxchars.get('id')
            self._chartype_id_to_pos_map[chartype_id] = len(self._chartype_id_to_pos_map)
            self._id_chartype_map[chartype_id] = col
            self._char_types.append(col)
            annotations = [i for i in nxchars.findall_annotations()]
            for annotation in annotations:
                self._parse_annotations(col, annotation)

    def create_standard_character_alphabet(self, char_matrix, symbol_list=None):
        """
        Returns a standard character state alphabet based on symbol_list.
        Defaults to '0' - '9' if not specified.
        """
        if symbol_list is None:
            symbol_list = [str(i) for i in range(10)]
        state_alphabet = StateAlphabet()
        for s in symbol_list:
            state_alphabet.append(StateIdentity(symbol=s))
        char_matrix.state_alphabets.append(state_alphabet)
        char_matrix.default_state_alphabet = state_alphabet
