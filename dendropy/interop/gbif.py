#! /usr/bin/env python

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
Wrappers for interacting with GBIF.
"""

import sys
if sys.version_info.major < 3:
    from urllib import urlencode
    from urllib import urlopen
else:
    from urllib.parse import urlencode
    from urllib.request import urlopen
from dendropy.datamodel import basemodel
from dendropy.dataio import xmlprocessing

class GbifXmlElement(xmlprocessing.XmlElement):

    GBIF_NAMESPACE = "http://portal.gbif.org/ws/response/gbif"
    TAXON_OCCURRENCE_NAMESPACE = "http://rs.tdwg.org/ontology/voc/TaxonOccurrence#"
    TAXON_CONCEPT_NAMESPACE = "http://rs.tdwg.org/ontology/voc/TaxonConcept#"
    TAXON_NAME_NAMESPACE = "http://rs.tdwg.org/ontology/voc/TaxonName#"
    RDF_NAMESPACE = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"

    def __init__(self, element, default_namespace=None):
        if default_namespace is None:
            default_namespace = GbifXmlElement.GBIF_NAMESPACE
        xmlprocessing.XmlElement.__init__(self,
                element=element,
                default_namespace=default_namespace)

    def get_about_attr(self):
        return self._element.get("{%s}about" % self.RDF_NAMESPACE)

    # /gbifResponse/dataProviders/dataProvider/dataResources/dataResource/occurrenceRecords/occurrenceRecord
    def iter_taxon_occurrence(self):
        return self.namespaced_getiterator("TaxonOccurrence",
                namespace=self.TAXON_OCCURRENCE_NAMESPACE)

    def _process_ret_val(self, element, text_only=False):
        if text_only:
            if element:
                return element.text
            else:
                return None
        else:
            return element

    def find_institution_code(self, text_only=False):
        e = self.namespaced_find("institutionCode", namespace=self.TAXON_OCCURRENCE_NAMESPACE)
        return self._process_ret_val(e, text_only)

    def find_collection_code(self, text_only=False):
        e = self.namespaced_find("collectionCode", namespace=self.TAXON_OCCURRENCE_NAMESPACE)
        return self._process_ret_val(e, text_only)

    def find_catalog_number(self, text_only=False):
        e = self.namespaced_find("catalogNumber", namespace=self.TAXON_OCCURRENCE_NAMESPACE)
        return self._process_ret_val(e, text_only)

    def find_longitude(self, text_only=False):
        e = self.namespaced_find("decimalLongitude", namespace=self.TAXON_OCCURRENCE_NAMESPACE)
        return self._process_ret_val(e, text_only)

    def find_latitude(self, text_only=False):
        e = self.namespaced_find("decimalLatitude", namespace=self.TAXON_OCCURRENCE_NAMESPACE)
        return self._process_ret_val(e, text_only)

    def find_taxon_name(self, text_only=False):
        # path = "{%(ns)s}identifiedTo/{%(ns)s}Identification/{%(ns)s}taxon_name" % {"ns": self.TAXON_OCCURRENCE_NAMESPACE}
        path = ["identifiedTo", "Identification", "taxonName"]
        e = self.namespaced_find(path, namespace=self.TAXON_OCCURRENCE_NAMESPACE)
        return self._process_ret_val(e, text_only)

class GbifDataProvenance(object):

    def __init__(self, xml=None):
        self.name = None
        self.gbif_key = None
        self.uri = None
        self.rights = None
        self.citation = None
        if xml:
            self.parse_xml(xml)

    def parse_xml(self, xml):
        self.gbif_key = xml.get("gbifKey")
        self.uri = xml.get_about_attr()
        self.name = xml.namespaced_find("name").text
        self.rights = xml.namespaced_find("rights").text
        self.citation = xml.namespaced_find("citation").text

class GbifOccurrenceRecord(object):

    def parse_from_stream(stream):
        xml_doc = xmlprocessing.XmlDocument(file_obj=stream,
                subelement_factory=GbifXmlElement)
        gb_recs = []
        for txo in xml_doc.root.iter_taxon_occurrence():
            gbo = GbifOccurrenceRecord()
            gbo.parse_taxon_occurrence_xml(txo)
            gb_recs.append(gbo)
        return gb_recs
    parse_from_stream = staticmethod(parse_from_stream)

    def __init__(self):
        self.gbif_key = None
        self.uri = None
        self.institution_code = None
        self.collection_code = None
        self.catalog_number = None
        self.taxon_name = None
        self.data_provider = None
        self.data_resource = None
        self._longitude = None
        self._latitude = None

    def subelement_factory(self, element):
        return GbifXmlElement(element)

    def parse_taxon_occurrence_xml(self, txo):
        self.gbif_key = txo.get("gbifKey")
        self.uri = txo.get_about_attr()
        self.institution_code = txo.find_institution_code(text_only=True)
        self.collection_code = txo.find_collection_code(text_only=True)
        self.catalog_number = txo.find_catalog_number(text_only=True)
        self.longitude = txo.find_longitude(text_only=True)
        self.latitude = txo.find_latitude(text_only=True)
        self.taxon_name = txo.find_taxon_name(text_only=True)

    def _get_longitude(self):
        return self._longitude
    def _set_longitude(self, value):
        if value is not None:
            try:
                self._longitude = float(value)
            except ValueError:
                self._longitude = value
    longitude = property(_get_longitude, _set_longitude)

    def _get_latitude(self):
        return self._latitude
    def _set_latitude(self, value):
        if value is not None:
            try:
                self._latitude = float(value)
            except ValueError:
                self._latitude = value
    latitude = property(_get_latitude, _set_latitude)

    def __str__(self):
        return "%s (%s) %s: %s [%s %s]" % (
                self.institution_code,
                self.collection_code,
                self.catalog_number,
                self.taxon_name,
                self.longitude,
                self.latitude)

    def _get_coordinates_as_string(self, sep=" "):
        return "%s%s%s" % (self.longitude, sep, self.latitude)
    coordinates_as_string = property(_get_coordinates_as_string)

    # def as_coordinate_annotation(self,
    #         name=None,
    #         name_prefix=None,
    #         namespace=None,
    #         name_is_prefixed=False,
    #         include_gbif_reference=True,
    #         include_metadata=True,
    #         dynamic=False):
    #     if name is None:
    #         name = "coordinates"
    #     if name_prefix is None or namespace is None:
    #         # name_prefix = "kml"
    #         # namespace = "http://earth.google.com/kml/2.2"
    #         # name_prefix = "ogckml"
    #         # namespace = "http://www.opengis.net/kml/2.2"
    #         name_prefix = "gml"
    #         namespace = "http://www.opengis.net/gml"
    #     if dynamic:
    #         is_attribute = True
    #         value = (self, "coordinates_as_string")
    #     else:
    #         is_attribute = False
    #         value = self.coordinates_as_string
    #     annote = basemodel.Annotation(
    #             name="pos",
    #             value=value,
    #             name_prefix=name_prefix,
    #             namespace=namespace,
    #             name_is_prefixed=name_is_prefixed,
    #             is_attribute=is_attribute,
    #             annotate_as_reference=False,
    #             )
    #     if include_gbif_reference:
    #         if dynamic:
    #             value = (self, "uri")
    #         else:
    #             value = self.uri
    #         subannote = basemodel.Annotation(
    #                 name="source",
    #                 value=value,
    #                 # name_prefix="dc",
    #                 # namespace="http://purl.org/dc/elements/1.1/",
    #                 name_prefix="dcterms",
    #                 namespace="http://purl.org/dc/terms/",
    #                 name_is_prefixed=False,
    #                 is_attribute=is_attribute,
    #                 annotate_as_reference=True,
    #                 )
    #         annote.annotations.add(subannote)
    #     if include_metadata:
    #         for attr in [
    #             ("institution_code", "institutionCode"),
    #             ("collection_code", "collectionCode"),
    #             ("catalog_number", "catalogNumber"),
    #             ("taxon_name", "scientificName"),
    #             ]:
    #             if dynamic:
    #                 value = (self, attr[0])
    #             else:
    #                 value = getattr(self, attr[0])
    #             subannote = basemodel.Annotation(
    #                     name=attr[1],
    #                     value=value,
    #                     name_prefix="dwc",
    #                     namespace="http://rs.tdwg.org/dwc/terms/",
    #                     name_is_prefixed=False,
    #                     is_attribute=is_attribute,
    #                     annotate_as_reference=False,
    #                     )
    #             annote.annotations.add(subannote)
    #     return annote

    def as_annotation(self,
            name="TaxonOccurrence",
            name_prefix="to",
            namespace="http://rs.tdwg.org/ontology/voc/TaxonOccurrence#",
            include_gbif_reference=True,
            dynamic=False):
        """
        Sample output (NeXML)::

            <meta xsi:type="nex:ResourceMeta" rel="to:TaxonOccurrence" id="d4324014736" >
                <meta xsi:type="nex:ResourceMeta" rel="dcterms:source" href="http://data.gbif.org/ws/rest/occurrence/get/44726287" id="d4324014800" />
                <meta xsi:type="nex:LiteralMeta" property="to:decimalLongitude" content="-116.02004" datatype="xsd:float" id="d4324014928" />
                <meta xsi:type="nex:LiteralMeta" property="to:decimalLatitude" content="34.67338" datatype="xsd:float" id="d4324014992" />
                <meta xsi:type="nex:LiteralMeta" property="to:institutionCode" content="ROM" datatype="xsd:string" id="d4324015056" />
                <meta xsi:type="nex:LiteralMeta" property="to:collectionCode" content="Herps" datatype="xsd:string" id="d4324015120" />
                <meta xsi:type="nex:LiteralMeta" property="to:catalogNumber" content="14584" datatype="xsd:string" id="d4324015184" />
                <meta xsi:type="nex:LiteralMeta" property="to:scientificName" content="Crotaphytus bicinctores" datatype="xsd:string" id="d4324015248" />
            </meta>

        """
        # name_prefix="dwc",
        # namespace="http://rs.tdwg.org/dwc/terms/",
        top_node = basemodel.Annotation(
                name=name,
                value=None,
                name_prefix=name_prefix,
                namespace=namespace,
                name_is_prefixed=False,
                is_attribute=False,
                annotate_as_reference=True,
                )
        if dynamic:
            is_attribute=True
        else:
            is_attribute=False
        if include_gbif_reference:
            if dynamic:
                value = (self, "uri")
            else:
                value = self.uri
            subannote = basemodel.Annotation(
                    name="source",
                    value=value,
                    # name_prefix="dc",
                    # namespace="http://purl.org/dc/elements/1.1/",
                    name_prefix="dcterms",
                    namespace="http://purl.org/dc/terms/",
                    name_is_prefixed=False,
                    is_attribute=is_attribute,
                    annotate_as_reference=True,
                    )
            top_node.annotations.add(subannote)
        for attr in [
            ("longitude", "decimalLongitude", "xsd:float"),
            ("latitude", "decimalLatitude", "xsd:float"),
            ("institution_code", "institutionCode", "xsd:string"),
            ("collection_code", "collectionCode", "xsd:string"),
            ("catalog_number", "catalogNumber", "xsd:string"),
            ("taxon_name", "scientificName", "xsd:string"),
            ]:
            if dynamic:
                value = (self, attr[0])
                is_attribute=True
            else:
                value = getattr(self, attr[0])
                is_attribute=False
            subannote = basemodel.Annotation(
                    name=attr[1],
                    value=value,
                    datatype_hint=attr[2],
                    name_prefix=name_prefix,
                    namespace=namespace,
                    name_is_prefixed=False,
                    is_attribute=is_attribute,
                    annotate_as_reference=False,
                    )
            top_node.annotations.add(subannote)
        return top_node

class GbifDb(object):

    def __init__(self):
        self.base_url = None

    def compose_query_url(self, action, query_dict):
        parts = []
        for k, v in query_dict.items():
            parts.append("%s=%s" % (k, v))
        query_url = self.base_url + action + "?" + "&".join(parts)
        return query_url

class GbifOccurrenceDb(GbifDb):

    def __init__(self):
        GbifDb.__init__(self)
        self.base_url = "http://data.gbif.org/ws/rest/occurrence/"

    def fetch_keys(self, **kwargs):
        url = self.compose_query_url(action="list",
                query_dict=kwargs)
        response = urlopen(url)
        return self.parse_list_keys(response)

    def fetch_occurrences(self, **kwargs):
        keys = self.fetch_keys(**kwargs)
        occurrences = []
        for key in keys:
            url = self.compose_query_url(action="get",
                    query_dict={"key": key})
            response = urlopen(url)
            occurrences.extend(self.parse_occurrence_records(response))
        return occurrences
        # url = self.compose_query_url(action="list",
        #         query_dict=kwargs)
        # response = urlopen(url)
        # return self.parse_occurrence_records(response)

    def parse_list_keys(self, stream):
        keys = []
        xml_doc = xmlprocessing.XmlDocument(file_obj=stream,
                subelement_factory=GbifXmlElement)
        xml_root = xml_doc.root
        for txml in xml_root.iter_taxon_occurrence():
            keys.append(txml.get("gbifKey"))
        return keys

    def parse_occurrence_records(self, stream):
        occurrences = GbifOccurrenceRecord.parse_from_stream(stream)
        return occurrences
    #     xml_doc = xmlprocessing.XmlDocument(file_obj=stream,
    #             subelement_factory=GbifXmlElement)
    #     xml_root = xml_doc.root
    #     gbif_recs = []
    #     for dps in xml_root.namespaced_findall("dataProviders", namespace=self.GBIF_NAMESPACE):
    #         for dp in dps.namespaced_findall("dataProvider", namespace=self.GBIF_NAMESPACE):
    #             data_provider = GbifDataProvenance(dp)
    #             for drs in dp.namespaced_findall("dataResources", namespace=self.GBIF_NAMESPACE):
    #                 for dr in drs.namespaced_findall("dataResource", namespace=self.GBIF_NAMESPACE):
    #                     data_resource = GbifDataProvenance(dr)
    #                     for occurs in dr.namedspaced_findall("occurrenceRecords", namespace=self.GBIF_NAMESPACE):
    #                         for occur in occurs.namespaced_findall("TaxonOccurrence", namespace=self.TAXON_OCCURRENCE_NAMESPACE):
    #                             pass

