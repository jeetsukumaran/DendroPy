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
Wrappers for interacting with GBIF.
"""

from dendropy.dataio import xmlparser

class GbifXmlElement(xmlparser.XmlElement):

    GBIF_NAMESPACE = "http://portal.gbif.org/ws/response/gbif"
    TAXON_OCCURRENCE_NAMESPACE = "http://rs.tdwg.org/ontology/voc/TaxonOccurrence#"
    TAXON_CONCEPT_NAMESPACE = "http://rs.tdwg.org/ontology/voc/TaxonConcept#"
    TAXON_NAME_NAMESPACE = "http://rs.tdwg.org/ontology/voc/TaxonName#"
    RDF_NAMESPACE = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"

    def __init__(self, element, default_namespace=None):
        if default_namespace is None:
            default_namespace = GbifXmlElement.GBIF_NAMESPACE
        xmlparser.XmlElement.__init__(self,
                element=element,
                default_namespace=default_namespace)

    def get_about_attr(self):
        return self._element.get("{%s}about" % self.RDF_NAMESPACE)

    # /gbifResponse/dataProviders/dataProvider/dataResources/dataResource/occurrenceRecords/occurrenceRecord
    def iter_taxon_occurrence(self):
        return self.namespaced_getiterator("TaxonOccurrence",
                namespace=self.TAXON_OCCURRENCE_NAMESPACE)

    def find_catalog_number(self):
        return self.namespaced_find("catalogNumber",
                namespace=self.TAXON_OCCURRENCE_NAMESPACE)

    def find_longitude(self):
        return self.namespaced_find("decimalLongitude",
                namespace=self.TAXON_OCCURRENCE_NAMESPACE)

    def find_latitude(self):
        return self.namespaced_find("decimalLatitude",
                namespace=self.TAXON_OCCURRENCE_NAMESPACE)

    def find_taxon_name(self):
        # path = "{%(ns)s}identifiedTo/{%(ns)s}Identification/{%(ns)s}taxon_name" % {"ns": self.TAXON_OCCURRENCE_NAMESPACE}
        path = ["identifiedTo", "Identification", "taxonName"]
        return self.namespaced_find(path, namespace=self.TAXON_OCCURRENCE_NAMESPACE)

class GbifOccurrence(object):

    def parse_from_stream(stream):
        xml_doc = xmlparser.XmlDocument(file_obj=stream,
                subelement_factory=GbifXmlElement)
        gb_recs = []
        for txo in xml_doc.root.iter_taxon_occurrence():
            gbo = GbifOccurrence()
            gbo.parse_taxon_occurrence_xml(txo)
            gb_recs.append(gbo)
        return gb_recs
    parse_from_stream = staticmethod(parse_from_stream)

    def __init__(self):
        self.gbif_key = None
        self.url = None
        self.catalog_number = None
        self.longitude = None
        self.latitude = None
        self.taxon_name = None

    def subelement_factory(self, element):
        return GbifXmlElement(element)

    def parse_taxon_occurrence_xml(self, txo):
        self.gbif_key = txo.get("gbifKey")
        self.url = txo.get_about_attr()
        self.catalog_number = txo.find_catalog_number().text
        self.longitude = txo.find_longitude().text
        self.latitude = txo.find_latitude().text
        self.taxon_name = txo.find_taxon_name().text

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
        for key in keys:
            url = self.compose_query_url(action="get",
                    query_dict={"key": key})
            response = urlopen(url)
            print response.read()

    def parse_list_keys(self, stream):
        keys = []
        xml_doc = xmlparser.XmlDocument(file_obj=stream,
                subelement_factory=GbifXmlElement)
        xml_root = xml_doc.root
        for txml in xml_root.iter_taxon_occurrence():
            keys.append(txml.get("gbifKey"))
        return keys
