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
Annotations testing.
"""

import sys
import os
import unittest
import tempfile
import copy

from dendropy.test.support import pathmap
from dendropy.test.support import datagen
from dendropy.test.support import datatest
import dendropy
from dendropy.dataio import nexml
from dendropy import treesplit
from dendropy import treecalc


class NexmlAnnotations(datatest.AnnotatedDataObjectVerificationTestCase):

    def setUp(self):
        self.tree_src_path = pathmap.tree_source_path("treebase_s373.xml")
        self.prefix_to_namespace = {
                "nex"      :  "http://www.nexml.org/2009",
                ""         :  "http://www.nexml.org/2009",
                "dc"       :  "http://purl.org/dc/elements/1.1/",
                "dcterms"  :  "http://purl.org/dc/terms/",
                "prism"    :  "http://prismstandard.org/namespaces/1.2/basic/",
                "rdf"      :  "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
                "rdfs"     :  "http://www.w3.org/2000/01/rdf-schema#",
                "skos"     :  "http://www.w3.org/2004/02/skos/core#",
                "tb"       :  "http://purl.org/phylo/treebase/2.0/terms#",
                "xsd"      :  "http://www.w3.org/2001/XMLSchema#",
                }
        self.meta = {}
        self.meta["dataset"] = [
                {"content": "Generated on Sat Jun 09 22:14:00 EDT 2012", "datatype": "xsd:string", "id": "meta4928", "property": "skos:changeNote", "type": "nex:LiteralMeta",},
                {"content": "Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlDocumentWriter@5a4b3e1d $Rev: 1060 $", "datatype": "xsd:string", "id": "meta4927", "property": "skos:historyNote", "type": "nex:LiteralMeta",},
                {"content": "109", "datatype": "xsd:string", "id": "meta4926", "property": "prism:volume", "type": "nex:LiteralMeta",},
                {"content": "Zoological Journal of the Linnean Society", "datatype": "xsd:string", "id": "meta4925", "property": "dc:publisher", "type": "nex:LiteralMeta",},
                {"content": "Zoological Journal of the Linnean Society", "datatype": "xsd:string", "id": "meta4924", "property": "prism:publicationName", "type": "nex:LiteralMeta",},
                {"content": "275-299", "datatype": "xsd:string", "id": "meta4923", "property": "prism:pageRange", "type": "nex:LiteralMeta",},
                {"content": "299", "datatype": "xsd:string", "id": "meta4922", "property": "prism:endingPage", "type": "nex:LiteralMeta",},
                {"content": "275", "datatype": "xsd:string", "id": "meta4921", "property": "prism:startingPage", "type": "nex:LiteralMeta",},
                {"content": "1993", "datatype": "xsd:string", "id": "meta4920", "property": "prism:publicationDate", "type": "nex:LiteralMeta",},
                {"content": "Rossman D.", "datatype": "xsd:string", "id": "meta4919", "property": "dc:contributor", "type": "nex:LiteralMeta",},
                {"content": "Wallach V.", "datatype": "xsd:string", "id": "meta4918", "property": "dc:contributor", "type": "nex:LiteralMeta",},
                {"content": "Cundall D.", "datatype": "xsd:string", "id": "meta4917", "property": "dc:contributor", "type": "nex:LiteralMeta",},
                {"content": "Cundall D., Wallach V., & Rossman D.", "datatype": "xsd:string", "id": "meta4916", "property": "dc:creator", "type": "nex:LiteralMeta",},
                {"content": "The systematic relationships of the snake genus Anomochilus.", "datatype": "xsd:string", "id": "meta4915", "property": "dc:title", "type": "nex:LiteralMeta",},
                {"content": "Cundall D., Wallach V., & Rossman D. 1993. The systematic relationships of the snake genus Anomochilus. Zoological Journal of the Linnean Society, 109: 275-299.", "datatype": "xsd:string", "id": "meta4914", "property": "dcterms:bibliographicCitation", "type": "nex:LiteralMeta",},
                {"content": "1998-09-22", "datatype": "xsd:string", "id": "meta4913", "property": "prism:creationDate", "type": "nex:LiteralMeta",},
                {"content": "1998-09-22", "datatype": "xsd:string", "id": "meta4912", "property": "prism:modificationDate", "type": "nex:LiteralMeta",},
                {"content": "1998-09-22", "datatype": "xsd:string", "id": "meta4911", "property": "dc:date", "type": "nex:LiteralMeta",},
                {"content": "S309", "datatype": "xsd:string", "id": "meta4910", "property": "tb:identifier.study.tb1", "type": "nex:LiteralMeta",},
                {"content": "373", "datatype": "xsd:string", "id": "meta4909", "property": "tb:identifier.study", "type": "nex:LiteralMeta",},
                {"content": "Study", "datatype": "xsd:string", "id": "meta4907", "property": "prism:section", "type": "nex:LiteralMeta",},
                ]
        self.meta["taxon_sets"] = {}
        self.meta["taxon_sets"]["Tls9816"] = [
                {"content": "Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlOTUWriter@62f5ae81 $Rev: 1040 $", "datatype": "xsd:string", "id": "meta4930", "property": "skos:historyNote", "type": "nex:LiteralMeta",},
                ]
        self.meta["taxon_sets"]["Tls9817"] = [
                {"content": "Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlOTUWriter@62f5ae81 $Rev: 1040 $", "datatype": "xsd:string", "id": "meta5040", "property": "skos:historyNote", "type": "nex:LiteralMeta",},
                ]
        self.meta["taxon_sets"]["Tls9818"] = [
                {"content": "Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlOTUWriter@62f5ae81 $Rev: 1040 $", "datatype": "xsd:string", "id": "meta5150", "property": "skos:historyNote", "type": "nex:LiteralMeta",}
                ]
        self.meta["taxon"] = {}
        self.meta["taxon"]["Tl52311"] = [
                {"content": "6757", "datatype": "xsd:long", "id": "meta4936", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "16387", "datatype": "xsd:long", "id": "meta4935", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/34989", "id": "meta4934", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:5434416", "id": "meta4933", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4932", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52310"] = [
                {"content": "343", "datatype": "xsd:long", "id": "meta4942", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "846", "datatype": "xsd:long", "id": "meta4941", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/42164", "id": "meta4940", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549759", "id": "meta4939", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4938", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52322"] = [
                {"content": "30007", "datatype": "xsd:long", "id": "meta4948", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "70126", "datatype": "xsd:long", "id": "meta4947", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/39698", "id": "meta4946", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549765", "id": "meta4945", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4944", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52318"] = [
                {"content": "3702", "datatype": "xsd:long", "id": "meta4954", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "8851", "datatype": "xsd:long", "id": "meta4953", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/51855", "id": "meta4952", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549821", "id": "meta4951", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4950", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52308"] = [
                {"content": "10453", "datatype": "xsd:long", "id": "meta4960", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "25017", "datatype": "xsd:long", "id": "meta4959", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/196245", "id": "meta4958", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2757603", "id": "meta4957", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4956", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52315"] = [
                {"content": "3652", "datatype": "xsd:long", "id": "meta4966", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "8705", "datatype": "xsd:long", "id": "meta4965", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/196244", "id": "meta4964", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:1770023", "id": "meta4963", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4962", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52317"] = [
                {"content": "24690", "datatype": "xsd:long", "id": "meta4972", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "57823", "datatype": "xsd:long", "id": "meta4971", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/34984", "id": "meta4970", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2757602", "id": "meta4969", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4968", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52314"] = [
                {"content": "16385", "datatype": "xsd:long", "id": "meta4978", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "38388", "datatype": "xsd:long", "id": "meta4977", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/39076", "id": "meta4976", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549764", "id": "meta4975", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4974", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52319"] = [
                {"content": "31032", "datatype": "xsd:long", "id": "meta4984", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "72453", "datatype": "xsd:long", "id": "meta4983", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/196251", "id": "meta4982", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549767", "id": "meta4981", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4980", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52320"] = [
                {"content": "1768", "datatype": "xsd:long", "id": "meta4990", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "4330", "datatype": "xsd:long", "id": "meta4989", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/327153", "id": "meta4988", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2546805", "id": "meta4987", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4986", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52312"] = [
                {"content": "30325", "datatype": "xsd:long", "id": "meta4995", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "70769", "datatype": "xsd:long", "id": "meta4994", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:5572245", "id": "meta4993", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4992", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52316"] = [
                {"content": "7969", "datatype": "xsd:long", "id": "meta5002", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "19155", "datatype": "xsd:long", "id": "meta5001", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/305692", "id": "meta5000", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"content": "Cylindrophiidae", "datatype": "xsd:string", "id": "meta4999", "property": "skos:altLabel", "type": "nex:LiteralMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549763", "id": "meta4998", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta4997", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52321"] = [
                {"content": "1642", "datatype": "xsd:long", "id": "meta5008", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "4102", "datatype": "xsd:long", "id": "meta5007", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/51842", "id": "meta5006", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549760", "id": "meta5005", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5004", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52313"] = [
                {"content": "15729", "datatype": "xsd:long", "id": "meta5014", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "36870", "datatype": "xsd:long", "id": "meta5013", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/34977", "id": "meta5012", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549783", "id": "meta5011", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5010", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl112723"] = [
                {"content": "30141", "datatype": "xsd:long", "id": "meta5020", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "70385", "datatype": "xsd:long", "id": "meta5019", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/34978", "id": "meta5018", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549784", "id": "meta5017", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5016", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52323"] = [
                {"content": "1760", "datatype": "xsd:long", "id": "meta5026", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "4314", "datatype": "xsd:long", "id": "meta5025", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/42186", "id": "meta5024", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549820", "id": "meta5023", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5022", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl112732"] = [
                {"content": "8926", "datatype": "xsd:long", "id": "meta5032", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "21467", "datatype": "xsd:long", "id": "meta5031", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/261508", "id": "meta5030", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:2549756", "id": "meta5029", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5028", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["taxon"]["Tl52309"] = [
                {"content": "1624", "datatype": "xsd:long", "id": "meta5038", "property": "tb:identifier.taxon", "type": "nex:LiteralMeta",},
                {"content": "4068", "datatype": "xsd:long", "id": "meta5037", "property": "tb:identifier.taxonVariant", "type": "nex:LiteralMeta",},
                {"href": "http://purl.uniprot.org/taxonomy/8548", "id": "meta5036", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://www.ubio.org/authority/metadata.php?lsid=urn:lsid:ubio.org:namebank:5952711", "id": "meta5035", "rel": "skos:closeMatch", "type": "nex:ResourceMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5034", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["tree_lists"] = {}
        self.meta["tree_lists"]["Tb5169"] = [
                {"content": "Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlTreeBlockWriter@5ace1e59 $Rev: 1040 $", "datatype": "xsd:string", "id": "meta5474", "property": "skos:historyNote", "type": "nex:LiteralMeta",},
                {"href": "S373", "id": "meta5473", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["tree_lists"]["Tb5168"] = [
                {"content": "Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlTreeBlockWriter@5ace1e59 $Rev: 1040 $", "datatype": "xsd:string", "id": "meta5474", "property": "skos:historyNote", "type": "nex:LiteralMeta",},
                {"href": "S373", "id": "meta5473", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["tree_lists"]["Tb5167"] = [
                {"content": "Mapped from TreeBASE schema using org.cipres.treebase.domain.nexus.nexml.NexmlTreeBlockWriter@5ace1e59 $Rev: 1040 $", "datatype": "xsd:string", "id": "meta5474", "property": "skos:historyNote", "type": "nex:LiteralMeta",},
                {"href": "S373", "id": "meta5473", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["tree"] = {}
        self.meta["tree"]["Tr3260"] = [
                {"content": "18", "datatype": "xsd:integer", "id": "meta5480", "property": "tb:ntax.tree", "type": "nex:LiteralMeta",},
                {"content": "Unrated", "datatype": "xsd:string", "id": "meta5479", "property": "tb:quality.tree", "type": "nex:LiteralMeta",},
                {"content": "Consensus", "datatype": "xsd:string", "id": "meta5478", "property": "tb:type.tree", "type": "nex:LiteralMeta",},
                {"content": "Species Tree", "datatype": "xsd:string", "id": "meta5477", "property": "tb:kind.tree", "type": "nex:LiteralMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5476", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["tree"]["Tr3258"] = [
                {"content": "18", "datatype": "xsd:integer", "id": "meta5548", "property": "tb:ntax.tree", "type": "nex:LiteralMeta",},
                {"content": "Unrated", "datatype": "xsd:string", "id": "meta5547", "property": "tb:quality.tree", "type": "nex:LiteralMeta",},
                {"content": "Consensus", "datatype": "xsd:string", "id": "meta5546", "property": "tb:type.tree", "type": "nex:LiteralMeta",},
                {"content": "Species Tree", "datatype": "xsd:string", "id": "meta5545", "property": "tb:kind.tree", "type": "nex:LiteralMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5544", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]
        self.meta["tree"]["Tr3259"] = [
                {"content": "18", "datatype": "xsd:integer", "id": "meta5624", "property": "tb:ntax.tree", "type": "nex:LiteralMeta",},
                {"content": "Unrated", "datatype": "xsd:string", "id": "meta5623", "property": "tb:quality.tree", "type": "nex:LiteralMeta",},
                {"content": "Consensus", "datatype": "xsd:string", "id": "meta5622", "property": "tb:type.tree", "type": "nex:LiteralMeta",},
                {"content": "Species Tree", "datatype": "xsd:string", "id": "meta5621", "property": "tb:kind.tree", "type": "nex:LiteralMeta",},
                {"href": "http://purl.org/phylo/treebase/phylows/study/TB2:S373", "id": "meta5620", "rel": "rdfs:isDefinedBy", "type": "nex:ResourceMeta",},
                ]

    def collect_by_qname(self, meta_elements):
        qname_to_content = {}
        for meta in meta_elements:
            if "href" in meta:
                value = meta["href"]
                qname = meta["rel"]
            else:
                value = meta["content"]
                qname = meta["property"]
            try:
                vtype = meta["datatype"]
            except KeyError:
                vtype = None
            try:
                qname_to_content[qname].add((value, vtype,))
            except KeyError:
                qname_to_content[qname] = set([(value, vtype,)])
        return qname_to_content

    def verify_metadata(self, observed, expected):
        self.assertEqual(len(observed), len(expected))
        qname_to_content = self.collect_by_qname(expected)
        for qname in qname_to_content:
            meta_set = qname_to_content[qname]
            annotes = observed.findall(prefixed_name=qname)
            self.assertEqual(len(meta_set), len(annotes))
            obs_set = set()
            for a in annotes:
                obs_set.add((a.value, a.datatype_hint,))
            self.assertEqual(meta_set, obs_set)

    def verify_dataset(self, dataset):
        annotations = dataset.annotations
        self.verify_metadata(
                observed=annotations,
                expected=self.meta["dataset"])
        for taxon_set in dataset.taxon_sets:
            annotations = taxon_set.annotations
            expected = self.meta["taxon_sets"][taxon_set.oid]
            self.verify_metadata(
                    observed=annotations,
                    expected=expected)
            for taxon in taxon_set:
                annotations = taxon.annotations
                expected = self.meta["taxon"][taxon.oid]
                self.verify_metadata(
                        observed=annotations,
                        expected=expected)
        for tree_list in dataset.tree_lists:
            annotations = tree_list.annotations
            expected = self.meta["tree_lists"][tree_list.oid]
            self.verify_metadata(
                    observed=annotations,
                    expected=expected)
            for tree in tree_list:
                annotations = tree.annotations
                expected = self.meta["tree"][tree.oid]
                self.verify_metadata(
                        observed=annotations,
                        expected=expected)

    def testParseNexmlAnnotations(self):
        s = pathmap.tree_source_stream("treebase_s373.xml")
        dataset = dendropy.DataSet(stream=s, schema="nexml")
        self.verify_dataset(dataset)

    # def testClonedNexmlAnnotations(self):
    #     s = pathmap.tree_source_stream("treebase_s373.xml")
    #     dataset1 = dendropy.DataSet(stream=s, schema="nexml")
    #     dataset2 = dendropy.DataSet(dataset1)
    #     self.verify_dataset(dataset2)

if __name__ == "__main__":
    unittest.main()


