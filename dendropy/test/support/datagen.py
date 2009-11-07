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
Data and value generation
"""

from random import Random
from dendropy.test.support import pathmap
from dendropy.test.support.framework import NodeRelationship
import dendropy

def four_taxon_tree1():
    taxa = dendropy.TaxonSet(['A', 'B', 'C', 'D'])
    tree = dendropy.Tree(taxon_set=taxa)
    assert tree.taxon_set == taxa
    tree.seed_node.oid = 'root'
    tree.seed_node.label = 'root'
    tree.seed_node.edge.length = 2.0
    i1 = tree.seed_node.new_child(oid='i1', label='i1')
    i1.edge.length = 2.5
    a = i1.new_child(oid='a', taxon=taxa.require_taxon(label='A'))
    a.edge.length = 3.5
    b = i1.new_child(oid='b', taxon=taxa.require_taxon(label='B'))
    b.edge.length = 3.5
    i2 = tree.seed_node.new_child(oid='i2', label='i2')
    i2.edge.length = 4.0
    c = i2.new_child(oid='c', taxon=taxa.require_taxon(label='C'))
    c.edge.length = 2.0
    d = i2.new_child(oid='d', taxon=taxa.require_taxon(label='D'))
    d.edge.length = 2.0
    return tree

def reference_tree_list_postorder_node_labels():
    return [
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Morelia boeleni','Node4837936','Node4837616','Morelia viridis','Node4837552','Antaresia maculosa','Node4837488','Antaresia stimsoni','Node4837424','Python timoriensis','Node4796016','Morelia oenpelliensis','Node4796304','Morelia bredli','Node4796336','Antaresia perthensis','Python brongersmai','Node4838576','Node4796272','Morelia carinata','Node4796144'],
        ['Aspidites ramsayi','Bothrochilus boa','Morelia boeleni','Morelia oenpelliensis','Node4839664','Node4795856','Morelia viridis','Node4839152','Antaresia perthensis','Node4839344','Liasis fuscus','Antaresia stimsoni','Python timoriensis','Morelia carinata','Node4840208','Python brongersmai','Node4840176','Antaresia maculosa','Node4840304','Node4840144','Node4839824','Morelia bredli','Node4839696','Node4839216'],
        ['Aspidites ramsayi','Bothrochilus boa','Antaresia perthensis','Node4841008','Python brongersmai','Node4840560','Python timoriensis','Antaresia maculosa','Node4841392','Morelia oenpelliensis','Node4841200','Node4840912','Liasis fuscus','Morelia viridis','Morelia bredli','Node4895152','Node4895088','Node4840720','Antaresia stimsoni','Node4840880','Morelia carinata','Morelia boeleni','Node4895344','Node4840752'],
        ['Aspidites ramsayi','Bothrochilus boa','Antaresia perthensis','Node4795920','Morelia bredli','Node4895632','Liasis fuscus','Morelia oenpelliensis','Node4896336','Morelia viridis','Node4896240','Python timoriensis','Morelia carinata','Node4896784','Antaresia maculosa','Node4896528','Node4896304','Python brongersmai','Node4896208','Antaresia stimsoni','Node4896016','Node4895824','Morelia boeleni','Node4895696'],
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Morelia carinata','Node4897680','Node4897488','Antaresia stimsoni','Antaresia perthensis','Node4896720','Antaresia maculosa','Node4840656','Morelia viridis','Morelia bredli','Morelia oenpelliensis','Node4898288','Node4898064','Python timoriensis','Node4897872','Node4897712','Python brongersmai','Morelia boeleni','Node4898608','Node4897520','Node4897360'],
        ['Aspidites ramsayi','Bothrochilus boa','Python brongersmai','Node4897168','Python timoriensis','Node5066960','Morelia bredli','Antaresia maculosa','Node5067120','Node4898480','Liasis fuscus','Morelia oenpelliensis','Node5067344','Antaresia stimsoni','Morelia carinata','Morelia boeleni','Node5068080','Node5067888','Antaresia perthensis','Node5067856','Node5067056','Morelia viridis','Node5067600','Node5066800'],
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Node5068592','Antaresia stimsoni','Morelia viridis','Python timoriensis','Morelia boeleni','Node5069392','Node5069296','Node5068976','Antaresia maculosa','Morelia oenpelliensis','Node5069424','Node5068912','Morelia bredli','Node5068656','Morelia carinata','Node5068880','Python brongersmai','Node5068816','Antaresia perthensis','Node4840688','Node5068464'],
        ['Aspidites ramsayi','Bothrochilus boa','Python timoriensis','Node5070448','Antaresia perthensis','Node5070384','Antaresia maculosa','Morelia oenpelliensis','Node5070640','Node5070320','Morelia carinata','Node5070256','Antaresia stimsoni','Node5069808','Morelia boeleni','Node5070160','Morelia viridis','Python brongersmai','Node5124176','Node5069968','Liasis fuscus','Node5070128','Morelia bredli','Node5070000'],
        ['Aspidites ramsayi','Bothrochilus boa','Python brongersmai','Node5125072','Liasis fuscus','Antaresia maculosa','Node5125456','Antaresia stimsoni','Node5125392','Morelia bredli','Antaresia perthensis','Python timoriensis','Node5126032','Node5125648','Node5125136','Morelia carinata','Node5125360','Morelia oenpelliensis','Node5125296','Morelia viridis','Morelia boeleni','Node5125584','Node4796112','Node5124944'],
        ['Aspidites ramsayi','Bothrochilus boa','Morelia bredli','Python timoriensis','Node5127120','Node5126384','Liasis fuscus','Morelia viridis','Node5126320','Morelia carinata','Node5069904','Node5126768','Antaresia perthensis','Antaresia maculosa','Morelia boeleni','Node5127280','Morelia oenpelliensis','Node5127568','Node5127472','Node5126256','Antaresia stimsoni','Node5126736','Python brongersmai','Node5126608'],
    ]

def reference_tree_list_newick_string():
    return """\
        ('Aspidites ramsayi':0.056823,(((((((('Bothrochilus boa':0.073578,('Liasis fuscus':0.074347,'Morelia boeleni':0.07628)Node4837936:0.003063)Node4837616:0.011155,'Morelia viridis':0.066566)Node4837552:0.010103,'Antaresia maculosa':0.080749)Node4837488:0.007681,'Antaresia stimsoni':0.061996)Node4837424:0.0,'Python timoriensis':0.113605)Node4796016:0.002264,'Morelia oenpelliensis':0.074047)Node4796304:0.007655,'Morelia bredli':0.065046)Node4796336:0.003495,('Antaresia perthensis':0.065004,'Python brongersmai':0.107706)Node4838576:0.018613)Node4796272:0.011159,'Morelia carinata':0.079321)Node4796144;
        ('Aspidites ramsayi':0.065297,((('Bothrochilus boa':0.073565,('Morelia boeleni':0.069885,'Morelia oenpelliensis':0.06247)Node4839664:0.007879)Node4795856:0.016022,'Morelia viridis':0.066015)Node4839152:0.010946,'Antaresia perthensis':0.079424)Node4839344:0.001444,(('Liasis fuscus':0.082049,('Antaresia stimsoni':0.056894,((('Python timoriensis':0.092077,'Morelia carinata':0.0579)Node4840208:0.015909,'Python brongersmai':0.114385)Node4840176:0.010752,'Antaresia maculosa':0.079792)Node4840304:0.004552)Node4840144:0.006907)Node4839824:0.008007,'Morelia bredli':0.063725)Node4839696:0.004077)Node4839216;
        ('Aspidites ramsayi':0.06111,((((('Bothrochilus boa':0.067908,'Antaresia perthensis':0.077196)Node4841008:0.008688,'Python brongersmai':0.1069)Node4840560:0.010367,(('Python timoriensis':0.109318,'Antaresia maculosa':0.081382)Node4841392:0.004129,'Morelia oenpelliensis':0.068424)Node4841200:0.007731)Node4840912:0.005167,('Liasis fuscus':0.079478,('Morelia viridis':0.062239,'Morelia bredli':0.060484)Node4895152:0.013155)Node4895088:0.00574)Node4840720:0.006988,'Antaresia stimsoni':0.063808)Node4840880:0.005408,('Morelia carinata':0.069501,'Morelia boeleni':0.073602)Node4895344:0.012121)Node4840752;
        ('Aspidites ramsayi':0.053028,((('Bothrochilus boa':0.072831,'Antaresia perthensis':0.073004)Node4795920:0.012413,'Morelia bredli':0.069933)Node4895632:0.002567,((((('Liasis fuscus':0.072339,'Morelia oenpelliensis':0.062058)Node4896336:0.011723,'Morelia viridis':0.065297)Node4896240:0.011519,(('Python timoriensis':0.093762,'Morelia carinata':0.055992)Node4896784:0.022912,'Antaresia maculosa':0.079945)Node4896528:0.007647)Node4896304:5.134e-05,'Python brongersmai':0.118976)Node4896208:0.004803,'Antaresia stimsoni':0.064721)Node4896016:0.005969)Node4895824:0.014376,'Morelia boeleni':0.076894)Node4895696;
        ('Aspidites ramsayi':0.06391,('Bothrochilus boa':0.071835,('Liasis fuscus':0.075788,'Morelia carinata':0.07116)Node4897680:0.010888)Node4897488:0.00451,(((('Antaresia stimsoni':0.048754,'Antaresia perthensis':0.069682)Node4896720:0.012356,'Antaresia maculosa':0.077058)Node4840656:0.013817,(('Morelia viridis':0.07026,('Morelia bredli':0.066182,'Morelia oenpelliensis':0.068166)Node4898288:0.001443)Node4898064:0.003147,'Python timoriensis':0.108243)Node4897872:0.012897)Node4897712:0.011944,('Python brongersmai':0.108072,'Morelia boeleni':0.073251)Node4898608:0.00612)Node4897520:0.005404)Node4897360;
        ('Aspidites ramsayi':0.06501,((('Bothrochilus boa':0.059984,'Python brongersmai':0.100555)Node4897168:0.017824,'Python timoriensis':0.102639)Node5066960:0.011722,('Morelia bredli':0.0538,'Antaresia maculosa':0.07843)Node5067120:0.014836)Node4898480:0.003644,((('Liasis fuscus':0.068438,'Morelia oenpelliensis':0.066815)Node5067344:0.015326,(('Antaresia stimsoni':0.057007,('Morelia carinata':0.062766,'Morelia boeleni':0.079145)Node5068080:0.013392)Node5067888:0.010508,'Antaresia perthensis':0.083409)Node5067856:0.005897)Node5067056:0.0,'Morelia viridis':0.072091)Node5067600:0.008523)Node5066800;
        ('Aspidites ramsayi':0.062297,('Bothrochilus boa':0.069165,'Liasis fuscus':0.073495)Node5068592:0.009899,(((((('Antaresia stimsoni':0.061022,('Morelia viridis':0.071291,('Python timoriensis':0.098776,'Morelia boeleni':0.075069)Node5069392:0.017245)Node5069296:0.004698)Node5068976:0.003603,('Antaresia maculosa':0.075494,'Morelia oenpelliensis':0.063757)Node5069424:0.017571)Node5068912:0.000386,'Morelia bredli':0.070614)Node5068656:0.00713,'Morelia carinata':0.079534)Node5068880:0.010249,'Python brongersmai':0.114381)Node5068816:0.003672,'Antaresia perthensis':0.080896)Node4840688:0.003212)Node5068464;
        ('Aspidites ramsayi':0.05678,(((((((('Bothrochilus boa':0.065928,'Python timoriensis':0.10077)Node5070448:0.013813,'Antaresia perthensis':0.083335)Node5070384:0.004377,('Antaresia maculosa':0.073392,'Morelia oenpelliensis':0.065748)Node5070640:0.016433)Node5070320:0.005104,'Morelia carinata':0.07646)Node5070256:0.001866,'Antaresia stimsoni':0.06129)Node5069808:0.012314,'Morelia boeleni':0.081982)Node5070160:0.002833,('Morelia viridis':0.064242,'Python brongersmai':0.109434)Node5124176:0.010734)Node5069968:0.00182,'Liasis fuscus':0.079156)Node5070128:0.014479,'Morelia bredli':0.060374)Node5070000;
        ('Aspidites ramsayi':0.063004,('Bothrochilus boa':0.061413,'Python brongersmai':0.099507)Node5125072:0.018059,(((((('Liasis fuscus':0.075877,'Antaresia maculosa':0.071375)Node5125456:0.014459,'Antaresia stimsoni':0.054698)Node5125392:0.002182,('Morelia bredli':0.062291,('Antaresia perthensis':0.074463,'Python timoriensis':0.103887)Node5126032:0.00521)Node5125648:0.011197)Node5125136:0.010763,'Morelia carinata':0.070762)Node5125360:0.016533,'Morelia oenpelliensis':0.068818)Node5125296:0.003714,('Morelia viridis':0.071307,'Morelia boeleni':0.073392)Node5125584:0.010478)Node4796112:0.00238)Node5124944;
        ('Aspidites ramsayi':0.053901,(((('Bothrochilus boa':0.07927,('Morelia bredli':0.053956,'Python timoriensis':0.098874)Node5127120:0.016769)Node5126384:0.003041,(('Liasis fuscus':0.075223,'Morelia viridis':0.065805)Node5126320:0.010712,'Morelia carinata':0.075934)Node5069904:0.008083)Node5126768:0.00777,('Antaresia perthensis':0.082425,(('Antaresia maculosa':0.075429,'Morelia boeleni':0.0787)Node5127280:0.010583,'Morelia oenpelliensis':0.070289)Node5127568:0.008508)Node5127472:0.00373)Node5126256:0.002034,'Antaresia stimsoni':0.064652)Node5126736:0.013947,'Python brongersmai':0.107868)Node5126608;
    """

def reference_tree_list_node_relationships():
    treelist_node_references = [
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node4796144', child_labels=[], edge_length=0.056823, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4837616', child_labels=[], edge_length=0.073578, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4837936', child_labels=[], edge_length=0.074347, taxon_label='Liasis fuscus'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4837936', child_labels=[], edge_length=0.07628, taxon_label='Morelia boeleni'),
            'Node4837936' : NodeRelationship(parent_label='Node4837616', child_labels=['Liasis fuscus','Morelia boeleni'], edge_length=0.003063, taxon_label=None),
            'Node4837616' : NodeRelationship(parent_label='Node4837552', child_labels=['Bothrochilus boa','Node4837936'], edge_length=0.011155, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node4837552', child_labels=[], edge_length=0.066566, taxon_label='Morelia viridis'),
            'Node4837552' : NodeRelationship(parent_label='Node4837488', child_labels=['Node4837616','Morelia viridis'], edge_length=0.010103, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4837488', child_labels=[], edge_length=0.080749, taxon_label='Antaresia maculosa'),
            'Node4837488' : NodeRelationship(parent_label='Node4837424', child_labels=['Node4837552','Antaresia maculosa'], edge_length=0.007681, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4837424', child_labels=[], edge_length=0.061996, taxon_label='Antaresia stimsoni'),
            'Node4837424' : NodeRelationship(parent_label='Node4796016', child_labels=['Node4837488','Antaresia stimsoni'], edge_length=0.0, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4796016', child_labels=[], edge_length=0.113605, taxon_label='Python timoriensis'),
            'Node4796016' : NodeRelationship(parent_label='Node4796304', child_labels=['Node4837424','Python timoriensis'], edge_length=0.002264, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4796304', child_labels=[], edge_length=0.074047, taxon_label='Morelia oenpelliensis'),
            'Node4796304' : NodeRelationship(parent_label='Node4796336', child_labels=['Node4796016','Morelia oenpelliensis'], edge_length=0.007655, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node4796336', child_labels=[], edge_length=0.065046, taxon_label='Morelia bredli'),
            'Node4796336' : NodeRelationship(parent_label='Node4796272', child_labels=['Node4796304','Morelia bredli'], edge_length=0.003495, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4838576', child_labels=[], edge_length=0.065004, taxon_label='Antaresia perthensis'),
            'Python brongersmai' : NodeRelationship(parent_label='Node4838576', child_labels=[], edge_length=0.107706, taxon_label='Python brongersmai'),
            'Node4838576' : NodeRelationship(parent_label='Node4796272', child_labels=['Antaresia perthensis','Python brongersmai'], edge_length=0.018613, taxon_label=None),
            'Node4796272' : NodeRelationship(parent_label='Node4796144', child_labels=['Node4796336','Node4838576'], edge_length=0.011159, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4796144', child_labels=[], edge_length=0.079321, taxon_label='Morelia carinata'),
            'Node4796144' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node4796272','Morelia carinata'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node4839216', child_labels=[], edge_length=0.065297, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4795856', child_labels=[], edge_length=0.073565, taxon_label='Bothrochilus boa'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4839664', child_labels=[], edge_length=0.069885, taxon_label='Morelia boeleni'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4839664', child_labels=[], edge_length=0.06247, taxon_label='Morelia oenpelliensis'),
            'Node4839664' : NodeRelationship(parent_label='Node4795856', child_labels=['Morelia boeleni','Morelia oenpelliensis'], edge_length=0.007879, taxon_label=None),
            'Node4795856' : NodeRelationship(parent_label='Node4839152', child_labels=['Bothrochilus boa','Node4839664'], edge_length=0.016022, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node4839152', child_labels=[], edge_length=0.066015, taxon_label='Morelia viridis'),
            'Node4839152' : NodeRelationship(parent_label='Node4839344', child_labels=['Node4795856','Morelia viridis'], edge_length=0.010946, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4839344', child_labels=[], edge_length=0.079424, taxon_label='Antaresia perthensis'),
            'Node4839344' : NodeRelationship(parent_label='Node4839216', child_labels=['Node4839152','Antaresia perthensis'], edge_length=0.001444, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4839824', child_labels=[], edge_length=0.082049, taxon_label='Liasis fuscus'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4840144', child_labels=[], edge_length=0.056894, taxon_label='Antaresia stimsoni'),
            'Python timoriensis' : NodeRelationship(parent_label='Node4840208', child_labels=[], edge_length=0.092077, taxon_label='Python timoriensis'),
            'Morelia carinata' : NodeRelationship(parent_label='Node4840208', child_labels=[], edge_length=0.0579, taxon_label='Morelia carinata'),
            'Node4840208' : NodeRelationship(parent_label='Node4840176', child_labels=['Python timoriensis','Morelia carinata'], edge_length=0.015909, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node4840176', child_labels=[], edge_length=0.114385, taxon_label='Python brongersmai'),
            'Node4840176' : NodeRelationship(parent_label='Node4840304', child_labels=['Node4840208','Python brongersmai'], edge_length=0.010752, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4840304', child_labels=[], edge_length=0.079792, taxon_label='Antaresia maculosa'),
            'Node4840304' : NodeRelationship(parent_label='Node4840144', child_labels=['Node4840176','Antaresia maculosa'], edge_length=0.004552, taxon_label=None),
            'Node4840144' : NodeRelationship(parent_label='Node4839824', child_labels=['Antaresia stimsoni','Node4840304'], edge_length=0.006907, taxon_label=None),
            'Node4839824' : NodeRelationship(parent_label='Node4839696', child_labels=['Liasis fuscus','Node4840144'], edge_length=0.008007, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node4839696', child_labels=[], edge_length=0.063725, taxon_label='Morelia bredli'),
            'Node4839696' : NodeRelationship(parent_label='Node4839216', child_labels=['Node4839824','Morelia bredli'], edge_length=0.004077, taxon_label=None),
            'Node4839216' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node4839344','Node4839696'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node4840752', child_labels=[], edge_length=0.06111, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4841008', child_labels=[], edge_length=0.067908, taxon_label='Bothrochilus boa'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4841008', child_labels=[], edge_length=0.077196, taxon_label='Antaresia perthensis'),
            'Node4841008' : NodeRelationship(parent_label='Node4840560', child_labels=['Bothrochilus boa','Antaresia perthensis'], edge_length=0.008688, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node4840560', child_labels=[], edge_length=0.1069, taxon_label='Python brongersmai'),
            'Node4840560' : NodeRelationship(parent_label='Node4840912', child_labels=['Node4841008','Python brongersmai'], edge_length=0.010367, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4841392', child_labels=[], edge_length=0.109318, taxon_label='Python timoriensis'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4841392', child_labels=[], edge_length=0.081382, taxon_label='Antaresia maculosa'),
            'Node4841392' : NodeRelationship(parent_label='Node4841200', child_labels=['Python timoriensis','Antaresia maculosa'], edge_length=0.004129, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4841200', child_labels=[], edge_length=0.068424, taxon_label='Morelia oenpelliensis'),
            'Node4841200' : NodeRelationship(parent_label='Node4840912', child_labels=['Node4841392','Morelia oenpelliensis'], edge_length=0.007731, taxon_label=None),
            'Node4840912' : NodeRelationship(parent_label='Node4840720', child_labels=['Node4840560','Node4841200'], edge_length=0.005167, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4895088', child_labels=[], edge_length=0.079478, taxon_label='Liasis fuscus'),
            'Morelia viridis' : NodeRelationship(parent_label='Node4895152', child_labels=[], edge_length=0.062239, taxon_label='Morelia viridis'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4895152', child_labels=[], edge_length=0.060484, taxon_label='Morelia bredli'),
            'Node4895152' : NodeRelationship(parent_label='Node4895088', child_labels=['Morelia viridis','Morelia bredli'], edge_length=0.013155, taxon_label=None),
            'Node4895088' : NodeRelationship(parent_label='Node4840720', child_labels=['Liasis fuscus','Node4895152'], edge_length=0.00574, taxon_label=None),
            'Node4840720' : NodeRelationship(parent_label='Node4840880', child_labels=['Node4840912','Node4895088'], edge_length=0.006988, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4840880', child_labels=[], edge_length=0.063808, taxon_label='Antaresia stimsoni'),
            'Node4840880' : NodeRelationship(parent_label='Node4840752', child_labels=['Node4840720','Antaresia stimsoni'], edge_length=0.005408, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4895344', child_labels=[], edge_length=0.069501, taxon_label='Morelia carinata'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4895344', child_labels=[], edge_length=0.073602, taxon_label='Morelia boeleni'),
            'Node4895344' : NodeRelationship(parent_label='Node4840752', child_labels=['Morelia carinata','Morelia boeleni'], edge_length=0.012121, taxon_label=None),
            'Node4840752' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node4840880','Node4895344'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node4895696', child_labels=[], edge_length=0.053028, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4795920', child_labels=[], edge_length=0.072831, taxon_label='Bothrochilus boa'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4795920', child_labels=[], edge_length=0.073004, taxon_label='Antaresia perthensis'),
            'Node4795920' : NodeRelationship(parent_label='Node4895632', child_labels=['Bothrochilus boa','Antaresia perthensis'], edge_length=0.012413, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node4895632', child_labels=[], edge_length=0.069933, taxon_label='Morelia bredli'),
            'Node4895632' : NodeRelationship(parent_label='Node4895824', child_labels=['Node4795920','Morelia bredli'], edge_length=0.002567, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4896336', child_labels=[], edge_length=0.072339, taxon_label='Liasis fuscus'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4896336', child_labels=[], edge_length=0.062058, taxon_label='Morelia oenpelliensis'),
            'Node4896336' : NodeRelationship(parent_label='Node4896240', child_labels=['Liasis fuscus','Morelia oenpelliensis'], edge_length=0.011723, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node4896240', child_labels=[], edge_length=0.065297, taxon_label='Morelia viridis'),
            'Node4896240' : NodeRelationship(parent_label='Node4896304', child_labels=['Node4896336','Morelia viridis'], edge_length=0.011519, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4896784', child_labels=[], edge_length=0.093762, taxon_label='Python timoriensis'),
            'Morelia carinata' : NodeRelationship(parent_label='Node4896784', child_labels=[], edge_length=0.055992, taxon_label='Morelia carinata'),
            'Node4896784' : NodeRelationship(parent_label='Node4896528', child_labels=['Python timoriensis','Morelia carinata'], edge_length=0.022912, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4896528', child_labels=[], edge_length=0.079945, taxon_label='Antaresia maculosa'),
            'Node4896528' : NodeRelationship(parent_label='Node4896304', child_labels=['Node4896784','Antaresia maculosa'], edge_length=0.007647, taxon_label=None),
            'Node4896304' : NodeRelationship(parent_label='Node4896208', child_labels=['Node4896240','Node4896528'], edge_length=5.134e-05, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node4896208', child_labels=[], edge_length=0.118976, taxon_label='Python brongersmai'),
            'Node4896208' : NodeRelationship(parent_label='Node4896016', child_labels=['Node4896304','Python brongersmai'], edge_length=0.004803, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4896016', child_labels=[], edge_length=0.064721, taxon_label='Antaresia stimsoni'),
            'Node4896016' : NodeRelationship(parent_label='Node4895824', child_labels=['Node4896208','Antaresia stimsoni'], edge_length=0.005969, taxon_label=None),
            'Node4895824' : NodeRelationship(parent_label='Node4895696', child_labels=['Node4895632','Node4896016'], edge_length=0.014376, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4895696', child_labels=[], edge_length=0.076894, taxon_label='Morelia boeleni'),
            'Node4895696' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node4895824','Morelia boeleni'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node4897360', child_labels=[], edge_length=0.06391, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4897488', child_labels=[], edge_length=0.071835, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4897680', child_labels=[], edge_length=0.075788, taxon_label='Liasis fuscus'),
            'Morelia carinata' : NodeRelationship(parent_label='Node4897680', child_labels=[], edge_length=0.07116, taxon_label='Morelia carinata'),
            'Node4897680' : NodeRelationship(parent_label='Node4897488', child_labels=['Liasis fuscus','Morelia carinata'], edge_length=0.010888, taxon_label=None),
            'Node4897488' : NodeRelationship(parent_label='Node4897360', child_labels=['Bothrochilus boa','Node4897680'], edge_length=0.00451, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4896720', child_labels=[], edge_length=0.048754, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4896720', child_labels=[], edge_length=0.069682, taxon_label='Antaresia perthensis'),
            'Node4896720' : NodeRelationship(parent_label='Node4840656', child_labels=['Antaresia stimsoni','Antaresia perthensis'], edge_length=0.012356, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4840656', child_labels=[], edge_length=0.077058, taxon_label='Antaresia maculosa'),
            'Node4840656' : NodeRelationship(parent_label='Node4897712', child_labels=['Node4896720','Antaresia maculosa'], edge_length=0.013817, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node4898064', child_labels=[], edge_length=0.07026, taxon_label='Morelia viridis'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4898288', child_labels=[], edge_length=0.066182, taxon_label='Morelia bredli'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4898288', child_labels=[], edge_length=0.068166, taxon_label='Morelia oenpelliensis'),
            'Node4898288' : NodeRelationship(parent_label='Node4898064', child_labels=['Morelia bredli','Morelia oenpelliensis'], edge_length=0.001443, taxon_label=None),
            'Node4898064' : NodeRelationship(parent_label='Node4897872', child_labels=['Morelia viridis','Node4898288'], edge_length=0.003147, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4897872', child_labels=[], edge_length=0.108243, taxon_label='Python timoriensis'),
            'Node4897872' : NodeRelationship(parent_label='Node4897712', child_labels=['Node4898064','Python timoriensis'], edge_length=0.012897, taxon_label=None),
            'Node4897712' : NodeRelationship(parent_label='Node4897520', child_labels=['Node4840656','Node4897872'], edge_length=0.011944, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node4898608', child_labels=[], edge_length=0.108072, taxon_label='Python brongersmai'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4898608', child_labels=[], edge_length=0.073251, taxon_label='Morelia boeleni'),
            'Node4898608' : NodeRelationship(parent_label='Node4897520', child_labels=['Python brongersmai','Morelia boeleni'], edge_length=0.00612, taxon_label=None),
            'Node4897520' : NodeRelationship(parent_label='Node4897360', child_labels=['Node4897712','Node4898608'], edge_length=0.005404, taxon_label=None),
            'Node4897360' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node4897488','Node4897520'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node5066800', child_labels=[], edge_length=0.06501, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4897168', child_labels=[], edge_length=0.059984, taxon_label='Bothrochilus boa'),
            'Python brongersmai' : NodeRelationship(parent_label='Node4897168', child_labels=[], edge_length=0.100555, taxon_label='Python brongersmai'),
            'Node4897168' : NodeRelationship(parent_label='Node5066960', child_labels=['Bothrochilus boa','Python brongersmai'], edge_length=0.017824, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node5066960', child_labels=[], edge_length=0.102639, taxon_label='Python timoriensis'),
            'Node5066960' : NodeRelationship(parent_label='Node4898480', child_labels=['Node4897168','Python timoriensis'], edge_length=0.011722, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node5067120', child_labels=[], edge_length=0.0538, taxon_label='Morelia bredli'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5067120', child_labels=[], edge_length=0.07843, taxon_label='Antaresia maculosa'),
            'Node5067120' : NodeRelationship(parent_label='Node4898480', child_labels=['Morelia bredli','Antaresia maculosa'], edge_length=0.014836, taxon_label=None),
            'Node4898480' : NodeRelationship(parent_label='Node5066800', child_labels=['Node5066960','Node5067120'], edge_length=0.003644, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node5067344', child_labels=[], edge_length=0.068438, taxon_label='Liasis fuscus'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node5067344', child_labels=[], edge_length=0.066815, taxon_label='Morelia oenpelliensis'),
            'Node5067344' : NodeRelationship(parent_label='Node5067056', child_labels=['Liasis fuscus','Morelia oenpelliensis'], edge_length=0.015326, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node5067888', child_labels=[], edge_length=0.057007, taxon_label='Antaresia stimsoni'),
            'Morelia carinata' : NodeRelationship(parent_label='Node5068080', child_labels=[], edge_length=0.062766, taxon_label='Morelia carinata'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node5068080', child_labels=[], edge_length=0.079145, taxon_label='Morelia boeleni'),
            'Node5068080' : NodeRelationship(parent_label='Node5067888', child_labels=['Morelia carinata','Morelia boeleni'], edge_length=0.013392, taxon_label=None),
            'Node5067888' : NodeRelationship(parent_label='Node5067856', child_labels=['Antaresia stimsoni','Node5068080'], edge_length=0.010508, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node5067856', child_labels=[], edge_length=0.083409, taxon_label='Antaresia perthensis'),
            'Node5067856' : NodeRelationship(parent_label='Node5067056', child_labels=['Node5067888','Antaresia perthensis'], edge_length=0.005897, taxon_label=None),
            'Node5067056' : NodeRelationship(parent_label='Node5067600', child_labels=['Node5067344','Node5067856'], edge_length=0.0, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node5067600', child_labels=[], edge_length=0.072091, taxon_label='Morelia viridis'),
            'Node5067600' : NodeRelationship(parent_label='Node5066800', child_labels=['Node5067056','Morelia viridis'], edge_length=0.008523, taxon_label=None),
            'Node5066800' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node4898480','Node5067600'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node5068464', child_labels=[], edge_length=0.062297, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node5068592', child_labels=[], edge_length=0.069165, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node5068592', child_labels=[], edge_length=0.073495, taxon_label='Liasis fuscus'),
            'Node5068592' : NodeRelationship(parent_label='Node5068464', child_labels=['Bothrochilus boa','Liasis fuscus'], edge_length=0.009899, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node5068976', child_labels=[], edge_length=0.061022, taxon_label='Antaresia stimsoni'),
            'Morelia viridis' : NodeRelationship(parent_label='Node5069296', child_labels=[], edge_length=0.071291, taxon_label='Morelia viridis'),
            'Python timoriensis' : NodeRelationship(parent_label='Node5069392', child_labels=[], edge_length=0.098776, taxon_label='Python timoriensis'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node5069392', child_labels=[], edge_length=0.075069, taxon_label='Morelia boeleni'),
            'Node5069392' : NodeRelationship(parent_label='Node5069296', child_labels=['Python timoriensis','Morelia boeleni'], edge_length=0.017245, taxon_label=None),
            'Node5069296' : NodeRelationship(parent_label='Node5068976', child_labels=['Morelia viridis','Node5069392'], edge_length=0.004698, taxon_label=None),
            'Node5068976' : NodeRelationship(parent_label='Node5068912', child_labels=['Antaresia stimsoni','Node5069296'], edge_length=0.003603, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5069424', child_labels=[], edge_length=0.075494, taxon_label='Antaresia maculosa'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node5069424', child_labels=[], edge_length=0.063757, taxon_label='Morelia oenpelliensis'),
            'Node5069424' : NodeRelationship(parent_label='Node5068912', child_labels=['Antaresia maculosa','Morelia oenpelliensis'], edge_length=0.017571, taxon_label=None),
            'Node5068912' : NodeRelationship(parent_label='Node5068656', child_labels=['Node5068976','Node5069424'], edge_length=0.000386, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node5068656', child_labels=[], edge_length=0.070614, taxon_label='Morelia bredli'),
            'Node5068656' : NodeRelationship(parent_label='Node5068880', child_labels=['Node5068912','Morelia bredli'], edge_length=0.00713, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node5068880', child_labels=[], edge_length=0.079534, taxon_label='Morelia carinata'),
            'Node5068880' : NodeRelationship(parent_label='Node5068816', child_labels=['Node5068656','Morelia carinata'], edge_length=0.010249, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node5068816', child_labels=[], edge_length=0.114381, taxon_label='Python brongersmai'),
            'Node5068816' : NodeRelationship(parent_label='Node4840688', child_labels=['Node5068880','Python brongersmai'], edge_length=0.003672, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4840688', child_labels=[], edge_length=0.080896, taxon_label='Antaresia perthensis'),
            'Node4840688' : NodeRelationship(parent_label='Node5068464', child_labels=['Node5068816','Antaresia perthensis'], edge_length=0.003212, taxon_label=None),
            'Node5068464' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node5068592','Node4840688'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node5070000', child_labels=[], edge_length=0.05678, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node5070448', child_labels=[], edge_length=0.065928, taxon_label='Bothrochilus boa'),
            'Python timoriensis' : NodeRelationship(parent_label='Node5070448', child_labels=[], edge_length=0.10077, taxon_label='Python timoriensis'),
            'Node5070448' : NodeRelationship(parent_label='Node5070384', child_labels=['Bothrochilus boa','Python timoriensis'], edge_length=0.013813, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node5070384', child_labels=[], edge_length=0.083335, taxon_label='Antaresia perthensis'),
            'Node5070384' : NodeRelationship(parent_label='Node5070320', child_labels=['Node5070448','Antaresia perthensis'], edge_length=0.004377, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5070640', child_labels=[], edge_length=0.073392, taxon_label='Antaresia maculosa'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node5070640', child_labels=[], edge_length=0.065748, taxon_label='Morelia oenpelliensis'),
            'Node5070640' : NodeRelationship(parent_label='Node5070320', child_labels=['Antaresia maculosa','Morelia oenpelliensis'], edge_length=0.016433, taxon_label=None),
            'Node5070320' : NodeRelationship(parent_label='Node5070256', child_labels=['Node5070384','Node5070640'], edge_length=0.005104, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node5070256', child_labels=[], edge_length=0.07646, taxon_label='Morelia carinata'),
            'Node5070256' : NodeRelationship(parent_label='Node5069808', child_labels=['Node5070320','Morelia carinata'], edge_length=0.001866, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node5069808', child_labels=[], edge_length=0.06129, taxon_label='Antaresia stimsoni'),
            'Node5069808' : NodeRelationship(parent_label='Node5070160', child_labels=['Node5070256','Antaresia stimsoni'], edge_length=0.012314, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node5070160', child_labels=[], edge_length=0.081982, taxon_label='Morelia boeleni'),
            'Node5070160' : NodeRelationship(parent_label='Node5069968', child_labels=['Node5069808','Morelia boeleni'], edge_length=0.002833, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node5124176', child_labels=[], edge_length=0.064242, taxon_label='Morelia viridis'),
            'Python brongersmai' : NodeRelationship(parent_label='Node5124176', child_labels=[], edge_length=0.109434, taxon_label='Python brongersmai'),
            'Node5124176' : NodeRelationship(parent_label='Node5069968', child_labels=['Morelia viridis','Python brongersmai'], edge_length=0.010734, taxon_label=None),
            'Node5069968' : NodeRelationship(parent_label='Node5070128', child_labels=['Node5070160','Node5124176'], edge_length=0.00182, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node5070128', child_labels=[], edge_length=0.079156, taxon_label='Liasis fuscus'),
            'Node5070128' : NodeRelationship(parent_label='Node5070000', child_labels=['Node5069968','Liasis fuscus'], edge_length=0.014479, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node5070000', child_labels=[], edge_length=0.060374, taxon_label='Morelia bredli'),
            'Node5070000' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node5070128','Morelia bredli'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node5124944', child_labels=[], edge_length=0.063004, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node5125072', child_labels=[], edge_length=0.061413, taxon_label='Bothrochilus boa'),
            'Python brongersmai' : NodeRelationship(parent_label='Node5125072', child_labels=[], edge_length=0.099507, taxon_label='Python brongersmai'),
            'Node5125072' : NodeRelationship(parent_label='Node5124944', child_labels=['Bothrochilus boa','Python brongersmai'], edge_length=0.018059, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node5125456', child_labels=[], edge_length=0.075877, taxon_label='Liasis fuscus'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5125456', child_labels=[], edge_length=0.071375, taxon_label='Antaresia maculosa'),
            'Node5125456' : NodeRelationship(parent_label='Node5125392', child_labels=['Liasis fuscus','Antaresia maculosa'], edge_length=0.014459, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node5125392', child_labels=[], edge_length=0.054698, taxon_label='Antaresia stimsoni'),
            'Node5125392' : NodeRelationship(parent_label='Node5125136', child_labels=['Node5125456','Antaresia stimsoni'], edge_length=0.002182, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node5125648', child_labels=[], edge_length=0.062291, taxon_label='Morelia bredli'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node5126032', child_labels=[], edge_length=0.074463, taxon_label='Antaresia perthensis'),
            'Python timoriensis' : NodeRelationship(parent_label='Node5126032', child_labels=[], edge_length=0.103887, taxon_label='Python timoriensis'),
            'Node5126032' : NodeRelationship(parent_label='Node5125648', child_labels=['Antaresia perthensis','Python timoriensis'], edge_length=0.00521, taxon_label=None),
            'Node5125648' : NodeRelationship(parent_label='Node5125136', child_labels=['Morelia bredli','Node5126032'], edge_length=0.011197, taxon_label=None),
            'Node5125136' : NodeRelationship(parent_label='Node5125360', child_labels=['Node5125392','Node5125648'], edge_length=0.010763, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node5125360', child_labels=[], edge_length=0.070762, taxon_label='Morelia carinata'),
            'Node5125360' : NodeRelationship(parent_label='Node5125296', child_labels=['Node5125136','Morelia carinata'], edge_length=0.016533, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node5125296', child_labels=[], edge_length=0.068818, taxon_label='Morelia oenpelliensis'),
            'Node5125296' : NodeRelationship(parent_label='Node4796112', child_labels=['Node5125360','Morelia oenpelliensis'], edge_length=0.003714, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node5125584', child_labels=[], edge_length=0.071307, taxon_label='Morelia viridis'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node5125584', child_labels=[], edge_length=0.073392, taxon_label='Morelia boeleni'),
            'Node5125584' : NodeRelationship(parent_label='Node4796112', child_labels=['Morelia viridis','Morelia boeleni'], edge_length=0.010478, taxon_label=None),
            'Node4796112' : NodeRelationship(parent_label='Node5124944', child_labels=['Node5125296','Node5125584'], edge_length=0.00238, taxon_label=None),
            'Node5124944' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node5125072','Node4796112'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node5126608', child_labels=[], edge_length=0.053901, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node5126384', child_labels=[], edge_length=0.07927, taxon_label='Bothrochilus boa'),
            'Morelia bredli' : NodeRelationship(parent_label='Node5127120', child_labels=[], edge_length=0.053956, taxon_label='Morelia bredli'),
            'Python timoriensis' : NodeRelationship(parent_label='Node5127120', child_labels=[], edge_length=0.098874, taxon_label='Python timoriensis'),
            'Node5127120' : NodeRelationship(parent_label='Node5126384', child_labels=['Morelia bredli','Python timoriensis'], edge_length=0.016769, taxon_label=None),
            'Node5126384' : NodeRelationship(parent_label='Node5126768', child_labels=['Bothrochilus boa','Node5127120'], edge_length=0.003041, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node5126320', child_labels=[], edge_length=0.075223, taxon_label='Liasis fuscus'),
            'Morelia viridis' : NodeRelationship(parent_label='Node5126320', child_labels=[], edge_length=0.065805, taxon_label='Morelia viridis'),
            'Node5126320' : NodeRelationship(parent_label='Node5069904', child_labels=['Liasis fuscus','Morelia viridis'], edge_length=0.010712, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node5069904', child_labels=[], edge_length=0.075934, taxon_label='Morelia carinata'),
            'Node5069904' : NodeRelationship(parent_label='Node5126768', child_labels=['Node5126320','Morelia carinata'], edge_length=0.008083, taxon_label=None),
            'Node5126768' : NodeRelationship(parent_label='Node5126256', child_labels=['Node5126384','Node5069904'], edge_length=0.00777, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node5127472', child_labels=[], edge_length=0.082425, taxon_label='Antaresia perthensis'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5127280', child_labels=[], edge_length=0.075429, taxon_label='Antaresia maculosa'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node5127280', child_labels=[], edge_length=0.0787, taxon_label='Morelia boeleni'),
            'Node5127280' : NodeRelationship(parent_label='Node5127568', child_labels=['Antaresia maculosa','Morelia boeleni'], edge_length=0.010583, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node5127568', child_labels=[], edge_length=0.070289, taxon_label='Morelia oenpelliensis'),
            'Node5127568' : NodeRelationship(parent_label='Node5127472', child_labels=['Node5127280','Morelia oenpelliensis'], edge_length=0.008508, taxon_label=None),
            'Node5127472' : NodeRelationship(parent_label='Node5126256', child_labels=['Antaresia perthensis','Node5127568'], edge_length=0.00373, taxon_label=None),
            'Node5126256' : NodeRelationship(parent_label='Node5126736', child_labels=['Node5126768','Node5127472'], edge_length=0.002034, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node5126736', child_labels=[], edge_length=0.064652, taxon_label='Antaresia stimsoni'),
            'Node5126736' : NodeRelationship(parent_label='Node5126608', child_labels=['Node5126256','Antaresia stimsoni'], edge_length=0.013947, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node5126608', child_labels=[], edge_length=0.107868, taxon_label='Python brongersmai'),
            'Node5126608' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node5126736','Python brongersmai'], edge_length=None, taxon_label=None),
        },
    ]
    return treelist_node_references

def reference_tree_list():
    reference_tree_list = dendropy.TreeList(label=None, oid="TreeList4250528")
    tax_4795760 = reference_tree_list.taxon_set.require_taxon(label="Aspidites ramsayi", oid="Taxon4795760")
    tax_4837744 = reference_tree_list.taxon_set.require_taxon(label="Bothrochilus boa", oid="Taxon4837744")
    tax_4796368 = reference_tree_list.taxon_set.require_taxon(label="Liasis fuscus", oid="Taxon4796368")
    tax_4838128 = reference_tree_list.taxon_set.require_taxon(label="Morelia boeleni", oid="Taxon4838128")
    tax_4838224 = reference_tree_list.taxon_set.require_taxon(label="Morelia viridis", oid="Taxon4838224")
    tax_4838288 = reference_tree_list.taxon_set.require_taxon(label="Antaresia maculosa", oid="Taxon4838288")
    tax_4838160 = reference_tree_list.taxon_set.require_taxon(label="Antaresia stimsoni", oid="Taxon4838160")
    tax_4838448 = reference_tree_list.taxon_set.require_taxon(label="Python timoriensis", oid="Taxon4838448")
    tax_4838512 = reference_tree_list.taxon_set.require_taxon(label="Morelia oenpelliensis", oid="Taxon4838512")
    tax_4838672 = reference_tree_list.taxon_set.require_taxon(label="Morelia bredli", oid="Taxon4838672")
    tax_4838832 = reference_tree_list.taxon_set.require_taxon(label="Antaresia perthensis", oid="Taxon4838832")
    tax_4838960 = reference_tree_list.taxon_set.require_taxon(label="Python brongersmai", oid="Taxon4838960")
    tax_4839056 = reference_tree_list.taxon_set.require_taxon(label="Morelia carinata", oid="Taxon4839056")
    tree_4796048 = dendropy.Tree(label="Tree01", taxon_set=reference_tree_list.taxon_set, oid="Tree4796048")
    reference_tree_list.append(tree_4796048, reindex_taxa=False)
    nd_4796208 = tree_4796048.seed_node.new_child(label="Node4796144", taxon=tax_4795760, edge_length=0.056823, oid="Node4796208")
    nd_4796208.edge.oid = "Edge4796240"
    nd_4796272 = tree_4796048.seed_node.new_child(label="Node4796144", taxon=None, edge_length=0.011159, oid="Node4796272")
    nd_4796272.edge.oid = "Edge4795888"
    nd_4839088 = tree_4796048.seed_node.new_child(label="Node4796144", taxon=tax_4839056, edge_length=0.079321, oid="Node4839088")
    nd_4839088.edge.oid = "Edge4838704"
    nd_4796336 = nd_4796272.new_child(label="Node4796272", taxon=None, edge_length=0.003495, oid="Node4796336")
    nd_4796336.edge.oid = "Edge4796400"
    nd_4838576 = nd_4796272.new_child(label="Node4796272", taxon=None, edge_length=0.018613, oid="Node4838576")
    nd_4838576.edge.oid = "Edge4838736"
    nd_4796304 = nd_4796336.new_child(label="Node4796336", taxon=None, edge_length=0.007655, oid="Node4796304")
    nd_4796304.edge.oid = "Edge4837456"
    nd_4838640 = nd_4796336.new_child(label="Node4796336", taxon=tax_4838672, edge_length=0.065046, oid="Node4838640")
    nd_4838640.edge.oid = "Edge4838608"
    nd_4796016 = nd_4796304.new_child(label="Node4796304", taxon=None, edge_length=0.002264, oid="Node4796016")
    nd_4796016.edge.oid = "Edge4837520"
    nd_4838384 = nd_4796304.new_child(label="Node4796304", taxon=tax_4838512, edge_length=0.074047, oid="Node4838384")
    nd_4838384.edge.oid = "Edge4838544"
    nd_4837424 = nd_4796016.new_child(label="Node4796016", taxon=None, edge_length=0.0, oid="Node4837424")
    nd_4837424.edge.oid = "Edge4837584"
    nd_4838320 = nd_4796016.new_child(label="Node4796016", taxon=tax_4838448, edge_length=0.113605, oid="Node4838320")
    nd_4838320.edge.oid = "Edge4838480"
    nd_4837488 = nd_4837424.new_child(label="Node4837424", taxon=None, edge_length=0.007681, oid="Node4837488")
    nd_4837488.edge.oid = "Edge4837648"
    nd_4838192 = nd_4837424.new_child(label="Node4837424", taxon=tax_4838160, edge_length=0.061996, oid="Node4838192")
    nd_4838192.edge.oid = "Edge4838352"
    nd_4837552 = nd_4837488.new_child(label="Node4837488", taxon=None, edge_length=0.010103, oid="Node4837552")
    nd_4837552.edge.oid = "Edge4837712"
    nd_4838032 = nd_4837488.new_child(label="Node4837488", taxon=tax_4838288, edge_length=0.080749, oid="Node4838032")
    nd_4838032.edge.oid = "Edge4838096"
    nd_4837616 = nd_4837552.new_child(label="Node4837552", taxon=None, edge_length=0.011155, oid="Node4837616")
    nd_4837616.edge.oid = "Edge4837776"
    nd_4838256 = nd_4837552.new_child(label="Node4837552", taxon=tax_4838224, edge_length=0.066566, oid="Node4838256")
    nd_4838256.edge.oid = "Edge4838000"
    nd_4837680 = nd_4837616.new_child(label="Node4837616", taxon=tax_4837744, edge_length=0.073578, oid="Node4837680")
    nd_4837680.edge.oid = "Edge4837840"
    nd_4837936 = nd_4837616.new_child(label="Node4837616", taxon=None, edge_length=0.003063, oid="Node4837936")
    nd_4837936.edge.oid = "Edge4837808"
    nd_4837968 = nd_4837936.new_child(label="Node4837936", taxon=tax_4796368, edge_length=0.074347, oid="Node4837968")
    nd_4837968.edge.oid = "Edge4837872"
    nd_4838064 = nd_4837936.new_child(label="Node4837936", taxon=tax_4838128, edge_length=0.07628, oid="Node4838064")
    nd_4838064.edge.oid = "Edge4837904"
    nd_4838416 = nd_4838576.new_child(label="Node4838576", taxon=tax_4838832, edge_length=0.065004, oid="Node4838416")
    nd_4838416.edge.oid = "Edge4838768"
    nd_4838928 = nd_4838576.new_child(label="Node4838576", taxon=tax_4838960, edge_length=0.107706, oid="Node4838928")
    nd_4838928.edge.oid = "Edge4838800"
    tree_4838896 = dendropy.Tree(label="Tree02", taxon_set=reference_tree_list.taxon_set, oid="Tree4838896")
    reference_tree_list.append(tree_4838896, reindex_taxa=False)
    nd_4839280 = tree_4838896.seed_node.new_child(label="Node4839216", taxon=tax_4795760, edge_length=0.065297, oid="Node4839280")
    nd_4839280.edge.oid = "Edge4839312"
    nd_4839344 = tree_4838896.seed_node.new_child(label="Node4839216", taxon=None, edge_length=0.001444, oid="Node4839344")
    nd_4839344.edge.oid = "Edge4839184"
    nd_4839696 = tree_4838896.seed_node.new_child(label="Node4839216", taxon=None, edge_length=0.004077, oid="Node4839696")
    nd_4839696.edge.oid = "Edge4839632"
    nd_4839152 = nd_4839344.new_child(label="Node4839344", taxon=None, edge_length=0.010946, oid="Node4839152")
    nd_4839152.edge.oid = "Edge4839440"
    nd_4839408 = nd_4839344.new_child(label="Node4839344", taxon=tax_4838832, edge_length=0.079424, oid="Node4839408")
    nd_4839408.edge.oid = "Edge4839856"
    nd_4795856 = nd_4839152.new_child(label="Node4839152", taxon=None, edge_length=0.016022, oid="Node4795856")
    nd_4795856.edge.oid = "Edge4839504"
    nd_4839888 = nd_4839152.new_child(label="Node4839152", taxon=tax_4838224, edge_length=0.066015, oid="Node4839888")
    nd_4839888.edge.oid = "Edge4839792"
    nd_4839376 = nd_4795856.new_child(label="Node4795856", taxon=tax_4837744, edge_length=0.073565, oid="Node4839376")
    nd_4839376.edge.oid = "Edge4839568"
    nd_4839664 = nd_4795856.new_child(label="Node4795856", taxon=None, edge_length=0.007879, oid="Node4839664")
    nd_4839664.edge.oid = "Edge4839600"
    nd_4839536 = nd_4839664.new_child(label="Node4839664", taxon=tax_4838128, edge_length=0.069885, oid="Node4839536")
    nd_4839536.edge.oid = "Edge4839472"
    nd_4839760 = nd_4839664.new_child(label="Node4839664", taxon=tax_4838512, edge_length=0.06247, oid="Node4839760")
    nd_4839760.edge.oid = "Edge4839728"
    nd_4839824 = nd_4839696.new_child(label="Node4839696", taxon=None, edge_length=0.008007, oid="Node4839824")
    nd_4839824.edge.oid = "Edge4839984"
    nd_4840496 = nd_4839696.new_child(label="Node4839696", taxon=tax_4838672, edge_length=0.063725, oid="Node4840496")
    nd_4840496.edge.oid = "Edge4840528"
    nd_4840016 = nd_4839824.new_child(label="Node4839824", taxon=tax_4796368, edge_length=0.082049, oid="Node4840016")
    nd_4840016.edge.oid = "Edge4840080"
    nd_4840144 = nd_4839824.new_child(label="Node4839824", taxon=None, edge_length=0.006907, oid="Node4840144")
    nd_4840144.edge.oid = "Edge4840112"
    nd_4840048 = nd_4840144.new_child(label="Node4840144", taxon=tax_4838160, edge_length=0.056894, oid="Node4840048")
    nd_4840048.edge.oid = "Edge4839952"
    nd_4840304 = nd_4840144.new_child(label="Node4840144", taxon=None, edge_length=0.004552, oid="Node4840304")
    nd_4840304.edge.oid = "Edge4796080"
    nd_4840176 = nd_4840304.new_child(label="Node4840304", taxon=None, edge_length=0.010752, oid="Node4840176")
    nd_4840176.edge.oid = "Edge4838992"
    nd_4840272 = nd_4840304.new_child(label="Node4840304", taxon=tax_4838288, edge_length=0.079792, oid="Node4840272")
    nd_4840272.edge.oid = "Edge4840336"
    nd_4840208 = nd_4840176.new_child(label="Node4840176", taxon=None, edge_length=0.015909, oid="Node4840208")
    nd_4840208.edge.oid = "Edge4840240"
    nd_4840464 = nd_4840176.new_child(label="Node4840176", taxon=tax_4838960, edge_length=0.114385, oid="Node4840464")
    nd_4840464.edge.oid = "Edge4840400"
    nd_4795984 = nd_4840208.new_child(label="Node4840208", taxon=tax_4838448, edge_length=0.092077, oid="Node4795984")
    nd_4795984.edge.oid = "Edge4839920"
    nd_4840432 = nd_4840208.new_child(label="Node4840208", taxon=tax_4839056, edge_length=0.0579, oid="Node4840432")
    nd_4840432.edge.oid = "Edge4840368"
    tree_4840592 = dendropy.Tree(label="Tree03", taxon_set=reference_tree_list.taxon_set, oid="Tree4840592")
    reference_tree_list.append(tree_4840592, reindex_taxa=False)
    nd_4840816 = tree_4840592.seed_node.new_child(label="Node4840752", taxon=tax_4795760, edge_length=0.06111, oid="Node4840816")
    nd_4840816.edge.oid = "Edge4840848"
    nd_4840880 = tree_4840592.seed_node.new_child(label="Node4840752", taxon=None, edge_length=0.005408, oid="Node4840880")
    nd_4840880.edge.oid = "Edge4838864"
    nd_4895344 = tree_4840592.seed_node.new_child(label="Node4840752", taxon=None, edge_length=0.012121, oid="Node4895344")
    nd_4895344.edge.oid = "Edge4895216"
    nd_4840720 = nd_4840880.new_child(label="Node4840880", taxon=None, edge_length=0.006988, oid="Node4840720")
    nd_4840720.edge.oid = "Edge4840976"
    nd_4895184 = nd_4840880.new_child(label="Node4840880", taxon=tax_4838160, edge_length=0.063808, oid="Node4895184")
    nd_4895184.edge.oid = "Edge4895312"
    nd_4840912 = nd_4840720.new_child(label="Node4840720", taxon=None, edge_length=0.005167, oid="Node4840912")
    nd_4840912.edge.oid = "Edge4841040"
    nd_4895088 = nd_4840720.new_child(label="Node4840720", taxon=None, edge_length=0.00574, oid="Node4895088")
    nd_4895088.edge.oid = "Edge4894864"
    nd_4840560 = nd_4840912.new_child(label="Node4840912", taxon=None, edge_length=0.010367, oid="Node4840560")
    nd_4840560.edge.oid = "Edge4841104"
    nd_4841200 = nd_4840912.new_child(label="Node4840912", taxon=None, edge_length=0.007731, oid="Node4841200")
    nd_4841200.edge.oid = "Edge4840944"
    nd_4841008 = nd_4840560.new_child(label="Node4840560", taxon=None, edge_length=0.008688, oid="Node4841008")
    nd_4841008.edge.oid = "Edge4841168"
    nd_4841360 = nd_4840560.new_child(label="Node4840560", taxon=tax_4838960, edge_length=0.1069, oid="Node4841360")
    nd_4841360.edge.oid = "Edge4841296"
    nd_4841072 = nd_4841008.new_child(label="Node4841008", taxon=tax_4837744, edge_length=0.067908, oid="Node4841072")
    nd_4841072.edge.oid = "Edge4841232"
    nd_4841328 = nd_4841008.new_child(label="Node4841008", taxon=tax_4838832, edge_length=0.077196, oid="Node4841328")
    nd_4841328.edge.oid = "Edge4841264"
    nd_4841392 = nd_4841200.new_child(label="Node4841200", taxon=None, edge_length=0.004129, oid="Node4841392")
    nd_4841392.edge.oid = "Edge4894768"
    nd_4841136 = nd_4841200.new_child(label="Node4841200", taxon=tax_4838512, edge_length=0.068424, oid="Node4841136")
    nd_4841136.edge.oid = "Edge4894800"
    nd_4841424 = nd_4841392.new_child(label="Node4841392", taxon=tax_4838448, edge_length=0.109318, oid="Node4841424")
    nd_4841424.edge.oid = "Edge4894832"
    nd_4841456 = nd_4841392.new_child(label="Node4841392", taxon=tax_4838288, edge_length=0.081382, oid="Node4841456")
    nd_4841456.edge.oid = "Edge4894928"
    nd_4895024 = nd_4895088.new_child(label="Node4895088", taxon=tax_4796368, edge_length=0.079478, oid="Node4895024")
    nd_4895024.edge.oid = "Edge4894896"
    nd_4895152 = nd_4895088.new_child(label="Node4895088", taxon=None, edge_length=0.013155, oid="Node4895152")
    nd_4895152.edge.oid = "Edge4895120"
    nd_4895056 = nd_4895152.new_child(label="Node4895152", taxon=tax_4838224, edge_length=0.062239, oid="Node4895056")
    nd_4895056.edge.oid = "Edge4894960"
    nd_4895280 = nd_4895152.new_child(label="Node4895152", taxon=tax_4838672, edge_length=0.060484, oid="Node4895280")
    nd_4895280.edge.oid = "Edge4895248"
    nd_4895408 = nd_4895344.new_child(label="Node4895344", taxon=tax_4839056, edge_length=0.069501, oid="Node4895408")
    nd_4895408.edge.oid = "Edge4895440"
    nd_4895568 = nd_4895344.new_child(label="Node4895344", taxon=tax_4838128, edge_length=0.073602, oid="Node4895568")
    nd_4895568.edge.oid = "Edge4895504"
    tree_4895536 = dendropy.Tree(label="Tree04", taxon_set=reference_tree_list.taxon_set, oid="Tree4895536")
    reference_tree_list.append(tree_4895536, reindex_taxa=False)
    nd_4895760 = tree_4895536.seed_node.new_child(label="Node4895696", taxon=tax_4795760, edge_length=0.053028, oid="Node4895760")
    nd_4895760.edge.oid = "Edge4895792"
    nd_4895824 = tree_4895536.seed_node.new_child(label="Node4895696", taxon=None, edge_length=0.014376, oid="Node4895824")
    nd_4895824.edge.oid = "Edge4895376"
    nd_4897264 = tree_4895536.seed_node.new_child(label="Node4895696", taxon=tax_4838128, edge_length=0.076894, oid="Node4897264")
    nd_4897264.edge.oid = "Edge4896752"
    nd_4895632 = nd_4895824.new_child(label="Node4895824", taxon=None, edge_length=0.002567, oid="Node4895632")
    nd_4895632.edge.oid = "Edge4895920"
    nd_4896016 = nd_4895824.new_child(label="Node4895824", taxon=None, edge_length=0.005969, oid="Node4896016")
    nd_4896016.edge.oid = "Edge4895888"
    nd_4795920 = nd_4895632.new_child(label="Node4895632", taxon=None, edge_length=0.012413, oid="Node4795920")
    nd_4795920.edge.oid = "Edge4895984"
    nd_4896176 = nd_4895632.new_child(label="Node4895632", taxon=tax_4838672, edge_length=0.069933, oid="Node4896176")
    nd_4896176.edge.oid = "Edge4896112"
    nd_4895856 = nd_4795920.new_child(label="Node4795920", taxon=tax_4837744, edge_length=0.072831, oid="Node4895856")
    nd_4895856.edge.oid = "Edge4896048"
    nd_4896144 = nd_4795920.new_child(label="Node4795920", taxon=tax_4838832, edge_length=0.073004, oid="Node4896144")
    nd_4896144.edge.oid = "Edge4896080"
    nd_4896208 = nd_4896016.new_child(label="Node4896016", taxon=None, edge_length=0.004803, oid="Node4896208")
    nd_4896208.edge.oid = "Edge4896272"
    nd_4896848 = nd_4896016.new_child(label="Node4896016", taxon=tax_4838160, edge_length=0.064721, oid="Node4896848")
    nd_4896848.edge.oid = "Edge4897040"
    nd_4896304 = nd_4896208.new_child(label="Node4896208", taxon=None, edge_length=5.134e-05, oid="Node4896304")
    nd_4896304.edge.oid = "Edge4896368"
    nd_4897136 = nd_4896208.new_child(label="Node4896208", taxon=tax_4838960, edge_length=0.118976, oid="Node4897136")
    nd_4897136.edge.oid = "Edge4897072"
    nd_4896240 = nd_4896304.new_child(label="Node4896304", taxon=None, edge_length=0.011519, oid="Node4896240")
    nd_4896240.edge.oid = "Edge4896432"
    nd_4896528 = nd_4896304.new_child(label="Node4896304", taxon=None, edge_length=0.007647, oid="Node4896528")
    nd_4896528.edge.oid = "Edge4895952"
    nd_4896336 = nd_4896240.new_child(label="Node4896240", taxon=None, edge_length=0.011723, oid="Node4896336")
    nd_4896336.edge.oid = "Edge4896496"
    nd_4896688 = nd_4896240.new_child(label="Node4896240", taxon=tax_4838224, edge_length=0.065297, oid="Node4896688")
    nd_4896688.edge.oid = "Edge4896656"
    nd_4896400 = nd_4896336.new_child(label="Node4896336", taxon=tax_4796368, edge_length=0.072339, oid="Node4896400")
    nd_4896400.edge.oid = "Edge4896560"
    nd_4896624 = nd_4896336.new_child(label="Node4896336", taxon=tax_4838512, edge_length=0.062058, oid="Node4896624")
    nd_4896624.edge.oid = "Edge4896592"
    nd_4896784 = nd_4896528.new_child(label="Node4896528", taxon=None, edge_length=0.022912, oid="Node4896784")
    nd_4896784.edge.oid = "Edge4896464"
    nd_4897008 = nd_4896528.new_child(label="Node4896528", taxon=tax_4838288, edge_length=0.079945, oid="Node4897008")
    nd_4897008.edge.oid = "Edge4896944"
    nd_4896816 = nd_4896784.new_child(label="Node4896784", taxon=tax_4838448, edge_length=0.093762, oid="Node4896816")
    nd_4896816.edge.oid = "Edge4896880"
    nd_4896976 = nd_4896784.new_child(label="Node4896784", taxon=tax_4839056, edge_length=0.055992, oid="Node4896976")
    nd_4896976.edge.oid = "Edge4896912"
    tree_4897104 = dendropy.Tree(label="Tree05", taxon_set=reference_tree_list.taxon_set, oid="Tree4897104")
    reference_tree_list.append(tree_4897104, reindex_taxa=False)
    nd_4897424 = tree_4897104.seed_node.new_child(label="Node4897360", taxon=tax_4795760, edge_length=0.06391, oid="Node4897424")
    nd_4897424.edge.oid = "Edge4897456"
    nd_4897488 = tree_4897104.seed_node.new_child(label="Node4897360", taxon=None, edge_length=0.00451, oid="Node4897488")
    nd_4897488.edge.oid = "Edge4895472"
    nd_4897520 = tree_4897104.seed_node.new_child(label="Node4897360", taxon=None, edge_length=0.005404, oid="Node4897520")
    nd_4897520.edge.oid = "Edge4897648"
    nd_4897232 = nd_4897488.new_child(label="Node4897488", taxon=tax_4837744, edge_length=0.071835, oid="Node4897232")
    nd_4897232.edge.oid = "Edge4897584"
    nd_4897680 = nd_4897488.new_child(label="Node4897488", taxon=None, edge_length=0.010888, oid="Node4897680")
    nd_4897680.edge.oid = "Edge4897616"
    nd_4897328 = nd_4897680.new_child(label="Node4897680", taxon=tax_4796368, edge_length=0.075788, oid="Node4897328")
    nd_4897328.edge.oid = "Edge4840624"
    nd_4897296 = nd_4897680.new_child(label="Node4897680", taxon=tax_4839056, edge_length=0.07116, oid="Node4897296")
    nd_4897296.edge.oid = "Edge4897552"
    nd_4897712 = nd_4897520.new_child(label="Node4897520", taxon=None, edge_length=0.011944, oid="Node4897712")
    nd_4897712.edge.oid = "Edge4894992"
    nd_4898608 = nd_4897520.new_child(label="Node4897520", taxon=None, edge_length=0.00612, oid="Node4898608")
    nd_4898608.edge.oid = "Edge4898512"
    nd_4840656 = nd_4897712.new_child(label="Node4897712", taxon=None, edge_length=0.013817, oid="Node4840656")
    nd_4840656.edge.oid = "Edge4897776"
    nd_4897872 = nd_4897712.new_child(label="Node4897712", taxon=None, edge_length=0.012897, oid="Node4897872")
    nd_4897872.edge.oid = "Edge4895664"
    nd_4896720 = nd_4840656.new_child(label="Node4840656", taxon=None, edge_length=0.012356, oid="Node4896720")
    nd_4896720.edge.oid = "Edge4897840"
    nd_4898032 = nd_4840656.new_child(label="Node4840656", taxon=tax_4838288, edge_length=0.077058, oid="Node4898032")
    nd_4898032.edge.oid = "Edge4897968"
    nd_4897744 = nd_4896720.new_child(label="Node4896720", taxon=tax_4838160, edge_length=0.048754, oid="Node4897744")
    nd_4897744.edge.oid = "Edge4897904"
    nd_4898000 = nd_4896720.new_child(label="Node4896720", taxon=tax_4838832, edge_length=0.069682, oid="Node4898000")
    nd_4898000.edge.oid = "Edge4897936"
    nd_4898064 = nd_4897872.new_child(label="Node4897872", taxon=None, edge_length=0.003147, oid="Node4898064")
    nd_4898064.edge.oid = "Edge4898128"
    nd_4898544 = nd_4897872.new_child(label="Node4897872", taxon=tax_4838448, edge_length=0.108243, oid="Node4898544")
    nd_4898544.edge.oid = "Edge4898416"
    nd_4898160 = nd_4898064.new_child(label="Node4898064", taxon=tax_4838224, edge_length=0.07026, oid="Node4898160")
    nd_4898160.edge.oid = "Edge4898224"
    nd_4898288 = nd_4898064.new_child(label="Node4898064", taxon=None, edge_length=0.001443, oid="Node4898288")
    nd_4898288.edge.oid = "Edge4898256"
    nd_4898192 = nd_4898288.new_child(label="Node4898288", taxon=tax_4838672, edge_length=0.066182, oid="Node4898192")
    nd_4898192.edge.oid = "Edge4898096"
    nd_4898448 = nd_4898288.new_child(label="Node4898288", taxon=tax_4838512, edge_length=0.068166, oid="Node4898448")
    nd_4898448.edge.oid = "Edge4898384"
    nd_4898320 = nd_4898608.new_child(label="Node4898608", taxon=tax_4838960, edge_length=0.108072, oid="Node4898320")
    nd_4898320.edge.oid = "Edge4897808"
    nd_4898704 = nd_4898608.new_child(label="Node4898608", taxon=tax_4838128, edge_length=0.073251, oid="Node4898704")
    nd_4898704.edge.oid = "Edge4898640"
    tree_4898672 = dendropy.Tree(label="Tree06", taxon_set=reference_tree_list.taxon_set, oid="Tree4898672")
    reference_tree_list.append(tree_4898672, reindex_taxa=False)
    nd_5066864 = tree_4898672.seed_node.new_child(label="Node5066800", taxon=tax_4795760, edge_length=0.06501, oid="Node5066864")
    nd_5066864.edge.oid = "Edge5066896"
    nd_4898480 = tree_4898672.seed_node.new_child(label="Node5066800", taxon=None, edge_length=0.003644, oid="Node4898480")
    nd_4898480.edge.oid = "Edge4898352"
    nd_5067600 = tree_4898672.seed_node.new_child(label="Node5066800", taxon=None, edge_length=0.008523, oid="Node5067600")
    nd_5067600.edge.oid = "Edge5067472"
    nd_5066960 = nd_4898480.new_child(label="Node4898480", taxon=None, edge_length=0.011722, oid="Node5066960")
    nd_5066960.edge.oid = "Edge5067024"
    nd_5067120 = nd_4898480.new_child(label="Node4898480", taxon=None, edge_length=0.014836, oid="Node5067120")
    nd_5067120.edge.oid = "Edge5066992"
    nd_4897168 = nd_5066960.new_child(label="Node5066960", taxon=None, edge_length=0.017824, oid="Node4897168")
    nd_4897168.edge.oid = "Edge5067088"
    nd_5067280 = nd_5066960.new_child(label="Node5066960", taxon=tax_4838448, edge_length=0.102639, oid="Node5067280")
    nd_5067280.edge.oid = "Edge5067216"
    nd_5066928 = nd_4897168.new_child(label="Node4897168", taxon=tax_4837744, edge_length=0.059984, oid="Node5066928")
    nd_5066928.edge.oid = "Edge5067152"
    nd_5067248 = nd_4897168.new_child(label="Node4897168", taxon=tax_4838960, edge_length=0.100555, oid="Node5067248")
    nd_5067248.edge.oid = "Edge5067184"
    nd_5067312 = nd_5067120.new_child(label="Node5067120", taxon=tax_4838672, edge_length=0.0538, oid="Node5067312")
    nd_5067312.edge.oid = "Edge5067376"
    nd_5067504 = nd_5067120.new_child(label="Node5067120", taxon=tax_4838288, edge_length=0.07843, oid="Node5067504")
    nd_5067504.edge.oid = "Edge5067440"
    nd_5067056 = nd_5067600.new_child(label="Node5067600", taxon=None, edge_length=0.0, oid="Node5067056")
    nd_5067056.edge.oid = "Edge5067536"
    nd_5068368 = nd_5067600.new_child(label="Node5067600", taxon=tax_4838224, edge_length=0.072091, oid="Node5068368")
    nd_5068368.edge.oid = "Edge5067696"
    nd_5067344 = nd_5067056.new_child(label="Node5067056", taxon=None, edge_length=0.015326, oid="Node5067344")
    nd_5067344.edge.oid = "Edge5067664"
    nd_5067856 = nd_5067056.new_child(label="Node5067056", taxon=None, edge_length=0.005897, oid="Node5067856")
    nd_5067856.edge.oid = "Edge5067824"
    nd_5067408 = nd_5067344.new_child(label="Node5067344", taxon=tax_4796368, edge_length=0.068438, oid="Node5067408")
    nd_5067408.edge.oid = "Edge5067728"
    nd_5067792 = nd_5067344.new_child(label="Node5067344", taxon=tax_4838512, edge_length=0.066815, oid="Node5067792")
    nd_5067792.edge.oid = "Edge5067760"
    nd_5067888 = nd_5067856.new_child(label="Node5067856", taxon=None, edge_length=0.010508, oid="Node5067888")
    nd_5067888.edge.oid = "Edge5067632"
    nd_5068304 = nd_5067856.new_child(label="Node5067856", taxon=tax_4838832, edge_length=0.083409, oid="Node5068304")
    nd_5068304.edge.oid = "Edge5068176"
    nd_5067920 = nd_5067888.new_child(label="Node5067888", taxon=tax_4838160, edge_length=0.057007, oid="Node5067920")
    nd_5067920.edge.oid = "Edge5067984"
    nd_5068080 = nd_5067888.new_child(label="Node5067888", taxon=None, edge_length=0.013392, oid="Node5068080")
    nd_5068080.edge.oid = "Edge5068016"
    nd_5067952 = nd_5068080.new_child(label="Node5068080", taxon=tax_4839056, edge_length=0.062766, oid="Node5067952")
    nd_5067952.edge.oid = "Edge5067568"
    nd_5068208 = nd_5068080.new_child(label="Node5068080", taxon=tax_4838128, edge_length=0.079145, oid="Node5068208")
    nd_5068208.edge.oid = "Edge5068144"
    tree_5068240 = dendropy.Tree(label="Tree07", taxon_set=reference_tree_list.taxon_set, oid="Tree5068240")
    reference_tree_list.append(tree_5068240, reindex_taxa=False)
    nd_5068528 = tree_5068240.seed_node.new_child(label="Node5068464", taxon=tax_4795760, edge_length=0.062297, oid="Node5068528")
    nd_5068528.edge.oid = "Edge5068560"
    nd_5068592 = tree_5068240.seed_node.new_child(label="Node5068464", taxon=None, edge_length=0.009899, oid="Node5068592")
    nd_5068592.edge.oid = "Edge5068432"
    nd_4840688 = tree_5068240.seed_node.new_child(label="Node5068464", taxon=None, edge_length=0.003212, oid="Node4840688")
    nd_4840688.edge.oid = "Edge5068624"
    nd_5068272 = nd_5068592.new_child(label="Node5068592", taxon=tax_4837744, edge_length=0.069165, oid="Node5068272")
    nd_5068272.edge.oid = "Edge5068688"
    nd_5068784 = nd_5068592.new_child(label="Node5068592", taxon=tax_4796368, edge_length=0.073495, oid="Node5068784")
    nd_5068784.edge.oid = "Edge5068720"
    nd_5068816 = nd_4840688.new_child(label="Node4840688", taxon=None, edge_length=0.003672, oid="Node5068816")
    nd_5068816.edge.oid = "Edge5068752"
    nd_5069456 = nd_4840688.new_child(label="Node4840688", taxon=tax_4838832, edge_length=0.080896, oid="Node5069456")
    nd_5069456.edge.oid = "Edge5068400"
    nd_5068880 = nd_5068816.new_child(label="Node5068816", taxon=None, edge_length=0.010249, oid="Node5068880")
    nd_5068880.edge.oid = "Edge5068944"
    nd_4898576 = nd_5068816.new_child(label="Node5068816", taxon=tax_4838960, edge_length=0.114381, oid="Node4898576")
    nd_4898576.edge.oid = "Edge5069328"
    nd_5068656 = nd_5068880.new_child(label="Node5068880", taxon=None, edge_length=0.00713, oid="Node5068656")
    nd_5068656.edge.oid = "Edge5069008"
    nd_4898768 = nd_5068880.new_child(label="Node5068880", taxon=tax_4839056, edge_length=0.079534, oid="Node4898768")
    nd_4898768.edge.oid = "Edge5069584"
    nd_5068912 = nd_5068656.new_child(label="Node5068656", taxon=None, edge_length=0.000386, oid="Node5068912")
    nd_5068912.edge.oid = "Edge5069072"
    nd_5069840 = nd_5068656.new_child(label="Node5068656", taxon=tax_4838672, edge_length=0.070614, oid="Node5069840")
    nd_5069840.edge.oid = "Edge5069712"
    nd_5068976 = nd_5068912.new_child(label="Node5068912", taxon=None, edge_length=0.003603, oid="Node5068976")
    nd_5068976.edge.oid = "Edge5069136"
    nd_5069424 = nd_5068912.new_child(label="Node5068912", taxon=None, edge_length=0.017571, oid="Node5069424")
    nd_5069424.edge.oid = "Edge5069520"
    nd_5069040 = nd_5068976.new_child(label="Node5068976", taxon=tax_4838160, edge_length=0.061022, oid="Node5069040")
    nd_5069040.edge.oid = "Edge5069200"
    nd_5069296 = nd_5068976.new_child(label="Node5068976", taxon=None, edge_length=0.004698, oid="Node5069296")
    nd_5069296.edge.oid = "Edge5069232"
    nd_5069168 = nd_5069296.new_child(label="Node5069296", taxon=tax_4838224, edge_length=0.071291, oid="Node5069168")
    nd_5069168.edge.oid = "Edge5069104"
    nd_5069392 = nd_5069296.new_child(label="Node5069296", taxon=None, edge_length=0.017245, oid="Node5069392")
    nd_5069392.edge.oid = "Edge5069360"
    nd_5069264 = nd_5069392.new_child(label="Node5069392", taxon=tax_4838448, edge_length=0.098776, oid="Node5069264")
    nd_5069264.edge.oid = "Edge5068848"
    nd_5069552 = nd_5069392.new_child(label="Node5069392", taxon=tax_4838128, edge_length=0.075069, oid="Node5069552")
    nd_5069552.edge.oid = "Edge5069488"
    nd_5069616 = nd_5069424.new_child(label="Node5069424", taxon=tax_4838288, edge_length=0.075494, oid="Node5069616")
    nd_5069616.edge.oid = "Edge5069648"
    nd_5069744 = nd_5069424.new_child(label="Node5069424", taxon=tax_4838512, edge_length=0.063757, oid="Node5069744")
    nd_5069744.edge.oid = "Edge5069680"
    tree_5069776 = dendropy.Tree(label="Tree08", taxon_set=reference_tree_list.taxon_set, oid="Tree5069776")
    reference_tree_list.append(tree_5069776, reindex_taxa=False)
    nd_5070064 = tree_5069776.seed_node.new_child(label="Node5070000", taxon=tax_4795760, edge_length=0.05678, oid="Node5070064")
    nd_5070064.edge.oid = "Edge5070096"
    nd_5070128 = tree_5069776.seed_node.new_child(label="Node5070000", taxon=None, edge_length=0.014479, oid="Node5070128")
    nd_5070128.edge.oid = "Edge5068048"
    nd_5124752 = tree_5069776.seed_node.new_child(label="Node5070000", taxon=tax_4838672, edge_length=0.060374, oid="Node5124752")
    nd_5124752.edge.oid = "Edge5124592"
    nd_5069968 = nd_5070128.new_child(label="Node5070128", taxon=None, edge_length=0.00182, oid="Node5069968")
    nd_5069968.edge.oid = "Edge5070224"
    nd_5124784 = nd_5070128.new_child(label="Node5070128", taxon=tax_4796368, edge_length=0.079156, oid="Node5124784")
    nd_5124784.edge.oid = "Edge5124688"
    nd_5070160 = nd_5069968.new_child(label="Node5069968", taxon=None, edge_length=0.002833, oid="Node5070160")
    nd_5070160.edge.oid = "Edge5070288"
    nd_5124176 = nd_5069968.new_child(label="Node5069968", taxon=None, edge_length=0.010734, oid="Node5124176")
    nd_5124176.edge.oid = "Edge5124464"
    nd_5069808 = nd_5070160.new_child(label="Node5070160", taxon=None, edge_length=0.012314, oid="Node5069808")
    nd_5069808.edge.oid = "Edge5070352"
    nd_5124336 = nd_5070160.new_child(label="Node5070160", taxon=tax_4838128, edge_length=0.081982, oid="Node5124336")
    nd_5124336.edge.oid = "Edge5124368"
    nd_5070256 = nd_5069808.new_child(label="Node5069808", taxon=None, edge_length=0.001866, oid="Node5070256")
    nd_5070256.edge.oid = "Edge5070416"
    nd_5070192 = nd_5069808.new_child(label="Node5069808", taxon=tax_4838160, edge_length=0.06129, oid="Node5070192")
    nd_5070192.edge.oid = "Edge5124272"
    nd_5070320 = nd_5070256.new_child(label="Node5070256", taxon=None, edge_length=0.005104, oid="Node5070320")
    nd_5070320.edge.oid = "Edge5070480"
    nd_5124400 = nd_5070256.new_child(label="Node5070256", taxon=tax_4839056, edge_length=0.07646, oid="Node5124400")
    nd_5124400.edge.oid = "Edge5124240"
    nd_5070384 = nd_5070320.new_child(label="Node5070320", taxon=None, edge_length=0.004377, oid="Node5070384")
    nd_5070384.edge.oid = "Edge5070544"
    nd_5070640 = nd_5070320.new_child(label="Node5070320", taxon=None, edge_length=0.016433, oid="Node5070640")
    nd_5070640.edge.oid = "Edge5070576"
    nd_5070448 = nd_5070384.new_child(label="Node5070384", taxon=None, edge_length=0.013813, oid="Node5070448")
    nd_5070448.edge.oid = "Edge5070608"
    nd_5070800 = nd_5070384.new_child(label="Node5070384", taxon=tax_4838832, edge_length=0.083335, oid="Node5070800")
    nd_5070800.edge.oid = "Edge5070736"
    nd_5070512 = nd_5070448.new_child(label="Node5070448", taxon=tax_4837744, edge_length=0.065928, oid="Node5070512")
    nd_5070512.edge.oid = "Edge5070672"
    nd_5070768 = nd_5070448.new_child(label="Node5070448", taxon=tax_4838448, edge_length=0.10077, oid="Node5070768")
    nd_5070768.edge.oid = "Edge5070704"
    nd_5124144 = nd_5070640.new_child(label="Node5070640", taxon=tax_4838288, edge_length=0.073392, oid="Node5124144")
    nd_5124144.edge.oid = "Edge5124208"
    nd_5070832 = nd_5070640.new_child(label="Node5070640", taxon=tax_4838512, edge_length=0.065748, oid="Node5070832")
    nd_5070832.edge.oid = "Edge5124304"
    nd_5124560 = nd_5124176.new_child(label="Node5124176", taxon=tax_4838224, edge_length=0.064242, oid="Node5124560")
    nd_5124560.edge.oid = "Edge5124496"
    nd_5124656 = nd_5124176.new_child(label="Node5124176", taxon=tax_4838960, edge_length=0.109434, oid="Node5124656")
    nd_5124656.edge.oid = "Edge5124624"
    tree_5124528 = dendropy.Tree(label="Tree09", taxon_set=reference_tree_list.taxon_set, oid="Tree5124528")
    reference_tree_list.append(tree_5124528, reindex_taxa=False)
    nd_5125008 = tree_5124528.seed_node.new_child(label="Node5124944", taxon=tax_4795760, edge_length=0.063004, oid="Node5125008")
    nd_5125008.edge.oid = "Edge5125040"
    nd_5125072 = tree_5124528.seed_node.new_child(label="Node5124944", taxon=None, edge_length=0.018059, oid="Node5125072")
    nd_5125072.edge.oid = "Edge5124912"
    nd_4796112 = tree_5124528.seed_node.new_child(label="Node5124944", taxon=None, edge_length=0.00238, oid="Node4796112")
    nd_4796112.edge.oid = "Edge5125104"
    nd_5124848 = nd_5125072.new_child(label="Node5125072", taxon=tax_4837744, edge_length=0.061413, oid="Node5124848")
    nd_5124848.edge.oid = "Edge5125168"
    nd_5125264 = nd_5125072.new_child(label="Node5125072", taxon=tax_4838960, edge_length=0.099507, oid="Node5125264")
    nd_5125264.edge.oid = "Edge5125200"
    nd_5125296 = nd_4796112.new_child(label="Node4796112", taxon=None, edge_length=0.003714, oid="Node5125296")
    nd_5125296.edge.oid = "Edge5125328"
    nd_5125584 = nd_4796112.new_child(label="Node4796112", taxon=None, edge_length=0.010478, oid="Node5125584")
    nd_5125584.edge.oid = "Edge5126224"
    nd_5125360 = nd_5125296.new_child(label="Node5125296", taxon=None, edge_length=0.016533, oid="Node5125360")
    nd_5125360.edge.oid = "Edge5125424"
    nd_5126192 = nd_5125296.new_child(label="Node5125296", taxon=tax_4838512, edge_length=0.068818, oid="Node5126192")
    nd_5126192.edge.oid = "Edge5126064"
    nd_5125136 = nd_5125360.new_child(label="Node5125360", taxon=None, edge_length=0.010763, oid="Node5125136")
    nd_5125136.edge.oid = "Edge5125488"
    nd_5126000 = nd_5125360.new_child(label="Node5125360", taxon=tax_4839056, edge_length=0.070762, oid="Node5126000")
    nd_5126000.edge.oid = "Edge5126128"
    nd_5125392 = nd_5125136.new_child(label="Node5125136", taxon=None, edge_length=0.002182, oid="Node5125392")
    nd_5125392.edge.oid = "Edge5125552"
    nd_5125648 = nd_5125136.new_child(label="Node5125136", taxon=None, edge_length=0.011197, oid="Node5125648")
    nd_5125648.edge.oid = "Edge5125232"
    nd_5125456 = nd_5125392.new_child(label="Node5125392", taxon=None, edge_length=0.014459, oid="Node5125456")
    nd_5125456.edge.oid = "Edge5125616"
    nd_5125808 = nd_5125392.new_child(label="Node5125392", taxon=tax_4838160, edge_length=0.054698, oid="Node5125808")
    nd_5125808.edge.oid = "Edge5125776"
    nd_5125520 = nd_5125456.new_child(label="Node5125456", taxon=tax_4796368, edge_length=0.075877, oid="Node5125520")
    nd_5125520.edge.oid = "Edge5125680"
    nd_5125744 = nd_5125456.new_child(label="Node5125456", taxon=tax_4838288, edge_length=0.071375, oid="Node5125744")
    nd_5125744.edge.oid = "Edge5125712"
    nd_5125840 = nd_5125648.new_child(label="Node5125648", taxon=tax_4838672, edge_length=0.062291, oid="Node5125840")
    nd_5125840.edge.oid = "Edge5125904"
    nd_5126032 = nd_5125648.new_child(label="Node5125648", taxon=None, edge_length=0.00521, oid="Node5126032")
    nd_5126032.edge.oid = "Edge5125968"
    nd_5125872 = nd_5126032.new_child(label="Node5126032", taxon=tax_4838832, edge_length=0.074463, oid="Node5125872")
    nd_5125872.edge.oid = "Edge5125936"
    nd_5126160 = nd_5126032.new_child(label="Node5126032", taxon=tax_4838448, edge_length=0.103887, oid="Node5126160")
    nd_5126160.edge.oid = "Edge5126096"
    nd_5126352 = nd_5125584.new_child(label="Node5125584", taxon=tax_4838224, edge_length=0.071307, oid="Node5126352")
    nd_5126352.edge.oid = "Edge5126288"
    nd_5126448 = nd_5125584.new_child(label="Node5125584", taxon=tax_4838128, edge_length=0.073392, oid="Node5126448")
    nd_5126448.edge.oid = "Edge5126416"
    tree_5126480 = dendropy.Tree(label="Tree10", taxon_set=reference_tree_list.taxon_set, oid="Tree5126480")
    reference_tree_list.append(tree_5126480, reindex_taxa=False)
    nd_5126672 = tree_5126480.seed_node.new_child(label="Node5126608", taxon=tax_4795760, edge_length=0.053901, oid="Node5126672")
    nd_5126672.edge.oid = "Edge5126704"
    nd_5126736 = tree_5126480.seed_node.new_child(label="Node5126608", taxon=None, edge_length=0.013947, oid="Node5126736")
    nd_5126736.edge.oid = "Edge5124816"
    nd_5127632 = tree_5126480.seed_node.new_child(label="Node5126608", taxon=tax_4838960, edge_length=0.107868, oid="Node5127632")
    nd_5127632.edge.oid = "Edge5127856"
    nd_5126256 = nd_5126736.new_child(label="Node5126736", taxon=None, edge_length=0.002034, oid="Node5126256")
    nd_5126256.edge.oid = "Edge5126832"
    nd_5127824 = nd_5126736.new_child(label="Node5126736", taxon=tax_4838160, edge_length=0.064652, oid="Node5127824")
    nd_5127824.edge.oid = "Edge5127600"
    nd_5126768 = nd_5126256.new_child(label="Node5126256", taxon=None, edge_length=0.00777, oid="Node5126768")
    nd_5126768.edge.oid = "Edge5126896"
    nd_5127472 = nd_5126256.new_child(label="Node5126256", taxon=None, edge_length=0.00373, oid="Node5127472")
    nd_5127472.edge.oid = "Edge5127152"
    nd_5126384 = nd_5126768.new_child(label="Node5126768", taxon=None, edge_length=0.003041, oid="Node5126384")
    nd_5126384.edge.oid = "Edge5126960"
    nd_5069904 = nd_5126768.new_child(label="Node5126768", taxon=None, edge_length=0.008083, oid="Node5069904")
    nd_5069904.edge.oid = "Edge5126544"
    nd_5126864 = nd_5126384.new_child(label="Node5126384", taxon=tax_4837744, edge_length=0.07927, oid="Node5126864")
    nd_5126864.edge.oid = "Edge5127024"
    nd_5127120 = nd_5126384.new_child(label="Node5126384", taxon=None, edge_length=0.016769, oid="Node5127120")
    nd_5127120.edge.oid = "Edge5127056"
    nd_5126992 = nd_5127120.new_child(label="Node5127120", taxon=tax_4838672, edge_length=0.053956, oid="Node5126992")
    nd_5126992.edge.oid = "Edge5126928"
    nd_5127248 = nd_5127120.new_child(label="Node5127120", taxon=tax_4838448, edge_length=0.098874, oid="Node5127248")
    nd_5127248.edge.oid = "Edge5127184"
    nd_5126320 = nd_5069904.new_child(label="Node5069904", taxon=None, edge_length=0.010712, oid="Node5126320")
    nd_5126320.edge.oid = "Edge5127088"
    nd_5127344 = nd_5069904.new_child(label="Node5069904", taxon=tax_4839056, edge_length=0.075934, oid="Node5127344")
    nd_5127344.edge.oid = "Edge5127312"
    nd_5069872 = nd_5126320.new_child(label="Node5126320", taxon=tax_4796368, edge_length=0.075223, oid="Node5069872")
    nd_5069872.edge.oid = "Edge5124432"
    nd_5126800 = nd_5126320.new_child(label="Node5126320", taxon=tax_4838224, edge_length=0.065805, oid="Node5126800")
    nd_5126800.edge.oid = "Edge5127216"
    nd_5127408 = nd_5127472.new_child(label="Node5127472", taxon=tax_4838832, edge_length=0.082425, oid="Node5127408")
    nd_5127408.edge.oid = "Edge5124720"
    nd_5127568 = nd_5127472.new_child(label="Node5127472", taxon=None, edge_length=0.008508, oid="Node5127568")
    nd_5127568.edge.oid = "Edge5127504"
    nd_5127280 = nd_5127568.new_child(label="Node5127568", taxon=None, edge_length=0.010583, oid="Node5127280")
    nd_5127280.edge.oid = "Edge5127376"
    nd_5127792 = nd_5127568.new_child(label="Node5127568", taxon=tax_4838512, edge_length=0.070289, oid="Node5127792")
    nd_5127792.edge.oid = "Edge5127728"
    nd_5127440 = nd_5127280.new_child(label="Node5127280", taxon=tax_4838288, edge_length=0.075429, oid="Node5127440")
    nd_5127440.edge.oid = "Edge5127664"
    nd_5127760 = nd_5127280.new_child(label="Node5127280", taxon=tax_4838128, edge_length=0.0787, oid="Node5127760")
    nd_5127760.edge.oid = "Edge5127696"

    # set labels of nodes with taxa to taxon label, else oid (for consistent
    # identification in debugging)
    for t in reference_tree_list:
        t.assign_node_labels_from_taxon_or_oid()

    return reference_tree_list





class RepeatedRandom(Random):
    """
    An overload of Random that returns numbers from a circular list of 1000
    numbers. This is useful for testing.
    """
    def __init__(self, x=None):
        """
        Initialize an instance. Optional argument x controls seeding,
        as for Random.seed().
        """
        self.period = 1000
        self.index = 0
        self.rand_nums = [
0.31443191705262208, 0.2271212707069612, 0.49820681873991857,
0.44722415062740772, 0.92936218264281234, 0.59873171836809469,
0.19982411078265749, 0.98903094700215477, 0.026814421675923739,
0.24553294072971865, 0.90473066455214113, 0.63500099590551295,
0.062485878633049552, 0.40324970479428146, 0.75995556949919885,
0.98363548217857744, 0.54357505325205402, 0.85219514190867551,
0.61210275315338891, 0.74896253672556812, 0.70873501751854029,
0.87114422189514784, 0.92753457916157112, 0.60789956222627695,
0.34263919254268715, 0.33226005584940799, 0.13578382005587253,
0.23424699675198624, 0.50146300494834883, 0.77817054412039099,
0.15334683515195979, 0.38000295395822037, 0.80658234859686051,
0.52900928702132466, 0.16707757903062137, 0.9906812734146403,
0.068901748752614456, 0.49099641208134737, 0.47428802125713176,
0.043623746806798258, 0.27436313661767664, 0.6062541793704701,
0.46702497674258392, 0.010313222742326711, 0.045727337925076994,
0.096940257042402833, 0.8451024639104916, 0.10813661350545056,
0.36447286770120435, 0.17255941523459573, 0.4458754882710313,
0.56730670677709205, 0.66242313381279672, 0.2411633825018501,
0.64646607162710512, 0.51004102559263365, 0.14276617568106564,
0.49980274479883058, 0.56829545095899059, 0.14666139536796241,
0.27122315425029742, 0.062076706517107283, 0.29784232388014031,
0.46937153843131396, 0.39662305708213441, 0.77502963788221912,
0.34431246165606266, 0.12406976774146794, 0.77954233508627691,
0.98442357618953502, 0.65870863138879365, 0.76317430726323587,
0.32193320354888499, 0.41567072312058628, 0.49826478349801484,
0.94425296232063916, 0.40038568591413637, 0.83149782296288244,
0.50820834756269906, 0.45158625835932231, 0.87137420135576282,
0.18969443418561149, 0.70471518687224011, 0.43044950270947491,
0.14105841021622767, 0.46819020426556379, 0.28429748379509456,
0.33490273723410269, 0.92074700534564446, 0.85633367206412803,
0.44757299380703963, 0.45672655371946225, 0.70731920284944627,
0.75300537008067314, 0.2232676076061364, 0.61512018215155906,
0.62469408147020833, 0.70676695570246395, 0.90236199735102374,
0.079786365777333557, 0.13031703283041907, 0.69722911855322822,
0.065641265674688865, 0.14489805312547432, 0.50224362737950645,
0.44567510180376657, 0.89500463323088286, 0.90910794360813263,
0.100405084345489, 0.55357619376108902, 0.2631932106857483,
0.5368157483119983, 0.99717971859556964, 0.90080123405742396,
0.8919749597060378, 0.96319408512128946, 0.72104763380068249,
0.33794570782212341, 0.92887336720989366, 0.57094598529869378,
0.81251580554905201, 0.49833566254350214, 0.074578046896970562,
0.15742891672534376, 0.56030040624723298, 0.62341890462128202,
0.41169999530912615, 0.80360387569467229, 0.24478587276981478,
0.53034050960413304, 0.14555402088453773, 0.57197267277344321,
0.7492744117968001, 0.40784211694021266, 0.11232415851698097,
0.29253785687755374, 0.20417762476978685, 0.45634845668850621,
0.97792963901394625, 0.7152084884227089, 0.27861419201158755,
0.62803675215958377, 0.85168636397984088, 0.6597634578414231,
0.078721939023294496, 0.1349959738362746, 0.27468049112529636,
0.02946495780617886, 0.42358091820270471, 0.77729283466797494,
0.02988444242940469, 0.83123858689441787, 0.76934867134442064,
0.64741431523640336, 0.3912084129095641, 0.6045302455094459,
0.34111190327675178, 0.55496514485688919, 0.60933955006653029,
0.54722840485533641, 0.10093773352434698, 0.58381620951367252,
0.076014095573963547, 0.4796220011714295, 0.10066813035566546,
0.9886687525572585, 0.73671841735259536, 0.77680789701691333,
0.50699813617100875, 0.16008125326001144, 0.87823152126474469,
0.14575661058543921, 0.60250564670920836, 0.61902906767558685,
0.049351327842071746, 0.49953084340237608, 0.41093474250555473,
0.97552953639011253, 0.51676060749834307, 0.93071762821570492,
0.24874925056574937, 0.060860377679168187, 0.39232281732286711,
0.80356313595940954, 0.786440293404014, 0.64625142562190385,
0.94871068219852805, 0.84212491538939049, 0.36449282817053852,
0.82234210696581844, 0.61187529740994917, 0.94733352117304892,
0.92772114792736216, 0.26248960149311207, 0.41887757379152035,
0.12936694229983947, 0.57001710049713061, 0.010005635615027653,
0.68070258603920941, 0.55882732694910564, 0.87295001985704523,
0.36658461453172919, 0.23237288683316459, 0.64910872189673208,
0.084974460323642309, 0.34096916028967794, 0.47612236632740101,
0.24552896768330368, 0.042646747131294127, 0.020766228915284457,
0.35707004990778979, 0.40859001679046225, 0.48306480955720965,
0.64214193409924802, 0.81819075869930569, 0.69944147181886729,
0.31408848949179724, 0.50536430567600688, 0.40477779512142853,
0.59493090819332872, 0.0089839382967883408, 0.67524592540041495,
0.63002984647419125, 0.177094541575886, 0.1749310289571564,
0.1932982649036219, 0.21408296727424447, 0.20929754851582516,
0.46529022624963767, 0.78707697700412727, 0.41333186638156294,
0.11574974641989877, 0.53867238158603425, 0.82175328149463944,
0.67639453631223012, 0.25014369054860197, 0.49854609920360005,
0.24030238095944334, 0.48457535655726824, 0.25746386922489906,
0.46106628357880619, 0.81808647586913241, 0.15491536835125597,
0.79565090323486443, 0.72113371275143234, 0.15892353894428957,
0.83542074414500689, 0.45886961733146725, 0.43162455397970179,
0.50755071650515615, 0.37972028616229669, 0.69097183736373469,
0.3392880552007751, 0.56391954383553367, 0.51596972483937542,
0.3350085000360179, 0.76262188690414423, 0.04059566986953278,
0.72096512147672664, 0.73804031506128298, 0.68867017527258689,
0.4839824851362271, 0.56147196547531775, 0.8322849786127503,
0.41930393927188214, 0.95129787073939975, 0.077389437920532544,
0.73312136710456399, 0.85968301610141051, 0.84910535174694934,
0.67672577614868568, 0.39840657520227629, 0.94849298075019639,
0.73238765873246292, 0.65963864561479935, 0.83040338914322942,
0.84313679569585709, 0.74823189764524778, 0.26818361235787402,
0.64491449933958134, 0.29926105221964849, 0.40622193728965883,
0.43894970523034738, 0.37205376178932681, 0.68288316787791026,
0.51152834706969119, 0.040094666318732486, 0.51835551855297335,
0.69914133860734295, 0.64041327950776628, 0.53911859102535886,
0.47419203100048946, 0.98833247945952285, 0.45942391149428241,
0.48638367210938549, 0.29548147882844189, 0.5516776265496679,
0.80321163510972327, 0.86442523514958547, 0.47263593924626335,
0.30319463616137066, 0.36245407246198802, 0.75558505880342719,
0.95739479621470969, 0.82587358532163613, 0.11247960960382175,
0.46419576976880894, 0.90535923412978114, 0.18505075271874682,
0.065084859736668332, 0.41088413144481051, 0.1163853084536014,
0.71752303978381093, 0.9696211810020936, 0.82196703884289235,
0.72049774195929861, 0.7752354427755731, 0.09388893717118818,
0.1325673610022654, 0.25435966713715885, 0.61212919784239284,
0.99337146179390368, 0.93141574366900326, 0.82812399481728349,
0.4428708918637918, 0.98219708766421743, 0.59186175672774433,
0.21553895655929745, 0.18919035141314622, 0.37778541717248304,
0.70150248593439224, 0.038493068549859011, 0.45017106304425603,
0.012349355068150825, 0.057210716141382956, 0.36131218967977441,
0.30368749532763406, 0.76748032875773708, 0.074477506286847683,
0.36360791444330609, 0.9623640161286835, 0.80248013595119039,
0.64840527058023545, 0.27315029782770062, 0.26570186631332293,
0.032536083555708806, 0.48705081387400317, 0.59687133594383202,
0.62350625415636296, 0.28500931577707378, 0.75038812771746921,
0.63642564096896859, 0.50344381173183284, 0.046246962255789281,
0.65096776559704583, 0.99335669174451846, 0.60052902908479966,
0.91310341608001355, 0.26914409939868555, 0.98510339278263348,
0.46710139889858215, 0.39325219794139132, 0.52623165641816971,
0.96037264125943966, 0.79579351771376217, 0.25476456445703821,
0.10940016989945811, 0.99020083999370367, 0.22552307457113296,
0.71042251786619393, 0.039664945075574831, 0.96216550041626125,
0.10016594114658017, 0.96625136816262613, 0.040562665307392276,
0.75347472744879929, 0.56515679140674746, 0.27644225575309123,
0.89878952868547379, 0.83261330046450299, 0.15034835204416519,
0.11408844903704218, 0.9576585991414176, 0.5366198880109444,
0.35420905430064598, 0.68357000407099699, 0.94642589611259054,
0.56475741360759735, 0.77754449030481632, 0.88616730803987898,
0.65641106985679387, 0.67437950243977707, 0.47221622506846839,
0.50062117520556204, 0.91970430648981094, 0.18131053303862033,
0.21272268806612415, 0.05377388641792813, 0.31007890951373962,
0.32618319789265249, 0.34290429703351155, 0.21549459205241117,
0.2524280170570633, 0.81709409505752828, 0.75815705258870691,
0.47234205496003512, 0.78357212776826923, 0.53932071107791935,
0.85077682027667501, 0.83311361991293376, 0.90786178925224947,
0.15239185151498524, 0.55034674406410433, 0.75795968586040596,
0.65221986813104249, 0.65040096790719815, 0.83824965080922464,
0.43799507441117524, 0.87074085705683202, 0.94866130510988389,
0.52427380088194253, 0.75466867079311706, 0.59072719583099642,
0.75196616934637517, 0.86124226022283668, 0.4083157725740546,
0.18255991491456602, 0.8094909841226301, 0.9238223349425172,
0.3489557609572872, 0.63121018395963568, 0.17545841707657228,
0.49821878413082843, 0.80327811868126764, 0.28922419662100818,
0.096050290301955776, 0.538241792459762, 0.1713976104609598,
0.42123504272249834, 0.070501770640446049, 0.10049193448270621,
0.40017124813453819, 0.49753518158427923, 0.9886569977137456,
0.0036223759804981936, 0.72003196594546082, 0.49450652913802473,
0.85446043489120027, 0.27576601630834741, 0.69832406474399744,
0.84560364778003083, 0.51772138467128959, 0.39288177031969562,
0.96205954073664657, 0.01036595300701304, 0.20813795619564024,
0.58211793695016312, 0.093444683710641629, 0.56200191781191788,
0.23670382060870676, 0.11486440233716622, 0.95144899622178924,
0.16548179609792213, 0.749604411358611, 0.8489642046684549,
0.74481664127095337, 0.14323986526663668, 0.57678405290507173,
0.70737322386084223, 0.67618648905068834, 0.2465839168756111,
0.78710019275683885, 0.28787543298114859, 0.097255946101287516,
0.52633018343137772, 0.13397988165453845, 0.34982242740349179,
0.30083568024580221, 0.23002499396511533, 0.39455145880055487,
0.53909282369214884, 0.5098939390012095, 0.99519630939119719,
0.65245304515657243, 0.6278994145755239, 0.28974657946211091,
0.60721277396019713, 0.89345895195698755, 0.29467240439427766,
0.85481899776248127, 0.52193388182858658, 0.9133801554519243,
0.50427500446909879, 0.85002896645153558, 0.87105707616743488,
0.88274737164353945, 0.39826295447590065, 0.97178570093923,
0.59557258786570377, 0.85176449380549812, 0.63033913649146966,
0.85138887004600872, 0.84230334700094767, 0.54577129632427634,
0.066205916782848928, 0.055465695420991001, 0.34186500465807157,
0.83648536922433414, 0.50566045826603923, 0.075635038536681964,
0.83140814268128949, 0.93468963034722885, 0.59507556482872892,
0.99365469502662762, 0.058149318098988489, 0.54942554831926116,
0.1649351536556477, 0.36156873463203754, 0.48927371055917002,
0.46630723876444657, 0.43494631989319743, 0.068712919541944251,
0.017189475299067558, 0.0088556842001139557, 0.36141417615913152,
0.85308467014593914, 0.62582981485046285, 0.9581126836643491,
0.61919083481355885, 0.27660431371386407, 0.69629612158347465,
0.1882410189585707, 0.88515805020035998, 0.22963170832643232,
0.26232029298479531, 0.37470837626189657, 0.23675339810004337,
0.69090831156040311, 0.80140513633505306, 0.92279857466770476,
0.6065036854306598, 0.018178829690600251, 0.093450588123038192,
0.76349911808481519, 0.31488524475447477, 0.22338899889386088,
0.34791728969913149, 0.28557610843267833, 0.35156310782500766,
0.92534773820401084, 0.82482329374533125, 0.21731554294414435,
0.36144672624854923, 0.65796790552131279, 0.19540075868896012,
0.16884193855701868, 0.36038857167145644, 0.59040341628371806,
0.88059940221166688, 0.96497956855411582, 0.081138766376561589,
0.21706711129437584, 0.090823322263668516, 0.2088945616604938,
0.073581231114019929, 0.46847688801420184, 0.075993449653287581,
0.61216656870467434, 0.76075674317013886, 0.61140551802038356,
0.16549265935398116, 0.96911925303347668, 0.60625163307811458,
0.72399365784065461, 0.002644877214834751, 0.4420448631751952,
0.27100053245915956, 0.01049476392458959, 0.34914797964595035,
0.62620571336031605, 0.99174582921261834, 0.10282266856087319,
0.42527709634381805, 0.61827337279979699, 0.042628563903215788,
0.99425805906104536, 0.18835802817121261, 0.51663287801309388,
0.00071687243225870834, 0.19224356026637945, 0.35424929691079809,
0.27771613021877772, 0.1702926412836574, 0.57744611720777628,
0.09411402810943903, 0.55800959037312126, 0.6755403793946978,
0.63971514472872804, 0.84208699726913872, 0.73924273979392963,
0.31797310932926437, 0.43560559324304082, 0.72822006078137302,
0.34021010712720468, 0.65800890008991064, 0.92506411624434082,
0.87027885934255056, 0.65165487411908385, 0.57590949220775711,
0.5305811800645196, 0.21300563938823236, 0.85319355429992927,
0.41813336013349511, 0.6402346817441521, 0.12696860712654501,
0.47007546629474051, 0.48542408063649056, 0.075530125640996482,
0.51893066398146626, 0.81844064972408193, 0.10747498923675491,
0.1096282809806246, 0.25375564153884489, 0.36633088159110827,
0.15801256597931457, 0.922151809927651, 0.031768569931787893,
0.43576325573288455, 0.21113431542672856, 0.31596674335427177,
0.85599318804265878, 0.81208478273596529, 0.10816056162310417,
0.90119858816078713, 0.44899529012597006, 0.33724510161922039,
0.88408157989409231, 0.087297142667446925, 0.45339378304422251,
0.92291152162924073, 0.37146410762513915, 0.28635300349987958,
0.52531368082180052, 0.19987533551229164, 0.28515195833401197,
0.78696780334000849, 0.18409735841751462, 0.63115068875856151,
0.014052055220890813, 0.75537970046662439, 0.69667760101543752,
0.53849726798980924, 0.69966987907192613, 0.68409265434583921,
0.96233996652121068, 0.33707239023588242, 0.78097691869862418,
0.77797849881511727, 0.80270691051387077, 0.48887933213516466,
0.66342940395001104, 0.75499746294870507, 0.19781780665223792,
0.18272687155761635, 0.96623704587636516, 0.5238656866056709,
0.44386011541475057, 0.64522237359183865, 0.50810450980814414,
0.76823725686412425, 0.47406197139443873, 0.41873386804615276,
0.47922922274530522, 0.31067765642017786, 0.59285344557631647,
0.40805366505854601, 0.81430345633987966, 0.71662142693747621,
0.85685183873929738, 0.14406177373290485, 0.51104814471267757,
0.14071252290705238, 0.06969316390364011, 0.2583719584556573,
0.68451057254140613, 0.0016165163630846857, 0.080538867172593731,
0.82231364544818153, 0.16851979870461298, 0.25747974558536035,
0.45300989528105839, 0.98526738767557143, 0.97451169000984994,
0.038117904404271652, 0.01765330749662497, 0.6744420693012253,
0.71078268782386889, 0.55552919512325083, 0.19525732457142908,
0.98498115572050793, 0.94424033466702706, 0.30238482785473619,
0.51735289069888513, 0.0218844360276651, 0.37095525152876008,
0.98856409873656281, 0.81171161518835877, 0.62205790570841402,
0.13903995589268614, 0.24317363346222531, 0.28691664626948732,
0.41529758456162613, 0.58102138159149475, 0.29511103490158153,
0.81032771243549684, 0.54664215124727655, 0.86806884248043081,
0.40148973967126578, 0.72785143234314909, 0.83093687762801249,
0.68138259848769756, 0.76590081936459165, 0.50823725006607534,
0.67456323322438305, 0.39077640138345837, 0.42605548372487245,
0.27751254873155762, 0.19227397148004932, 0.43513135333712738,
0.012639800127961398, 0.24034640152670483, 0.10027447128694145,
0.89318562867730877, 0.50849414601968046, 0.20657343439508324,
0.57147540212723991, 0.47361978176002362, 0.11753188924212987,
0.98670021096046046, 0.8207811703836605, 0.66086365421202642,
0.60966298584013634, 0.94363899598878753, 0.70628481771581986,
0.17426291736596278, 0.024325687803500751, 0.40140066436409716,
0.89813827762266019, 0.23245614268034809, 0.41721013902193649,
0.74346705142425296, 0.69053314604711236, 0.55492823893072951,
0.87243520608233738, 0.8312181578062765, 0.97684515195591803,
0.35216590061944664, 0.69012293976323458, 0.66421283923840491,
0.88350728730414396, 0.76583235834404084, 0.88512324584587998,
0.28611466957123011, 0.64601337868076381, 0.14372944791838049,
0.78288264166083665, 0.2487441999079687, 0.4718699501149034,
0.52975394857724545, 0.33705125203762321, 0.090787594259293392,
0.31240428763858863, 0.90097506501517788, 0.59462802131684955,
0.49988656761918837, 0.96370898732399146, 0.56268934292939077,
0.36938414960144983, 0.7883504205258377, 0.25721869603089698,
0.28997277985430103, 0.84515931522936061, 0.81404163063999102,
0.10469687422346452, 0.4395925152686998, 0.87762561700081865,
0.68571872122065192, 0.15240218006329653, 0.080152884004887737,
0.29408742629437801, 0.43230665934130696, 0.95930050828503155,
0.88164080303001868, 0.014228459746536415, 0.54117421597118331,
0.3093599606166112, 0.2253169883201902, 0.91886214568858338,
0.1911888563036076, 0.95893282911333255, 0.42194833819579569,
0.64958524906175685, 0.27438310864894144, 0.47728534096729958,
0.72632268640363284, 0.44452695053483104, 0.43552220277519649,
0.40980432172030312, 0.77230248591055173, 0.9077817668261553,
0.61445585506231115, 0.82788913948136988, 0.36586382522685101,
0.13940582072259522, 0.66770701082795725, 0.62362384669896942,
0.15927430003372178, 0.90459467071989419, 0.13333361766042984,
0.11466869457171014, 0.28276573967969287, 0.21175722591814683,
0.73761579576599035, 0.46162293964555479, 0.53931182885661688,
0.19346800949975318, 0.59273833997939807, 0.85593387242226937,
0.20636310091981991, 0.74740172028853169, 0.0097987635504059867,
0.53220187213375414, 0.10477473988968611, 0.42588178132469767,
0.49116443200629045, 0.99047240022140914, 0.68964358689657168,
0.036585644664419381, 0.38343278241796397, 0.78422102834621799,
0.53000667898091136, 0.56832676501516, 0.1618813939818381,
0.78120478334302312, 0.81318259290902628, 0.40539491724822418,
0.55822463931191346, 0.48382947242366969, 0.013999510513232449,
0.30319124448812262, 0.36263529307915177, 0.97471886280882758,
0.89616494921964651, 0.84141953467994779, 0.27375941926765657,
0.76617869775906056, 0.96900963369293314, 0.43640341590511755,
0.037322202213532885, 0.64673899008553493, 0.35767838217965664,
0.7889272154535828, 0.20109932743575909, 0.68153577205564464,
0.91221429179121627, 0.086009200734357472, 0.59721147029656663,
0.27523447803309387, 0.73045391968310958, 0.024866069256412437,
0.86107232375125931, 0.91463723977350431, 0.87186138748263264,
0.39324191303855216, 0.096161379593851071, 0.32068370878413976,
0.78897571163642655, 0.78474613759719869, 0.34331601538487666,
0.40828309056141754, 0.15694516839134898, 0.83789401420436682,
0.80883177334456757, 0.21466083788682122, 0.24325789921874352,
0.83206267445430881, 0.38798636193955949, 0.95060347757178409,
0.24085962895484037, 0.27990366953894874, 0.56592755582224441,
0.18271630447568943, 0.33717455928632301, 0.9592910233932247,
0.69853867366898292, 0.45801598729692572, 0.87473822819048719,
0.20057204411755836, 0.41238641408676235, 0.48833064608189392,
0.84049561182553301, 0.42519968407784081, 0.43566093025099217,
0.041244638889872287, 0.97655130365177889, 0.28173859355824082,
0.98771197197012328, 0.81858850686964801, 0.26580819922534471,
0.45944722778299152, 0.3291016001586401, 0.31426388220315016,
0.13855305849387545, 0.74478342127839636, 0.43227250519402238,
0.94387295682438122, 0.45541032020373573, 0.77009137570947128,
0.31092801905533718, 0.66941813407075434, 0.4268626820379563,
0.20020363748681624, 0.92365265242403805, 0.23553797096803619,
0.4799416350221184, 0.94229554199687893, 0.8205115632610912,
0.15081497950722322, 0.056775689990967493, 0.52364515059006689,
0.61167866114344005, 0.49131714402095361, 0.27172462716061663,
0.34401473726073106, 0.1309776400644278, 0.1155305645196063,
0.64795318728951179, 0.05890399840832139, 0.3826763630244947,
0.72087369958143899, 0.27041817171916704, 0.6326836251042649,
0.75987299662928076, 0.52755368073763942, 0.87262870339772391,
0.6556925798696458, 0.67540950390279975, 0.14720672733677387,
0.81844445101944563, 0.59474364657015155, 0.16749016121926463,
0.62166156198457123, 0.4558918955138721, 0.17287013915801952,
0.56204353906568527, 0.94645178075956382, 0.83732842723959022,
0.98291695415637659, 0.29899190475031368, 0.98472724100866804,
0.17415887677241315, 0.36897659959880857, 0.48239959673737964,
0.62276551493865662, 0.030010365540096617, 0.7689900566409642,
0.032718717678784426, 0.062869321963432312, 0.11209046647623322,
0.84948915240213352, 0.83212575464788829, 0.061793679140391355,
0.22412336143745482, 0.55303755867428062, 0.79750889129943847,
0.80870433937874264, 0.62110916094845803, 0.23896504325793966,
0.20101301051309206, 0.68067301304891681, 0.30134264414403189,
0.91596006365773341, 0.89156690577779907, 0.6733373469244831,
0.4752275437610407, 0.37783271284483444, 0.45649144239182682,
0.33465933618787669, 0.18984314534253122, 0.31178645703848173,
0.35169927993508443, 0.34634536877091138, 0.56400514652008393,
0.1061010450319515, 0.78804117314026834, 0.29742453358636112,
0.99729169101350168, 0.14889988220740147, 0.57392391117992947,
0.18823074515724048]

        self.seed(x)
        self.gauss_next = None

    def seed(self, a=None):
        """Initialize internal state."""
        if a is None:
            a = 1
        self.index = a % self.period
        super(Random, self).seed(a)
        self.gauss_next = None

    def getstate(self):
        """Return internal state; can be passed to setstate() later."""
        gs = list(Random.getstate(self)[1])
        gs.append(self.index)
        return -1, tuple(gs), self.gauss_next

    def setstate(self, state):
        """Restore internal state from object returned by getstate()."""
        version = state[0]
        if version == -1:
            version, internalstate, self.gauss_next = state
            l_is = list(internalstate)
            self.index = l_is.pop(-1)
            t = (Random.VERSION, tuple(l_is), self.gauss_next)
            Random.setstate(self, t)
        else:
            raise ValueError("state with version %s passed to "
                             "Random.setstate() of version %s" %
                             (version, -1))
    def jumpahead(self, n):
        self.index = (self.index + n) % self.period
    def random(self):
        r = self.rand_nums[self.index]
        self.index += 1
        if self.index >= self.period:
            self.index = 0
        return r
