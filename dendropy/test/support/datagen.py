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
from dendropy.dataobject.tree import NodeRelationship
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
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Morelia boeleni','Node5594256','Node5593936','Morelia viridis','Node5593872','Antaresia maculosa','Node5593808','Antaresia stimsoni','Node5593744','Python timoriensis','Node5593328','Morelia oenpelliensis','Node5593616','Morelia bredli','Node5593648','Antaresia perthensis','Python brongersmai','Node5594896','Node5593584','Morelia carinata','Node5593456'],
        ['Aspidites ramsayi','Bothrochilus boa','Morelia boeleni','Morelia oenpelliensis','Node7353200','Node5593168','Morelia viridis','Node7352688','Antaresia perthensis','Node7352880','Liasis fuscus','Antaresia stimsoni','Python timoriensis','Morelia carinata','Node7353744','Python brongersmai','Node7353712','Antaresia maculosa','Node7353840','Node7353680','Node7353360','Morelia bredli','Node7353232','Node7352752'],
        ['Aspidites ramsayi','Bothrochilus boa','Antaresia perthensis','Node7354544','Python brongersmai','Node7354096','Python timoriensis','Antaresia maculosa','Node7354928','Morelia oenpelliensis','Node7354736','Node7354448','Liasis fuscus','Morelia viridis','Morelia bredli','Node7355408','Node7355344','Node7354256','Antaresia stimsoni','Node7354416','Morelia carinata','Morelia boeleni','Node7355600','Node7354288'],
        ['Aspidites ramsayi','Bothrochilus boa','Antaresia perthensis','Node7356112','Morelia bredli','Node7355920','Liasis fuscus','Morelia oenpelliensis','Node7409872','Morelia viridis','Node7409840','Python timoriensis','Morelia carinata','Node7410320','Antaresia maculosa','Node7410064','Node7356208','Python brongersmai','Node7409776','Antaresia stimsoni','Node7356144','Node7356080','Morelia boeleni','Node7355952'],
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Morelia carinata','Node7411216','Node7411024','Antaresia stimsoni','Antaresia perthensis','Node7411344','Antaresia maculosa','Node5593424','Morelia viridis','Morelia bredli','Morelia oenpelliensis','Node7411824','Node7411600','Python timoriensis','Node7411504','Node7411248','Python brongersmai','Morelia boeleni','Node7412144','Node7411440','Node7410896'],
        ['Aspidites ramsayi','Bothrochilus boa','Python brongersmai','Node7412528','Python timoriensis','Node7412016','Morelia bredli','Antaresia maculosa','Node7412688','Node7412496','Liasis fuscus','Morelia oenpelliensis','Node7412912','Antaresia stimsoni','Morelia carinata','Morelia boeleni','Node7413648','Node7413456','Antaresia perthensis','Node7413424','Node7412624','Morelia viridis','Node7413168','Node7412368'],
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Node7504304','Antaresia stimsoni','Morelia viridis','Python timoriensis','Morelia boeleni','Node7505104','Node7505008','Node7504688','Antaresia maculosa','Morelia oenpelliensis','Node7505136','Node7504624','Morelia bredli','Node7504368','Morelia carinata','Node7504592','Python brongersmai','Node7504528','Antaresia perthensis','Node5593232','Node7504176'],
        ['Aspidites ramsayi','Bothrochilus boa','Python timoriensis','Node7506160','Antaresia perthensis','Node7506096','Antaresia maculosa','Morelia oenpelliensis','Node7506352','Node7506032','Morelia carinata','Node7505968','Antaresia stimsoni','Node7505872','Morelia boeleni','Node7412272','Morelia viridis','Python brongersmai','Node7506768','Node7505584','Liasis fuscus','Node7505840','Morelia bredli','Node7505680'],
        ['Aspidites ramsayi','Bothrochilus boa','Python brongersmai','Node7507504','Liasis fuscus','Antaresia maculosa','Node7507888','Antaresia stimsoni','Node7507824','Morelia bredli','Antaresia perthensis','Python timoriensis','Node7561744','Node7561296','Node7507568','Morelia carinata','Node7507792','Morelia oenpelliensis','Node7507536','Morelia viridis','Morelia boeleni','Node7561520','Node7507728','Node7507376'],
        ['Aspidites ramsayi','Bothrochilus boa','Morelia bredli','Python timoriensis','Node7562832','Node7562480','Liasis fuscus','Morelia viridis','Node7562864','Morelia carinata','Node7563056','Node7355728','Antaresia perthensis','Antaresia maculosa','Morelia boeleni','Node7563024','Morelia oenpelliensis','Node7563280','Node7563184','Node7562032','Antaresia stimsoni','Node7562448','Python brongersmai','Node7562320'],
    ]

def reference_tree_list_newick_string():
    return """\
        ('Aspidites ramsayi':0.056823,(((((((('Bothrochilus boa':0.073578,('Liasis fuscus':0.074347,'Morelia boeleni':0.07628)Node5594256:0.003063)Node5593936:0.011155,'Morelia viridis':0.066566)Node5593872:0.010103,'Antaresia maculosa':0.080749)Node5593808:0.007681,'Antaresia stimsoni':0.061996)Node5593744:0.0,'Python timoriensis':0.113605)Node5593328:0.002264,'Morelia oenpelliensis':0.074047)Node5593616:0.007655,'Morelia bredli':0.065046)Node5593648:0.003495,('Antaresia perthensis':0.065004,'Python brongersmai':0.107706)Node5594896:0.018613)Node5593584:0.011159,'Morelia carinata':0.079321)Node5593456;
        ('Aspidites ramsayi':0.065297,((('Bothrochilus boa':0.073565,('Morelia boeleni':0.069885,'Morelia oenpelliensis':0.06247)Node7353200:0.007879)Node5593168:0.016022,'Morelia viridis':0.066015)Node7352688:0.010946,'Antaresia perthensis':0.079424)Node7352880:0.001444,(('Liasis fuscus':0.082049,('Antaresia stimsoni':0.056894,((('Python timoriensis':0.092077,'Morelia carinata':0.0579)Node7353744:0.015909,'Python brongersmai':0.114385)Node7353712:0.010752,'Antaresia maculosa':0.079792)Node7353840:0.004552)Node7353680:0.006907)Node7353360:0.008007,'Morelia bredli':0.063725)Node7353232:0.004077)Node7352752;
        ('Aspidites ramsayi':0.06111,((((('Bothrochilus boa':0.067908,'Antaresia perthensis':0.077196)Node7354544:0.008688,'Python brongersmai':0.1069)Node7354096:0.010367,(('Python timoriensis':0.109318,'Antaresia maculosa':0.081382)Node7354928:0.004129,'Morelia oenpelliensis':0.068424)Node7354736:0.007731)Node7354448:0.005167,('Liasis fuscus':0.079478,('Morelia viridis':0.062239,'Morelia bredli':0.060484)Node7355408:0.013155)Node7355344:0.00574)Node7354256:0.006988,'Antaresia stimsoni':0.063808)Node7354416:0.005408,('Morelia carinata':0.069501,'Morelia boeleni':0.073602)Node7355600:0.012121)Node7354288;
        ('Aspidites ramsayi':0.053028,((('Bothrochilus boa':0.072831,'Antaresia perthensis':0.073004)Node7356112:0.012413,'Morelia bredli':0.069933)Node7355920:0.002567,((((('Liasis fuscus':0.072339,'Morelia oenpelliensis':0.062058)Node7409872:0.011723,'Morelia viridis':0.065297)Node7409840:0.011519,(('Python timoriensis':0.093762,'Morelia carinata':0.055992)Node7410320:0.022912,'Antaresia maculosa':0.079945)Node7410064:0.007647)Node7356208:5.134e-05,'Python brongersmai':0.118976)Node7409776:0.004803,'Antaresia stimsoni':0.064721)Node7356144:0.005969)Node7356080:0.014376,'Morelia boeleni':0.076894)Node7355952;
        ('Aspidites ramsayi':0.06391,('Bothrochilus boa':0.071835,('Liasis fuscus':0.075788,'Morelia carinata':0.07116)Node7411216:0.010888)Node7411024:0.00451,(((('Antaresia stimsoni':0.048754,'Antaresia perthensis':0.069682)Node7411344:0.012356,'Antaresia maculosa':0.077058)Node5593424:0.013817,(('Morelia viridis':0.07026,('Morelia bredli':0.066182,'Morelia oenpelliensis':0.068166)Node7411824:0.001443)Node7411600:0.003147,'Python timoriensis':0.108243)Node7411504:0.012897)Node7411248:0.011944,('Python brongersmai':0.108072,'Morelia boeleni':0.073251)Node7412144:0.00612)Node7411440:0.005404)Node7410896;
        ('Aspidites ramsayi':0.06501,((('Bothrochilus boa':0.059984,'Python brongersmai':0.100555)Node7412528:0.017824,'Python timoriensis':0.102639)Node7412016:0.011722,('Morelia bredli':0.0538,'Antaresia maculosa':0.07843)Node7412688:0.014836)Node7412496:0.003644,((('Liasis fuscus':0.068438,'Morelia oenpelliensis':0.066815)Node7412912:0.015326,(('Antaresia stimsoni':0.057007,('Morelia carinata':0.062766,'Morelia boeleni':0.079145)Node7413648:0.013392)Node7413456:0.010508,'Antaresia perthensis':0.083409)Node7413424:0.005897)Node7412624:0.0,'Morelia viridis':0.072091)Node7413168:0.008523)Node7412368;
        ('Aspidites ramsayi':0.062297,('Bothrochilus boa':0.069165,'Liasis fuscus':0.073495)Node7504304:0.009899,(((((('Antaresia stimsoni':0.061022,('Morelia viridis':0.071291,('Python timoriensis':0.098776,'Morelia boeleni':0.075069)Node7505104:0.017245)Node7505008:0.004698)Node7504688:0.003603,('Antaresia maculosa':0.075494,'Morelia oenpelliensis':0.063757)Node7505136:0.017571)Node7504624:0.000386,'Morelia bredli':0.070614)Node7504368:0.00713,'Morelia carinata':0.079534)Node7504592:0.010249,'Python brongersmai':0.114381)Node7504528:0.003672,'Antaresia perthensis':0.080896)Node5593232:0.003212)Node7504176;
        ('Aspidites ramsayi':0.05678,(((((((('Bothrochilus boa':0.065928,'Python timoriensis':0.10077)Node7506160:0.013813,'Antaresia perthensis':0.083335)Node7506096:0.004377,('Antaresia maculosa':0.073392,'Morelia oenpelliensis':0.065748)Node7506352:0.016433)Node7506032:0.005104,'Morelia carinata':0.07646)Node7505968:0.001866,'Antaresia stimsoni':0.06129)Node7505872:0.012314,'Morelia boeleni':0.081982)Node7412272:0.002833,('Morelia viridis':0.064242,'Python brongersmai':0.109434)Node7506768:0.010734)Node7505584:0.00182,'Liasis fuscus':0.079156)Node7505840:0.014479,'Morelia bredli':0.060374)Node7505680;
        ('Aspidites ramsayi':0.063004,('Bothrochilus boa':0.061413,'Python brongersmai':0.099507)Node7507504:0.018059,(((((('Liasis fuscus':0.075877,'Antaresia maculosa':0.071375)Node7507888:0.014459,'Antaresia stimsoni':0.054698)Node7507824:0.002182,('Morelia bredli':0.062291,('Antaresia perthensis':0.074463,'Python timoriensis':0.103887)Node7561744:0.00521)Node7561296:0.011197)Node7507568:0.010763,'Morelia carinata':0.070762)Node7507792:0.016533,'Morelia oenpelliensis':0.068818)Node7507536:0.003714,('Morelia viridis':0.071307,'Morelia boeleni':0.073392)Node7561520:0.010478)Node7507728:0.00238)Node7507376;
        ('Aspidites ramsayi':0.053901,(((('Bothrochilus boa':0.07927,('Morelia bredli':0.053956,'Python timoriensis':0.098874)Node7562832:0.016769)Node7562480:0.003041,(('Liasis fuscus':0.075223,'Morelia viridis':0.065805)Node7562864:0.010712,'Morelia carinata':0.075934)Node7563056:0.008083)Node7355728:0.00777,('Antaresia perthensis':0.082425,(('Antaresia maculosa':0.075429,'Morelia boeleni':0.0787)Node7563024:0.010583,'Morelia oenpelliensis':0.070289)Node7563280:0.008508)Node7563184:0.00373)Node7562032:0.002034,'Antaresia stimsoni':0.064652)Node7562448:0.013947,'Python brongersmai':0.107868)Node7562320;
    """

def reference_tree_list_node_relationships():
    treelist_node_references = [
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node5593456', child_labels=[], edge_length=0.056823, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node5593936', child_labels=[], edge_length=0.073578, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node5594256', child_labels=[], edge_length=0.074347, taxon_label='Liasis fuscus'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node5594256', child_labels=[], edge_length=0.07628, taxon_label='Morelia boeleni'),
            'Node5594256' : NodeRelationship(parent_label='Node5593936', child_labels=['Liasis fuscus','Morelia boeleni'], edge_length=0.003063, taxon_label=None),
            'Node5593936' : NodeRelationship(parent_label='Node5593872', child_labels=['Bothrochilus boa','Node5594256'], edge_length=0.011155, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node5593872', child_labels=[], edge_length=0.066566, taxon_label='Morelia viridis'),
            'Node5593872' : NodeRelationship(parent_label='Node5593808', child_labels=['Node5593936','Morelia viridis'], edge_length=0.010103, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5593808', child_labels=[], edge_length=0.080749, taxon_label='Antaresia maculosa'),
            'Node5593808' : NodeRelationship(parent_label='Node5593744', child_labels=['Node5593872','Antaresia maculosa'], edge_length=0.007681, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node5593744', child_labels=[], edge_length=0.061996, taxon_label='Antaresia stimsoni'),
            'Node5593744' : NodeRelationship(parent_label='Node5593328', child_labels=['Node5593808','Antaresia stimsoni'], edge_length=0.0, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node5593328', child_labels=[], edge_length=0.113605, taxon_label='Python timoriensis'),
            'Node5593328' : NodeRelationship(parent_label='Node5593616', child_labels=['Node5593744','Python timoriensis'], edge_length=0.002264, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node5593616', child_labels=[], edge_length=0.074047, taxon_label='Morelia oenpelliensis'),
            'Node5593616' : NodeRelationship(parent_label='Node5593648', child_labels=['Node5593328','Morelia oenpelliensis'], edge_length=0.007655, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node5593648', child_labels=[], edge_length=0.065046, taxon_label='Morelia bredli'),
            'Node5593648' : NodeRelationship(parent_label='Node5593584', child_labels=['Node5593616','Morelia bredli'], edge_length=0.003495, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node5594896', child_labels=[], edge_length=0.065004, taxon_label='Antaresia perthensis'),
            'Python brongersmai' : NodeRelationship(parent_label='Node5594896', child_labels=[], edge_length=0.107706, taxon_label='Python brongersmai'),
            'Node5594896' : NodeRelationship(parent_label='Node5593584', child_labels=['Antaresia perthensis','Python brongersmai'], edge_length=0.018613, taxon_label=None),
            'Node5593584' : NodeRelationship(parent_label='Node5593456', child_labels=['Node5593648','Node5594896'], edge_length=0.011159, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node5593456', child_labels=[], edge_length=0.079321, taxon_label='Morelia carinata'),
            'Node5593456' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node5593584','Morelia carinata'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7352752', child_labels=[], edge_length=0.065297, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node5593168', child_labels=[], edge_length=0.073565, taxon_label='Bothrochilus boa'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7353200', child_labels=[], edge_length=0.069885, taxon_label='Morelia boeleni'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7353200', child_labels=[], edge_length=0.06247, taxon_label='Morelia oenpelliensis'),
            'Node7353200' : NodeRelationship(parent_label='Node5593168', child_labels=['Morelia boeleni','Morelia oenpelliensis'], edge_length=0.007879, taxon_label=None),
            'Node5593168' : NodeRelationship(parent_label='Node7352688', child_labels=['Bothrochilus boa','Node7353200'], edge_length=0.016022, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7352688', child_labels=[], edge_length=0.066015, taxon_label='Morelia viridis'),
            'Node7352688' : NodeRelationship(parent_label='Node7352880', child_labels=['Node5593168','Morelia viridis'], edge_length=0.010946, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7352880', child_labels=[], edge_length=0.079424, taxon_label='Antaresia perthensis'),
            'Node7352880' : NodeRelationship(parent_label='Node7352752', child_labels=['Node7352688','Antaresia perthensis'], edge_length=0.001444, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7353360', child_labels=[], edge_length=0.082049, taxon_label='Liasis fuscus'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7353680', child_labels=[], edge_length=0.056894, taxon_label='Antaresia stimsoni'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7353744', child_labels=[], edge_length=0.092077, taxon_label='Python timoriensis'),
            'Morelia carinata' : NodeRelationship(parent_label='Node7353744', child_labels=[], edge_length=0.0579, taxon_label='Morelia carinata'),
            'Node7353744' : NodeRelationship(parent_label='Node7353712', child_labels=['Python timoriensis','Morelia carinata'], edge_length=0.015909, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7353712', child_labels=[], edge_length=0.114385, taxon_label='Python brongersmai'),
            'Node7353712' : NodeRelationship(parent_label='Node7353840', child_labels=['Node7353744','Python brongersmai'], edge_length=0.010752, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7353840', child_labels=[], edge_length=0.079792, taxon_label='Antaresia maculosa'),
            'Node7353840' : NodeRelationship(parent_label='Node7353680', child_labels=['Node7353712','Antaresia maculosa'], edge_length=0.004552, taxon_label=None),
            'Node7353680' : NodeRelationship(parent_label='Node7353360', child_labels=['Antaresia stimsoni','Node7353840'], edge_length=0.006907, taxon_label=None),
            'Node7353360' : NodeRelationship(parent_label='Node7353232', child_labels=['Liasis fuscus','Node7353680'], edge_length=0.008007, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7353232', child_labels=[], edge_length=0.063725, taxon_label='Morelia bredli'),
            'Node7353232' : NodeRelationship(parent_label='Node7352752', child_labels=['Node7353360','Morelia bredli'], edge_length=0.004077, taxon_label=None),
            'Node7352752' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7352880','Node7353232'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7354288', child_labels=[], edge_length=0.06111, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7354544', child_labels=[], edge_length=0.067908, taxon_label='Bothrochilus boa'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7354544', child_labels=[], edge_length=0.077196, taxon_label='Antaresia perthensis'),
            'Node7354544' : NodeRelationship(parent_label='Node7354096', child_labels=['Bothrochilus boa','Antaresia perthensis'], edge_length=0.008688, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7354096', child_labels=[], edge_length=0.1069, taxon_label='Python brongersmai'),
            'Node7354096' : NodeRelationship(parent_label='Node7354448', child_labels=['Node7354544','Python brongersmai'], edge_length=0.010367, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node7354928', child_labels=[], edge_length=0.109318, taxon_label='Python timoriensis'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7354928', child_labels=[], edge_length=0.081382, taxon_label='Antaresia maculosa'),
            'Node7354928' : NodeRelationship(parent_label='Node7354736', child_labels=['Python timoriensis','Antaresia maculosa'], edge_length=0.004129, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7354736', child_labels=[], edge_length=0.068424, taxon_label='Morelia oenpelliensis'),
            'Node7354736' : NodeRelationship(parent_label='Node7354448', child_labels=['Node7354928','Morelia oenpelliensis'], edge_length=0.007731, taxon_label=None),
            'Node7354448' : NodeRelationship(parent_label='Node7354256', child_labels=['Node7354096','Node7354736'], edge_length=0.005167, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7355344', child_labels=[], edge_length=0.079478, taxon_label='Liasis fuscus'),
            'Morelia viridis' : NodeRelationship(parent_label='Node7355408', child_labels=[], edge_length=0.062239, taxon_label='Morelia viridis'),
            'Morelia bredli' : NodeRelationship(parent_label='Node7355408', child_labels=[], edge_length=0.060484, taxon_label='Morelia bredli'),
            'Node7355408' : NodeRelationship(parent_label='Node7355344', child_labels=['Morelia viridis','Morelia bredli'], edge_length=0.013155, taxon_label=None),
            'Node7355344' : NodeRelationship(parent_label='Node7354256', child_labels=['Liasis fuscus','Node7355408'], edge_length=0.00574, taxon_label=None),
            'Node7354256' : NodeRelationship(parent_label='Node7354416', child_labels=['Node7354448','Node7355344'], edge_length=0.006988, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7354416', child_labels=[], edge_length=0.063808, taxon_label='Antaresia stimsoni'),
            'Node7354416' : NodeRelationship(parent_label='Node7354288', child_labels=['Node7354256','Antaresia stimsoni'], edge_length=0.005408, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7355600', child_labels=[], edge_length=0.069501, taxon_label='Morelia carinata'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7355600', child_labels=[], edge_length=0.073602, taxon_label='Morelia boeleni'),
            'Node7355600' : NodeRelationship(parent_label='Node7354288', child_labels=['Morelia carinata','Morelia boeleni'], edge_length=0.012121, taxon_label=None),
            'Node7354288' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7354416','Node7355600'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7355952', child_labels=[], edge_length=0.053028, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7356112', child_labels=[], edge_length=0.072831, taxon_label='Bothrochilus boa'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7356112', child_labels=[], edge_length=0.073004, taxon_label='Antaresia perthensis'),
            'Node7356112' : NodeRelationship(parent_label='Node7355920', child_labels=['Bothrochilus boa','Antaresia perthensis'], edge_length=0.012413, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7355920', child_labels=[], edge_length=0.069933, taxon_label='Morelia bredli'),
            'Node7355920' : NodeRelationship(parent_label='Node7356080', child_labels=['Node7356112','Morelia bredli'], edge_length=0.002567, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7409872', child_labels=[], edge_length=0.072339, taxon_label='Liasis fuscus'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7409872', child_labels=[], edge_length=0.062058, taxon_label='Morelia oenpelliensis'),
            'Node7409872' : NodeRelationship(parent_label='Node7409840', child_labels=['Liasis fuscus','Morelia oenpelliensis'], edge_length=0.011723, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7409840', child_labels=[], edge_length=0.065297, taxon_label='Morelia viridis'),
            'Node7409840' : NodeRelationship(parent_label='Node7356208', child_labels=['Node7409872','Morelia viridis'], edge_length=0.011519, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node7410320', child_labels=[], edge_length=0.093762, taxon_label='Python timoriensis'),
            'Morelia carinata' : NodeRelationship(parent_label='Node7410320', child_labels=[], edge_length=0.055992, taxon_label='Morelia carinata'),
            'Node7410320' : NodeRelationship(parent_label='Node7410064', child_labels=['Python timoriensis','Morelia carinata'], edge_length=0.022912, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7410064', child_labels=[], edge_length=0.079945, taxon_label='Antaresia maculosa'),
            'Node7410064' : NodeRelationship(parent_label='Node7356208', child_labels=['Node7410320','Antaresia maculosa'], edge_length=0.007647, taxon_label=None),
            'Node7356208' : NodeRelationship(parent_label='Node7409776', child_labels=['Node7409840','Node7410064'], edge_length=5.134e-05, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7409776', child_labels=[], edge_length=0.118976, taxon_label='Python brongersmai'),
            'Node7409776' : NodeRelationship(parent_label='Node7356144', child_labels=['Node7356208','Python brongersmai'], edge_length=0.004803, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7356144', child_labels=[], edge_length=0.064721, taxon_label='Antaresia stimsoni'),
            'Node7356144' : NodeRelationship(parent_label='Node7356080', child_labels=['Node7409776','Antaresia stimsoni'], edge_length=0.005969, taxon_label=None),
            'Node7356080' : NodeRelationship(parent_label='Node7355952', child_labels=['Node7355920','Node7356144'], edge_length=0.014376, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7355952', child_labels=[], edge_length=0.076894, taxon_label='Morelia boeleni'),
            'Node7355952' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7356080','Morelia boeleni'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7410896', child_labels=[], edge_length=0.06391, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7411024', child_labels=[], edge_length=0.071835, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7411216', child_labels=[], edge_length=0.075788, taxon_label='Liasis fuscus'),
            'Morelia carinata' : NodeRelationship(parent_label='Node7411216', child_labels=[], edge_length=0.07116, taxon_label='Morelia carinata'),
            'Node7411216' : NodeRelationship(parent_label='Node7411024', child_labels=['Liasis fuscus','Morelia carinata'], edge_length=0.010888, taxon_label=None),
            'Node7411024' : NodeRelationship(parent_label='Node7410896', child_labels=['Bothrochilus boa','Node7411216'], edge_length=0.00451, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7411344', child_labels=[], edge_length=0.048754, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7411344', child_labels=[], edge_length=0.069682, taxon_label='Antaresia perthensis'),
            'Node7411344' : NodeRelationship(parent_label='Node5593424', child_labels=['Antaresia stimsoni','Antaresia perthensis'], edge_length=0.012356, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5593424', child_labels=[], edge_length=0.077058, taxon_label='Antaresia maculosa'),
            'Node5593424' : NodeRelationship(parent_label='Node7411248', child_labels=['Node7411344','Antaresia maculosa'], edge_length=0.013817, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7411600', child_labels=[], edge_length=0.07026, taxon_label='Morelia viridis'),
            'Morelia bredli' : NodeRelationship(parent_label='Node7411824', child_labels=[], edge_length=0.066182, taxon_label='Morelia bredli'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7411824', child_labels=[], edge_length=0.068166, taxon_label='Morelia oenpelliensis'),
            'Node7411824' : NodeRelationship(parent_label='Node7411600', child_labels=['Morelia bredli','Morelia oenpelliensis'], edge_length=0.001443, taxon_label=None),
            'Node7411600' : NodeRelationship(parent_label='Node7411504', child_labels=['Morelia viridis','Node7411824'], edge_length=0.003147, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node7411504', child_labels=[], edge_length=0.108243, taxon_label='Python timoriensis'),
            'Node7411504' : NodeRelationship(parent_label='Node7411248', child_labels=['Node7411600','Python timoriensis'], edge_length=0.012897, taxon_label=None),
            'Node7411248' : NodeRelationship(parent_label='Node7411440', child_labels=['Node5593424','Node7411504'], edge_length=0.011944, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7412144', child_labels=[], edge_length=0.108072, taxon_label='Python brongersmai'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7412144', child_labels=[], edge_length=0.073251, taxon_label='Morelia boeleni'),
            'Node7412144' : NodeRelationship(parent_label='Node7411440', child_labels=['Python brongersmai','Morelia boeleni'], edge_length=0.00612, taxon_label=None),
            'Node7411440' : NodeRelationship(parent_label='Node7410896', child_labels=['Node7411248','Node7412144'], edge_length=0.005404, taxon_label=None),
            'Node7410896' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7411024','Node7411440'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7412368', child_labels=[], edge_length=0.06501, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7412528', child_labels=[], edge_length=0.059984, taxon_label='Bothrochilus boa'),
            'Python brongersmai' : NodeRelationship(parent_label='Node7412528', child_labels=[], edge_length=0.100555, taxon_label='Python brongersmai'),
            'Node7412528' : NodeRelationship(parent_label='Node7412016', child_labels=['Bothrochilus boa','Python brongersmai'], edge_length=0.017824, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node7412016', child_labels=[], edge_length=0.102639, taxon_label='Python timoriensis'),
            'Node7412016' : NodeRelationship(parent_label='Node7412496', child_labels=['Node7412528','Python timoriensis'], edge_length=0.011722, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7412688', child_labels=[], edge_length=0.0538, taxon_label='Morelia bredli'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7412688', child_labels=[], edge_length=0.07843, taxon_label='Antaresia maculosa'),
            'Node7412688' : NodeRelationship(parent_label='Node7412496', child_labels=['Morelia bredli','Antaresia maculosa'], edge_length=0.014836, taxon_label=None),
            'Node7412496' : NodeRelationship(parent_label='Node7412368', child_labels=['Node7412016','Node7412688'], edge_length=0.003644, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7412912', child_labels=[], edge_length=0.068438, taxon_label='Liasis fuscus'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7412912', child_labels=[], edge_length=0.066815, taxon_label='Morelia oenpelliensis'),
            'Node7412912' : NodeRelationship(parent_label='Node7412624', child_labels=['Liasis fuscus','Morelia oenpelliensis'], edge_length=0.015326, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7413456', child_labels=[], edge_length=0.057007, taxon_label='Antaresia stimsoni'),
            'Morelia carinata' : NodeRelationship(parent_label='Node7413648', child_labels=[], edge_length=0.062766, taxon_label='Morelia carinata'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7413648', child_labels=[], edge_length=0.079145, taxon_label='Morelia boeleni'),
            'Node7413648' : NodeRelationship(parent_label='Node7413456', child_labels=['Morelia carinata','Morelia boeleni'], edge_length=0.013392, taxon_label=None),
            'Node7413456' : NodeRelationship(parent_label='Node7413424', child_labels=['Antaresia stimsoni','Node7413648'], edge_length=0.010508, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7413424', child_labels=[], edge_length=0.083409, taxon_label='Antaresia perthensis'),
            'Node7413424' : NodeRelationship(parent_label='Node7412624', child_labels=['Node7413456','Antaresia perthensis'], edge_length=0.005897, taxon_label=None),
            'Node7412624' : NodeRelationship(parent_label='Node7413168', child_labels=['Node7412912','Node7413424'], edge_length=0.0, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7413168', child_labels=[], edge_length=0.072091, taxon_label='Morelia viridis'),
            'Node7413168' : NodeRelationship(parent_label='Node7412368', child_labels=['Node7412624','Morelia viridis'], edge_length=0.008523, taxon_label=None),
            'Node7412368' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7412496','Node7413168'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7504176', child_labels=[], edge_length=0.062297, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7504304', child_labels=[], edge_length=0.069165, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7504304', child_labels=[], edge_length=0.073495, taxon_label='Liasis fuscus'),
            'Node7504304' : NodeRelationship(parent_label='Node7504176', child_labels=['Bothrochilus boa','Liasis fuscus'], edge_length=0.009899, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7504688', child_labels=[], edge_length=0.061022, taxon_label='Antaresia stimsoni'),
            'Morelia viridis' : NodeRelationship(parent_label='Node7505008', child_labels=[], edge_length=0.071291, taxon_label='Morelia viridis'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7505104', child_labels=[], edge_length=0.098776, taxon_label='Python timoriensis'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7505104', child_labels=[], edge_length=0.075069, taxon_label='Morelia boeleni'),
            'Node7505104' : NodeRelationship(parent_label='Node7505008', child_labels=['Python timoriensis','Morelia boeleni'], edge_length=0.017245, taxon_label=None),
            'Node7505008' : NodeRelationship(parent_label='Node7504688', child_labels=['Morelia viridis','Node7505104'], edge_length=0.004698, taxon_label=None),
            'Node7504688' : NodeRelationship(parent_label='Node7504624', child_labels=['Antaresia stimsoni','Node7505008'], edge_length=0.003603, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7505136', child_labels=[], edge_length=0.075494, taxon_label='Antaresia maculosa'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7505136', child_labels=[], edge_length=0.063757, taxon_label='Morelia oenpelliensis'),
            'Node7505136' : NodeRelationship(parent_label='Node7504624', child_labels=['Antaresia maculosa','Morelia oenpelliensis'], edge_length=0.017571, taxon_label=None),
            'Node7504624' : NodeRelationship(parent_label='Node7504368', child_labels=['Node7504688','Node7505136'], edge_length=0.000386, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7504368', child_labels=[], edge_length=0.070614, taxon_label='Morelia bredli'),
            'Node7504368' : NodeRelationship(parent_label='Node7504592', child_labels=['Node7504624','Morelia bredli'], edge_length=0.00713, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7504592', child_labels=[], edge_length=0.079534, taxon_label='Morelia carinata'),
            'Node7504592' : NodeRelationship(parent_label='Node7504528', child_labels=['Node7504368','Morelia carinata'], edge_length=0.010249, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7504528', child_labels=[], edge_length=0.114381, taxon_label='Python brongersmai'),
            'Node7504528' : NodeRelationship(parent_label='Node5593232', child_labels=['Node7504592','Python brongersmai'], edge_length=0.003672, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node5593232', child_labels=[], edge_length=0.080896, taxon_label='Antaresia perthensis'),
            'Node5593232' : NodeRelationship(parent_label='Node7504176', child_labels=['Node7504528','Antaresia perthensis'], edge_length=0.003212, taxon_label=None),
            'Node7504176' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7504304','Node5593232'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7505680', child_labels=[], edge_length=0.05678, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7506160', child_labels=[], edge_length=0.065928, taxon_label='Bothrochilus boa'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7506160', child_labels=[], edge_length=0.10077, taxon_label='Python timoriensis'),
            'Node7506160' : NodeRelationship(parent_label='Node7506096', child_labels=['Bothrochilus boa','Python timoriensis'], edge_length=0.013813, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7506096', child_labels=[], edge_length=0.083335, taxon_label='Antaresia perthensis'),
            'Node7506096' : NodeRelationship(parent_label='Node7506032', child_labels=['Node7506160','Antaresia perthensis'], edge_length=0.004377, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7506352', child_labels=[], edge_length=0.073392, taxon_label='Antaresia maculosa'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7506352', child_labels=[], edge_length=0.065748, taxon_label='Morelia oenpelliensis'),
            'Node7506352' : NodeRelationship(parent_label='Node7506032', child_labels=['Antaresia maculosa','Morelia oenpelliensis'], edge_length=0.016433, taxon_label=None),
            'Node7506032' : NodeRelationship(parent_label='Node7505968', child_labels=['Node7506096','Node7506352'], edge_length=0.005104, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7505968', child_labels=[], edge_length=0.07646, taxon_label='Morelia carinata'),
            'Node7505968' : NodeRelationship(parent_label='Node7505872', child_labels=['Node7506032','Morelia carinata'], edge_length=0.001866, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7505872', child_labels=[], edge_length=0.06129, taxon_label='Antaresia stimsoni'),
            'Node7505872' : NodeRelationship(parent_label='Node7412272', child_labels=['Node7505968','Antaresia stimsoni'], edge_length=0.012314, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7412272', child_labels=[], edge_length=0.081982, taxon_label='Morelia boeleni'),
            'Node7412272' : NodeRelationship(parent_label='Node7505584', child_labels=['Node7505872','Morelia boeleni'], edge_length=0.002833, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7506768', child_labels=[], edge_length=0.064242, taxon_label='Morelia viridis'),
            'Python brongersmai' : NodeRelationship(parent_label='Node7506768', child_labels=[], edge_length=0.109434, taxon_label='Python brongersmai'),
            'Node7506768' : NodeRelationship(parent_label='Node7505584', child_labels=['Morelia viridis','Python brongersmai'], edge_length=0.010734, taxon_label=None),
            'Node7505584' : NodeRelationship(parent_label='Node7505840', child_labels=['Node7412272','Node7506768'], edge_length=0.00182, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7505840', child_labels=[], edge_length=0.079156, taxon_label='Liasis fuscus'),
            'Node7505840' : NodeRelationship(parent_label='Node7505680', child_labels=['Node7505584','Liasis fuscus'], edge_length=0.014479, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7505680', child_labels=[], edge_length=0.060374, taxon_label='Morelia bredli'),
            'Node7505680' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7505840','Morelia bredli'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7507376', child_labels=[], edge_length=0.063004, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7507504', child_labels=[], edge_length=0.061413, taxon_label='Bothrochilus boa'),
            'Python brongersmai' : NodeRelationship(parent_label='Node7507504', child_labels=[], edge_length=0.099507, taxon_label='Python brongersmai'),
            'Node7507504' : NodeRelationship(parent_label='Node7507376', child_labels=['Bothrochilus boa','Python brongersmai'], edge_length=0.018059, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7507888', child_labels=[], edge_length=0.075877, taxon_label='Liasis fuscus'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7507888', child_labels=[], edge_length=0.071375, taxon_label='Antaresia maculosa'),
            'Node7507888' : NodeRelationship(parent_label='Node7507824', child_labels=['Liasis fuscus','Antaresia maculosa'], edge_length=0.014459, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7507824', child_labels=[], edge_length=0.054698, taxon_label='Antaresia stimsoni'),
            'Node7507824' : NodeRelationship(parent_label='Node7507568', child_labels=['Node7507888','Antaresia stimsoni'], edge_length=0.002182, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7561296', child_labels=[], edge_length=0.062291, taxon_label='Morelia bredli'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7561744', child_labels=[], edge_length=0.074463, taxon_label='Antaresia perthensis'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7561744', child_labels=[], edge_length=0.103887, taxon_label='Python timoriensis'),
            'Node7561744' : NodeRelationship(parent_label='Node7561296', child_labels=['Antaresia perthensis','Python timoriensis'], edge_length=0.00521, taxon_label=None),
            'Node7561296' : NodeRelationship(parent_label='Node7507568', child_labels=['Morelia bredli','Node7561744'], edge_length=0.011197, taxon_label=None),
            'Node7507568' : NodeRelationship(parent_label='Node7507792', child_labels=['Node7507824','Node7561296'], edge_length=0.010763, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7507792', child_labels=[], edge_length=0.070762, taxon_label='Morelia carinata'),
            'Node7507792' : NodeRelationship(parent_label='Node7507536', child_labels=['Node7507568','Morelia carinata'], edge_length=0.016533, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7507536', child_labels=[], edge_length=0.068818, taxon_label='Morelia oenpelliensis'),
            'Node7507536' : NodeRelationship(parent_label='Node7507728', child_labels=['Node7507792','Morelia oenpelliensis'], edge_length=0.003714, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7561520', child_labels=[], edge_length=0.071307, taxon_label='Morelia viridis'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7561520', child_labels=[], edge_length=0.073392, taxon_label='Morelia boeleni'),
            'Node7561520' : NodeRelationship(parent_label='Node7507728', child_labels=['Morelia viridis','Morelia boeleni'], edge_length=0.010478, taxon_label=None),
            'Node7507728' : NodeRelationship(parent_label='Node7507376', child_labels=['Node7507536','Node7561520'], edge_length=0.00238, taxon_label=None),
            'Node7507376' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7507504','Node7507728'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7562320', child_labels=[], edge_length=0.053901, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7562480', child_labels=[], edge_length=0.07927, taxon_label='Bothrochilus boa'),
            'Morelia bredli' : NodeRelationship(parent_label='Node7562832', child_labels=[], edge_length=0.053956, taxon_label='Morelia bredli'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7562832', child_labels=[], edge_length=0.098874, taxon_label='Python timoriensis'),
            'Node7562832' : NodeRelationship(parent_label='Node7562480', child_labels=['Morelia bredli','Python timoriensis'], edge_length=0.016769, taxon_label=None),
            'Node7562480' : NodeRelationship(parent_label='Node7355728', child_labels=['Bothrochilus boa','Node7562832'], edge_length=0.003041, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7562864', child_labels=[], edge_length=0.075223, taxon_label='Liasis fuscus'),
            'Morelia viridis' : NodeRelationship(parent_label='Node7562864', child_labels=[], edge_length=0.065805, taxon_label='Morelia viridis'),
            'Node7562864' : NodeRelationship(parent_label='Node7563056', child_labels=['Liasis fuscus','Morelia viridis'], edge_length=0.010712, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7563056', child_labels=[], edge_length=0.075934, taxon_label='Morelia carinata'),
            'Node7563056' : NodeRelationship(parent_label='Node7355728', child_labels=['Node7562864','Morelia carinata'], edge_length=0.008083, taxon_label=None),
            'Node7355728' : NodeRelationship(parent_label='Node7562032', child_labels=['Node7562480','Node7563056'], edge_length=0.00777, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7563184', child_labels=[], edge_length=0.082425, taxon_label='Antaresia perthensis'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7563024', child_labels=[], edge_length=0.075429, taxon_label='Antaresia maculosa'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7563024', child_labels=[], edge_length=0.0787, taxon_label='Morelia boeleni'),
            'Node7563024' : NodeRelationship(parent_label='Node7563280', child_labels=['Antaresia maculosa','Morelia boeleni'], edge_length=0.010583, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7563280', child_labels=[], edge_length=0.070289, taxon_label='Morelia oenpelliensis'),
            'Node7563280' : NodeRelationship(parent_label='Node7563184', child_labels=['Node7563024','Morelia oenpelliensis'], edge_length=0.008508, taxon_label=None),
            'Node7563184' : NodeRelationship(parent_label='Node7562032', child_labels=['Antaresia perthensis','Node7563280'], edge_length=0.00373, taxon_label=None),
            'Node7562032' : NodeRelationship(parent_label='Node7562448', child_labels=['Node7355728','Node7563184'], edge_length=0.002034, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7562448', child_labels=[], edge_length=0.064652, taxon_label='Antaresia stimsoni'),
            'Node7562448' : NodeRelationship(parent_label='Node7562320', child_labels=['Node7562032','Antaresia stimsoni'], edge_length=0.013947, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7562320', child_labels=[], edge_length=0.107868, taxon_label='Python brongersmai'),
            'Node7562320' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7562448','Python brongersmai'], edge_length=None, taxon_label=None),
        },
    ]
    return treelist_node_references

def reference_taxon_set():
    taxon_set = dendropy.TaxonSet()
    taxon_set.new_taxon(label="Aspidites ramsayi", oid="Taxon5593072")
    taxon_set.new_taxon(label="Bothrochilus boa", oid="Taxon5594064")
    taxon_set.new_taxon(label="Liasis fuscus", oid="Taxon5593680")
    taxon_set.new_taxon(label="Morelia boeleni", oid="Taxon5594448")
    taxon_set.new_taxon(label="Morelia viridis", oid="Taxon5594544")
    taxon_set.new_taxon(label="Antaresia maculosa", oid="Taxon5594608")
    taxon_set.new_taxon(label="Antaresia stimsoni", oid="Taxon5594480")
    taxon_set.new_taxon(label="Python timoriensis", oid="Taxon5594768")
    taxon_set.new_taxon(label="Morelia oenpelliensis", oid="Taxon5594832")
    taxon_set.new_taxon(label="Morelia bredli", oid="Taxon5594992")
    taxon_set.new_taxon(label="Antaresia perthensis", oid="Taxon5595120")
    taxon_set.new_taxon(label="Python brongersmai", oid="Taxon5595024")
    taxon_set.new_taxon(label="Morelia carinata", oid="Taxon7352592")
    return taxon_set

def reference_tree_list():
    taxon_set = reference_taxon_set()
    tree_list = dendropy.TreeList(label=None, oid="TreeList4398560", taxon_set=taxon_set)
    assert tree_list.taxon_set is taxon_set
    assert len(tree_list.taxon_set) == 13
    assert tree_list.taxon_set.has_taxa(labels=["Aspidites ramsayi",
            "Bothrochilus boa", "Liasis fuscus", "Morelia boeleni", "Morelia viridis",
            "Antaresia maculosa", "Antaresia stimsoni", "Python timoriensis",
            "Morelia oenpelliensis", "Morelia bredli", "Antaresia perthensis",
            "Python brongersmai", "Morelia carinata"]), \
            [t.label for t in tree_list.taxon_set]
    tax_5593072 = taxon_set.require_taxon(label="Aspidites ramsayi", oid="Taxon5593072")
    tax_5594064 = taxon_set.require_taxon(label="Bothrochilus boa", oid="Taxon5594064")
    tax_5593680 = taxon_set.require_taxon(label="Liasis fuscus", oid="Taxon5593680")
    tax_5594448 = taxon_set.require_taxon(label="Morelia boeleni", oid="Taxon5594448")
    tax_5594544 = taxon_set.require_taxon(label="Morelia viridis", oid="Taxon5594544")
    tax_5594608 = taxon_set.require_taxon(label="Antaresia maculosa", oid="Taxon5594608")
    tax_5594480 = taxon_set.require_taxon(label="Antaresia stimsoni", oid="Taxon5594480")
    tax_5594768 = taxon_set.require_taxon(label="Python timoriensis", oid="Taxon5594768")
    tax_5594832 = taxon_set.require_taxon(label="Morelia oenpelliensis", oid="Taxon5594832")
    tax_5594992 = taxon_set.require_taxon(label="Morelia bredli", oid="Taxon5594992")
    tax_5595120 = taxon_set.require_taxon(label="Antaresia perthensis", oid="Taxon5595120")
    tax_5595024 = taxon_set.require_taxon(label="Python brongersmai", oid="Taxon5595024")
    tax_7352592 = taxon_set.require_taxon(label="Morelia carinata", oid="Taxon7352592")
    assert len(tree_list.taxon_set) == 13
    tree_5593360 = dendropy.Tree(label="Tree01", taxon_set=tree_list.taxon_set, oid="Tree5593360")
    tree_list.append(tree_5593360, reindex_taxa=False)
    tree_5593360.seed_node.oid = 'Node5593456'
    nd_5593520 = tree_5593360.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.056823, oid="Node5593520")
    nd_5593520.edge.oid = "Edge5593552"
    nd_5593584 = tree_5593360.seed_node.new_child(label="Node5593584", taxon=None, edge_length=0.011159, oid="Node5593584")
    nd_5593584.edge.oid = "Edge5593200"
    nd_7352624 = tree_5593360.seed_node.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.079321, oid="Node7352624")
    nd_7352624.edge.oid = "Edge7352400"
    nd_5593648 = nd_5593584.new_child(label="Node5593648", taxon=None, edge_length=0.003495, oid="Node5593648")
    nd_5593648.edge.oid = "Edge5593712"
    nd_5594896 = nd_5593584.new_child(label="Node5594896", taxon=None, edge_length=0.018613, oid="Node5594896")
    nd_5594896.edge.oid = "Edge5595056"
    nd_5593616 = nd_5593648.new_child(label="Node5593616", taxon=None, edge_length=0.007655, oid="Node5593616")
    nd_5593616.edge.oid = "Edge5593776"
    nd_5594960 = nd_5593648.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.065046, oid="Node5594960")
    nd_5594960.edge.oid = "Edge5594928"
    nd_5593328 = nd_5593616.new_child(label="Node5593328", taxon=None, edge_length=0.002264, oid="Node5593328")
    nd_5593328.edge.oid = "Edge5593840"
    nd_5594704 = nd_5593616.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.074047, oid="Node5594704")
    nd_5594704.edge.oid = "Edge5594864"
    nd_5593744 = nd_5593328.new_child(label="Node5593744", taxon=None, edge_length=0.0, oid="Node5593744")
    nd_5593744.edge.oid = "Edge5593904"
    nd_5594640 = nd_5593328.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.113605, oid="Node5594640")
    nd_5594640.edge.oid = "Edge5594800"
    nd_5593808 = nd_5593744.new_child(label="Node5593808", taxon=None, edge_length=0.007681, oid="Node5593808")
    nd_5593808.edge.oid = "Edge5593968"
    nd_5594512 = nd_5593744.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.061996, oid="Node5594512")
    nd_5594512.edge.oid = "Edge5594672"
    nd_5593872 = nd_5593808.new_child(label="Node5593872", taxon=None, edge_length=0.010103, oid="Node5593872")
    nd_5593872.edge.oid = "Edge5594032"
    nd_5594352 = nd_5593808.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.080749, oid="Node5594352")
    nd_5594352.edge.oid = "Edge5594416"
    nd_5593936 = nd_5593872.new_child(label="Node5593936", taxon=None, edge_length=0.011155, oid="Node5593936")
    nd_5593936.edge.oid = "Edge5594096"
    nd_5594576 = nd_5593872.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.066566, oid="Node5594576")
    nd_5594576.edge.oid = "Edge5594320"
    nd_5594000 = nd_5593936.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.073578, oid="Node5594000")
    nd_5594000.edge.oid = "Edge5594160"
    nd_5594256 = nd_5593936.new_child(label="Node5594256", taxon=None, edge_length=0.003063, oid="Node5594256")
    nd_5594256.edge.oid = "Edge5594128"
    nd_5594288 = nd_5594256.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.074347, oid="Node5594288")
    nd_5594288.edge.oid = "Edge5594192"
    nd_5594384 = nd_5594256.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.07628, oid="Node5594384")
    nd_5594384.edge.oid = "Edge5594224"
    nd_5594736 = nd_5594896.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.065004, oid="Node5594736")
    nd_5594736.edge.oid = "Edge7352368"
    nd_5595088 = nd_5594896.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.107706, oid="Node5595088")
    nd_5595088.edge.oid = "Edge7352464"
    tree_7352496 = dendropy.Tree(label="Tree02", taxon_set=tree_list.taxon_set, oid="Tree7352496")
    tree_list.append(tree_7352496, reindex_taxa=False)
    tree_7352496.seed_node.oid = 'Node7352752'
    nd_7352816 = tree_7352496.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.065297, oid="Node7352816")
    nd_7352816.edge.oid = "Edge7352848"
    nd_7352880 = tree_7352496.seed_node.new_child(label="Node7352880", taxon=None, edge_length=0.001444, oid="Node7352880")
    nd_7352880.edge.oid = "Edge7352720"
    nd_7353232 = tree_7352496.seed_node.new_child(label="Node7353232", taxon=None, edge_length=0.004077, oid="Node7353232")
    nd_7353232.edge.oid = "Edge7353168"
    nd_7352688 = nd_7352880.new_child(label="Node7352688", taxon=None, edge_length=0.010946, oid="Node7352688")
    nd_7352688.edge.oid = "Edge7352976"
    nd_7352944 = nd_7352880.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.079424, oid="Node7352944")
    nd_7352944.edge.oid = "Edge7353392"
    nd_5593168 = nd_7352688.new_child(label="Node5593168", taxon=None, edge_length=0.016022, oid="Node5593168")
    nd_5593168.edge.oid = "Edge7353040"
    nd_7353424 = nd_7352688.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.066015, oid="Node7353424")
    nd_7353424.edge.oid = "Edge7353328"
    nd_7352912 = nd_5593168.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.073565, oid="Node7352912")
    nd_7352912.edge.oid = "Edge7353104"
    nd_7353200 = nd_5593168.new_child(label="Node7353200", taxon=None, edge_length=0.007879, oid="Node7353200")
    nd_7353200.edge.oid = "Edge7353136"
    nd_7353072 = nd_7353200.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.069885, oid="Node7353072")
    nd_7353072.edge.oid = "Edge7353008"
    nd_7353296 = nd_7353200.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.06247, oid="Node7353296")
    nd_7353296.edge.oid = "Edge7353264"
    nd_7353360 = nd_7353232.new_child(label="Node7353360", taxon=None, edge_length=0.008007, oid="Node7353360")
    nd_7353360.edge.oid = "Edge7353520"
    nd_7354032 = nd_7353232.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.063725, oid="Node7354032")
    nd_7354032.edge.oid = "Edge7354064"
    nd_7353552 = nd_7353360.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.082049, oid="Node7353552")
    nd_7353552.edge.oid = "Edge7353616"
    nd_7353680 = nd_7353360.new_child(label="Node7353680", taxon=None, edge_length=0.006907, oid="Node7353680")
    nd_7353680.edge.oid = "Edge7353648"
    nd_7353584 = nd_7353680.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.056894, oid="Node7353584")
    nd_7353584.edge.oid = "Edge7353488"
    nd_7353840 = nd_7353680.new_child(label="Node7353840", taxon=None, edge_length=0.004552, oid="Node7353840")
    nd_7353840.edge.oid = "Edge7353776"
    nd_7353712 = nd_7353840.new_child(label="Node7353712", taxon=None, edge_length=0.010752, oid="Node7353712")
    nd_7353712.edge.oid = "Edge7353456"
    nd_7352528 = nd_7353840.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.079792, oid="Node7352528")
    nd_7352528.edge.oid = "Edge7353872"
    nd_7353744 = nd_7353712.new_child(label="Node7353744", taxon=None, edge_length=0.015909, oid="Node7353744")
    nd_7353744.edge.oid = "Edge7353936"
    nd_5593296 = nd_7353712.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.114385, oid="Node5593296")
    nd_5593296.edge.oid = "Edge7353968"
    nd_7353808 = nd_7353744.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.092077, oid="Node7353808")
    nd_7353808.edge.oid = "Edge5593392"
    nd_7352560 = nd_7353744.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.0579, oid="Node7352560")
    nd_7352560.edge.oid = "Edge7353904"
    tree_7354128 = dendropy.Tree(label="Tree03", taxon_set=tree_list.taxon_set, oid="Tree7354128")
    tree_list.append(tree_7354128, reindex_taxa=False)
    tree_7354128.seed_node.oid = 'Node7354288'
    nd_7354352 = tree_7354128.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.06111, oid="Node7354352")
    nd_7354352.edge.oid = "Edge7354384"
    nd_7354416 = tree_7354128.seed_node.new_child(label="Node7354416", taxon=None, edge_length=0.005408, oid="Node7354416")
    nd_7354416.edge.oid = "Edge7352432"
    nd_7355600 = tree_7354128.seed_node.new_child(label="Node7355600", taxon=None, edge_length=0.012121, oid="Node7355600")
    nd_7355600.edge.oid = "Edge7355472"
    nd_7354256 = nd_7354416.new_child(label="Node7354256", taxon=None, edge_length=0.006988, oid="Node7354256")
    nd_7354256.edge.oid = "Edge7354512"
    nd_7355440 = nd_7354416.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.063808, oid="Node7355440")
    nd_7355440.edge.oid = "Edge7355568"
    nd_7354448 = nd_7354256.new_child(label="Node7354448", taxon=None, edge_length=0.005167, oid="Node7354448")
    nd_7354448.edge.oid = "Edge7354576"
    nd_7355344 = nd_7354256.new_child(label="Node7355344", taxon=None, edge_length=0.00574, oid="Node7355344")
    nd_7355344.edge.oid = "Edge7354672"
    nd_7354096 = nd_7354448.new_child(label="Node7354096", taxon=None, edge_length=0.010367, oid="Node7354096")
    nd_7354096.edge.oid = "Edge7354640"
    nd_7354736 = nd_7354448.new_child(label="Node7354736", taxon=None, edge_length=0.007731, oid="Node7354736")
    nd_7354736.edge.oid = "Edge7354480"
    nd_7354544 = nd_7354096.new_child(label="Node7354544", taxon=None, edge_length=0.008688, oid="Node7354544")
    nd_7354544.edge.oid = "Edge7354704"
    nd_7354896 = nd_7354096.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.1069, oid="Node7354896")
    nd_7354896.edge.oid = "Edge7354832"
    nd_7354608 = nd_7354544.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.067908, oid="Node7354608")
    nd_7354608.edge.oid = "Edge7354768"
    nd_7354864 = nd_7354544.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.077196, oid="Node7354864")
    nd_7354864.edge.oid = "Edge7354800"
    nd_7354928 = nd_7354736.new_child(label="Node7354928", taxon=None, edge_length=0.004129, oid="Node7354928")
    nd_7354928.edge.oid = "Edge7354992"
    nd_7355216 = nd_7354736.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.068424, oid="Node7355216")
    nd_7355216.edge.oid = "Edge7355152"
    nd_7355024 = nd_7354928.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.109318, oid="Node7355024")
    nd_7355024.edge.oid = "Edge7355088"
    nd_7355184 = nd_7354928.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.081382, oid="Node7355184")
    nd_7355184.edge.oid = "Edge7355120"
    nd_7355280 = nd_7355344.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.079478, oid="Node7355280")
    nd_7355280.edge.oid = "Edge7355056"
    nd_7355408 = nd_7355344.new_child(label="Node7355408", taxon=None, edge_length=0.013155, oid="Node7355408")
    nd_7355408.edge.oid = "Edge7355376"
    nd_7355312 = nd_7355408.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.062239, oid="Node7355312")
    nd_7355312.edge.oid = "Edge7354960"
    nd_7355536 = nd_7355408.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.060484, oid="Node7355536")
    nd_7355536.edge.oid = "Edge7355504"
    nd_7355664 = nd_7355600.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.069501, oid="Node7355664")
    nd_7355664.edge.oid = "Edge7355696"
    nd_7355824 = nd_7355600.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.073602, oid="Node7355824")
    nd_7355824.edge.oid = "Edge7355760"
    tree_7355792 = dendropy.Tree(label="Tree04", taxon_set=tree_list.taxon_set, oid="Tree7355792")
    tree_list.append(tree_7355792, reindex_taxa=False)
    tree_7355792.seed_node.oid = 'Node7355952'
    nd_7356016 = tree_7355792.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.053028, oid="Node7356016")
    nd_7356016.edge.oid = "Edge7356048"
    nd_7356080 = tree_7355792.seed_node.new_child(label="Node7356080", taxon=None, edge_length=0.014376, oid="Node7356080")
    nd_7356080.edge.oid = "Edge7354224"
    nd_7410800 = tree_7355792.seed_node.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.076894, oid="Node7410800")
    nd_7410800.edge.oid = "Edge7410288"
    nd_7355920 = nd_7356080.new_child(label="Node7355920", taxon=None, edge_length=0.002567, oid="Node7355920")
    nd_7355920.edge.oid = "Edge7356176"
    nd_7356144 = nd_7356080.new_child(label="Node7356144", taxon=None, edge_length=0.005969, oid="Node7356144")
    nd_7356144.edge.oid = "Edge7409808"
    nd_7356112 = nd_7355920.new_child(label="Node7356112", taxon=None, edge_length=0.012413, oid="Node7356112")
    nd_7356112.edge.oid = "Edge7356240"
    nd_7356368 = nd_7355920.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.069933, oid="Node7356368")
    nd_7356368.edge.oid = "Edge7356272"
    nd_7355632 = nd_7356112.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.072831, oid="Node7355632")
    nd_7355632.edge.oid = "Edge7356304"
    nd_7356400 = nd_7356112.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.073004, oid="Node7356400")
    nd_7356400.edge.oid = "Edge7356336"
    nd_7409776 = nd_7356144.new_child(label="Node7409776", taxon=None, edge_length=0.004803, oid="Node7409776")
    nd_7409776.edge.oid = "Edge7409744"
    nd_7410384 = nd_7356144.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.064721, oid="Node7410384")
    nd_7410384.edge.oid = "Edge7410576"
    nd_7356208 = nd_7409776.new_child(label="Node7356208", taxon=None, edge_length=5.134e-05, oid="Node7356208")
    nd_7356208.edge.oid = "Edge7409904"
    nd_7410672 = nd_7409776.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.118976, oid="Node7410672")
    nd_7410672.edge.oid = "Edge7410608"
    nd_7409840 = nd_7356208.new_child(label="Node7409840", taxon=None, edge_length=0.011519, oid="Node7409840")
    nd_7409840.edge.oid = "Edge7409968"
    nd_7410064 = nd_7356208.new_child(label="Node7410064", taxon=None, edge_length=0.007647, oid="Node7410064")
    nd_7410064.edge.oid = "Edge7409712"
    nd_7409872 = nd_7409840.new_child(label="Node7409872", taxon=None, edge_length=0.011723, oid="Node7409872")
    nd_7409872.edge.oid = "Edge7410032"
    nd_7410224 = nd_7409840.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.065297, oid="Node7410224")
    nd_7410224.edge.oid = "Edge7410192"
    nd_7409936 = nd_7409872.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.072339, oid="Node7409936")
    nd_7409936.edge.oid = "Edge7410096"
    nd_7410160 = nd_7409872.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.062058, oid="Node7410160")
    nd_7410160.edge.oid = "Edge7410128"
    nd_7410320 = nd_7410064.new_child(label="Node7410320", taxon=None, edge_length=0.022912, oid="Node7410320")
    nd_7410320.edge.oid = "Edge7410000"
    nd_7410544 = nd_7410064.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.079945, oid="Node7410544")
    nd_7410544.edge.oid = "Edge7410480"
    nd_7410352 = nd_7410320.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.093762, oid="Node7410352")
    nd_7410352.edge.oid = "Edge7410416"
    nd_7410512 = nd_7410320.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.055992, oid="Node7410512")
    nd_7410512.edge.oid = "Edge7410448"
    tree_7410640 = dendropy.Tree(label="Tree05", taxon_set=tree_list.taxon_set, oid="Tree7410640")
    tree_list.append(tree_7410640, reindex_taxa=False)
    tree_7410640.seed_node.oid = 'Node7410896'
    nd_7410960 = tree_7410640.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.06391, oid="Node7410960")
    nd_7410960.edge.oid = "Edge7410992"
    nd_7411024 = tree_7410640.seed_node.new_child(label="Node7411024", taxon=None, edge_length=0.00451, oid="Node7411024")
    nd_7411024.edge.oid = "Edge7410864"
    nd_7411440 = tree_7410640.seed_node.new_child(label="Node7411440", taxon=None, edge_length=0.005404, oid="Node7411440")
    nd_7411440.edge.oid = "Edge7411088"
    nd_7410256 = nd_7411024.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.071835, oid="Node7410256")
    nd_7410256.edge.oid = "Edge7411120"
    nd_7411216 = nd_7411024.new_child(label="Node7411216", taxon=None, edge_length=0.010888, oid="Node7411216")
    nd_7411216.edge.oid = "Edge7411152"
    nd_7411056 = nd_7411216.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.075788, oid="Node7411056")
    nd_7411056.edge.oid = "Edge7411184"
    nd_7411312 = nd_7411216.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.07116, oid="Node7411312")
    nd_7411312.edge.oid = "Edge7411280"
    nd_7411248 = nd_7411440.new_child(label="Node7411248", taxon=None, edge_length=0.011944, oid="Node7411248")
    nd_7411248.edge.oid = "Edge7355248"
    nd_7412144 = nd_7411440.new_child(label="Node7412144", taxon=None, edge_length=0.00612, oid="Node7412144")
    nd_7412144.edge.oid = "Edge7412048"
    nd_5593424 = nd_7411248.new_child(label="Node5593424", taxon=None, edge_length=0.013817, oid="Node5593424")
    nd_5593424.edge.oid = "Edge7354160"
    nd_7411504 = nd_7411248.new_child(label="Node7411504", taxon=None, edge_length=0.012897, oid="Node7411504")
    nd_7411504.edge.oid = "Edge7411408"
    nd_7411344 = nd_5593424.new_child(label="Node7411344", taxon=None, edge_length=0.012356, oid="Node7411344")
    nd_7411344.edge.oid = "Edge7410832"
    nd_7354192 = nd_5593424.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.077058, oid="Node7354192")
    nd_7354192.edge.oid = "Edge7410768"
    nd_7355888 = nd_7411344.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.048754, oid="Node7355888")
    nd_7355888.edge.oid = "Edge7411376"
    nd_7411536 = nd_7411344.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.069682, oid="Node7411536")
    nd_7411536.edge.oid = "Edge7411472"
    nd_7411600 = nd_7411504.new_child(label="Node7411600", taxon=None, edge_length=0.003147, oid="Node7411600")
    nd_7411600.edge.oid = "Edge7411664"
    nd_7412080 = nd_7411504.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.108243, oid="Node7412080")
    nd_7412080.edge.oid = "Edge7411952"
    nd_7411696 = nd_7411600.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.07026, oid="Node7411696")
    nd_7411696.edge.oid = "Edge7411760"
    nd_7411824 = nd_7411600.new_child(label="Node7411824", taxon=None, edge_length=0.001443, oid="Node7411824")
    nd_7411824.edge.oid = "Edge7411792"
    nd_7411728 = nd_7411824.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.066182, oid="Node7411728")
    nd_7411728.edge.oid = "Edge7411632"
    nd_7411984 = nd_7411824.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.068166, oid="Node7411984")
    nd_7411984.edge.oid = "Edge7411920"
    nd_7411856 = nd_7412144.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.108072, oid="Node7411856")
    nd_7411856.edge.oid = "Edge7411568"
    nd_7412240 = nd_7412144.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.073251, oid="Node7412240")
    nd_7412240.edge.oid = "Edge7412176"
    tree_7412208 = dendropy.Tree(label="Tree06", taxon_set=tree_list.taxon_set, oid="Tree7412208")
    tree_list.append(tree_7412208, reindex_taxa=False)
    tree_7412208.seed_node.oid = 'Node7412368'
    nd_7412432 = tree_7412208.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.06501, oid="Node7412432")
    nd_7412432.edge.oid = "Edge7412464"
    nd_7412496 = tree_7412208.seed_node.new_child(label="Node7412496", taxon=None, edge_length=0.003644, oid="Node7412496")
    nd_7412496.edge.oid = "Edge7410704"
    nd_7413168 = tree_7412208.seed_node.new_child(label="Node7413168", taxon=None, edge_length=0.008523, oid="Node7413168")
    nd_7413168.edge.oid = "Edge7413040"
    nd_7412016 = nd_7412496.new_child(label="Node7412016", taxon=None, edge_length=0.011722, oid="Node7412016")
    nd_7412016.edge.oid = "Edge7412592"
    nd_7412688 = nd_7412496.new_child(label="Node7412688", taxon=None, edge_length=0.014836, oid="Node7412688")
    nd_7412688.edge.oid = "Edge7412560"
    nd_7412528 = nd_7412016.new_child(label="Node7412528", taxon=None, edge_length=0.017824, oid="Node7412528")
    nd_7412528.edge.oid = "Edge7412656"
    nd_7412848 = nd_7412016.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.102639, oid="Node7412848")
    nd_7412848.edge.oid = "Edge7412784"
    nd_7411888 = nd_7412528.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.059984, oid="Node7411888")
    nd_7411888.edge.oid = "Edge7412720"
    nd_7412816 = nd_7412528.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.100555, oid="Node7412816")
    nd_7412816.edge.oid = "Edge7412752"
    nd_7412880 = nd_7412688.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.0538, oid="Node7412880")
    nd_7412880.edge.oid = "Edge7412944"
    nd_7413072 = nd_7412688.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.07843, oid="Node7413072")
    nd_7413072.edge.oid = "Edge7413008"
    nd_7412624 = nd_7413168.new_child(label="Node7412624", taxon=None, edge_length=0.0, oid="Node7412624")
    nd_7412624.edge.oid = "Edge7413104"
    nd_7413680 = nd_7413168.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.072091, oid="Node7413680")
    nd_7413680.edge.oid = "Edge7504080"
    nd_7412912 = nd_7412624.new_child(label="Node7412912", taxon=None, edge_length=0.015326, oid="Node7412912")
    nd_7412912.edge.oid = "Edge7413232"
    nd_7413424 = nd_7412624.new_child(label="Node7413424", taxon=None, edge_length=0.005897, oid="Node7413424")
    nd_7413424.edge.oid = "Edge7413392"
    nd_7412976 = nd_7412912.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.068438, oid="Node7412976")
    nd_7412976.edge.oid = "Edge7413296"
    nd_7413360 = nd_7412912.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.066815, oid="Node7413360")
    nd_7413360.edge.oid = "Edge7413328"
    nd_7413456 = nd_7413424.new_child(label="Node7413456", taxon=None, edge_length=0.010508, oid="Node7413456")
    nd_7413456.edge.oid = "Edge7413200"
    nd_7413712 = nd_7413424.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.083409, oid="Node7413712")
    nd_7413712.edge.oid = "Edge7413264"
    nd_7413488 = nd_7413456.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.057007, oid="Node7413488")
    nd_7413488.edge.oid = "Edge7413552"
    nd_7413648 = nd_7413456.new_child(label="Node7413648", taxon=None, edge_length=0.013392, oid="Node7413648")
    nd_7413648.edge.oid = "Edge7413584"
    nd_7413520 = nd_7413648.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.062766, oid="Node7413520")
    nd_7413520.edge.oid = "Edge7413136"
    nd_7413744 = nd_7413648.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.079145, oid="Node7413744")
    nd_7413744.edge.oid = "Edge7413616"
    tree_7504016 = dendropy.Tree(label="Tree07", taxon_set=tree_list.taxon_set, oid="Tree7504016")
    tree_list.append(tree_7504016, reindex_taxa=False)
    tree_7504016.seed_node.oid = 'Node7504176'
    nd_7504240 = tree_7504016.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.062297, oid="Node7504240")
    nd_7504240.edge.oid = "Edge7504272"
    nd_7504304 = tree_7504016.seed_node.new_child(label="Node7504304", taxon=None, edge_length=0.009899, oid="Node7504304")
    nd_7504304.edge.oid = "Edge7504144"
    nd_5593232 = tree_7504016.seed_node.new_child(label="Node5593232", taxon=None, edge_length=0.003212, oid="Node5593232")
    nd_5593232.edge.oid = "Edge7504336"
    nd_7503920 = nd_7504304.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.069165, oid="Node7503920")
    nd_7503920.edge.oid = "Edge7504400"
    nd_7504496 = nd_7504304.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.073495, oid="Node7504496")
    nd_7504496.edge.oid = "Edge7504432"
    nd_7504528 = nd_5593232.new_child(label="Node7504528", taxon=None, edge_length=0.003672, oid="Node7504528")
    nd_7504528.edge.oid = "Edge7504464"
    nd_7505488 = nd_5593232.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.080896, oid="Node7505488")
    nd_7505488.edge.oid = "Edge7412304"
    nd_7504592 = nd_7504528.new_child(label="Node7504592", taxon=None, edge_length=0.010249, oid="Node7504592")
    nd_7504592.edge.oid = "Edge7504656"
    nd_7505296 = nd_7504528.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.114381, oid="Node7505296")
    nd_7505296.edge.oid = "Edge7505168"
    nd_7504368 = nd_7504592.new_child(label="Node7504368", taxon=None, edge_length=0.00713, oid="Node7504368")
    nd_7504368.edge.oid = "Edge7504720"
    nd_7505040 = nd_7504592.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.079534, oid="Node7505040")
    nd_7505040.edge.oid = "Edge7505520"
    nd_7504624 = nd_7504368.new_child(label="Node7504624", taxon=None, edge_length=0.000386, oid="Node7504624")
    nd_7504624.edge.oid = "Edge7504784"
    nd_7505552 = nd_7504368.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.070614, oid="Node7505552")
    nd_7505552.edge.oid = "Edge7505424"
    nd_7504688 = nd_7504624.new_child(label="Node7504688", taxon=None, edge_length=0.003603, oid="Node7504688")
    nd_7504688.edge.oid = "Edge7504848"
    nd_7505136 = nd_7504624.new_child(label="Node7505136", taxon=None, edge_length=0.017571, oid="Node7505136")
    nd_7505136.edge.oid = "Edge7505232"
    nd_7504752 = nd_7504688.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.061022, oid="Node7504752")
    nd_7504752.edge.oid = "Edge7504912"
    nd_7505008 = nd_7504688.new_child(label="Node7505008", taxon=None, edge_length=0.004698, oid="Node7505008")
    nd_7505008.edge.oid = "Edge7504944"
    nd_7504880 = nd_7505008.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.071291, oid="Node7504880")
    nd_7504880.edge.oid = "Edge7504816"
    nd_7505104 = nd_7505008.new_child(label="Node7505104", taxon=None, edge_length=0.017245, oid="Node7505104")
    nd_7505104.edge.oid = "Edge7505072"
    nd_7504976 = nd_7505104.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.098776, oid="Node7504976")
    nd_7504976.edge.oid = "Edge7504560"
    nd_7505264 = nd_7505104.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.075069, oid="Node7505264")
    nd_7505264.edge.oid = "Edge7505200"
    nd_7505328 = nd_7505136.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.075494, oid="Node7505328")
    nd_7505328.edge.oid = "Edge7505360"
    nd_7505456 = nd_7505136.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.063757, oid="Node7505456")
    nd_7505456.edge.oid = "Edge7505392"
    tree_7505648 = dendropy.Tree(label="Tree08", taxon_set=tree_list.taxon_set, oid="Tree7505648")
    tree_list.append(tree_7505648, reindex_taxa=False)
    tree_7505648.seed_node.oid = 'Node7505680'
    nd_7505776 = tree_7505648.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.05678, oid="Node7505776")
    nd_7505776.edge.oid = "Edge7505808"
    nd_7505840 = tree_7505648.seed_node.new_child(label="Node7505840", taxon=None, edge_length=0.014479, oid="Node7505840")
    nd_7505840.edge.oid = "Edge7503984"
    nd_7507184 = tree_7505648.seed_node.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.060374, oid="Node7507184")
    nd_7507184.edge.oid = "Edge7507024"
    nd_7505584 = nd_7505840.new_child(label="Node7505584", taxon=None, edge_length=0.00182, oid="Node7505584")
    nd_7505584.edge.oid = "Edge7505936"
    nd_7507216 = nd_7505840.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.079156, oid="Node7507216")
    nd_7507216.edge.oid = "Edge7507120"
    nd_7412272 = nd_7505584.new_child(label="Node7412272", taxon=None, edge_length=0.002833, oid="Node7412272")
    nd_7412272.edge.oid = "Edge7506000"
    nd_7506768 = nd_7505584.new_child(label="Node7506768", taxon=None, edge_length=0.010734, oid="Node7506768")
    nd_7506768.edge.oid = "Edge7506896"
    nd_7505872 = nd_7412272.new_child(label="Node7505872", taxon=None, edge_length=0.012314, oid="Node7505872")
    nd_7505872.edge.oid = "Edge7506064"
    nd_7506288 = nd_7412272.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.081982, oid="Node7506288")
    nd_7506288.edge.oid = "Edge7506576"
    nd_7505968 = nd_7505872.new_child(label="Node7505968", taxon=None, edge_length=0.001866, oid="Node7505968")
    nd_7505968.edge.oid = "Edge7506128"
    nd_7506640 = nd_7505872.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.06129, oid="Node7506640")
    nd_7506640.edge.oid = "Edge7506800"
    nd_7506032 = nd_7505968.new_child(label="Node7506032", taxon=None, edge_length=0.005104, oid="Node7506032")
    nd_7506032.edge.oid = "Edge7506192"
    nd_7506832 = nd_7505968.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.07646, oid="Node7506832")
    nd_7506832.edge.oid = "Edge7506704"
    nd_7506096 = nd_7506032.new_child(label="Node7506096", taxon=None, edge_length=0.004377, oid="Node7506096")
    nd_7506096.edge.oid = "Edge7506256"
    nd_7506352 = nd_7506032.new_child(label="Node7506352", taxon=None, edge_length=0.016433, oid="Node7506352")
    nd_7506352.edge.oid = "Edge7505904"
    nd_7506160 = nd_7506096.new_child(label="Node7506160", taxon=None, edge_length=0.013813, oid="Node7506160")
    nd_7506160.edge.oid = "Edge7506320"
    nd_7506512 = nd_7506096.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.083335, oid="Node7506512")
    nd_7506512.edge.oid = "Edge7506448"
    nd_7506224 = nd_7506160.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.065928, oid="Node7506224")
    nd_7506224.edge.oid = "Edge7506384"
    nd_7506480 = nd_7506160.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.10077, oid="Node7506480")
    nd_7506480.edge.oid = "Edge7506416"
    nd_7506544 = nd_7506352.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.073392, oid="Node7506544")
    nd_7506544.edge.oid = "Edge7506608"
    nd_7506736 = nd_7506352.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.065748, oid="Node7506736")
    nd_7506736.edge.oid = "Edge7506672"
    nd_7506992 = nd_7506768.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.064242, oid="Node7506992")
    nd_7506992.edge.oid = "Edge7506928"
    nd_7507088 = nd_7506768.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.109434, oid="Node7507088")
    nd_7507088.edge.oid = "Edge7507056"
    tree_7506960 = dendropy.Tree(label="Tree09", taxon_set=tree_list.taxon_set, oid="Tree7506960")
    tree_list.append(tree_7506960, reindex_taxa=False)
    tree_7506960.seed_node.oid = 'Node7507376'
    nd_7507440 = tree_7506960.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.063004, oid="Node7507440")
    nd_7507440.edge.oid = "Edge7507472"
    nd_7507504 = tree_7506960.seed_node.new_child(label="Node7507504", taxon=None, edge_length=0.018059, oid="Node7507504")
    nd_7507504.edge.oid = "Edge7504112"
    nd_7507728 = tree_7506960.seed_node.new_child(label="Node7507728", taxon=None, edge_length=0.00238, oid="Node7507728")
    nd_7507728.edge.oid = "Edge7507664"
    nd_7506864 = nd_7507504.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.061413, oid="Node7506864")
    nd_7506864.edge.oid = "Edge7507600"
    nd_7507696 = nd_7507504.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.099507, oid="Node7507696")
    nd_7507696.edge.oid = "Edge7507632"
    nd_7507536 = nd_7507728.new_child(label="Node7507536", taxon=None, edge_length=0.003714, oid="Node7507536")
    nd_7507536.edge.oid = "Edge7507760"
    nd_7561520 = nd_7507728.new_child(label="Node7561520", taxon=None, edge_length=0.010478, oid="Node7561520")
    nd_7561520.edge.oid = "Edge7561936"
    nd_7507792 = nd_7507536.new_child(label="Node7507792", taxon=None, edge_length=0.016533, oid="Node7507792")
    nd_7507792.edge.oid = "Edge7507856"
    nd_7561904 = nd_7507536.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.068818, oid="Node7561904")
    nd_7561904.edge.oid = "Edge7561776"
    nd_7507568 = nd_7507792.new_child(label="Node7507568", taxon=None, edge_length=0.010763, oid="Node7507568")
    nd_7507568.edge.oid = "Edge7507920"
    nd_7561712 = nd_7507792.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.070762, oid="Node7561712")
    nd_7561712.edge.oid = "Edge7561840"
    nd_7507824 = nd_7507568.new_child(label="Node7507824", taxon=None, edge_length=0.002182, oid="Node7507824")
    nd_7507824.edge.oid = "Edge7561264"
    nd_7561296 = nd_7507568.new_child(label="Node7561296", taxon=None, edge_length=0.011197, oid="Node7561296")
    nd_7561296.edge.oid = "Edge7561488"
    nd_7507888 = nd_7507824.new_child(label="Node7507888", taxon=None, edge_length=0.014459, oid="Node7507888")
    nd_7507888.edge.oid = "Edge7561328"
    nd_7507344 = nd_7507824.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.054698, oid="Node7507344")
    nd_7507344.edge.oid = "Edge7561360"
    nd_7507952 = nd_7507888.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.075877, oid="Node7507952")
    nd_7507952.edge.oid = "Edge7561392"
    nd_7561456 = nd_7507888.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.071375, oid="Node7561456")
    nd_7561456.edge.oid = "Edge7561424"
    nd_7561552 = nd_7561296.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.062291, oid="Node7561552")
    nd_7561552.edge.oid = "Edge7561616"
    nd_7561744 = nd_7561296.new_child(label="Node7561744", taxon=None, edge_length=0.00521, oid="Node7561744")
    nd_7561744.edge.oid = "Edge7561680"
    nd_7561584 = nd_7561744.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.074463, oid="Node7561584")
    nd_7561584.edge.oid = "Edge7561648"
    nd_7561872 = nd_7561744.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.103887, oid="Node7561872")
    nd_7561872.edge.oid = "Edge7561808"
    nd_7562064 = nd_7561520.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.071307, oid="Node7562064")
    nd_7562064.edge.oid = "Edge7562000"
    nd_7562160 = nd_7561520.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.073392, oid="Node7562160")
    nd_7562160.edge.oid = "Edge7562128"
    tree_7562192 = dendropy.Tree(label="Tree10", taxon_set=tree_list.taxon_set, oid="Tree7562192")
    tree_list.append(tree_7562192, reindex_taxa=False)
    tree_7562192.seed_node.oid = 'Node7562320'
    nd_7562384 = tree_7562192.seed_node.new_child(label="Aspidites ramsayi", taxon=tax_5593072, edge_length=0.053901, oid="Node7562384")
    nd_7562384.edge.oid = "Edge7562416"
    nd_7562448 = tree_7562192.seed_node.new_child(label="Node7562448", taxon=None, edge_length=0.013947, oid="Node7562448")
    nd_7562448.edge.oid = "Edge7562096"
    nd_7563344 = tree_7562192.seed_node.new_child(label="Python brongersmai", taxon=tax_5595024, edge_length=0.107868, oid="Node7563344")
    nd_7563344.edge.oid = "Edge7563568"
    nd_7562032 = nd_7562448.new_child(label="Node7562032", taxon=None, edge_length=0.002034, oid="Node7562032")
    nd_7562032.edge.oid = "Edge7562544"
    nd_7563536 = nd_7562448.new_child(label="Antaresia stimsoni", taxon=tax_5594480, edge_length=0.064652, oid="Node7563536")
    nd_7563536.edge.oid = "Edge7563312"
    nd_7355728 = nd_7562032.new_child(label="Node7355728", taxon=None, edge_length=0.00777, oid="Node7355728")
    nd_7355728.edge.oid = "Edge7562608"
    nd_7563184 = nd_7562032.new_child(label="Node7563184", taxon=None, edge_length=0.00373, oid="Node7563184")
    nd_7563184.edge.oid = "Edge7562512"
    nd_7562480 = nd_7355728.new_child(label="Node7562480", taxon=None, edge_length=0.003041, oid="Node7562480")
    nd_7562480.edge.oid = "Edge7562672"
    nd_7563056 = nd_7355728.new_child(label="Node7563056", taxon=None, edge_length=0.008083, oid="Node7563056")
    nd_7563056.edge.oid = "Edge7562928"
    nd_7562576 = nd_7562480.new_child(label="Bothrochilus boa", taxon=tax_5594064, edge_length=0.07927, oid="Node7562576")
    nd_7562576.edge.oid = "Edge7562736"
    nd_7562832 = nd_7562480.new_child(label="Node7562832", taxon=None, edge_length=0.016769, oid="Node7562832")
    nd_7562832.edge.oid = "Edge7562768"
    nd_7562704 = nd_7562832.new_child(label="Morelia bredli", taxon=tax_5594992, edge_length=0.053956, oid="Node7562704")
    nd_7562704.edge.oid = "Edge7562640"
    nd_7562960 = nd_7562832.new_child(label="Python timoriensis", taxon=tax_5594768, edge_length=0.098874, oid="Node7562960")
    nd_7562960.edge.oid = "Edge7562896"
    nd_7562864 = nd_7563056.new_child(label="Node7562864", taxon=None, edge_length=0.010712, oid="Node7562864")
    nd_7562864.edge.oid = "Edge7562992"
    nd_7503952 = nd_7563056.new_child(label="Morelia carinata", taxon=tax_7352592, edge_length=0.075934, oid="Node7503952")
    nd_7503952.edge.oid = "Edge7563088"
    nd_7562800 = nd_7562864.new_child(label="Liasis fuscus", taxon=tax_5593680, edge_length=0.075223, oid="Node7562800")
    nd_7562800.edge.oid = "Edge7507152"
    nd_7505744 = nd_7562864.new_child(label="Morelia viridis", taxon=tax_5594544, edge_length=0.065805, oid="Node7505744")
    nd_7505744.edge.oid = "Edge7507280"
    nd_7563120 = nd_7563184.new_child(label="Antaresia perthensis", taxon=tax_5595120, edge_length=0.082425, oid="Node7563120")
    nd_7563120.edge.oid = "Edge7561968"
    nd_7563280 = nd_7563184.new_child(label="Node7563280", taxon=None, edge_length=0.008508, oid="Node7563280")
    nd_7563280.edge.oid = "Edge7563216"
    nd_7563024 = nd_7563280.new_child(label="Node7563024", taxon=None, edge_length=0.010583, oid="Node7563024")
    nd_7563024.edge.oid = "Edge7562256"
    nd_7563504 = nd_7563280.new_child(label="Morelia oenpelliensis", taxon=tax_5594832, edge_length=0.070289, oid="Node7563504")
    nd_7563504.edge.oid = "Edge7563440"
    nd_7563152 = nd_7563024.new_child(label="Antaresia maculosa", taxon=tax_5594608, edge_length=0.075429, oid="Node7563152")
    nd_7563152.edge.oid = "Edge7563376"
    nd_7563472 = nd_7563024.new_child(label="Morelia boeleni", taxon=tax_5594448, edge_length=0.0787, oid="Node7563472")
    nd_7563472.edge.oid = "Edge7563408"

    # set labels of nodes with taxa to taxon label, else oid (for consistent
    # identification in debugging)
    for t in tree_list:
        t.assign_node_labels_from_taxon_or_oid()

    return tree_list

def reference_tree_list_postorder_node_labels():
    return [
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Morelia boeleni','Node5594256','Node5593936','Morelia viridis','Node5593872','Antaresia maculosa','Node5593808','Antaresia stimsoni','Node5593744','Python timoriensis','Node5593328','Morelia oenpelliensis','Node5593616','Morelia bredli','Node5593648','Antaresia perthensis','Python brongersmai','Node5594896','Node5593584','Morelia carinata','Node5593456'],
        ['Aspidites ramsayi','Bothrochilus boa','Morelia boeleni','Morelia oenpelliensis','Node7353200','Node5593168','Morelia viridis','Node7352688','Antaresia perthensis','Node7352880','Liasis fuscus','Antaresia stimsoni','Python timoriensis','Morelia carinata','Node7353744','Python brongersmai','Node7353712','Antaresia maculosa','Node7353840','Node7353680','Node7353360','Morelia bredli','Node7353232','Node7352752'],
        ['Aspidites ramsayi','Bothrochilus boa','Antaresia perthensis','Node7354544','Python brongersmai','Node7354096','Python timoriensis','Antaresia maculosa','Node7354928','Morelia oenpelliensis','Node7354736','Node7354448','Liasis fuscus','Morelia viridis','Morelia bredli','Node7355408','Node7355344','Node7354256','Antaresia stimsoni','Node7354416','Morelia carinata','Morelia boeleni','Node7355600','Node7354288'],
        ['Aspidites ramsayi','Bothrochilus boa','Antaresia perthensis','Node7356112','Morelia bredli','Node7355920','Liasis fuscus','Morelia oenpelliensis','Node7409872','Morelia viridis','Node7409840','Python timoriensis','Morelia carinata','Node7410320','Antaresia maculosa','Node7410064','Node7356208','Python brongersmai','Node7409776','Antaresia stimsoni','Node7356144','Node7356080','Morelia boeleni','Node7355952'],
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Morelia carinata','Node7411216','Node7411024','Antaresia stimsoni','Antaresia perthensis','Node7411344','Antaresia maculosa','Node5593424','Morelia viridis','Morelia bredli','Morelia oenpelliensis','Node7411824','Node7411600','Python timoriensis','Node7411504','Node7411248','Python brongersmai','Morelia boeleni','Node7412144','Node7411440','Node7410896'],
        ['Aspidites ramsayi','Bothrochilus boa','Python brongersmai','Node7412528','Python timoriensis','Node7412016','Morelia bredli','Antaresia maculosa','Node7412688','Node7412496','Liasis fuscus','Morelia oenpelliensis','Node7412912','Antaresia stimsoni','Morelia carinata','Morelia boeleni','Node7413648','Node7413456','Antaresia perthensis','Node7413424','Node7412624','Morelia viridis','Node7413168','Node7412368'],
        ['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Node7504304','Antaresia stimsoni','Morelia viridis','Python timoriensis','Morelia boeleni','Node7505104','Node7505008','Node7504688','Antaresia maculosa','Morelia oenpelliensis','Node7505136','Node7504624','Morelia bredli','Node7504368','Morelia carinata','Node7504592','Python brongersmai','Node7504528','Antaresia perthensis','Node5593232','Node7504176'],
        ['Aspidites ramsayi','Bothrochilus boa','Python timoriensis','Node7506160','Antaresia perthensis','Node7506096','Antaresia maculosa','Morelia oenpelliensis','Node7506352','Node7506032','Morelia carinata','Node7505968','Antaresia stimsoni','Node7505872','Morelia boeleni','Node7412272','Morelia viridis','Python brongersmai','Node7506768','Node7505584','Liasis fuscus','Node7505840','Morelia bredli','Node7505680'],
        ['Aspidites ramsayi','Bothrochilus boa','Python brongersmai','Node7507504','Liasis fuscus','Antaresia maculosa','Node7507888','Antaresia stimsoni','Node7507824','Morelia bredli','Antaresia perthensis','Python timoriensis','Node7561744','Node7561296','Node7507568','Morelia carinata','Node7507792','Morelia oenpelliensis','Node7507536','Morelia viridis','Morelia boeleni','Node7561520','Node7507728','Node7507376'],
        ['Aspidites ramsayi','Bothrochilus boa','Morelia bredli','Python timoriensis','Node7562832','Node7562480','Liasis fuscus','Morelia viridis','Node7562864','Morelia carinata','Node7563056','Node7355728','Antaresia perthensis','Antaresia maculosa','Morelia boeleni','Node7563024','Morelia oenpelliensis','Node7563280','Node7563184','Node7562032','Antaresia stimsoni','Node7562448','Python brongersmai','Node7562320'],
    ]

def reference_tree_list_newick_string():
    return """\
        ('Aspidites ramsayi':0.056823,(((((((('Bothrochilus boa':0.073578,('Liasis fuscus':0.074347,'Morelia boeleni':0.07628)Node5594256:0.003063)Node5593936:0.011155,'Morelia viridis':0.066566)Node5593872:0.010103,'Antaresia maculosa':0.080749)Node5593808:0.007681,'Antaresia stimsoni':0.061996)Node5593744:0.0,'Python timoriensis':0.113605)Node5593328:0.002264,'Morelia oenpelliensis':0.074047)Node5593616:0.007655,'Morelia bredli':0.065046)Node5593648:0.003495,('Antaresia perthensis':0.065004,'Python brongersmai':0.107706)Node5594896:0.018613)Node5593584:0.011159,'Morelia carinata':0.079321)Node5593456;
        ('Aspidites ramsayi':0.065297,((('Bothrochilus boa':0.073565,('Morelia boeleni':0.069885,'Morelia oenpelliensis':0.06247)Node7353200:0.007879)Node5593168:0.016022,'Morelia viridis':0.066015)Node7352688:0.010946,'Antaresia perthensis':0.079424)Node7352880:0.001444,(('Liasis fuscus':0.082049,('Antaresia stimsoni':0.056894,((('Python timoriensis':0.092077,'Morelia carinata':0.0579)Node7353744:0.015909,'Python brongersmai':0.114385)Node7353712:0.010752,'Antaresia maculosa':0.079792)Node7353840:0.004552)Node7353680:0.006907)Node7353360:0.008007,'Morelia bredli':0.063725)Node7353232:0.004077)Node7352752;
        ('Aspidites ramsayi':0.06111,((((('Bothrochilus boa':0.067908,'Antaresia perthensis':0.077196)Node7354544:0.008688,'Python brongersmai':0.1069)Node7354096:0.010367,(('Python timoriensis':0.109318,'Antaresia maculosa':0.081382)Node7354928:0.004129,'Morelia oenpelliensis':0.068424)Node7354736:0.007731)Node7354448:0.005167,('Liasis fuscus':0.079478,('Morelia viridis':0.062239,'Morelia bredli':0.060484)Node7355408:0.013155)Node7355344:0.00574)Node7354256:0.006988,'Antaresia stimsoni':0.063808)Node7354416:0.005408,('Morelia carinata':0.069501,'Morelia boeleni':0.073602)Node7355600:0.012121)Node7354288;
        ('Aspidites ramsayi':0.053028,((('Bothrochilus boa':0.072831,'Antaresia perthensis':0.073004)Node7356112:0.012413,'Morelia bredli':0.069933)Node7355920:0.002567,((((('Liasis fuscus':0.072339,'Morelia oenpelliensis':0.062058)Node7409872:0.011723,'Morelia viridis':0.065297)Node7409840:0.011519,(('Python timoriensis':0.093762,'Morelia carinata':0.055992)Node7410320:0.022912,'Antaresia maculosa':0.079945)Node7410064:0.007647)Node7356208:5.134e-05,'Python brongersmai':0.118976)Node7409776:0.004803,'Antaresia stimsoni':0.064721)Node7356144:0.005969)Node7356080:0.014376,'Morelia boeleni':0.076894)Node7355952;
        ('Aspidites ramsayi':0.06391,('Bothrochilus boa':0.071835,('Liasis fuscus':0.075788,'Morelia carinata':0.07116)Node7411216:0.010888)Node7411024:0.00451,(((('Antaresia stimsoni':0.048754,'Antaresia perthensis':0.069682)Node7411344:0.012356,'Antaresia maculosa':0.077058)Node5593424:0.013817,(('Morelia viridis':0.07026,('Morelia bredli':0.066182,'Morelia oenpelliensis':0.068166)Node7411824:0.001443)Node7411600:0.003147,'Python timoriensis':0.108243)Node7411504:0.012897)Node7411248:0.011944,('Python brongersmai':0.108072,'Morelia boeleni':0.073251)Node7412144:0.00612)Node7411440:0.005404)Node7410896;
        ('Aspidites ramsayi':0.06501,((('Bothrochilus boa':0.059984,'Python brongersmai':0.100555)Node7412528:0.017824,'Python timoriensis':0.102639)Node7412016:0.011722,('Morelia bredli':0.0538,'Antaresia maculosa':0.07843)Node7412688:0.014836)Node7412496:0.003644,((('Liasis fuscus':0.068438,'Morelia oenpelliensis':0.066815)Node7412912:0.015326,(('Antaresia stimsoni':0.057007,('Morelia carinata':0.062766,'Morelia boeleni':0.079145)Node7413648:0.013392)Node7413456:0.010508,'Antaresia perthensis':0.083409)Node7413424:0.005897)Node7412624:0.0,'Morelia viridis':0.072091)Node7413168:0.008523)Node7412368;
        ('Aspidites ramsayi':0.062297,('Bothrochilus boa':0.069165,'Liasis fuscus':0.073495)Node7504304:0.009899,(((((('Antaresia stimsoni':0.061022,('Morelia viridis':0.071291,('Python timoriensis':0.098776,'Morelia boeleni':0.075069)Node7505104:0.017245)Node7505008:0.004698)Node7504688:0.003603,('Antaresia maculosa':0.075494,'Morelia oenpelliensis':0.063757)Node7505136:0.017571)Node7504624:0.000386,'Morelia bredli':0.070614)Node7504368:0.00713,'Morelia carinata':0.079534)Node7504592:0.010249,'Python brongersmai':0.114381)Node7504528:0.003672,'Antaresia perthensis':0.080896)Node5593232:0.003212)Node7504176;
        ('Aspidites ramsayi':0.05678,(((((((('Bothrochilus boa':0.065928,'Python timoriensis':0.10077)Node7506160:0.013813,'Antaresia perthensis':0.083335)Node7506096:0.004377,('Antaresia maculosa':0.073392,'Morelia oenpelliensis':0.065748)Node7506352:0.016433)Node7506032:0.005104,'Morelia carinata':0.07646)Node7505968:0.001866,'Antaresia stimsoni':0.06129)Node7505872:0.012314,'Morelia boeleni':0.081982)Node7412272:0.002833,('Morelia viridis':0.064242,'Python brongersmai':0.109434)Node7506768:0.010734)Node7505584:0.00182,'Liasis fuscus':0.079156)Node7505840:0.014479,'Morelia bredli':0.060374)Node7505680;
        ('Aspidites ramsayi':0.063004,('Bothrochilus boa':0.061413,'Python brongersmai':0.099507)Node7507504:0.018059,(((((('Liasis fuscus':0.075877,'Antaresia maculosa':0.071375)Node7507888:0.014459,'Antaresia stimsoni':0.054698)Node7507824:0.002182,('Morelia bredli':0.062291,('Antaresia perthensis':0.074463,'Python timoriensis':0.103887)Node7561744:0.00521)Node7561296:0.011197)Node7507568:0.010763,'Morelia carinata':0.070762)Node7507792:0.016533,'Morelia oenpelliensis':0.068818)Node7507536:0.003714,('Morelia viridis':0.071307,'Morelia boeleni':0.073392)Node7561520:0.010478)Node7507728:0.00238)Node7507376;
        ('Aspidites ramsayi':0.053901,(((('Bothrochilus boa':0.07927,('Morelia bredli':0.053956,'Python timoriensis':0.098874)Node7562832:0.016769)Node7562480:0.003041,(('Liasis fuscus':0.075223,'Morelia viridis':0.065805)Node7562864:0.010712,'Morelia carinata':0.075934)Node7563056:0.008083)Node7355728:0.00777,('Antaresia perthensis':0.082425,(('Antaresia maculosa':0.075429,'Morelia boeleni':0.0787)Node7563024:0.010583,'Morelia oenpelliensis':0.070289)Node7563280:0.008508)Node7563184:0.00373)Node7562032:0.002034,'Antaresia stimsoni':0.064652)Node7562448:0.013947,'Python brongersmai':0.107868)Node7562320;
    """

def reference_tree_list_node_relationships():
    treelist_node_references = [
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node5593456', child_labels=[], edge_length=0.056823, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node5593936', child_labels=[], edge_length=0.073578, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node5594256', child_labels=[], edge_length=0.074347, taxon_label='Liasis fuscus'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node5594256', child_labels=[], edge_length=0.07628, taxon_label='Morelia boeleni'),
            'Node5594256' : NodeRelationship(parent_label='Node5593936', child_labels=['Liasis fuscus','Morelia boeleni'], edge_length=0.003063, taxon_label=None),
            'Node5593936' : NodeRelationship(parent_label='Node5593872', child_labels=['Bothrochilus boa','Node5594256'], edge_length=0.011155, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node5593872', child_labels=[], edge_length=0.066566, taxon_label='Morelia viridis'),
            'Node5593872' : NodeRelationship(parent_label='Node5593808', child_labels=['Node5593936','Morelia viridis'], edge_length=0.010103, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5593808', child_labels=[], edge_length=0.080749, taxon_label='Antaresia maculosa'),
            'Node5593808' : NodeRelationship(parent_label='Node5593744', child_labels=['Node5593872','Antaresia maculosa'], edge_length=0.007681, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node5593744', child_labels=[], edge_length=0.061996, taxon_label='Antaresia stimsoni'),
            'Node5593744' : NodeRelationship(parent_label='Node5593328', child_labels=['Node5593808','Antaresia stimsoni'], edge_length=0.0, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node5593328', child_labels=[], edge_length=0.113605, taxon_label='Python timoriensis'),
            'Node5593328' : NodeRelationship(parent_label='Node5593616', child_labels=['Node5593744','Python timoriensis'], edge_length=0.002264, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node5593616', child_labels=[], edge_length=0.074047, taxon_label='Morelia oenpelliensis'),
            'Node5593616' : NodeRelationship(parent_label='Node5593648', child_labels=['Node5593328','Morelia oenpelliensis'], edge_length=0.007655, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node5593648', child_labels=[], edge_length=0.065046, taxon_label='Morelia bredli'),
            'Node5593648' : NodeRelationship(parent_label='Node5593584', child_labels=['Node5593616','Morelia bredli'], edge_length=0.003495, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node5594896', child_labels=[], edge_length=0.065004, taxon_label='Antaresia perthensis'),
            'Python brongersmai' : NodeRelationship(parent_label='Node5594896', child_labels=[], edge_length=0.107706, taxon_label='Python brongersmai'),
            'Node5594896' : NodeRelationship(parent_label='Node5593584', child_labels=['Antaresia perthensis','Python brongersmai'], edge_length=0.018613, taxon_label=None),
            'Node5593584' : NodeRelationship(parent_label='Node5593456', child_labels=['Node5593648','Node5594896'], edge_length=0.011159, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node5593456', child_labels=[], edge_length=0.079321, taxon_label='Morelia carinata'),
            'Node5593456' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node5593584','Morelia carinata'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7352752', child_labels=[], edge_length=0.065297, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node5593168', child_labels=[], edge_length=0.073565, taxon_label='Bothrochilus boa'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7353200', child_labels=[], edge_length=0.069885, taxon_label='Morelia boeleni'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7353200', child_labels=[], edge_length=0.06247, taxon_label='Morelia oenpelliensis'),
            'Node7353200' : NodeRelationship(parent_label='Node5593168', child_labels=['Morelia boeleni','Morelia oenpelliensis'], edge_length=0.007879, taxon_label=None),
            'Node5593168' : NodeRelationship(parent_label='Node7352688', child_labels=['Bothrochilus boa','Node7353200'], edge_length=0.016022, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7352688', child_labels=[], edge_length=0.066015, taxon_label='Morelia viridis'),
            'Node7352688' : NodeRelationship(parent_label='Node7352880', child_labels=['Node5593168','Morelia viridis'], edge_length=0.010946, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7352880', child_labels=[], edge_length=0.079424, taxon_label='Antaresia perthensis'),
            'Node7352880' : NodeRelationship(parent_label='Node7352752', child_labels=['Node7352688','Antaresia perthensis'], edge_length=0.001444, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7353360', child_labels=[], edge_length=0.082049, taxon_label='Liasis fuscus'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7353680', child_labels=[], edge_length=0.056894, taxon_label='Antaresia stimsoni'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7353744', child_labels=[], edge_length=0.092077, taxon_label='Python timoriensis'),
            'Morelia carinata' : NodeRelationship(parent_label='Node7353744', child_labels=[], edge_length=0.0579, taxon_label='Morelia carinata'),
            'Node7353744' : NodeRelationship(parent_label='Node7353712', child_labels=['Python timoriensis','Morelia carinata'], edge_length=0.015909, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7353712', child_labels=[], edge_length=0.114385, taxon_label='Python brongersmai'),
            'Node7353712' : NodeRelationship(parent_label='Node7353840', child_labels=['Node7353744','Python brongersmai'], edge_length=0.010752, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7353840', child_labels=[], edge_length=0.079792, taxon_label='Antaresia maculosa'),
            'Node7353840' : NodeRelationship(parent_label='Node7353680', child_labels=['Node7353712','Antaresia maculosa'], edge_length=0.004552, taxon_label=None),
            'Node7353680' : NodeRelationship(parent_label='Node7353360', child_labels=['Antaresia stimsoni','Node7353840'], edge_length=0.006907, taxon_label=None),
            'Node7353360' : NodeRelationship(parent_label='Node7353232', child_labels=['Liasis fuscus','Node7353680'], edge_length=0.008007, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7353232', child_labels=[], edge_length=0.063725, taxon_label='Morelia bredli'),
            'Node7353232' : NodeRelationship(parent_label='Node7352752', child_labels=['Node7353360','Morelia bredli'], edge_length=0.004077, taxon_label=None),
            'Node7352752' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7352880','Node7353232'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7354288', child_labels=[], edge_length=0.06111, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7354544', child_labels=[], edge_length=0.067908, taxon_label='Bothrochilus boa'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7354544', child_labels=[], edge_length=0.077196, taxon_label='Antaresia perthensis'),
            'Node7354544' : NodeRelationship(parent_label='Node7354096', child_labels=['Bothrochilus boa','Antaresia perthensis'], edge_length=0.008688, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7354096', child_labels=[], edge_length=0.1069, taxon_label='Python brongersmai'),
            'Node7354096' : NodeRelationship(parent_label='Node7354448', child_labels=['Node7354544','Python brongersmai'], edge_length=0.010367, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node7354928', child_labels=[], edge_length=0.109318, taxon_label='Python timoriensis'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7354928', child_labels=[], edge_length=0.081382, taxon_label='Antaresia maculosa'),
            'Node7354928' : NodeRelationship(parent_label='Node7354736', child_labels=['Python timoriensis','Antaresia maculosa'], edge_length=0.004129, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7354736', child_labels=[], edge_length=0.068424, taxon_label='Morelia oenpelliensis'),
            'Node7354736' : NodeRelationship(parent_label='Node7354448', child_labels=['Node7354928','Morelia oenpelliensis'], edge_length=0.007731, taxon_label=None),
            'Node7354448' : NodeRelationship(parent_label='Node7354256', child_labels=['Node7354096','Node7354736'], edge_length=0.005167, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7355344', child_labels=[], edge_length=0.079478, taxon_label='Liasis fuscus'),
            'Morelia viridis' : NodeRelationship(parent_label='Node7355408', child_labels=[], edge_length=0.062239, taxon_label='Morelia viridis'),
            'Morelia bredli' : NodeRelationship(parent_label='Node7355408', child_labels=[], edge_length=0.060484, taxon_label='Morelia bredli'),
            'Node7355408' : NodeRelationship(parent_label='Node7355344', child_labels=['Morelia viridis','Morelia bredli'], edge_length=0.013155, taxon_label=None),
            'Node7355344' : NodeRelationship(parent_label='Node7354256', child_labels=['Liasis fuscus','Node7355408'], edge_length=0.00574, taxon_label=None),
            'Node7354256' : NodeRelationship(parent_label='Node7354416', child_labels=['Node7354448','Node7355344'], edge_length=0.006988, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7354416', child_labels=[], edge_length=0.063808, taxon_label='Antaresia stimsoni'),
            'Node7354416' : NodeRelationship(parent_label='Node7354288', child_labels=['Node7354256','Antaresia stimsoni'], edge_length=0.005408, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7355600', child_labels=[], edge_length=0.069501, taxon_label='Morelia carinata'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7355600', child_labels=[], edge_length=0.073602, taxon_label='Morelia boeleni'),
            'Node7355600' : NodeRelationship(parent_label='Node7354288', child_labels=['Morelia carinata','Morelia boeleni'], edge_length=0.012121, taxon_label=None),
            'Node7354288' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7354416','Node7355600'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7355952', child_labels=[], edge_length=0.053028, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7356112', child_labels=[], edge_length=0.072831, taxon_label='Bothrochilus boa'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7356112', child_labels=[], edge_length=0.073004, taxon_label='Antaresia perthensis'),
            'Node7356112' : NodeRelationship(parent_label='Node7355920', child_labels=['Bothrochilus boa','Antaresia perthensis'], edge_length=0.012413, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7355920', child_labels=[], edge_length=0.069933, taxon_label='Morelia bredli'),
            'Node7355920' : NodeRelationship(parent_label='Node7356080', child_labels=['Node7356112','Morelia bredli'], edge_length=0.002567, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7409872', child_labels=[], edge_length=0.072339, taxon_label='Liasis fuscus'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7409872', child_labels=[], edge_length=0.062058, taxon_label='Morelia oenpelliensis'),
            'Node7409872' : NodeRelationship(parent_label='Node7409840', child_labels=['Liasis fuscus','Morelia oenpelliensis'], edge_length=0.011723, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7409840', child_labels=[], edge_length=0.065297, taxon_label='Morelia viridis'),
            'Node7409840' : NodeRelationship(parent_label='Node7356208', child_labels=['Node7409872','Morelia viridis'], edge_length=0.011519, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node7410320', child_labels=[], edge_length=0.093762, taxon_label='Python timoriensis'),
            'Morelia carinata' : NodeRelationship(parent_label='Node7410320', child_labels=[], edge_length=0.055992, taxon_label='Morelia carinata'),
            'Node7410320' : NodeRelationship(parent_label='Node7410064', child_labels=['Python timoriensis','Morelia carinata'], edge_length=0.022912, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7410064', child_labels=[], edge_length=0.079945, taxon_label='Antaresia maculosa'),
            'Node7410064' : NodeRelationship(parent_label='Node7356208', child_labels=['Node7410320','Antaresia maculosa'], edge_length=0.007647, taxon_label=None),
            'Node7356208' : NodeRelationship(parent_label='Node7409776', child_labels=['Node7409840','Node7410064'], edge_length=5.134e-05, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7409776', child_labels=[], edge_length=0.118976, taxon_label='Python brongersmai'),
            'Node7409776' : NodeRelationship(parent_label='Node7356144', child_labels=['Node7356208','Python brongersmai'], edge_length=0.004803, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7356144', child_labels=[], edge_length=0.064721, taxon_label='Antaresia stimsoni'),
            'Node7356144' : NodeRelationship(parent_label='Node7356080', child_labels=['Node7409776','Antaresia stimsoni'], edge_length=0.005969, taxon_label=None),
            'Node7356080' : NodeRelationship(parent_label='Node7355952', child_labels=['Node7355920','Node7356144'], edge_length=0.014376, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7355952', child_labels=[], edge_length=0.076894, taxon_label='Morelia boeleni'),
            'Node7355952' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7356080','Morelia boeleni'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7410896', child_labels=[], edge_length=0.06391, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7411024', child_labels=[], edge_length=0.071835, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7411216', child_labels=[], edge_length=0.075788, taxon_label='Liasis fuscus'),
            'Morelia carinata' : NodeRelationship(parent_label='Node7411216', child_labels=[], edge_length=0.07116, taxon_label='Morelia carinata'),
            'Node7411216' : NodeRelationship(parent_label='Node7411024', child_labels=['Liasis fuscus','Morelia carinata'], edge_length=0.010888, taxon_label=None),
            'Node7411024' : NodeRelationship(parent_label='Node7410896', child_labels=['Bothrochilus boa','Node7411216'], edge_length=0.00451, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7411344', child_labels=[], edge_length=0.048754, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7411344', child_labels=[], edge_length=0.069682, taxon_label='Antaresia perthensis'),
            'Node7411344' : NodeRelationship(parent_label='Node5593424', child_labels=['Antaresia stimsoni','Antaresia perthensis'], edge_length=0.012356, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node5593424', child_labels=[], edge_length=0.077058, taxon_label='Antaresia maculosa'),
            'Node5593424' : NodeRelationship(parent_label='Node7411248', child_labels=['Node7411344','Antaresia maculosa'], edge_length=0.013817, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7411600', child_labels=[], edge_length=0.07026, taxon_label='Morelia viridis'),
            'Morelia bredli' : NodeRelationship(parent_label='Node7411824', child_labels=[], edge_length=0.066182, taxon_label='Morelia bredli'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7411824', child_labels=[], edge_length=0.068166, taxon_label='Morelia oenpelliensis'),
            'Node7411824' : NodeRelationship(parent_label='Node7411600', child_labels=['Morelia bredli','Morelia oenpelliensis'], edge_length=0.001443, taxon_label=None),
            'Node7411600' : NodeRelationship(parent_label='Node7411504', child_labels=['Morelia viridis','Node7411824'], edge_length=0.003147, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node7411504', child_labels=[], edge_length=0.108243, taxon_label='Python timoriensis'),
            'Node7411504' : NodeRelationship(parent_label='Node7411248', child_labels=['Node7411600','Python timoriensis'], edge_length=0.012897, taxon_label=None),
            'Node7411248' : NodeRelationship(parent_label='Node7411440', child_labels=['Node5593424','Node7411504'], edge_length=0.011944, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7412144', child_labels=[], edge_length=0.108072, taxon_label='Python brongersmai'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7412144', child_labels=[], edge_length=0.073251, taxon_label='Morelia boeleni'),
            'Node7412144' : NodeRelationship(parent_label='Node7411440', child_labels=['Python brongersmai','Morelia boeleni'], edge_length=0.00612, taxon_label=None),
            'Node7411440' : NodeRelationship(parent_label='Node7410896', child_labels=['Node7411248','Node7412144'], edge_length=0.005404, taxon_label=None),
            'Node7410896' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7411024','Node7411440'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7412368', child_labels=[], edge_length=0.06501, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7412528', child_labels=[], edge_length=0.059984, taxon_label='Bothrochilus boa'),
            'Python brongersmai' : NodeRelationship(parent_label='Node7412528', child_labels=[], edge_length=0.100555, taxon_label='Python brongersmai'),
            'Node7412528' : NodeRelationship(parent_label='Node7412016', child_labels=['Bothrochilus boa','Python brongersmai'], edge_length=0.017824, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node7412016', child_labels=[], edge_length=0.102639, taxon_label='Python timoriensis'),
            'Node7412016' : NodeRelationship(parent_label='Node7412496', child_labels=['Node7412528','Python timoriensis'], edge_length=0.011722, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7412688', child_labels=[], edge_length=0.0538, taxon_label='Morelia bredli'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7412688', child_labels=[], edge_length=0.07843, taxon_label='Antaresia maculosa'),
            'Node7412688' : NodeRelationship(parent_label='Node7412496', child_labels=['Morelia bredli','Antaresia maculosa'], edge_length=0.014836, taxon_label=None),
            'Node7412496' : NodeRelationship(parent_label='Node7412368', child_labels=['Node7412016','Node7412688'], edge_length=0.003644, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7412912', child_labels=[], edge_length=0.068438, taxon_label='Liasis fuscus'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7412912', child_labels=[], edge_length=0.066815, taxon_label='Morelia oenpelliensis'),
            'Node7412912' : NodeRelationship(parent_label='Node7412624', child_labels=['Liasis fuscus','Morelia oenpelliensis'], edge_length=0.015326, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7413456', child_labels=[], edge_length=0.057007, taxon_label='Antaresia stimsoni'),
            'Morelia carinata' : NodeRelationship(parent_label='Node7413648', child_labels=[], edge_length=0.062766, taxon_label='Morelia carinata'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7413648', child_labels=[], edge_length=0.079145, taxon_label='Morelia boeleni'),
            'Node7413648' : NodeRelationship(parent_label='Node7413456', child_labels=['Morelia carinata','Morelia boeleni'], edge_length=0.013392, taxon_label=None),
            'Node7413456' : NodeRelationship(parent_label='Node7413424', child_labels=['Antaresia stimsoni','Node7413648'], edge_length=0.010508, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7413424', child_labels=[], edge_length=0.083409, taxon_label='Antaresia perthensis'),
            'Node7413424' : NodeRelationship(parent_label='Node7412624', child_labels=['Node7413456','Antaresia perthensis'], edge_length=0.005897, taxon_label=None),
            'Node7412624' : NodeRelationship(parent_label='Node7413168', child_labels=['Node7412912','Node7413424'], edge_length=0.0, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7413168', child_labels=[], edge_length=0.072091, taxon_label='Morelia viridis'),
            'Node7413168' : NodeRelationship(parent_label='Node7412368', child_labels=['Node7412624','Morelia viridis'], edge_length=0.008523, taxon_label=None),
            'Node7412368' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7412496','Node7413168'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7504176', child_labels=[], edge_length=0.062297, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7504304', child_labels=[], edge_length=0.069165, taxon_label='Bothrochilus boa'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7504304', child_labels=[], edge_length=0.073495, taxon_label='Liasis fuscus'),
            'Node7504304' : NodeRelationship(parent_label='Node7504176', child_labels=['Bothrochilus boa','Liasis fuscus'], edge_length=0.009899, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7504688', child_labels=[], edge_length=0.061022, taxon_label='Antaresia stimsoni'),
            'Morelia viridis' : NodeRelationship(parent_label='Node7505008', child_labels=[], edge_length=0.071291, taxon_label='Morelia viridis'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7505104', child_labels=[], edge_length=0.098776, taxon_label='Python timoriensis'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7505104', child_labels=[], edge_length=0.075069, taxon_label='Morelia boeleni'),
            'Node7505104' : NodeRelationship(parent_label='Node7505008', child_labels=['Python timoriensis','Morelia boeleni'], edge_length=0.017245, taxon_label=None),
            'Node7505008' : NodeRelationship(parent_label='Node7504688', child_labels=['Morelia viridis','Node7505104'], edge_length=0.004698, taxon_label=None),
            'Node7504688' : NodeRelationship(parent_label='Node7504624', child_labels=['Antaresia stimsoni','Node7505008'], edge_length=0.003603, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7505136', child_labels=[], edge_length=0.075494, taxon_label='Antaresia maculosa'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7505136', child_labels=[], edge_length=0.063757, taxon_label='Morelia oenpelliensis'),
            'Node7505136' : NodeRelationship(parent_label='Node7504624', child_labels=['Antaresia maculosa','Morelia oenpelliensis'], edge_length=0.017571, taxon_label=None),
            'Node7504624' : NodeRelationship(parent_label='Node7504368', child_labels=['Node7504688','Node7505136'], edge_length=0.000386, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7504368', child_labels=[], edge_length=0.070614, taxon_label='Morelia bredli'),
            'Node7504368' : NodeRelationship(parent_label='Node7504592', child_labels=['Node7504624','Morelia bredli'], edge_length=0.00713, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7504592', child_labels=[], edge_length=0.079534, taxon_label='Morelia carinata'),
            'Node7504592' : NodeRelationship(parent_label='Node7504528', child_labels=['Node7504368','Morelia carinata'], edge_length=0.010249, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7504528', child_labels=[], edge_length=0.114381, taxon_label='Python brongersmai'),
            'Node7504528' : NodeRelationship(parent_label='Node5593232', child_labels=['Node7504592','Python brongersmai'], edge_length=0.003672, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node5593232', child_labels=[], edge_length=0.080896, taxon_label='Antaresia perthensis'),
            'Node5593232' : NodeRelationship(parent_label='Node7504176', child_labels=['Node7504528','Antaresia perthensis'], edge_length=0.003212, taxon_label=None),
            'Node7504176' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7504304','Node5593232'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7505680', child_labels=[], edge_length=0.05678, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7506160', child_labels=[], edge_length=0.065928, taxon_label='Bothrochilus boa'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7506160', child_labels=[], edge_length=0.10077, taxon_label='Python timoriensis'),
            'Node7506160' : NodeRelationship(parent_label='Node7506096', child_labels=['Bothrochilus boa','Python timoriensis'], edge_length=0.013813, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7506096', child_labels=[], edge_length=0.083335, taxon_label='Antaresia perthensis'),
            'Node7506096' : NodeRelationship(parent_label='Node7506032', child_labels=['Node7506160','Antaresia perthensis'], edge_length=0.004377, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7506352', child_labels=[], edge_length=0.073392, taxon_label='Antaresia maculosa'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7506352', child_labels=[], edge_length=0.065748, taxon_label='Morelia oenpelliensis'),
            'Node7506352' : NodeRelationship(parent_label='Node7506032', child_labels=['Antaresia maculosa','Morelia oenpelliensis'], edge_length=0.016433, taxon_label=None),
            'Node7506032' : NodeRelationship(parent_label='Node7505968', child_labels=['Node7506096','Node7506352'], edge_length=0.005104, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7505968', child_labels=[], edge_length=0.07646, taxon_label='Morelia carinata'),
            'Node7505968' : NodeRelationship(parent_label='Node7505872', child_labels=['Node7506032','Morelia carinata'], edge_length=0.001866, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7505872', child_labels=[], edge_length=0.06129, taxon_label='Antaresia stimsoni'),
            'Node7505872' : NodeRelationship(parent_label='Node7412272', child_labels=['Node7505968','Antaresia stimsoni'], edge_length=0.012314, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7412272', child_labels=[], edge_length=0.081982, taxon_label='Morelia boeleni'),
            'Node7412272' : NodeRelationship(parent_label='Node7505584', child_labels=['Node7505872','Morelia boeleni'], edge_length=0.002833, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7506768', child_labels=[], edge_length=0.064242, taxon_label='Morelia viridis'),
            'Python brongersmai' : NodeRelationship(parent_label='Node7506768', child_labels=[], edge_length=0.109434, taxon_label='Python brongersmai'),
            'Node7506768' : NodeRelationship(parent_label='Node7505584', child_labels=['Morelia viridis','Python brongersmai'], edge_length=0.010734, taxon_label=None),
            'Node7505584' : NodeRelationship(parent_label='Node7505840', child_labels=['Node7412272','Node7506768'], edge_length=0.00182, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7505840', child_labels=[], edge_length=0.079156, taxon_label='Liasis fuscus'),
            'Node7505840' : NodeRelationship(parent_label='Node7505680', child_labels=['Node7505584','Liasis fuscus'], edge_length=0.014479, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7505680', child_labels=[], edge_length=0.060374, taxon_label='Morelia bredli'),
            'Node7505680' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7505840','Morelia bredli'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7507376', child_labels=[], edge_length=0.063004, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7507504', child_labels=[], edge_length=0.061413, taxon_label='Bothrochilus boa'),
            'Python brongersmai' : NodeRelationship(parent_label='Node7507504', child_labels=[], edge_length=0.099507, taxon_label='Python brongersmai'),
            'Node7507504' : NodeRelationship(parent_label='Node7507376', child_labels=['Bothrochilus boa','Python brongersmai'], edge_length=0.018059, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7507888', child_labels=[], edge_length=0.075877, taxon_label='Liasis fuscus'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7507888', child_labels=[], edge_length=0.071375, taxon_label='Antaresia maculosa'),
            'Node7507888' : NodeRelationship(parent_label='Node7507824', child_labels=['Liasis fuscus','Antaresia maculosa'], edge_length=0.014459, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7507824', child_labels=[], edge_length=0.054698, taxon_label='Antaresia stimsoni'),
            'Node7507824' : NodeRelationship(parent_label='Node7507568', child_labels=['Node7507888','Antaresia stimsoni'], edge_length=0.002182, taxon_label=None),
            'Morelia bredli' : NodeRelationship(parent_label='Node7561296', child_labels=[], edge_length=0.062291, taxon_label='Morelia bredli'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7561744', child_labels=[], edge_length=0.074463, taxon_label='Antaresia perthensis'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7561744', child_labels=[], edge_length=0.103887, taxon_label='Python timoriensis'),
            'Node7561744' : NodeRelationship(parent_label='Node7561296', child_labels=['Antaresia perthensis','Python timoriensis'], edge_length=0.00521, taxon_label=None),
            'Node7561296' : NodeRelationship(parent_label='Node7507568', child_labels=['Morelia bredli','Node7561744'], edge_length=0.011197, taxon_label=None),
            'Node7507568' : NodeRelationship(parent_label='Node7507792', child_labels=['Node7507824','Node7561296'], edge_length=0.010763, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7507792', child_labels=[], edge_length=0.070762, taxon_label='Morelia carinata'),
            'Node7507792' : NodeRelationship(parent_label='Node7507536', child_labels=['Node7507568','Morelia carinata'], edge_length=0.016533, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7507536', child_labels=[], edge_length=0.068818, taxon_label='Morelia oenpelliensis'),
            'Node7507536' : NodeRelationship(parent_label='Node7507728', child_labels=['Node7507792','Morelia oenpelliensis'], edge_length=0.003714, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node7561520', child_labels=[], edge_length=0.071307, taxon_label='Morelia viridis'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7561520', child_labels=[], edge_length=0.073392, taxon_label='Morelia boeleni'),
            'Node7561520' : NodeRelationship(parent_label='Node7507728', child_labels=['Morelia viridis','Morelia boeleni'], edge_length=0.010478, taxon_label=None),
            'Node7507728' : NodeRelationship(parent_label='Node7507376', child_labels=['Node7507536','Node7561520'], edge_length=0.00238, taxon_label=None),
            'Node7507376' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7507504','Node7507728'], edge_length=None, taxon_label=None),
        },
        {
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node7562320', child_labels=[], edge_length=0.053901, taxon_label='Aspidites ramsayi'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node7562480', child_labels=[], edge_length=0.07927, taxon_label='Bothrochilus boa'),
            'Morelia bredli' : NodeRelationship(parent_label='Node7562832', child_labels=[], edge_length=0.053956, taxon_label='Morelia bredli'),
            'Python timoriensis' : NodeRelationship(parent_label='Node7562832', child_labels=[], edge_length=0.098874, taxon_label='Python timoriensis'),
            'Node7562832' : NodeRelationship(parent_label='Node7562480', child_labels=['Morelia bredli','Python timoriensis'], edge_length=0.016769, taxon_label=None),
            'Node7562480' : NodeRelationship(parent_label='Node7355728', child_labels=['Bothrochilus boa','Node7562832'], edge_length=0.003041, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node7562864', child_labels=[], edge_length=0.075223, taxon_label='Liasis fuscus'),
            'Morelia viridis' : NodeRelationship(parent_label='Node7562864', child_labels=[], edge_length=0.065805, taxon_label='Morelia viridis'),
            'Node7562864' : NodeRelationship(parent_label='Node7563056', child_labels=['Liasis fuscus','Morelia viridis'], edge_length=0.010712, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node7563056', child_labels=[], edge_length=0.075934, taxon_label='Morelia carinata'),
            'Node7563056' : NodeRelationship(parent_label='Node7355728', child_labels=['Node7562864','Morelia carinata'], edge_length=0.008083, taxon_label=None),
            'Node7355728' : NodeRelationship(parent_label='Node7562032', child_labels=['Node7562480','Node7563056'], edge_length=0.00777, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node7563184', child_labels=[], edge_length=0.082425, taxon_label='Antaresia perthensis'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node7563024', child_labels=[], edge_length=0.075429, taxon_label='Antaresia maculosa'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node7563024', child_labels=[], edge_length=0.0787, taxon_label='Morelia boeleni'),
            'Node7563024' : NodeRelationship(parent_label='Node7563280', child_labels=['Antaresia maculosa','Morelia boeleni'], edge_length=0.010583, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node7563280', child_labels=[], edge_length=0.070289, taxon_label='Morelia oenpelliensis'),
            'Node7563280' : NodeRelationship(parent_label='Node7563184', child_labels=['Node7563024','Morelia oenpelliensis'], edge_length=0.008508, taxon_label=None),
            'Node7563184' : NodeRelationship(parent_label='Node7562032', child_labels=['Antaresia perthensis','Node7563280'], edge_length=0.00373, taxon_label=None),
            'Node7562032' : NodeRelationship(parent_label='Node7562448', child_labels=['Node7355728','Node7563184'], edge_length=0.002034, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node7562448', child_labels=[], edge_length=0.064652, taxon_label='Antaresia stimsoni'),
            'Node7562448' : NodeRelationship(parent_label='Node7562320', child_labels=['Node7562032','Antaresia stimsoni'], edge_length=0.013947, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node7562320', child_labels=[], edge_length=0.107868, taxon_label='Python brongersmai'),
            'Node7562320' : NodeRelationship(parent_label=None, child_labels=['Aspidites ramsayi','Node7562448','Python brongersmai'], edge_length=None, taxon_label=None),
        },
    ]
    return treelist_node_references

def reference_dna_dict():
    dna_dict = {}
    dna_dict['Aspidites ramsayi'] = "------------------------TTCGGCTCAATACTACTAACATGCTTAGGTYTACAAGTACTAACCGGCTTYTTCCTAGCCGTCCACTACACCGCAAACATTAACCTGGCATTCTCATCTATCGTTCACATCACCCGAGATGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGTATCTACATTCACATTGCACGAGGATTATACTACGGATCCTACCTTAACAAAGAAACCTGAATATCGGGTATTACATTACTCATTACACTAATAGCAACCGCCTTCTTCGGATATGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTCACCGCCGTACCATACCTAGGTACATCTCTAACAACCTGACTGTGAGGAGGATTCGCAATCAATGATCCCACCCTAACACGATTCTTTGCGCTACACTTCATCCTACCATTCGCAATCATCTCCTTGTCCTCACTACACATTATCTTACTTCACGAAGAAGGCTCCAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCGTACCACTCCCACAAAGACCTCCTCCTCCTAACGCTAATAATTATATCCTTATTTATTATCACCTCGTTCTTCCCAGACATCTTTAACGACCCCGACAATTTCTCAAAAGCCAACCCCCTAGTAACACCACAGCACATTAAACCAGAGTGGTACTTCCTATTTGCCTACGGCATCCTACGATCCATCCCCAACAAACTAGGAGGCGCACTAGCCCTAGTAATATCAATCATAATCCTATTTACCATTCCATTCACACACACAGCCTACCTTCGCCCTATAACCTTCCGCCCCCTATCACAACTAATATTCTGAACACTAGTCTCAACATTCATCACCATTACATGAGCCGCCACAAAACCAGTAGAACCACCATTCATTATAATTAGCCAAGTAACCTCAACACTATACTTCATATTCTTCTTATCAACACCCATCCTAGGATGAATAGAAAATAAAATAATAAACATTTCAT"
    dna_dict['Bothrochilus boa'] = "------------------------TTTGGCTCAATATTATTAACATGCCTGGCCCTACAAGTACTAACCGGCTTCTTCCTGGCCGTCCACTACACAGCAAACATCAACCTGGCATTCTCATCCATTATTCACATCACCCGAGATGTCCCATATGGCTGAATAATACAAAACCTGCACGCCATCGGAGCATCCATATTCTTCATTTGCGTATACATTCACATCGCACGAGGACTATACTACGGGTCATACCTAAACAAAGAAACCTGAATATCTGGCATTACCCTGCTCATCACACTAATAGCGACCGCCTTCTTTGGATATGTCCTCCCGTGAGGACAGATATCATTCTGAGCCGCAACCGTAATTACAAACCTGTTAACAGCAGTACCCTACCTGGGCACATCACTAACAACCTGGTTGTGAGGCGGGTTCGCAATCAATGACCCCACCCTAACACGATTCTTCGCACTACACTTCATCCTACCATTCGCAATCATCTCCCTATCCTCACTACACATCATCCTACTTCACGAGGAGGGCTCCAGCAACCCACTAGGAACCAACCCAGACATCGACAAAATCCCTTTCCACCCCTACCACTCCCACAAAGACTTTCTTCTTCTAACACTAATAACCCTATCCTTACTCATCATCGTCTTATTCTTCCCAGACATCTTTAACGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGATACTTCCTATTCGCCTACGGCATCCTACGTTCAATCCCCAATAAACTAGGAGGCGCACTAGCCTTAGTAATATCAATCATAATCCTATTCACCATCCCATTCACACACACAGCCCACCTTCGCCCCATAACCTTCCGTCCACTATCACAACTAATATTCTGAACATTAGTGTCAACATTTATCACTATCACATGGGCCGCCACAAAACCAGTAGAACCACCATTTATCACTATCAGCCAAACAACCTCAACACTATATTTTACATTCTTTTTACTTACCCCCATCCTAGGCTGAATAGAAAACAAAATAATAAAAACTCCCT"
    dna_dict['Liasis fuscus'] = "ACCAATATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTTTAGCCCTACAAGTATTAACCGGATTCTTCCTGGCTGTCCACTATACAGCAAATATTGACCTGGCATTCTCATCCATCATCCACATCACTCGAGACGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCAATATTCTTCATTTGTATCTACATCCACATCGCCCGAGGCCTATACTACGGATCATACCTCAACAAAGAAACCTGAATATCCGGCATCACCCTACTTATCACACTAATAGCAACCGCCTTCTTCGGGTACGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTTCTTACCGCCGTACCCTACCTAGGCACATCCTTGACAACCTGATTATGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGATTCTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCTTTATCCTCACTACACGTCATCCTCCTCCACGAGGAGGGGTCTAGCAACCCACTAGGGACTAACCCAGACATCGACAAAATCCCATTCCACCCTTACCACTCCCACAAAGACCTTCTCCTACTAACACTAATAATAATATCCCTACTCATTATTGTTTCCTTTTTCCCAGACATCTTCAACGACCCAGACAACTTCTCAAAAGCCAACCCCTTGGTAACACCACAACACATTAAACCAGAATGATACTTCCTGTTCGCCTACGGCATCCTACGATCTATTCCCAACAAACTTGGAGGAGCATTAGCTCTAGTAATATCAATCATAATCTTATTTTCTACCCCATTCACACACACAGCCCACCTCCGCCCTATAACTTTCCGCCCACTATCCCAACTAATATTCTGAACACTAATCTCAACATTCATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTCATCATCATCAGCCAAGTAACCTCAATACTATACTTCACATTTTTCCTATCCATCCCCATTCTAGGATGGGTAGAGAACAAAATTATAAACACCCCCT"
    dna_dict['Antaresia stimsoni'] = "------------------------TTCGGCTCAATACTATTAACATGTCTAGCCCTACAAGTATTAACCGGCTTTTTCCTAGCCGTTCATTATACAGCAAACATTAACCTAGCATTTTCATCCATCGTTCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTTATTTGTATTTATATTCACATCGCACGCGGACTATACTATGGATCCTACCTCAACAAAGAAACCTGAATATCCGGTATCACCCTGCTCATCACACTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTACCATACCTAGGCACATCGCTAACAACCTGACTGTGAGGGGGGTTCGCAATCAATGACCCCACCCTAACACGATTCTTCGCACTACACTTTATTCTACCATTCGCAATCATCTCCTTATCCTCCCTACACATTATCTTACTACACGAAGAAGGCTCAAGCAACCCACTAGGAACTAACCCAGACATCGACAAAATCCCATTCCACCCATACCACACCCACAAAGATCTACTCTTATTAACACTAATAGTTTTACTCCTATTCATTATCATTTCATTTTCCCCAGATATCTTTAATGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGGTACTTCCTATTTGCCTACGGTATCCTACGATCCATTCCCAACAAACTAGGAGGCGCACTAGCCTTAGTAATATCTATCATAATCCTATTCACCATCCCATTCACACACACAGCCCACCTACGCCCAATAACCTTCCGCCCACTCTCACAACTAACATTCTGAACACTGGTCTCAACATTCGCCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTTATCACCATCAGTCAAGTAACCTCAACACTATACTTCGCATTCTTCCTATCTATCCCAATTCTCGGATGAGTAGAAAACAAAATAATAAACATTTCAT"
    dna_dict['Morelia viridis'] = "------------------------TTCGGYTCAATACTATTAACATGCCTAGCCCTACAAGTATTAACCGGCTTCTTCCTAGCCGTTCACTACACAGSAAACATTAATCTAGCATTCTCATCCATCATCCACATCTCCCGAGATGTTCCATACGGTTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATCTACATCCATATTGCACGAGGATTATACTATGGATCCTACCTCAACAAAGAAACCTGAATATCCGGTATTACCCTACTCATCACACTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTTACCGCTGTACCCTACCTGGGTACATCACTAACAACCTGATTATGAGGCGGATTCGCAATCAACGACCCAACCCTAACACGATTTTTTGCACTGCACTTCATCCTACCATTCGCAATCATCTCTCTATCCTCACTTCATGTTATCCTACTCCACGAAGAAGGCTCCAGCAATCCACTAGGAACTAATCCAGATATCGATAAAATCCCATTTCATCCATACCACTCCTACAAAGACCTACTCCTACTAACACTAATAATCCTATTCTTATTCATCATCGTTTCATTCTTCCCAGACATTTTTAACGATCCGGACAACTTCTCAAAAGCTAACCCATTAGTAACACCACAACACATCAAACCAGAATGGTATTTCCTATTCGCCTACGGCATTCTACGATCCATCCCCAACAAACTAGGAGGCGCATTAGCCTTAGTAATATCAATCATAATCCTATTTACCATCCCATTTACACACACAGCCTACCTCCGCCCCATAACCTTTCGTCCACTATCACAATTAATATTTTGAACATTGGTTGCAACATTCGCCACTATTACATGGGCTGCCACAAAACCAGTAGAACCCCCATTTATCCTCATTAGCCAAGTAACTTCAACACTATATTTCACATTCTTCTTATCCATCCCAATTTTAGGATGAATGGAAAATAAAATAATAAACATCCCCT"
    dna_dict['Morelia bredli'] = "------------------------TTCGGCTCAATACTATTAACATGCCTAGCCCTGCAAATCCTAACCGGCTTCTTTTTAGCGGTCCACTACACAGCAAACATCAACCTAGCATTCTCATCCATCATCCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAASCTACACGCCATCGGAGCATCCCTATTCTTCATCTGCATCTACATCCATATCGCACGTGGGTTATACTACGGATCCTATCTCAACAAAGAAACCTGAATATCCGGTATTACCCTACTCATCACACTAATAGCAACTGCCTTCTTCGGTTATGTCCTTCCATGAGGACAAATATCATTCTRRGCCGCAACTGTAATTACAAATCTACTCACCGCCGTACCATACCTGGGCACATCACTAACAACCTGACTATGAGGCGGATTCGCAATCAATGACCCCACCTTAACACGATTCTTCGCGCTACACTTCATCCTACCATTCGCAATCATCTCTCTCTCTTCACTACACATTATCTTACTTCACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACACCCACAAAGACCTGCTCCTACTAACGCTCATAATCCTGTTCCTATTCATTATCGTCTCATTCCTCCCAGATATCTTCAATGACCCAGACAACTTCTCAAAAGCTAACCCCTTGGTAACACCACAACACATTAAACCAGAGTGGTACTTCCTATTTGCCTATGGCATTCTACGATCCATCCCCAATAAACTAGGAGGCGCACTAGCCCTAATAATATCGATCCTAATTCTATTCACGATCCCATTCATACACACAGCCTATCTCCGCCCTATAACCTTCCGCCCCCTGTCACAACTTATATTTTGAACACTAATCTCAACATTCGCCACCATTACATGAGCTGCCACAAAACCAGTAGAACCCCCATTCATCATCATTAGCCAAGTAACCTCAACACTATACTTCACATTCTTCCTAACCATCCCAATTCTAGGGTGAATAGAAAACAAAATAATAAACATCTCCT"
    dna_dict['Antaresia perthensis'] = "------------------------TTCGGMTCAATACTACTAACATGTTTAGCCTTACAAGTACTAACCGGCTTTTTTTTGGCCGTCCACTACACAGCAAACATCAACCTGGCATTTTCATCCATCATTCACATTACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTGCACGCCATCGGAGCATCCATATTCTTCATTTGCATCTACATTCACATTGCACGCGGACTCTACTACGGATCCTACCTCAACAAAGAAACCTGGATATCGGGAATTACCCTCCTCATCACACTGATAGCTACCGCCTTCTTCGGCTACGTCCTCCCATGAGGACAGATATCATTCTGAGCCGCAACAGTAATCACCAACCTACTCACCGCTGTACCCTACCTAGGCACATCACTAACAACCTGACTATGAGGGGGGTTCGCAATCAACGATCCCACCCTGACGCGATTCTTCGCACTACACTTCATCCTACCATTCGCAATCATTTCTCTATCATCCTTACACATTATCTTACTACACGAAGAGGGCTCCAGCAACCCACTAGGAACCAACCCAGACATTGACAAAATCCCATTCCACCCATACCACTCCCACAAAGACCTGCTCCTACTAACACTTATAATTTTATTCTTATTCATTATCCTATCGTTCTCCCCGGACATCTTCAACGACCCAGATAACTTCTCAAAAGCTAACCCCCTAGTAACACCACAACACATTAAACCAGAATGGTACTTCCTGTTCGCCTACGGAATTCTACGATCCATTCCAAATAAATTAGGAGGCGCACTAGCCCTTGTGATATCAATCATAATCCTATTTACCATCCCATTCATACACACCGCCCACCTACGCCCAATAACCTTCCGCCCACTTTCACAACTTATATTCTGAACACTAGTCTCAACATTTGCCACCATTACATGAGCCGCCACAAAACCGGTAGAGCCCCCATACATCCTCATTAGCCAAGTGACCGCAACACTATACTTCACATTCTTTCTATCCATCCCAATCCTAGGATGAATAGAAAACAAAATAATAAACACCTCCT"
    dna_dict['Python timoriensis'] = "------------------------TTCGGCTCACTACTATTAACATGTCTAGCCCTACAAGTATTAACTGGTTTTTTCCTAGCCGTTCACTACACAGCAAACATTAACCTGGCATTTTCATCCATCATTCACATCACCCGAGACGTCCCATACGGTTGAATGATACAAAACCTCCACGCCATCGGAGCATCCATATTTTTCATTTGTATTTACATCCACATCGCACGAGGCCTATACTACGGATCATATYTTAACAAAGAAACTTGAATATCAGGCATCACCCTACTCATCACATTAATAGCTACTGCTTTCTTCGGATATGTTCTTCCATGAGGACAAATATCATTCTGRGCCGCAACTGTAATTACAAACCTACTTACAGCCGTACCATACCTGGGCACATCATTAACAACCTGACTCTGAGGCGGATTTGCAATCAACGACCCAACTCTAACACGATTCTTCGCACTACACTTTATCCTACCATTCGCAATCATCTCACTATCCTCACTACACATTATCTTACTCCATGAAGAAGGTTCTAGCAACCCCCTAGGAACTAACCCAGACATCGATAAAATCCCATTCCACCCCTATCATTCCCACAAAGACTTCCTCTTACTAATACTAATAATTCTATTTTTATTCATTATCGTTTCATTCTTCCCAGATATTTTCAACGACCCAGACAATTTCTCAAAAGCTAACCCACTAGTGACACCACAACACATTAAACCAGAATGATACTTCCTATTTGCCTATGGCATCCTACGATCTATTCCCAATAAACTAGGAGGAGCCCTAGCCCTAGTAATATCTATTATAATTCTATTCACCATCCCATTCACACACACAGCCTATCTTCGTCCAATAACCTTTCGCCCTTTCTCACAATTCATATTCTGAACACTAATCACTACATTCATCACCATCACATGAGCCGCTACAAAACCTGTAGAACCACCATTCATTATTATCAGCCAAGCGACATCAACACTATATTTCACCTTCTTCATTTCAATCCCTCTTCTAGGCTGAATAGAAAACAAAATAATACACCTCAATT"
    dna_dict['Antaresia maculosa'] = "------------------------TTCGGCTCAATATTATTAACATGSCTGGSCCTGCAAGTATTAACCGGATTCTTCTTGGCCGTCCATTACACAGCAAATATCAACCTAGCATTTTCATCCATTATTCACATCACCCGAGATGTCCCATACGGCTGAATAATACAAAACCTACACGCCATYGGAGCCTCCATATTCTTCATTTGCATTTATATCCACATTKYACGAGGACTATACTACGGATCTTACCTCAATAAAGAAACCTGAATATCTGGTATCACCCTTCTTATCACACTAATAGCAACAGCCTTCTTCGGTTACGTTCTCCCATGAGGACARATATCATTCTGAGCCGCAACCGTAATCACAAACTTACTTACCGCCGTCCCATATCTAGGCACATCACTAACAACATGATTATGAGGGGGCTTTGCAATCAATGATCCCACCCTGACACGATTCTTCGCACTACACTTTATTCTACCATTCGCAATTATCTCCTTATCCTCACTACACATTATTTTACTTCACGAAGAGGGCTCCAGCAACCCACTAGGAACCAACCCAGACATTGACAAAATTCCATTTCACCCCTACCACTCTCATAAAGACCTGCTCCTACTCACACTAATAATTCTACTCTTACTCACTATCGTCTCATTTCTCCCAGACATCTTCAATGACCCAGACAACTTCTCAAAAGCTAATCCCCTAGTAACACCACAACACATCAAACCAGAGTGATATTTCCTATTCGCCTATGGCATTTTACGATCCATCCCTAATAAACTGGGAGGCGCACTAGCCCTAGTAATATCAATCATAATCCTATTTACCATTCCATTCACACACACAGCCCGCCTACGCCCCATAACCTTCCGCCCACTATCACAACTAATATTTTGAACATTAGTATCAACATTCGCCACCATTACATGAGCCGCCACAAAACCAGTAGAACCCCCATTTATCATCATCAGCCAAGCAACTTCAATTCTATACTTCACATTCTTCCTATCTACCCCAATTCTAGGGTGGGTGGAAAACAAAATAATAAATATCTCCT"
    dna_dict['Morelia carinata'] = "-----------------------CTTCGGCTCGATACTATTAACATGTTTAGCCCTACAAGTATTAACCGGCTTCTTCTTAGCTGTTCACTACACAGCAAACATTAACCTAGCATTCTCATCCATCATTCACATCACCCGAGACGTCCCATACGGCTGAATAATACAAAATCTGCACGCCATCGGAGCATCCATATTCTTCATCTGCATTTACATTCATATTGCACGAGGACTATACTATGGGTCTTACCTCAACAAAGAAACCTGAATATCTGGTATCACCCTGCTAATTATCCTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTACCCTACTTAGGCACATCACTAACAACCTGGCTTTGAGGCGGATTCGCAATCAATGACCCAACTCTAACACGATTCTTCGCACTACACTTCATCCTACCATTCGCTATTATCTCCCTATCCTCACTACACATTATCTTGCTTCACGAAGAAGGTTCTAGCAACCCCTTAGGAACCAACCCGGACATCGACAAAATCCCATTCCACCCCTATCACACCTACAAAGATCTTCTTCTACTAACAGTAATAATCCTATTTTTATTCATTATCGTTTCATTCTTCCCAGACATTTTTAATGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGATACTTTCTATTTGCCTACGGAATTCTACGATCCATCCCCAATAAACTGGGAGGAGCATTAGCCCTAGTAATATCAATTATAATCCTATTCACCATTCCATTTATACACACAGCCCATCTTCGCCCAATAACCTTCCGCCCACTATCMCAACTAATATTTTGAACACTAATCTCAACATTTATCACTATCACATGAGCTGCCACAAAACCAGTAGAACCCCCATTCATTACCATCAGCCAAGCAACTTCAGCGCTATACTTCACGTTCTTCCTAACCACCCCAATTCTAGGATGAGTAGAAAATAAAATARTAAACATTCCCT"
    dna_dict['Python brongersmai'] = "------------------------TTCGGTTCAATATTACTCACTTGCCTAGTCCTACAAGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACAGCAAACATCAACCTAGCATTTTCCTCTATTATACACATCACCCGCGACGTCCCATACGGCTGAATAATACAAAACTTACACGCYATCGGCGCATCTATATTTTTCATCTGCATCTATATCCACATCGCACGAGGACTATATTACGGCTCCTATCTCAATAAAGAAACCTGAATGTCTGGCATTACACTCCTCATCACACTAATAGCAACCGCTTTTTTCGGATATGTCCTCCCATGAGGACAGATGTCATTCTGAGCCGCAACCGTAATCACCAATCTACTAACTGCTGTACCATACCTAGGCACAACCCTAACAACCTGATTATGAGGAGGGTTCGCAATCAACGACCCCACCCTCACACGATTCTTTGCACTACACTTCATCCTACCTTTCGCAATCATCTCTTTATCATCACTACACATTATTCTCCTTCATGAAGAAGGATCTAGCAATCCACTAGGAACCAACCCCGACATCGACAAAATCCCATTCCACCCATACCACTCTCACAAAGACTTCCTCCTACTCACACTATTAATCCTTTTCCTATTTATCATTGTCTCCTTCTTCCCAGACATTTTTAATGATCCAGATAACTTCTCAAAAGCTAACCCCCTTGTCACACCCCAACACATTAAACCAGAATGATACTTCTTATTCGCTTACGGAATCCTACGATCCATCCCAAACAAACTAGGTGGCGCATTAGCATTAGTAATATCAATCATAATCTTATTTACCATCCCATTCACACACACAGCCCATCTTCGCCCTATAACCTTCCGACCGTTCTCACAACTAATATTCTGAACACTAATCTCAACATTCATCACCATCACATGGGCCGCCACAAAACCAGTAGAACCCCCATATATTATCATCAGCCAAGCAACTGCAGCATTATACTTCACCTTCTTCATCTCTACACCCCTCCTGGGCTGAATAGAAAATAAAATAACAAACACTCCCT"
    dna_dict['Morelia boeleni'] = "------------------------TTCGGCTCTATACTATTAACATGCTTAGGCTTACAAGTAATAACCGGCTTCTTCCTAGCCGTACACTACACAGCAAACATCAACTTAGCATTCTCATCCATCATCCACATCACCCGAGACGTCCCATACGGCTGAATAATACAAAACTTGCACGCTATCGGAGCATCTATATTCTTCATCTGCATTTACATCCACATCGCACGAGGGTTGTACTACGGATCATACCTTAACAAAGAAACCTGAATATCTGGCATTACCCTACTTATCACATTAATAGCAACTGCCTTCTTTGGATACGTTCTCCCATGAGGACAAATATCATTCTGAGSSGCWRCMGTWATCACAAACCTACTCACTGCCATCCCTTATCTAGGCACATCACTAACAACTTGACTATGAGGCGGATTCGCAATCAATGATCCTACACTAACACGATTTTTCGCACTACACTTCATCCTTCCATTCGCAATCATTTCCTTATCCTCACTACACATCATCCTACTCCACGAAGAAGGTTCCAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCATACCACTCCTACAAAGATCTTCTCCTACTAACATTAATAACCCTGTTCTTATTTATCATCGTCTCATTCTTCCCAGATATTTTTAACGACCCAGACAACTTTTCAAAAGCTAATCCCCTAGTAACACCACAACACATCAAACCTGAGTGATACTTTCTATTCGCCTATGGCATCCTACGATCCATCCCCAACAAACTAGGAGGTGCATTAGCCCTAGTAATATCAATCATGATCCTGTTTACCATCCCGTTTACACATACAGCCCACCTCCGTCCTATAACCTTCCGTCCACTCTCACAACTAATATTCTGAATATTAGTATCAACATTCATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCACCATTCATTACCATCAGTCAAGTAACCTCAACACTTTACTTCACATTCTTCTTATCCATCCCCATCCTAGGATGGATAGAAAACAAAATAATAGACATTCCAT"
    dna_dict['Morelia oenpelliensis'] = "------------------------TTCGGCTCAATACTATTAACATGCCTAGCCCTACAAGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACAGCAAATATCAACCTAGCATTTTCATCCATTATCCACATCACCCGTGACGTCCCATACGGTTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATTTATATTCACATCGCTCGAGGACTATACTATGGGTCATACCTTAACAAAGAAACCTGAATATCCGGTATCACCCTACTCATTACACTAATAGCAACCGCCTTCTTCGGATATGTTCTTCCATGAGGACAGATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCCGTACCATACCTAGGCACATCACTAACAACATGACTATGAGGCGGATTCGCAATCAATGACCCAACCCTAACCCGATTCTTTGCATTACACTTCATCCTACCATTTGCAATTATCTCTTTATCCTCACTACATATCATCCTACTCCATGAAGAAGGTTCCAGCAATCCATTAGGAACCAATCCTGACATTGACAAAATCCCATTCCACCCATACCACTCCCACAAAGACCTACTCCTATTAACACTAATAACCCTACTCCTATTCATTATTGTCTCATTCTTCCCAGATATTTTCAACGACCCAGACAACTTCTCAAAAGCCAATCCCATAGTAACACCACAACACATCAAACCAGAATGGTACTTCCTATTTGCCTACGGCATCCTACGATCCATCCCAAACAAACTAGGAGGCGCATTAGCCCTAGTAATATCAATTATAATTCTATTCACCGCCCCATTCACACATACAGCCTACCTACGCCCTATAACCTTTCGCCCACTTTCACAACTAATATTCTGAGCACTAGTATCAACATTCATCACCATTACATGAGCCGCCACAAAACCAGTAGAACCACCATTTATCATCATCAGCCAAACAACTGCAACACTATACTTCACATTCTTCTTATCCATCCCCATCACAGGATGAATTGAAAACAAAATAATAAACACCCACT"
    return dna_dict

def reference_dna_array():
    taxon_set = reference_taxon_set()
    dna = dendropy.DnaCharacterArray(taxon_set=taxon_set)
    sa = dna.default_state_alphabet
    assert len(dna.taxon_set) == 13
    assert dna.taxon_set.has_taxa(labels=["Aspidites ramsayi",
            "Bothrochilus boa", "Liasis fuscus", "Morelia boeleni", "Morelia viridis",
            "Antaresia maculosa", "Antaresia stimsoni", "Python timoriensis",
            "Morelia oenpelliensis", "Morelia bredli", "Antaresia perthensis",
            "Python brongersmai", "Morelia carinata"]), \
            [t.label for t in dna.taxon_set]

    dna_dict = reference_dna_dict()
    for t in dna.taxon_set:
        dna[t] = sa.get_states_as_cells(symbols=dna_dict[t.label])
    return dna

def _get_standard_state_alphabet(symbols):
    sa = dendropy.StateAlphabet()
    for symbol in symbols:
        sa.append(dendropy.StateAlphabetElement(symbol=symbol))
    sa.append(dendropy.StateAlphabetElement(symbol="?",
                                       multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                                       member_states=sa.get_states(symbols=symbols)))
    sa.append(dendropy.StateAlphabetElement(symbol="-",
                                       multistate=dendropy.StateAlphabetElement.AMBIGUOUS_STATE,
                                       member_states=sa.get_states(symbols=symbols)))
    return sa

def _get_standard_cells(col_type, symbols):
    cells = [dendropy.CharacterDataCell(value=s, column_type=col_type) for s in col_type.state_alphabet.get_states(symbols=symbols)]
    return cells

def reference_standard_array():
    taxon_set = reference_taxon_set()
    ca1 = dendropy.StandardCharacterArray(taxon_set=taxon_set)
    assert len(ca1.taxon_set) == 13
    sa1 = _get_standard_state_alphabet("012")
    sa2 = _get_standard_state_alphabet("XYZ")
    sa3 = _get_standard_state_alphabet("JKL")
    ca1.state_alphabets = [sa1, sa2, sa3]
    col_012 = dendropy.ColumnType(state_alphabet=sa1, label="COL_012")
    col_xyz = dendropy.ColumnType(state_alphabet=sa2, label="COL_XYZ")
    col_jkl = dendropy.ColumnType(state_alphabet=sa3, label="COL_JKL")
    ca1.column_types = [col_012, col_xyz, col_jkl]
    for t in taxon_set:
        ca1[t] = dendropy.CharacterDataVector(_get_standard_cells(col_012, "001122-??012")) \
               + dendropy.CharacterDataVector(_get_standard_cells(col_xyz, "XYZXYZ??-XXZ")) \
               + dendropy.CharacterDataVector(_get_standard_cells(col_jkl, "JKJLKL-??KJJ"))
    return ca1

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
