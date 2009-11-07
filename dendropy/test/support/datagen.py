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

def reference_tree_list_newick_source_path():
    return pathmap.tree_source_path("pythonidae_cytb.nexus.tre")

def reference_tree_list_newick_string():
    return """\
('Aspidites ramsayi':0.056823,(((((((('Bothrochilus boa':0.073578,('Liasis fuscus':0.074347,'Morelia boeleni':0.07628)Node100:0.003063)Node95:0.011155,'Morelia viridis':0.066566)Node93:0.010103,'Antaresia maculosa':0.080749)Node91:0.007681,'Antaresia stimsoni':0.061996)Node89:0.0,'Python timoriensis':0.113605)Node87:0.002264,'Morelia oenpelliensis':0.074047)Node85:0.007655,'Morelia bredli':0.065046)Node83:0.003495,('Antaresia perthensis':0.065004,'Python brongersmai':0.107706)Node126:0.018613)Node81:0.011159,'Morelia carinata':0.079321)Node87;
('Aspidites ramsayi':0.065297,((('Bothrochilus boa':0.073565,('Morelia boeleni':0.069885,'Morelia oenpelliensis':0.06247)Node153:0.007879)Node149:0.016022,'Morelia viridis':0.066015)Node147:0.010946,'Antaresia perthensis':0.079424)Node145:0.001444,(('Liasis fuscus':0.082049,('Antaresia stimsoni':0.056894,((('Python timoriensis':0.092077,'Morelia carinata':0.0579)Node177:0.015909,'Python brongersmai':0.114385)Node175:0.010752,'Antaresia maculosa':0.079792)Node173:0.004552)Node169:0.006907)Node165:0.008007,'Morelia bredli':0.063725)Node163:0.004077)Node137;
('Aspidites ramsayi':0.06111,((((('Bothrochilus boa':0.067908,'Antaresia perthensis':0.077196)Node205:0.008688,'Python brongersmai':0.1069)Node203:0.010367,(('Python timoriensis':0.109318,'Antaresia maculosa':0.081382)Node215:0.004129,'Morelia oenpelliensis':0.068424)Node213:0.007731)Node201:0.005167,('Liasis fuscus':0.079478,('Morelia viridis':0.062239,'Morelia bredli':0.060484)Node227:0.013155)Node223:0.00574)Node199:0.006988,'Antaresia stimsoni':0.063808)Node197:0.005408,('Morelia carinata':0.069501,'Morelia boeleni':0.073602)Node235:0.012121)Node187;
('Aspidites ramsayi':0.053028,((('Bothrochilus boa':0.072831,'Antaresia perthensis':0.073004)Node253:0.012413,'Morelia bredli':0.069933)Node251:0.002567,((((('Liasis fuscus':0.072339,'Morelia oenpelliensis':0.062058)Node269:0.011723,'Morelia viridis':0.065297)Node267:0.011519,(('Python timoriensis':0.093762,'Morelia carinata':0.055992)Node279:0.022912,'Antaresia maculosa':0.079945)Node277:0.007647)Node265:5.134e-05,'Python brongersmai':0.118976)Node263:0.004803,'Antaresia stimsoni':0.064721)Node261:0.005969)Node249:0.014376,'Morelia boeleni':0.076894)Node237;
('Aspidites ramsayi':0.06391,('Bothrochilus boa':0.071835,('Liasis fuscus':0.075788,'Morelia carinata':0.07116)Node305:0.010888)Node301:0.00451,(((('Antaresia stimsoni':0.048754,'Antaresia perthensis':0.069682)Node317:0.012356,'Antaresia maculosa':0.077058)Node315:0.013817,(('Morelia viridis':0.07026,('Morelia bredli':0.066182,'Morelia oenpelliensis':0.068166)Node331:0.001443)Node327:0.003147,'Python timoriensis':0.108243)Node325:0.012897)Node313:0.011944,('Python brongersmai':0.108072,'Morelia boeleni':0.073251)Node339:0.00612)Node311:0.005404)Node287;
('Aspidites ramsayi':0.06501,((('Bothrochilus boa':0.059984,'Python brongersmai':0.100555)Node357:0.017824,'Python timoriensis':0.102639)Node355:0.011722,('Morelia bredli':0.0538,'Antaresia maculosa':0.07843)Node365:0.014836)Node353:0.003644,((('Liasis fuscus':0.068438,'Morelia oenpelliensis':0.066815)Node375:0.015326,(('Antaresia stimsoni':0.057007,('Morelia carinata':0.062766,'Morelia boeleni':0.079145)Node387:0.013392)Node383:0.010508,'Antaresia perthensis':0.083409)Node381:0.005897)Node373:0.0,'Morelia viridis':0.072091)Node371:0.008523)Node337;
('Aspidites ramsayi':0.062297,('Bothrochilus boa':0.069165,'Liasis fuscus':0.073495)Node405:0.009899,(((((('Antaresia stimsoni':0.061022,('Morelia viridis':0.071291,('Python timoriensis':0.098776,'Morelia boeleni':0.075069)Node429:0.017245)Node425:0.004698)Node421:0.003603,('Antaresia maculosa':0.075494,'Morelia oenpelliensis':0.063757)Node435:0.017571)Node419:0.000386,'Morelia bredli':0.070614)Node417:0.00713,'Morelia carinata':0.079534)Node415:0.010249,'Python brongersmai':0.114381)Node413:0.003672,'Antaresia perthensis':0.080896)Node411:0.003212)Node387;
('Aspidites ramsayi':0.05678,(((((((('Bothrochilus boa':0.065928,'Python timoriensis':0.10077)Node471:0.013813,'Antaresia perthensis':0.083335)Node469:0.004377,('Antaresia maculosa':0.073392,'Morelia oenpelliensis':0.065748)Node479:0.016433)Node467:0.005104,'Morelia carinata':0.07646)Node465:0.001866,'Antaresia stimsoni':0.06129)Node463:0.012314,'Morelia boeleni':0.081982)Node461:0.002833,('Morelia viridis':0.064242,'Python brongersmai':0.109434)Node491:0.010734)Node459:0.00182,'Liasis fuscus':0.079156)Node457:0.014479,'Morelia bredli':0.060374)Node437;
('Aspidites ramsayi':0.063004,('Bothrochilus boa':0.061413,'Python brongersmai':0.099507)Node509:0.018059,(((((('Liasis fuscus':0.075877,'Antaresia maculosa':0.071375)Node525:0.014459,'Antaresia stimsoni':0.054698)Node523:0.002182,('Morelia bredli':0.062291,('Antaresia perthensis':0.074463,'Python timoriensis':0.103887)Node537:0.00521)Node533:0.011197)Node521:0.010763,'Morelia carinata':0.070762)Node519:0.016533,'Morelia oenpelliensis':0.068818)Node517:0.003714,('Morelia viridis':0.071307,'Morelia boeleni':0.073392)Node547:0.010478)Node515:0.00238)Node487;
('Aspidites ramsayi':0.053901,(((('Bothrochilus boa':0.07927,('Morelia bredli':0.053956,'Python timoriensis':0.098874)Node571:0.016769)Node567:0.003041,(('Liasis fuscus':0.075223,'Morelia viridis':0.065805)Node579:0.010712,'Morelia carinata':0.075934)Node577:0.008083)Node565:0.00777,('Antaresia perthensis':0.082425,(('Antaresia maculosa':0.075429,'Morelia boeleni':0.0787)Node593:0.010583,'Morelia oenpelliensis':0.070289)Node591:0.008508)Node587:0.00373)Node563:0.002034,'Antaresia stimsoni':0.064652)Node561:0.013947,'Python brongersmai':0.107868)Node537;
"""

def reference_tree_list_postorder_node_labels():
    return """\
['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Morelia boeleni','Node100','Node95','Morelia viridis','Node93','Antaresia maculosa','Node91','Antaresia stimsoni','Node89','Python timoriensis','Node87','Morelia oenpelliensis','Node85','Morelia bredli','Node83','Antaresia perthensis','Python brongersmai','Node126','Node81','Morelia carinata','Node87']
['Aspidites ramsayi','Bothrochilus boa','Morelia boeleni','Morelia oenpelliensis','Node153','Node149','Morelia viridis','Node147','Antaresia perthensis','Node145','Liasis fuscus','Antaresia stimsoni','Python timoriensis','Morelia carinata','Node177','Python brongersmai','Node175','Antaresia maculosa','Node173','Node169','Node165','Morelia bredli','Node163','Node137']
['Aspidites ramsayi','Bothrochilus boa','Antaresia perthensis','Node205','Python brongersmai','Node203','Python timoriensis','Antaresia maculosa','Node215','Morelia oenpelliensis','Node213','Node201','Liasis fuscus','Morelia viridis','Morelia bredli','Node227','Node223','Node199','Antaresia stimsoni','Node197','Morelia carinata','Morelia boeleni','Node235','Node187']
['Aspidites ramsayi','Bothrochilus boa','Antaresia perthensis','Node253','Morelia bredli','Node251','Liasis fuscus','Morelia oenpelliensis','Node269','Morelia viridis','Node267','Python timoriensis','Morelia carinata','Node279','Antaresia maculosa','Node277','Node265','Python brongersmai','Node263','Antaresia stimsoni','Node261','Node249','Morelia boeleni','Node237']
['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Morelia carinata','Node305','Node301','Antaresia stimsoni','Antaresia perthensis','Node317','Antaresia maculosa','Node315','Morelia viridis','Morelia bredli','Morelia oenpelliensis','Node331','Node327','Python timoriensis','Node325','Node313','Python brongersmai','Morelia boeleni','Node339','Node311','Node287']
['Aspidites ramsayi','Bothrochilus boa','Python brongersmai','Node357','Python timoriensis','Node355','Morelia bredli','Antaresia maculosa','Node365','Node353','Liasis fuscus','Morelia oenpelliensis','Node375','Antaresia stimsoni','Morelia carinata','Morelia boeleni','Node387','Node383','Antaresia perthensis','Node381','Node373','Morelia viridis','Node371','Node337']
['Aspidites ramsayi','Bothrochilus boa','Liasis fuscus','Node405','Antaresia stimsoni','Morelia viridis','Python timoriensis','Morelia boeleni','Node429','Node425','Node421','Antaresia maculosa','Morelia oenpelliensis','Node435','Node419','Morelia bredli','Node417','Morelia carinata','Node415','Python brongersmai','Node413','Antaresia perthensis','Node411','Node387']
['Aspidites ramsayi','Bothrochilus boa','Python timoriensis','Node471','Antaresia perthensis','Node469','Antaresia maculosa','Morelia oenpelliensis','Node479','Node467','Morelia carinata','Node465','Antaresia stimsoni','Node463','Morelia boeleni','Node461','Morelia viridis','Python brongersmai','Node491','Node459','Liasis fuscus','Node457','Morelia bredli','Node437']
['Aspidites ramsayi','Bothrochilus boa','Python brongersmai','Node509','Liasis fuscus','Antaresia maculosa','Node525','Antaresia stimsoni','Node523','Morelia bredli','Antaresia perthensis','Python timoriensis','Node537','Node533','Node521','Morelia carinata','Node519','Morelia oenpelliensis','Node517','Morelia viridis','Morelia boeleni','Node547','Node515','Node487']
['Aspidites ramsayi','Bothrochilus boa','Morelia bredli','Python timoriensis','Node571','Node567','Liasis fuscus','Morelia viridis','Node579','Morelia carinata','Node577','Node565','Antaresia perthensis','Antaresia maculosa','Morelia boeleni','Node593','Morelia oenpelliensis','Node591','Node587','Node563','Antaresia stimsoni','Node561','Python brongersmai','Node537']
"""

def reference_tree_list_relationships():
    """
    Returns list, with each element a dictionary describing the node
    relationships for the corresponding tree in the reference tree list.
    The keys for the dictionary are the node labels, and the values are a tuple
    consisting of: (<PARENT NODE LABEL>, [<CHILD NODE LABELS>], <EDGE LENGTH>).
    """

    return [
            {
                'Aspidites ramsayi' : ['Node87', [], 0.056823],
                'Bothrochilus boa' : ['Node95', [], 0.073578],
                'Liasis fuscus' : ['Node100', [], 0.074347],
                'Morelia boeleni' : ['Node100', [], 0.07628],
                'Node100' : ['Node95', ['Liasis fuscus','Morelia boeleni'], 0.003063],
                'Node95' : ['Node93', ['Bothrochilus boa','Node100'], 0.011155],
                'Morelia viridis' : ['Node93', [], 0.066566],
                'Node93' : ['Node91', ['Node95','Morelia viridis'], 0.010103],
                'Antaresia maculosa' : ['Node91', [], 0.080749],
                'Node91' : ['Node89', ['Node93','Antaresia maculosa'], 0.007681],
                'Antaresia stimsoni' : ['Node89', [], 0.061996],
                'Node89' : ['Node87', ['Node91','Antaresia stimsoni'], 0.0],
                'Python timoriensis' : ['Node87', [], 0.113605],
                'Node87' : ['Node85', ['Node89','Python timoriensis'], 0.002264],
                'Morelia oenpelliensis' : ['Node85', [], 0.074047],
                'Node85' : ['Node83', ['Node87','Morelia oenpelliensis'], 0.007655],
                'Morelia bredli' : ['Node83', [], 0.065046],
                'Node83' : ['Node81', ['Node85','Morelia bredli'], 0.003495],
                'Antaresia perthensis' : ['Node126', [], 0.065004],
                'Python brongersmai' : ['Node126', [], 0.107706],
                'Node126' : ['Node81', ['Antaresia perthensis','Python brongersmai'], 0.018613],
                'Node81' : ['Node87', ['Node83','Node126'], 0.011159],
                'Morelia carinata' : ['Node87', [], 0.079321],
                'Node87' : [None, ['Aspidites ramsayi','Node81','Morelia carinata'], None],
            },
            {
                'Aspidites ramsayi' : ['Node137', [], 0.065297],
                'Bothrochilus boa' : ['Node149', [], 0.073565],
                'Morelia boeleni' : ['Node153', [], 0.069885],
                'Morelia oenpelliensis' : ['Node153', [], 0.06247],
                'Node153' : ['Node149', ['Morelia boeleni','Morelia oenpelliensis'], 0.007879],
                'Node149' : ['Node147', ['Bothrochilus boa','Node153'], 0.016022],
                'Morelia viridis' : ['Node147', [], 0.066015],
                'Node147' : ['Node145', ['Node149','Morelia viridis'], 0.010946],
                'Antaresia perthensis' : ['Node145', [], 0.079424],
                'Node145' : ['Node137', ['Node147','Antaresia perthensis'], 0.001444],
                'Liasis fuscus' : ['Node165', [], 0.082049],
                'Antaresia stimsoni' : ['Node169', [], 0.056894],
                'Python timoriensis' : ['Node177', [], 0.092077],
                'Morelia carinata' : ['Node177', [], 0.0579],
                'Node177' : ['Node175', ['Python timoriensis','Morelia carinata'], 0.015909],
                'Python brongersmai' : ['Node175', [], 0.114385],
                'Node175' : ['Node173', ['Node177','Python brongersmai'], 0.010752],
                'Antaresia maculosa' : ['Node173', [], 0.079792],
                'Node173' : ['Node169', ['Node175','Antaresia maculosa'], 0.004552],
                'Node169' : ['Node165', ['Antaresia stimsoni','Node173'], 0.006907],
                'Node165' : ['Node163', ['Liasis fuscus','Node169'], 0.008007],
                'Morelia bredli' : ['Node163', [], 0.063725],
                'Node163' : ['Node137', ['Node165','Morelia bredli'], 0.004077],
                'Node137' : [None, ['Aspidites ramsayi','Node145','Node163'], None],
            },
            {
                'Aspidites ramsayi' : ['Node187', [], 0.06111],
                'Bothrochilus boa' : ['Node205', [], 0.067908],
                'Antaresia perthensis' : ['Node205', [], 0.077196],
                'Node205' : ['Node203', ['Bothrochilus boa','Antaresia perthensis'], 0.008688],
                'Python brongersmai' : ['Node203', [], 0.1069],
                'Node203' : ['Node201', ['Node205','Python brongersmai'], 0.010367],
                'Python timoriensis' : ['Node215', [], 0.109318],
                'Antaresia maculosa' : ['Node215', [], 0.081382],
                'Node215' : ['Node213', ['Python timoriensis','Antaresia maculosa'], 0.004129],
                'Morelia oenpelliensis' : ['Node213', [], 0.068424],
                'Node213' : ['Node201', ['Node215','Morelia oenpelliensis'], 0.007731],
                'Node201' : ['Node199', ['Node203','Node213'], 0.005167],
                'Liasis fuscus' : ['Node223', [], 0.079478],
                'Morelia viridis' : ['Node227', [], 0.062239],
                'Morelia bredli' : ['Node227', [], 0.060484],
                'Node227' : ['Node223', ['Morelia viridis','Morelia bredli'], 0.013155],
                'Node223' : ['Node199', ['Liasis fuscus','Node227'], 0.00574],
                'Node199' : ['Node197', ['Node201','Node223'], 0.006988],
                'Antaresia stimsoni' : ['Node197', [], 0.063808],
                'Node197' : ['Node187', ['Node199','Antaresia stimsoni'], 0.005408],
                'Morelia carinata' : ['Node235', [], 0.069501],
                'Morelia boeleni' : ['Node235', [], 0.073602],
                'Node235' : ['Node187', ['Morelia carinata','Morelia boeleni'], 0.012121],
                'Node187' : [None, ['Aspidites ramsayi','Node197','Node235'], None],
            },
            {
                'Aspidites ramsayi' : ['Node237', [], 0.053028],
                'Bothrochilus boa' : ['Node253', [], 0.072831],
                'Antaresia perthensis' : ['Node253', [], 0.073004],
                'Node253' : ['Node251', ['Bothrochilus boa','Antaresia perthensis'], 0.012413],
                'Morelia bredli' : ['Node251', [], 0.069933],
                'Node251' : ['Node249', ['Node253','Morelia bredli'], 0.002567],
                'Liasis fuscus' : ['Node269', [], 0.072339],
                'Morelia oenpelliensis' : ['Node269', [], 0.062058],
                'Node269' : ['Node267', ['Liasis fuscus','Morelia oenpelliensis'], 0.011723],
                'Morelia viridis' : ['Node267', [], 0.065297],
                'Node267' : ['Node265', ['Node269','Morelia viridis'], 0.011519],
                'Python timoriensis' : ['Node279', [], 0.093762],
                'Morelia carinata' : ['Node279', [], 0.055992],
                'Node279' : ['Node277', ['Python timoriensis','Morelia carinata'], 0.022912],
                'Antaresia maculosa' : ['Node277', [], 0.079945],
                'Node277' : ['Node265', ['Node279','Antaresia maculosa'], 0.007647],
                'Node265' : ['Node263', ['Node267','Node277'], 5.134e-05],
                'Python brongersmai' : ['Node263', [], 0.118976],
                'Node263' : ['Node261', ['Node265','Python brongersmai'], 0.004803],
                'Antaresia stimsoni' : ['Node261', [], 0.064721],
                'Node261' : ['Node249', ['Node263','Antaresia stimsoni'], 0.005969],
                'Node249' : ['Node237', ['Node251','Node261'], 0.014376],
                'Morelia boeleni' : ['Node237', [], 0.076894],
                'Node237' : [None, ['Aspidites ramsayi','Node249','Morelia boeleni'], None],
            },
            {
                'Aspidites ramsayi' : ['Node287', [], 0.06391],
                'Bothrochilus boa' : ['Node301', [], 0.071835],
                'Liasis fuscus' : ['Node305', [], 0.075788],
                'Morelia carinata' : ['Node305', [], 0.07116],
                'Node305' : ['Node301', ['Liasis fuscus','Morelia carinata'], 0.010888],
                'Node301' : ['Node287', ['Bothrochilus boa','Node305'], 0.00451],
                'Antaresia stimsoni' : ['Node317', [], 0.048754],
                'Antaresia perthensis' : ['Node317', [], 0.069682],
                'Node317' : ['Node315', ['Antaresia stimsoni','Antaresia perthensis'], 0.012356],
                'Antaresia maculosa' : ['Node315', [], 0.077058],
                'Node315' : ['Node313', ['Node317','Antaresia maculosa'], 0.013817],
                'Morelia viridis' : ['Node327', [], 0.07026],
                'Morelia bredli' : ['Node331', [], 0.066182],
                'Morelia oenpelliensis' : ['Node331', [], 0.068166],
                'Node331' : ['Node327', ['Morelia bredli','Morelia oenpelliensis'], 0.001443],
                'Node327' : ['Node325', ['Morelia viridis','Node331'], 0.003147],
                'Python timoriensis' : ['Node325', [], 0.108243],
                'Node325' : ['Node313', ['Node327','Python timoriensis'], 0.012897],
                'Node313' : ['Node311', ['Node315','Node325'], 0.011944],
                'Python brongersmai' : ['Node339', [], 0.108072],
                'Morelia boeleni' : ['Node339', [], 0.073251],
                'Node339' : ['Node311', ['Python brongersmai','Morelia boeleni'], 0.00612],
                'Node311' : ['Node287', ['Node313','Node339'], 0.005404],
                'Node287' : [None, ['Aspidites ramsayi','Node301','Node311'], None],
            },
            {
                'Aspidites ramsayi' : ['Node337', [], 0.06501],
                'Bothrochilus boa' : ['Node357', [], 0.059984],
                'Python brongersmai' : ['Node357', [], 0.100555],
                'Node357' : ['Node355', ['Bothrochilus boa','Python brongersmai'], 0.017824],
                'Python timoriensis' : ['Node355', [], 0.102639],
                'Node355' : ['Node353', ['Node357','Python timoriensis'], 0.011722],
                'Morelia bredli' : ['Node365', [], 0.0538],
                'Antaresia maculosa' : ['Node365', [], 0.07843],
                'Node365' : ['Node353', ['Morelia bredli','Antaresia maculosa'], 0.014836],
                'Node353' : ['Node337', ['Node355','Node365'], 0.003644],
                'Liasis fuscus' : ['Node375', [], 0.068438],
                'Morelia oenpelliensis' : ['Node375', [], 0.066815],
                'Node375' : ['Node373', ['Liasis fuscus','Morelia oenpelliensis'], 0.015326],
                'Antaresia stimsoni' : ['Node383', [], 0.057007],
                'Morelia carinata' : ['Node387', [], 0.062766],
                'Morelia boeleni' : ['Node387', [], 0.079145],
                'Node387' : ['Node383', ['Morelia carinata','Morelia boeleni'], 0.013392],
                'Node383' : ['Node381', ['Antaresia stimsoni','Node387'], 0.010508],
                'Antaresia perthensis' : ['Node381', [], 0.083409],
                'Node381' : ['Node373', ['Node383','Antaresia perthensis'], 0.005897],
                'Node373' : ['Node371', ['Node375','Node381'], 0.0],
                'Morelia viridis' : ['Node371', [], 0.072091],
                'Node371' : ['Node337', ['Node373','Morelia viridis'], 0.008523],
                'Node337' : [None, ['Aspidites ramsayi','Node353','Node371'], None],
            },
            {
                'Aspidites ramsayi' : ['Node387', [], 0.062297],
                'Bothrochilus boa' : ['Node405', [], 0.069165],
                'Liasis fuscus' : ['Node405', [], 0.073495],
                'Node405' : ['Node387', ['Bothrochilus boa','Liasis fuscus'], 0.009899],
                'Antaresia stimsoni' : ['Node421', [], 0.061022],
                'Morelia viridis' : ['Node425', [], 0.071291],
                'Python timoriensis' : ['Node429', [], 0.098776],
                'Morelia boeleni' : ['Node429', [], 0.075069],
                'Node429' : ['Node425', ['Python timoriensis','Morelia boeleni'], 0.017245],
                'Node425' : ['Node421', ['Morelia viridis','Node429'], 0.004698],
                'Node421' : ['Node419', ['Antaresia stimsoni','Node425'], 0.003603],
                'Antaresia maculosa' : ['Node435', [], 0.075494],
                'Morelia oenpelliensis' : ['Node435', [], 0.063757],
                'Node435' : ['Node419', ['Antaresia maculosa','Morelia oenpelliensis'], 0.017571],
                'Node419' : ['Node417', ['Node421','Node435'], 0.000386],
                'Morelia bredli' : ['Node417', [], 0.070614],
                'Node417' : ['Node415', ['Node419','Morelia bredli'], 0.00713],
                'Morelia carinata' : ['Node415', [], 0.079534],
                'Node415' : ['Node413', ['Node417','Morelia carinata'], 0.010249],
                'Python brongersmai' : ['Node413', [], 0.114381],
                'Node413' : ['Node411', ['Node415','Python brongersmai'], 0.003672],
                'Antaresia perthensis' : ['Node411', [], 0.080896],
                'Node411' : ['Node387', ['Node413','Antaresia perthensis'], 0.003212],
                'Node387' : [None, ['Aspidites ramsayi','Node405','Node411'], None],
            },
            {
                'Aspidites ramsayi' : ['Node437', [], 0.05678],
                'Bothrochilus boa' : ['Node471', [], 0.065928],
                'Python timoriensis' : ['Node471', [], 0.10077],
                'Node471' : ['Node469', ['Bothrochilus boa','Python timoriensis'], 0.013813],
                'Antaresia perthensis' : ['Node469', [], 0.083335],
                'Node469' : ['Node467', ['Node471','Antaresia perthensis'], 0.004377],
                'Antaresia maculosa' : ['Node479', [], 0.073392],
                'Morelia oenpelliensis' : ['Node479', [], 0.065748],
                'Node479' : ['Node467', ['Antaresia maculosa','Morelia oenpelliensis'], 0.016433],
                'Node467' : ['Node465', ['Node469','Node479'], 0.005104],
                'Morelia carinata' : ['Node465', [], 0.07646],
                'Node465' : ['Node463', ['Node467','Morelia carinata'], 0.001866],
                'Antaresia stimsoni' : ['Node463', [], 0.06129],
                'Node463' : ['Node461', ['Node465','Antaresia stimsoni'], 0.012314],
                'Morelia boeleni' : ['Node461', [], 0.081982],
                'Node461' : ['Node459', ['Node463','Morelia boeleni'], 0.002833],
                'Morelia viridis' : ['Node491', [], 0.064242],
                'Python brongersmai' : ['Node491', [], 0.109434],
                'Node491' : ['Node459', ['Morelia viridis','Python brongersmai'], 0.010734],
                'Node459' : ['Node457', ['Node461','Node491'], 0.00182],
                'Liasis fuscus' : ['Node457', [], 0.079156],
                'Node457' : ['Node437', ['Node459','Liasis fuscus'], 0.014479],
                'Morelia bredli' : ['Node437', [], 0.060374],
                'Node437' : [None, ['Aspidites ramsayi','Node457','Morelia bredli'], None],
            },
            {
                'Aspidites ramsayi' : ['Node487', [], 0.063004],
                'Bothrochilus boa' : ['Node509', [], 0.061413],
                'Python brongersmai' : ['Node509', [], 0.099507],
                'Node509' : ['Node487', ['Bothrochilus boa','Python brongersmai'], 0.018059],
                'Liasis fuscus' : ['Node525', [], 0.075877],
                'Antaresia maculosa' : ['Node525', [], 0.071375],
                'Node525' : ['Node523', ['Liasis fuscus','Antaresia maculosa'], 0.014459],
                'Antaresia stimsoni' : ['Node523', [], 0.054698],
                'Node523' : ['Node521', ['Node525','Antaresia stimsoni'], 0.002182],
                'Morelia bredli' : ['Node533', [], 0.062291],
                'Antaresia perthensis' : ['Node537', [], 0.074463],
                'Python timoriensis' : ['Node537', [], 0.103887],
                'Node537' : ['Node533', ['Antaresia perthensis','Python timoriensis'], 0.00521],
                'Node533' : ['Node521', ['Morelia bredli','Node537'], 0.011197],
                'Node521' : ['Node519', ['Node523','Node533'], 0.010763],
                'Morelia carinata' : ['Node519', [], 0.070762],
                'Node519' : ['Node517', ['Node521','Morelia carinata'], 0.016533],
                'Morelia oenpelliensis' : ['Node517', [], 0.068818],
                'Node517' : ['Node515', ['Node519','Morelia oenpelliensis'], 0.003714],
                'Morelia viridis' : ['Node547', [], 0.071307],
                'Morelia boeleni' : ['Node547', [], 0.073392],
                'Node547' : ['Node515', ['Morelia viridis','Morelia boeleni'], 0.010478],
                'Node515' : ['Node487', ['Node517','Node547'], 0.00238],
                'Node487' : [None, ['Aspidites ramsayi','Node509','Node515'], None],
            },
            {
                'Aspidites ramsayi' : ['Node537', [], 0.053901],
                'Bothrochilus boa' : ['Node567', [], 0.07927],
                'Morelia bredli' : ['Node571', [], 0.053956],
                'Python timoriensis' : ['Node571', [], 0.098874],
                'Node571' : ['Node567', ['Morelia bredli','Python timoriensis'], 0.016769],
                'Node567' : ['Node565', ['Bothrochilus boa','Node571'], 0.003041],
                'Liasis fuscus' : ['Node579', [], 0.075223],
                'Morelia viridis' : ['Node579', [], 0.065805],
                'Node579' : ['Node577', ['Liasis fuscus','Morelia viridis'], 0.010712],
                'Morelia carinata' : ['Node577', [], 0.075934],
                'Node577' : ['Node565', ['Node579','Morelia carinata'], 0.008083],
                'Node565' : ['Node563', ['Node567','Node577'], 0.00777],
                'Antaresia perthensis' : ['Node587', [], 0.082425],
                'Antaresia maculosa' : ['Node593', [], 0.075429],
                'Morelia boeleni' : ['Node593', [], 0.0787],
                'Node593' : ['Node591', ['Antaresia maculosa','Morelia boeleni'], 0.010583],
                'Morelia oenpelliensis' : ['Node591', [], 0.070289],
                'Node591' : ['Node587', ['Node593','Morelia oenpelliensis'], 0.008508],
                'Node587' : ['Node563', ['Antaresia perthensis','Node591'], 0.00373],
                'Node563' : ['Node561', ['Node565','Node587'], 0.002034],
                'Antaresia stimsoni' : ['Node561', [], 0.064652],
                'Node561' : ['Node537', ['Node563','Antaresia stimsoni'], 0.013947],
                'Python brongersmai' : ['Node537', [], 0.107868],
                'Node537' : [None, ['Aspidites ramsayi','Node561','Python brongersmai'], None],
            },
        ]

def reference_tree_list():
    tree_list_4286432 = dendropy.TreeList(label=None, oid="DendroPy Reference Tree")
    tax_4787248 = tree_list_4286432.taxon_set.require_taxon(label="Aspidites ramsayi", oid="Taxon80")
    tax_4825808 = tree_list_4286432.taxon_set.require_taxon(label="Bothrochilus boa", oid="Taxon99")
    tax_4788112 = tree_list_4286432.taxon_set.require_taxon(label="Liasis fuscus", oid="Taxon104")
    tax_4826544 = tree_list_4286432.taxon_set.require_taxon(label="Morelia boeleni", oid="Taxon107")
    tax_4826736 = tree_list_4286432.taxon_set.require_taxon(label="Morelia viridis", oid="Taxon110")
    tax_4826864 = tree_list_4286432.taxon_set.require_taxon(label="Antaresia maculosa", oid="Taxon113")
    tax_4827024 = tree_list_4286432.taxon_set.require_taxon(label="Antaresia stimsoni", oid="Taxon116")
    tax_4827280 = tree_list_4286432.taxon_set.require_taxon(label="Python timoriensis", oid="Taxon119")
    tax_4827216 = tree_list_4286432.taxon_set.require_taxon(label="Morelia oenpelliensis", oid="Taxon122")
    tax_4827632 = tree_list_4286432.taxon_set.require_taxon(label="Morelia bredli", oid="Taxon125")
    tax_4827920 = tree_list_4286432.taxon_set.require_taxon(label="Antaresia perthensis", oid="Taxon130")
    tax_4828208 = tree_list_4286432.taxon_set.require_taxon(label="Python brongersmai", oid="Taxon133")
    tax_4828400 = tree_list_4286432.taxon_set.require_taxon(label="Morelia carinata", oid="Taxon136")
    tree_4787536 = dendropy.Tree(label="PAUP 1", taxon_set=tree_list_4286432.taxon_set, oid="Tree72")
    tree_list_4286432.append(tree_4787536, reindex_taxa=False)
    nd_4787824 = tree_4787536.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.056823, oid="Node78")
    nd_4787824.edge.oid = "Edge79"
    nd_4787952 = tree_4787536.seed_node.new_child(label=None, taxon=None, edge_length=0.011159, oid="Node81")
    nd_4787952.edge.oid = "Edge82"
    nd_4828432 = tree_4787536.seed_node.new_child(label=None, taxon=tax_4828400, edge_length=0.079321, oid="Node134")
    nd_4828432.edge.oid = "Edge135"
    nd_4788080 = nd_4787952.new_child(label=None, taxon=None, edge_length=0.003495, oid="Node83")
    nd_4788080.edge.oid = "Edge84"
    nd_4827408 = nd_4787952.new_child(label=None, taxon=None, edge_length=0.018613, oid="Node126")
    nd_4827408.edge.oid = "Edge127"
    nd_4788016 = nd_4788080.new_child(label=None, taxon=None, edge_length=0.007655, oid="Node85")
    nd_4788016.edge.oid = "Edge86"
    nd_4827536 = nd_4788080.new_child(label=None, taxon=tax_4827632, edge_length=0.065046, oid="Node123")
    nd_4827536.edge.oid = "Edge124"
    nd_4788144 = nd_4788016.new_child(label=None, taxon=None, edge_length=0.002264, oid="Node87")
    nd_4788144.edge.oid = "Edge88"
    nd_4827376 = nd_4788016.new_child(label=None, taxon=tax_4827216, edge_length=0.074047, oid="Node120")
    nd_4827376.edge.oid = "Edge121"
    nd_4825168 = nd_4788144.new_child(label=None, taxon=None, edge_length=0.0, oid="Node89")
    nd_4825168.edge.oid = "Edge90"
    nd_4826832 = nd_4788144.new_child(label=None, taxon=tax_4827280, edge_length=0.113605, oid="Node117")
    nd_4826832.edge.oid = "Edge118"
    nd_4825296 = nd_4825168.new_child(label=None, taxon=None, edge_length=0.007681, oid="Node91")
    nd_4825296.edge.oid = "Edge92"
    nd_4826576 = nd_4825168.new_child(label=None, taxon=tax_4827024, edge_length=0.061996, oid="Node114")
    nd_4826576.edge.oid = "Edge115"
    nd_4825424 = nd_4825296.new_child(label=None, taxon=None, edge_length=0.010103, oid="Node93")
    nd_4825424.edge.oid = "Edge94"
    nd_4826640 = nd_4825296.new_child(label=None, taxon=tax_4826864, edge_length=0.080749, oid="Node111")
    nd_4826640.edge.oid = "Edge112"
    nd_4825552 = nd_4825424.new_child(label=None, taxon=None, edge_length=0.011155, oid="Node95")
    nd_4825552.edge.oid = "Edge96"
    nd_4826768 = nd_4825424.new_child(label=None, taxon=tax_4826736, edge_length=0.066566, oid="Node108")
    nd_4826768.edge.oid = "Edge109"
    nd_4825680 = nd_4825552.new_child(label=None, taxon=tax_4825808, edge_length=0.073578, oid="Node97")
    nd_4825680.edge.oid = "Edge98"
    nd_4826160 = nd_4825552.new_child(label=None, taxon=None, edge_length=0.003063, oid="Node100")
    nd_4826160.edge.oid = "Edge101"
    nd_4826224 = nd_4826160.new_child(label=None, taxon=tax_4788112, edge_length=0.074347, oid="Node102")
    nd_4826224.edge.oid = "Edge103"
    nd_4826448 = nd_4826160.new_child(label=None, taxon=tax_4826544, edge_length=0.07628, oid="Node105")
    nd_4826448.edge.oid = "Edge106"
    nd_4827792 = nd_4827408.new_child(label=None, taxon=tax_4827920, edge_length=0.065004, oid="Node128")
    nd_4827792.edge.oid = "Edge129"
    nd_4828144 = nd_4827408.new_child(label=None, taxon=tax_4828208, edge_length=0.107706, oid="Node131")
    nd_4828144.edge.oid = "Edge132"
    tree_4828560 = dendropy.Tree(label="PAUP 2", taxon_set=tree_list_4286432.taxon_set, oid="Tree137")
    tree_list_4286432.append(tree_4828560, reindex_taxa=False)
    nd_4828848 = tree_4828560.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.065297, oid="Node143")
    nd_4828848.edge.oid = "Edge144"
    nd_4828976 = tree_4828560.seed_node.new_child(label=None, taxon=None, edge_length=0.001444, oid="Node145")
    nd_4828976.edge.oid = "Edge146"
    nd_4858800 = tree_4828560.seed_node.new_child(label=None, taxon=None, edge_length=0.004077, oid="Node163")
    nd_4858800.edge.oid = "Edge164"
    nd_4828592 = nd_4828976.new_child(label=None, taxon=None, edge_length=0.010946, oid="Node147")
    nd_4828592.edge.oid = "Edge148"
    nd_4829072 = nd_4828976.new_child(label=None, taxon=tax_4827920, edge_length=0.079424, oid="Node161")
    nd_4829072.edge.oid = "Edge162"
    nd_4787344 = nd_4828592.new_child(label=None, taxon=None, edge_length=0.016022, oid="Node149")
    nd_4787344.edge.oid = "Edge150"
    nd_4858704 = nd_4828592.new_child(label=None, taxon=tax_4826736, edge_length=0.066015, oid="Node159")
    nd_4858704.edge.oid = "Edge160"
    nd_4829104 = nd_4787344.new_child(label=None, taxon=tax_4825808, edge_length=0.073565, oid="Node151")
    nd_4829104.edge.oid = "Edge152"
    nd_4858288 = nd_4787344.new_child(label=None, taxon=None, edge_length=0.007879, oid="Node153")
    nd_4858288.edge.oid = "Edge154"
    nd_4858064 = nd_4858288.new_child(label=None, taxon=tax_4826544, edge_length=0.069885, oid="Node155")
    nd_4858064.edge.oid = "Edge156"
    nd_4858512 = nd_4858288.new_child(label=None, taxon=tax_4827216, edge_length=0.06247, oid="Node157")
    nd_4858512.edge.oid = "Edge158"
    nd_4858896 = nd_4858800.new_child(label=None, taxon=None, edge_length=0.008007, oid="Node165")
    nd_4858896.edge.oid = "Edge166"
    nd_4860048 = nd_4858800.new_child(label=None, taxon=tax_4827632, edge_length=0.063725, oid="Node187")
    nd_4860048.edge.oid = "Edge188"
    nd_4858960 = nd_4858896.new_child(label=None, taxon=tax_4788112, edge_length=0.082049, oid="Node167")
    nd_4858960.edge.oid = "Edge168"
    nd_4859280 = nd_4858896.new_child(label=None, taxon=None, edge_length=0.006907, oid="Node169")
    nd_4859280.edge.oid = "Edge170"
    nd_4859088 = nd_4859280.new_child(label=None, taxon=tax_4827024, edge_length=0.056894, oid="Node171")
    nd_4859088.edge.oid = "Edge172"
    nd_4859568 = nd_4859280.new_child(label=None, taxon=None, edge_length=0.004552, oid="Node173")
    nd_4859568.edge.oid = "Edge174"
    nd_4858928 = nd_4859568.new_child(label=None, taxon=None, edge_length=0.010752, oid="Node175")
    nd_4858928.edge.oid = "Edge176"
    nd_4859920 = nd_4859568.new_child(label=None, taxon=tax_4826864, edge_length=0.079792, oid="Node185")
    nd_4859920.edge.oid = "Edge186"
    nd_4828240 = nd_4858928.new_child(label=None, taxon=None, edge_length=0.015909, oid="Node177")
    nd_4828240.edge.oid = "Edge178"
    nd_4787632 = nd_4858928.new_child(label=None, taxon=tax_4828208, edge_length=0.114385, oid="Node183")
    nd_4787632.edge.oid = "Edge184"
    nd_4787696 = nd_4828240.new_child(label=None, taxon=tax_4827280, edge_length=0.092077, oid="Node179")
    nd_4787696.edge.oid = "Edge180"
    nd_4828656 = nd_4828240.new_child(label=None, taxon=tax_4828400, edge_length=0.0579, oid="Node181")
    nd_4828656.edge.oid = "Edge182"
    tree_4860240 = dendropy.Tree(label="PAUP 3", taxon_set=tree_list_4286432.taxon_set, oid="Tree189")
    tree_list_4286432.append(tree_4860240, reindex_taxa=False)
    nd_4860592 = tree_4860240.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.06111, oid="Node195")
    nd_4860592.edge.oid = "Edge196"
    nd_4860720 = tree_4860240.seed_node.new_child(label=None, taxon=None, edge_length=0.005408, oid="Node197")
    nd_4860720.edge.oid = "Edge198"
    nd_4883536 = tree_4860240.seed_node.new_child(label=None, taxon=None, edge_length=0.012121, oid="Node235")
    nd_4883536.edge.oid = "Edge236"
    nd_4860304 = nd_4860720.new_child(label=None, taxon=None, edge_length=0.006988, oid="Node199")
    nd_4860304.edge.oid = "Edge200"
    nd_4882768 = nd_4860720.new_child(label=None, taxon=tax_4827024, edge_length=0.063808, oid="Node233")
    nd_4882768.edge.oid = "Edge234"
    nd_4787664 = nd_4860304.new_child(label=None, taxon=None, edge_length=0.005167, oid="Node201")
    nd_4787664.edge.oid = "Edge202"
    nd_4883024 = nd_4860304.new_child(label=None, taxon=None, edge_length=0.00574, oid="Node223")
    nd_4883024.edge.oid = "Edge224"
    nd_4860848 = nd_4787664.new_child(label=None, taxon=None, edge_length=0.010367, oid="Node203")
    nd_4860848.edge.oid = "Edge204"
    nd_4861360 = nd_4787664.new_child(label=None, taxon=None, edge_length=0.007731, oid="Node213")
    nd_4861360.edge.oid = "Edge214"
    nd_4860976 = nd_4860848.new_child(label=None, taxon=None, edge_length=0.008688, oid="Node205")
    nd_4860976.edge.oid = "Edge206"
    nd_4861616 = nd_4860848.new_child(label=None, taxon=tax_4828208, edge_length=0.1069, oid="Node211")
    nd_4861616.edge.oid = "Edge212"
    nd_4861104 = nd_4860976.new_child(label=None, taxon=tax_4825808, edge_length=0.067908, oid="Node207")
    nd_4861104.edge.oid = "Edge208"
    nd_4861584 = nd_4860976.new_child(label=None, taxon=tax_4827920, edge_length=0.077196, oid="Node209")
    nd_4861584.edge.oid = "Edge210"
    nd_4861808 = nd_4861360.new_child(label=None, taxon=None, edge_length=0.004129, oid="Node215")
    nd_4861808.edge.oid = "Edge216"
    nd_4861680 = nd_4861360.new_child(label=None, taxon=tax_4827216, edge_length=0.068424, oid="Node221")
    nd_4861680.edge.oid = "Edge222"
    nd_4861872 = nd_4861808.new_child(label=None, taxon=tax_4827280, edge_length=0.109318, oid="Node217")
    nd_4861872.edge.oid = "Edge218"
    nd_4861840 = nd_4861808.new_child(label=None, taxon=tax_4826864, edge_length=0.081382, oid="Node219")
    nd_4861840.edge.oid = "Edge220"
    nd_4882896 = nd_4883024.new_child(label=None, taxon=tax_4788112, edge_length=0.079478, oid="Node225")
    nd_4882896.edge.oid = "Edge226"
    nd_4883216 = nd_4883024.new_child(label=None, taxon=None, edge_length=0.013155, oid="Node227")
    nd_4883216.edge.oid = "Edge228"
    nd_4882832 = nd_4883216.new_child(label=None, taxon=tax_4826736, edge_length=0.062239, oid="Node229")
    nd_4882832.edge.oid = "Edge230"
    nd_4883472 = nd_4883216.new_child(label=None, taxon=tax_4827632, edge_length=0.060484, oid="Node231")
    nd_4883472.edge.oid = "Edge232"
    nd_4883728 = nd_4883536.new_child(label=None, taxon=tax_4828400, edge_length=0.069501, oid="Node237")
    nd_4883728.edge.oid = "Edge238"
    nd_4884016 = nd_4883536.new_child(label=None, taxon=tax_4826544, edge_length=0.073602, oid="Node239")
    nd_4884016.edge.oid = "Edge240"
    tree_4883792 = dendropy.Tree(label="PAUP 4", taxon_set=tree_list_4286432.taxon_set, oid="Tree241")
    tree_list_4286432.append(tree_4883792, reindex_taxa=False)
    nd_4884400 = tree_4883792.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.053028, oid="Node247")
    nd_4884400.edge.oid = "Edge248"
    nd_4884528 = tree_4883792.seed_node.new_child(label=None, taxon=None, edge_length=0.014376, oid="Node249")
    nd_4884528.edge.oid = "Edge250"
    nd_5030736 = tree_4883792.seed_node.new_child(label=None, taxon=tax_4826544, edge_length=0.076894, oid="Node291")
    nd_5030736.edge.oid = "Edge292"
    nd_4884144 = nd_4884528.new_child(label=None, taxon=None, edge_length=0.002567, oid="Node251")
    nd_4884144.edge.oid = "Edge252"
    nd_4884912 = nd_4884528.new_child(label=None, taxon=None, edge_length=0.005969, oid="Node261")
    nd_4884912.edge.oid = "Edge262"
    nd_4828688 = nd_4884144.new_child(label=None, taxon=None, edge_length=0.012413, oid="Node253")
    nd_4828688.edge.oid = "Edge254"
    nd_4885168 = nd_4884144.new_child(label=None, taxon=tax_4827632, edge_length=0.069933, oid="Node259")
    nd_4885168.edge.oid = "Edge260"
    nd_4884656 = nd_4828688.new_child(label=None, taxon=tax_4825808, edge_length=0.072831, oid="Node255")
    nd_4884656.edge.oid = "Edge256"
    nd_4885136 = nd_4828688.new_child(label=None, taxon=tax_4827920, edge_length=0.073004, oid="Node257")
    nd_4885136.edge.oid = "Edge258"
    nd_4885360 = nd_4884912.new_child(label=None, taxon=None, edge_length=0.004803, oid="Node263")
    nd_4885360.edge.oid = "Edge264"
    nd_5030512 = nd_4884912.new_child(label=None, taxon=tax_4827024, edge_length=0.064721, oid="Node289")
    nd_5030512.edge.oid = "Edge290"
    nd_4885424 = nd_4885360.new_child(label=None, taxon=None, edge_length=5.134e-05, oid="Node265")
    nd_4885424.edge.oid = "Edge266"
    nd_5030480 = nd_4885360.new_child(label=None, taxon=tax_4828208, edge_length=0.118976, oid="Node287")
    nd_5030480.edge.oid = "Edge288"
    nd_4885392 = nd_4885424.new_child(label=None, taxon=None, edge_length=0.011519, oid="Node267")
    nd_4885392.edge.oid = "Edge268"
    nd_4885936 = nd_4885424.new_child(label=None, taxon=None, edge_length=0.007647, oid="Node277")
    nd_4885936.edge.oid = "Edge278"
    nd_4885552 = nd_4885392.new_child(label=None, taxon=None, edge_length=0.011723, oid="Node269")
    nd_4885552.edge.oid = "Edge270"
    nd_4886192 = nd_4885392.new_child(label=None, taxon=tax_4826736, edge_length=0.065297, oid="Node275")
    nd_4886192.edge.oid = "Edge276"
    nd_4885680 = nd_4885552.new_child(label=None, taxon=tax_4788112, edge_length=0.072339, oid="Node271")
    nd_4885680.edge.oid = "Edge272"
    nd_4886128 = nd_4885552.new_child(label=None, taxon=tax_4827216, edge_length=0.062058, oid="Node273")
    nd_4886128.edge.oid = "Edge274"
    nd_4886416 = nd_4885936.new_child(label=None, taxon=None, edge_length=0.022912, oid="Node279")
    nd_4886416.edge.oid = "Edge280"
    nd_4886288 = nd_4885936.new_child(label=None, taxon=tax_4826864, edge_length=0.079945, oid="Node285")
    nd_4886288.edge.oid = "Edge286"
    nd_4886448 = nd_4886416.new_child(label=None, taxon=tax_4827280, edge_length=0.093762, oid="Node281")
    nd_4886448.edge.oid = "Edge282"
    nd_4886384 = nd_4886416.new_child(label=None, taxon=tax_4828400, edge_length=0.055992, oid="Node283")
    nd_4886384.edge.oid = "Edge284"
    tree_5030608 = dendropy.Tree(label="PAUP 5", taxon_set=tree_list_4286432.taxon_set, oid="Tree293")
    tree_list_4286432.append(tree_5030608, reindex_taxa=False)
    nd_5031088 = tree_5030608.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.06391, oid="Node299")
    nd_5031088.edge.oid = "Edge300"
    nd_5031216 = tree_5030608.seed_node.new_child(label=None, taxon=None, edge_length=0.00451, oid="Node301")
    nd_5031216.edge.oid = "Edge302"
    nd_4884176 = tree_5030608.seed_node.new_child(label=None, taxon=None, edge_length=0.005404, oid="Node311")
    nd_4884176.edge.oid = "Edge312"
    nd_5030320 = nd_5031216.new_child(label=None, taxon=tax_4825808, edge_length=0.071835, oid="Node303")
    nd_5030320.edge.oid = "Edge304"
    nd_5031568 = nd_5031216.new_child(label=None, taxon=None, edge_length=0.010888, oid="Node305")
    nd_5031568.edge.oid = "Edge306"
    nd_5031344 = nd_5031568.new_child(label=None, taxon=tax_4788112, edge_length=0.075788, oid="Node307")
    nd_5031344.edge.oid = "Edge308"
    nd_5031792 = nd_5031568.new_child(label=None, taxon=tax_4828400, edge_length=0.07116, oid="Node309")
    nd_5031792.edge.oid = "Edge310"
    nd_4860432 = nd_4884176.new_child(label=None, taxon=None, edge_length=0.011944, oid="Node313")
    nd_4860432.edge.oid = "Edge314"
    nd_5033392 = nd_4884176.new_child(label=None, taxon=None, edge_length=0.00612, oid="Node339")
    nd_5033392.edge.oid = "Edge340"
    nd_4884208 = nd_4860432.new_child(label=None, taxon=None, edge_length=0.013817, oid="Node315")
    nd_4884208.edge.oid = "Edge316"
    nd_5032304 = nd_4860432.new_child(label=None, taxon=None, edge_length=0.012897, oid="Node325")
    nd_5032304.edge.oid = "Edge326"
    nd_4860336 = nd_4884208.new_child(label=None, taxon=None, edge_length=0.012356, oid="Node317")
    nd_4860336.edge.oid = "Edge318"
    nd_4860400 = nd_4884208.new_child(label=None, taxon=tax_4826864, edge_length=0.077058, oid="Node323")
    nd_4860400.edge.oid = "Edge324"
    nd_5031536 = nd_4860336.new_child(label=None, taxon=tax_4827024, edge_length=0.048754, oid="Node319")
    nd_5031536.edge.oid = "Edge320"
    nd_5032208 = nd_4860336.new_child(label=None, taxon=tax_4827920, edge_length=0.069682, oid="Node321")
    nd_5032208.edge.oid = "Edge322"
    nd_5032432 = nd_5032304.new_child(label=None, taxon=None, edge_length=0.003147, oid="Node327")
    nd_5032432.edge.oid = "Edge328"
    nd_5033264 = nd_5032304.new_child(label=None, taxon=tax_4827280, edge_length=0.108243, oid="Node337")
    nd_5033264.edge.oid = "Edge338"
    nd_5032496 = nd_5032432.new_child(label=None, taxon=tax_4826736, edge_length=0.07026, oid="Node329")
    nd_5032496.edge.oid = "Edge330"
    nd_5032816 = nd_5032432.new_child(label=None, taxon=None, edge_length=0.001443, oid="Node331")
    nd_5032816.edge.oid = "Edge332"
    nd_5032624 = nd_5032816.new_child(label=None, taxon=tax_4827632, edge_length=0.066182, oid="Node333")
    nd_5032624.edge.oid = "Edge334"
    nd_5033104 = nd_5032816.new_child(label=None, taxon=tax_4827216, edge_length=0.068166, oid="Node335")
    nd_5033104.edge.oid = "Edge336"
    nd_5032464 = nd_5033392.new_child(label=None, taxon=tax_4828208, edge_length=0.108072, oid="Node341")
    nd_5032464.edge.oid = "Edge342"
    nd_5033616 = nd_5033392.new_child(label=None, taxon=tax_4826544, edge_length=0.073251, oid="Node343")
    nd_5033616.edge.oid = "Edge344"
    tree_5032880 = dendropy.Tree(label="PAUP 6", taxon_set=tree_list_4286432.taxon_set, oid="Tree345")
    tree_list_4286432.append(tree_5032880, reindex_taxa=False)
    nd_5058608 = tree_5032880.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.06501, oid="Node351")
    nd_5058608.edge.oid = "Edge352"
    nd_5033200 = tree_5032880.seed_node.new_child(label=None, taxon=None, edge_length=0.003644, oid="Node353")
    nd_5033200.edge.oid = "Edge354"
    nd_5060016 = tree_5032880.seed_node.new_child(label=None, taxon=None, edge_length=0.008523, oid="Node371")
    nd_5060016.edge.oid = "Edge372"
    nd_5058800 = nd_5033200.new_child(label=None, taxon=None, edge_length=0.011722, oid="Node355")
    nd_5058800.edge.oid = "Edge356"
    nd_5059120 = nd_5033200.new_child(label=None, taxon=None, edge_length=0.014836, oid="Node365")
    nd_5059120.edge.oid = "Edge366"
    nd_5030928 = nd_5058800.new_child(label=None, taxon=None, edge_length=0.017824, oid="Node357")
    nd_5030928.edge.oid = "Edge358"
    nd_5059376 = nd_5058800.new_child(label=None, taxon=tax_4827280, edge_length=0.102639, oid="Node363")
    nd_5059376.edge.oid = "Edge364"
    nd_5058864 = nd_5030928.new_child(label=None, taxon=tax_4825808, edge_length=0.059984, oid="Node359")
    nd_5058864.edge.oid = "Edge360"
    nd_5059344 = nd_5030928.new_child(label=None, taxon=tax_4828208, edge_length=0.100555, oid="Node361")
    nd_5059344.edge.oid = "Edge362"
    nd_5059568 = nd_5059120.new_child(label=None, taxon=tax_4827632, edge_length=0.0538, oid="Node367")
    nd_5059568.edge.oid = "Edge368"
    nd_5059856 = nd_5059120.new_child(label=None, taxon=tax_4826864, edge_length=0.07843, oid="Node369")
    nd_5059856.edge.oid = "Edge370"
    nd_5059440 = nd_5060016.new_child(label=None, taxon=None, edge_length=0.0, oid="Node373")
    nd_5059440.edge.oid = "Edge374"
    nd_5061552 = nd_5060016.new_child(label=None, taxon=tax_4826736, edge_length=0.072091, oid="Node395")
    nd_5061552.edge.oid = "Edge396"
    nd_5059600 = nd_5059440.new_child(label=None, taxon=None, edge_length=0.015326, oid="Node375")
    nd_5059600.edge.oid = "Edge376"
    nd_5060528 = nd_5059440.new_child(label=None, taxon=None, edge_length=0.005897, oid="Node381")
    nd_5060528.edge.oid = "Edge382"
    nd_5059888 = nd_5059600.new_child(label=None, taxon=tax_4788112, edge_length=0.068438, oid="Node377")
    nd_5059888.edge.oid = "Edge378"
    nd_5060464 = nd_5059600.new_child(label=None, taxon=tax_4827216, edge_length=0.066815, oid="Node379")
    nd_5060464.edge.oid = "Edge380"
    nd_5060624 = nd_5060528.new_child(label=None, taxon=None, edge_length=0.010508, oid="Node383")
    nd_5060624.edge.oid = "Edge384"
    nd_5061424 = nd_5060528.new_child(label=None, taxon=tax_4827920, edge_length=0.083409, oid="Node393")
    nd_5061424.edge.oid = "Edge394"
    nd_5060656 = nd_5060624.new_child(label=None, taxon=tax_4827024, edge_length=0.057007, oid="Node385")
    nd_5060656.edge.oid = "Edge386"
    nd_5061008 = nd_5060624.new_child(label=None, taxon=None, edge_length=0.013392, oid="Node387")
    nd_5061008.edge.oid = "Edge388"
    nd_5060784 = nd_5061008.new_child(label=None, taxon=tax_4828400, edge_length=0.062766, oid="Node389")
    nd_5060784.edge.oid = "Edge390"
    nd_5061264 = nd_5061008.new_child(label=None, taxon=tax_4826544, edge_length=0.079145, oid="Node391")
    nd_5061264.edge.oid = "Edge392"
    tree_5060592 = dendropy.Tree(label="PAUP 7", taxon_set=tree_list_4286432.taxon_set, oid="Tree397")
    tree_list_4286432.append(tree_5060592, reindex_taxa=False)
    nd_5061904 = tree_5060592.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.062297, oid="Node403")
    nd_5061904.edge.oid = "Edge404"
    nd_5062032 = tree_5060592.seed_node.new_child(label=None, taxon=None, edge_length=0.009899, oid="Node405")
    nd_5062032.edge.oid = "Edge406"
    nd_4787408 = tree_5060592.seed_node.new_child(label=None, taxon=None, edge_length=0.003212, oid="Node411")
    nd_4787408.edge.oid = "Edge412"
    nd_5061296 = nd_5062032.new_child(label=None, taxon=tax_4825808, edge_length=0.069165, oid="Node407")
    nd_5061296.edge.oid = "Edge408"
    nd_5062384 = nd_5062032.new_child(label=None, taxon=tax_4788112, edge_length=0.073495, oid="Node409")
    nd_5062384.edge.oid = "Edge410"
    nd_5062416 = nd_4787408.new_child(label=None, taxon=None, edge_length=0.003672, oid="Node413")
    nd_5062416.edge.oid = "Edge414"
    nd_5033808 = nd_4787408.new_child(label=None, taxon=tax_4827920, edge_length=0.080896, oid="Node447")
    nd_5033808.edge.oid = "Edge448"
    nd_5062544 = nd_5062416.new_child(label=None, taxon=None, edge_length=0.010249, oid="Node415")
    nd_5062544.edge.oid = "Edge416"
    nd_5088560 = nd_5062416.new_child(label=None, taxon=tax_4828208, edge_length=0.114381, oid="Node445")
    nd_5088560.edge.oid = "Edge446"
    nd_5062480 = nd_5062544.new_child(label=None, taxon=None, edge_length=0.00713, oid="Node417")
    nd_5062480.edge.oid = "Edge418"
    nd_5089008 = nd_5062544.new_child(label=None, taxon=tax_4828400, edge_length=0.079534, oid="Node443")
    nd_5089008.edge.oid = "Edge444"
    nd_5087280 = nd_5062480.new_child(label=None, taxon=None, edge_length=0.000386, oid="Node419")
    nd_5087280.edge.oid = "Edge420"
    nd_5089072 = nd_5062480.new_child(label=None, taxon=tax_4827632, edge_length=0.070614, oid="Node441")
    nd_5089072.edge.oid = "Edge442"
    nd_5087408 = nd_5087280.new_child(label=None, taxon=None, edge_length=0.003603, oid="Node421")
    nd_5087408.edge.oid = "Edge422"
    nd_5062512 = nd_5087280.new_child(label=None, taxon=None, edge_length=0.017571, oid="Node435")
    nd_5062512.edge.oid = "Edge436"
    nd_5087536 = nd_5087408.new_child(label=None, taxon=tax_4827024, edge_length=0.061022, oid="Node423")
    nd_5087536.edge.oid = "Edge424"
    nd_5088016 = nd_5087408.new_child(label=None, taxon=None, edge_length=0.004698, oid="Node425")
    nd_5088016.edge.oid = "Edge426"
    nd_5087792 = nd_5088016.new_child(label=None, taxon=tax_4826736, edge_length=0.071291, oid="Node427")
    nd_5087792.edge.oid = "Edge428"
    nd_5088240 = nd_5088016.new_child(label=None, taxon=None, edge_length=0.017245, oid="Node429")
    nd_5088240.edge.oid = "Edge430"
    nd_5087664 = nd_5088240.new_child(label=None, taxon=tax_4827280, edge_length=0.098776, oid="Node431")
    nd_5087664.edge.oid = "Edge432"
    nd_5088528 = nd_5088240.new_child(label=None, taxon=tax_4826544, edge_length=0.075069, oid="Node433")
    nd_5088528.edge.oid = "Edge434"
    nd_5088656 = nd_5062512.new_child(label=None, taxon=tax_4826864, edge_length=0.075494, oid="Node437")
    nd_5088656.edge.oid = "Edge438"
    nd_5088912 = nd_5062512.new_child(label=None, taxon=tax_4827216, edge_length=0.063757, oid="Node439")
    nd_5088912.edge.oid = "Edge440"
    tree_5033360 = dendropy.Tree(label="PAUP 8", taxon_set=tree_list_4286432.taxon_set, oid="Tree449")
    tree_list_4286432.append(tree_5033360, reindex_taxa=False)
    nd_5089552 = tree_5033360.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.05678, oid="Node455")
    nd_5089552.edge.oid = "Edge456"
    nd_5033840 = tree_5033360.seed_node.new_child(label=None, taxon=None, edge_length=0.014479, oid="Node457")
    nd_5033840.edge.oid = "Edge458"
    nd_5116944 = tree_5033360.seed_node.new_child(label=None, taxon=tax_4827632, edge_length=0.060374, oid="Node499")
    nd_5116944.edge.oid = "Edge500"
    nd_5089136 = nd_5033840.new_child(label=None, taxon=None, edge_length=0.00182, oid="Node459")
    nd_5089136.edge.oid = "Edge460"
    nd_5116976 = nd_5033840.new_child(label=None, taxon=tax_4788112, edge_length=0.079156, oid="Node497")
    nd_5116976.edge.oid = "Edge498"
    nd_4884240 = nd_5089136.new_child(label=None, taxon=None, edge_length=0.002833, oid="Node461")
    nd_4884240.edge.oid = "Edge462"
    nd_5116144 = nd_5089136.new_child(label=None, taxon=None, edge_length=0.010734, oid="Node491")
    nd_5116144.edge.oid = "Edge492"
    nd_5089808 = nd_4884240.new_child(label=None, taxon=None, edge_length=0.012314, oid="Node463")
    nd_5089808.edge.oid = "Edge464"
    nd_5116272 = nd_4884240.new_child(label=None, taxon=tax_4826544, edge_length=0.081982, oid="Node489")
    nd_5116272.edge.oid = "Edge490"
    nd_5089936 = nd_5089808.new_child(label=None, taxon=None, edge_length=0.001866, oid="Node465")
    nd_5089936.edge.oid = "Edge466"
    nd_5091024 = nd_5089808.new_child(label=None, taxon=tax_4827024, edge_length=0.06129, oid="Node487")
    nd_5091024.edge.oid = "Edge488"
    nd_5090064 = nd_5089936.new_child(label=None, taxon=None, edge_length=0.005104, oid="Node467")
    nd_5090064.edge.oid = "Edge468"
    nd_5116208 = nd_5089936.new_child(label=None, taxon=tax_4828400, edge_length=0.07646, oid="Node485")
    nd_5116208.edge.oid = "Edge486"
    nd_5090192 = nd_5090064.new_child(label=None, taxon=None, edge_length=0.004377, oid="Node469")
    nd_5090192.edge.oid = "Edge470"
    nd_5090704 = nd_5090064.new_child(label=None, taxon=None, edge_length=0.016433, oid="Node479")
    nd_5090704.edge.oid = "Edge480"
    nd_5090320 = nd_5090192.new_child(label=None, taxon=None, edge_length=0.013813, oid="Node471")
    nd_5090320.edge.oid = "Edge472"
    nd_5090960 = nd_5090192.new_child(label=None, taxon=tax_4827920, edge_length=0.083335, oid="Node477")
    nd_5090960.edge.oid = "Edge478"
    nd_5090448 = nd_5090320.new_child(label=None, taxon=tax_4825808, edge_length=0.065928, oid="Node473")
    nd_5090448.edge.oid = "Edge474"
    nd_5090928 = nd_5090320.new_child(label=None, taxon=tax_4827280, edge_length=0.10077, oid="Node475")
    nd_5090928.edge.oid = "Edge476"
    nd_5091152 = nd_5090704.new_child(label=None, taxon=tax_4826864, edge_length=0.073392, oid="Node481")
    nd_5091152.edge.oid = "Edge482"
    nd_5091184 = nd_5090704.new_child(label=None, taxon=tax_4827216, edge_length=0.065748, oid="Node483")
    nd_5091184.edge.oid = "Edge484"
    nd_5116560 = nd_5116144.new_child(label=None, taxon=tax_4826736, edge_length=0.064242, oid="Node493")
    nd_5116560.edge.oid = "Edge494"
    nd_5116784 = nd_5116144.new_child(label=None, taxon=tax_4828208, edge_length=0.109434, oid="Node495")
    nd_5116784.edge.oid = "Edge496"
    tree_5117040 = dendropy.Tree(label="PAUP 9", taxon_set=tree_list_4286432.taxon_set, oid="Tree501")
    tree_list_4286432.append(tree_5117040, reindex_taxa=False)
    nd_5117456 = tree_5117040.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.063004, oid="Node507")
    nd_5117456.edge.oid = "Edge508"
    nd_5117584 = tree_5117040.seed_node.new_child(label=None, taxon=None, edge_length=0.018059, oid="Node509")
    nd_5117584.edge.oid = "Edge510"
    nd_5061776 = tree_5117040.seed_node.new_child(label=None, taxon=None, edge_length=0.00238, oid="Node515")
    nd_5061776.edge.oid = "Edge516"
    nd_5117104 = nd_5117584.new_child(label=None, taxon=tax_4825808, edge_length=0.061413, oid="Node511")
    nd_5117104.edge.oid = "Edge512"
    nd_5117936 = nd_5117584.new_child(label=None, taxon=tax_4828208, edge_length=0.099507, oid="Node513")
    nd_5117936.edge.oid = "Edge514"
    nd_5117968 = nd_5061776.new_child(label=None, taxon=None, edge_length=0.003714, oid="Node517")
    nd_5117968.edge.oid = "Edge518"
    nd_5119824 = nd_5061776.new_child(label=None, taxon=None, edge_length=0.010478, oid="Node547")
    nd_5119824.edge.oid = "Edge548"
    nd_5118096 = nd_5117968.new_child(label=None, taxon=None, edge_length=0.016533, oid="Node519")
    nd_5118096.edge.oid = "Edge520"
    nd_5119760 = nd_5117968.new_child(label=None, taxon=tax_4827216, edge_length=0.068818, oid="Node545")
    nd_5119760.edge.oid = "Edge546"
    nd_5118064 = nd_5118096.new_child(label=None, taxon=None, edge_length=0.010763, oid="Node521")
    nd_5118064.edge.oid = "Edge522"
    nd_5119248 = nd_5118096.new_child(label=None, taxon=tax_4828400, edge_length=0.070762, oid="Node543")
    nd_5119248.edge.oid = "Edge544"
    nd_5118224 = nd_5118064.new_child(label=None, taxon=None, edge_length=0.002182, oid="Node523")
    nd_5118224.edge.oid = "Edge524"
    nd_5118736 = nd_5118064.new_child(label=None, taxon=None, edge_length=0.011197, oid="Node533")
    nd_5118736.edge.oid = "Edge534"
    nd_5118352 = nd_5118224.new_child(label=None, taxon=None, edge_length=0.014459, oid="Node525")
    nd_5118352.edge.oid = "Edge526"
    nd_5118992 = nd_5118224.new_child(label=None, taxon=tax_4827024, edge_length=0.054698, oid="Node531")
    nd_5118992.edge.oid = "Edge532"
    nd_5118480 = nd_5118352.new_child(label=None, taxon=tax_4788112, edge_length=0.075877, oid="Node527")
    nd_5118480.edge.oid = "Edge528"
    nd_5118928 = nd_5118352.new_child(label=None, taxon=tax_4826864, edge_length=0.071375, oid="Node529")
    nd_5118928.edge.oid = "Edge530"
    nd_5119184 = nd_5118736.new_child(label=None, taxon=tax_4827632, edge_length=0.062291, oid="Node535")
    nd_5119184.edge.oid = "Edge536"
    nd_5119472 = nd_5118736.new_child(label=None, taxon=None, edge_length=0.00521, oid="Node537")
    nd_5119472.edge.oid = "Edge538"
    nd_5119216 = nd_5119472.new_child(label=None, taxon=tax_4827920, edge_length=0.074463, oid="Node539")
    nd_5119216.edge.oid = "Edge540"
    nd_5119728 = nd_5119472.new_child(label=None, taxon=tax_4827280, edge_length=0.103887, oid="Node541")
    nd_5119728.edge.oid = "Edge542"
    nd_5144688 = nd_5119824.new_child(label=None, taxon=tax_4826736, edge_length=0.071307, oid="Node549")
    nd_5144688.edge.oid = "Edge550"
    nd_5144944 = nd_5119824.new_child(label=None, taxon=tax_4826544, edge_length=0.073392, oid="Node551")
    nd_5144944.edge.oid = "Edge552"
    tree_5144752 = dendropy.Tree(label="PAUP 10", taxon_set=tree_list_4286432.taxon_set, oid="Tree553")
    tree_list_4286432.append(tree_5144752, reindex_taxa=False)
    nd_5145360 = tree_5144752.seed_node.new_child(label=None, taxon=tax_4787248, edge_length=0.053901, oid="Node559")
    nd_5145360.edge.oid = "Edge560"
    nd_5089392 = tree_5144752.seed_node.new_child(label=None, taxon=None, edge_length=0.013947, oid="Node561")
    nd_5089392.edge.oid = "Edge562"
    nd_5147280 = tree_5144752.seed_node.new_child(label=None, taxon=tax_4828208, edge_length=0.107868, oid="Node603")
    nd_5147280.edge.oid = "Edge604"
    nd_5145072 = nd_5089392.new_child(label=None, taxon=None, edge_length=0.002034, oid="Node563")
    nd_5145072.edge.oid = "Edge564"
    nd_5147632 = nd_5089392.new_child(label=None, taxon=tax_4827024, edge_length=0.064652, oid="Node601")
    nd_5147632.edge.oid = "Edge602"
    nd_5061744 = nd_5145072.new_child(label=None, taxon=None, edge_length=0.00777, oid="Node565")
    nd_5061744.edge.oid = "Edge566"
    nd_5146896 = nd_5145072.new_child(label=None, taxon=None, edge_length=0.00373, oid="Node587")
    nd_5146896.edge.oid = "Edge588"
    nd_5145616 = nd_5061744.new_child(label=None, taxon=None, edge_length=0.003041, oid="Node567")
    nd_5145616.edge.oid = "Edge568"
    nd_5146640 = nd_5061744.new_child(label=None, taxon=None, edge_length=0.008083, oid="Node577")
    nd_5146640.edge.oid = "Edge578"
    nd_5145744 = nd_5145616.new_child(label=None, taxon=tax_4825808, edge_length=0.07927, oid="Node569")
    nd_5145744.edge.oid = "Edge570"
    nd_5146224 = nd_5145616.new_child(label=None, taxon=None, edge_length=0.016769, oid="Node571")
    nd_5146224.edge.oid = "Edge572"
    nd_5146000 = nd_5146224.new_child(label=None, taxon=tax_4827632, edge_length=0.053956, oid="Node573")
    nd_5146000.edge.oid = "Edge574"
    nd_5146480 = nd_5146224.new_child(label=None, taxon=tax_4827280, edge_length=0.098874, oid="Node575")
    nd_5146480.edge.oid = "Edge576"
    nd_5146256 = nd_5146640.new_child(label=None, taxon=None, edge_length=0.010712, oid="Node579")
    nd_5146256.edge.oid = "Edge580"
    nd_5117328 = nd_5146640.new_child(label=None, taxon=tax_4828400, edge_length=0.075934, oid="Node585")
    nd_5117328.edge.oid = "Edge586"
    nd_5089360 = nd_5146256.new_child(label=None, taxon=tax_4788112, edge_length=0.075223, oid="Node581")
    nd_5089360.edge.oid = "Edge582"
    nd_5089200 = nd_5146256.new_child(label=None, taxon=tax_4826736, edge_length=0.065805, oid="Node583")
    nd_5089200.edge.oid = "Edge584"
    nd_5146768 = nd_5146896.new_child(label=None, taxon=tax_4827920, edge_length=0.082425, oid="Node589")
    nd_5146768.edge.oid = "Edge590"
    nd_5147120 = nd_5146896.new_child(label=None, taxon=None, edge_length=0.008508, oid="Node591")
    nd_5147120.edge.oid = "Edge592"
    nd_5146736 = nd_5147120.new_child(label=None, taxon=None, edge_length=0.010583, oid="Node593")
    nd_5146736.edge.oid = "Edge594"
    nd_5147536 = nd_5147120.new_child(label=None, taxon=tax_4827216, edge_length=0.070289, oid="Node599")
    nd_5147536.edge.oid = "Edge600"
    nd_5146864 = nd_5146736.new_child(label=None, taxon=tax_4826864, edge_length=0.075429, oid="Node595")
    nd_5146864.edge.oid = "Edge596"
    nd_5147504 = nd_5146736.new_child(label=None, taxon=tax_4826544, edge_length=0.0787, oid="Node597")
    nd_5147504.edge.oid = "Edge598"

    # set labels of nodes with taxa to taxon label, else oid (for consistent
    # identification in debugging)
    for t in tree_list_4286432:
        t.assign_node_labels_from_taxon_or_oid()
    return tree_list_4286432

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
