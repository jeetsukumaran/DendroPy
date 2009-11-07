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
import dendropy

def get_standard_four_taxon_tree():
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

def get_canonical_tree_list():
    tree_list_4286432 = dendropy.TreeList(label=None, oid="TreeList70")
    tax_4787280 = tree_list_4286432.taxon_set.require_taxon(label="Aspidites ramsayi", oid="Taxon80")
    tax_4825840 = tree_list_4286432.taxon_set.require_taxon(label="Bothrochilus boa", oid="Taxon99")
    tax_4788144 = tree_list_4286432.taxon_set.require_taxon(label="Liasis fuscus", oid="Taxon104")
    tax_4826576 = tree_list_4286432.taxon_set.require_taxon(label="Morelia boeleni", oid="Taxon107")
    tax_4826768 = tree_list_4286432.taxon_set.require_taxon(label="Morelia viridis", oid="Taxon110")
    tax_4826896 = tree_list_4286432.taxon_set.require_taxon(label="Antaresia maculosa", oid="Taxon113")
    tax_4827056 = tree_list_4286432.taxon_set.require_taxon(label="Antaresia stimsoni", oid="Taxon116")
    tax_4827312 = tree_list_4286432.taxon_set.require_taxon(label="Python timoriensis", oid="Taxon119")
    tax_4827248 = tree_list_4286432.taxon_set.require_taxon(label="Morelia oenpelliensis", oid="Taxon122")
    tax_4827664 = tree_list_4286432.taxon_set.require_taxon(label="Morelia bredli", oid="Taxon125")
    tax_4827952 = tree_list_4286432.taxon_set.require_taxon(label="Antaresia perthensis", oid="Taxon130")
    tax_4828240 = tree_list_4286432.taxon_set.require_taxon(label="Python brongersmai", oid="Taxon133")
    tax_4828432 = tree_list_4286432.taxon_set.require_taxon(label="Morelia carinata", oid="Taxon136")
    tree_4787568 = dendropy.Tree(label="PAUP 1", taxon_set=tree_list_4286432.taxon_set, oid="Tree72")
    tree_list_4286432.append(tree_4787568, reindex_taxa=False)
    nd_4787856 = tree_4787568.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.056823, oid="Node78")
    nd_4787984 = tree_4787568.seed_node.new_child(label=None, taxon=None, edge_length=0.011159, oid="Node81")
    nd_4828464 = tree_4787568.seed_node.new_child(label=None, taxon=tax_4828432, edge_length=0.079321, oid="Node134")
    nd_4828464.edge.oid = "Edge135"
    nd_4828464.edge.oid = "Edge135"
    nd_4788112 = nd_4787984.new_child(label=None, taxon=None, edge_length=0.003495, oid="Node83")
    nd_4827440 = nd_4787984.new_child(label=None, taxon=None, edge_length=0.018613, oid="Node126")
    nd_4827440.edge.oid = "Edge127"
    nd_4788048 = nd_4788112.new_child(label=None, taxon=None, edge_length=0.007655, oid="Node85")
    nd_4827568 = nd_4788112.new_child(label=None, taxon=tax_4827664, edge_length=0.065046, oid="Node123")
    nd_4827568.edge.oid = "Edge124"
    nd_4788176 = nd_4788048.new_child(label=None, taxon=None, edge_length=0.002264, oid="Node87")
    nd_4827408 = nd_4788048.new_child(label=None, taxon=tax_4827248, edge_length=0.074047, oid="Node120")
    nd_4827408.edge.oid = "Edge121"
    nd_4825200 = nd_4788176.new_child(label=None, taxon=None, edge_length=0.0, oid="Node89")
    nd_4826864 = nd_4788176.new_child(label=None, taxon=tax_4827312, edge_length=0.113605, oid="Node117")
    nd_4826864.edge.oid = "Edge118"
    nd_4825328 = nd_4825200.new_child(label=None, taxon=None, edge_length=0.007681, oid="Node91")
    nd_4826608 = nd_4825200.new_child(label=None, taxon=tax_4827056, edge_length=0.061996, oid="Node114")
    nd_4826608.edge.oid = "Edge115"
    nd_4825456 = nd_4825328.new_child(label=None, taxon=None, edge_length=0.010103, oid="Node93")
    nd_4826672 = nd_4825328.new_child(label=None, taxon=tax_4826896, edge_length=0.080749, oid="Node111")
    nd_4826672.edge.oid = "Edge112"
    nd_4825584 = nd_4825456.new_child(label=None, taxon=None, edge_length=0.011155, oid="Node95")
    nd_4826800 = nd_4825456.new_child(label=None, taxon=tax_4826768, edge_length=0.066566, oid="Node108")
    nd_4826800.edge.oid = "Edge109"
    nd_4825712 = nd_4825584.new_child(label=None, taxon=tax_4825840, edge_length=0.073578, oid="Node97")
    nd_4826192 = nd_4825584.new_child(label=None, taxon=None, edge_length=0.003063, oid="Node100")
    nd_4826192.edge.oid = "Edge101"
    nd_4826192.edge.oid = "Edge101"
    nd_4826256 = nd_4826192.new_child(label=None, taxon=tax_4788144, edge_length=0.074347, oid="Node102")
    nd_4826480 = nd_4826192.new_child(label=None, taxon=tax_4826576, edge_length=0.07628, oid="Node105")
    nd_4826480.edge.oid = "Edge106"
    nd_4826480.edge.oid = "Edge106"
    nd_4826480.edge.oid = "Edge106"
    nd_4826480.edge.oid = "Edge106"
    nd_4826480.edge.oid = "Edge106"
    nd_4826480.edge.oid = "Edge106"
    nd_4826480.edge.oid = "Edge106"
    nd_4826480.edge.oid = "Edge106"
    nd_4826480.edge.oid = "Edge106"
    nd_4827824 = nd_4827440.new_child(label=None, taxon=tax_4827952, edge_length=0.065004, oid="Node128")
    nd_4828176 = nd_4827440.new_child(label=None, taxon=tax_4828240, edge_length=0.107706, oid="Node131")
    nd_4828176.edge.oid = "Edge132"
    nd_4828176.edge.oid = "Edge132"
    nd_4828176.edge.oid = "Edge132"
    nd_4828176.edge.oid = "Edge132"
    tree_4828592 = dendropy.Tree(label="PAUP 2", taxon_set=tree_list_4286432.taxon_set, oid="Tree137")
    tree_list_4286432.append(tree_4828592, reindex_taxa=False)
    nd_4828880 = tree_4828592.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.065297, oid="Node143")
    nd_4829008 = tree_4828592.seed_node.new_child(label=None, taxon=None, edge_length=0.001444, oid="Node145")
    nd_4858832 = tree_4828592.seed_node.new_child(label=None, taxon=None, edge_length=0.004077, oid="Node163")
    nd_4858832.edge.oid = "Edge164"
    nd_4858832.edge.oid = "Edge164"
    nd_4828624 = nd_4829008.new_child(label=None, taxon=None, edge_length=0.010946, oid="Node147")
    nd_4829104 = nd_4829008.new_child(label=None, taxon=tax_4827952, edge_length=0.079424, oid="Node161")
    nd_4829104.edge.oid = "Edge162"
    nd_4787376 = nd_4828624.new_child(label=None, taxon=None, edge_length=0.016022, oid="Node149")
    nd_4858736 = nd_4828624.new_child(label=None, taxon=tax_4826768, edge_length=0.066015, oid="Node159")
    nd_4858736.edge.oid = "Edge160"
    nd_4829136 = nd_4787376.new_child(label=None, taxon=tax_4825840, edge_length=0.073565, oid="Node151")
    nd_4858320 = nd_4787376.new_child(label=None, taxon=None, edge_length=0.007879, oid="Node153")
    nd_4858320.edge.oid = "Edge154"
    nd_4858320.edge.oid = "Edge154"
    nd_4858096 = nd_4858320.new_child(label=None, taxon=tax_4826576, edge_length=0.069885, oid="Node155")
    nd_4858544 = nd_4858320.new_child(label=None, taxon=tax_4827248, edge_length=0.06247, oid="Node157")
    nd_4858544.edge.oid = "Edge158"
    nd_4858544.edge.oid = "Edge158"
    nd_4858544.edge.oid = "Edge158"
    nd_4858544.edge.oid = "Edge158"
    nd_4858544.edge.oid = "Edge158"
    nd_4858928 = nd_4858832.new_child(label=None, taxon=None, edge_length=0.008007, oid="Node165")
    nd_4860080 = nd_4858832.new_child(label=None, taxon=tax_4827664, edge_length=0.063725, oid="Node187")
    nd_4860080.edge.oid = "Edge188"
    nd_4858992 = nd_4858928.new_child(label=None, taxon=tax_4788144, edge_length=0.082049, oid="Node167")
    nd_4859312 = nd_4858928.new_child(label=None, taxon=None, edge_length=0.006907, oid="Node169")
    nd_4859312.edge.oid = "Edge170"
    nd_4859312.edge.oid = "Edge170"
    nd_4859120 = nd_4859312.new_child(label=None, taxon=tax_4827056, edge_length=0.056894, oid="Node171")
    nd_4859600 = nd_4859312.new_child(label=None, taxon=None, edge_length=0.004552, oid="Node173")
    nd_4859600.edge.oid = "Edge174"
    nd_4859600.edge.oid = "Edge174"
    nd_4858960 = nd_4859600.new_child(label=None, taxon=None, edge_length=0.010752, oid="Node175")
    nd_4859952 = nd_4859600.new_child(label=None, taxon=tax_4826896, edge_length=0.079792, oid="Node185")
    nd_4859952.edge.oid = "Edge186"
    nd_4828272 = nd_4858960.new_child(label=None, taxon=None, edge_length=0.015909, oid="Node177")
    nd_4787664 = nd_4858960.new_child(label=None, taxon=tax_4828240, edge_length=0.114385, oid="Node183")
    nd_4787664.edge.oid = "Edge184"
    nd_4787728 = nd_4828272.new_child(label=None, taxon=tax_4827312, edge_length=0.092077, oid="Node179")
    nd_4828688 = nd_4828272.new_child(label=None, taxon=tax_4828432, edge_length=0.0579, oid="Node181")
    nd_4828688.edge.oid = "Edge182"
    nd_4828688.edge.oid = "Edge182"
    nd_4828688.edge.oid = "Edge182"
    nd_4828688.edge.oid = "Edge182"
    nd_4828688.edge.oid = "Edge182"
    nd_4828688.edge.oid = "Edge182"
    tree_4860272 = dendropy.Tree(label="PAUP 3", taxon_set=tree_list_4286432.taxon_set, oid="Tree189")
    tree_list_4286432.append(tree_4860272, reindex_taxa=False)
    nd_4860624 = tree_4860272.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.06111, oid="Node195")
    nd_4860752 = tree_4860272.seed_node.new_child(label=None, taxon=None, edge_length=0.005408, oid="Node197")
    nd_4883568 = tree_4860272.seed_node.new_child(label=None, taxon=None, edge_length=0.012121, oid="Node235")
    nd_4883568.edge.oid = "Edge236"
    nd_4883568.edge.oid = "Edge236"
    nd_4860336 = nd_4860752.new_child(label=None, taxon=None, edge_length=0.006988, oid="Node199")
    nd_4882800 = nd_4860752.new_child(label=None, taxon=tax_4827056, edge_length=0.063808, oid="Node233")
    nd_4882800.edge.oid = "Edge234"
    nd_4787696 = nd_4860336.new_child(label=None, taxon=None, edge_length=0.005167, oid="Node201")
    nd_4883056 = nd_4860336.new_child(label=None, taxon=None, edge_length=0.00574, oid="Node223")
    nd_4883056.edge.oid = "Edge224"
    nd_4860880 = nd_4787696.new_child(label=None, taxon=None, edge_length=0.010367, oid="Node203")
    nd_4861392 = nd_4787696.new_child(label=None, taxon=None, edge_length=0.007731, oid="Node213")
    nd_4861392.edge.oid = "Edge214"
    nd_4861008 = nd_4860880.new_child(label=None, taxon=None, edge_length=0.008688, oid="Node205")
    nd_4861648 = nd_4860880.new_child(label=None, taxon=tax_4828240, edge_length=0.1069, oid="Node211")
    nd_4861648.edge.oid = "Edge212"
    nd_4861136 = nd_4861008.new_child(label=None, taxon=tax_4825840, edge_length=0.067908, oid="Node207")
    nd_4861616 = nd_4861008.new_child(label=None, taxon=tax_4827952, edge_length=0.077196, oid="Node209")
    nd_4861616.edge.oid = "Edge210"
    nd_4861616.edge.oid = "Edge210"
    nd_4861616.edge.oid = "Edge210"
    nd_4861616.edge.oid = "Edge210"
    nd_4861840 = nd_4861392.new_child(label=None, taxon=None, edge_length=0.004129, oid="Node215")
    nd_4861712 = nd_4861392.new_child(label=None, taxon=tax_4827248, edge_length=0.068424, oid="Node221")
    nd_4861712.edge.oid = "Edge222"
    nd_4861904 = nd_4861840.new_child(label=None, taxon=tax_4827312, edge_length=0.109318, oid="Node217")
    nd_4861872 = nd_4861840.new_child(label=None, taxon=tax_4826896, edge_length=0.081382, oid="Node219")
    nd_4861872.edge.oid = "Edge220"
    nd_4861872.edge.oid = "Edge220"
    nd_4861872.edge.oid = "Edge220"
    nd_4861872.edge.oid = "Edge220"
    nd_4882928 = nd_4883056.new_child(label=None, taxon=tax_4788144, edge_length=0.079478, oid="Node225")
    nd_4883248 = nd_4883056.new_child(label=None, taxon=None, edge_length=0.013155, oid="Node227")
    nd_4883248.edge.oid = "Edge228"
    nd_4883248.edge.oid = "Edge228"
    nd_4882864 = nd_4883248.new_child(label=None, taxon=tax_4826768, edge_length=0.062239, oid="Node229")
    nd_4883504 = nd_4883248.new_child(label=None, taxon=tax_4827664, edge_length=0.060484, oid="Node231")
    nd_4883504.edge.oid = "Edge232"
    nd_4883504.edge.oid = "Edge232"
    nd_4883504.edge.oid = "Edge232"
    nd_4883504.edge.oid = "Edge232"
    nd_4883760 = nd_4883568.new_child(label=None, taxon=tax_4828432, edge_length=0.069501, oid="Node237")
    nd_4884048 = nd_4883568.new_child(label=None, taxon=tax_4826576, edge_length=0.073602, oid="Node239")
    nd_4884048.edge.oid = "Edge240"
    nd_4884048.edge.oid = "Edge240"
    nd_4884048.edge.oid = "Edge240"
    tree_4883824 = dendropy.Tree(label="PAUP 4", taxon_set=tree_list_4286432.taxon_set, oid="Tree241")
    tree_list_4286432.append(tree_4883824, reindex_taxa=False)
    nd_4884432 = tree_4883824.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.053028, oid="Node247")
    nd_4884560 = tree_4883824.seed_node.new_child(label=None, taxon=None, edge_length=0.014376, oid="Node249")
    nd_5030768 = tree_4883824.seed_node.new_child(label=None, taxon=tax_4826576, edge_length=0.076894, oid="Node291")
    nd_5030768.edge.oid = "Edge292"
    nd_5030768.edge.oid = "Edge292"
    nd_4884176 = nd_4884560.new_child(label=None, taxon=None, edge_length=0.002567, oid="Node251")
    nd_4884944 = nd_4884560.new_child(label=None, taxon=None, edge_length=0.005969, oid="Node261")
    nd_4884944.edge.oid = "Edge262"
    nd_4828720 = nd_4884176.new_child(label=None, taxon=None, edge_length=0.012413, oid="Node253")
    nd_4885200 = nd_4884176.new_child(label=None, taxon=tax_4827664, edge_length=0.069933, oid="Node259")
    nd_4885200.edge.oid = "Edge260"
    nd_4884688 = nd_4828720.new_child(label=None, taxon=tax_4825840, edge_length=0.072831, oid="Node255")
    nd_4885168 = nd_4828720.new_child(label=None, taxon=tax_4827952, edge_length=0.073004, oid="Node257")
    nd_4885168.edge.oid = "Edge258"
    nd_4885168.edge.oid = "Edge258"
    nd_4885168.edge.oid = "Edge258"
    nd_4885168.edge.oid = "Edge258"
    nd_4885392 = nd_4884944.new_child(label=None, taxon=None, edge_length=0.004803, oid="Node263")
    nd_5030544 = nd_4884944.new_child(label=None, taxon=tax_4827056, edge_length=0.064721, oid="Node289")
    nd_5030544.edge.oid = "Edge290"
    nd_4885456 = nd_4885392.new_child(label=None, taxon=None, edge_length=5.134e-05, oid="Node265")
    nd_5030512 = nd_4885392.new_child(label=None, taxon=tax_4828240, edge_length=0.118976, oid="Node287")
    nd_5030512.edge.oid = "Edge288"
    nd_4885424 = nd_4885456.new_child(label=None, taxon=None, edge_length=0.011519, oid="Node267")
    nd_4885968 = nd_4885456.new_child(label=None, taxon=None, edge_length=0.007647, oid="Node277")
    nd_4885968.edge.oid = "Edge278"
    nd_4885584 = nd_4885424.new_child(label=None, taxon=None, edge_length=0.011723, oid="Node269")
    nd_4886224 = nd_4885424.new_child(label=None, taxon=tax_4826768, edge_length=0.065297, oid="Node275")
    nd_4886224.edge.oid = "Edge276"
    nd_4885712 = nd_4885584.new_child(label=None, taxon=tax_4788144, edge_length=0.072339, oid="Node271")
    nd_4886160 = nd_4885584.new_child(label=None, taxon=tax_4827248, edge_length=0.062058, oid="Node273")
    nd_4886160.edge.oid = "Edge274"
    nd_4886160.edge.oid = "Edge274"
    nd_4886160.edge.oid = "Edge274"
    nd_4886160.edge.oid = "Edge274"
    nd_4886448 = nd_4885968.new_child(label=None, taxon=None, edge_length=0.022912, oid="Node279")
    nd_4886320 = nd_4885968.new_child(label=None, taxon=tax_4826896, edge_length=0.079945, oid="Node285")
    nd_4886320.edge.oid = "Edge286"
    nd_4886480 = nd_4886448.new_child(label=None, taxon=tax_4827312, edge_length=0.093762, oid="Node281")
    nd_4886416 = nd_4886448.new_child(label=None, taxon=tax_4828432, edge_length=0.055992, oid="Node283")
    nd_4886416.edge.oid = "Edge284"
    nd_4886416.edge.oid = "Edge284"
    nd_4886416.edge.oid = "Edge284"
    nd_4886416.edge.oid = "Edge284"
    nd_4886416.edge.oid = "Edge284"
    nd_4886416.edge.oid = "Edge284"
    nd_4886416.edge.oid = "Edge284"
    tree_5030640 = dendropy.Tree(label="PAUP 5", taxon_set=tree_list_4286432.taxon_set, oid="Tree293")
    tree_list_4286432.append(tree_5030640, reindex_taxa=False)
    nd_5031120 = tree_5030640.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.06391, oid="Node299")
    nd_5031248 = tree_5030640.seed_node.new_child(label=None, taxon=None, edge_length=0.00451, oid="Node301")
    nd_4884208 = tree_5030640.seed_node.new_child(label=None, taxon=None, edge_length=0.005404, oid="Node311")
    nd_4884208.edge.oid = "Edge312"
    nd_4884208.edge.oid = "Edge312"
    nd_5030352 = nd_5031248.new_child(label=None, taxon=tax_4825840, edge_length=0.071835, oid="Node303")
    nd_5031600 = nd_5031248.new_child(label=None, taxon=None, edge_length=0.010888, oid="Node305")
    nd_5031600.edge.oid = "Edge306"
    nd_5031600.edge.oid = "Edge306"
    nd_5031376 = nd_5031600.new_child(label=None, taxon=tax_4788144, edge_length=0.075788, oid="Node307")
    nd_5031824 = nd_5031600.new_child(label=None, taxon=tax_4828432, edge_length=0.07116, oid="Node309")
    nd_5031824.edge.oid = "Edge310"
    nd_5031824.edge.oid = "Edge310"
    nd_5031824.edge.oid = "Edge310"
    nd_4860464 = nd_4884208.new_child(label=None, taxon=None, edge_length=0.011944, oid="Node313")
    nd_5033424 = nd_4884208.new_child(label=None, taxon=None, edge_length=0.00612, oid="Node339")
    nd_5033424.edge.oid = "Edge340"
    nd_4884240 = nd_4860464.new_child(label=None, taxon=None, edge_length=0.013817, oid="Node315")
    nd_5032336 = nd_4860464.new_child(label=None, taxon=None, edge_length=0.012897, oid="Node325")
    nd_5032336.edge.oid = "Edge326"
    nd_4860368 = nd_4884240.new_child(label=None, taxon=None, edge_length=0.012356, oid="Node317")
    nd_4860432 = nd_4884240.new_child(label=None, taxon=tax_4826896, edge_length=0.077058, oid="Node323")
    nd_4860432.edge.oid = "Edge324"
    nd_5031568 = nd_4860368.new_child(label=None, taxon=tax_4827056, edge_length=0.048754, oid="Node319")
    nd_5032240 = nd_4860368.new_child(label=None, taxon=tax_4827952, edge_length=0.069682, oid="Node321")
    nd_5032240.edge.oid = "Edge322"
    nd_5032240.edge.oid = "Edge322"
    nd_5032240.edge.oid = "Edge322"
    nd_5032240.edge.oid = "Edge322"
    nd_5032464 = nd_5032336.new_child(label=None, taxon=None, edge_length=0.003147, oid="Node327")
    nd_5033296 = nd_5032336.new_child(label=None, taxon=tax_4827312, edge_length=0.108243, oid="Node337")
    nd_5033296.edge.oid = "Edge338"
    nd_5032528 = nd_5032464.new_child(label=None, taxon=tax_4826768, edge_length=0.07026, oid="Node329")
    nd_5032848 = nd_5032464.new_child(label=None, taxon=None, edge_length=0.001443, oid="Node331")
    nd_5032848.edge.oid = "Edge332"
    nd_5032848.edge.oid = "Edge332"
    nd_5032656 = nd_5032848.new_child(label=None, taxon=tax_4827664, edge_length=0.066182, oid="Node333")
    nd_5033136 = nd_5032848.new_child(label=None, taxon=tax_4827248, edge_length=0.068166, oid="Node335")
    nd_5033136.edge.oid = "Edge336"
    nd_5033136.edge.oid = "Edge336"
    nd_5033136.edge.oid = "Edge336"
    nd_5033136.edge.oid = "Edge336"
    nd_5032496 = nd_5033424.new_child(label=None, taxon=tax_4828240, edge_length=0.108072, oid="Node341")
    nd_5033648 = nd_5033424.new_child(label=None, taxon=tax_4826576, edge_length=0.073251, oid="Node343")
    nd_5033648.edge.oid = "Edge344"
    nd_5033648.edge.oid = "Edge344"
    nd_5033648.edge.oid = "Edge344"
    tree_5032912 = dendropy.Tree(label="PAUP 6", taxon_set=tree_list_4286432.taxon_set, oid="Tree345")
    tree_list_4286432.append(tree_5032912, reindex_taxa=False)
    nd_5058640 = tree_5032912.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.06501, oid="Node351")
    nd_5033232 = tree_5032912.seed_node.new_child(label=None, taxon=None, edge_length=0.003644, oid="Node353")
    nd_5060048 = tree_5032912.seed_node.new_child(label=None, taxon=None, edge_length=0.008523, oid="Node371")
    nd_5060048.edge.oid = "Edge372"
    nd_5060048.edge.oid = "Edge372"
    nd_5058832 = nd_5033232.new_child(label=None, taxon=None, edge_length=0.011722, oid="Node355")
    nd_5059152 = nd_5033232.new_child(label=None, taxon=None, edge_length=0.014836, oid="Node365")
    nd_5059152.edge.oid = "Edge366"
    nd_5030960 = nd_5058832.new_child(label=None, taxon=None, edge_length=0.017824, oid="Node357")
    nd_5059408 = nd_5058832.new_child(label=None, taxon=tax_4827312, edge_length=0.102639, oid="Node363")
    nd_5059408.edge.oid = "Edge364"
    nd_5058896 = nd_5030960.new_child(label=None, taxon=tax_4825840, edge_length=0.059984, oid="Node359")
    nd_5059376 = nd_5030960.new_child(label=None, taxon=tax_4828240, edge_length=0.100555, oid="Node361")
    nd_5059376.edge.oid = "Edge362"
    nd_5059376.edge.oid = "Edge362"
    nd_5059376.edge.oid = "Edge362"
    nd_5059376.edge.oid = "Edge362"
    nd_5059600 = nd_5059152.new_child(label=None, taxon=tax_4827664, edge_length=0.0538, oid="Node367")
    nd_5059888 = nd_5059152.new_child(label=None, taxon=tax_4826896, edge_length=0.07843, oid="Node369")
    nd_5059888.edge.oid = "Edge370"
    nd_5059888.edge.oid = "Edge370"
    nd_5059888.edge.oid = "Edge370"
    nd_5059472 = nd_5060048.new_child(label=None, taxon=None, edge_length=0.0, oid="Node373")
    nd_5061584 = nd_5060048.new_child(label=None, taxon=tax_4826768, edge_length=0.072091, oid="Node395")
    nd_5061584.edge.oid = "Edge396"
    nd_5059632 = nd_5059472.new_child(label=None, taxon=None, edge_length=0.015326, oid="Node375")
    nd_5060560 = nd_5059472.new_child(label=None, taxon=None, edge_length=0.005897, oid="Node381")
    nd_5060560.edge.oid = "Edge382"
    nd_5059920 = nd_5059632.new_child(label=None, taxon=tax_4788144, edge_length=0.068438, oid="Node377")
    nd_5060496 = nd_5059632.new_child(label=None, taxon=tax_4827248, edge_length=0.066815, oid="Node379")
    nd_5060496.edge.oid = "Edge380"
    nd_5060496.edge.oid = "Edge380"
    nd_5060496.edge.oid = "Edge380"
    nd_5060656 = nd_5060560.new_child(label=None, taxon=None, edge_length=0.010508, oid="Node383")
    nd_5061456 = nd_5060560.new_child(label=None, taxon=tax_4827952, edge_length=0.083409, oid="Node393")
    nd_5061456.edge.oid = "Edge394"
    nd_5060688 = nd_5060656.new_child(label=None, taxon=tax_4827056, edge_length=0.057007, oid="Node385")
    nd_5061040 = nd_5060656.new_child(label=None, taxon=None, edge_length=0.013392, oid="Node387")
    nd_5061040.edge.oid = "Edge388"
    nd_5061040.edge.oid = "Edge388"
    nd_5060816 = nd_5061040.new_child(label=None, taxon=tax_4828432, edge_length=0.062766, oid="Node389")
    nd_5061296 = nd_5061040.new_child(label=None, taxon=tax_4826576, edge_length=0.079145, oid="Node391")
    nd_5061296.edge.oid = "Edge392"
    nd_5061296.edge.oid = "Edge392"
    nd_5061296.edge.oid = "Edge392"
    nd_5061296.edge.oid = "Edge392"
    nd_5061296.edge.oid = "Edge392"
    tree_5060624 = dendropy.Tree(label="PAUP 7", taxon_set=tree_list_4286432.taxon_set, oid="Tree397")
    tree_list_4286432.append(tree_5060624, reindex_taxa=False)
    nd_5061936 = tree_5060624.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.062297, oid="Node403")
    nd_5062064 = tree_5060624.seed_node.new_child(label=None, taxon=None, edge_length=0.009899, oid="Node405")
    nd_4787440 = tree_5060624.seed_node.new_child(label=None, taxon=None, edge_length=0.003212, oid="Node411")
    nd_4787440.edge.oid = "Edge412"
    nd_4787440.edge.oid = "Edge412"
    nd_5061328 = nd_5062064.new_child(label=None, taxon=tax_4825840, edge_length=0.069165, oid="Node407")
    nd_5062416 = nd_5062064.new_child(label=None, taxon=tax_4788144, edge_length=0.073495, oid="Node409")
    nd_5062416.edge.oid = "Edge410"
    nd_5062416.edge.oid = "Edge410"
    nd_5062416.edge.oid = "Edge410"
    nd_5062448 = nd_4787440.new_child(label=None, taxon=None, edge_length=0.003672, oid="Node413")
    nd_5033840 = nd_4787440.new_child(label=None, taxon=tax_4827952, edge_length=0.080896, oid="Node447")
    nd_5033840.edge.oid = "Edge448"
    nd_5062576 = nd_5062448.new_child(label=None, taxon=None, edge_length=0.010249, oid="Node415")
    nd_5088592 = nd_5062448.new_child(label=None, taxon=tax_4828240, edge_length=0.114381, oid="Node445")
    nd_5088592.edge.oid = "Edge446"
    nd_5062512 = nd_5062576.new_child(label=None, taxon=None, edge_length=0.00713, oid="Node417")
    nd_5089040 = nd_5062576.new_child(label=None, taxon=tax_4828432, edge_length=0.079534, oid="Node443")
    nd_5089040.edge.oid = "Edge444"
    nd_5087312 = nd_5062512.new_child(label=None, taxon=None, edge_length=0.000386, oid="Node419")
    nd_5089104 = nd_5062512.new_child(label=None, taxon=tax_4827664, edge_length=0.070614, oid="Node441")
    nd_5089104.edge.oid = "Edge442"
    nd_5087440 = nd_5087312.new_child(label=None, taxon=None, edge_length=0.003603, oid="Node421")
    nd_5062544 = nd_5087312.new_child(label=None, taxon=None, edge_length=0.017571, oid="Node435")
    nd_5062544.edge.oid = "Edge436"
    nd_5087568 = nd_5087440.new_child(label=None, taxon=tax_4827056, edge_length=0.061022, oid="Node423")
    nd_5088048 = nd_5087440.new_child(label=None, taxon=None, edge_length=0.004698, oid="Node425")
    nd_5088048.edge.oid = "Edge426"
    nd_5088048.edge.oid = "Edge426"
    nd_5087824 = nd_5088048.new_child(label=None, taxon=tax_4826768, edge_length=0.071291, oid="Node427")
    nd_5088272 = nd_5088048.new_child(label=None, taxon=None, edge_length=0.017245, oid="Node429")
    nd_5088272.edge.oid = "Edge430"
    nd_5088272.edge.oid = "Edge430"
    nd_5087696 = nd_5088272.new_child(label=None, taxon=tax_4827312, edge_length=0.098776, oid="Node431")
    nd_5088560 = nd_5088272.new_child(label=None, taxon=tax_4826576, edge_length=0.075069, oid="Node433")
    nd_5088560.edge.oid = "Edge434"
    nd_5088560.edge.oid = "Edge434"
    nd_5088560.edge.oid = "Edge434"
    nd_5088688 = nd_5062544.new_child(label=None, taxon=tax_4826896, edge_length=0.075494, oid="Node437")
    nd_5088944 = nd_5062544.new_child(label=None, taxon=tax_4827248, edge_length=0.063757, oid="Node439")
    nd_5088944.edge.oid = "Edge440"
    nd_5088944.edge.oid = "Edge440"
    nd_5088944.edge.oid = "Edge440"
    nd_5088944.edge.oid = "Edge440"
    nd_5088944.edge.oid = "Edge440"
    nd_5088944.edge.oid = "Edge440"
    nd_5088944.edge.oid = "Edge440"
    tree_5033392 = dendropy.Tree(label="PAUP 8", taxon_set=tree_list_4286432.taxon_set, oid="Tree449")
    tree_list_4286432.append(tree_5033392, reindex_taxa=False)
    nd_5089584 = tree_5033392.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.05678, oid="Node455")
    nd_5033872 = tree_5033392.seed_node.new_child(label=None, taxon=None, edge_length=0.014479, oid="Node457")
    nd_5116976 = tree_5033392.seed_node.new_child(label=None, taxon=tax_4827664, edge_length=0.060374, oid="Node499")
    nd_5116976.edge.oid = "Edge500"
    nd_5116976.edge.oid = "Edge500"
    nd_5089168 = nd_5033872.new_child(label=None, taxon=None, edge_length=0.00182, oid="Node459")
    nd_5117008 = nd_5033872.new_child(label=None, taxon=tax_4788144, edge_length=0.079156, oid="Node497")
    nd_5117008.edge.oid = "Edge498"
    nd_4884272 = nd_5089168.new_child(label=None, taxon=None, edge_length=0.002833, oid="Node461")
    nd_5116176 = nd_5089168.new_child(label=None, taxon=None, edge_length=0.010734, oid="Node491")
    nd_5116176.edge.oid = "Edge492"
    nd_5089840 = nd_4884272.new_child(label=None, taxon=None, edge_length=0.012314, oid="Node463")
    nd_5116304 = nd_4884272.new_child(label=None, taxon=tax_4826576, edge_length=0.081982, oid="Node489")
    nd_5116304.edge.oid = "Edge490"
    nd_5089968 = nd_5089840.new_child(label=None, taxon=None, edge_length=0.001866, oid="Node465")
    nd_5091056 = nd_5089840.new_child(label=None, taxon=tax_4827056, edge_length=0.06129, oid="Node487")
    nd_5091056.edge.oid = "Edge488"
    nd_5090096 = nd_5089968.new_child(label=None, taxon=None, edge_length=0.005104, oid="Node467")
    nd_5116240 = nd_5089968.new_child(label=None, taxon=tax_4828432, edge_length=0.07646, oid="Node485")
    nd_5116240.edge.oid = "Edge486"
    nd_5090224 = nd_5090096.new_child(label=None, taxon=None, edge_length=0.004377, oid="Node469")
    nd_5090736 = nd_5090096.new_child(label=None, taxon=None, edge_length=0.016433, oid="Node479")
    nd_5090736.edge.oid = "Edge480"
    nd_5090352 = nd_5090224.new_child(label=None, taxon=None, edge_length=0.013813, oid="Node471")
    nd_5090992 = nd_5090224.new_child(label=None, taxon=tax_4827952, edge_length=0.083335, oid="Node477")
    nd_5090992.edge.oid = "Edge478"
    nd_5090480 = nd_5090352.new_child(label=None, taxon=tax_4825840, edge_length=0.065928, oid="Node473")
    nd_5090960 = nd_5090352.new_child(label=None, taxon=tax_4827312, edge_length=0.10077, oid="Node475")
    nd_5090960.edge.oid = "Edge476"
    nd_5090960.edge.oid = "Edge476"
    nd_5090960.edge.oid = "Edge476"
    nd_5090960.edge.oid = "Edge476"
    nd_5091184 = nd_5090736.new_child(label=None, taxon=tax_4826896, edge_length=0.073392, oid="Node481")
    nd_5091216 = nd_5090736.new_child(label=None, taxon=tax_4827248, edge_length=0.065748, oid="Node483")
    nd_5091216.edge.oid = "Edge484"
    nd_5091216.edge.oid = "Edge484"
    nd_5091216.edge.oid = "Edge484"
    nd_5091216.edge.oid = "Edge484"
    nd_5091216.edge.oid = "Edge484"
    nd_5091216.edge.oid = "Edge484"
    nd_5116592 = nd_5116176.new_child(label=None, taxon=tax_4826768, edge_length=0.064242, oid="Node493")
    nd_5116816 = nd_5116176.new_child(label=None, taxon=tax_4828240, edge_length=0.109434, oid="Node495")
    nd_5116816.edge.oid = "Edge496"
    nd_5116816.edge.oid = "Edge496"
    nd_5116816.edge.oid = "Edge496"
    nd_5116816.edge.oid = "Edge496"
    nd_5116816.edge.oid = "Edge496"
    tree_5117072 = dendropy.Tree(label="PAUP 9", taxon_set=tree_list_4286432.taxon_set, oid="Tree501")
    tree_list_4286432.append(tree_5117072, reindex_taxa=False)
    nd_5117488 = tree_5117072.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.063004, oid="Node507")
    nd_5117616 = tree_5117072.seed_node.new_child(label=None, taxon=None, edge_length=0.018059, oid="Node509")
    nd_5061808 = tree_5117072.seed_node.new_child(label=None, taxon=None, edge_length=0.00238, oid="Node515")
    nd_5061808.edge.oid = "Edge516"
    nd_5061808.edge.oid = "Edge516"
    nd_5117136 = nd_5117616.new_child(label=None, taxon=tax_4825840, edge_length=0.061413, oid="Node511")
    nd_5117968 = nd_5117616.new_child(label=None, taxon=tax_4828240, edge_length=0.099507, oid="Node513")
    nd_5117968.edge.oid = "Edge514"
    nd_5117968.edge.oid = "Edge514"
    nd_5117968.edge.oid = "Edge514"
    nd_5118000 = nd_5061808.new_child(label=None, taxon=None, edge_length=0.003714, oid="Node517")
    nd_5119888 = nd_5061808.new_child(label=None, taxon=None, edge_length=0.010478, oid="Node547")
    nd_5119888.edge.oid = "Edge548"
    nd_5118128 = nd_5118000.new_child(label=None, taxon=None, edge_length=0.016533, oid="Node519")
    nd_5119792 = nd_5118000.new_child(label=None, taxon=tax_4827248, edge_length=0.068818, oid="Node545")
    nd_5119792.edge.oid = "Edge546"
    nd_5118096 = nd_5118128.new_child(label=None, taxon=None, edge_length=0.010763, oid="Node521")
    nd_5119280 = nd_5118128.new_child(label=None, taxon=tax_4828432, edge_length=0.070762, oid="Node543")
    nd_5119280.edge.oid = "Edge544"
    nd_5118256 = nd_5118096.new_child(label=None, taxon=None, edge_length=0.002182, oid="Node523")
    nd_5118768 = nd_5118096.new_child(label=None, taxon=None, edge_length=0.011197, oid="Node533")
    nd_5118768.edge.oid = "Edge534"
    nd_5118384 = nd_5118256.new_child(label=None, taxon=None, edge_length=0.014459, oid="Node525")
    nd_5119024 = nd_5118256.new_child(label=None, taxon=tax_4827056, edge_length=0.054698, oid="Node531")
    nd_5119024.edge.oid = "Edge532"
    nd_5118512 = nd_5118384.new_child(label=None, taxon=tax_4788144, edge_length=0.075877, oid="Node527")
    nd_5118960 = nd_5118384.new_child(label=None, taxon=tax_4826896, edge_length=0.071375, oid="Node529")
    nd_5118960.edge.oid = "Edge530"
    nd_5118960.edge.oid = "Edge530"
    nd_5118960.edge.oid = "Edge530"
    nd_5118960.edge.oid = "Edge530"
    nd_5119216 = nd_5118768.new_child(label=None, taxon=tax_4827664, edge_length=0.062291, oid="Node535")
    nd_5119504 = nd_5118768.new_child(label=None, taxon=None, edge_length=0.00521, oid="Node537")
    nd_5119504.edge.oid = "Edge538"
    nd_5119504.edge.oid = "Edge538"
    nd_5119248 = nd_5119504.new_child(label=None, taxon=tax_4827952, edge_length=0.074463, oid="Node539")
    nd_5119760 = nd_5119504.new_child(label=None, taxon=tax_4827312, edge_length=0.103887, oid="Node541")
    nd_5119760.edge.oid = "Edge542"
    nd_5119760.edge.oid = "Edge542"
    nd_5119760.edge.oid = "Edge542"
    nd_5119760.edge.oid = "Edge542"
    nd_5119760.edge.oid = "Edge542"
    nd_5144752 = nd_5119888.new_child(label=None, taxon=tax_4826768, edge_length=0.071307, oid="Node549")
    nd_5144976 = nd_5119888.new_child(label=None, taxon=tax_4826576, edge_length=0.073392, oid="Node551")
    nd_5144976.edge.oid = "Edge552"
    nd_5144976.edge.oid = "Edge552"
    nd_5144976.edge.oid = "Edge552"
    tree_5144656 = dendropy.Tree(label="PAUP 10", taxon_set=tree_list_4286432.taxon_set, oid="Tree553")
    tree_list_4286432.append(tree_5144656, reindex_taxa=False)
    nd_5145392 = tree_5144656.seed_node.new_child(label=None, taxon=tax_4787280, edge_length=0.053901, oid="Node559")
    nd_5089424 = tree_5144656.seed_node.new_child(label=None, taxon=None, edge_length=0.013947, oid="Node561")
    nd_5147312 = tree_5144656.seed_node.new_child(label=None, taxon=tax_4828240, edge_length=0.107868, oid="Node603")
    nd_5147312.edge.oid = "Edge604"
    nd_5147312.edge.oid = "Edge604"
    nd_5145104 = nd_5089424.new_child(label=None, taxon=None, edge_length=0.002034, oid="Node563")
    nd_5147664 = nd_5089424.new_child(label=None, taxon=tax_4827056, edge_length=0.064652, oid="Node601")
    nd_5147664.edge.oid = "Edge602"
    nd_5061776 = nd_5145104.new_child(label=None, taxon=None, edge_length=0.00777, oid="Node565")
    nd_5146928 = nd_5145104.new_child(label=None, taxon=None, edge_length=0.00373, oid="Node587")
    nd_5146928.edge.oid = "Edge588"
    nd_5145648 = nd_5061776.new_child(label=None, taxon=None, edge_length=0.003041, oid="Node567")
    nd_5146672 = nd_5061776.new_child(label=None, taxon=None, edge_length=0.008083, oid="Node577")
    nd_5146672.edge.oid = "Edge578"
    nd_5145776 = nd_5145648.new_child(label=None, taxon=tax_4825840, edge_length=0.07927, oid="Node569")
    nd_5146256 = nd_5145648.new_child(label=None, taxon=None, edge_length=0.016769, oid="Node571")
    nd_5146256.edge.oid = "Edge572"
    nd_5146256.edge.oid = "Edge572"
    nd_5146032 = nd_5146256.new_child(label=None, taxon=tax_4827664, edge_length=0.053956, oid="Node573")
    nd_5146512 = nd_5146256.new_child(label=None, taxon=tax_4827312, edge_length=0.098874, oid="Node575")
    nd_5146512.edge.oid = "Edge576"
    nd_5146512.edge.oid = "Edge576"
    nd_5146512.edge.oid = "Edge576"
    nd_5146288 = nd_5146672.new_child(label=None, taxon=None, edge_length=0.010712, oid="Node579")
    nd_5117360 = nd_5146672.new_child(label=None, taxon=tax_4828432, edge_length=0.075934, oid="Node585")
    nd_5117360.edge.oid = "Edge586"
    nd_5089392 = nd_5146288.new_child(label=None, taxon=tax_4788144, edge_length=0.075223, oid="Node581")
    nd_5089232 = nd_5146288.new_child(label=None, taxon=tax_4826768, edge_length=0.065805, oid="Node583")
    nd_5089232.edge.oid = "Edge584"
    nd_5089232.edge.oid = "Edge584"
    nd_5089232.edge.oid = "Edge584"
    nd_5089232.edge.oid = "Edge584"
    nd_5146800 = nd_5146928.new_child(label=None, taxon=tax_4827952, edge_length=0.082425, oid="Node589")
    nd_5147152 = nd_5146928.new_child(label=None, taxon=None, edge_length=0.008508, oid="Node591")
    nd_5147152.edge.oid = "Edge592"
    nd_5147152.edge.oid = "Edge592"
    nd_5146768 = nd_5147152.new_child(label=None, taxon=None, edge_length=0.010583, oid="Node593")
    nd_5147568 = nd_5147152.new_child(label=None, taxon=tax_4827248, edge_length=0.070289, oid="Node599")
    nd_5147568.edge.oid = "Edge600"
    nd_5146896 = nd_5146768.new_child(label=None, taxon=tax_4826896, edge_length=0.075429, oid="Node595")
    nd_5147536 = nd_5146768.new_child(label=None, taxon=tax_4826576, edge_length=0.0787, oid="Node597")
    nd_5147536.edge.oid = "Edge598"
    nd_5147536.edge.oid = "Edge598"
    nd_5147536.edge.oid = "Edge598"
    nd_5147536.edge.oid = "Edge598"
    nd_5147536.edge.oid = "Edge598"
    nd_5147536.edge.oid = "Edge598"
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
