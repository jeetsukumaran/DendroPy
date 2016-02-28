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
Tests of tree metrics.
"""

import random
import math
import unittest
from dendropy.test.support import dendropytest
from dendropy.test.support import pathmap

import dendropy
from dendropy.calculate import treemeasure
from dendropy.calculate import treecompare
from dendropy.utility.textprocessing import StringIO

def _get_reference_tree_list(taxon_namespace=None):
    tree_list = dendropy.TreeList(label=None, taxon_namespace=taxon_namespace)
    tax_4313741136 = tree_list.taxon_namespace.require_taxon(label="Antaresia childreni")
    tax_4313741328 = tree_list.taxon_namespace.require_taxon(label="Antaresia maculosa")
    tax_4313741456 = tree_list.taxon_namespace.require_taxon(label="Antaresia melanocephalus")
    tax_4313741584 = tree_list.taxon_namespace.require_taxon(label="Antaresia perthensis")
    tax_4313741712 = tree_list.taxon_namespace.require_taxon(label="Antaresia ramsayi")
    tax_4313741840 = tree_list.taxon_namespace.require_taxon(label="Antaresia stimsoni")
    tax_4313741904 = tree_list.taxon_namespace.require_taxon(label="Apodora papuana")
    tax_4313741968 = tree_list.taxon_namespace.require_taxon(label="Bothrochilus boa")
    tax_4313742032 = tree_list.taxon_namespace.require_taxon(label="Candoia aspera")
    tax_4313742160 = tree_list.taxon_namespace.require_taxon(label="Liasis albertisii")
    tax_4313742224 = tree_list.taxon_namespace.require_taxon(label="Liasis fuscus")
    tax_4313742288 = tree_list.taxon_namespace.require_taxon(label="Liasis mackloti")
    tax_4313742352 = tree_list.taxon_namespace.require_taxon(label="Liasis olivaceus")
    tax_4313742480 = tree_list.taxon_namespace.require_taxon(label="Loxocemus bicolor")
    tax_4313742608 = tree_list.taxon_namespace.require_taxon(label="Morelia amethistina")
    tax_4313742672 = tree_list.taxon_namespace.require_taxon(label="Morelia boeleni")
    tax_4313742736 = tree_list.taxon_namespace.require_taxon(label="Morelia bredli")
    tax_4313742800 = tree_list.taxon_namespace.require_taxon(label="Morelia carinata")
    tax_4313742928 = tree_list.taxon_namespace.require_taxon(label="Morelia clastolepis")
    tax_4313743056 = tree_list.taxon_namespace.require_taxon(label="Morelia kinghorni")
    tax_4313743120 = tree_list.taxon_namespace.require_taxon(label="Morelia nauta")
    tax_4313743248 = tree_list.taxon_namespace.require_taxon(label="Morelia oenpelliensis")
    tax_4313743312 = tree_list.taxon_namespace.require_taxon(label="Morelia spilota")
    tax_4313759824 = tree_list.taxon_namespace.require_taxon(label="Morelia tracyae")
    tax_4313759888 = tree_list.taxon_namespace.require_taxon(label="Morelia viridisN")
    tax_4313759952 = tree_list.taxon_namespace.require_taxon(label="Morelia viridisS")
    tax_4313760016 = tree_list.taxon_namespace.require_taxon(label="Python curtus")
    tax_4313760080 = tree_list.taxon_namespace.require_taxon(label="Python molurus")
    tax_4313760144 = tree_list.taxon_namespace.require_taxon(label="Python regius")
    tax_4313760272 = tree_list.taxon_namespace.require_taxon(label="Python reticulatus")
    tax_4313760336 = tree_list.taxon_namespace.require_taxon(label="Python sebae")
    tax_4313760464 = tree_list.taxon_namespace.require_taxon(label="Python timoriensis")
    tax_4313760592 = tree_list.taxon_namespace.require_taxon(label="Xenopeltis unicolor")
    tree_4313760848 = dendropy.Tree(label="Tree01", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4313760848, reindex_taxa=False)
    nd_4313761232 = tree_4313760848.seed_node.new_child(label="Node4313761232", taxon=None, edge_length=78.6266408419)
    nd_4313801936 = tree_4313760848.seed_node.new_child(label="Node4313801936", taxon=None, edge_length=229.880308935)
    nd_4313761360 = nd_4313761232.new_child(label="Node4313761360", taxon=None, edge_length=123.936295332)
    nd_4313801552 = nd_4313761232.new_child(label="Node4313801552", taxon=None, edge_length=202.422980755)
    nd_4313761488 = nd_4313761360.new_child(label="Node4313761488", taxon=None, edge_length=53.7887032791)
    nd_4313779728 = nd_4313761360.new_child(label="Node4313779728", taxon=None, edge_length=69.6091072013)
    nd_4313761616 = nd_4313761488.new_child(label="Node4313761616", taxon=None, edge_length=4.50960412336)
    nd_4313779088 = nd_4313761488.new_child(label="Node4313779088", taxon=None, edge_length=24.4333103788)
    nd_4313761744 = nd_4313761616.new_child(label="Node4313761744", taxon=None, edge_length=1.29166787247)
    nd_4313778704 = nd_4313761616.new_child(label="Node4313778704", taxon=None, edge_length=13.1562713928)
    nd_4313761872 = nd_4313761744.new_child(label="Node4313761872", taxon=None, edge_length=0.739270951321)
    nd_4313777808 = nd_4313761744.new_child(label="Node4313777808", taxon=None, edge_length=16.8726126715)
    nd_4313762000 = nd_4313761872.new_child(label="Node4313762000", taxon=None, edge_length=4.00932181352)
    nd_4313777552 = nd_4313761872.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=18.3098649933)
    nd_4313762128 = nd_4313762000.new_child(label="Node4313762128", taxon=None, edge_length=10.4795105126)
    nd_4313762704 = nd_4313762000.new_child(label="Node4313762704", taxon=None, edge_length=3.67431549248)
    nd_4313762256 = nd_4313762128.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=3.82103266723)
    nd_4313762448 = nd_4313762128.new_child(label="Node4313762448", taxon=None, edge_length=2.72293867732)
    nd_4313762576 = nd_4313762448.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=1.09809398991)
    nd_4313762384 = nd_4313762448.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=1.09809398991)
    nd_4313762960 = nd_4313762704.new_child(label="Node4313762960", taxon=None, edge_length=4.47413907662)
    nd_4313776720 = nd_4313762704.new_child(label="Node4313776720", taxon=None, edge_length=3.27824986702)
    nd_4313763088 = nd_4313762960.new_child(label="Node4313763088", taxon=None, edge_length=1.10594701364)
    nd_4313776208 = nd_4313762960.new_child(label="Node4313776208", taxon=None, edge_length=0.946531890581)
    nd_4313763216 = nd_4313763088.new_child(label="Node4313763216", taxon=None, edge_length=3.56964311758)
    nd_4313763600 = nd_4313763088.new_child(label="Node4313763600", taxon=None, edge_length=3.54960462938)
    nd_4313763344 = nd_4313763216.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=1.47649847948)
    nd_4313762832 = nd_4313763216.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=1.47649847948)
    nd_4313763728 = nd_4313763600.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=1.49653696768)
    nd_4313763472 = nd_4313763600.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=1.49653696768)
    nd_4313776464 = nd_4313776208.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=5.20555672012)
    nd_4313776336 = nd_4313776208.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=5.20555672012)
    nd_4313776912 = nd_4313776720.new_child(label="Node4313776912", taxon=None, edge_length=2.88197852526)
    nd_4313777296 = nd_4313776720.new_child(label="Node4313777296", taxon=None, edge_length=6.86415378064)
    nd_4313777040 = nd_4313776912.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=4.46599929505)
    nd_4313776784 = nd_4313776912.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=4.46599929505)
    nd_4313777424 = nd_4313777296.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=0.483824039664)
    nd_4313777168 = nd_4313777296.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=0.483824039664)
    nd_4313777936 = nd_4313777808.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=2.17652327312)
    nd_4313777680 = nd_4313777808.new_child(label="Node4313777680", taxon=None, edge_length=1.67230791531)
    nd_4313778192 = nd_4313777680.new_child(label="Node4313778192", taxon=None, edge_length=0.491713738136)
    nd_4313778576 = nd_4313777680.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=0.504215357803)
    nd_4313778320 = nd_4313778192.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=0.0125016196671)
    nd_4313778064 = nd_4313778192.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=0.0125016196671)
    nd_4313778832 = nd_4313778704.new_child(label="Node4313778832", taxon=None, edge_length=2.98661848623)
    nd_4313779216 = nd_4313778704.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=7.18453242432)
    nd_4313778960 = nd_4313778832.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=4.19791393809)
    nd_4313778448 = nd_4313778832.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=4.19791393809)
    nd_4313779472 = nd_4313779088.new_child(label="Node4313779472", taxon=None, edge_length=0.207889736001)
    nd_4313779856 = nd_4313779088.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=0.417097561686)
    nd_4313779600 = nd_4313779472.new_child(label="Python regius", taxon=tax_4313760144, edge_length=0.209207825685)
    nd_4313779344 = nd_4313779472.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=0.209207825685)
    nd_4313780112 = nd_4313779728.new_child(label="Node4313780112", taxon=None, edge_length=1.24643505521)
    nd_4313801424 = nd_4313779728.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=9.03000401821)
    nd_4313800784 = nd_4313780112.new_child(label="Node4313800784", taxon=None, edge_length=2.6364754365)
    nd_4313801168 = nd_4313780112.new_child(label="Node4313801168", taxon=None, edge_length=7.12573328141)
    nd_4313800912 = nd_4313800784.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=5.1470935265)
    nd_4313779984 = nd_4313800784.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=5.1470935265)
    nd_4313801296 = nd_4313801168.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=0.657835681585)
    nd_4313801104 = nd_4313801168.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=0.657835681585)
    nd_4313801808 = nd_4313801552.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=0.152425796315)
    nd_4313801680 = nd_4313801552.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=0.152425796315)
    nd_4313802192 = nd_4313801936.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=51.3217384582)
    nd_4313802064 = nd_4313801936.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=51.3217384582)
    tree_4313802320 = dendropy.Tree(label="Tree02", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4313802320, reindex_taxa=False)
    nd_4313802704 = tree_4313802320.seed_node.new_child(label="Node4313802704", taxon=None, edge_length=18.8917197007)
    nd_4316157840 = tree_4313802320.seed_node.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=100.189141925)
    nd_4313802832 = nd_4313802704.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=81.2974222246)
    nd_4313803024 = nd_4313802704.new_child(label="Node4313803024", taxon=None, edge_length=33.1565984398)
    nd_4313803152 = nd_4313803024.new_child(label="Node4313803152", taxon=None, edge_length=6.57324583185)
    nd_4316156752 = nd_4313803024.new_child(label="Node4316156752", taxon=None, edge_length=16.2594519516)
    nd_4313803280 = nd_4313803152.new_child(label="Node4313803280", taxon=None, edge_length=0.76583222117)
    nd_4313829328 = nd_4313803152.new_child(label="Node4313829328", taxon=None, edge_length=2.53377266123)
    nd_4313803408 = nd_4313803280.new_child(label="Node4313803408", taxon=None, edge_length=4.5936111676)
    nd_4313826384 = nd_4313803280.new_child(label="Node4313826384", taxon=None, edge_length=1.27553605821)
    nd_4313803536 = nd_4313803408.new_child(label="Node4313803536", taxon=None, edge_length=15.1180336863)
    nd_4313804432 = nd_4313803408.new_child(label="Node4313804432", taxon=None, edge_length=3.40951166184)
    nd_4313803664 = nd_4313803536.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=21.0901008779)
    nd_4313802960 = nd_4313803536.new_child(label="Node4313802960", taxon=None, edge_length=3.38653541663)
    nd_4313803920 = nd_4313802960.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=17.7035654613)
    nd_4313803792 = nd_4313802960.new_child(label="Node4313803792", taxon=None, edge_length=2.6244717729)
    nd_4313804112 = nd_4313803792.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=15.0790936884)
    nd_4313804304 = nd_4313803792.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=15.0790936884)
    nd_4313804560 = nd_4313804432.new_child(label="Node4313804560", taxon=None, edge_length=5.00375936144)
    nd_4313825616 = nd_4313804432.new_child(label="Node4313825616", taxon=None, edge_length=5.81119736053)
    nd_4313804688 = nd_4313804560.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=27.7948635409)
    nd_4313804240 = nd_4313804560.new_child(label="Node4313804240", taxon=None, edge_length=8.32618746237)
    nd_4313825488 = nd_4313804240.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=19.4686760786)
    nd_4313825360 = nd_4313804240.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=19.4686760786)
    nd_4313825872 = nd_4313825616.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=26.9874255418)
    nd_4313825744 = nd_4313825616.new_child(label="Node4313825744", taxon=None, edge_length=2.25683638168)
    nd_4313826128 = nd_4313825744.new_child(label="Node4313826128", taxon=None, edge_length=16.6530983052)
    nd_4313826512 = nd_4313825744.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=24.7305891602)
    nd_4313826256 = nd_4313826128.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=8.07749085501)
    nd_4313826000 = nd_4313826128.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=8.07749085501)
    nd_4313826768 = nd_4313826384.new_child(label="Node4313826768", taxon=None, edge_length=4.33214615343)
    nd_4313828048 = nd_4313826384.new_child(label="Node4313828048", taxon=None, edge_length=10.9652932592)
    nd_4313826896 = nd_4313826768.new_child(label="Node4313826896", taxon=None, edge_length=3.37363071467)
    nd_4313827664 = nd_4313826768.new_child(label="Node4313827664", taxon=None, edge_length=16.3762764593)
    nd_4313827024 = nd_4313826896.new_child(label="Node4313827024", taxon=None, edge_length=26.1365403684)
    nd_4313827408 = nd_4313826896.new_child(label="Node4313827408", taxon=None, edge_length=12.1064068345)
    nd_4313827152 = nd_4313827024.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=5.68389243709)
    nd_4313826640 = nd_4313827024.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=5.68389243709)
    nd_4313827536 = nd_4313827408.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=19.714025971)
    nd_4313827280 = nd_4313827408.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=19.714025971)
    nd_4313827920 = nd_4313827664.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=18.8177870609)
    nd_4313827792 = nd_4313827664.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=18.8177870609)
    nd_4313828304 = nd_4313828048.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=28.5609164144)
    nd_4313828176 = nd_4313828048.new_child(label="Node4313828176", taxon=None, edge_length=11.3916491298)
    nd_4313828560 = nd_4313828176.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=17.1692672846)
    nd_4313828432 = nd_4313828176.new_child(label="Node4313828432", taxon=None, edge_length=10.4784522084)
    nd_4313828816 = nd_4313828432.new_child(label="Node4313828816", taxon=None, edge_length=1.44575855569)
    nd_4313829200 = nd_4313828432.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=6.69081507619)
    nd_4313828944 = nd_4313828816.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=5.2450565205)
    nd_4313828688 = nd_4313828816.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=5.2450565205)
    nd_4316155984 = nd_4313829328.new_child(label="Node4316155984", taxon=None, edge_length=23.8923728386)
    nd_4313802512 = nd_4313829328.new_child(label="Node4313802512", taxon=None, edge_length=22.5856079922)
    nd_4316156112 = nd_4316155984.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=15.1414324532)
    nd_4316156368 = nd_4316155984.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=15.1414324532)
    nd_4316156496 = nd_4313802512.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=16.4481972995)
    nd_4313802448 = nd_4313802512.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=16.4481972995)
    nd_4316156880 = nd_4316156752.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=31.8813718333)
    nd_4316156624 = nd_4316156752.new_child(label="Node4316156624", taxon=None, edge_length=5.97313984611)
    nd_4316157136 = nd_4316156624.new_child(label="Node4316157136", taxon=None, edge_length=8.94343133576)
    nd_4316157904 = nd_4316156624.new_child(label="Python regius", taxon=tax_4313760144, edge_length=25.9082319872)
    nd_4316157264 = nd_4316157136.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=16.9648006514)
    nd_4316157008 = nd_4316157136.new_child(label="Node4316157008", taxon=None, edge_length=4.66373979181)
    nd_4316157584 = nd_4316157008.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=12.3010608596)
    nd_4316157456 = nd_4316157008.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=12.3010608596)
    tree_4316158160 = dendropy.Tree(label="Tree03", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316158160, reindex_taxa=False)
    nd_4316158416 = tree_4316158160.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=109.372663833)
    nd_4316158608 = tree_4316158160.seed_node.new_child(label="Node4316158608", taxon=None, edge_length=23.0231215792)
    nd_4316158736 = nd_4316158608.new_child(label="Node4316158736", taxon=None, edge_length=25.5450384687)
    nd_4316203024 = nd_4316158608.new_child(label="Node4316203024", taxon=None, edge_length=3.40243621372)
    nd_4316158864 = nd_4316158736.new_child(label="Node4316158864", taxon=None, edge_length=8.2214192007)
    nd_4316202128 = nd_4316158736.new_child(label="Node4316202128", taxon=None, edge_length=17.448321597)
    nd_4316158992 = nd_4316158864.new_child(label="Node4316158992", taxon=None, edge_length=6.47536868825)
    nd_4316201552 = nd_4316158864.new_child(label="Node4316201552", taxon=None, edge_length=22.2791730332)
    nd_4316159120 = nd_4316158992.new_child(label="Node4316159120", taxon=None, edge_length=2.53965409687)
    nd_4316179792 = nd_4316158992.new_child(label="Node4316179792", taxon=None, edge_length=7.82653683343)
    nd_4316159248 = nd_4316159120.new_child(label="Node4316159248", taxon=None, edge_length=2.41802497137)
    nd_4316178128 = nd_4316159120.new_child(label="Node4316178128", taxon=None, edge_length=5.90175129715)
    nd_4316159376 = nd_4316159248.new_child(label="Node4316159376", taxon=None, edge_length=6.39175712039)
    nd_4316177104 = nd_4316159248.new_child(label="Node4316177104", taxon=None, edge_length=9.11329596086)
    nd_4316159504 = nd_4316159376.new_child(label="Node4316159504", taxon=None, edge_length=4.82772953939)
    nd_4316176720 = nd_4316159376.new_child(label="Node4316176720", taxon=None, edge_length=17.7955396972)
    nd_4316159632 = nd_4316159504.new_child(label="Node4316159632", taxon=None, edge_length=18.9362130146)
    nd_4316176464 = nd_4316159504.new_child(label="Node4316176464", taxon=None, edge_length=8.55401069191)
    nd_4316159760 = nd_4316159632.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=10.9943371535)
    nd_4316158544 = nd_4316159632.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=10.9943371535)
    nd_4316176592 = nd_4316176464.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=21.3765394762)
    nd_4316159888 = nd_4316176464.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=21.3765394762)
    nd_4316176976 = nd_4316176720.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=16.9627400103)
    nd_4316176848 = nd_4316176720.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=16.9627400103)
    nd_4316177360 = nd_4316177104.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=32.036740867)
    nd_4316177232 = nd_4316177104.new_child(label="Node4316177232", taxon=None, edge_length=14.7791918926)
    nd_4316177616 = nd_4316177232.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=17.2575489744)
    nd_4316177488 = nd_4316177232.new_child(label="Node4316177488", taxon=None, edge_length=13.3651095585)
    nd_4316177872 = nd_4316177488.new_child(label="Node4316177872", taxon=None, edge_length=0.439451186875)
    nd_4316178256 = nd_4316177488.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=3.89243941597)
    nd_4316178000 = nd_4316177872.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=3.4529882291)
    nd_4316177744 = nd_4316177872.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=3.4529882291)
    nd_4316178512 = nd_4316178128.new_child(label="Node4316178512", taxon=None, edge_length=11.5242004723)
    nd_4316179280 = nd_4316178128.new_child(label="Node4316179280", taxon=None, edge_length=13.3915851761)
    nd_4316178640 = nd_4316178512.new_child(label="Node4316178640", taxon=None, edge_length=1.13984981318)
    nd_4316179024 = nd_4316178512.new_child(label="Node4316179024", taxon=None, edge_length=13.1543046788)
    nd_4316178768 = nd_4316178640.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=25.0022602166)
    nd_4316178384 = nd_4316178640.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=25.0022602166)
    nd_4316179152 = nd_4316179024.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=12.987805351)
    nd_4316178896 = nd_4316179024.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=12.987805351)
    nd_4316179536 = nd_4316179280.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=24.274725326)
    nd_4316179408 = nd_4316179280.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=24.274725326)
    nd_4316179920 = nd_4316179792.new_child(label="Node4316179920", taxon=None, edge_length=7.8146690322)
    nd_4316180432 = nd_4316179792.new_child(label="Node4316180432", taxon=None, edge_length=5.10842077756)
    nd_4316180048 = nd_4316179920.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=30.4665100305)
    nd_4316179664 = nd_4316179920.new_child(label="Node4316179664", taxon=None, edge_length=11.0043198537)
    nd_4316180304 = nd_4316179664.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=19.4621901768)
    nd_4316180176 = nd_4316179664.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=19.4621901768)
    nd_4316201104 = nd_4316180432.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=33.1727582851)
    nd_4316158224 = nd_4316180432.new_child(label="Node4316158224", taxon=None, edge_length=4.7141022378)
    nd_4316201296 = nd_4316158224.new_child(label="Node4316201296", taxon=None, edge_length=19.4308450954)
    nd_4316201744 = nd_4316158224.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=28.4586560473)
    nd_4316201424 = nd_4316201296.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=9.02781095195)
    nd_4316201616 = nd_4316201296.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=9.02781095195)
    nd_4316202000 = nd_4316201552.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=30.3039115511)
    nd_4316201872 = nd_4316201552.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=30.3039115511)
    nd_4316202384 = nd_4316202128.new_child(label="Node4316202384", taxon=None, edge_length=6.14068810478)
    nd_4316202896 = nd_4316202128.new_child(label="Python regius", taxon=tax_4313760144, edge_length=43.356182188)
    nd_4316202512 = nd_4316202384.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=37.2154940833)
    nd_4316202256 = nd_4316202384.new_child(label="Node4316202256", taxon=None, edge_length=15.2807163378)
    nd_4316202768 = nd_4316202256.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=21.9347777454)
    nd_4316202640 = nd_4316202256.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=21.9347777454)
    nd_4316203280 = nd_4316203024.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=82.94710604)
    nd_4316203152 = nd_4316203024.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=82.94710604)
    tree_4316203536 = dendropy.Tree(label="Tree04", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316203536, reindex_taxa=False)
    nd_4316203792 = tree_4316203536.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=157.750076773)
    nd_4316203984 = tree_4316203536.seed_node.new_child(label="Node4316203984", taxon=None, edge_length=44.9789688242)
    nd_4316204112 = nd_4316203984.new_child(label="Node4316204112", taxon=None, edge_length=19.5811677101)
    nd_4316252880 = nd_4316203984.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=112.771107949)
    nd_4316204240 = nd_4316204112.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=93.189940239)
    nd_4316203920 = nd_4316204112.new_child(label="Node4316203920", taxon=None, edge_length=27.291533515)
    nd_4316204496 = nd_4316203920.new_child(label="Node4316204496", taxon=None, edge_length=10.3875398007)
    nd_4316251856 = nd_4316203920.new_child(label="Node4316251856", taxon=None, edge_length=19.8895314814)
    nd_4316204624 = nd_4316204496.new_child(label="Node4316204624", taxon=None, edge_length=9.74044751021)
    nd_4316251600 = nd_4316204496.new_child(label="Node4316251600", taxon=None, edge_length=22.1202819118)
    nd_4316204752 = nd_4316204624.new_child(label="Node4316204752", taxon=None, edge_length=2.88538115771)
    nd_4316227472 = nd_4316204624.new_child(label="Node4316227472", taxon=None, edge_length=2.72676824396)
    nd_4316204880 = nd_4316204752.new_child(label="Node4316204880", taxon=None, edge_length=1.34047834167)
    nd_4316227216 = nd_4316204752.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=42.8850382554)
    nd_4316205008 = nd_4316204880.new_child(label="Node4316205008", taxon=None, edge_length=2.31767871982)
    nd_4316226768 = nd_4316204880.new_child(label="Node4316226768", taxon=None, edge_length=18.2547475502)
    nd_4316225680 = nd_4316205008.new_child(label="Node4316225680", taxon=None, edge_length=6.3930928479)
    nd_4316226448 = nd_4316205008.new_child(label="Node4316226448", taxon=None, edge_length=23.2404397828)
    nd_4316225808 = nd_4316225680.new_child(label="Node4316225808", taxon=None, edge_length=24.6792964706)
    nd_4316226192 = nd_4316225680.new_child(label="Node4316226192", taxon=None, edge_length=8.52936801714)
    nd_4316225936 = nd_4316225808.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=8.15449187544)
    nd_4316204368 = nd_4316225808.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=8.15449187544)
    nd_4316226320 = nd_4316226192.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=24.3044203289)
    nd_4316226128 = nd_4316226192.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=24.3044203289)
    nd_4316226640 = nd_4316226448.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=15.9864414111)
    nd_4316226832 = nd_4316226448.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=15.9864414111)
    nd_4316227088 = nd_4316226768.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=23.2898123636)
    nd_4316226960 = nd_4316226768.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=23.2898123636)
    nd_4316227600 = nd_4316227472.new_child(label="Node4316227600", taxon=None, edge_length=14.2175774566)
    nd_4316229392 = nd_4316227472.new_child(label="Node4316229392", taxon=None, edge_length=3.9347409374)
    nd_4316227728 = nd_4316227600.new_child(label="Node4316227728", taxon=None, edge_length=12.5474231006)
    nd_4316228112 = nd_4316227600.new_child(label="Node4316228112", taxon=None, edge_length=1.26678175478)
    nd_4316227856 = nd_4316227728.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=16.278650612)
    nd_4316227344 = nd_4316227728.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=16.278650612)
    nd_4316228240 = nd_4316228112.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=27.5592919578)
    nd_4316227984 = nd_4316228112.new_child(label="Node4316227984", taxon=None, edge_length=13.3039152583)
    nd_4316228496 = nd_4316227984.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=14.2553766995)
    nd_4316228368 = nd_4316227984.new_child(label="Node4316228368", taxon=None, edge_length=1.66170791236)
    nd_4316228752 = nd_4316228368.new_child(label="Node4316228752", taxon=None, edge_length=8.89489387836)
    nd_4316229264 = nd_4316228368.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=12.5936687872)
    nd_4316228880 = nd_4316228752.new_child(label="Node4316228880", taxon=None, edge_length=0.230110019205)
    nd_4316229136 = nd_4316228752.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=3.69877490882)
    nd_4316229008 = nd_4316228880.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=3.46866488962)
    nd_4316203408 = nd_4316228880.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=3.46866488962)
    nd_4316229520 = nd_4316229392.new_child(label="Node4316229520", taxon=None, edge_length=5.88218975316)
    nd_4316250576 = nd_4316229392.new_child(label="Node4316250576", taxon=None, edge_length=7.11128547149)
    nd_4316250192 = nd_4316229520.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=33.2267204786)
    nd_4316203600 = nd_4316229520.new_child(label="Node4316203600", taxon=None, edge_length=13.4306199458)
    nd_4316250448 = nd_4316203600.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=19.7961005329)
    nd_4316250320 = nd_4316203600.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=19.7961005329)
    nd_4316250832 = nd_4316250576.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=31.9976247603)
    nd_4316250704 = nd_4316250576.new_child(label="Node4316250704", taxon=None, edge_length=4.71436528425)
    nd_4316251088 = nd_4316250704.new_child(label="Node4316251088", taxon=None, edge_length=15.9285543528)
    nd_4316251472 = nd_4316250704.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=27.2832594761)
    nd_4316251216 = nd_4316251088.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=11.3547051232)
    nd_4316250960 = nd_4316251088.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=11.3547051232)
    nd_4316251728 = nd_4316251600.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=33.3905850116)
    nd_4316251344 = nd_4316251600.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=33.3905850116)
    nd_4316252112 = nd_4316251856.new_child(label="Node4316252112", taxon=None, edge_length=11.1907198563)
    nd_4316252624 = nd_4316251856.new_child(label="Python regius", taxon=tax_4313760144, edge_length=46.0088752426)
    nd_4316252240 = nd_4316252112.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=34.8181553863)
    nd_4316251984 = nd_4316252112.new_child(label="Node4316251984", taxon=None, edge_length=7.89583224277)
    nd_4316252496 = nd_4316251984.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=26.9223231435)
    nd_4316252368 = nd_4316251984.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=26.9223231435)
    tree_4316252752 = dendropy.Tree(label="Tree05", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316252752, reindex_taxa=False)
    nd_4316253264 = tree_4316252752.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=117.405482731)
    nd_4316253456 = tree_4316252752.seed_node.new_child(label="Node4316253456", taxon=None, edge_length=26.8517813278)
    nd_4316253584 = nd_4316253456.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=90.5537014032)
    nd_4316253392 = nd_4316253456.new_child(label="Node4316253392", taxon=None, edge_length=16.5508683851)
    nd_4316253776 = nd_4316253392.new_child(label="Node4316253776", taxon=None, edge_length=16.6353899137)
    nd_4316302224 = nd_4316253392.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=74.0028330181)
    nd_4316253904 = nd_4316253776.new_child(label="Node4316253904", taxon=None, edge_length=6.56321044476)
    nd_4316301328 = nd_4316253776.new_child(label="Node4316301328", taxon=None, edge_length=18.1852453647)
    nd_4316254032 = nd_4316253904.new_child(label="Node4316254032", taxon=None, edge_length=12.3993667148)
    nd_4316301072 = nd_4316253904.new_child(label="Node4316301072", taxon=None, edge_length=23.3069353709)
    nd_4316254160 = nd_4316254032.new_child(label="Node4316254160", taxon=None, edge_length=1.85501747057)
    nd_4316276944 = nd_4316254032.new_child(label="Node4316276944", taxon=None, edge_length=2.5597754437)
    nd_4316274832 = nd_4316254160.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=36.5498484742)
    nd_4316275024 = nd_4316254160.new_child(label="Node4316275024", taxon=None, edge_length=1.04661210215)
    nd_4316275152 = nd_4316275024.new_child(label="Node4316275152", taxon=None, edge_length=6.27408700912)
    nd_4316276432 = nd_4316275024.new_child(label="Node4316276432", taxon=None, edge_length=15.4075337774)
    nd_4316275280 = nd_4316275152.new_child(label="Node4316275280", taxon=None, edge_length=4.14502379891)
    nd_4316276048 = nd_4316275152.new_child(label="Node4316276048", taxon=None, edge_length=17.4625828861)
    nd_4316275408 = nd_4316275280.new_child(label="Node4316275408", taxon=None, edge_length=18.7911640437)
    nd_4316275792 = nd_4316275280.new_child(label="Node4316275792", taxon=None, edge_length=4.65165886356)
    nd_4316275536 = nd_4316275408.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=6.29296152027)
    nd_4316274960 = nd_4316275408.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=6.29296152027)
    nd_4316275920 = nd_4316275792.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=20.4324667004)
    nd_4316275664 = nd_4316275792.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=20.4324667004)
    nd_4316276304 = nd_4316276048.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=11.7665664768)
    nd_4316276176 = nd_4316276048.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=11.7665664768)
    nd_4316276688 = nd_4316276432.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=20.0957025946)
    nd_4316276560 = nd_4316276432.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=20.0957025946)
    nd_4316277072 = nd_4316276944.new_child(label="Node4316277072", taxon=None, edge_length=10.9999278199)
    nd_4316299408 = nd_4316276944.new_child(label="Node4316299408", taxon=None, edge_length=4.77611939702)
    nd_4316277200 = nd_4316277072.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=24.8451626811)
    nd_4316276816 = nd_4316277072.new_child(label="Node4316276816", taxon=None, edge_length=2.78884378851)
    nd_4316277392 = nd_4316276816.new_child(label="Node4316277392", taxon=None, edge_length=7.72979171285)
    nd_4316278352 = nd_4316276816.new_child(label="Node4316278352", taxon=None, edge_length=12.4767981898)
    nd_4316277520 = nd_4316277392.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=14.3265271797)
    nd_4316277584 = nd_4316277392.new_child(label="Node4316277584", taxon=None, edge_length=3.38770784397)
    nd_4316277712 = nd_4316277584.new_child(label="Node4316277712", taxon=None, edge_length=5.65307392619)
    nd_4316278224 = nd_4316277584.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=10.9388193358)
    nd_4316277840 = nd_4316277712.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=5.28574540958)
    nd_4316253008 = nd_4316277712.new_child(label="Node4316253008", taxon=None, edge_length=2.55413552294)
    nd_4316278096 = nd_4316253008.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.73160988664)
    nd_4316277968 = nd_4316253008.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=2.73160988664)
    nd_4316278608 = nd_4316278352.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=9.57952070277)
    nd_4316278480 = nd_4316278352.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=9.57952070277)
    nd_4316299536 = nd_4316299408.new_child(label="Node4316299536", taxon=None, edge_length=6.62946104565)
    nd_4316300048 = nd_4316299408.new_child(label="Node4316300048", taxon=None, edge_length=3.07031045323)
    nd_4316299664 = nd_4316299536.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=24.4395100584)
    nd_4316299344 = nd_4316299536.new_child(label="Node4316299344", taxon=None, edge_length=7.74158436759)
    nd_4316299920 = nd_4316299344.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=16.6979256908)
    nd_4316299792 = nd_4316299344.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=16.6979256908)
    nd_4316300304 = nd_4316300048.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=27.9986606508)
    nd_4316300176 = nd_4316300048.new_child(label="Node4316300176", taxon=None, edge_length=4.88214307652)
    nd_4316300560 = nd_4316300176.new_child(label="Node4316300560", taxon=None, edge_length=16.2064963991)
    nd_4316300944 = nd_4316300176.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=23.1165175743)
    nd_4316300688 = nd_4316300560.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=6.91002117518)
    nd_4316300432 = nd_4316300560.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=6.91002117518)
    nd_4316301200 = nd_4316301072.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=27.4972972887)
    nd_4316300816 = nd_4316301072.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=27.4972972887)
    nd_4316301584 = nd_4316301328.new_child(label="Node4316301584", taxon=None, edge_length=6.08258726043)
    nd_4316302096 = nd_4316301328.new_child(label="Python regius", taxon=tax_4313760144, edge_length=39.1821977396)
    nd_4316301712 = nd_4316301584.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=33.0996104792)
    nd_4316301456 = nd_4316301584.new_child(label="Node4316301456", taxon=None, edge_length=14.293345062)
    nd_4316301968 = nd_4316301456.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=18.8062654172)
    nd_4316301840 = nd_4316301456.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=18.8062654172)
    tree_4316302352 = dendropy.Tree(label="Tree06", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316302352, reindex_taxa=False)
    nd_4316302736 = tree_4316302352.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=126.159749307)
    nd_4316302928 = tree_4316302352.seed_node.new_child(label="Node4316302928", taxon=None, edge_length=29.5182899091)
    nd_4316303056 = nd_4316302928.new_child(label="Node4316303056", taxon=None, edge_length=20.4272017262)
    nd_4316347728 = nd_4316302928.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=96.6414593981)
    nd_4316303184 = nd_4316303056.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=76.2142576718)
    nd_4316302864 = nd_4316303056.new_child(label="Node4316302864", taxon=None, edge_length=22.3340694246)
    nd_4316323920 = nd_4316302864.new_child(label="Node4316323920", taxon=None, edge_length=7.33856679933)
    nd_4316346704 = nd_4316302864.new_child(label="Node4316346704", taxon=None, edge_length=17.9081848887)
    nd_4316324048 = nd_4316323920.new_child(label="Node4316324048", taxon=None, edge_length=8.24675869634)
    nd_4316346448 = nd_4316323920.new_child(label="Node4316346448", taxon=None, edge_length=23.556392099)
    nd_4316324176 = nd_4316324048.new_child(label="Node4316324176", taxon=None, edge_length=1.78294733026)
    nd_4316326288 = nd_4316324048.new_child(label="Node4316326288", taxon=None, edge_length=3.26747919544)
    nd_4316324304 = nd_4316324176.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=36.5119154213)
    nd_4316324496 = nd_4316324176.new_child(label="Node4316324496", taxon=None, edge_length=1.38475127065)
    nd_4316324624 = nd_4316324496.new_child(label="Node4316324624", taxon=None, edge_length=2.24740276648)
    nd_4316326032 = nd_4316324496.new_child(label="Node4316326032", taxon=None, edge_length=12.9799500845)
    nd_4316324752 = nd_4316324624.new_child(label="Node4316324752", taxon=None, edge_length=17.772328432)
    nd_4316325136 = nd_4316324624.new_child(label="Node4316325136", taxon=None, edge_length=4.4885389798)
    nd_4316324880 = nd_4316324752.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=15.1074329521)
    nd_4316324432 = nd_4316324752.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=15.1074329521)
    nd_4316325264 = nd_4316325136.new_child(label="Node4316325264", taxon=None, edge_length=20.8200876951)
    nd_4316325648 = nd_4316325136.new_child(label="Node4316325648", taxon=None, edge_length=4.37289319177)
    nd_4316325392 = nd_4316325264.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=7.5711347092)
    nd_4316325008 = nd_4316325264.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=7.5711347092)
    nd_4316325776 = nd_4316325648.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=24.0183292126)
    nd_4316325520 = nd_4316325648.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=24.0183292126)
    nd_4316326160 = nd_4316326032.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=22.1472140662)
    nd_4316325904 = nd_4316326032.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=22.1472140662)
    nd_4316326416 = nd_4316326288.new_child(label="Node4316326416", taxon=None, edge_length=12.2960670728)
    nd_4316344784 = nd_4316326288.new_child(label="Node4316344784", taxon=None, edge_length=4.17834235089)
    nd_4316326544 = nd_4316326416.new_child(label="Node4316326544", taxon=None, edge_length=1.62908011033)
    nd_4316327184 = nd_4316326416.new_child(label="Node4316327184", taxon=None, edge_length=6.20414166284)
    nd_4316326672 = nd_4316326544.new_child(label="Node4316326672", taxon=None, edge_length=11.7610480778)
    nd_4316327056 = nd_4316326544.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=21.1022363729)
    nd_4316326800 = nd_4316326672.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=9.34118829512)
    nd_4316302480 = nd_4316326672.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=9.34118829512)
    nd_4316327312 = nd_4316327184.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=16.5271748204)
    nd_4316326928 = nd_4316327184.new_child(label="Node4316326928", taxon=None, edge_length=4.27339647744)
    nd_4316327568 = nd_4316326928.new_child(label="Node4316327568", taxon=None, edge_length=7.91289243511)
    nd_4316344656 = nd_4316326928.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=12.253778343)
    nd_4316327696 = nd_4316327568.new_child(label="Node4316327696", taxon=None, edge_length=1.89524872622)
    nd_4316344528 = nd_4316327568.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=4.34088590785)
    nd_4316327824 = nd_4316327696.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=2.44563718163)
    nd_4316327440 = nd_4316327696.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.44563718163)
    nd_4316344912 = nd_4316344784.new_child(label="Node4316344912", taxon=None, edge_length=5.03710806168)
    nd_4316345424 = nd_4316344784.new_child(label="Node4316345424", taxon=None, edge_length=4.09757269601)
    nd_4316345040 = nd_4316344912.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=25.8119331435)
    nd_4316344400 = nd_4316344912.new_child(label="Node4316344400", taxon=None, edge_length=11.4337686931)
    nd_4316345296 = nd_4316344400.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=14.3781644505)
    nd_4316345168 = nd_4316344400.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=14.3781644505)
    nd_4316345680 = nd_4316345424.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=26.7514685092)
    nd_4316345552 = nd_4316345424.new_child(label="Node4316345552", taxon=None, edge_length=2.33638652753)
    nd_4316345936 = nd_4316345552.new_child(label="Node4316345936", taxon=None, edge_length=17.4436719435)
    nd_4316346320 = nd_4316345552.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=24.4150819817)
    nd_4316346064 = nd_4316345936.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=6.97141003814)
    nd_4316345808 = nd_4316345936.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=6.97141003814)
    nd_4316346576 = nd_4316346448.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=22.9852293488)
    nd_4316346192 = nd_4316346448.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=22.9852293488)
    nd_4316346960 = nd_4316346704.new_child(label="Node4316346960", taxon=None, edge_length=6.75664224871)
    nd_4316347472 = nd_4316346704.new_child(label="Python regius", taxon=tax_4313760144, edge_length=35.9720033585)
    nd_4316347088 = nd_4316346960.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=29.2153611098)
    nd_4316346832 = nd_4316346960.new_child(label="Node4316346832", taxon=None, edge_length=8.15978945225)
    nd_4316347344 = nd_4316346832.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=21.0555716576)
    nd_4316347216 = nd_4316346832.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=21.0555716576)
    tree_4316347600 = dendropy.Tree(label="Tree07", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316347600, reindex_taxa=False)
    nd_4316348112 = tree_4316347600.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=124.564186516)
    nd_4316348304 = tree_4316347600.seed_node.new_child(label="Node4316348304", taxon=None, edge_length=36.3676780441)
    nd_4316368976 = nd_4316348304.new_child(label="Node4316368976", taxon=None, edge_length=11.1789504571)
    nd_4316397200 = nd_4316348304.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=88.196508472)
    nd_4316369104 = nd_4316368976.new_child(label="Node4316369104", taxon=None, edge_length=20.3346663059)
    nd_4316396944 = nd_4316368976.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=77.0175580149)
    nd_4316369232 = nd_4316369104.new_child(label="Node4316369232", taxon=None, edge_length=10.7040090023)
    nd_4316396048 = nd_4316369104.new_child(label="Node4316396048", taxon=None, edge_length=19.1477140595)
    nd_4316369360 = nd_4316369232.new_child(label="Node4316369360", taxon=None, edge_length=7.28944403695)
    nd_4316395664 = nd_4316369232.new_child(label="Node4316395664", taxon=None, edge_length=20.035415736)
    nd_4316369488 = nd_4316369360.new_child(label="Node4316369488", taxon=None, edge_length=2.76745367097)
    nd_4316372752 = nd_4316369360.new_child(label="Node4316372752", taxon=None, edge_length=1.09335779327)
    nd_4316369616 = nd_4316369488.new_child(label="Node4316369616", taxon=None, edge_length=12.5842551962)
    nd_4316371280 = nd_4316369488.new_child(label="Node4316371280", taxon=None, edge_length=7.66239308607)
    nd_4316369744 = nd_4316369616.new_child(label="Node4316369744", taxon=None, edge_length=12.0191169658)
    nd_4316370128 = nd_4316369616.new_child(label="Node4316370128", taxon=None, edge_length=2.25737379321)
    nd_4316369872 = nd_4316369744.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=11.3186128368)
    nd_4316348240 = nd_4316369744.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=11.3186128368)
    nd_4316370256 = nd_4316370128.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=21.0803560094)
    nd_4316370000 = nd_4316370128.new_child(label="Node4316370000", taxon=None, edge_length=6.93019779087)
    nd_4316370512 = nd_4316370000.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=14.1501582185)
    nd_4316370384 = nd_4316370000.new_child(label="Node4316370384", taxon=None, edge_length=3.68709033509)
    nd_4316347856 = nd_4316370384.new_child(label="Node4316347856", taxon=None, edge_length=4.34576715349)
    nd_4316371152 = nd_4316370384.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=10.4630678834)
    nd_4316370768 = nd_4316347856.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=6.11730072991)
    nd_4316370640 = nd_4316347856.new_child(label="Node4316370640", taxon=None, edge_length=2.07374196069)
    nd_4316371024 = nd_4316370640.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=4.04355876922)
    nd_4316370896 = nd_4316370640.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=4.04355876922)
    nd_4316371536 = nd_4316371280.new_child(label="Node4316371536", taxon=None, edge_length=5.79725116231)
    nd_4316372048 = nd_4316371280.new_child(label="Node4316372048", taxon=None, edge_length=3.71811515781)
    nd_4316371664 = nd_4316371536.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=22.4623407504)
    nd_4316371408 = nd_4316371536.new_child(label="Node4316371408", taxon=None, edge_length=5.91246822929)
    nd_4316371920 = nd_4316371408.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=16.5498725211)
    nd_4316371792 = nd_4316371408.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=16.5498725211)
    nd_4316372304 = nd_4316372048.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=24.5414767549)
    nd_4316372176 = nd_4316372048.new_child(label="Node4316372176", taxon=None, edge_length=5.68085904285)
    nd_4316372496 = nd_4316372176.new_child(label="Node4316372496", taxon=None, edge_length=11.551508813)
    nd_4316372944 = nd_4316372176.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=18.8606177121)
    nd_4316372624 = nd_4316372496.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=7.30910889905)
    nd_4316372816 = nd_4316372496.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=7.30910889905)
    nd_4316393744 = nd_4316372752.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=37.5960808765)
    nd_4316393616 = nd_4316372752.new_child(label="Node4316393616", taxon=None, edge_length=2.73294846613)
    nd_4316393936 = nd_4316393616.new_child(label="Node4316393936", taxon=None, edge_length=2.22916081797)
    nd_4316395408 = nd_4316393616.new_child(label="Node4316395408", taxon=None, edge_length=13.2144705479)
    nd_4316394064 = nd_4316393936.new_child(label="Node4316394064", taxon=None, edge_length=19.1660901297)
    nd_4316394512 = nd_4316393936.new_child(label="Node4316394512", taxon=None, edge_length=7.55440984409)
    nd_4316394192 = nd_4316394064.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=13.4678814627)
    nd_4316394384 = nd_4316394064.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=13.4678814627)
    nd_4316394640 = nd_4316394512.new_child(label="Node4316394640", taxon=None, edge_length=18.9933631628)
    nd_4316395024 = nd_4316394512.new_child(label="Node4316395024", taxon=None, edge_length=5.77339164759)
    nd_4316394768 = nd_4316394640.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=6.08619858554)
    nd_4316394320 = nd_4316394640.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=6.08619858554)
    nd_4316395152 = nd_4316395024.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=19.3061701007)
    nd_4316394896 = nd_4316395024.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=19.3061701007)
    nd_4316395536 = nd_4316395408.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=21.6486618625)
    nd_4316395280 = nd_4316395408.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=21.6486618625)
    nd_4316395920 = nd_4316395664.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=25.9434669707)
    nd_4316395792 = nd_4316395664.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=25.9434669707)
    nd_4316396304 = nd_4316396048.new_child(label="Node4316396304", taxon=None, edge_length=7.07244422644)
    nd_4316396816 = nd_4316396048.new_child(label="Python regius", taxon=tax_4313760144, edge_length=37.5351776494)
    nd_4316396432 = nd_4316396304.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=30.462733423)
    nd_4316396176 = nd_4316396304.new_child(label="Node4316396176", taxon=None, edge_length=14.4432592042)
    nd_4316396688 = nd_4316396176.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=16.0194742188)
    nd_4316396560 = nd_4316396176.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=16.0194742188)
    tree_4316397072 = dendropy.Tree(label="Tree08", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316397072, reindex_taxa=False)
    nd_4316418128 = tree_4316397072.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=95.8502441646)
    nd_4316418320 = tree_4316397072.seed_node.new_child(label="Node4316418320", taxon=None, edge_length=21.8741644934)
    nd_4316418448 = nd_4316418320.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=73.9760796713)
    nd_4316418256 = nd_4316418320.new_child(label="Node4316418256", taxon=None, edge_length=9.52951598189)
    nd_4316418704 = nd_4316418256.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=64.4465636894)
    nd_4316418576 = nd_4316418256.new_child(label="Node4316418576", taxon=None, edge_length=22.9882151659)
    nd_4316418960 = nd_4316418576.new_child(label="Node4316418960", taxon=None, edge_length=6.43697452247)
    nd_4316445776 = nd_4316418576.new_child(label="Node4316445776", taxon=None, edge_length=13.8002436509)
    nd_4316419088 = nd_4316418960.new_child(label="Node4316419088", taxon=None, edge_length=6.7231540226)
    nd_4316445392 = nd_4316418960.new_child(label="Node4316445392", taxon=None, edge_length=15.1871703268)
    nd_4316397392 = nd_4316419088.new_child(label="Node4316397392", taxon=None, edge_length=1.03257897787)
    nd_4316443088 = nd_4316419088.new_child(label="Node4316443088", taxon=None, edge_length=1.21568464049)
    nd_4316419216 = nd_4316397392.new_child(label="Node4316419216", taxon=None, edge_length=7.52518318022)
    nd_4316421136 = nd_4316397392.new_child(label="Node4316421136", taxon=None, edge_length=2.39608223465)
    nd_4316419344 = nd_4316419216.new_child(label="Node4316419344", taxon=None, edge_length=2.42060556338)
    nd_4316421008 = nd_4316419216.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=19.7404578203)
    nd_4316419472 = nd_4316419344.new_child(label="Node4316419472", taxon=None, edge_length=8.29264517113)
    nd_4316419856 = nd_4316419344.new_child(label="Node4316419856", taxon=None, edge_length=6.70163613113)
    nd_4316419600 = nd_4316419472.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=9.02720708579)
    nd_4316418832 = nd_4316419472.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=9.02720708579)
    nd_4316419984 = nd_4316419856.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=10.6182161258)
    nd_4316419728 = nd_4316419856.new_child(label="Node4316419728", taxon=None, edge_length=3.47880840545)
    nd_4316420240 = nd_4316419728.new_child(label="Node4316420240", taxon=None, edge_length=4.27223311967)
    nd_4316420752 = nd_4316419728.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=7.13940772034)
    nd_4316420368 = nd_4316420240.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=2.86717460067)
    nd_4316420112 = nd_4316420240.new_child(label="Node4316420112", taxon=None, edge_length=0.371215520464)
    nd_4316420624 = nd_4316420112.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=2.49595908021)
    nd_4316420496 = nd_4316420112.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.49595908021)
    nd_4316421264 = nd_4316421136.new_child(label="Node4316421264", taxon=None, edge_length=3.24355281662)
    nd_4316421776 = nd_4316421136.new_child(label="Node4316421776", taxon=None, edge_length=3.2027013644)
    nd_4316421392 = nd_4316421264.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=21.6260059492)
    nd_4316420880 = nd_4316421264.new_child(label="Node4316420880", taxon=None, edge_length=9.1533399099)
    nd_4316421648 = nd_4316420880.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=12.4726660393)
    nd_4316421520 = nd_4316420880.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=12.4726660393)
    nd_4316422032 = nd_4316421776.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=21.6668574015)
    nd_4316421904 = nd_4316421776.new_child(label="Node4316421904", taxon=None, edge_length=4.33602073619)
    nd_4316442832 = nd_4316421904.new_child(label="Node4316442832", taxon=None, edge_length=11.5233229214)
    nd_4316443216 = nd_4316421904.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=17.3308366653)
    nd_4316442960 = nd_4316442832.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=5.8075137439)
    nd_4316442704 = nd_4316442832.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=5.8075137439)
    nd_4316443472 = nd_4316443088.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=27.0825353379)
    nd_4316443344 = nd_4316443088.new_child(label="Node4316443344", taxon=None, edge_length=1.37229916019)
    nd_4316443728 = nd_4316443344.new_child(label="Node4316443728", taxon=None, edge_length=2.64946637554)
    nd_4316445136 = nd_4316443344.new_child(label="Node4316445136", taxon=None, edge_length=11.4050202795)
    nd_4316443856 = nd_4316443728.new_child(label="Node4316443856", taxon=None, edge_length=13.5545767859)
    nd_4316444240 = nd_4316443728.new_child(label="Node4316444240", taxon=None, edge_length=4.67390676307)
    nd_4316443984 = nd_4316443856.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=9.50619301624)
    nd_4316443600 = nd_4316443856.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=9.50619301624)
    nd_4316444368 = nd_4316444240.new_child(label="Node4316444368", taxon=None, edge_length=12.8995814401)
    nd_4316444752 = nd_4316444240.new_child(label="Node4316444752", taxon=None, edge_length=1.38849394051)
    nd_4316444496 = nd_4316444368.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=5.48728159904)
    nd_4316444112 = nd_4316444368.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=5.48728159904)
    nd_4316444880 = nd_4316444752.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=16.9983690986)
    nd_4316444624 = nd_4316444752.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=16.9983690986)
    nd_4316445264 = nd_4316445136.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=14.3052158982)
    nd_4316445008 = nd_4316445136.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=14.3052158982)
    nd_4316445648 = nd_4316445392.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=19.8342036742)
    nd_4316445520 = nd_4316445392.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=19.8342036742)
    nd_4316446032 = nd_4316445776.new_child(label="Node4316446032", taxon=None, edge_length=3.93234771773)
    nd_4316446544 = nd_4316445776.new_child(label="Python regius", taxon=tax_4313760144, edge_length=27.6581048725)
    nd_4316446160 = nd_4316446032.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=23.7257571548)
    nd_4316445904 = nd_4316446032.new_child(label="Node4316445904", taxon=None, edge_length=7.63118798191)
    nd_4316446416 = nd_4316445904.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=16.0945691729)
    nd_4316446288 = nd_4316445904.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=16.0945691729)
    tree_4316446672 = dendropy.Tree(label="Tree09", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316446672, reindex_taxa=False)
    nd_4316467600 = tree_4316446672.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=126.147943419)
    nd_4316467344 = tree_4316446672.seed_node.new_child(label="Node4316467344", taxon=None, edge_length=40.5157695372)
    nd_4316467792 = nd_4316467344.new_child(label="Node4316467792", taxon=None, edge_length=13.6608326978)
    nd_4316512592 = nd_4316467344.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=85.6321738818)
    nd_4316467920 = nd_4316467792.new_child(label="Node4316467920", taxon=None, edge_length=17.1918011574)
    nd_4316512336 = nd_4316467792.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=71.971341184)
    nd_4316468048 = nd_4316467920.new_child(label="Node4316468048", taxon=None, edge_length=5.77299961102)
    nd_4316490896 = nd_4316467920.new_child(label="Node4316490896", taxon=None, edge_length=21.166651853)
    nd_4316468176 = nd_4316468048.new_child(label="Node4316468176", taxon=None, edge_length=13.7010200713)
    nd_4316490576 = nd_4316468048.new_child(label="Node4316490576", taxon=None, edge_length=24.6978541753)
    nd_4316468304 = nd_4316468176.new_child(label="Node4316468304", taxon=None, edge_length=0.617658990864)
    nd_4316470480 = nd_4316468176.new_child(label="Node4316470480", taxon=None, edge_length=1.35670146604)
    nd_4316468432 = nd_4316468304.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=34.6878613533)
    nd_4316467408 = nd_4316468304.new_child(label="Node4316467408", taxon=None, edge_length=2.23641376685)
    nd_4316468688 = nd_4316467408.new_child(label="Node4316468688", taxon=None, edge_length=10.5685372807)
    nd_4316469072 = nd_4316467408.new_child(label="Node4316469072", taxon=None, edge_length=2.11685229318)
    nd_4316468816 = nd_4316468688.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=21.8829103058)
    nd_4316468560 = nd_4316468688.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=21.8829103058)
    nd_4316469200 = nd_4316469072.new_child(label="Node4316469200", taxon=None, edge_length=17.1757724208)
    nd_4316469584 = nd_4316469072.new_child(label="Node4316469584", taxon=None, edge_length=3.96372676423)
    nd_4316469328 = nd_4316469200.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=13.1588228725)
    nd_4316468944 = nd_4316469200.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=13.1588228725)
    nd_4316469712 = nd_4316469584.new_child(label="Node4316469712", taxon=None, edge_length=19.5683902852)
    nd_4316470096 = nd_4316469584.new_child(label="Node4316470096", taxon=None, edge_length=4.82785669688)
    nd_4316469840 = nd_4316469712.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=6.80247824391)
    nd_4316469456 = nd_4316469712.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=6.80247824391)
    nd_4316470224 = nd_4316470096.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=21.5430118322)
    nd_4316469968 = nd_4316470096.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=21.5430118322)
    nd_4316470608 = nd_4316470480.new_child(label="Node4316470608", taxon=None, edge_length=10.3401216686)
    nd_4316488848 = nd_4316470480.new_child(label="Node4316488848", taxon=None, edge_length=7.03488295134)
    nd_4316470736 = nd_4316470608.new_child(label="Node4316470736", taxon=None, edge_length=2.31259127041)
    nd_4316488592 = nd_4316470608.new_child(label="Node4316488592", taxon=None, edge_length=13.9826784484)
    nd_4316470864 = nd_4316470736.new_child(label="Node4316470864", taxon=None, edge_length=7.79402846351)
    nd_4316488336 = nd_4316470736.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=21.2961059391)
    nd_4316470992 = nd_4316470864.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=13.5020774756)
    nd_4316470352 = nd_4316470864.new_child(label="Node4316470352", taxon=None, edge_length=3.5072599877)
    nd_4316471248 = nd_4316470352.new_child(label="Node4316471248", taxon=None, edge_length=6.20487512451)
    nd_4316488208 = nd_4316470352.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=9.99481748791)
    nd_4316487824 = nd_4316471248.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=3.78994236341)
    nd_4316471120 = nd_4316471248.new_child(label="Node4316471120", taxon=None, edge_length=1.25204312348)
    nd_4316488080 = nd_4316471120.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=2.53789923993)
    nd_4316487952 = nd_4316471120.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.53789923993)
    nd_4316488720 = nd_4316488592.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=9.62601876119)
    nd_4316488464 = nd_4316488592.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=9.62601876119)
    nd_4316489104 = nd_4316488848.new_child(label="Node4316489104", taxon=None, edge_length=4.39672345568)
    nd_4316489616 = nd_4316488848.new_child(label="Node4316489616", taxon=None, edge_length=3.94778865925)
    nd_4316489232 = nd_4316489104.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=22.5172124711)
    nd_4316488976 = nd_4316489104.new_child(label="Node4316488976", taxon=None, edge_length=11.9061835258)
    nd_4316489488 = nd_4316488976.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=10.6110289453)
    nd_4316489360 = nd_4316488976.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=10.6110289453)
    nd_4316489872 = nd_4316489616.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=22.9661472676)
    nd_4316489744 = nd_4316489616.new_child(label="Node4316489744", taxon=None, edge_length=3.7713704432)
    nd_4316490128 = nd_4316489744.new_child(label="Node4316490128", taxon=None, edge_length=11.9409692012)
    nd_4316490512 = nd_4316489744.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=19.1947768244)
    nd_4316490256 = nd_4316490128.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=7.25380762314)
    nd_4316490000 = nd_4316490128.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=7.25380762314)
    nd_4316490768 = nd_4316490576.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=24.3086862402)
    nd_4316490640 = nd_4316490576.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=24.3086862402)
    nd_4316491152 = nd_4316490896.new_child(label="Node4316491152", taxon=None, edge_length=5.80983351578)
    nd_4316491664 = nd_4316490896.new_child(label="Python regius", taxon=tax_4313760144, edge_length=33.6128881735)
    nd_4316491280 = nd_4316491152.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=27.8030546577)
    nd_4316491024 = nd_4316491152.new_child(label="Node4316491024", taxon=None, edge_length=10.4395768579)
    nd_4316491536 = nd_4316491024.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=17.3634777999)
    nd_4316491408 = nd_4316491024.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=17.3634777999)
    tree_4316512464 = dendropy.Tree(label="Tree10", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316512464, reindex_taxa=False)
    nd_4316512976 = tree_4316512464.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=146.770054852)
    nd_4316513168 = tree_4316512464.seed_node.new_child(label="Node4316513168", taxon=None, edge_length=49.9930471528)
    nd_4316513296 = nd_4316513168.new_child(label="Node4316513296", taxon=None, edge_length=13.4634525107)
    nd_4316558096 = nd_4316513168.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=96.7770076997)
    nd_4316513424 = nd_4316513296.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=83.313555189)
    nd_4316513104 = nd_4316513296.new_child(label="Node4316513104", taxon=None, edge_length=30.4776798161)
    nd_4316513680 = nd_4316513104.new_child(label="Node4316513680", taxon=None, edge_length=8.15039409252)
    nd_4316536464 = nd_4316513104.new_child(label="Node4316536464", taxon=None, edge_length=14.4274535052)
    nd_4316513808 = nd_4316513680.new_child(label="Node4316513808", taxon=None, edge_length=9.51023675819)
    nd_4316512784 = nd_4316513680.new_child(label="Node4316512784", taxon=None, edge_length=25.3895116055)
    nd_4316513936 = nd_4316513808.new_child(label="Node4316513936", taxon=None, edge_length=1.94845077373)
    nd_4316534096 = nd_4316513808.new_child(label="Node4316534096", taxon=None, edge_length=1.20705242999)
    nd_4316514064 = nd_4316513936.new_child(label="Node4316514064", taxon=None, edge_length=3.25905094181)
    nd_4316515728 = nd_4316513936.new_child(label="Node4316515728", taxon=None, edge_length=10.0348009179)
    nd_4316514192 = nd_4316514064.new_child(label="Node4316514192", taxon=None, edge_length=4.82232736143)
    nd_4316514704 = nd_4316514064.new_child(label="Node4316514704", taxon=None, edge_length=4.58004825968)
    nd_4316514320 = nd_4316514192.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=25.1454154452)
    nd_4316513552 = nd_4316514192.new_child(label="Node4316513552", taxon=None, edge_length=10.5788836702)
    nd_4316514576 = nd_4316513552.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=14.566531775)
    nd_4316514448 = nd_4316513552.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=14.566531775)
    nd_4316514960 = nd_4316514704.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=25.387694547)
    nd_4316514832 = nd_4316514704.new_child(label="Node4316514832", taxon=None, edge_length=5.75635440071)
    nd_4316515216 = nd_4316514832.new_child(label="Node4316515216", taxon=None, edge_length=12.7115868187)
    nd_4316515600 = nd_4316514832.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=19.6313401463)
    nd_4316515344 = nd_4316515216.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=6.91975332756)
    nd_4316515088 = nd_4316515216.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=6.91975332756)
    nd_4316515856 = nd_4316515728.new_child(label="Node4316515856", taxon=None, edge_length=14.4811675935)
    nd_4316516240 = nd_4316515728.new_child(label="Node4316516240", taxon=None, edge_length=0.795030531054)
    nd_4316515984 = nd_4316515856.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=8.71082523711)
    nd_4316515472 = nd_4316515856.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=8.71082523711)
    nd_4316532816 = nd_4316516240.new_child(label="Node4316532816", taxon=None, edge_length=7.95009719313)
    nd_4316533840 = nd_4316516240.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=22.3969622995)
    nd_4316532944 = nd_4316532816.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=14.4468651064)
    nd_4316516112 = nd_4316532816.new_child(label="Node4316516112", taxon=None, edge_length=2.27822479328)
    nd_4316533200 = nd_4316516112.new_child(label="Node4316533200", taxon=None, edge_length=6.98040063458)
    nd_4316533648 = nd_4316516112.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=12.1686403131)
    nd_4316533328 = nd_4316533200.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=5.18823967855)
    nd_4316533072 = nd_4316533200.new_child(label="Node4316533072", taxon=None, edge_length=1.71064757594)
    nd_4316533520 = nd_4316533072.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=3.47759210261)
    nd_4316533712 = nd_4316533072.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=3.47759210261)
    nd_4316534224 = nd_4316534096.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=33.9681920922)
    nd_4316533968 = nd_4316534096.new_child(label="Node4316533968", taxon=None, edge_length=1.28486469944)
    nd_4316534480 = nd_4316533968.new_child(label="Node4316534480", taxon=None, edge_length=12.4520799939)
    nd_4316534864 = nd_4316533968.new_child(label="Node4316534864", taxon=None, edge_length=1.68023264943)
    nd_4316534608 = nd_4316534480.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=20.2312473989)
    nd_4316534352 = nd_4316534480.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=20.2312473989)
    nd_4316534992 = nd_4316534864.new_child(label="Node4316534992", taxon=None, edge_length=19.0383478987)
    nd_4316535376 = nd_4316534864.new_child(label="Node4316535376", taxon=None, edge_length=4.75943584051)
    nd_4316535120 = nd_4316534992.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=11.9647468446)
    nd_4316534736 = nd_4316534992.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=11.9647468446)
    nd_4316535504 = nd_4316535376.new_child(label="Node4316535504", taxon=None, edge_length=17.1180008393)
    nd_4316535888 = nd_4316535376.new_child(label="Node4316535888", taxon=None, edge_length=2.6817531508)
    nd_4316535632 = nd_4316535504.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=9.12565806357)
    nd_4316535248 = nd_4316535504.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=9.12565806357)
    nd_4316536016 = nd_4316535888.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=23.561905752)
    nd_4316512720 = nd_4316535888.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=23.561905752)
    nd_4316536336 = nd_4316512784.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=19.2959696748)
    nd_4316536208 = nd_4316512784.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=19.2959696748)
    nd_4316536656 = nd_4316536464.new_child(label="Node4316536656", taxon=None, edge_length=6.97033769)
    nd_4316557776 = nd_4316536464.new_child(label="Python regius", taxon=tax_4313760144, edge_length=38.4084218677)
    nd_4316536784 = nd_4316536656.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=31.4380841777)
    nd_4316557520 = nd_4316536656.new_child(label="Node4316557520", taxon=None, edge_length=13.5557299793)
    nd_4316557648 = nd_4316557520.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=17.8823541984)
    nd_4316557456 = nd_4316557520.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=17.8823541984)
    tree_4316557904 = dendropy.Tree(label="Tree11", taxon_namespace=tree_list.taxon_namespace)
    tree_list.append(tree_4316557904, reindex_taxa=False)
    nd_4316558480 = tree_4316557904.seed_node.new_child(label="Node4316558480", taxon=None, edge_length=71.3865451194)
    nd_4316607440 = tree_4316557904.seed_node.new_child(label="Candoia aspera", taxon=tax_4313742032, edge_length=211.46817309)
    nd_4316558608 = nd_4316558480.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=140.081627971)
    nd_4316558800 = nd_4316558480.new_child(label="Node4316558800", taxon=None, edge_length=26.0234610563)
    nd_4316558928 = nd_4316558800.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=114.058166915)
    nd_4316558736 = nd_4316558800.new_child(label="Node4316558736", taxon=None, edge_length=40.1699176918)
    nd_4316559184 = nd_4316558736.new_child(label="Node4316559184", taxon=None, edge_length=12.7558559915)
    nd_4316606544 = nd_4316558736.new_child(label="Node4316606544", taxon=None, edge_length=17.7722054357)
    nd_4316559312 = nd_4316559184.new_child(label="Node4316559312", taxon=None, edge_length=12.1691499629)
    nd_4316585744 = nd_4316559184.new_child(label="Node4316585744", taxon=None, edge_length=28.055060848)
    nd_4316559440 = nd_4316559312.new_child(label="Node4316559440", taxon=None, edge_length=4.4227398576)
    nd_4316583696 = nd_4316559312.new_child(label="Node4316583696", taxon=None, edge_length=3.26334313874)
    nd_4316559568 = nd_4316559440.new_child(label="Node4316559568", taxon=None, edge_length=5.98387013923)
    nd_4316561232 = nd_4316559440.new_child(label="Node4316561232", taxon=None, edge_length=17.9619739892)
    nd_4316559696 = nd_4316559568.new_child(label="Node4316559696", taxon=None, edge_length=7.53311752135)
    nd_4316560208 = nd_4316559568.new_child(label="Node4316560208", taxon=None, edge_length=4.18557870369)
    nd_4316559824 = nd_4316559696.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=31.0235157502)
    nd_4316559056 = nd_4316559696.new_child(label="Node4316559056", taxon=None, edge_length=12.7654788618)
    nd_4316560080 = nd_4316559056.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=18.2580368884)
    nd_4316559952 = nd_4316559056.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=18.2580368884)
    nd_4316560464 = nd_4316560208.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=34.3710545678)
    nd_4316560336 = nd_4316560208.new_child(label="Node4316560336", taxon=None, edge_length=7.81168349566)
    nd_4316560720 = nd_4316560336.new_child(label="Node4316560720", taxon=None, edge_length=15.350690271)
    nd_4316561104 = nd_4316560336.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=26.5593710722)
    nd_4316560848 = nd_4316560720.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=11.2086808012)
    nd_4316560592 = nd_4316560720.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=11.2086808012)
    nd_4316561360 = nd_4316561232.new_child(label="Node4316561360", taxon=None, edge_length=1.39777518086)
    nd_4316583312 = nd_4316561232.new_child(label="Node4316583312", taxon=None, edge_length=13.6008570384)
    nd_4316582032 = nd_4316561360.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=25.1807542407)
    nd_4316560976 = nd_4316561360.new_child(label="Node4316560976", taxon=None, edge_length=7.6242060025)
    nd_4316582288 = nd_4316560976.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=17.5565482382)
    nd_4316582160 = nd_4316560976.new_child(label="Node4316582160", taxon=None, edge_length=3.73213849687)
    nd_4316582544 = nd_4316582160.new_child(label="Node4316582544", taxon=None, edge_length=8.62088071739)
    nd_4316583056 = nd_4316582160.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=13.8244097414)
    nd_4316582672 = nd_4316582544.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=5.20352902397)
    nd_4316582416 = nd_4316582544.new_child(label="Node4316582416", taxon=None, edge_length=2.83199057731)
    nd_4316582928 = nd_4316582416.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.37153844665)
    nd_4316582800 = nd_4316582416.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=2.37153844665)
    nd_4316583440 = nd_4316583312.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=12.9776723832)
    nd_4316583184 = nd_4316583312.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=12.9776723832)
    nd_4316583824 = nd_4316583696.new_child(label="Node4316583824", taxon=None, edge_length=3.5026210108)
    nd_4316585104 = nd_4316583696.new_child(label="Node4316585104", taxon=None, edge_length=0.85376044557)
    nd_4316583952 = nd_4316583824.new_child(label="Node4316583952", taxon=None, edge_length=23.3327851176)
    nd_4316584336 = nd_4316583824.new_child(label="Node4316584336", taxon=None, edge_length=12.8697548268)
    nd_4316584080 = nd_4316583952.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=18.8644940012)
    nd_4316583568 = nd_4316583952.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=18.8644940012)
    nd_4316584464 = nd_4316584336.new_child(label="Node4316584464", taxon=None, edge_length=22.0000981488)
    nd_4316558288 = nd_4316584336.new_child(label="Node4316558288", taxon=None, edge_length=5.24037108992)
    nd_4316584592 = nd_4316584464.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=7.32742614319)
    nd_4316584208 = nd_4316584464.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=7.32742614319)
    nd_4316584848 = nd_4316558288.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=24.0871532021)
    nd_4316558224 = nd_4316558288.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=24.0871532021)
    nd_4316585232 = nd_4316585104.new_child(label="Node4316585232", taxon=None, edge_length=19.3585494438)
    nd_4316585616 = nd_4316585104.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=44.846139684)
    nd_4316585360 = nd_4316585232.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=25.4875902402)
    nd_4316584976 = nd_4316585232.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=25.4875902402)
    nd_4316585872 = nd_4316585744.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=33.0773323832)
    nd_4316585488 = nd_4316585744.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=33.0773323832)
    nd_4316606800 = nd_4316606544.new_child(label="Node4316606800", taxon=None, edge_length=13.6364666177)
    nd_4316607312 = nd_4316606544.new_child(label="Python regius", taxon=tax_4313760144, edge_length=56.1160437871)
    nd_4316606928 = nd_4316606800.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=42.4795771694)
    nd_4316606672 = nd_4316606800.new_child(label="Node4316606672", taxon=None, edge_length=16.6495052056)
    nd_4316607184 = nd_4316606672.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=25.8300719638)
    nd_4316607056 = nd_4316606672.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=25.8300719638)

    # set labels of nodes with taxa to taxon label, else oid (for consistent
    # identification in debugging)
    # for t in tree_list:
    #     t.assign_node_labels_from_taxon_or_oid()

    return tree_list

class TreeUnaryMetricsTest(unittest.TestCase):

    def testNBar(self):
        trees = _get_reference_tree_list()
        # trees = dendropy.TreeList.get_from_path(
        #         src=pathmap.tree_source_path("pythonidae.beast.mcmc.trees"),
        #         schema='nexus')
        expected_values = [
            7.818181818181818,
            7.515151515151516,
            7.666666666666667,
            8.727272727272727,
            8.757575757575758,
            8.636363636363637,
            8.727272727272727,
            8.757575757575758,
            8.727272727272727,
            8.727272727272727,
            8.575757575757576,
            ]
        for idx, tree in enumerate(trees):
            observed = treemeasure.N_bar(tree)
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def test_colless_tree_imbalance(self):
        trees = _get_reference_tree_list()
        # for tree in trees:
        #     print tree.colless_tree_imbalance()
        expected_values = [
            0.3024193548387097,
            0.2540322580645161,
            0.2762096774193548,
            0.3548387096774194,
            0.35685483870967744,
            0.344758064516129,
            0.3548387096774194,
            0.35685483870967744,
            0.3548387096774194,
            0.3548387096774194,
            0.3407258064516129,
            ]
        for idx, tree in enumerate(trees):
            observed = treemeasure.colless_tree_imbalance(tree, normalize="max")
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def test_colless_tree_imbalance2(self):
        # library(apTreeshape)
        # data(hivtree.treeshape)
        # print(paste("colless, raw: ", colless(hivtree.treeshape), sep=""))
        # print(paste("colless, pda: ", colless(hivtree.treeshape, "pda"), sep=""))
        # print(paste("colless, yule: ", colless(hivtree.treeshape, "yule"), sep=""))
        # # [1] "colless, raw: 992"
        # # [1] "colless, pda: 0.369977836654251"
        # # [1] "colless, yule: 0.993137704712054"
        tree = dendropy.Tree.get_from_path(
                src=pathmap.tree_source_path("hiv1.nexus"),
                schema='nexus')
        self.assertAlmostEqual(treemeasure.colless_tree_imbalance(tree, normalize=None), 992)
        self.assertAlmostEqual(treemeasure.colless_tree_imbalance(tree, normalize="pda"), 0.3699778366542512686443)
        self.assertAlmostEqual(treemeasure.colless_tree_imbalance(tree, normalize="yule"), 0.9931377047120540924041)

    def test_sackin_index(self):
        # library(apTreeshape)
        # data(hivtree.treeshape)
        # print(paste("sackin, raw: ", sackin(hivtree.treeshape), sep=""))
        # print(paste("sackin, pda: ", sackin(hivtree.treeshape, "pda"), sep=""))
        # print(paste("sackin, yule: ", sackin(hivtree.treeshape, "yule"), sep=""))
        # # [1] "sackin, raw: 2028"
        # # [1] "sackin, pda: 0.756365980579457"
        # # [1] "sackin, yule: 0.822783440343329"
        tree = dendropy.Tree.get_from_path(
                src=pathmap.tree_source_path("hiv1.nexus"),
                schema='nexus')
        self.assertAlmostEqual(treemeasure.sackin_index(tree, normalize=None), 2028)
        self.assertAlmostEqual(treemeasure.sackin_index(tree, normalize="pda"), 0.756365980579457)
        self.assertAlmostEqual(treemeasure.sackin_index(tree, normalize="yule"), 0.822783440343329)

    def test_b1(self):
        trees = _get_reference_tree_list()
        # for tree in trees:
        #     print tree.B1()
        expected_values = [
            18.686544011544008,
            16.862301587301587,
            18.012301587301586,
            15.803210678210679,
            15.803210678210679,
            16.219877344877347,
            15.80321067821068,
            15.80321067821068,
            15.803210678210679,
            15.80321067821068,
            16.10321067821068,
            ]
        for idx, tree in enumerate(trees):
            observed = treemeasure.B1(tree)
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def test_treeness(self):
        trees = _get_reference_tree_list()
        # for tree in trees:
        #     print tree.treeness()
        expected_values = [
            0.82043976304486,
            0.30678033634423607,
            0.2686940663128338,
            0.2674702980152253,
            0.2731856127080352,
            0.26942308963183575,
            0.2764640737121644,
            0.26096444220828763,
            0.2846852453916621,
            0.2791363657987356,
            0.28304948441090816,
            ]
        for idx, tree in enumerate(trees):
            observed = treemeasure.treeness(tree)
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def test_gamma2(self):
        trees = _get_reference_tree_list()
        # for tree in trees:
        #     print tree.pybus_harvey_gamma()
        expected_values = [
            6.690070011342222,
            -2.1016546214332665,
            -2.2071830302961493,
            -0.9868763184862083,
            -1.1223514055125514,
            -1.0914035287339103,
            -0.9432772103480326,
            -0.9855794349340775,
            -0.7566110136514949,
            -0.4693672063234924,
            0.08314644690264045,
            ]
        for idx, tree in enumerate(trees):
            observed = treemeasure.pybus_harvey_gamma(tree)
            expected = expected_values[idx]
            self.assertAlmostEqual(expected, observed)

    def testPHGamma(self):
        newick_str = "((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);"
        tree = dendropy.Tree.get_from_stream(StringIO(newick_str), schema="newick")
        g = treemeasure.pybus_harvey_gamma(tree)
        self.assertAlmostEqual(g, 0.546276, 4)

class TreeEuclideanDistTest(unittest.TestCase):

    def runTest(self):
         tree_list = dendropy.TreeList.get_from_stream(
            StringIO("""((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);
                        ((t5:2.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):1):0.065840,t3:0.170221):0.383247);
                        ((t5:0.161175,t6:0.161175):0.392293,((t2:0.075411,(t4:0.104381,t1:0.075411):1):0.065840,t3:0.170221):0.383247);
                        ((t5:0.161175,t6:0.161175):0.392293,((t4:0.104381,(t2:0.075411,t1:0.075411):0.028969):0.065840,t3:0.170221):0.383247);
                        """),
            schema="newick")
         for t in tree_list:
             t.encode_bipartitions()
         self.assertAlmostEqual(treecompare.euclidean_distance(tree_list[0], tree_list[1]), 2.0)
         self.assertAlmostEqual(treecompare.euclidean_distance(tree_list[0], tree_list[2]), math.sqrt(2.0))
         self.assertAlmostEqual(treecompare.euclidean_distance(tree_list[0], tree_list[3]), 0.97103099999999998)
         self.assertAlmostEqual(treecompare.euclidean_distance(tree_list[1], tree_list[2]), math.sqrt(6.0))
         self.assertAlmostEqual(treecompare.euclidean_distance(tree_list[1], tree_list[3]), 2.2232636377544162)
         self.assertAlmostEqual(treecompare.euclidean_distance(tree_list[2], tree_list[3]), 1.000419513484718)

class TreeSymmetricDistTest(unittest.TestCase):

    def runTest(self):
         ref = dendropy.Tree.get_from_stream(StringIO("((t5,t6),((t4,(t2,t1)),t3));"), schema="newick")
         taxon_namespace = ref.taxon_namespace
         ref.encode_bipartitions()
         o_tree = dendropy.Tree.get_from_stream(StringIO("((t1,t2),((t4,(t5,t6)),t3));"), schema="newick", taxon_namespace=taxon_namespace)
         o_tree.encode_bipartitions()
         self.assertEqual(treecompare.symmetric_difference(o_tree, ref), 2)

class TreeCompareTests(dendropytest.ExtendedTestCase):

    def setUp(self):
        tns = dendropy.TaxonNamespace()
        self.tree_list1 = _get_reference_tree_list(taxon_namespace=tns)
        self.tree_list2 = _get_reference_tree_list(taxon_namespace=tns)

    # def testNonMutatingDistinctTaxonNamespaceSameStructComparisons(self):
    #     tl1_ts = self.tree_list1.taxon_namespace
    #     tl2_ts = self.tree_list2.taxon_namespace
    #     self.assertIsNot(tl1_ts, tl2_ts)
    #     for i, t1 in enumerate(self.tree_list1):
    #         t2 = self.tree_list2[i]
    #         t1_ts = t1.taxon_namespace
    #         t2_ts = t2.taxon_namespace
    #         self.assertIsNot(t1_ts, t2_ts)
    #         self.assertEqual(t1.symmetric_difference(t2), 0)
    #         self.assertAlmostEqual(t1.euclidean_distance(t2), 0)
    #         self.assertAlmostEqual(t1.robinson_foulds_distance(t2), 0)
    #         self.assertIs(t1.taxon_namespace, t1_ts)
    #         self.assertIs(t2.taxon_namespace, t2_ts)
    #     self.assertIs(self.tree_list1.taxon_namespace, tl1_ts)
    #     self.assertIs(self.tree_list2.taxon_namespace, tl2_ts)

    def testSymmetricDifferences(self):
        expected = {
            (0,1):60, (0,2):60, (0,3):60, (0,4):60, (0,5):60, (0,6):60, (0,7):60, (0,8):60,
            (0,9):60, (0,10):60, (1,2):14, (1,3):24, (1,4):22, (1,5):20, (1,6):24, (1,7):22, (1,8):24,
            (1,9):22, (1,10):22, (2,3):18, (2,4):16, (2,5):16, (2,6):18, (2,7):16, (2,8):18, (2,9):16, (2,10):16,
            (3,4):4, (3,5):4, (3,6):0, (3,7):4, (3,8):0, (3,9):2, (3,10):4, (4,5):2, (4,6):4, (4,7):0, (4,8):4,
            (4,9):2, (4,10):4, (5,6):4, (5,7):2, (5,8):4, (5,9):2, (5,10):4, (6,7):4, (6,8):0, (6,9):2, (6,10):4,
            (7,8):4, (7,9):2, (7,10):4, (8,9):2, (8,10):4, (9,10):2,
        }
        for i, t1 in enumerate(self.tree_list1[:-1]):
            for j, t2 in enumerate(self.tree_list2[i+1:]):
                v = treecompare.symmetric_difference(t1, t2)
                self.assertEqual(expected[(i, i+j+1)], v)
#                print "(%d,%d):%d," % (i, i+j+1, v),
#                if (i * i+j+1) % 6 == 0:
#                    print

    def testEuclideanDistances(self):
        expected = {
            (0,1):442.518379997, (0,2):458.269219125, (0,3):492.707662859, (0,4):457.731995932, (0,5):463.419798784, (0,6):462.181969494,
            (0,7):439.865064545, (0,8):462.3054297, (0,9):479.06569226, (0,10):544.720324057, (1,2):105.534825723, (1,3):168.86739068, (1,4):119.287056085, (1,5):127.221894919, (1,6):125.918517173,
            (1,7):102.290062347, (1,8):130.5296198, (1,9):154.336066685, (1,10):247.555999428, (2,3):89.1098950842, (2,4):45.5124918081,
            (2,5):52.2607244547, (2,6):53.0477320261, (2,7):62.1391636266, (2,8):59.898883066, (2,9):79.3921379438, (2,10):172.187021923,
            (3,4):73.4046806483, (3,5):61.7211889655, (3,6):63.308525227,
            (3,7):113.043429355, (3,8):64.9098905352, (3,9):43.9926843558, (3,10):91.395044825, (4,5):22.881252195, (4,6):24.686671743,
            (4,7):47.14854215, (4,8):30.4425119229, (4,9):58.4893274048, (4,10):158.948156946, (5,6):24.7029660833, (5,7):56.9022982438, (5,8):25.0745838358, (5,9):45.9638357231, (5,10):146.364107049,
            (6,7):56.1301333366, (6,8):20.3469798051, (6,9):43.429825221, (6,10):145.712937469, (7,8):58.1647873304, (7,9):89.4537113125, (7,10):197.098347126, (8,9):40.5187846693, (8,10):145.393476072,
            (9,10):111.210401924,
        }
        for i, t1 in enumerate(self.tree_list1[:-1]):
            for j, t2 in enumerate(self.tree_list2[i+1:]):
                v = treecompare.euclidean_distance(t1, t2)
                self.assertAlmostEqual(expected[(i, i+j+1)], v)
#                print "(%d,%d):%s," % (i, i+j+1, v),
#                if (i * i+j+1) % 6 == 0:
#                    print

    def testRobinsonFouldsDistances(self):
        expected = {
            (0,1):1849.2928245, (0,2):2058.49072588, (0,3):2196.0995614, (0,4):1953.16064964, (0,5):1984.76411566, (0,6):1943.24487014,
            (0,7):1723.09194669, (0,8):1920.18504491, (0,9):1998.04696628, (0,10):2406.42091465, (1,2):508.212960297, (1,3):702.092000773, (1,4):579.45550447, (1,5):577.047914047, (1,6):596.881857714,
            (1,7):535.123132276, (1,8):611.28408319, (1,9):632.852687475, (1,10):857.759045631, (2,3):364.804588356, (2,4):283.907134148,
            (2,5):305.534136399, (2,6):318.128572842, (2,7):424.71989186, (2,8):351.751319705, (2,9):358.03680072, (2,10):531.731219604,
            (3,4):315.556017395, (3,5):271.016494089, (3,6):314.906668504,
            (3,7):517.444273417, (3,8):343.014958112, (3,9):266.498405531, (3,10):278.870282525, (4,5):133.642994808, (4,6):134.649689854,
            (4,7):260.010627711, (4,8):148.901405649, (4,9):187.954728978, (4,10):477.315325085, (5,6):150.970483084, (5,7):277.342311245, (5,8):144.704886539, (5,9):159.326519241, (5,10):449.629738145,
            (6,7):256.511541887, (6,8):119.487158128, (6,9):182.878241583, (6,10):493.201642403, (7,8):237.16728985, (7,9):296.353239488, (7,10):709.696300851, (8,9):171.021015022, (8,10):522.572965967,
            (9,10):435.439226227,
        }
        for i, t1 in enumerate(self.tree_list1[:-1]):
            for j, t2 in enumerate(self.tree_list2[i+1:]):
                v = treecompare.robinson_foulds_distance(t1, t2)
                self.assertAlmostEqual(expected[(i, i+j+1)], v)
#                print "(%d,%d):%s," % (i, i+j+1, v),
#                if (i * i+j+1) % 6 == 0:
#                    print

class FrequencyOfBipartitionsTests(unittest.TestCase):

    def testCount1(self):
        trees = dendropy.TreeList.get_from_path(
                src=pathmap.tree_source_path('pythonidae.random.bd0301.tre'),
                schema='nexus')
        bipartition_leaves = ['Python regius', 'Apodora papuana']
        f = trees.frequency_of_bipartition(labels=bipartition_leaves)
        self.assertAlmostEqual(f, 0.04)

    def testRaisesIndexError(self):
        trees = dendropy.TreeList.get_from_path(
                src=pathmap.tree_source_path('pythonidae.random.bd0301.tre'),
                schema='nexus')
        bipartition_leaves = ['Bad Taxon', 'Apodora papuana']
        self.assertRaises(IndexError, trees.frequency_of_bipartition, labels=bipartition_leaves)

    def test_freqs2(self):
        trees = dendropy.TreeList.get(
                path=pathmap.tree_source_path("pythonidae.mb.run1.t"),
                schema='nexus')
        test_sets = [
                # labels, split bitmask, frequency
                ( ["Python molurus", "Python regius"],                                24576, 0.00990099009901 ),
                ( ["Morelia clastolepis", "Morelia nauta", "Morelia kinghorni"],  536877056,   0.990099009901 ),
                ( ["Liasis olivaceus", "Apodora papuana"],                          4194432,   0.346534653465 ),
                ]
        for labels, split_bitmask, exp_freq in test_sets:
            self.assertAlmostEqual(trees.frequency_of_bipartition(labels=labels), exp_freq)
            taxa = trees.taxon_namespace.get_taxa(labels=labels)
            self.assertAlmostEqual(trees.frequency_of_bipartition(taxa=taxa), exp_freq)
            split_bitmask = trees.taxon_namespace.taxa_bitmask(labels=labels)
            self.assertEqual(split_bitmask, split_bitmask)
            self.assertAlmostEqual(trees.frequency_of_bipartition(split_bitmask=split_bitmask), exp_freq)
            bipartition = dendropy.Bipartition(bitmask=split_bitmask)
            self.assertAlmostEqual(trees.frequency_of_bipartition(bipartition=bipartition), exp_freq)

    def test_tree_bipartitions_encoding(self):
        trees = dendropy.TreeList.get(
                path=pathmap.tree_source_path("pythonidae.mb.run1.t"),
                schema='nexus')
        labels = ["Liasis olivaceus", "Apodora papuana"]
        freq = trees.frequency_of_bipartition(
                labels=labels,
                is_bipartitions_updated=True) # this will be ignore as the member `Tree.bipartition_encoding` are not populated
        self.assertAlmostEqual(freq, 0.346534653465)

if __name__ == "__main__":
    unittest.main()
