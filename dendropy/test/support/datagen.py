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

def reference_trees_filename(schema):
    if schema == "nexus":
        return "pythonidae.reference-trees.nexus"
    elif schema == "newick":
        return "pythonidae.reference-trees.newick"
    else:
        raise ValueError("Reference trees not available in '%s' format" % schema)

def reference_tree_list_postorder_node_labels():
    return [
        ['Antaresia childreni','Morelia clastolepis','Loxocemus bicolor','Node4313762448','Node4313762128','Antaresia melanocephalus','Antaresia perthensis','Node4313763216','Python timoriensis','Liasis albertisii','Node4313763600','Node4313763088','Python sebae','Antaresia maculosa','Node4313776208','Node4313762960','Python curtus','Bothrochilus boa','Node4313776912','Xenopeltis unicolor','Liasis fuscus','Node4313777296','Node4313776720','Node4313762704','Node4313762000','Morelia boeleni','Node4313761872','Morelia tracyae','Morelia viridisN','Python reticulatus','Node4313778192','Morelia bredli','Node4313777680','Node4313777808','Node4313761744','Liasis mackloti','Antaresia ramsayi','Node4313778832','Morelia kinghorni','Node4313778704','Node4313761616','Python regius','Candola aspera','Node4313779472','Antaresia stimsoni','Node4313779088','Node4313761488','Python molurus','Morelia amethistina','Node4313800784','Morelia viridisS','Liasis olivaceus','Node4313801168','Node4313780112','Morelia nauta','Node4313779728','Node4313761360','Apodora papuana','Morelia carinata','Node4313801552','Node4313761232','Morelia spilota','Morelia oenpelliensis','Node4313801936','Node4313761104'],
        ['Candola aspera','Morelia amethistina','Morelia spilota','Morelia oenpelliensis','Morelia bredli','Node4313803792','Node4313802960','Node4313803536','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4313804240','Node4313804560','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4313826128','Antaresia perthensis','Node4313825744','Node4313825616','Node4313804432','Node4313803408','Liasis fuscus','Liasis mackloti','Node4313827024','Apodora papuana','Liasis olivaceus','Node4313827408','Node4313826896','Antaresia melanocephalus','Antaresia ramsayi','Node4313827664','Node4313826768','Morelia boeleni','Morelia tracyae','Morelia kinghorni','Morelia nauta','Node4313828816','Morelia clastolepis','Node4313828432','Node4313828176','Node4313828048','Node4313826384','Node4313803280','Python timoriensis','Python reticulatus','Node4316155984','Bothrochilus boa','Liasis albertisii','Node4313802512','Node4313829328','Node4313803152','Xenopeltis unicolor','Python curtus','Python sebae','Python molurus','Node4316157008','Node4316157136','Python regius','Node4316156624','Node4316156752','Node4313803024','Node4313802704','Loxocemus bicolor','Node4313802576'],
        ['Candola aspera','Liasis fuscus','Liasis mackloti','Node4316159632','Apodora papuana','Liasis olivaceus','Node4316176464','Node4316159504','Antaresia melanocephalus','Antaresia ramsayi','Node4316176720','Node4316159376','Morelia boeleni','Morelia tracyae','Morelia kinghorni','Morelia nauta','Node4316177872','Morelia clastolepis','Node4316177488','Node4316177232','Node4316177104','Node4316159248','Morelia oenpelliensis','Morelia amethistina','Node4316178640','Morelia spilota','Morelia bredli','Node4316179024','Node4316178512','Bothrochilus boa','Liasis albertisii','Node4316179280','Node4316178128','Node4316159120','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316179664','Node4316179920','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316201296','Antaresia perthensis','Node4316158224','Node4316180432','Node4316179792','Node4316158992','Python timoriensis','Python reticulatus','Node4316201552','Node4316158864','Python curtus','Python sebae','Python molurus','Node4316202256','Node4316202384','Python regius','Node4316202128','Node4316158736','Loxocemus bicolor','Xenopeltis unicolor','Node4316203024','Node4316158608','Node4316158288'],
        ['Candola aspera','Loxocemus bicolor','Liasis fuscus','Liasis mackloti','Node4316225808','Apodora papuana','Liasis olivaceus','Node4316226192','Node4316225680','Antaresia melanocephalus','Antaresia ramsayi','Node4316226448','Node4316205008','Bothrochilus boa','Liasis albertisii','Node4316226768','Node4316204880','Morelia boeleni','Node4316204752','Morelia spilota','Morelia bredli','Node4316227728','Morelia oenpelliensis','Morelia tracyae','Morelia clastolepis','Morelia kinghorni','Node4316228880','Morelia nauta','Node4316228752','Morelia amethistina','Node4316228368','Node4316227984','Node4316228112','Node4316227600','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316203600','Node4316229520','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316251088','Antaresia perthensis','Node4316250704','Node4316250576','Node4316229392','Node4316227472','Node4316204624','Python timoriensis','Python reticulatus','Node4316251600','Node4316204496','Python curtus','Python sebae','Python molurus','Node4316251984','Node4316252112','Python regius','Node4316251856','Node4316203920','Node4316204112','Xenopeltis unicolor','Node4316203984','Node4316203664'],
        ['Candola aspera','Xenopeltis unicolor','Morelia boeleni','Liasis fuscus','Liasis mackloti','Node4316275408','Apodora papuana','Liasis olivaceus','Node4316275792','Node4316275280','Antaresia melanocephalus','Antaresia ramsayi','Node4316276048','Node4316275152','Bothrochilus boa','Liasis albertisii','Node4316276432','Node4316275024','Node4316254160','Morelia oenpelliensis','Morelia tracyae','Morelia clastolepis','Morelia kinghorni','Morelia nauta','Node4316253008','Node4316277712','Morelia amethistina','Node4316277584','Node4316277392','Morelia spilota','Morelia bredli','Node4316278352','Node4316276816','Node4316277072','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316299344','Node4316299536','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316300560','Antaresia perthensis','Node4316300176','Node4316300048','Node4316299408','Node4316276944','Node4316254032','Python timoriensis','Python reticulatus','Node4316301072','Node4316253904','Python curtus','Python sebae','Python molurus','Node4316301456','Node4316301584','Python regius','Node4316301328','Node4316253776','Loxocemus bicolor','Node4316253392','Node4316253456','Node4316253136'],
        ['Candola aspera','Loxocemus bicolor','Morelia boeleni','Antaresia melanocephalus','Antaresia ramsayi','Node4316324752','Liasis fuscus','Liasis mackloti','Node4316325264','Apodora papuana','Liasis olivaceus','Node4316325648','Node4316325136','Node4316324624','Bothrochilus boa','Liasis albertisii','Node4316326032','Node4316324496','Node4316324176','Morelia spilota','Morelia bredli','Node4316326672','Morelia oenpelliensis','Node4316326544','Morelia tracyae','Morelia nauta','Morelia kinghorni','Node4316327696','Morelia clastolepis','Node4316327568','Morelia amethistina','Node4316326928','Node4316327184','Node4316326416','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316344400','Node4316344912','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316345936','Antaresia perthensis','Node4316345552','Node4316345424','Node4316344784','Node4316326288','Node4316324048','Python timoriensis','Python reticulatus','Node4316346448','Node4316323920','Python curtus','Python sebae','Python molurus','Node4316346832','Node4316346960','Python regius','Node4316346704','Node4316302864','Node4316303056','Xenopeltis unicolor','Node4316302928','Node4316302608'],
        ['Candola aspera','Morelia spilota','Morelia bredli','Node4316369744','Morelia oenpelliensis','Morelia tracyae','Morelia nauta','Morelia kinghorni','Morelia clastolepis','Node4316370640','Node4316347856','Morelia amethistina','Node4316370384','Node4316370000','Node4316370128','Node4316369616','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316371408','Node4316371536','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316372496','Antaresia perthensis','Node4316372176','Node4316372048','Node4316371280','Node4316369488','Morelia boeleni','Antaresia melanocephalus','Antaresia ramsayi','Node4316394064','Liasis fuscus','Liasis mackloti','Node4316394640','Apodora papuana','Liasis olivaceus','Node4316395024','Node4316394512','Node4316393936','Bothrochilus boa','Liasis albertisii','Node4316395408','Node4316393616','Node4316372752','Node4316369360','Python timoriensis','Python reticulatus','Node4316395664','Node4316369232','Python curtus','Python sebae','Python molurus','Node4316396176','Node4316396304','Python regius','Node4316396048','Node4316369104','Loxocemus bicolor','Node4316368976','Xenopeltis unicolor','Node4316348304','Node4316347984'],
        ['Candola aspera','Xenopeltis unicolor','Loxocemus bicolor','Morelia spilota','Morelia bredli','Node4316419472','Morelia tracyae','Morelia clastolepis','Morelia nauta','Morelia kinghorni','Node4316420112','Node4316420240','Morelia amethistina','Node4316419728','Node4316419856','Node4316419344','Morelia oenpelliensis','Node4316419216','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316420880','Node4316421264','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316442832','Antaresia perthensis','Node4316421904','Node4316421776','Node4316421136','Node4316397392','Morelia boeleni','Antaresia melanocephalus','Antaresia ramsayi','Node4316443856','Liasis fuscus','Liasis mackloti','Node4316444368','Apodora papuana','Liasis olivaceus','Node4316444752','Node4316444240','Node4316443728','Bothrochilus boa','Liasis albertisii','Node4316445136','Node4316443344','Node4316443088','Node4316419088','Python timoriensis','Python reticulatus','Node4316445392','Node4316418960','Python curtus','Python sebae','Python molurus','Node4316445904','Node4316446032','Python regius','Node4316445776','Node4316418576','Node4316418256','Node4316418320','Node4316397456'],
        ['Candola aspera','Morelia boeleni','Bothrochilus boa','Liasis albertisii','Node4316468688','Antaresia melanocephalus','Antaresia ramsayi','Node4316469200','Liasis fuscus','Liasis mackloti','Node4316469712','Apodora papuana','Liasis olivaceus','Node4316470096','Node4316469584','Node4316469072','Node4316467408','Node4316468304','Morelia tracyae','Morelia nauta','Morelia clastolepis','Morelia kinghorni','Node4316471120','Node4316471248','Morelia amethistina','Node4316470352','Node4316470864','Morelia oenpelliensis','Node4316470736','Morelia spilota','Morelia bredli','Node4316488592','Node4316470608','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316488976','Node4316489104','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316490128','Antaresia perthensis','Node4316489744','Node4316489616','Node4316488848','Node4316470480','Node4316468176','Python timoriensis','Python reticulatus','Node4316490576','Node4316468048','Python curtus','Python sebae','Python molurus','Node4316491024','Node4316491152','Python regius','Node4316490896','Node4316467920','Loxocemus bicolor','Node4316467792','Xenopeltis unicolor','Node4316467344','Node4316467472'],
        ['Candola aspera','Loxocemus bicolor','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316513552','Node4316514192','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316515216','Antaresia perthensis','Node4316514832','Node4316514704','Node4316514064','Morelia spilota','Morelia bredli','Node4316515856','Morelia tracyae','Morelia clastolepis','Morelia kinghorni','Morelia nauta','Node4316533072','Node4316533200','Morelia amethistina','Node4316516112','Node4316532816','Morelia oenpelliensis','Node4316516240','Node4316515728','Node4316513936','Morelia boeleni','Bothrochilus boa','Liasis albertisii','Node4316534480','Antaresia melanocephalus','Antaresia ramsayi','Node4316534992','Liasis fuscus','Liasis mackloti','Node4316535504','Apodora papuana','Liasis olivaceus','Node4316535888','Node4316535376','Node4316534864','Node4316533968','Node4316534096','Node4316513808','Python timoriensis','Python reticulatus','Node4316512784','Node4316513680','Python curtus','Python sebae','Python molurus','Node4316557520','Node4316536656','Python regius','Node4316536464','Node4316513104','Node4316513296','Xenopeltis unicolor','Node4316513168','Node4316512848'],
        ['Xenopeltis unicolor','Loxocemus bicolor','Morelia carinata','Morelia viridisS','Morelia viridisN','Node4316559056','Node4316559696','Antaresia maculosa','Antaresia childreni','Antaresia stimsoni','Node4316560720','Antaresia perthensis','Node4316560336','Node4316560208','Node4316559568','Morelia oenpelliensis','Morelia tracyae','Morelia clastolepis','Morelia kinghorni','Morelia nauta','Node4316582416','Node4316582544','Morelia amethistina','Node4316582160','Node4316560976','Node4316561360','Morelia spilota','Morelia bredli','Node4316583312','Node4316561232','Node4316559440','Antaresia melanocephalus','Antaresia ramsayi','Node4316583952','Liasis fuscus','Liasis mackloti','Node4316584464','Apodora papuana','Liasis olivaceus','Node4316558288','Node4316584336','Node4316583824','Bothrochilus boa','Liasis albertisii','Node4316585232','Morelia boeleni','Node4316585104','Node4316583696','Node4316559312','Python timoriensis','Python reticulatus','Node4316585744','Node4316559184','Python curtus','Python sebae','Python molurus','Node4316606672','Node4316606800','Python regius','Node4316606544','Node4316558736','Node4316558800','Node4316558480','Candola aspera','Node4316558352'],
    ]

def reference_tree_list_newick_string(taxon_set=None):
    return """\
    ((((((((('Antaresia childreni':3.82103266723,('Morelia clastolepis':1.09809398991,'Loxocemus bicolor':1.09809398991)Node4313762448:2.72293867732)Node4313762128:10.4795105126,(((('Antaresia melanocephalus':1.47649847948,'Antaresia perthensis':1.47649847948)Node4313763216:3.56964311758,('Python timoriensis':1.49653696768,'Liasis albertisii':1.49653696768)Node4313763600:3.54960462938)Node4313763088:1.10594701364,('Python sebae':5.20555672012,'Antaresia maculosa':5.20555672012)Node4313776208:0.946531890581)Node4313762960:4.47413907662,(('Python curtus':4.46599929505,'Bothrochilus boa':4.46599929505)Node4313776912:2.88197852526,('Xenopeltis unicolor':0.483824039664,'Liasis fuscus':0.483824039664)Node4313777296:6.86415378064)Node4313776720:3.27824986702)Node4313762704:3.67431549248)Node4313762000:4.00932181352,'Morelia boeleni':18.3098649933)Node4313761872:0.739270951321,('Morelia tracyae':2.17652327312,(('Morelia viridisN':0.0125016196671,'Python reticulatus':0.0125016196671)Node4313778192:0.491713738136,'Morelia bredli':0.504215357803)Node4313777680:1.67230791531)Node4313777808:16.8726126715)Node4313761744:1.29166787247,(('Liasis mackloti':4.19791393809,'Antaresia ramsayi':4.19791393809)Node4313778832:2.98661848623,'Morelia kinghorni':7.18453242432)Node4313778704:13.1562713928)Node4313761616:4.50960412336,(('Python regius':0.209207825685,'Candola aspera':0.209207825685)Node4313779472:0.207889736001,'Antaresia stimsoni':0.417097561686)Node4313779088:24.4333103788)Node4313761488:53.7887032791,((('Python molurus':5.1470935265,'Morelia amethistina':5.1470935265)Node4313800784:2.6364754365,('Morelia viridisS':0.657835681585,'Liasis olivaceus':0.657835681585)Node4313801168:7.12573328141)Node4313780112:1.24643505521,'Morelia nauta':9.03000401821)Node4313779728:69.6091072013)Node4313761360:123.936295332,('Apodora papuana':0.152425796315,'Morelia carinata':0.152425796315)Node4313801552:202.422980755)Node4313761232:78.6266408419,('Morelia spilota':51.3217384582,'Morelia oenpelliensis':51.3217384582)Node4313801936:229.880308935)Node4313761104;
    (('Candola aspera':81.2974222246,((((('Morelia amethistina':21.0901008779,('Morelia spilota':17.7035654613,('Morelia oenpelliensis':15.0790936884,'Morelia bredli':15.0790936884)Node4313803792:2.6244717729)Node4313802960:3.38653541663)Node4313803536:15.1180336863,(('Morelia carinata':27.7948635409,('Morelia viridisS':19.4686760786,'Morelia viridisN':19.4686760786)Node4313804240:8.32618746237)Node4313804560:5.00375936144,('Antaresia maculosa':26.9874255418,(('Antaresia childreni':8.07749085501,'Antaresia stimsoni':8.07749085501)Node4313826128:16.6530983052,'Antaresia perthensis':24.7305891602)Node4313825744:2.25683638168)Node4313825616:5.81119736053)Node4313804432:3.40951166184)Node4313803408:4.5936111676,(((('Liasis fuscus':5.68389243709,'Liasis mackloti':5.68389243709)Node4313827024:26.1365403684,('Apodora papuana':19.714025971,'Liasis olivaceus':19.714025971)Node4313827408:12.1064068345)Node4313826896:3.37363071467,('Antaresia melanocephalus':18.8177870609,'Antaresia ramsayi':18.8177870609)Node4313827664:16.3762764593)Node4313826768:4.33214615343,('Morelia boeleni':28.5609164144,('Morelia tracyae':17.1692672846,(('Morelia kinghorni':5.2450565205,'Morelia nauta':5.2450565205)Node4313828816:1.44575855569,'Morelia clastolepis':6.69081507619)Node4313828432:10.4784522084)Node4313828176:11.3916491298)Node4313828048:10.9652932592)Node4313826384:1.27553605821)Node4313803280:0.76583222117,(('Python timoriensis':15.1414324532,'Python reticulatus':15.1414324532)Node4316155984:23.8923728386,('Bothrochilus boa':16.4481972995,'Liasis albertisii':16.4481972995)Node4313802512:22.5856079922)Node4313829328:2.53377266123)Node4313803152:6.57324583185,('Xenopeltis unicolor':31.8813718333,(('Python curtus':16.9648006514,('Python sebae':12.3010608596,'Python molurus':12.3010608596)Node4316157008:4.66373979181)Node4316157136:8.94343133576,'Python regius':25.9082319872)Node4316156624:5.97313984611)Node4316156752:16.2594519516)Node4313803024:33.1565984398)Node4313802704:18.8917197007,'Loxocemus bicolor':100.189141925)Node4313802576;
    ('Candola aspera':109.372663833,((((((((('Liasis fuscus':10.9943371535,'Liasis mackloti':10.9943371535)Node4316159632:18.9362130146,('Apodora papuana':21.3765394762,'Liasis olivaceus':21.3765394762)Node4316176464:8.55401069191)Node4316159504:4.82772953939,('Antaresia melanocephalus':16.9627400103,'Antaresia ramsayi':16.9627400103)Node4316176720:17.7955396972)Node4316159376:6.39175712039,('Morelia boeleni':32.036740867,('Morelia tracyae':17.2575489744,(('Morelia kinghorni':3.4529882291,'Morelia nauta':3.4529882291)Node4316177872:0.439451186875,'Morelia clastolepis':3.89243941597)Node4316177488:13.3651095585)Node4316177232:14.7791918926)Node4316177104:9.11329596086)Node4316159248:2.41802497137,((('Morelia oenpelliensis':25.0022602166,'Morelia amethistina':25.0022602166)Node4316178640:1.13984981318,('Morelia spilota':12.987805351,'Morelia bredli':12.987805351)Node4316179024:13.1543046788)Node4316178512:11.5242004723,('Bothrochilus boa':24.274725326,'Liasis albertisii':24.274725326)Node4316179280:13.3915851761)Node4316178128:5.90175129715)Node4316159120:2.53965409687,(('Morelia carinata':30.4665100305,('Morelia viridisS':19.4621901768,'Morelia viridisN':19.4621901768)Node4316179664:11.0043198537)Node4316179920:7.8146690322,('Antaresia maculosa':33.1727582851,(('Antaresia childreni':9.02781095195,'Antaresia stimsoni':9.02781095195)Node4316201296:19.4308450954,'Antaresia perthensis':28.4586560473)Node4316158224:4.7141022378)Node4316180432:5.10842077756)Node4316179792:7.82653683343)Node4316158992:6.47536868825,('Python timoriensis':30.3039115511,'Python reticulatus':30.3039115511)Node4316201552:22.2791730332)Node4316158864:8.2214192007,(('Python curtus':37.2154940833,('Python sebae':21.9347777454,'Python molurus':21.9347777454)Node4316202256:15.2807163378)Node4316202384:6.14068810478,'Python regius':43.356182188)Node4316202128:17.448321597)Node4316158736:25.5450384687,('Loxocemus bicolor':82.94710604,'Xenopeltis unicolor':82.94710604)Node4316203024:3.40243621372)Node4316158608:23.0231215792)Node4316158288;
    ('Candola aspera':157.750076773,(('Loxocemus bicolor':93.189940239,(((((((('Liasis fuscus':8.15449187544,'Liasis mackloti':8.15449187544)Node4316225808:24.6792964706,('Apodora papuana':24.3044203289,'Liasis olivaceus':24.3044203289)Node4316226192:8.52936801714)Node4316225680:6.3930928479,('Antaresia melanocephalus':15.9864414111,'Antaresia ramsayi':15.9864414111)Node4316226448:23.2404397828)Node4316205008:2.31767871982,('Bothrochilus boa':23.2898123636,'Liasis albertisii':23.2898123636)Node4316226768:18.2547475502)Node4316204880:1.34047834167,'Morelia boeleni':42.8850382554)Node4316204752:2.88538115771,((('Morelia spilota':16.278650612,'Morelia bredli':16.278650612)Node4316227728:12.5474231006,('Morelia oenpelliensis':27.5592919578,('Morelia tracyae':14.2553766995,((('Morelia clastolepis':3.46866488962,'Morelia kinghorni':3.46866488962)Node4316228880:0.230110019205,'Morelia nauta':3.69877490882)Node4316228752:8.89489387836,'Morelia amethistina':12.5936687872)Node4316228368:1.66170791236)Node4316227984:13.3039152583)Node4316228112:1.26678175478)Node4316227600:14.2175774566,(('Morelia carinata':33.2267204786,('Morelia viridisS':19.7961005329,'Morelia viridisN':19.7961005329)Node4316203600:13.4306199458)Node4316229520:5.88218975316,('Antaresia maculosa':31.9976247603,(('Antaresia childreni':11.3547051232,'Antaresia stimsoni':11.3547051232)Node4316251088:15.9285543528,'Antaresia perthensis':27.2832594761)Node4316250704:4.71436528425)Node4316250576:7.11128547149)Node4316229392:3.9347409374)Node4316227472:2.72676824396)Node4316204624:9.74044751021,('Python timoriensis':33.3905850116,'Python reticulatus':33.3905850116)Node4316251600:22.1202819118)Node4316204496:10.3875398007,(('Python curtus':34.8181553863,('Python sebae':26.9223231435,'Python molurus':26.9223231435)Node4316251984:7.89583224277)Node4316252112:11.1907198563,'Python regius':46.0088752426)Node4316251856:19.8895314814)Node4316203920:27.291533515)Node4316204112:19.5811677101,'Xenopeltis unicolor':112.771107949)Node4316203984:44.9789688242)Node4316203664;
    ('Candola aspera':117.405482731,('Xenopeltis unicolor':90.5537014032,((((('Morelia boeleni':36.5498484742,(((('Liasis fuscus':6.29296152027,'Liasis mackloti':6.29296152027)Node4316275408:18.7911640437,('Apodora papuana':20.4324667004,'Liasis olivaceus':20.4324667004)Node4316275792:4.65165886356)Node4316275280:4.14502379891,('Antaresia melanocephalus':11.7665664768,'Antaresia ramsayi':11.7665664768)Node4316276048:17.4625828861)Node4316275152:6.27408700912,('Bothrochilus boa':20.0957025946,'Liasis albertisii':20.0957025946)Node4316276432:15.4075337774)Node4316275024:1.04661210215)Node4316254160:1.85501747057,(('Morelia oenpelliensis':24.8451626811,(('Morelia tracyae':14.3265271797,(('Morelia clastolepis':5.28574540958,('Morelia kinghorni':2.73160988664,'Morelia nauta':2.73160988664)Node4316253008:2.55413552294)Node4316277712:5.65307392619,'Morelia amethistina':10.9388193358)Node4316277584:3.38770784397)Node4316277392:7.72979171285,('Morelia spilota':9.57952070277,'Morelia bredli':9.57952070277)Node4316278352:12.4767981898)Node4316276816:2.78884378851)Node4316277072:10.9999278199,(('Morelia carinata':24.4395100584,('Morelia viridisS':16.6979256908,'Morelia viridisN':16.6979256908)Node4316299344:7.74158436759)Node4316299536:6.62946104565,('Antaresia maculosa':27.9986606508,(('Antaresia childreni':6.91002117518,'Antaresia stimsoni':6.91002117518)Node4316300560:16.2064963991,'Antaresia perthensis':23.1165175743)Node4316300176:4.88214307652)Node4316300048:3.07031045323)Node4316299408:4.77611939702)Node4316276944:2.5597754437)Node4316254032:12.3993667148,('Python timoriensis':27.4972972887,'Python reticulatus':27.4972972887)Node4316301072:23.3069353709)Node4316253904:6.56321044476,(('Python curtus':33.0996104792,('Python sebae':18.8062654172,'Python molurus':18.8062654172)Node4316301456:14.293345062)Node4316301584:6.08258726043,'Python regius':39.1821977396)Node4316301328:18.1852453647)Node4316253776:16.6353899137,'Loxocemus bicolor':74.0028330181)Node4316253392:16.5508683851)Node4316253456:26.8517813278)Node4316253136;
    ('Candola aspera':126.159749307,(('Loxocemus bicolor':76.2142576718,(((('Morelia boeleni':36.5119154213,((('Antaresia melanocephalus':15.1074329521,'Antaresia ramsayi':15.1074329521)Node4316324752:17.772328432,(('Liasis fuscus':7.5711347092,'Liasis mackloti':7.5711347092)Node4316325264:20.8200876951,('Apodora papuana':24.0183292126,'Liasis olivaceus':24.0183292126)Node4316325648:4.37289319177)Node4316325136:4.4885389798)Node4316324624:2.24740276648,('Bothrochilus boa':22.1472140662,'Liasis albertisii':22.1472140662)Node4316326032:12.9799500845)Node4316324496:1.38475127065)Node4316324176:1.78294733026,(((('Morelia spilota':9.34118829512,'Morelia bredli':9.34118829512)Node4316326672:11.7610480778,'Morelia oenpelliensis':21.1022363729)Node4316326544:1.62908011033,('Morelia tracyae':16.5271748204,((('Morelia nauta':2.44563718163,'Morelia kinghorni':2.44563718163)Node4316327696:1.89524872622,'Morelia clastolepis':4.34088590785)Node4316327568:7.91289243511,'Morelia amethistina':12.253778343)Node4316326928:4.27339647744)Node4316327184:6.20414166284)Node4316326416:12.2960670728,(('Morelia carinata':25.8119331435,('Morelia viridisS':14.3781644505,'Morelia viridisN':14.3781644505)Node4316344400:11.4337686931)Node4316344912:5.03710806168,('Antaresia maculosa':26.7514685092,(('Antaresia childreni':6.97141003814,'Antaresia stimsoni':6.97141003814)Node4316345936:17.4436719435,'Antaresia perthensis':24.4150819817)Node4316345552:2.33638652753)Node4316345424:4.09757269601)Node4316344784:4.17834235089)Node4316326288:3.26747919544)Node4316324048:8.24675869634,('Python timoriensis':22.9852293488,'Python reticulatus':22.9852293488)Node4316346448:23.556392099)Node4316323920:7.33856679933,(('Python curtus':29.2153611098,('Python sebae':21.0555716576,'Python molurus':21.0555716576)Node4316346832:8.15978945225)Node4316346960:6.75664224871,'Python regius':35.9720033585)Node4316346704:17.9081848887)Node4316302864:22.3340694246)Node4316303056:20.4272017262,'Xenopeltis unicolor':96.6414593981)Node4316302928:29.5182899091)Node4316302608;
    ('Candola aspera':124.564186516,(((((((('Morelia spilota':11.3186128368,'Morelia bredli':11.3186128368)Node4316369744:12.0191169658,('Morelia oenpelliensis':21.0803560094,('Morelia tracyae':14.1501582185,(('Morelia nauta':6.11730072991,('Morelia kinghorni':4.04355876922,'Morelia clastolepis':4.04355876922)Node4316370640:2.07374196069)Node4316347856:4.34576715349,'Morelia amethistina':10.4630678834)Node4316370384:3.68709033509)Node4316370000:6.93019779087)Node4316370128:2.25737379321)Node4316369616:12.5842551962,(('Morelia carinata':22.4623407504,('Morelia viridisS':16.5498725211,'Morelia viridisN':16.5498725211)Node4316371408:5.91246822929)Node4316371536:5.79725116231,('Antaresia maculosa':24.5414767549,(('Antaresia childreni':7.30910889905,'Antaresia stimsoni':7.30910889905)Node4316372496:11.551508813,'Antaresia perthensis':18.8606177121)Node4316372176:5.68085904285)Node4316372048:3.71811515781)Node4316371280:7.66239308607)Node4316369488:2.76745367097,('Morelia boeleni':37.5960808765,((('Antaresia melanocephalus':13.4678814627,'Antaresia ramsayi':13.4678814627)Node4316394064:19.1660901297,(('Liasis fuscus':6.08619858554,'Liasis mackloti':6.08619858554)Node4316394640:18.9933631628,('Apodora papuana':19.3061701007,'Liasis olivaceus':19.3061701007)Node4316395024:5.77339164759)Node4316394512:7.55440984409)Node4316393936:2.22916081797,('Bothrochilus boa':21.6486618625,'Liasis albertisii':21.6486618625)Node4316395408:13.2144705479)Node4316393616:2.73294846613)Node4316372752:1.09335779327)Node4316369360:7.28944403695,('Python timoriensis':25.9434669707,'Python reticulatus':25.9434669707)Node4316395664:20.035415736)Node4316369232:10.7040090023,(('Python curtus':30.462733423,('Python sebae':16.0194742188,'Python molurus':16.0194742188)Node4316396176:14.4432592042)Node4316396304:7.07244422644,'Python regius':37.5351776494)Node4316396048:19.1477140595)Node4316369104:20.3346663059,'Loxocemus bicolor':77.0175580149)Node4316368976:11.1789504571,'Xenopeltis unicolor':88.196508472)Node4316348304:36.3676780441)Node4316347984;
    ('Candola aspera':95.8502441646,('Xenopeltis unicolor':73.9760796713,('Loxocemus bicolor':64.4465636894,((((((('Morelia spilota':9.02720708579,'Morelia bredli':9.02720708579)Node4316419472:8.29264517113,('Morelia tracyae':10.6182161258,(('Morelia clastolepis':2.86717460067,('Morelia nauta':2.49595908021,'Morelia kinghorni':2.49595908021)Node4316420112:0.371215520464)Node4316420240:4.27223311967,'Morelia amethistina':7.13940772034)Node4316419728:3.47880840545)Node4316419856:6.70163613113)Node4316419344:2.42060556338,'Morelia oenpelliensis':19.7404578203)Node4316419216:7.52518318022,(('Morelia carinata':21.6260059492,('Morelia viridisS':12.4726660393,'Morelia viridisN':12.4726660393)Node4316420880:9.1533399099)Node4316421264:3.24355281662,('Antaresia maculosa':21.6668574015,(('Antaresia childreni':5.8075137439,'Antaresia stimsoni':5.8075137439)Node4316442832:11.5233229214,'Antaresia perthensis':17.3308366653)Node4316421904:4.33602073619)Node4316421776:3.2027013644)Node4316421136:2.39608223465)Node4316397392:1.03257897787,('Morelia boeleni':27.0825353379,((('Antaresia melanocephalus':9.50619301624,'Antaresia ramsayi':9.50619301624)Node4316443856:13.5545767859,(('Liasis fuscus':5.48728159904,'Liasis mackloti':5.48728159904)Node4316444368:12.8995814401,('Apodora papuana':16.9983690986,'Liasis olivaceus':16.9983690986)Node4316444752:1.38849394051)Node4316444240:4.67390676307)Node4316443728:2.64946637554,('Bothrochilus boa':14.3052158982,'Liasis albertisii':14.3052158982)Node4316445136:11.4050202795)Node4316443344:1.37229916019)Node4316443088:1.21568464049)Node4316419088:6.7231540226,('Python timoriensis':19.8342036742,'Python reticulatus':19.8342036742)Node4316445392:15.1871703268)Node4316418960:6.43697452247,(('Python curtus':23.7257571548,('Python sebae':16.0945691729,'Python molurus':16.0945691729)Node4316445904:7.63118798191)Node4316446032:3.93234771773,'Python regius':27.6581048725)Node4316445776:13.8002436509)Node4316418576:22.9882151659)Node4316418256:9.52951598189)Node4316418320:21.8741644934)Node4316397456;
    ('Candola aspera':126.147943419,(((((('Morelia boeleni':34.6878613533,(('Bothrochilus boa':21.8829103058,'Liasis albertisii':21.8829103058)Node4316468688:10.5685372807,(('Antaresia melanocephalus':13.1588228725,'Antaresia ramsayi':13.1588228725)Node4316469200:17.1757724208,(('Liasis fuscus':6.80247824391,'Liasis mackloti':6.80247824391)Node4316469712:19.5683902852,('Apodora papuana':21.5430118322,'Liasis olivaceus':21.5430118322)Node4316470096:4.82785669688)Node4316469584:3.96372676423)Node4316469072:2.11685229318)Node4316467408:2.23641376685)Node4316468304:0.617658990864,(((('Morelia tracyae':13.5020774756,(('Morelia nauta':3.78994236341,('Morelia clastolepis':2.53789923993,'Morelia kinghorni':2.53789923993)Node4316471120:1.25204312348)Node4316471248:6.20487512451,'Morelia amethistina':9.99481748791)Node4316470352:3.5072599877)Node4316470864:7.79402846351,'Morelia oenpelliensis':21.2961059391)Node4316470736:2.31259127041,('Morelia spilota':9.62601876119,'Morelia bredli':9.62601876119)Node4316488592:13.9826784484)Node4316470608:10.3401216686,(('Morelia carinata':22.5172124711,('Morelia viridisS':10.6110289453,'Morelia viridisN':10.6110289453)Node4316488976:11.9061835258)Node4316489104:4.39672345568,('Antaresia maculosa':22.9661472676,(('Antaresia childreni':7.25380762314,'Antaresia stimsoni':7.25380762314)Node4316490128:11.9409692012,'Antaresia perthensis':19.1947768244)Node4316489744:3.7713704432)Node4316489616:3.94778865925)Node4316488848:7.03488295134)Node4316470480:1.35670146604)Node4316468176:13.7010200713,('Python timoriensis':24.3086862402,'Python reticulatus':24.3086862402)Node4316490576:24.6978541753)Node4316468048:5.77299961102,(('Python curtus':27.8030546577,('Python sebae':17.3634777999,'Python molurus':17.3634777999)Node4316491024:10.4395768579)Node4316491152:5.80983351578,'Python regius':33.6128881735)Node4316490896:21.166651853)Node4316467920:17.1918011574,'Loxocemus bicolor':71.971341184)Node4316467792:13.6608326978,'Xenopeltis unicolor':85.6321738818)Node4316467344:40.5157695372)Node4316467472;
    ('Candola aspera':146.770054852,(('Loxocemus bicolor':83.313555189,(((((('Morelia carinata':25.1454154452,('Morelia viridisS':14.566531775,'Morelia viridisN':14.566531775)Node4316513552:10.5788836702)Node4316514192:4.82232736143,('Antaresia maculosa':25.387694547,(('Antaresia childreni':6.91975332756,'Antaresia stimsoni':6.91975332756)Node4316515216:12.7115868187,'Antaresia perthensis':19.6313401463)Node4316514832:5.75635440071)Node4316514704:4.58004825968)Node4316514064:3.25905094181,(('Morelia spilota':8.71082523711,'Morelia bredli':8.71082523711)Node4316515856:14.4811675935,(('Morelia tracyae':14.4468651064,(('Morelia clastolepis':5.18823967855,('Morelia kinghorni':3.47759210261,'Morelia nauta':3.47759210261)Node4316533072:1.71064757594)Node4316533200:6.98040063458,'Morelia amethistina':12.1686403131)Node4316516112:2.27822479328)Node4316532816:7.95009719313,'Morelia oenpelliensis':22.3969622995)Node4316516240:0.795030531054)Node4316515728:10.0348009179)Node4316513936:1.94845077373,('Morelia boeleni':33.9681920922,(('Bothrochilus boa':20.2312473989,'Liasis albertisii':20.2312473989)Node4316534480:12.4520799939,(('Antaresia melanocephalus':11.9647468446,'Antaresia ramsayi':11.9647468446)Node4316534992:19.0383478987,(('Liasis fuscus':9.12565806357,'Liasis mackloti':9.12565806357)Node4316535504:17.1180008393,('Apodora papuana':23.561905752,'Liasis olivaceus':23.561905752)Node4316535888:2.6817531508)Node4316535376:4.75943584051)Node4316534864:1.68023264943)Node4316533968:1.28486469944)Node4316534096:1.20705242999)Node4316513808:9.51023675819,('Python timoriensis':19.2959696748,'Python reticulatus':19.2959696748)Node4316512784:25.3895116055)Node4316513680:8.15039409252,(('Python curtus':31.4380841777,('Python sebae':17.8823541984,'Python molurus':17.8823541984)Node4316557520:13.5557299793)Node4316536656:6.97033769,'Python regius':38.4084218677)Node4316536464:14.4274535052)Node4316513104:30.4776798161)Node4316513296:13.4634525107,'Xenopeltis unicolor':96.7770076997)Node4316513168:49.9930471528)Node4316512848;
    (('Xenopeltis unicolor':140.081627971,('Loxocemus bicolor':114.058166915,(((((('Morelia carinata':31.0235157502,('Morelia viridisS':18.2580368884,'Morelia viridisN':18.2580368884)Node4316559056:12.7654788618)Node4316559696:7.53311752135,('Antaresia maculosa':34.3710545678,(('Antaresia childreni':11.2086808012,'Antaresia stimsoni':11.2086808012)Node4316560720:15.350690271,'Antaresia perthensis':26.5593710722)Node4316560336:7.81168349566)Node4316560208:4.18557870369)Node4316559568:5.98387013923,(('Morelia oenpelliensis':25.1807542407,('Morelia tracyae':17.5565482382,(('Morelia clastolepis':5.20352902397,('Morelia kinghorni':2.37153844665,'Morelia nauta':2.37153844665)Node4316582416:2.83199057731)Node4316582544:8.62088071739,'Morelia amethistina':13.8244097414)Node4316582160:3.73213849687)Node4316560976:7.6242060025)Node4316561360:1.39777518086,('Morelia spilota':12.9776723832,'Morelia bredli':12.9776723832)Node4316583312:13.6008570384)Node4316561232:17.9619739892)Node4316559440:4.4227398576,((('Antaresia melanocephalus':18.8644940012,'Antaresia ramsayi':18.8644940012)Node4316583952:23.3327851176,(('Liasis fuscus':7.32742614319,'Liasis mackloti':7.32742614319)Node4316584464:22.0000981488,('Apodora papuana':24.0871532021,'Liasis olivaceus':24.0871532021)Node4316558288:5.24037108992)Node4316584336:12.8697548268)Node4316583824:3.5026210108,(('Bothrochilus boa':25.4875902402,'Liasis albertisii':25.4875902402)Node4316585232:19.3585494438,'Morelia boeleni':44.846139684)Node4316585104:0.85376044557)Node4316583696:3.26334313874)Node4316559312:12.1691499629,('Python timoriensis':33.0773323832,'Python reticulatus':33.0773323832)Node4316585744:28.055060848)Node4316559184:12.7558559915,(('Python curtus':42.4795771694,('Python sebae':25.8300719638,'Python molurus':25.8300719638)Node4316606672:16.6495052056)Node4316606800:13.6364666177,'Python regius':56.1160437871)Node4316606544:17.7722054357)Node4316558736:40.1699176918)Node4316558800:26.0234610563)Node4316558480:71.3865451194,'Candola aspera':211.46817309)Node4316558352;
"""

def reference_tree_list_node_relationships():
    treelist_node_references = [
        {
            'Node4313761104' : NodeRelationship(parent_label=None, child_labels=['Node4313761232','Node4313801936'], edge_length=None, taxon_label=None),
            'Node4313761232' : NodeRelationship(parent_label='Node4313761104', child_labels=['Node4313761360','Node4313801552'], edge_length=78.6266408419, taxon_label=None),
            'Node4313761360' : NodeRelationship(parent_label='Node4313761232', child_labels=['Node4313761488','Node4313779728'], edge_length=123.936295332, taxon_label=None),
            'Node4313761488' : NodeRelationship(parent_label='Node4313761360', child_labels=['Node4313761616','Node4313779088'], edge_length=53.7887032791, taxon_label=None),
            'Node4313761616' : NodeRelationship(parent_label='Node4313761488', child_labels=['Node4313761744','Node4313778704'], edge_length=4.50960412336, taxon_label=None),
            'Node4313761744' : NodeRelationship(parent_label='Node4313761616', child_labels=['Node4313761872','Node4313777808'], edge_length=1.29166787247, taxon_label=None),
            'Node4313761872' : NodeRelationship(parent_label='Node4313761744', child_labels=['Node4313762000','Morelia boeleni'], edge_length=0.739270951321, taxon_label=None),
            'Node4313762000' : NodeRelationship(parent_label='Node4313761872', child_labels=['Node4313762128','Node4313762704'], edge_length=4.00932181352, taxon_label=None),
            'Node4313762128' : NodeRelationship(parent_label='Node4313762000', child_labels=['Antaresia childreni','Node4313762448'], edge_length=10.4795105126, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4313762128', child_labels=[], edge_length=3.82103266723, taxon_label='Antaresia childreni'),
            'Node4313762448' : NodeRelationship(parent_label='Node4313762128', child_labels=['Morelia clastolepis','Loxocemus bicolor'], edge_length=2.72293867732, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4313762448', child_labels=[], edge_length=1.09809398991, taxon_label='Morelia clastolepis'),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4313762448', child_labels=[], edge_length=1.09809398991, taxon_label='Loxocemus bicolor'),
            'Node4313762704' : NodeRelationship(parent_label='Node4313762000', child_labels=['Node4313762960','Node4313776720'], edge_length=3.67431549248, taxon_label=None),
            'Node4313762960' : NodeRelationship(parent_label='Node4313762704', child_labels=['Node4313763088','Node4313776208'], edge_length=4.47413907662, taxon_label=None),
            'Node4313763088' : NodeRelationship(parent_label='Node4313762960', child_labels=['Node4313763216','Node4313763600'], edge_length=1.10594701364, taxon_label=None),
            'Node4313763216' : NodeRelationship(parent_label='Node4313763088', child_labels=['Antaresia melanocephalus','Antaresia perthensis'], edge_length=3.56964311758, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4313763216', child_labels=[], edge_length=1.47649847948, taxon_label='Antaresia melanocephalus'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4313763216', child_labels=[], edge_length=1.47649847948, taxon_label='Antaresia perthensis'),
            'Node4313763600' : NodeRelationship(parent_label='Node4313763088', child_labels=['Python timoriensis','Liasis albertisii'], edge_length=3.54960462938, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4313763600', child_labels=[], edge_length=1.49653696768, taxon_label='Python timoriensis'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4313763600', child_labels=[], edge_length=1.49653696768, taxon_label='Liasis albertisii'),
            'Node4313776208' : NodeRelationship(parent_label='Node4313762960', child_labels=['Python sebae','Antaresia maculosa'], edge_length=0.946531890581, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4313776208', child_labels=[], edge_length=5.20555672012, taxon_label='Python sebae'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4313776208', child_labels=[], edge_length=5.20555672012, taxon_label='Antaresia maculosa'),
            'Node4313776720' : NodeRelationship(parent_label='Node4313762704', child_labels=['Node4313776912','Node4313777296'], edge_length=3.27824986702, taxon_label=None),
            'Node4313776912' : NodeRelationship(parent_label='Node4313776720', child_labels=['Python curtus','Bothrochilus boa'], edge_length=2.88197852526, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4313776912', child_labels=[], edge_length=4.46599929505, taxon_label='Python curtus'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4313776912', child_labels=[], edge_length=4.46599929505, taxon_label='Bothrochilus boa'),
            'Node4313777296' : NodeRelationship(parent_label='Node4313776720', child_labels=['Xenopeltis unicolor','Liasis fuscus'], edge_length=6.86415378064, taxon_label=None),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4313777296', child_labels=[], edge_length=0.483824039664, taxon_label='Xenopeltis unicolor'),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4313777296', child_labels=[], edge_length=0.483824039664, taxon_label='Liasis fuscus'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4313761872', child_labels=[], edge_length=18.3098649933, taxon_label='Morelia boeleni'),
            'Node4313777808' : NodeRelationship(parent_label='Node4313761744', child_labels=['Morelia tracyae','Node4313777680'], edge_length=16.8726126715, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4313777808', child_labels=[], edge_length=2.17652327312, taxon_label='Morelia tracyae'),
            'Node4313777680' : NodeRelationship(parent_label='Node4313777808', child_labels=['Node4313778192','Morelia bredli'], edge_length=1.67230791531, taxon_label=None),
            'Node4313778192' : NodeRelationship(parent_label='Node4313777680', child_labels=['Morelia viridisN','Python reticulatus'], edge_length=0.491713738136, taxon_label=None),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4313778192', child_labels=[], edge_length=0.0125016196671, taxon_label='Morelia viridisN'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4313778192', child_labels=[], edge_length=0.0125016196671, taxon_label='Python reticulatus'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4313777680', child_labels=[], edge_length=0.504215357803, taxon_label='Morelia bredli'),
            'Node4313778704' : NodeRelationship(parent_label='Node4313761616', child_labels=['Node4313778832','Morelia kinghorni'], edge_length=13.1562713928, taxon_label=None),
            'Node4313778832' : NodeRelationship(parent_label='Node4313778704', child_labels=['Liasis mackloti','Antaresia ramsayi'], edge_length=2.98661848623, taxon_label=None),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4313778832', child_labels=[], edge_length=4.19791393809, taxon_label='Liasis mackloti'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4313778832', child_labels=[], edge_length=4.19791393809, taxon_label='Antaresia ramsayi'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4313778704', child_labels=[], edge_length=7.18453242432, taxon_label='Morelia kinghorni'),
            'Node4313779088' : NodeRelationship(parent_label='Node4313761488', child_labels=['Node4313779472','Antaresia stimsoni'], edge_length=24.4333103788, taxon_label=None),
            'Node4313779472' : NodeRelationship(parent_label='Node4313779088', child_labels=['Python regius','Candola aspera'], edge_length=0.207889736001, taxon_label=None),
            'Python regius' : NodeRelationship(parent_label='Node4313779472', child_labels=[], edge_length=0.209207825685, taxon_label='Python regius'),
            'Candola aspera' : NodeRelationship(parent_label='Node4313779472', child_labels=[], edge_length=0.209207825685, taxon_label='Candola aspera'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4313779088', child_labels=[], edge_length=0.417097561686, taxon_label='Antaresia stimsoni'),
            'Node4313779728' : NodeRelationship(parent_label='Node4313761360', child_labels=['Node4313780112','Morelia nauta'], edge_length=69.6091072013, taxon_label=None),
            'Node4313780112' : NodeRelationship(parent_label='Node4313779728', child_labels=['Node4313800784','Node4313801168'], edge_length=1.24643505521, taxon_label=None),
            'Node4313800784' : NodeRelationship(parent_label='Node4313780112', child_labels=['Python molurus','Morelia amethistina'], edge_length=2.6364754365, taxon_label=None),
            'Python molurus' : NodeRelationship(parent_label='Node4313800784', child_labels=[], edge_length=5.1470935265, taxon_label='Python molurus'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4313800784', child_labels=[], edge_length=5.1470935265, taxon_label='Morelia amethistina'),
            'Node4313801168' : NodeRelationship(parent_label='Node4313780112', child_labels=['Morelia viridisS','Liasis olivaceus'], edge_length=7.12573328141, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4313801168', child_labels=[], edge_length=0.657835681585, taxon_label='Morelia viridisS'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4313801168', child_labels=[], edge_length=0.657835681585, taxon_label='Liasis olivaceus'),
            'Morelia nauta' : NodeRelationship(parent_label='Node4313779728', child_labels=[], edge_length=9.03000401821, taxon_label='Morelia nauta'),
            'Node4313801552' : NodeRelationship(parent_label='Node4313761232', child_labels=['Apodora papuana','Morelia carinata'], edge_length=202.422980755, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4313801552', child_labels=[], edge_length=0.152425796315, taxon_label='Apodora papuana'),
            'Morelia carinata' : NodeRelationship(parent_label='Node4313801552', child_labels=[], edge_length=0.152425796315, taxon_label='Morelia carinata'),
            'Node4313801936' : NodeRelationship(parent_label='Node4313761104', child_labels=['Morelia spilota','Morelia oenpelliensis'], edge_length=229.880308935, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4313801936', child_labels=[], edge_length=51.3217384582, taxon_label='Morelia spilota'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4313801936', child_labels=[], edge_length=51.3217384582, taxon_label='Morelia oenpelliensis'),
        },
        {
            'Node4313802576' : NodeRelationship(parent_label=None, child_labels=['Node4313802704','Loxocemus bicolor'], edge_length=None, taxon_label=None),
            'Node4313802704' : NodeRelationship(parent_label='Node4313802576', child_labels=['Candola aspera','Node4313803024'], edge_length=18.8917197007, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4313802704', child_labels=[], edge_length=81.2974222246, taxon_label='Candola aspera'),
            'Node4313803024' : NodeRelationship(parent_label='Node4313802704', child_labels=['Node4313803152','Node4316156752'], edge_length=33.1565984398, taxon_label=None),
            'Node4313803152' : NodeRelationship(parent_label='Node4313803024', child_labels=['Node4313803280','Node4313829328'], edge_length=6.57324583185, taxon_label=None),
            'Node4313803280' : NodeRelationship(parent_label='Node4313803152', child_labels=['Node4313803408','Node4313826384'], edge_length=0.76583222117, taxon_label=None),
            'Node4313803408' : NodeRelationship(parent_label='Node4313803280', child_labels=['Node4313803536','Node4313804432'], edge_length=4.5936111676, taxon_label=None),
            'Node4313803536' : NodeRelationship(parent_label='Node4313803408', child_labels=['Morelia amethistina','Node4313802960'], edge_length=15.1180336863, taxon_label=None),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4313803536', child_labels=[], edge_length=21.0901008779, taxon_label='Morelia amethistina'),
            'Node4313802960' : NodeRelationship(parent_label='Node4313803536', child_labels=['Morelia spilota','Node4313803792'], edge_length=3.38653541663, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4313802960', child_labels=[], edge_length=17.7035654613, taxon_label='Morelia spilota'),
            'Node4313803792' : NodeRelationship(parent_label='Node4313802960', child_labels=['Morelia oenpelliensis','Morelia bredli'], edge_length=2.6244717729, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4313803792', child_labels=[], edge_length=15.0790936884, taxon_label='Morelia oenpelliensis'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4313803792', child_labels=[], edge_length=15.0790936884, taxon_label='Morelia bredli'),
            'Node4313804432' : NodeRelationship(parent_label='Node4313803408', child_labels=['Node4313804560','Node4313825616'], edge_length=3.40951166184, taxon_label=None),
            'Node4313804560' : NodeRelationship(parent_label='Node4313804432', child_labels=['Morelia carinata','Node4313804240'], edge_length=5.00375936144, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4313804560', child_labels=[], edge_length=27.7948635409, taxon_label='Morelia carinata'),
            'Node4313804240' : NodeRelationship(parent_label='Node4313804560', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=8.32618746237, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4313804240', child_labels=[], edge_length=19.4686760786, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4313804240', child_labels=[], edge_length=19.4686760786, taxon_label='Morelia viridisN'),
            'Node4313825616' : NodeRelationship(parent_label='Node4313804432', child_labels=['Antaresia maculosa','Node4313825744'], edge_length=5.81119736053, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4313825616', child_labels=[], edge_length=26.9874255418, taxon_label='Antaresia maculosa'),
            'Node4313825744' : NodeRelationship(parent_label='Node4313825616', child_labels=['Node4313826128','Antaresia perthensis'], edge_length=2.25683638168, taxon_label=None),
            'Node4313826128' : NodeRelationship(parent_label='Node4313825744', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=16.6530983052, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4313826128', child_labels=[], edge_length=8.07749085501, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4313826128', child_labels=[], edge_length=8.07749085501, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4313825744', child_labels=[], edge_length=24.7305891602, taxon_label='Antaresia perthensis'),
            'Node4313826384' : NodeRelationship(parent_label='Node4313803280', child_labels=['Node4313826768','Node4313828048'], edge_length=1.27553605821, taxon_label=None),
            'Node4313826768' : NodeRelationship(parent_label='Node4313826384', child_labels=['Node4313826896','Node4313827664'], edge_length=4.33214615343, taxon_label=None),
            'Node4313826896' : NodeRelationship(parent_label='Node4313826768', child_labels=['Node4313827024','Node4313827408'], edge_length=3.37363071467, taxon_label=None),
            'Node4313827024' : NodeRelationship(parent_label='Node4313826896', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=26.1365403684, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4313827024', child_labels=[], edge_length=5.68389243709, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4313827024', child_labels=[], edge_length=5.68389243709, taxon_label='Liasis mackloti'),
            'Node4313827408' : NodeRelationship(parent_label='Node4313826896', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=12.1064068345, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4313827408', child_labels=[], edge_length=19.714025971, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4313827408', child_labels=[], edge_length=19.714025971, taxon_label='Liasis olivaceus'),
            'Node4313827664' : NodeRelationship(parent_label='Node4313826768', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=16.3762764593, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4313827664', child_labels=[], edge_length=18.8177870609, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4313827664', child_labels=[], edge_length=18.8177870609, taxon_label='Antaresia ramsayi'),
            'Node4313828048' : NodeRelationship(parent_label='Node4313826384', child_labels=['Morelia boeleni','Node4313828176'], edge_length=10.9652932592, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4313828048', child_labels=[], edge_length=28.5609164144, taxon_label='Morelia boeleni'),
            'Node4313828176' : NodeRelationship(parent_label='Node4313828048', child_labels=['Morelia tracyae','Node4313828432'], edge_length=11.3916491298, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4313828176', child_labels=[], edge_length=17.1692672846, taxon_label='Morelia tracyae'),
            'Node4313828432' : NodeRelationship(parent_label='Node4313828176', child_labels=['Node4313828816','Morelia clastolepis'], edge_length=10.4784522084, taxon_label=None),
            'Node4313828816' : NodeRelationship(parent_label='Node4313828432', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=1.44575855569, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4313828816', child_labels=[], edge_length=5.2450565205, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node4313828816', child_labels=[], edge_length=5.2450565205, taxon_label='Morelia nauta'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4313828432', child_labels=[], edge_length=6.69081507619, taxon_label='Morelia clastolepis'),
            'Node4313829328' : NodeRelationship(parent_label='Node4313803152', child_labels=['Node4316155984','Node4313802512'], edge_length=2.53377266123, taxon_label=None),
            'Node4316155984' : NodeRelationship(parent_label='Node4313829328', child_labels=['Python timoriensis','Python reticulatus'], edge_length=23.8923728386, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316155984', child_labels=[], edge_length=15.1414324532, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316155984', child_labels=[], edge_length=15.1414324532, taxon_label='Python reticulatus'),
            'Node4313802512' : NodeRelationship(parent_label='Node4313829328', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=22.5856079922, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4313802512', child_labels=[], edge_length=16.4481972995, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4313802512', child_labels=[], edge_length=16.4481972995, taxon_label='Liasis albertisii'),
            'Node4316156752' : NodeRelationship(parent_label='Node4313803024', child_labels=['Xenopeltis unicolor','Node4316156624'], edge_length=16.2594519516, taxon_label=None),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316156752', child_labels=[], edge_length=31.8813718333, taxon_label='Xenopeltis unicolor'),
            'Node4316156624' : NodeRelationship(parent_label='Node4316156752', child_labels=['Node4316157136','Python regius'], edge_length=5.97313984611, taxon_label=None),
            'Node4316157136' : NodeRelationship(parent_label='Node4316156624', child_labels=['Python curtus','Node4316157008'], edge_length=8.94343133576, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316157136', child_labels=[], edge_length=16.9648006514, taxon_label='Python curtus'),
            'Node4316157008' : NodeRelationship(parent_label='Node4316157136', child_labels=['Python sebae','Python molurus'], edge_length=4.66373979181, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316157008', child_labels=[], edge_length=12.3010608596, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316157008', child_labels=[], edge_length=12.3010608596, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316156624', child_labels=[], edge_length=25.9082319872, taxon_label='Python regius'),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4313802576', child_labels=[], edge_length=100.189141925, taxon_label='Loxocemus bicolor'),
        },
        {
            'Node4316158288' : NodeRelationship(parent_label=None, child_labels=['Candola aspera','Node4316158608'], edge_length=None, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4316158288', child_labels=[], edge_length=109.372663833, taxon_label='Candola aspera'),
            'Node4316158608' : NodeRelationship(parent_label='Node4316158288', child_labels=['Node4316158736','Node4316203024'], edge_length=23.0231215792, taxon_label=None),
            'Node4316158736' : NodeRelationship(parent_label='Node4316158608', child_labels=['Node4316158864','Node4316202128'], edge_length=25.5450384687, taxon_label=None),
            'Node4316158864' : NodeRelationship(parent_label='Node4316158736', child_labels=['Node4316158992','Node4316201552'], edge_length=8.2214192007, taxon_label=None),
            'Node4316158992' : NodeRelationship(parent_label='Node4316158864', child_labels=['Node4316159120','Node4316179792'], edge_length=6.47536868825, taxon_label=None),
            'Node4316159120' : NodeRelationship(parent_label='Node4316158992', child_labels=['Node4316159248','Node4316178128'], edge_length=2.53965409687, taxon_label=None),
            'Node4316159248' : NodeRelationship(parent_label='Node4316159120', child_labels=['Node4316159376','Node4316177104'], edge_length=2.41802497137, taxon_label=None),
            'Node4316159376' : NodeRelationship(parent_label='Node4316159248', child_labels=['Node4316159504','Node4316176720'], edge_length=6.39175712039, taxon_label=None),
            'Node4316159504' : NodeRelationship(parent_label='Node4316159376', child_labels=['Node4316159632','Node4316176464'], edge_length=4.82772953939, taxon_label=None),
            'Node4316159632' : NodeRelationship(parent_label='Node4316159504', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=18.9362130146, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316159632', child_labels=[], edge_length=10.9943371535, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316159632', child_labels=[], edge_length=10.9943371535, taxon_label='Liasis mackloti'),
            'Node4316176464' : NodeRelationship(parent_label='Node4316159504', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=8.55401069191, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316176464', child_labels=[], edge_length=21.3765394762, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316176464', child_labels=[], edge_length=21.3765394762, taxon_label='Liasis olivaceus'),
            'Node4316176720' : NodeRelationship(parent_label='Node4316159376', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=17.7955396972, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316176720', child_labels=[], edge_length=16.9627400103, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316176720', child_labels=[], edge_length=16.9627400103, taxon_label='Antaresia ramsayi'),
            'Node4316177104' : NodeRelationship(parent_label='Node4316159248', child_labels=['Morelia boeleni','Node4316177232'], edge_length=9.11329596086, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316177104', child_labels=[], edge_length=32.036740867, taxon_label='Morelia boeleni'),
            'Node4316177232' : NodeRelationship(parent_label='Node4316177104', child_labels=['Morelia tracyae','Node4316177488'], edge_length=14.7791918926, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316177232', child_labels=[], edge_length=17.2575489744, taxon_label='Morelia tracyae'),
            'Node4316177488' : NodeRelationship(parent_label='Node4316177232', child_labels=['Node4316177872','Morelia clastolepis'], edge_length=13.3651095585, taxon_label=None),
            'Node4316177872' : NodeRelationship(parent_label='Node4316177488', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=0.439451186875, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316177872', child_labels=[], edge_length=3.4529882291, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316177872', child_labels=[], edge_length=3.4529882291, taxon_label='Morelia nauta'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316177488', child_labels=[], edge_length=3.89243941597, taxon_label='Morelia clastolepis'),
            'Node4316178128' : NodeRelationship(parent_label='Node4316159120', child_labels=['Node4316178512','Node4316179280'], edge_length=5.90175129715, taxon_label=None),
            'Node4316178512' : NodeRelationship(parent_label='Node4316178128', child_labels=['Node4316178640','Node4316179024'], edge_length=11.5242004723, taxon_label=None),
            'Node4316178640' : NodeRelationship(parent_label='Node4316178512', child_labels=['Morelia oenpelliensis','Morelia amethistina'], edge_length=1.13984981318, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316178640', child_labels=[], edge_length=25.0022602166, taxon_label='Morelia oenpelliensis'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316178640', child_labels=[], edge_length=25.0022602166, taxon_label='Morelia amethistina'),
            'Node4316179024' : NodeRelationship(parent_label='Node4316178512', child_labels=['Morelia spilota','Morelia bredli'], edge_length=13.1543046788, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316179024', child_labels=[], edge_length=12.987805351, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316179024', child_labels=[], edge_length=12.987805351, taxon_label='Morelia bredli'),
            'Node4316179280' : NodeRelationship(parent_label='Node4316178128', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=13.3915851761, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316179280', child_labels=[], edge_length=24.274725326, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316179280', child_labels=[], edge_length=24.274725326, taxon_label='Liasis albertisii'),
            'Node4316179792' : NodeRelationship(parent_label='Node4316158992', child_labels=['Node4316179920','Node4316180432'], edge_length=7.82653683343, taxon_label=None),
            'Node4316179920' : NodeRelationship(parent_label='Node4316179792', child_labels=['Morelia carinata','Node4316179664'], edge_length=7.8146690322, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316179920', child_labels=[], edge_length=30.4665100305, taxon_label='Morelia carinata'),
            'Node4316179664' : NodeRelationship(parent_label='Node4316179920', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=11.0043198537, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316179664', child_labels=[], edge_length=19.4621901768, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316179664', child_labels=[], edge_length=19.4621901768, taxon_label='Morelia viridisN'),
            'Node4316180432' : NodeRelationship(parent_label='Node4316179792', child_labels=['Antaresia maculosa','Node4316158224'], edge_length=5.10842077756, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316180432', child_labels=[], edge_length=33.1727582851, taxon_label='Antaresia maculosa'),
            'Node4316158224' : NodeRelationship(parent_label='Node4316180432', child_labels=['Node4316201296','Antaresia perthensis'], edge_length=4.7141022378, taxon_label=None),
            'Node4316201296' : NodeRelationship(parent_label='Node4316158224', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=19.4308450954, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316201296', child_labels=[], edge_length=9.02781095195, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316201296', child_labels=[], edge_length=9.02781095195, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316158224', child_labels=[], edge_length=28.4586560473, taxon_label='Antaresia perthensis'),
            'Node4316201552' : NodeRelationship(parent_label='Node4316158864', child_labels=['Python timoriensis','Python reticulatus'], edge_length=22.2791730332, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316201552', child_labels=[], edge_length=30.3039115511, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316201552', child_labels=[], edge_length=30.3039115511, taxon_label='Python reticulatus'),
            'Node4316202128' : NodeRelationship(parent_label='Node4316158736', child_labels=['Node4316202384','Python regius'], edge_length=17.448321597, taxon_label=None),
            'Node4316202384' : NodeRelationship(parent_label='Node4316202128', child_labels=['Python curtus','Node4316202256'], edge_length=6.14068810478, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316202384', child_labels=[], edge_length=37.2154940833, taxon_label='Python curtus'),
            'Node4316202256' : NodeRelationship(parent_label='Node4316202384', child_labels=['Python sebae','Python molurus'], edge_length=15.2807163378, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316202256', child_labels=[], edge_length=21.9347777454, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316202256', child_labels=[], edge_length=21.9347777454, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316202128', child_labels=[], edge_length=43.356182188, taxon_label='Python regius'),
            'Node4316203024' : NodeRelationship(parent_label='Node4316158608', child_labels=['Loxocemus bicolor','Xenopeltis unicolor'], edge_length=3.40243621372, taxon_label=None),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316203024', child_labels=[], edge_length=82.94710604, taxon_label='Loxocemus bicolor'),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316203024', child_labels=[], edge_length=82.94710604, taxon_label='Xenopeltis unicolor'),
        },
        {
            'Node4316203664' : NodeRelationship(parent_label=None, child_labels=['Candola aspera','Node4316203984'], edge_length=None, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4316203664', child_labels=[], edge_length=157.750076773, taxon_label='Candola aspera'),
            'Node4316203984' : NodeRelationship(parent_label='Node4316203664', child_labels=['Node4316204112','Xenopeltis unicolor'], edge_length=44.9789688242, taxon_label=None),
            'Node4316204112' : NodeRelationship(parent_label='Node4316203984', child_labels=['Loxocemus bicolor','Node4316203920'], edge_length=19.5811677101, taxon_label=None),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316204112', child_labels=[], edge_length=93.189940239, taxon_label='Loxocemus bicolor'),
            'Node4316203920' : NodeRelationship(parent_label='Node4316204112', child_labels=['Node4316204496','Node4316251856'], edge_length=27.291533515, taxon_label=None),
            'Node4316204496' : NodeRelationship(parent_label='Node4316203920', child_labels=['Node4316204624','Node4316251600'], edge_length=10.3875398007, taxon_label=None),
            'Node4316204624' : NodeRelationship(parent_label='Node4316204496', child_labels=['Node4316204752','Node4316227472'], edge_length=9.74044751021, taxon_label=None),
            'Node4316204752' : NodeRelationship(parent_label='Node4316204624', child_labels=['Node4316204880','Morelia boeleni'], edge_length=2.88538115771, taxon_label=None),
            'Node4316204880' : NodeRelationship(parent_label='Node4316204752', child_labels=['Node4316205008','Node4316226768'], edge_length=1.34047834167, taxon_label=None),
            'Node4316205008' : NodeRelationship(parent_label='Node4316204880', child_labels=['Node4316225680','Node4316226448'], edge_length=2.31767871982, taxon_label=None),
            'Node4316225680' : NodeRelationship(parent_label='Node4316205008', child_labels=['Node4316225808','Node4316226192'], edge_length=6.3930928479, taxon_label=None),
            'Node4316225808' : NodeRelationship(parent_label='Node4316225680', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=24.6792964706, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316225808', child_labels=[], edge_length=8.15449187544, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316225808', child_labels=[], edge_length=8.15449187544, taxon_label='Liasis mackloti'),
            'Node4316226192' : NodeRelationship(parent_label='Node4316225680', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=8.52936801714, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316226192', child_labels=[], edge_length=24.3044203289, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316226192', child_labels=[], edge_length=24.3044203289, taxon_label='Liasis olivaceus'),
            'Node4316226448' : NodeRelationship(parent_label='Node4316205008', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=23.2404397828, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316226448', child_labels=[], edge_length=15.9864414111, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316226448', child_labels=[], edge_length=15.9864414111, taxon_label='Antaresia ramsayi'),
            'Node4316226768' : NodeRelationship(parent_label='Node4316204880', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=18.2547475502, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316226768', child_labels=[], edge_length=23.2898123636, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316226768', child_labels=[], edge_length=23.2898123636, taxon_label='Liasis albertisii'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316204752', child_labels=[], edge_length=42.8850382554, taxon_label='Morelia boeleni'),
            'Node4316227472' : NodeRelationship(parent_label='Node4316204624', child_labels=['Node4316227600','Node4316229392'], edge_length=2.72676824396, taxon_label=None),
            'Node4316227600' : NodeRelationship(parent_label='Node4316227472', child_labels=['Node4316227728','Node4316228112'], edge_length=14.2175774566, taxon_label=None),
            'Node4316227728' : NodeRelationship(parent_label='Node4316227600', child_labels=['Morelia spilota','Morelia bredli'], edge_length=12.5474231006, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316227728', child_labels=[], edge_length=16.278650612, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316227728', child_labels=[], edge_length=16.278650612, taxon_label='Morelia bredli'),
            'Node4316228112' : NodeRelationship(parent_label='Node4316227600', child_labels=['Morelia oenpelliensis','Node4316227984'], edge_length=1.26678175478, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316228112', child_labels=[], edge_length=27.5592919578, taxon_label='Morelia oenpelliensis'),
            'Node4316227984' : NodeRelationship(parent_label='Node4316228112', child_labels=['Morelia tracyae','Node4316228368'], edge_length=13.3039152583, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316227984', child_labels=[], edge_length=14.2553766995, taxon_label='Morelia tracyae'),
            'Node4316228368' : NodeRelationship(parent_label='Node4316227984', child_labels=['Node4316228752','Morelia amethistina'], edge_length=1.66170791236, taxon_label=None),
            'Node4316228752' : NodeRelationship(parent_label='Node4316228368', child_labels=['Node4316228880','Morelia nauta'], edge_length=8.89489387836, taxon_label=None),
            'Node4316228880' : NodeRelationship(parent_label='Node4316228752', child_labels=['Morelia clastolepis','Morelia kinghorni'], edge_length=0.230110019205, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316228880', child_labels=[], edge_length=3.46866488962, taxon_label='Morelia clastolepis'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316228880', child_labels=[], edge_length=3.46866488962, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316228752', child_labels=[], edge_length=3.69877490882, taxon_label='Morelia nauta'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316228368', child_labels=[], edge_length=12.5936687872, taxon_label='Morelia amethistina'),
            'Node4316229392' : NodeRelationship(parent_label='Node4316227472', child_labels=['Node4316229520','Node4316250576'], edge_length=3.9347409374, taxon_label=None),
            'Node4316229520' : NodeRelationship(parent_label='Node4316229392', child_labels=['Morelia carinata','Node4316203600'], edge_length=5.88218975316, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316229520', child_labels=[], edge_length=33.2267204786, taxon_label='Morelia carinata'),
            'Node4316203600' : NodeRelationship(parent_label='Node4316229520', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=13.4306199458, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316203600', child_labels=[], edge_length=19.7961005329, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316203600', child_labels=[], edge_length=19.7961005329, taxon_label='Morelia viridisN'),
            'Node4316250576' : NodeRelationship(parent_label='Node4316229392', child_labels=['Antaresia maculosa','Node4316250704'], edge_length=7.11128547149, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316250576', child_labels=[], edge_length=31.9976247603, taxon_label='Antaresia maculosa'),
            'Node4316250704' : NodeRelationship(parent_label='Node4316250576', child_labels=['Node4316251088','Antaresia perthensis'], edge_length=4.71436528425, taxon_label=None),
            'Node4316251088' : NodeRelationship(parent_label='Node4316250704', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=15.9285543528, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316251088', child_labels=[], edge_length=11.3547051232, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316251088', child_labels=[], edge_length=11.3547051232, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316250704', child_labels=[], edge_length=27.2832594761, taxon_label='Antaresia perthensis'),
            'Node4316251600' : NodeRelationship(parent_label='Node4316204496', child_labels=['Python timoriensis','Python reticulatus'], edge_length=22.1202819118, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316251600', child_labels=[], edge_length=33.3905850116, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316251600', child_labels=[], edge_length=33.3905850116, taxon_label='Python reticulatus'),
            'Node4316251856' : NodeRelationship(parent_label='Node4316203920', child_labels=['Node4316252112','Python regius'], edge_length=19.8895314814, taxon_label=None),
            'Node4316252112' : NodeRelationship(parent_label='Node4316251856', child_labels=['Python curtus','Node4316251984'], edge_length=11.1907198563, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316252112', child_labels=[], edge_length=34.8181553863, taxon_label='Python curtus'),
            'Node4316251984' : NodeRelationship(parent_label='Node4316252112', child_labels=['Python sebae','Python molurus'], edge_length=7.89583224277, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316251984', child_labels=[], edge_length=26.9223231435, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316251984', child_labels=[], edge_length=26.9223231435, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316251856', child_labels=[], edge_length=46.0088752426, taxon_label='Python regius'),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316203984', child_labels=[], edge_length=112.771107949, taxon_label='Xenopeltis unicolor'),
        },
        {
            'Node4316253136' : NodeRelationship(parent_label=None, child_labels=['Candola aspera','Node4316253456'], edge_length=None, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4316253136', child_labels=[], edge_length=117.405482731, taxon_label='Candola aspera'),
            'Node4316253456' : NodeRelationship(parent_label='Node4316253136', child_labels=['Xenopeltis unicolor','Node4316253392'], edge_length=26.8517813278, taxon_label=None),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316253456', child_labels=[], edge_length=90.5537014032, taxon_label='Xenopeltis unicolor'),
            'Node4316253392' : NodeRelationship(parent_label='Node4316253456', child_labels=['Node4316253776','Loxocemus bicolor'], edge_length=16.5508683851, taxon_label=None),
            'Node4316253776' : NodeRelationship(parent_label='Node4316253392', child_labels=['Node4316253904','Node4316301328'], edge_length=16.6353899137, taxon_label=None),
            'Node4316253904' : NodeRelationship(parent_label='Node4316253776', child_labels=['Node4316254032','Node4316301072'], edge_length=6.56321044476, taxon_label=None),
            'Node4316254032' : NodeRelationship(parent_label='Node4316253904', child_labels=['Node4316254160','Node4316276944'], edge_length=12.3993667148, taxon_label=None),
            'Node4316254160' : NodeRelationship(parent_label='Node4316254032', child_labels=['Morelia boeleni','Node4316275024'], edge_length=1.85501747057, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316254160', child_labels=[], edge_length=36.5498484742, taxon_label='Morelia boeleni'),
            'Node4316275024' : NodeRelationship(parent_label='Node4316254160', child_labels=['Node4316275152','Node4316276432'], edge_length=1.04661210215, taxon_label=None),
            'Node4316275152' : NodeRelationship(parent_label='Node4316275024', child_labels=['Node4316275280','Node4316276048'], edge_length=6.27408700912, taxon_label=None),
            'Node4316275280' : NodeRelationship(parent_label='Node4316275152', child_labels=['Node4316275408','Node4316275792'], edge_length=4.14502379891, taxon_label=None),
            'Node4316275408' : NodeRelationship(parent_label='Node4316275280', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=18.7911640437, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316275408', child_labels=[], edge_length=6.29296152027, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316275408', child_labels=[], edge_length=6.29296152027, taxon_label='Liasis mackloti'),
            'Node4316275792' : NodeRelationship(parent_label='Node4316275280', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=4.65165886356, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316275792', child_labels=[], edge_length=20.4324667004, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316275792', child_labels=[], edge_length=20.4324667004, taxon_label='Liasis olivaceus'),
            'Node4316276048' : NodeRelationship(parent_label='Node4316275152', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=17.4625828861, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316276048', child_labels=[], edge_length=11.7665664768, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316276048', child_labels=[], edge_length=11.7665664768, taxon_label='Antaresia ramsayi'),
            'Node4316276432' : NodeRelationship(parent_label='Node4316275024', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=15.4075337774, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316276432', child_labels=[], edge_length=20.0957025946, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316276432', child_labels=[], edge_length=20.0957025946, taxon_label='Liasis albertisii'),
            'Node4316276944' : NodeRelationship(parent_label='Node4316254032', child_labels=['Node4316277072','Node4316299408'], edge_length=2.5597754437, taxon_label=None),
            'Node4316277072' : NodeRelationship(parent_label='Node4316276944', child_labels=['Morelia oenpelliensis','Node4316276816'], edge_length=10.9999278199, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316277072', child_labels=[], edge_length=24.8451626811, taxon_label='Morelia oenpelliensis'),
            'Node4316276816' : NodeRelationship(parent_label='Node4316277072', child_labels=['Node4316277392','Node4316278352'], edge_length=2.78884378851, taxon_label=None),
            'Node4316277392' : NodeRelationship(parent_label='Node4316276816', child_labels=['Morelia tracyae','Node4316277584'], edge_length=7.72979171285, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316277392', child_labels=[], edge_length=14.3265271797, taxon_label='Morelia tracyae'),
            'Node4316277584' : NodeRelationship(parent_label='Node4316277392', child_labels=['Node4316277712','Morelia amethistina'], edge_length=3.38770784397, taxon_label=None),
            'Node4316277712' : NodeRelationship(parent_label='Node4316277584', child_labels=['Morelia clastolepis','Node4316253008'], edge_length=5.65307392619, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316277712', child_labels=[], edge_length=5.28574540958, taxon_label='Morelia clastolepis'),
            'Node4316253008' : NodeRelationship(parent_label='Node4316277712', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=2.55413552294, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316253008', child_labels=[], edge_length=2.73160988664, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316253008', child_labels=[], edge_length=2.73160988664, taxon_label='Morelia nauta'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316277584', child_labels=[], edge_length=10.9388193358, taxon_label='Morelia amethistina'),
            'Node4316278352' : NodeRelationship(parent_label='Node4316276816', child_labels=['Morelia spilota','Morelia bredli'], edge_length=12.4767981898, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316278352', child_labels=[], edge_length=9.57952070277, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316278352', child_labels=[], edge_length=9.57952070277, taxon_label='Morelia bredli'),
            'Node4316299408' : NodeRelationship(parent_label='Node4316276944', child_labels=['Node4316299536','Node4316300048'], edge_length=4.77611939702, taxon_label=None),
            'Node4316299536' : NodeRelationship(parent_label='Node4316299408', child_labels=['Morelia carinata','Node4316299344'], edge_length=6.62946104565, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316299536', child_labels=[], edge_length=24.4395100584, taxon_label='Morelia carinata'),
            'Node4316299344' : NodeRelationship(parent_label='Node4316299536', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=7.74158436759, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316299344', child_labels=[], edge_length=16.6979256908, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316299344', child_labels=[], edge_length=16.6979256908, taxon_label='Morelia viridisN'),
            'Node4316300048' : NodeRelationship(parent_label='Node4316299408', child_labels=['Antaresia maculosa','Node4316300176'], edge_length=3.07031045323, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316300048', child_labels=[], edge_length=27.9986606508, taxon_label='Antaresia maculosa'),
            'Node4316300176' : NodeRelationship(parent_label='Node4316300048', child_labels=['Node4316300560','Antaresia perthensis'], edge_length=4.88214307652, taxon_label=None),
            'Node4316300560' : NodeRelationship(parent_label='Node4316300176', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=16.2064963991, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316300560', child_labels=[], edge_length=6.91002117518, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316300560', child_labels=[], edge_length=6.91002117518, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316300176', child_labels=[], edge_length=23.1165175743, taxon_label='Antaresia perthensis'),
            'Node4316301072' : NodeRelationship(parent_label='Node4316253904', child_labels=['Python timoriensis','Python reticulatus'], edge_length=23.3069353709, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316301072', child_labels=[], edge_length=27.4972972887, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316301072', child_labels=[], edge_length=27.4972972887, taxon_label='Python reticulatus'),
            'Node4316301328' : NodeRelationship(parent_label='Node4316253776', child_labels=['Node4316301584','Python regius'], edge_length=18.1852453647, taxon_label=None),
            'Node4316301584' : NodeRelationship(parent_label='Node4316301328', child_labels=['Python curtus','Node4316301456'], edge_length=6.08258726043, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316301584', child_labels=[], edge_length=33.0996104792, taxon_label='Python curtus'),
            'Node4316301456' : NodeRelationship(parent_label='Node4316301584', child_labels=['Python sebae','Python molurus'], edge_length=14.293345062, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316301456', child_labels=[], edge_length=18.8062654172, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316301456', child_labels=[], edge_length=18.8062654172, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316301328', child_labels=[], edge_length=39.1821977396, taxon_label='Python regius'),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316253392', child_labels=[], edge_length=74.0028330181, taxon_label='Loxocemus bicolor'),
        },
        {
            'Node4316302608' : NodeRelationship(parent_label=None, child_labels=['Candola aspera','Node4316302928'], edge_length=None, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4316302608', child_labels=[], edge_length=126.159749307, taxon_label='Candola aspera'),
            'Node4316302928' : NodeRelationship(parent_label='Node4316302608', child_labels=['Node4316303056','Xenopeltis unicolor'], edge_length=29.5182899091, taxon_label=None),
            'Node4316303056' : NodeRelationship(parent_label='Node4316302928', child_labels=['Loxocemus bicolor','Node4316302864'], edge_length=20.4272017262, taxon_label=None),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316303056', child_labels=[], edge_length=76.2142576718, taxon_label='Loxocemus bicolor'),
            'Node4316302864' : NodeRelationship(parent_label='Node4316303056', child_labels=['Node4316323920','Node4316346704'], edge_length=22.3340694246, taxon_label=None),
            'Node4316323920' : NodeRelationship(parent_label='Node4316302864', child_labels=['Node4316324048','Node4316346448'], edge_length=7.33856679933, taxon_label=None),
            'Node4316324048' : NodeRelationship(parent_label='Node4316323920', child_labels=['Node4316324176','Node4316326288'], edge_length=8.24675869634, taxon_label=None),
            'Node4316324176' : NodeRelationship(parent_label='Node4316324048', child_labels=['Morelia boeleni','Node4316324496'], edge_length=1.78294733026, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316324176', child_labels=[], edge_length=36.5119154213, taxon_label='Morelia boeleni'),
            'Node4316324496' : NodeRelationship(parent_label='Node4316324176', child_labels=['Node4316324624','Node4316326032'], edge_length=1.38475127065, taxon_label=None),
            'Node4316324624' : NodeRelationship(parent_label='Node4316324496', child_labels=['Node4316324752','Node4316325136'], edge_length=2.24740276648, taxon_label=None),
            'Node4316324752' : NodeRelationship(parent_label='Node4316324624', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=17.772328432, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316324752', child_labels=[], edge_length=15.1074329521, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316324752', child_labels=[], edge_length=15.1074329521, taxon_label='Antaresia ramsayi'),
            'Node4316325136' : NodeRelationship(parent_label='Node4316324624', child_labels=['Node4316325264','Node4316325648'], edge_length=4.4885389798, taxon_label=None),
            'Node4316325264' : NodeRelationship(parent_label='Node4316325136', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=20.8200876951, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316325264', child_labels=[], edge_length=7.5711347092, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316325264', child_labels=[], edge_length=7.5711347092, taxon_label='Liasis mackloti'),
            'Node4316325648' : NodeRelationship(parent_label='Node4316325136', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=4.37289319177, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316325648', child_labels=[], edge_length=24.0183292126, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316325648', child_labels=[], edge_length=24.0183292126, taxon_label='Liasis olivaceus'),
            'Node4316326032' : NodeRelationship(parent_label='Node4316324496', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=12.9799500845, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316326032', child_labels=[], edge_length=22.1472140662, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316326032', child_labels=[], edge_length=22.1472140662, taxon_label='Liasis albertisii'),
            'Node4316326288' : NodeRelationship(parent_label='Node4316324048', child_labels=['Node4316326416','Node4316344784'], edge_length=3.26747919544, taxon_label=None),
            'Node4316326416' : NodeRelationship(parent_label='Node4316326288', child_labels=['Node4316326544','Node4316327184'], edge_length=12.2960670728, taxon_label=None),
            'Node4316326544' : NodeRelationship(parent_label='Node4316326416', child_labels=['Node4316326672','Morelia oenpelliensis'], edge_length=1.62908011033, taxon_label=None),
            'Node4316326672' : NodeRelationship(parent_label='Node4316326544', child_labels=['Morelia spilota','Morelia bredli'], edge_length=11.7610480778, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316326672', child_labels=[], edge_length=9.34118829512, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316326672', child_labels=[], edge_length=9.34118829512, taxon_label='Morelia bredli'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316326544', child_labels=[], edge_length=21.1022363729, taxon_label='Morelia oenpelliensis'),
            'Node4316327184' : NodeRelationship(parent_label='Node4316326416', child_labels=['Morelia tracyae','Node4316326928'], edge_length=6.20414166284, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316327184', child_labels=[], edge_length=16.5271748204, taxon_label='Morelia tracyae'),
            'Node4316326928' : NodeRelationship(parent_label='Node4316327184', child_labels=['Node4316327568','Morelia amethistina'], edge_length=4.27339647744, taxon_label=None),
            'Node4316327568' : NodeRelationship(parent_label='Node4316326928', child_labels=['Node4316327696','Morelia clastolepis'], edge_length=7.91289243511, taxon_label=None),
            'Node4316327696' : NodeRelationship(parent_label='Node4316327568', child_labels=['Morelia nauta','Morelia kinghorni'], edge_length=1.89524872622, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316327696', child_labels=[], edge_length=2.44563718163, taxon_label='Morelia nauta'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316327696', child_labels=[], edge_length=2.44563718163, taxon_label='Morelia kinghorni'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316327568', child_labels=[], edge_length=4.34088590785, taxon_label='Morelia clastolepis'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316326928', child_labels=[], edge_length=12.253778343, taxon_label='Morelia amethistina'),
            'Node4316344784' : NodeRelationship(parent_label='Node4316326288', child_labels=['Node4316344912','Node4316345424'], edge_length=4.17834235089, taxon_label=None),
            'Node4316344912' : NodeRelationship(parent_label='Node4316344784', child_labels=['Morelia carinata','Node4316344400'], edge_length=5.03710806168, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316344912', child_labels=[], edge_length=25.8119331435, taxon_label='Morelia carinata'),
            'Node4316344400' : NodeRelationship(parent_label='Node4316344912', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=11.4337686931, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316344400', child_labels=[], edge_length=14.3781644505, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316344400', child_labels=[], edge_length=14.3781644505, taxon_label='Morelia viridisN'),
            'Node4316345424' : NodeRelationship(parent_label='Node4316344784', child_labels=['Antaresia maculosa','Node4316345552'], edge_length=4.09757269601, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316345424', child_labels=[], edge_length=26.7514685092, taxon_label='Antaresia maculosa'),
            'Node4316345552' : NodeRelationship(parent_label='Node4316345424', child_labels=['Node4316345936','Antaresia perthensis'], edge_length=2.33638652753, taxon_label=None),
            'Node4316345936' : NodeRelationship(parent_label='Node4316345552', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=17.4436719435, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316345936', child_labels=[], edge_length=6.97141003814, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316345936', child_labels=[], edge_length=6.97141003814, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316345552', child_labels=[], edge_length=24.4150819817, taxon_label='Antaresia perthensis'),
            'Node4316346448' : NodeRelationship(parent_label='Node4316323920', child_labels=['Python timoriensis','Python reticulatus'], edge_length=23.556392099, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316346448', child_labels=[], edge_length=22.9852293488, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316346448', child_labels=[], edge_length=22.9852293488, taxon_label='Python reticulatus'),
            'Node4316346704' : NodeRelationship(parent_label='Node4316302864', child_labels=['Node4316346960','Python regius'], edge_length=17.9081848887, taxon_label=None),
            'Node4316346960' : NodeRelationship(parent_label='Node4316346704', child_labels=['Python curtus','Node4316346832'], edge_length=6.75664224871, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316346960', child_labels=[], edge_length=29.2153611098, taxon_label='Python curtus'),
            'Node4316346832' : NodeRelationship(parent_label='Node4316346960', child_labels=['Python sebae','Python molurus'], edge_length=8.15978945225, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316346832', child_labels=[], edge_length=21.0555716576, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316346832', child_labels=[], edge_length=21.0555716576, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316346704', child_labels=[], edge_length=35.9720033585, taxon_label='Python regius'),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316302928', child_labels=[], edge_length=96.6414593981, taxon_label='Xenopeltis unicolor'),
        },
        {
            'Node4316347984' : NodeRelationship(parent_label=None, child_labels=['Candola aspera','Node4316348304'], edge_length=None, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4316347984', child_labels=[], edge_length=124.564186516, taxon_label='Candola aspera'),
            'Node4316348304' : NodeRelationship(parent_label='Node4316347984', child_labels=['Node4316368976','Xenopeltis unicolor'], edge_length=36.3676780441, taxon_label=None),
            'Node4316368976' : NodeRelationship(parent_label='Node4316348304', child_labels=['Node4316369104','Loxocemus bicolor'], edge_length=11.1789504571, taxon_label=None),
            'Node4316369104' : NodeRelationship(parent_label='Node4316368976', child_labels=['Node4316369232','Node4316396048'], edge_length=20.3346663059, taxon_label=None),
            'Node4316369232' : NodeRelationship(parent_label='Node4316369104', child_labels=['Node4316369360','Node4316395664'], edge_length=10.7040090023, taxon_label=None),
            'Node4316369360' : NodeRelationship(parent_label='Node4316369232', child_labels=['Node4316369488','Node4316372752'], edge_length=7.28944403695, taxon_label=None),
            'Node4316369488' : NodeRelationship(parent_label='Node4316369360', child_labels=['Node4316369616','Node4316371280'], edge_length=2.76745367097, taxon_label=None),
            'Node4316369616' : NodeRelationship(parent_label='Node4316369488', child_labels=['Node4316369744','Node4316370128'], edge_length=12.5842551962, taxon_label=None),
            'Node4316369744' : NodeRelationship(parent_label='Node4316369616', child_labels=['Morelia spilota','Morelia bredli'], edge_length=12.0191169658, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316369744', child_labels=[], edge_length=11.3186128368, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316369744', child_labels=[], edge_length=11.3186128368, taxon_label='Morelia bredli'),
            'Node4316370128' : NodeRelationship(parent_label='Node4316369616', child_labels=['Morelia oenpelliensis','Node4316370000'], edge_length=2.25737379321, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316370128', child_labels=[], edge_length=21.0803560094, taxon_label='Morelia oenpelliensis'),
            'Node4316370000' : NodeRelationship(parent_label='Node4316370128', child_labels=['Morelia tracyae','Node4316370384'], edge_length=6.93019779087, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316370000', child_labels=[], edge_length=14.1501582185, taxon_label='Morelia tracyae'),
            'Node4316370384' : NodeRelationship(parent_label='Node4316370000', child_labels=['Node4316347856','Morelia amethistina'], edge_length=3.68709033509, taxon_label=None),
            'Node4316347856' : NodeRelationship(parent_label='Node4316370384', child_labels=['Morelia nauta','Node4316370640'], edge_length=4.34576715349, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316347856', child_labels=[], edge_length=6.11730072991, taxon_label='Morelia nauta'),
            'Node4316370640' : NodeRelationship(parent_label='Node4316347856', child_labels=['Morelia kinghorni','Morelia clastolepis'], edge_length=2.07374196069, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316370640', child_labels=[], edge_length=4.04355876922, taxon_label='Morelia kinghorni'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316370640', child_labels=[], edge_length=4.04355876922, taxon_label='Morelia clastolepis'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316370384', child_labels=[], edge_length=10.4630678834, taxon_label='Morelia amethistina'),
            'Node4316371280' : NodeRelationship(parent_label='Node4316369488', child_labels=['Node4316371536','Node4316372048'], edge_length=7.66239308607, taxon_label=None),
            'Node4316371536' : NodeRelationship(parent_label='Node4316371280', child_labels=['Morelia carinata','Node4316371408'], edge_length=5.79725116231, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316371536', child_labels=[], edge_length=22.4623407504, taxon_label='Morelia carinata'),
            'Node4316371408' : NodeRelationship(parent_label='Node4316371536', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=5.91246822929, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316371408', child_labels=[], edge_length=16.5498725211, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316371408', child_labels=[], edge_length=16.5498725211, taxon_label='Morelia viridisN'),
            'Node4316372048' : NodeRelationship(parent_label='Node4316371280', child_labels=['Antaresia maculosa','Node4316372176'], edge_length=3.71811515781, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316372048', child_labels=[], edge_length=24.5414767549, taxon_label='Antaresia maculosa'),
            'Node4316372176' : NodeRelationship(parent_label='Node4316372048', child_labels=['Node4316372496','Antaresia perthensis'], edge_length=5.68085904285, taxon_label=None),
            'Node4316372496' : NodeRelationship(parent_label='Node4316372176', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=11.551508813, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316372496', child_labels=[], edge_length=7.30910889905, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316372496', child_labels=[], edge_length=7.30910889905, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316372176', child_labels=[], edge_length=18.8606177121, taxon_label='Antaresia perthensis'),
            'Node4316372752' : NodeRelationship(parent_label='Node4316369360', child_labels=['Morelia boeleni','Node4316393616'], edge_length=1.09335779327, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316372752', child_labels=[], edge_length=37.5960808765, taxon_label='Morelia boeleni'),
            'Node4316393616' : NodeRelationship(parent_label='Node4316372752', child_labels=['Node4316393936','Node4316395408'], edge_length=2.73294846613, taxon_label=None),
            'Node4316393936' : NodeRelationship(parent_label='Node4316393616', child_labels=['Node4316394064','Node4316394512'], edge_length=2.22916081797, taxon_label=None),
            'Node4316394064' : NodeRelationship(parent_label='Node4316393936', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=19.1660901297, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316394064', child_labels=[], edge_length=13.4678814627, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316394064', child_labels=[], edge_length=13.4678814627, taxon_label='Antaresia ramsayi'),
            'Node4316394512' : NodeRelationship(parent_label='Node4316393936', child_labels=['Node4316394640','Node4316395024'], edge_length=7.55440984409, taxon_label=None),
            'Node4316394640' : NodeRelationship(parent_label='Node4316394512', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=18.9933631628, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316394640', child_labels=[], edge_length=6.08619858554, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316394640', child_labels=[], edge_length=6.08619858554, taxon_label='Liasis mackloti'),
            'Node4316395024' : NodeRelationship(parent_label='Node4316394512', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=5.77339164759, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316395024', child_labels=[], edge_length=19.3061701007, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316395024', child_labels=[], edge_length=19.3061701007, taxon_label='Liasis olivaceus'),
            'Node4316395408' : NodeRelationship(parent_label='Node4316393616', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=13.2144705479, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316395408', child_labels=[], edge_length=21.6486618625, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316395408', child_labels=[], edge_length=21.6486618625, taxon_label='Liasis albertisii'),
            'Node4316395664' : NodeRelationship(parent_label='Node4316369232', child_labels=['Python timoriensis','Python reticulatus'], edge_length=20.035415736, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316395664', child_labels=[], edge_length=25.9434669707, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316395664', child_labels=[], edge_length=25.9434669707, taxon_label='Python reticulatus'),
            'Node4316396048' : NodeRelationship(parent_label='Node4316369104', child_labels=['Node4316396304','Python regius'], edge_length=19.1477140595, taxon_label=None),
            'Node4316396304' : NodeRelationship(parent_label='Node4316396048', child_labels=['Python curtus','Node4316396176'], edge_length=7.07244422644, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316396304', child_labels=[], edge_length=30.462733423, taxon_label='Python curtus'),
            'Node4316396176' : NodeRelationship(parent_label='Node4316396304', child_labels=['Python sebae','Python molurus'], edge_length=14.4432592042, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316396176', child_labels=[], edge_length=16.0194742188, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316396176', child_labels=[], edge_length=16.0194742188, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316396048', child_labels=[], edge_length=37.5351776494, taxon_label='Python regius'),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316368976', child_labels=[], edge_length=77.0175580149, taxon_label='Loxocemus bicolor'),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316348304', child_labels=[], edge_length=88.196508472, taxon_label='Xenopeltis unicolor'),
        },
        {
            'Node4316397456' : NodeRelationship(parent_label=None, child_labels=['Candola aspera','Node4316418320'], edge_length=None, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4316397456', child_labels=[], edge_length=95.8502441646, taxon_label='Candola aspera'),
            'Node4316418320' : NodeRelationship(parent_label='Node4316397456', child_labels=['Xenopeltis unicolor','Node4316418256'], edge_length=21.8741644934, taxon_label=None),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316418320', child_labels=[], edge_length=73.9760796713, taxon_label='Xenopeltis unicolor'),
            'Node4316418256' : NodeRelationship(parent_label='Node4316418320', child_labels=['Loxocemus bicolor','Node4316418576'], edge_length=9.52951598189, taxon_label=None),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316418256', child_labels=[], edge_length=64.4465636894, taxon_label='Loxocemus bicolor'),
            'Node4316418576' : NodeRelationship(parent_label='Node4316418256', child_labels=['Node4316418960','Node4316445776'], edge_length=22.9882151659, taxon_label=None),
            'Node4316418960' : NodeRelationship(parent_label='Node4316418576', child_labels=['Node4316419088','Node4316445392'], edge_length=6.43697452247, taxon_label=None),
            'Node4316419088' : NodeRelationship(parent_label='Node4316418960', child_labels=['Node4316397392','Node4316443088'], edge_length=6.7231540226, taxon_label=None),
            'Node4316397392' : NodeRelationship(parent_label='Node4316419088', child_labels=['Node4316419216','Node4316421136'], edge_length=1.03257897787, taxon_label=None),
            'Node4316419216' : NodeRelationship(parent_label='Node4316397392', child_labels=['Node4316419344','Morelia oenpelliensis'], edge_length=7.52518318022, taxon_label=None),
            'Node4316419344' : NodeRelationship(parent_label='Node4316419216', child_labels=['Node4316419472','Node4316419856'], edge_length=2.42060556338, taxon_label=None),
            'Node4316419472' : NodeRelationship(parent_label='Node4316419344', child_labels=['Morelia spilota','Morelia bredli'], edge_length=8.29264517113, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316419472', child_labels=[], edge_length=9.02720708579, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316419472', child_labels=[], edge_length=9.02720708579, taxon_label='Morelia bredli'),
            'Node4316419856' : NodeRelationship(parent_label='Node4316419344', child_labels=['Morelia tracyae','Node4316419728'], edge_length=6.70163613113, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316419856', child_labels=[], edge_length=10.6182161258, taxon_label='Morelia tracyae'),
            'Node4316419728' : NodeRelationship(parent_label='Node4316419856', child_labels=['Node4316420240','Morelia amethistina'], edge_length=3.47880840545, taxon_label=None),
            'Node4316420240' : NodeRelationship(parent_label='Node4316419728', child_labels=['Morelia clastolepis','Node4316420112'], edge_length=4.27223311967, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316420240', child_labels=[], edge_length=2.86717460067, taxon_label='Morelia clastolepis'),
            'Node4316420112' : NodeRelationship(parent_label='Node4316420240', child_labels=['Morelia nauta','Morelia kinghorni'], edge_length=0.371215520464, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316420112', child_labels=[], edge_length=2.49595908021, taxon_label='Morelia nauta'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316420112', child_labels=[], edge_length=2.49595908021, taxon_label='Morelia kinghorni'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316419728', child_labels=[], edge_length=7.13940772034, taxon_label='Morelia amethistina'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316419216', child_labels=[], edge_length=19.7404578203, taxon_label='Morelia oenpelliensis'),
            'Node4316421136' : NodeRelationship(parent_label='Node4316397392', child_labels=['Node4316421264','Node4316421776'], edge_length=2.39608223465, taxon_label=None),
            'Node4316421264' : NodeRelationship(parent_label='Node4316421136', child_labels=['Morelia carinata','Node4316420880'], edge_length=3.24355281662, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316421264', child_labels=[], edge_length=21.6260059492, taxon_label='Morelia carinata'),
            'Node4316420880' : NodeRelationship(parent_label='Node4316421264', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=9.1533399099, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316420880', child_labels=[], edge_length=12.4726660393, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316420880', child_labels=[], edge_length=12.4726660393, taxon_label='Morelia viridisN'),
            'Node4316421776' : NodeRelationship(parent_label='Node4316421136', child_labels=['Antaresia maculosa','Node4316421904'], edge_length=3.2027013644, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316421776', child_labels=[], edge_length=21.6668574015, taxon_label='Antaresia maculosa'),
            'Node4316421904' : NodeRelationship(parent_label='Node4316421776', child_labels=['Node4316442832','Antaresia perthensis'], edge_length=4.33602073619, taxon_label=None),
            'Node4316442832' : NodeRelationship(parent_label='Node4316421904', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=11.5233229214, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316442832', child_labels=[], edge_length=5.8075137439, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316442832', child_labels=[], edge_length=5.8075137439, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316421904', child_labels=[], edge_length=17.3308366653, taxon_label='Antaresia perthensis'),
            'Node4316443088' : NodeRelationship(parent_label='Node4316419088', child_labels=['Morelia boeleni','Node4316443344'], edge_length=1.21568464049, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316443088', child_labels=[], edge_length=27.0825353379, taxon_label='Morelia boeleni'),
            'Node4316443344' : NodeRelationship(parent_label='Node4316443088', child_labels=['Node4316443728','Node4316445136'], edge_length=1.37229916019, taxon_label=None),
            'Node4316443728' : NodeRelationship(parent_label='Node4316443344', child_labels=['Node4316443856','Node4316444240'], edge_length=2.64946637554, taxon_label=None),
            'Node4316443856' : NodeRelationship(parent_label='Node4316443728', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=13.5545767859, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316443856', child_labels=[], edge_length=9.50619301624, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316443856', child_labels=[], edge_length=9.50619301624, taxon_label='Antaresia ramsayi'),
            'Node4316444240' : NodeRelationship(parent_label='Node4316443728', child_labels=['Node4316444368','Node4316444752'], edge_length=4.67390676307, taxon_label=None),
            'Node4316444368' : NodeRelationship(parent_label='Node4316444240', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=12.8995814401, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316444368', child_labels=[], edge_length=5.48728159904, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316444368', child_labels=[], edge_length=5.48728159904, taxon_label='Liasis mackloti'),
            'Node4316444752' : NodeRelationship(parent_label='Node4316444240', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=1.38849394051, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316444752', child_labels=[], edge_length=16.9983690986, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316444752', child_labels=[], edge_length=16.9983690986, taxon_label='Liasis olivaceus'),
            'Node4316445136' : NodeRelationship(parent_label='Node4316443344', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=11.4050202795, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316445136', child_labels=[], edge_length=14.3052158982, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316445136', child_labels=[], edge_length=14.3052158982, taxon_label='Liasis albertisii'),
            'Node4316445392' : NodeRelationship(parent_label='Node4316418960', child_labels=['Python timoriensis','Python reticulatus'], edge_length=15.1871703268, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316445392', child_labels=[], edge_length=19.8342036742, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316445392', child_labels=[], edge_length=19.8342036742, taxon_label='Python reticulatus'),
            'Node4316445776' : NodeRelationship(parent_label='Node4316418576', child_labels=['Node4316446032','Python regius'], edge_length=13.8002436509, taxon_label=None),
            'Node4316446032' : NodeRelationship(parent_label='Node4316445776', child_labels=['Python curtus','Node4316445904'], edge_length=3.93234771773, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316446032', child_labels=[], edge_length=23.7257571548, taxon_label='Python curtus'),
            'Node4316445904' : NodeRelationship(parent_label='Node4316446032', child_labels=['Python sebae','Python molurus'], edge_length=7.63118798191, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316445904', child_labels=[], edge_length=16.0945691729, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316445904', child_labels=[], edge_length=16.0945691729, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316445776', child_labels=[], edge_length=27.6581048725, taxon_label='Python regius'),
        },
        {
            'Node4316467472' : NodeRelationship(parent_label=None, child_labels=['Candola aspera','Node4316467344'], edge_length=None, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4316467472', child_labels=[], edge_length=126.147943419, taxon_label='Candola aspera'),
            'Node4316467344' : NodeRelationship(parent_label='Node4316467472', child_labels=['Node4316467792','Xenopeltis unicolor'], edge_length=40.5157695372, taxon_label=None),
            'Node4316467792' : NodeRelationship(parent_label='Node4316467344', child_labels=['Node4316467920','Loxocemus bicolor'], edge_length=13.6608326978, taxon_label=None),
            'Node4316467920' : NodeRelationship(parent_label='Node4316467792', child_labels=['Node4316468048','Node4316490896'], edge_length=17.1918011574, taxon_label=None),
            'Node4316468048' : NodeRelationship(parent_label='Node4316467920', child_labels=['Node4316468176','Node4316490576'], edge_length=5.77299961102, taxon_label=None),
            'Node4316468176' : NodeRelationship(parent_label='Node4316468048', child_labels=['Node4316468304','Node4316470480'], edge_length=13.7010200713, taxon_label=None),
            'Node4316468304' : NodeRelationship(parent_label='Node4316468176', child_labels=['Morelia boeleni','Node4316467408'], edge_length=0.617658990864, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316468304', child_labels=[], edge_length=34.6878613533, taxon_label='Morelia boeleni'),
            'Node4316467408' : NodeRelationship(parent_label='Node4316468304', child_labels=['Node4316468688','Node4316469072'], edge_length=2.23641376685, taxon_label=None),
            'Node4316468688' : NodeRelationship(parent_label='Node4316467408', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=10.5685372807, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316468688', child_labels=[], edge_length=21.8829103058, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316468688', child_labels=[], edge_length=21.8829103058, taxon_label='Liasis albertisii'),
            'Node4316469072' : NodeRelationship(parent_label='Node4316467408', child_labels=['Node4316469200','Node4316469584'], edge_length=2.11685229318, taxon_label=None),
            'Node4316469200' : NodeRelationship(parent_label='Node4316469072', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=17.1757724208, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316469200', child_labels=[], edge_length=13.1588228725, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316469200', child_labels=[], edge_length=13.1588228725, taxon_label='Antaresia ramsayi'),
            'Node4316469584' : NodeRelationship(parent_label='Node4316469072', child_labels=['Node4316469712','Node4316470096'], edge_length=3.96372676423, taxon_label=None),
            'Node4316469712' : NodeRelationship(parent_label='Node4316469584', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=19.5683902852, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316469712', child_labels=[], edge_length=6.80247824391, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316469712', child_labels=[], edge_length=6.80247824391, taxon_label='Liasis mackloti'),
            'Node4316470096' : NodeRelationship(parent_label='Node4316469584', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=4.82785669688, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316470096', child_labels=[], edge_length=21.5430118322, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316470096', child_labels=[], edge_length=21.5430118322, taxon_label='Liasis olivaceus'),
            'Node4316470480' : NodeRelationship(parent_label='Node4316468176', child_labels=['Node4316470608','Node4316488848'], edge_length=1.35670146604, taxon_label=None),
            'Node4316470608' : NodeRelationship(parent_label='Node4316470480', child_labels=['Node4316470736','Node4316488592'], edge_length=10.3401216686, taxon_label=None),
            'Node4316470736' : NodeRelationship(parent_label='Node4316470608', child_labels=['Node4316470864','Morelia oenpelliensis'], edge_length=2.31259127041, taxon_label=None),
            'Node4316470864' : NodeRelationship(parent_label='Node4316470736', child_labels=['Morelia tracyae','Node4316470352'], edge_length=7.79402846351, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316470864', child_labels=[], edge_length=13.5020774756, taxon_label='Morelia tracyae'),
            'Node4316470352' : NodeRelationship(parent_label='Node4316470864', child_labels=['Node4316471248','Morelia amethistina'], edge_length=3.5072599877, taxon_label=None),
            'Node4316471248' : NodeRelationship(parent_label='Node4316470352', child_labels=['Morelia nauta','Node4316471120'], edge_length=6.20487512451, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316471248', child_labels=[], edge_length=3.78994236341, taxon_label='Morelia nauta'),
            'Node4316471120' : NodeRelationship(parent_label='Node4316471248', child_labels=['Morelia clastolepis','Morelia kinghorni'], edge_length=1.25204312348, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316471120', child_labels=[], edge_length=2.53789923993, taxon_label='Morelia clastolepis'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316471120', child_labels=[], edge_length=2.53789923993, taxon_label='Morelia kinghorni'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316470352', child_labels=[], edge_length=9.99481748791, taxon_label='Morelia amethistina'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316470736', child_labels=[], edge_length=21.2961059391, taxon_label='Morelia oenpelliensis'),
            'Node4316488592' : NodeRelationship(parent_label='Node4316470608', child_labels=['Morelia spilota','Morelia bredli'], edge_length=13.9826784484, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316488592', child_labels=[], edge_length=9.62601876119, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316488592', child_labels=[], edge_length=9.62601876119, taxon_label='Morelia bredli'),
            'Node4316488848' : NodeRelationship(parent_label='Node4316470480', child_labels=['Node4316489104','Node4316489616'], edge_length=7.03488295134, taxon_label=None),
            'Node4316489104' : NodeRelationship(parent_label='Node4316488848', child_labels=['Morelia carinata','Node4316488976'], edge_length=4.39672345568, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316489104', child_labels=[], edge_length=22.5172124711, taxon_label='Morelia carinata'),
            'Node4316488976' : NodeRelationship(parent_label='Node4316489104', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=11.9061835258, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316488976', child_labels=[], edge_length=10.6110289453, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316488976', child_labels=[], edge_length=10.6110289453, taxon_label='Morelia viridisN'),
            'Node4316489616' : NodeRelationship(parent_label='Node4316488848', child_labels=['Antaresia maculosa','Node4316489744'], edge_length=3.94778865925, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316489616', child_labels=[], edge_length=22.9661472676, taxon_label='Antaresia maculosa'),
            'Node4316489744' : NodeRelationship(parent_label='Node4316489616', child_labels=['Node4316490128','Antaresia perthensis'], edge_length=3.7713704432, taxon_label=None),
            'Node4316490128' : NodeRelationship(parent_label='Node4316489744', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=11.9409692012, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316490128', child_labels=[], edge_length=7.25380762314, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316490128', child_labels=[], edge_length=7.25380762314, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316489744', child_labels=[], edge_length=19.1947768244, taxon_label='Antaresia perthensis'),
            'Node4316490576' : NodeRelationship(parent_label='Node4316468048', child_labels=['Python timoriensis','Python reticulatus'], edge_length=24.6978541753, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316490576', child_labels=[], edge_length=24.3086862402, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316490576', child_labels=[], edge_length=24.3086862402, taxon_label='Python reticulatus'),
            'Node4316490896' : NodeRelationship(parent_label='Node4316467920', child_labels=['Node4316491152','Python regius'], edge_length=21.166651853, taxon_label=None),
            'Node4316491152' : NodeRelationship(parent_label='Node4316490896', child_labels=['Python curtus','Node4316491024'], edge_length=5.80983351578, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316491152', child_labels=[], edge_length=27.8030546577, taxon_label='Python curtus'),
            'Node4316491024' : NodeRelationship(parent_label='Node4316491152', child_labels=['Python sebae','Python molurus'], edge_length=10.4395768579, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316491024', child_labels=[], edge_length=17.3634777999, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316491024', child_labels=[], edge_length=17.3634777999, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316490896', child_labels=[], edge_length=33.6128881735, taxon_label='Python regius'),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316467792', child_labels=[], edge_length=71.971341184, taxon_label='Loxocemus bicolor'),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316467344', child_labels=[], edge_length=85.6321738818, taxon_label='Xenopeltis unicolor'),
        },
        {
            'Node4316512848' : NodeRelationship(parent_label=None, child_labels=['Candola aspera','Node4316513168'], edge_length=None, taxon_label=None),
            'Candola aspera' : NodeRelationship(parent_label='Node4316512848', child_labels=[], edge_length=146.770054852, taxon_label='Candola aspera'),
            'Node4316513168' : NodeRelationship(parent_label='Node4316512848', child_labels=['Node4316513296','Xenopeltis unicolor'], edge_length=49.9930471528, taxon_label=None),
            'Node4316513296' : NodeRelationship(parent_label='Node4316513168', child_labels=['Loxocemus bicolor','Node4316513104'], edge_length=13.4634525107, taxon_label=None),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316513296', child_labels=[], edge_length=83.313555189, taxon_label='Loxocemus bicolor'),
            'Node4316513104' : NodeRelationship(parent_label='Node4316513296', child_labels=['Node4316513680','Node4316536464'], edge_length=30.4776798161, taxon_label=None),
            'Node4316513680' : NodeRelationship(parent_label='Node4316513104', child_labels=['Node4316513808','Node4316512784'], edge_length=8.15039409252, taxon_label=None),
            'Node4316513808' : NodeRelationship(parent_label='Node4316513680', child_labels=['Node4316513936','Node4316534096'], edge_length=9.51023675819, taxon_label=None),
            'Node4316513936' : NodeRelationship(parent_label='Node4316513808', child_labels=['Node4316514064','Node4316515728'], edge_length=1.94845077373, taxon_label=None),
            'Node4316514064' : NodeRelationship(parent_label='Node4316513936', child_labels=['Node4316514192','Node4316514704'], edge_length=3.25905094181, taxon_label=None),
            'Node4316514192' : NodeRelationship(parent_label='Node4316514064', child_labels=['Morelia carinata','Node4316513552'], edge_length=4.82232736143, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316514192', child_labels=[], edge_length=25.1454154452, taxon_label='Morelia carinata'),
            'Node4316513552' : NodeRelationship(parent_label='Node4316514192', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=10.5788836702, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316513552', child_labels=[], edge_length=14.566531775, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316513552', child_labels=[], edge_length=14.566531775, taxon_label='Morelia viridisN'),
            'Node4316514704' : NodeRelationship(parent_label='Node4316514064', child_labels=['Antaresia maculosa','Node4316514832'], edge_length=4.58004825968, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316514704', child_labels=[], edge_length=25.387694547, taxon_label='Antaresia maculosa'),
            'Node4316514832' : NodeRelationship(parent_label='Node4316514704', child_labels=['Node4316515216','Antaresia perthensis'], edge_length=5.75635440071, taxon_label=None),
            'Node4316515216' : NodeRelationship(parent_label='Node4316514832', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=12.7115868187, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316515216', child_labels=[], edge_length=6.91975332756, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316515216', child_labels=[], edge_length=6.91975332756, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316514832', child_labels=[], edge_length=19.6313401463, taxon_label='Antaresia perthensis'),
            'Node4316515728' : NodeRelationship(parent_label='Node4316513936', child_labels=['Node4316515856','Node4316516240'], edge_length=10.0348009179, taxon_label=None),
            'Node4316515856' : NodeRelationship(parent_label='Node4316515728', child_labels=['Morelia spilota','Morelia bredli'], edge_length=14.4811675935, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316515856', child_labels=[], edge_length=8.71082523711, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316515856', child_labels=[], edge_length=8.71082523711, taxon_label='Morelia bredli'),
            'Node4316516240' : NodeRelationship(parent_label='Node4316515728', child_labels=['Node4316532816','Morelia oenpelliensis'], edge_length=0.795030531054, taxon_label=None),
            'Node4316532816' : NodeRelationship(parent_label='Node4316516240', child_labels=['Morelia tracyae','Node4316516112'], edge_length=7.95009719313, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316532816', child_labels=[], edge_length=14.4468651064, taxon_label='Morelia tracyae'),
            'Node4316516112' : NodeRelationship(parent_label='Node4316532816', child_labels=['Node4316533200','Morelia amethistina'], edge_length=2.27822479328, taxon_label=None),
            'Node4316533200' : NodeRelationship(parent_label='Node4316516112', child_labels=['Morelia clastolepis','Node4316533072'], edge_length=6.98040063458, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316533200', child_labels=[], edge_length=5.18823967855, taxon_label='Morelia clastolepis'),
            'Node4316533072' : NodeRelationship(parent_label='Node4316533200', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=1.71064757594, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316533072', child_labels=[], edge_length=3.47759210261, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316533072', child_labels=[], edge_length=3.47759210261, taxon_label='Morelia nauta'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316516112', child_labels=[], edge_length=12.1686403131, taxon_label='Morelia amethistina'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316516240', child_labels=[], edge_length=22.3969622995, taxon_label='Morelia oenpelliensis'),
            'Node4316534096' : NodeRelationship(parent_label='Node4316513808', child_labels=['Morelia boeleni','Node4316533968'], edge_length=1.20705242999, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316534096', child_labels=[], edge_length=33.9681920922, taxon_label='Morelia boeleni'),
            'Node4316533968' : NodeRelationship(parent_label='Node4316534096', child_labels=['Node4316534480','Node4316534864'], edge_length=1.28486469944, taxon_label=None),
            'Node4316534480' : NodeRelationship(parent_label='Node4316533968', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=12.4520799939, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316534480', child_labels=[], edge_length=20.2312473989, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316534480', child_labels=[], edge_length=20.2312473989, taxon_label='Liasis albertisii'),
            'Node4316534864' : NodeRelationship(parent_label='Node4316533968', child_labels=['Node4316534992','Node4316535376'], edge_length=1.68023264943, taxon_label=None),
            'Node4316534992' : NodeRelationship(parent_label='Node4316534864', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=19.0383478987, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316534992', child_labels=[], edge_length=11.9647468446, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316534992', child_labels=[], edge_length=11.9647468446, taxon_label='Antaresia ramsayi'),
            'Node4316535376' : NodeRelationship(parent_label='Node4316534864', child_labels=['Node4316535504','Node4316535888'], edge_length=4.75943584051, taxon_label=None),
            'Node4316535504' : NodeRelationship(parent_label='Node4316535376', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=17.1180008393, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316535504', child_labels=[], edge_length=9.12565806357, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316535504', child_labels=[], edge_length=9.12565806357, taxon_label='Liasis mackloti'),
            'Node4316535888' : NodeRelationship(parent_label='Node4316535376', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=2.6817531508, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316535888', child_labels=[], edge_length=23.561905752, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316535888', child_labels=[], edge_length=23.561905752, taxon_label='Liasis olivaceus'),
            'Node4316512784' : NodeRelationship(parent_label='Node4316513680', child_labels=['Python timoriensis','Python reticulatus'], edge_length=25.3895116055, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316512784', child_labels=[], edge_length=19.2959696748, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316512784', child_labels=[], edge_length=19.2959696748, taxon_label='Python reticulatus'),
            'Node4316536464' : NodeRelationship(parent_label='Node4316513104', child_labels=['Node4316536656','Python regius'], edge_length=14.4274535052, taxon_label=None),
            'Node4316536656' : NodeRelationship(parent_label='Node4316536464', child_labels=['Python curtus','Node4316557520'], edge_length=6.97033769, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316536656', child_labels=[], edge_length=31.4380841777, taxon_label='Python curtus'),
            'Node4316557520' : NodeRelationship(parent_label='Node4316536656', child_labels=['Python sebae','Python molurus'], edge_length=13.5557299793, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316557520', child_labels=[], edge_length=17.8823541984, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316557520', child_labels=[], edge_length=17.8823541984, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316536464', child_labels=[], edge_length=38.4084218677, taxon_label='Python regius'),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316513168', child_labels=[], edge_length=96.7770076997, taxon_label='Xenopeltis unicolor'),
        },
        {
            'Node4316558352' : NodeRelationship(parent_label=None, child_labels=['Node4316558480','Candola aspera'], edge_length=None, taxon_label=None),
            'Node4316558480' : NodeRelationship(parent_label='Node4316558352', child_labels=['Xenopeltis unicolor','Node4316558800'], edge_length=71.3865451194, taxon_label=None),
            'Xenopeltis unicolor' : NodeRelationship(parent_label='Node4316558480', child_labels=[], edge_length=140.081627971, taxon_label='Xenopeltis unicolor'),
            'Node4316558800' : NodeRelationship(parent_label='Node4316558480', child_labels=['Loxocemus bicolor','Node4316558736'], edge_length=26.0234610563, taxon_label=None),
            'Loxocemus bicolor' : NodeRelationship(parent_label='Node4316558800', child_labels=[], edge_length=114.058166915, taxon_label='Loxocemus bicolor'),
            'Node4316558736' : NodeRelationship(parent_label='Node4316558800', child_labels=['Node4316559184','Node4316606544'], edge_length=40.1699176918, taxon_label=None),
            'Node4316559184' : NodeRelationship(parent_label='Node4316558736', child_labels=['Node4316559312','Node4316585744'], edge_length=12.7558559915, taxon_label=None),
            'Node4316559312' : NodeRelationship(parent_label='Node4316559184', child_labels=['Node4316559440','Node4316583696'], edge_length=12.1691499629, taxon_label=None),
            'Node4316559440' : NodeRelationship(parent_label='Node4316559312', child_labels=['Node4316559568','Node4316561232'], edge_length=4.4227398576, taxon_label=None),
            'Node4316559568' : NodeRelationship(parent_label='Node4316559440', child_labels=['Node4316559696','Node4316560208'], edge_length=5.98387013923, taxon_label=None),
            'Node4316559696' : NodeRelationship(parent_label='Node4316559568', child_labels=['Morelia carinata','Node4316559056'], edge_length=7.53311752135, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node4316559696', child_labels=[], edge_length=31.0235157502, taxon_label='Morelia carinata'),
            'Node4316559056' : NodeRelationship(parent_label='Node4316559696', child_labels=['Morelia viridisS','Morelia viridisN'], edge_length=12.7654788618, taxon_label=None),
            'Morelia viridisS' : NodeRelationship(parent_label='Node4316559056', child_labels=[], edge_length=18.2580368884, taxon_label='Morelia viridisS'),
            'Morelia viridisN' : NodeRelationship(parent_label='Node4316559056', child_labels=[], edge_length=18.2580368884, taxon_label='Morelia viridisN'),
            'Node4316560208' : NodeRelationship(parent_label='Node4316559568', child_labels=['Antaresia maculosa','Node4316560336'], edge_length=4.18557870369, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node4316560208', child_labels=[], edge_length=34.3710545678, taxon_label='Antaresia maculosa'),
            'Node4316560336' : NodeRelationship(parent_label='Node4316560208', child_labels=['Node4316560720','Antaresia perthensis'], edge_length=7.81168349566, taxon_label=None),
            'Node4316560720' : NodeRelationship(parent_label='Node4316560336', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=15.350690271, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node4316560720', child_labels=[], edge_length=11.2086808012, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node4316560720', child_labels=[], edge_length=11.2086808012, taxon_label='Antaresia stimsoni'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node4316560336', child_labels=[], edge_length=26.5593710722, taxon_label='Antaresia perthensis'),
            'Node4316561232' : NodeRelationship(parent_label='Node4316559440', child_labels=['Node4316561360','Node4316583312'], edge_length=17.9619739892, taxon_label=None),
            'Node4316561360' : NodeRelationship(parent_label='Node4316561232', child_labels=['Morelia oenpelliensis','Node4316560976'], edge_length=1.39777518086, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node4316561360', child_labels=[], edge_length=25.1807542407, taxon_label='Morelia oenpelliensis'),
            'Node4316560976' : NodeRelationship(parent_label='Node4316561360', child_labels=['Morelia tracyae','Node4316582160'], edge_length=7.6242060025, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node4316560976', child_labels=[], edge_length=17.5565482382, taxon_label='Morelia tracyae'),
            'Node4316582160' : NodeRelationship(parent_label='Node4316560976', child_labels=['Node4316582544','Morelia amethistina'], edge_length=3.73213849687, taxon_label=None),
            'Node4316582544' : NodeRelationship(parent_label='Node4316582160', child_labels=['Morelia clastolepis','Node4316582416'], edge_length=8.62088071739, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node4316582544', child_labels=[], edge_length=5.20352902397, taxon_label='Morelia clastolepis'),
            'Node4316582416' : NodeRelationship(parent_label='Node4316582544', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=2.83199057731, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node4316582416', child_labels=[], edge_length=2.37153844665, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node4316582416', child_labels=[], edge_length=2.37153844665, taxon_label='Morelia nauta'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node4316582160', child_labels=[], edge_length=13.8244097414, taxon_label='Morelia amethistina'),
            'Node4316583312' : NodeRelationship(parent_label='Node4316561232', child_labels=['Morelia spilota','Morelia bredli'], edge_length=13.6008570384, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node4316583312', child_labels=[], edge_length=12.9776723832, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node4316583312', child_labels=[], edge_length=12.9776723832, taxon_label='Morelia bredli'),
            'Node4316583696' : NodeRelationship(parent_label='Node4316559312', child_labels=['Node4316583824','Node4316585104'], edge_length=3.26334313874, taxon_label=None),
            'Node4316583824' : NodeRelationship(parent_label='Node4316583696', child_labels=['Node4316583952','Node4316584336'], edge_length=3.5026210108, taxon_label=None),
            'Node4316583952' : NodeRelationship(parent_label='Node4316583824', child_labels=['Antaresia melanocephalus','Antaresia ramsayi'], edge_length=23.3327851176, taxon_label=None),
            'Antaresia melanocephalus' : NodeRelationship(parent_label='Node4316583952', child_labels=[], edge_length=18.8644940012, taxon_label='Antaresia melanocephalus'),
            'Antaresia ramsayi' : NodeRelationship(parent_label='Node4316583952', child_labels=[], edge_length=18.8644940012, taxon_label='Antaresia ramsayi'),
            'Node4316584336' : NodeRelationship(parent_label='Node4316583824', child_labels=['Node4316584464','Node4316558288'], edge_length=12.8697548268, taxon_label=None),
            'Node4316584464' : NodeRelationship(parent_label='Node4316584336', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=22.0000981488, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node4316584464', child_labels=[], edge_length=7.32742614319, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node4316584464', child_labels=[], edge_length=7.32742614319, taxon_label='Liasis mackloti'),
            'Node4316558288' : NodeRelationship(parent_label='Node4316584336', child_labels=['Apodora papuana','Liasis olivaceus'], edge_length=5.24037108992, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node4316558288', child_labels=[], edge_length=24.0871532021, taxon_label='Apodora papuana'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node4316558288', child_labels=[], edge_length=24.0871532021, taxon_label='Liasis olivaceus'),
            'Node4316585104' : NodeRelationship(parent_label='Node4316583696', child_labels=['Node4316585232','Morelia boeleni'], edge_length=0.85376044557, taxon_label=None),
            'Node4316585232' : NodeRelationship(parent_label='Node4316585104', child_labels=['Bothrochilus boa','Liasis albertisii'], edge_length=19.3585494438, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node4316585232', child_labels=[], edge_length=25.4875902402, taxon_label='Bothrochilus boa'),
            'Liasis albertisii' : NodeRelationship(parent_label='Node4316585232', child_labels=[], edge_length=25.4875902402, taxon_label='Liasis albertisii'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node4316585104', child_labels=[], edge_length=44.846139684, taxon_label='Morelia boeleni'),
            'Node4316585744' : NodeRelationship(parent_label='Node4316559184', child_labels=['Python timoriensis','Python reticulatus'], edge_length=28.055060848, taxon_label=None),
            'Python timoriensis' : NodeRelationship(parent_label='Node4316585744', child_labels=[], edge_length=33.0773323832, taxon_label='Python timoriensis'),
            'Python reticulatus' : NodeRelationship(parent_label='Node4316585744', child_labels=[], edge_length=33.0773323832, taxon_label='Python reticulatus'),
            'Node4316606544' : NodeRelationship(parent_label='Node4316558736', child_labels=['Node4316606800','Python regius'], edge_length=17.7722054357, taxon_label=None),
            'Node4316606800' : NodeRelationship(parent_label='Node4316606544', child_labels=['Python curtus','Node4316606672'], edge_length=13.6364666177, taxon_label=None),
            'Python curtus' : NodeRelationship(parent_label='Node4316606800', child_labels=[], edge_length=42.4795771694, taxon_label='Python curtus'),
            'Node4316606672' : NodeRelationship(parent_label='Node4316606800', child_labels=['Python sebae','Python molurus'], edge_length=16.6495052056, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node4316606672', child_labels=[], edge_length=25.8300719638, taxon_label='Python sebae'),
            'Python molurus' : NodeRelationship(parent_label='Node4316606672', child_labels=[], edge_length=25.8300719638, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node4316606544', child_labels=[], edge_length=56.1160437871, taxon_label='Python regius'),
            'Candola aspera' : NodeRelationship(parent_label='Node4316558352', child_labels=[], edge_length=211.46817309, taxon_label='Candola aspera'),
        },
    ]
    return treelist_node_references

def reference_tree_list(taxon_set=None):
    tree_list = dendropy.TreeList(label=None, oid="TreeList4303123040", taxon_set=taxon_set)
    tax_4313741136 = tree_list.taxon_set.require_taxon(label="Antaresia childreni", oid="Taxon4313741136")
    tax_4313741328 = tree_list.taxon_set.require_taxon(label="Antaresia maculosa", oid="Taxon4313741328")
    tax_4313741456 = tree_list.taxon_set.require_taxon(label="Antaresia melanocephalus", oid="Taxon4313741456")
    tax_4313741584 = tree_list.taxon_set.require_taxon(label="Antaresia perthensis", oid="Taxon4313741584")
    tax_4313741712 = tree_list.taxon_set.require_taxon(label="Antaresia ramsayi", oid="Taxon4313741712")
    tax_4313741840 = tree_list.taxon_set.require_taxon(label="Antaresia stimsoni", oid="Taxon4313741840")
    tax_4313741904 = tree_list.taxon_set.require_taxon(label="Apodora papuana", oid="Taxon4313741904")
    tax_4313741968 = tree_list.taxon_set.require_taxon(label="Bothrochilus boa", oid="Taxon4313741968")
    tax_4313742032 = tree_list.taxon_set.require_taxon(label="Candola aspera", oid="Taxon4313742032")
    tax_4313742160 = tree_list.taxon_set.require_taxon(label="Liasis albertisii", oid="Taxon4313742160")
    tax_4313742224 = tree_list.taxon_set.require_taxon(label="Liasis fuscus", oid="Taxon4313742224")
    tax_4313742288 = tree_list.taxon_set.require_taxon(label="Liasis mackloti", oid="Taxon4313742288")
    tax_4313742352 = tree_list.taxon_set.require_taxon(label="Liasis olivaceus", oid="Taxon4313742352")
    tax_4313742480 = tree_list.taxon_set.require_taxon(label="Loxocemus bicolor", oid="Taxon4313742480")
    tax_4313742608 = tree_list.taxon_set.require_taxon(label="Morelia amethistina", oid="Taxon4313742608")
    tax_4313742672 = tree_list.taxon_set.require_taxon(label="Morelia boeleni", oid="Taxon4313742672")
    tax_4313742736 = tree_list.taxon_set.require_taxon(label="Morelia bredli", oid="Taxon4313742736")
    tax_4313742800 = tree_list.taxon_set.require_taxon(label="Morelia carinata", oid="Taxon4313742800")
    tax_4313742928 = tree_list.taxon_set.require_taxon(label="Morelia clastolepis", oid="Taxon4313742928")
    tax_4313743056 = tree_list.taxon_set.require_taxon(label="Morelia kinghorni", oid="Taxon4313743056")
    tax_4313743120 = tree_list.taxon_set.require_taxon(label="Morelia nauta", oid="Taxon4313743120")
    tax_4313743248 = tree_list.taxon_set.require_taxon(label="Morelia oenpelliensis", oid="Taxon4313743248")
    tax_4313743312 = tree_list.taxon_set.require_taxon(label="Morelia spilota", oid="Taxon4313743312")
    tax_4313759824 = tree_list.taxon_set.require_taxon(label="Morelia tracyae", oid="Taxon4313759824")
    tax_4313759888 = tree_list.taxon_set.require_taxon(label="Morelia viridisN", oid="Taxon4313759888")
    tax_4313759952 = tree_list.taxon_set.require_taxon(label="Morelia viridisS", oid="Taxon4313759952")
    tax_4313760016 = tree_list.taxon_set.require_taxon(label="Python curtus", oid="Taxon4313760016")
    tax_4313760080 = tree_list.taxon_set.require_taxon(label="Python molurus", oid="Taxon4313760080")
    tax_4313760144 = tree_list.taxon_set.require_taxon(label="Python regius", oid="Taxon4313760144")
    tax_4313760272 = tree_list.taxon_set.require_taxon(label="Python reticulatus", oid="Taxon4313760272")
    tax_4313760336 = tree_list.taxon_set.require_taxon(label="Python sebae", oid="Taxon4313760336")
    tax_4313760464 = tree_list.taxon_set.require_taxon(label="Python timoriensis", oid="Taxon4313760464")
    tax_4313760592 = tree_list.taxon_set.require_taxon(label="Xenopeltis unicolor", oid="Taxon4313760592")
    tree_4313760848 = dendropy.Tree(label="Tree01", taxon_set=tree_list.taxon_set, oid="Tree4313760848")
    tree_list.append(tree_4313760848, reindex_taxa=False)
    tree_4313760848.seed_node.oid = 'Node4313761104'
    nd_4313761232 = tree_4313760848.seed_node.new_child(label="Node4313761232", taxon=None, edge_length=78.6266408419, oid="Node4313761232")
    nd_4313761232.edge.oid = "Edge4313761296"
    nd_4313801936 = tree_4313760848.seed_node.new_child(label="Node4313801936", taxon=None, edge_length=229.880308935, oid="Node4313801936")
    nd_4313801936.edge.oid = "Edge4313802128"
    nd_4313761360 = nd_4313761232.new_child(label="Node4313761360", taxon=None, edge_length=123.936295332, oid="Node4313761360")
    nd_4313761360.edge.oid = "Edge4313761424"
    nd_4313801552 = nd_4313761232.new_child(label="Node4313801552", taxon=None, edge_length=202.422980755, oid="Node4313801552")
    nd_4313801552.edge.oid = "Edge4313801744"
    nd_4313761488 = nd_4313761360.new_child(label="Node4313761488", taxon=None, edge_length=53.7887032791, oid="Node4313761488")
    nd_4313761488.edge.oid = "Edge4313761552"
    nd_4313779728 = nd_4313761360.new_child(label="Node4313779728", taxon=None, edge_length=69.6091072013, oid="Node4313779728")
    nd_4313779728.edge.oid = "Edge4313780048"
    nd_4313761616 = nd_4313761488.new_child(label="Node4313761616", taxon=None, edge_length=4.50960412336, oid="Node4313761616")
    nd_4313761616.edge.oid = "Edge4313761680"
    nd_4313779088 = nd_4313761488.new_child(label="Node4313779088", taxon=None, edge_length=24.4333103788, oid="Node4313779088")
    nd_4313779088.edge.oid = "Edge4313779408"
    nd_4313761744 = nd_4313761616.new_child(label="Node4313761744", taxon=None, edge_length=1.29166787247, oid="Node4313761744")
    nd_4313761744.edge.oid = "Edge4313761808"
    nd_4313778704 = nd_4313761616.new_child(label="Node4313778704", taxon=None, edge_length=13.1562713928, oid="Node4313778704")
    nd_4313778704.edge.oid = "Edge4313778768"
    nd_4313761872 = nd_4313761744.new_child(label="Node4313761872", taxon=None, edge_length=0.739270951321, oid="Node4313761872")
    nd_4313761872.edge.oid = "Edge4313761936"
    nd_4313777808 = nd_4313761744.new_child(label="Node4313777808", taxon=None, edge_length=16.8726126715, oid="Node4313777808")
    nd_4313777808.edge.oid = "Edge4313777872"
    nd_4313762000 = nd_4313761872.new_child(label="Node4313762000", taxon=None, edge_length=4.00932181352, oid="Node4313762000")
    nd_4313762000.edge.oid = "Edge4313762064"
    nd_4313777552 = nd_4313761872.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=18.3098649933, oid="Node4313777552")
    nd_4313777552.edge.oid = "Edge4313777744"
    nd_4313762128 = nd_4313762000.new_child(label="Node4313762128", taxon=None, edge_length=10.4795105126, oid="Node4313762128")
    nd_4313762128.edge.oid = "Edge4313762192"
    nd_4313762704 = nd_4313762000.new_child(label="Node4313762704", taxon=None, edge_length=3.67431549248, oid="Node4313762704")
    nd_4313762704.edge.oid = "Edge4313762896"
    nd_4313762256 = nd_4313762128.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=3.82103266723, oid="Node4313762256")
    nd_4313762256.edge.oid = "Edge4313762320"
    nd_4313762448 = nd_4313762128.new_child(label="Node4313762448", taxon=None, edge_length=2.72293867732, oid="Node4313762448")
    nd_4313762448.edge.oid = "Edge4313762512"
    nd_4313762576 = nd_4313762448.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=1.09809398991, oid="Node4313762576")
    nd_4313762576.edge.oid = "Edge4313762640"
    nd_4313762384 = nd_4313762448.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=1.09809398991, oid="Node4313762384")
    nd_4313762384.edge.oid = "Edge4313762768"
    nd_4313762960 = nd_4313762704.new_child(label="Node4313762960", taxon=None, edge_length=4.47413907662, oid="Node4313762960")
    nd_4313762960.edge.oid = "Edge4313763024"
    nd_4313776720 = nd_4313762704.new_child(label="Node4313776720", taxon=None, edge_length=3.27824986702, oid="Node4313776720")
    nd_4313776720.edge.oid = "Edge4313776848"
    nd_4313763088 = nd_4313762960.new_child(label="Node4313763088", taxon=None, edge_length=1.10594701364, oid="Node4313763088")
    nd_4313763088.edge.oid = "Edge4313763152"
    nd_4313776208 = nd_4313762960.new_child(label="Node4313776208", taxon=None, edge_length=0.946531890581, oid="Node4313776208")
    nd_4313776208.edge.oid = "Edge4313776400"
    nd_4313763216 = nd_4313763088.new_child(label="Node4313763216", taxon=None, edge_length=3.56964311758, oid="Node4313763216")
    nd_4313763216.edge.oid = "Edge4313763280"
    nd_4313763600 = nd_4313763088.new_child(label="Node4313763600", taxon=None, edge_length=3.54960462938, oid="Node4313763600")
    nd_4313763600.edge.oid = "Edge4313763664"
    nd_4313763344 = nd_4313763216.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=1.47649847948, oid="Node4313763344")
    nd_4313763344.edge.oid = "Edge4313763408"
    nd_4313762832 = nd_4313763216.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=1.47649847948, oid="Node4313762832")
    nd_4313762832.edge.oid = "Edge4313763536"
    nd_4313763728 = nd_4313763600.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=1.49653696768, oid="Node4313763728")
    nd_4313763728.edge.oid = "Edge4313763792"
    nd_4313763472 = nd_4313763600.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=1.49653696768, oid="Node4313763472")
    nd_4313763472.edge.oid = "Edge4313776272"
    nd_4313776464 = nd_4313776208.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=5.20555672012, oid="Node4313776464")
    nd_4313776464.edge.oid = "Edge4313776528"
    nd_4313776336 = nd_4313776208.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=5.20555672012, oid="Node4313776336")
    nd_4313776336.edge.oid = "Edge4313776656"
    nd_4313776912 = nd_4313776720.new_child(label="Node4313776912", taxon=None, edge_length=2.88197852526, oid="Node4313776912")
    nd_4313776912.edge.oid = "Edge4313776976"
    nd_4313777296 = nd_4313776720.new_child(label="Node4313777296", taxon=None, edge_length=6.86415378064, oid="Node4313777296")
    nd_4313777296.edge.oid = "Edge4313777360"
    nd_4313777040 = nd_4313776912.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=4.46599929505, oid="Node4313777040")
    nd_4313777040.edge.oid = "Edge4313777104"
    nd_4313776784 = nd_4313776912.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=4.46599929505, oid="Node4313776784")
    nd_4313776784.edge.oid = "Edge4313777232"
    nd_4313777424 = nd_4313777296.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=0.483824039664, oid="Node4313777424")
    nd_4313777424.edge.oid = "Edge4313777488"
    nd_4313777168 = nd_4313777296.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=0.483824039664, oid="Node4313777168")
    nd_4313777168.edge.oid = "Edge4313777616"
    nd_4313777936 = nd_4313777808.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=2.17652327312, oid="Node4313777936")
    nd_4313777936.edge.oid = "Edge4313778000"
    nd_4313777680 = nd_4313777808.new_child(label="Node4313777680", taxon=None, edge_length=1.67230791531, oid="Node4313777680")
    nd_4313777680.edge.oid = "Edge4313778128"
    nd_4313778192 = nd_4313777680.new_child(label="Node4313778192", taxon=None, edge_length=0.491713738136, oid="Node4313778192")
    nd_4313778192.edge.oid = "Edge4313778256"
    nd_4313778576 = nd_4313777680.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=0.504215357803, oid="Node4313778576")
    nd_4313778576.edge.oid = "Edge4313778640"
    nd_4313778320 = nd_4313778192.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=0.0125016196671, oid="Node4313778320")
    nd_4313778320.edge.oid = "Edge4313778384"
    nd_4313778064 = nd_4313778192.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=0.0125016196671, oid="Node4313778064")
    nd_4313778064.edge.oid = "Edge4313778512"
    nd_4313778832 = nd_4313778704.new_child(label="Node4313778832", taxon=None, edge_length=2.98661848623, oid="Node4313778832")
    nd_4313778832.edge.oid = "Edge4313778896"
    nd_4313779216 = nd_4313778704.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=7.18453242432, oid="Node4313779216")
    nd_4313779216.edge.oid = "Edge4313779280"
    nd_4313778960 = nd_4313778832.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=4.19791393809, oid="Node4313778960")
    nd_4313778960.edge.oid = "Edge4313779024"
    nd_4313778448 = nd_4313778832.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=4.19791393809, oid="Node4313778448")
    nd_4313778448.edge.oid = "Edge4313779152"
    nd_4313779472 = nd_4313779088.new_child(label="Node4313779472", taxon=None, edge_length=0.207889736001, oid="Node4313779472")
    nd_4313779472.edge.oid = "Edge4313779536"
    nd_4313779856 = nd_4313779088.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=0.417097561686, oid="Node4313779856")
    nd_4313779856.edge.oid = "Edge4313779920"
    nd_4313779600 = nd_4313779472.new_child(label="Python regius", taxon=tax_4313760144, edge_length=0.209207825685, oid="Node4313779600")
    nd_4313779600.edge.oid = "Edge4313779664"
    nd_4313779344 = nd_4313779472.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=0.209207825685, oid="Node4313779344")
    nd_4313779344.edge.oid = "Edge4313779792"
    nd_4313780112 = nd_4313779728.new_child(label="Node4313780112", taxon=None, edge_length=1.24643505521, oid="Node4313780112")
    nd_4313780112.edge.oid = "Edge4313780176"
    nd_4313801424 = nd_4313779728.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=9.03000401821, oid="Node4313801424")
    nd_4313801424.edge.oid = "Edge4313801616"
    nd_4313800784 = nd_4313780112.new_child(label="Node4313800784", taxon=None, edge_length=2.6364754365, oid="Node4313800784")
    nd_4313800784.edge.oid = "Edge4313800848"
    nd_4313801168 = nd_4313780112.new_child(label="Node4313801168", taxon=None, edge_length=7.12573328141, oid="Node4313801168")
    nd_4313801168.edge.oid = "Edge4313801232"
    nd_4313800912 = nd_4313800784.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=5.1470935265, oid="Node4313800912")
    nd_4313800912.edge.oid = "Edge4313800976"
    nd_4313779984 = nd_4313800784.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=5.1470935265, oid="Node4313779984")
    nd_4313779984.edge.oid = "Edge4313801040"
    nd_4313801296 = nd_4313801168.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=0.657835681585, oid="Node4313801296")
    nd_4313801296.edge.oid = "Edge4313801360"
    nd_4313801104 = nd_4313801168.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=0.657835681585, oid="Node4313801104")
    nd_4313801104.edge.oid = "Edge4313801488"
    nd_4313801808 = nd_4313801552.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=0.152425796315, oid="Node4313801808")
    nd_4313801808.edge.oid = "Edge4313801872"
    nd_4313801680 = nd_4313801552.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=0.152425796315, oid="Node4313801680")
    nd_4313801680.edge.oid = "Edge4313802000"
    nd_4313802192 = nd_4313801936.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=51.3217384582, oid="Node4313802192")
    nd_4313802192.edge.oid = "Edge4313802256"
    nd_4313802064 = nd_4313801936.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=51.3217384582, oid="Node4313802064")
    nd_4313802064.edge.oid = "Edge4313802384"
    tree_4313802320 = dendropy.Tree(label="Tree02", taxon_set=tree_list.taxon_set, oid="Tree4313802320")
    tree_list.append(tree_4313802320, reindex_taxa=False)
    tree_4313802320.seed_node.oid = 'Node4313802576'
    nd_4313802704 = tree_4313802320.seed_node.new_child(label="Node4313802704", taxon=None, edge_length=18.8917197007, oid="Node4313802704")
    nd_4313802704.edge.oid = "Edge4313802768"
    nd_4316157840 = tree_4313802320.seed_node.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=100.189141925, oid="Node4316157840")
    nd_4316157840.edge.oid = "Edge4316158096"
    nd_4313802832 = nd_4313802704.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=81.2974222246, oid="Node4313802832")
    nd_4313802832.edge.oid = "Edge4313802896"
    nd_4313803024 = nd_4313802704.new_child(label="Node4313803024", taxon=None, edge_length=33.1565984398, oid="Node4313803024")
    nd_4313803024.edge.oid = "Edge4313803088"
    nd_4313803152 = nd_4313803024.new_child(label="Node4313803152", taxon=None, edge_length=6.57324583185, oid="Node4313803152")
    nd_4313803152.edge.oid = "Edge4313803216"
    nd_4316156752 = nd_4313803024.new_child(label="Node4316156752", taxon=None, edge_length=16.2594519516, oid="Node4316156752")
    nd_4316156752.edge.oid = "Edge4316156816"
    nd_4313803280 = nd_4313803152.new_child(label="Node4313803280", taxon=None, edge_length=0.76583222117, oid="Node4313803280")
    nd_4313803280.edge.oid = "Edge4313803344"
    nd_4313829328 = nd_4313803152.new_child(label="Node4313829328", taxon=None, edge_length=2.53377266123, oid="Node4313829328")
    nd_4313829328.edge.oid = "Edge4313829072"
    nd_4313803408 = nd_4313803280.new_child(label="Node4313803408", taxon=None, edge_length=4.5936111676, oid="Node4313803408")
    nd_4313803408.edge.oid = "Edge4313803472"
    nd_4313826384 = nd_4313803280.new_child(label="Node4313826384", taxon=None, edge_length=1.27553605821, oid="Node4313826384")
    nd_4313826384.edge.oid = "Edge4313826704"
    nd_4313803536 = nd_4313803408.new_child(label="Node4313803536", taxon=None, edge_length=15.1180336863, oid="Node4313803536")
    nd_4313803536.edge.oid = "Edge4313803600"
    nd_4313804432 = nd_4313803408.new_child(label="Node4313804432", taxon=None, edge_length=3.40951166184, oid="Node4313804432")
    nd_4313804432.edge.oid = "Edge4313804496"
    nd_4313803664 = nd_4313803536.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=21.0901008779, oid="Node4313803664")
    nd_4313803664.edge.oid = "Edge4313803728"
    nd_4313802960 = nd_4313803536.new_child(label="Node4313802960", taxon=None, edge_length=3.38653541663, oid="Node4313802960")
    nd_4313802960.edge.oid = "Edge4313803856"
    nd_4313803920 = nd_4313802960.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=17.7035654613, oid="Node4313803920")
    nd_4313803920.edge.oid = "Edge4313803984"
    nd_4313803792 = nd_4313802960.new_child(label="Node4313803792", taxon=None, edge_length=2.6244717729, oid="Node4313803792")
    nd_4313803792.edge.oid = "Edge4313804048"
    nd_4313804112 = nd_4313803792.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=15.0790936884, oid="Node4313804112")
    nd_4313804112.edge.oid = "Edge4313804176"
    nd_4313804304 = nd_4313803792.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=15.0790936884, oid="Node4313804304")
    nd_4313804304.edge.oid = "Edge4313804368"
    nd_4313804560 = nd_4313804432.new_child(label="Node4313804560", taxon=None, edge_length=5.00375936144, oid="Node4313804560")
    nd_4313804560.edge.oid = "Edge4313804624"
    nd_4313825616 = nd_4313804432.new_child(label="Node4313825616", taxon=None, edge_length=5.81119736053, oid="Node4313825616")
    nd_4313825616.edge.oid = "Edge4313825808"
    nd_4313804688 = nd_4313804560.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=27.7948635409, oid="Node4313804688")
    nd_4313804688.edge.oid = "Edge4313804752"
    nd_4313804240 = nd_4313804560.new_child(label="Node4313804240", taxon=None, edge_length=8.32618746237, oid="Node4313804240")
    nd_4313804240.edge.oid = "Edge4313825424"
    nd_4313825488 = nd_4313804240.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=19.4686760786, oid="Node4313825488")
    nd_4313825488.edge.oid = "Edge4313825552"
    nd_4313825360 = nd_4313804240.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=19.4686760786, oid="Node4313825360")
    nd_4313825360.edge.oid = "Edge4313825680"
    nd_4313825872 = nd_4313825616.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=26.9874255418, oid="Node4313825872")
    nd_4313825872.edge.oid = "Edge4313825936"
    nd_4313825744 = nd_4313825616.new_child(label="Node4313825744", taxon=None, edge_length=2.25683638168, oid="Node4313825744")
    nd_4313825744.edge.oid = "Edge4313826064"
    nd_4313826128 = nd_4313825744.new_child(label="Node4313826128", taxon=None, edge_length=16.6530983052, oid="Node4313826128")
    nd_4313826128.edge.oid = "Edge4313826192"
    nd_4313826512 = nd_4313825744.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=24.7305891602, oid="Node4313826512")
    nd_4313826512.edge.oid = "Edge4313826576"
    nd_4313826256 = nd_4313826128.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=8.07749085501, oid="Node4313826256")
    nd_4313826256.edge.oid = "Edge4313826320"
    nd_4313826000 = nd_4313826128.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=8.07749085501, oid="Node4313826000")
    nd_4313826000.edge.oid = "Edge4313826448"
    nd_4313826768 = nd_4313826384.new_child(label="Node4313826768", taxon=None, edge_length=4.33214615343, oid="Node4313826768")
    nd_4313826768.edge.oid = "Edge4313826832"
    nd_4313828048 = nd_4313826384.new_child(label="Node4313828048", taxon=None, edge_length=10.9652932592, oid="Node4313828048")
    nd_4313828048.edge.oid = "Edge4313828240"
    nd_4313826896 = nd_4313826768.new_child(label="Node4313826896", taxon=None, edge_length=3.37363071467, oid="Node4313826896")
    nd_4313826896.edge.oid = "Edge4313826960"
    nd_4313827664 = nd_4313826768.new_child(label="Node4313827664", taxon=None, edge_length=16.3762764593, oid="Node4313827664")
    nd_4313827664.edge.oid = "Edge4313827856"
    nd_4313827024 = nd_4313826896.new_child(label="Node4313827024", taxon=None, edge_length=26.1365403684, oid="Node4313827024")
    nd_4313827024.edge.oid = "Edge4313827088"
    nd_4313827408 = nd_4313826896.new_child(label="Node4313827408", taxon=None, edge_length=12.1064068345, oid="Node4313827408")
    nd_4313827408.edge.oid = "Edge4313827472"
    nd_4313827152 = nd_4313827024.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=5.68389243709, oid="Node4313827152")
    nd_4313827152.edge.oid = "Edge4313827216"
    nd_4313826640 = nd_4313827024.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=5.68389243709, oid="Node4313826640")
    nd_4313826640.edge.oid = "Edge4313827344"
    nd_4313827536 = nd_4313827408.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=19.714025971, oid="Node4313827536")
    nd_4313827536.edge.oid = "Edge4313827600"
    nd_4313827280 = nd_4313827408.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=19.714025971, oid="Node4313827280")
    nd_4313827280.edge.oid = "Edge4313827728"
    nd_4313827920 = nd_4313827664.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=18.8177870609, oid="Node4313827920")
    nd_4313827920.edge.oid = "Edge4313827984"
    nd_4313827792 = nd_4313827664.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=18.8177870609, oid="Node4313827792")
    nd_4313827792.edge.oid = "Edge4313828112"
    nd_4313828304 = nd_4313828048.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=28.5609164144, oid="Node4313828304")
    nd_4313828304.edge.oid = "Edge4313828368"
    nd_4313828176 = nd_4313828048.new_child(label="Node4313828176", taxon=None, edge_length=11.3916491298, oid="Node4313828176")
    nd_4313828176.edge.oid = "Edge4313828496"
    nd_4313828560 = nd_4313828176.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=17.1692672846, oid="Node4313828560")
    nd_4313828560.edge.oid = "Edge4313828624"
    nd_4313828432 = nd_4313828176.new_child(label="Node4313828432", taxon=None, edge_length=10.4784522084, oid="Node4313828432")
    nd_4313828432.edge.oid = "Edge4313828752"
    nd_4313828816 = nd_4313828432.new_child(label="Node4313828816", taxon=None, edge_length=1.44575855569, oid="Node4313828816")
    nd_4313828816.edge.oid = "Edge4313828880"
    nd_4313829200 = nd_4313828432.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=6.69081507619, oid="Node4313829200")
    nd_4313829200.edge.oid = "Edge4313829264"
    nd_4313828944 = nd_4313828816.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=5.2450565205, oid="Node4313828944")
    nd_4313828944.edge.oid = "Edge4313829008"
    nd_4313828688 = nd_4313828816.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=5.2450565205, oid="Node4313828688")
    nd_4313828688.edge.oid = "Edge4313829136"
    nd_4316155984 = nd_4313829328.new_child(label="Node4316155984", taxon=None, edge_length=23.8923728386, oid="Node4316155984")
    nd_4316155984.edge.oid = "Edge4316156048"
    nd_4313802512 = nd_4313829328.new_child(label="Node4313802512", taxon=None, edge_length=22.5856079922, oid="Node4313802512")
    nd_4313802512.edge.oid = "Edge4316156304"
    nd_4316156112 = nd_4316155984.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=15.1414324532, oid="Node4316156112")
    nd_4316156112.edge.oid = "Edge4316156176"
    nd_4316156368 = nd_4316155984.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=15.1414324532, oid="Node4316156368")
    nd_4316156368.edge.oid = "Edge4316156432"
    nd_4316156496 = nd_4313802512.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=16.4481972995, oid="Node4316156496")
    nd_4316156496.edge.oid = "Edge4316156560"
    nd_4313802448 = nd_4313802512.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=16.4481972995, oid="Node4313802448")
    nd_4313802448.edge.oid = "Edge4316156688"
    nd_4316156880 = nd_4316156752.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=31.8813718333, oid="Node4316156880")
    nd_4316156880.edge.oid = "Edge4316156944"
    nd_4316156624 = nd_4316156752.new_child(label="Node4316156624", taxon=None, edge_length=5.97313984611, oid="Node4316156624")
    nd_4316156624.edge.oid = "Edge4316157072"
    nd_4316157136 = nd_4316156624.new_child(label="Node4316157136", taxon=None, edge_length=8.94343133576, oid="Node4316157136")
    nd_4316157136.edge.oid = "Edge4316157200"
    nd_4316157904 = nd_4316156624.new_child(label="Python regius", taxon=tax_4313760144, edge_length=25.9082319872, oid="Node4316157904")
    nd_4316157904.edge.oid = "Edge4316157968"
    nd_4316157264 = nd_4316157136.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=16.9648006514, oid="Node4316157264")
    nd_4316157264.edge.oid = "Edge4316157328"
    nd_4316157008 = nd_4316157136.new_child(label="Node4316157008", taxon=None, edge_length=4.66373979181, oid="Node4316157008")
    nd_4316157008.edge.oid = "Edge4316157520"
    nd_4316157584 = nd_4316157008.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=12.3010608596, oid="Node4316157584")
    nd_4316157584.edge.oid = "Edge4316157648"
    nd_4316157456 = nd_4316157008.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=12.3010608596, oid="Node4316157456")
    nd_4316157456.edge.oid = "Edge4316157776"
    tree_4316158160 = dendropy.Tree(label="Tree03", taxon_set=tree_list.taxon_set, oid="Tree4316158160")
    tree_list.append(tree_4316158160, reindex_taxa=False)
    tree_4316158160.seed_node.oid = 'Node4316158288'
    nd_4316158416 = tree_4316158160.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=109.372663833, oid="Node4316158416")
    nd_4316158416.edge.oid = "Edge4316158480"
    nd_4316158608 = tree_4316158160.seed_node.new_child(label="Node4316158608", taxon=None, edge_length=23.0231215792, oid="Node4316158608")
    nd_4316158608.edge.oid = "Edge4316158672"
    nd_4316158736 = nd_4316158608.new_child(label="Node4316158736", taxon=None, edge_length=25.5450384687, oid="Node4316158736")
    nd_4316158736.edge.oid = "Edge4316158800"
    nd_4316203024 = nd_4316158608.new_child(label="Node4316203024", taxon=None, edge_length=3.40243621372, oid="Node4316203024")
    nd_4316203024.edge.oid = "Edge4316203216"
    nd_4316158864 = nd_4316158736.new_child(label="Node4316158864", taxon=None, edge_length=8.2214192007, oid="Node4316158864")
    nd_4316158864.edge.oid = "Edge4316158928"
    nd_4316202128 = nd_4316158736.new_child(label="Node4316202128", taxon=None, edge_length=17.448321597, oid="Node4316202128")
    nd_4316202128.edge.oid = "Edge4316202320"
    nd_4316158992 = nd_4316158864.new_child(label="Node4316158992", taxon=None, edge_length=6.47536868825, oid="Node4316158992")
    nd_4316158992.edge.oid = "Edge4316159056"
    nd_4316201552 = nd_4316158864.new_child(label="Node4316201552", taxon=None, edge_length=22.2791730332, oid="Node4316201552")
    nd_4316201552.edge.oid = "Edge4316201936"
    nd_4316159120 = nd_4316158992.new_child(label="Node4316159120", taxon=None, edge_length=2.53965409687, oid="Node4316159120")
    nd_4316159120.edge.oid = "Edge4316159184"
    nd_4316179792 = nd_4316158992.new_child(label="Node4316179792", taxon=None, edge_length=7.82653683343, oid="Node4316179792")
    nd_4316179792.edge.oid = "Edge4316179856"
    nd_4316159248 = nd_4316159120.new_child(label="Node4316159248", taxon=None, edge_length=2.41802497137, oid="Node4316159248")
    nd_4316159248.edge.oid = "Edge4316159312"
    nd_4316178128 = nd_4316159120.new_child(label="Node4316178128", taxon=None, edge_length=5.90175129715, oid="Node4316178128")
    nd_4316178128.edge.oid = "Edge4316178448"
    nd_4316159376 = nd_4316159248.new_child(label="Node4316159376", taxon=None, edge_length=6.39175712039, oid="Node4316159376")
    nd_4316159376.edge.oid = "Edge4316159440"
    nd_4316177104 = nd_4316159248.new_child(label="Node4316177104", taxon=None, edge_length=9.11329596086, oid="Node4316177104")
    nd_4316177104.edge.oid = "Edge4316177296"
    nd_4316159504 = nd_4316159376.new_child(label="Node4316159504", taxon=None, edge_length=4.82772953939, oid="Node4316159504")
    nd_4316159504.edge.oid = "Edge4316159568"
    nd_4316176720 = nd_4316159376.new_child(label="Node4316176720", taxon=None, edge_length=17.7955396972, oid="Node4316176720")
    nd_4316176720.edge.oid = "Edge4316176912"
    nd_4316159632 = nd_4316159504.new_child(label="Node4316159632", taxon=None, edge_length=18.9362130146, oid="Node4316159632")
    nd_4316159632.edge.oid = "Edge4316159696"
    nd_4316176464 = nd_4316159504.new_child(label="Node4316176464", taxon=None, edge_length=8.55401069191, oid="Node4316176464")
    nd_4316176464.edge.oid = "Edge4316176528"
    nd_4316159760 = nd_4316159632.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=10.9943371535, oid="Node4316159760")
    nd_4316159760.edge.oid = "Edge4316159824"
    nd_4316158544 = nd_4316159632.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=10.9943371535, oid="Node4316158544")
    nd_4316158544.edge.oid = "Edge4316159952"
    nd_4316176592 = nd_4316176464.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=21.3765394762, oid="Node4316176592")
    nd_4316176592.edge.oid = "Edge4316176656"
    nd_4316159888 = nd_4316176464.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=21.3765394762, oid="Node4316159888")
    nd_4316159888.edge.oid = "Edge4316176784"
    nd_4316176976 = nd_4316176720.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=16.9627400103, oid="Node4316176976")
    nd_4316176976.edge.oid = "Edge4316177040"
    nd_4316176848 = nd_4316176720.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=16.9627400103, oid="Node4316176848")
    nd_4316176848.edge.oid = "Edge4316177168"
    nd_4316177360 = nd_4316177104.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=32.036740867, oid="Node4316177360")
    nd_4316177360.edge.oid = "Edge4316177424"
    nd_4316177232 = nd_4316177104.new_child(label="Node4316177232", taxon=None, edge_length=14.7791918926, oid="Node4316177232")
    nd_4316177232.edge.oid = "Edge4316177552"
    nd_4316177616 = nd_4316177232.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=17.2575489744, oid="Node4316177616")
    nd_4316177616.edge.oid = "Edge4316177680"
    nd_4316177488 = nd_4316177232.new_child(label="Node4316177488", taxon=None, edge_length=13.3651095585, oid="Node4316177488")
    nd_4316177488.edge.oid = "Edge4316177808"
    nd_4316177872 = nd_4316177488.new_child(label="Node4316177872", taxon=None, edge_length=0.439451186875, oid="Node4316177872")
    nd_4316177872.edge.oid = "Edge4316177936"
    nd_4316178256 = nd_4316177488.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=3.89243941597, oid="Node4316178256")
    nd_4316178256.edge.oid = "Edge4316178320"
    nd_4316178000 = nd_4316177872.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=3.4529882291, oid="Node4316178000")
    nd_4316178000.edge.oid = "Edge4316178064"
    nd_4316177744 = nd_4316177872.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=3.4529882291, oid="Node4316177744")
    nd_4316177744.edge.oid = "Edge4316178192"
    nd_4316178512 = nd_4316178128.new_child(label="Node4316178512", taxon=None, edge_length=11.5242004723, oid="Node4316178512")
    nd_4316178512.edge.oid = "Edge4316178576"
    nd_4316179280 = nd_4316178128.new_child(label="Node4316179280", taxon=None, edge_length=13.3915851761, oid="Node4316179280")
    nd_4316179280.edge.oid = "Edge4316179472"
    nd_4316178640 = nd_4316178512.new_child(label="Node4316178640", taxon=None, edge_length=1.13984981318, oid="Node4316178640")
    nd_4316178640.edge.oid = "Edge4316178704"
    nd_4316179024 = nd_4316178512.new_child(label="Node4316179024", taxon=None, edge_length=13.1543046788, oid="Node4316179024")
    nd_4316179024.edge.oid = "Edge4316179088"
    nd_4316178768 = nd_4316178640.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=25.0022602166, oid="Node4316178768")
    nd_4316178768.edge.oid = "Edge4316178832"
    nd_4316178384 = nd_4316178640.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=25.0022602166, oid="Node4316178384")
    nd_4316178384.edge.oid = "Edge4316178960"
    nd_4316179152 = nd_4316179024.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=12.987805351, oid="Node4316179152")
    nd_4316179152.edge.oid = "Edge4316179216"
    nd_4316178896 = nd_4316179024.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=12.987805351, oid="Node4316178896")
    nd_4316178896.edge.oid = "Edge4316179344"
    nd_4316179536 = nd_4316179280.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=24.274725326, oid="Node4316179536")
    nd_4316179536.edge.oid = "Edge4316179600"
    nd_4316179408 = nd_4316179280.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=24.274725326, oid="Node4316179408")
    nd_4316179408.edge.oid = "Edge4316179728"
    nd_4316179920 = nd_4316179792.new_child(label="Node4316179920", taxon=None, edge_length=7.8146690322, oid="Node4316179920")
    nd_4316179920.edge.oid = "Edge4316179984"
    nd_4316180432 = nd_4316179792.new_child(label="Node4316180432", taxon=None, edge_length=5.10842077756, oid="Node4316180432")
    nd_4316180432.edge.oid = "Edge4316158032"
    nd_4316180048 = nd_4316179920.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=30.4665100305, oid="Node4316180048")
    nd_4316180048.edge.oid = "Edge4316180112"
    nd_4316179664 = nd_4316179920.new_child(label="Node4316179664", taxon=None, edge_length=11.0043198537, oid="Node4316179664")
    nd_4316179664.edge.oid = "Edge4316180240"
    nd_4316180304 = nd_4316179664.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=19.4621901768, oid="Node4316180304")
    nd_4316180304.edge.oid = "Edge4316180368"
    nd_4316180176 = nd_4316179664.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=19.4621901768, oid="Node4316180176")
    nd_4316180176.edge.oid = "Edge4316201040"
    nd_4316201104 = nd_4316180432.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=33.1727582851, oid="Node4316201104")
    nd_4316201104.edge.oid = "Edge4316201168"
    nd_4316158224 = nd_4316180432.new_child(label="Node4316158224", taxon=None, edge_length=4.7141022378, oid="Node4316158224")
    nd_4316158224.edge.oid = "Edge4316201232"
    nd_4316201296 = nd_4316158224.new_child(label="Node4316201296", taxon=None, edge_length=19.4308450954, oid="Node4316201296")
    nd_4316201296.edge.oid = "Edge4316201360"
    nd_4316201744 = nd_4316158224.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=28.4586560473, oid="Node4316201744")
    nd_4316201744.edge.oid = "Edge4316201808"
    nd_4316201424 = nd_4316201296.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=9.02781095195, oid="Node4316201424")
    nd_4316201424.edge.oid = "Edge4316201488"
    nd_4316201616 = nd_4316201296.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=9.02781095195, oid="Node4316201616")
    nd_4316201616.edge.oid = "Edge4316201680"
    nd_4316202000 = nd_4316201552.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=30.3039115511, oid="Node4316202000")
    nd_4316202000.edge.oid = "Edge4316202064"
    nd_4316201872 = nd_4316201552.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=30.3039115511, oid="Node4316201872")
    nd_4316201872.edge.oid = "Edge4316202192"
    nd_4316202384 = nd_4316202128.new_child(label="Node4316202384", taxon=None, edge_length=6.14068810478, oid="Node4316202384")
    nd_4316202384.edge.oid = "Edge4316202448"
    nd_4316202896 = nd_4316202128.new_child(label="Python regius", taxon=tax_4313760144, edge_length=43.356182188, oid="Node4316202896")
    nd_4316202896.edge.oid = "Edge4316203088"
    nd_4316202512 = nd_4316202384.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=37.2154940833, oid="Node4316202512")
    nd_4316202512.edge.oid = "Edge4316202576"
    nd_4316202256 = nd_4316202384.new_child(label="Node4316202256", taxon=None, edge_length=15.2807163378, oid="Node4316202256")
    nd_4316202256.edge.oid = "Edge4316202704"
    nd_4316202768 = nd_4316202256.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=21.9347777454, oid="Node4316202768")
    nd_4316202768.edge.oid = "Edge4316202832"
    nd_4316202640 = nd_4316202256.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=21.9347777454, oid="Node4316202640")
    nd_4316202640.edge.oid = "Edge4316202960"
    nd_4316203280 = nd_4316203024.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=82.94710604, oid="Node4316203280")
    nd_4316203280.edge.oid = "Edge4316203344"
    nd_4316203152 = nd_4316203024.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=82.94710604, oid="Node4316203152")
    nd_4316203152.edge.oid = "Edge4316203472"
    tree_4316203536 = dendropy.Tree(label="Tree04", taxon_set=tree_list.taxon_set, oid="Tree4316203536")
    tree_list.append(tree_4316203536, reindex_taxa=False)
    tree_4316203536.seed_node.oid = 'Node4316203664'
    nd_4316203792 = tree_4316203536.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=157.750076773, oid="Node4316203792")
    nd_4316203792.edge.oid = "Edge4316203856"
    nd_4316203984 = tree_4316203536.seed_node.new_child(label="Node4316203984", taxon=None, edge_length=44.9789688242, oid="Node4316203984")
    nd_4316203984.edge.oid = "Edge4316204048"
    nd_4316204112 = nd_4316203984.new_child(label="Node4316204112", taxon=None, edge_length=19.5811677101, oid="Node4316204112")
    nd_4316204112.edge.oid = "Edge4316204176"
    nd_4316252880 = nd_4316203984.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=112.771107949, oid="Node4316252880")
    nd_4316252880.edge.oid = "Edge4316252944"
    nd_4316204240 = nd_4316204112.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=93.189940239, oid="Node4316204240")
    nd_4316204240.edge.oid = "Edge4316204304"
    nd_4316203920 = nd_4316204112.new_child(label="Node4316203920", taxon=None, edge_length=27.291533515, oid="Node4316203920")
    nd_4316203920.edge.oid = "Edge4316204432"
    nd_4316204496 = nd_4316203920.new_child(label="Node4316204496", taxon=None, edge_length=10.3875398007, oid="Node4316204496")
    nd_4316204496.edge.oid = "Edge4316204560"
    nd_4316251856 = nd_4316203920.new_child(label="Node4316251856", taxon=None, edge_length=19.8895314814, oid="Node4316251856")
    nd_4316251856.edge.oid = "Edge4316252048"
    nd_4316204624 = nd_4316204496.new_child(label="Node4316204624", taxon=None, edge_length=9.74044751021, oid="Node4316204624")
    nd_4316204624.edge.oid = "Edge4316204688"
    nd_4316251600 = nd_4316204496.new_child(label="Node4316251600", taxon=None, edge_length=22.1202819118, oid="Node4316251600")
    nd_4316251600.edge.oid = "Edge4316251664"
    nd_4316204752 = nd_4316204624.new_child(label="Node4316204752", taxon=None, edge_length=2.88538115771, oid="Node4316204752")
    nd_4316204752.edge.oid = "Edge4316204816"
    nd_4316227472 = nd_4316204624.new_child(label="Node4316227472", taxon=None, edge_length=2.72676824396, oid="Node4316227472")
    nd_4316227472.edge.oid = "Edge4316227536"
    nd_4316204880 = nd_4316204752.new_child(label="Node4316204880", taxon=None, edge_length=1.34047834167, oid="Node4316204880")
    nd_4316204880.edge.oid = "Edge4316204944"
    nd_4316227216 = nd_4316204752.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=42.8850382554, oid="Node4316227216")
    nd_4316227216.edge.oid = "Edge4316227408"
    nd_4316205008 = nd_4316204880.new_child(label="Node4316205008", taxon=None, edge_length=2.31767871982, oid="Node4316205008")
    nd_4316205008.edge.oid = "Edge4316225616"
    nd_4316226768 = nd_4316204880.new_child(label="Node4316226768", taxon=None, edge_length=18.2547475502, oid="Node4316226768")
    nd_4316226768.edge.oid = "Edge4316227024"
    nd_4316225680 = nd_4316205008.new_child(label="Node4316225680", taxon=None, edge_length=6.3930928479, oid="Node4316225680")
    nd_4316225680.edge.oid = "Edge4316225744"
    nd_4316226448 = nd_4316205008.new_child(label="Node4316226448", taxon=None, edge_length=23.2404397828, oid="Node4316226448")
    nd_4316226448.edge.oid = "Edge4316226576"
    nd_4316225808 = nd_4316225680.new_child(label="Node4316225808", taxon=None, edge_length=24.6792964706, oid="Node4316225808")
    nd_4316225808.edge.oid = "Edge4316225872"
    nd_4316226192 = nd_4316225680.new_child(label="Node4316226192", taxon=None, edge_length=8.52936801714, oid="Node4316226192")
    nd_4316226192.edge.oid = "Edge4316226256"
    nd_4316225936 = nd_4316225808.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=8.15449187544, oid="Node4316225936")
    nd_4316225936.edge.oid = "Edge4316226000"
    nd_4316204368 = nd_4316225808.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=8.15449187544, oid="Node4316204368")
    nd_4316204368.edge.oid = "Edge4316226064"
    nd_4316226320 = nd_4316226192.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=24.3044203289, oid="Node4316226320")
    nd_4316226320.edge.oid = "Edge4316226384"
    nd_4316226128 = nd_4316226192.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=24.3044203289, oid="Node4316226128")
    nd_4316226128.edge.oid = "Edge4316226512"
    nd_4316226640 = nd_4316226448.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=15.9864414111, oid="Node4316226640")
    nd_4316226640.edge.oid = "Edge4316226704"
    nd_4316226832 = nd_4316226448.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=15.9864414111, oid="Node4316226832")
    nd_4316226832.edge.oid = "Edge4316226896"
    nd_4316227088 = nd_4316226768.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=23.2898123636, oid="Node4316227088")
    nd_4316227088.edge.oid = "Edge4316227152"
    nd_4316226960 = nd_4316226768.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=23.2898123636, oid="Node4316226960")
    nd_4316226960.edge.oid = "Edge4316227280"
    nd_4316227600 = nd_4316227472.new_child(label="Node4316227600", taxon=None, edge_length=14.2175774566, oid="Node4316227600")
    nd_4316227600.edge.oid = "Edge4316227664"
    nd_4316229392 = nd_4316227472.new_child(label="Node4316229392", taxon=None, edge_length=3.9347409374, oid="Node4316229392")
    nd_4316229392.edge.oid = "Edge4316229456"
    nd_4316227728 = nd_4316227600.new_child(label="Node4316227728", taxon=None, edge_length=12.5474231006, oid="Node4316227728")
    nd_4316227728.edge.oid = "Edge4316227792"
    nd_4316228112 = nd_4316227600.new_child(label="Node4316228112", taxon=None, edge_length=1.26678175478, oid="Node4316228112")
    nd_4316228112.edge.oid = "Edge4316228176"
    nd_4316227856 = nd_4316227728.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=16.278650612, oid="Node4316227856")
    nd_4316227856.edge.oid = "Edge4316227920"
    nd_4316227344 = nd_4316227728.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=16.278650612, oid="Node4316227344")
    nd_4316227344.edge.oid = "Edge4316228048"
    nd_4316228240 = nd_4316228112.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=27.5592919578, oid="Node4316228240")
    nd_4316228240.edge.oid = "Edge4316228304"
    nd_4316227984 = nd_4316228112.new_child(label="Node4316227984", taxon=None, edge_length=13.3039152583, oid="Node4316227984")
    nd_4316227984.edge.oid = "Edge4316228432"
    nd_4316228496 = nd_4316227984.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=14.2553766995, oid="Node4316228496")
    nd_4316228496.edge.oid = "Edge4316228560"
    nd_4316228368 = nd_4316227984.new_child(label="Node4316228368", taxon=None, edge_length=1.66170791236, oid="Node4316228368")
    nd_4316228368.edge.oid = "Edge4316228688"
    nd_4316228752 = nd_4316228368.new_child(label="Node4316228752", taxon=None, edge_length=8.89489387836, oid="Node4316228752")
    nd_4316228752.edge.oid = "Edge4316228816"
    nd_4316229264 = nd_4316228368.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=12.5936687872, oid="Node4316229264")
    nd_4316229264.edge.oid = "Edge4316229328"
    nd_4316228880 = nd_4316228752.new_child(label="Node4316228880", taxon=None, edge_length=0.230110019205, oid="Node4316228880")
    nd_4316228880.edge.oid = "Edge4316228944"
    nd_4316229136 = nd_4316228752.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=3.69877490882, oid="Node4316229136")
    nd_4316229136.edge.oid = "Edge4316229200"
    nd_4316229008 = nd_4316228880.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=3.46866488962, oid="Node4316229008")
    nd_4316229008.edge.oid = "Edge4316229072"
    nd_4316203408 = nd_4316228880.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=3.46866488962, oid="Node4316203408")
    nd_4316203408.edge.oid = "Edge4316228624"
    nd_4316229520 = nd_4316229392.new_child(label="Node4316229520", taxon=None, edge_length=5.88218975316, oid="Node4316229520")
    nd_4316229520.edge.oid = "Edge4316229584"
    nd_4316250576 = nd_4316229392.new_child(label="Node4316250576", taxon=None, edge_length=7.11128547149, oid="Node4316250576")
    nd_4316250576.edge.oid = "Edge4316250768"
    nd_4316250192 = nd_4316229520.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=33.2267204786, oid="Node4316250192")
    nd_4316250192.edge.oid = "Edge4316250256"
    nd_4316203600 = nd_4316229520.new_child(label="Node4316203600", taxon=None, edge_length=13.4306199458, oid="Node4316203600")
    nd_4316203600.edge.oid = "Edge4316250384"
    nd_4316250448 = nd_4316203600.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=19.7961005329, oid="Node4316250448")
    nd_4316250448.edge.oid = "Edge4316250512"
    nd_4316250320 = nd_4316203600.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=19.7961005329, oid="Node4316250320")
    nd_4316250320.edge.oid = "Edge4316250640"
    nd_4316250832 = nd_4316250576.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=31.9976247603, oid="Node4316250832")
    nd_4316250832.edge.oid = "Edge4316250896"
    nd_4316250704 = nd_4316250576.new_child(label="Node4316250704", taxon=None, edge_length=4.71436528425, oid="Node4316250704")
    nd_4316250704.edge.oid = "Edge4316251024"
    nd_4316251088 = nd_4316250704.new_child(label="Node4316251088", taxon=None, edge_length=15.9285543528, oid="Node4316251088")
    nd_4316251088.edge.oid = "Edge4316251152"
    nd_4316251472 = nd_4316250704.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=27.2832594761, oid="Node4316251472")
    nd_4316251472.edge.oid = "Edge4316251536"
    nd_4316251216 = nd_4316251088.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=11.3547051232, oid="Node4316251216")
    nd_4316251216.edge.oid = "Edge4316251280"
    nd_4316250960 = nd_4316251088.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=11.3547051232, oid="Node4316250960")
    nd_4316250960.edge.oid = "Edge4316251408"
    nd_4316251728 = nd_4316251600.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=33.3905850116, oid="Node4316251728")
    nd_4316251728.edge.oid = "Edge4316251792"
    nd_4316251344 = nd_4316251600.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=33.3905850116, oid="Node4316251344")
    nd_4316251344.edge.oid = "Edge4316251920"
    nd_4316252112 = nd_4316251856.new_child(label="Node4316252112", taxon=None, edge_length=11.1907198563, oid="Node4316252112")
    nd_4316252112.edge.oid = "Edge4316252176"
    nd_4316252624 = nd_4316251856.new_child(label="Python regius", taxon=tax_4313760144, edge_length=46.0088752426, oid="Node4316252624")
    nd_4316252624.edge.oid = "Edge4316252816"
    nd_4316252240 = nd_4316252112.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=34.8181553863, oid="Node4316252240")
    nd_4316252240.edge.oid = "Edge4316252304"
    nd_4316251984 = nd_4316252112.new_child(label="Node4316251984", taxon=None, edge_length=7.89583224277, oid="Node4316251984")
    nd_4316251984.edge.oid = "Edge4316252432"
    nd_4316252496 = nd_4316251984.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=26.9223231435, oid="Node4316252496")
    nd_4316252496.edge.oid = "Edge4316252560"
    nd_4316252368 = nd_4316251984.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=26.9223231435, oid="Node4316252368")
    nd_4316252368.edge.oid = "Edge4316252688"
    tree_4316252752 = dendropy.Tree(label="Tree05", taxon_set=tree_list.taxon_set, oid="Tree4316252752")
    tree_list.append(tree_4316252752, reindex_taxa=False)
    tree_4316252752.seed_node.oid = 'Node4316253136'
    nd_4316253264 = tree_4316252752.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=117.405482731, oid="Node4316253264")
    nd_4316253264.edge.oid = "Edge4316253328"
    nd_4316253456 = tree_4316252752.seed_node.new_child(label="Node4316253456", taxon=None, edge_length=26.8517813278, oid="Node4316253456")
    nd_4316253456.edge.oid = "Edge4316253520"
    nd_4316253584 = nd_4316253456.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=90.5537014032, oid="Node4316253584")
    nd_4316253584.edge.oid = "Edge4316253648"
    nd_4316253392 = nd_4316253456.new_child(label="Node4316253392", taxon=None, edge_length=16.5508683851, oid="Node4316253392")
    nd_4316253392.edge.oid = "Edge4316253712"
    nd_4316253776 = nd_4316253392.new_child(label="Node4316253776", taxon=None, edge_length=16.6353899137, oid="Node4316253776")
    nd_4316253776.edge.oid = "Edge4316253840"
    nd_4316302224 = nd_4316253392.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=74.0028330181, oid="Node4316302224")
    nd_4316302224.edge.oid = "Edge4316302416"
    nd_4316253904 = nd_4316253776.new_child(label="Node4316253904", taxon=None, edge_length=6.56321044476, oid="Node4316253904")
    nd_4316253904.edge.oid = "Edge4316253968"
    nd_4316301328 = nd_4316253776.new_child(label="Node4316301328", taxon=None, edge_length=18.1852453647, oid="Node4316301328")
    nd_4316301328.edge.oid = "Edge4316301520"
    nd_4316254032 = nd_4316253904.new_child(label="Node4316254032", taxon=None, edge_length=12.3993667148, oid="Node4316254032")
    nd_4316254032.edge.oid = "Edge4316254096"
    nd_4316301072 = nd_4316253904.new_child(label="Node4316301072", taxon=None, edge_length=23.3069353709, oid="Node4316301072")
    nd_4316301072.edge.oid = "Edge4316301136"
    nd_4316254160 = nd_4316254032.new_child(label="Node4316254160", taxon=None, edge_length=1.85501747057, oid="Node4316254160")
    nd_4316254160.edge.oid = "Edge4316274768"
    nd_4316276944 = nd_4316254032.new_child(label="Node4316276944", taxon=None, edge_length=2.5597754437, oid="Node4316276944")
    nd_4316276944.edge.oid = "Edge4316277008"
    nd_4316274832 = nd_4316254160.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=36.5498484742, oid="Node4316274832")
    nd_4316274832.edge.oid = "Edge4316274896"
    nd_4316275024 = nd_4316254160.new_child(label="Node4316275024", taxon=None, edge_length=1.04661210215, oid="Node4316275024")
    nd_4316275024.edge.oid = "Edge4316275088"
    nd_4316275152 = nd_4316275024.new_child(label="Node4316275152", taxon=None, edge_length=6.27408700912, oid="Node4316275152")
    nd_4316275152.edge.oid = "Edge4316275216"
    nd_4316276432 = nd_4316275024.new_child(label="Node4316276432", taxon=None, edge_length=15.4075337774, oid="Node4316276432")
    nd_4316276432.edge.oid = "Edge4316276624"
    nd_4316275280 = nd_4316275152.new_child(label="Node4316275280", taxon=None, edge_length=4.14502379891, oid="Node4316275280")
    nd_4316275280.edge.oid = "Edge4316275344"
    nd_4316276048 = nd_4316275152.new_child(label="Node4316276048", taxon=None, edge_length=17.4625828861, oid="Node4316276048")
    nd_4316276048.edge.oid = "Edge4316276240"
    nd_4316275408 = nd_4316275280.new_child(label="Node4316275408", taxon=None, edge_length=18.7911640437, oid="Node4316275408")
    nd_4316275408.edge.oid = "Edge4316275472"
    nd_4316275792 = nd_4316275280.new_child(label="Node4316275792", taxon=None, edge_length=4.65165886356, oid="Node4316275792")
    nd_4316275792.edge.oid = "Edge4316275856"
    nd_4316275536 = nd_4316275408.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=6.29296152027, oid="Node4316275536")
    nd_4316275536.edge.oid = "Edge4316275600"
    nd_4316274960 = nd_4316275408.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=6.29296152027, oid="Node4316274960")
    nd_4316274960.edge.oid = "Edge4316275728"
    nd_4316275920 = nd_4316275792.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=20.4324667004, oid="Node4316275920")
    nd_4316275920.edge.oid = "Edge4316275984"
    nd_4316275664 = nd_4316275792.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=20.4324667004, oid="Node4316275664")
    nd_4316275664.edge.oid = "Edge4316276112"
    nd_4316276304 = nd_4316276048.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=11.7665664768, oid="Node4316276304")
    nd_4316276304.edge.oid = "Edge4316276368"
    nd_4316276176 = nd_4316276048.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=11.7665664768, oid="Node4316276176")
    nd_4316276176.edge.oid = "Edge4316276496"
    nd_4316276688 = nd_4316276432.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=20.0957025946, oid="Node4316276688")
    nd_4316276688.edge.oid = "Edge4316276752"
    nd_4316276560 = nd_4316276432.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=20.0957025946, oid="Node4316276560")
    nd_4316276560.edge.oid = "Edge4316276880"
    nd_4316277072 = nd_4316276944.new_child(label="Node4316277072", taxon=None, edge_length=10.9999278199, oid="Node4316277072")
    nd_4316277072.edge.oid = "Edge4316277136"
    nd_4316299408 = nd_4316276944.new_child(label="Node4316299408", taxon=None, edge_length=4.77611939702, oid="Node4316299408")
    nd_4316299408.edge.oid = "Edge4316299472"
    nd_4316277200 = nd_4316277072.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=24.8451626811, oid="Node4316277200")
    nd_4316277200.edge.oid = "Edge4316277264"
    nd_4316276816 = nd_4316277072.new_child(label="Node4316276816", taxon=None, edge_length=2.78884378851, oid="Node4316276816")
    nd_4316276816.edge.oid = "Edge4316277328"
    nd_4316277392 = nd_4316276816.new_child(label="Node4316277392", taxon=None, edge_length=7.72979171285, oid="Node4316277392")
    nd_4316277392.edge.oid = "Edge4316277456"
    nd_4316278352 = nd_4316276816.new_child(label="Node4316278352", taxon=None, edge_length=12.4767981898, oid="Node4316278352")
    nd_4316278352.edge.oid = "Edge4316278544"
    nd_4316277520 = nd_4316277392.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=14.3265271797, oid="Node4316277520")
    nd_4316277520.edge.oid = "Edge4316253072"
    nd_4316277584 = nd_4316277392.new_child(label="Node4316277584", taxon=None, edge_length=3.38770784397, oid="Node4316277584")
    nd_4316277584.edge.oid = "Edge4316277648"
    nd_4316277712 = nd_4316277584.new_child(label="Node4316277712", taxon=None, edge_length=5.65307392619, oid="Node4316277712")
    nd_4316277712.edge.oid = "Edge4316277776"
    nd_4316278224 = nd_4316277584.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=10.9388193358, oid="Node4316278224")
    nd_4316278224.edge.oid = "Edge4316278416"
    nd_4316277840 = nd_4316277712.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=5.28574540958, oid="Node4316277840")
    nd_4316277840.edge.oid = "Edge4316277904"
    nd_4316253008 = nd_4316277712.new_child(label="Node4316253008", taxon=None, edge_length=2.55413552294, oid="Node4316253008")
    nd_4316253008.edge.oid = "Edge4316278032"
    nd_4316278096 = nd_4316253008.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.73160988664, oid="Node4316278096")
    nd_4316278096.edge.oid = "Edge4316278160"
    nd_4316277968 = nd_4316253008.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=2.73160988664, oid="Node4316277968")
    nd_4316277968.edge.oid = "Edge4316278288"
    nd_4316278608 = nd_4316278352.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=9.57952070277, oid="Node4316278608")
    nd_4316278608.edge.oid = "Edge4316278672"
    nd_4316278480 = nd_4316278352.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=9.57952070277, oid="Node4316278480")
    nd_4316278480.edge.oid = "Edge4316278736"
    nd_4316299536 = nd_4316299408.new_child(label="Node4316299536", taxon=None, edge_length=6.62946104565, oid="Node4316299536")
    nd_4316299536.edge.oid = "Edge4316299600"
    nd_4316300048 = nd_4316299408.new_child(label="Node4316300048", taxon=None, edge_length=3.07031045323, oid="Node4316300048")
    nd_4316300048.edge.oid = "Edge4316300240"
    nd_4316299664 = nd_4316299536.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=24.4395100584, oid="Node4316299664")
    nd_4316299664.edge.oid = "Edge4316299728"
    nd_4316299344 = nd_4316299536.new_child(label="Node4316299344", taxon=None, edge_length=7.74158436759, oid="Node4316299344")
    nd_4316299344.edge.oid = "Edge4316299856"
    nd_4316299920 = nd_4316299344.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=16.6979256908, oid="Node4316299920")
    nd_4316299920.edge.oid = "Edge4316299984"
    nd_4316299792 = nd_4316299344.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=16.6979256908, oid="Node4316299792")
    nd_4316299792.edge.oid = "Edge4316300112"
    nd_4316300304 = nd_4316300048.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=27.9986606508, oid="Node4316300304")
    nd_4316300304.edge.oid = "Edge4316300368"
    nd_4316300176 = nd_4316300048.new_child(label="Node4316300176", taxon=None, edge_length=4.88214307652, oid="Node4316300176")
    nd_4316300176.edge.oid = "Edge4316300496"
    nd_4316300560 = nd_4316300176.new_child(label="Node4316300560", taxon=None, edge_length=16.2064963991, oid="Node4316300560")
    nd_4316300560.edge.oid = "Edge4316300624"
    nd_4316300944 = nd_4316300176.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=23.1165175743, oid="Node4316300944")
    nd_4316300944.edge.oid = "Edge4316301008"
    nd_4316300688 = nd_4316300560.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=6.91002117518, oid="Node4316300688")
    nd_4316300688.edge.oid = "Edge4316300752"
    nd_4316300432 = nd_4316300560.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=6.91002117518, oid="Node4316300432")
    nd_4316300432.edge.oid = "Edge4316300880"
    nd_4316301200 = nd_4316301072.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=27.4972972887, oid="Node4316301200")
    nd_4316301200.edge.oid = "Edge4316301264"
    nd_4316300816 = nd_4316301072.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=27.4972972887, oid="Node4316300816")
    nd_4316300816.edge.oid = "Edge4316301392"
    nd_4316301584 = nd_4316301328.new_child(label="Node4316301584", taxon=None, edge_length=6.08258726043, oid="Node4316301584")
    nd_4316301584.edge.oid = "Edge4316301648"
    nd_4316302096 = nd_4316301328.new_child(label="Python regius", taxon=tax_4313760144, edge_length=39.1821977396, oid="Node4316302096")
    nd_4316302096.edge.oid = "Edge4316302288"
    nd_4316301712 = nd_4316301584.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=33.0996104792, oid="Node4316301712")
    nd_4316301712.edge.oid = "Edge4316301776"
    nd_4316301456 = nd_4316301584.new_child(label="Node4316301456", taxon=None, edge_length=14.293345062, oid="Node4316301456")
    nd_4316301456.edge.oid = "Edge4316301904"
    nd_4316301968 = nd_4316301456.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=18.8062654172, oid="Node4316301968")
    nd_4316301968.edge.oid = "Edge4316302032"
    nd_4316301840 = nd_4316301456.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=18.8062654172, oid="Node4316301840")
    nd_4316301840.edge.oid = "Edge4316302160"
    tree_4316302352 = dendropy.Tree(label="Tree06", taxon_set=tree_list.taxon_set, oid="Tree4316302352")
    tree_list.append(tree_4316302352, reindex_taxa=False)
    tree_4316302352.seed_node.oid = 'Node4316302608'
    nd_4316302736 = tree_4316302352.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=126.159749307, oid="Node4316302736")
    nd_4316302736.edge.oid = "Edge4316302800"
    nd_4316302928 = tree_4316302352.seed_node.new_child(label="Node4316302928", taxon=None, edge_length=29.5182899091, oid="Node4316302928")
    nd_4316302928.edge.oid = "Edge4316302992"
    nd_4316303056 = nd_4316302928.new_child(label="Node4316303056", taxon=None, edge_length=20.4272017262, oid="Node4316303056")
    nd_4316303056.edge.oid = "Edge4316303120"
    nd_4316347728 = nd_4316302928.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=96.6414593981, oid="Node4316347728")
    nd_4316347728.edge.oid = "Edge4316347792"
    nd_4316303184 = nd_4316303056.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=76.2142576718, oid="Node4316303184")
    nd_4316303184.edge.oid = "Edge4316303248"
    nd_4316302864 = nd_4316303056.new_child(label="Node4316302864", taxon=None, edge_length=22.3340694246, oid="Node4316302864")
    nd_4316302864.edge.oid = "Edge4316303312"
    nd_4316323920 = nd_4316302864.new_child(label="Node4316323920", taxon=None, edge_length=7.33856679933, oid="Node4316323920")
    nd_4316323920.edge.oid = "Edge4316323984"
    nd_4316346704 = nd_4316302864.new_child(label="Node4316346704", taxon=None, edge_length=17.9081848887, oid="Node4316346704")
    nd_4316346704.edge.oid = "Edge4316346896"
    nd_4316324048 = nd_4316323920.new_child(label="Node4316324048", taxon=None, edge_length=8.24675869634, oid="Node4316324048")
    nd_4316324048.edge.oid = "Edge4316324112"
    nd_4316346448 = nd_4316323920.new_child(label="Node4316346448", taxon=None, edge_length=23.556392099, oid="Node4316346448")
    nd_4316346448.edge.oid = "Edge4316346512"
    nd_4316324176 = nd_4316324048.new_child(label="Node4316324176", taxon=None, edge_length=1.78294733026, oid="Node4316324176")
    nd_4316324176.edge.oid = "Edge4316324240"
    nd_4316326288 = nd_4316324048.new_child(label="Node4316326288", taxon=None, edge_length=3.26747919544, oid="Node4316326288")
    nd_4316326288.edge.oid = "Edge4316326352"
    nd_4316324304 = nd_4316324176.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=36.5119154213, oid="Node4316324304")
    nd_4316324304.edge.oid = "Edge4316324368"
    nd_4316324496 = nd_4316324176.new_child(label="Node4316324496", taxon=None, edge_length=1.38475127065, oid="Node4316324496")
    nd_4316324496.edge.oid = "Edge4316324560"
    nd_4316324624 = nd_4316324496.new_child(label="Node4316324624", taxon=None, edge_length=2.24740276648, oid="Node4316324624")
    nd_4316324624.edge.oid = "Edge4316324688"
    nd_4316326032 = nd_4316324496.new_child(label="Node4316326032", taxon=None, edge_length=12.9799500845, oid="Node4316326032")
    nd_4316326032.edge.oid = "Edge4316326096"
    nd_4316324752 = nd_4316324624.new_child(label="Node4316324752", taxon=None, edge_length=17.772328432, oid="Node4316324752")
    nd_4316324752.edge.oid = "Edge4316324816"
    nd_4316325136 = nd_4316324624.new_child(label="Node4316325136", taxon=None, edge_length=4.4885389798, oid="Node4316325136")
    nd_4316325136.edge.oid = "Edge4316325200"
    nd_4316324880 = nd_4316324752.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=15.1074329521, oid="Node4316324880")
    nd_4316324880.edge.oid = "Edge4316324944"
    nd_4316324432 = nd_4316324752.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=15.1074329521, oid="Node4316324432")
    nd_4316324432.edge.oid = "Edge4316325072"
    nd_4316325264 = nd_4316325136.new_child(label="Node4316325264", taxon=None, edge_length=20.8200876951, oid="Node4316325264")
    nd_4316325264.edge.oid = "Edge4316325328"
    nd_4316325648 = nd_4316325136.new_child(label="Node4316325648", taxon=None, edge_length=4.37289319177, oid="Node4316325648")
    nd_4316325648.edge.oid = "Edge4316325712"
    nd_4316325392 = nd_4316325264.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=7.5711347092, oid="Node4316325392")
    nd_4316325392.edge.oid = "Edge4316325456"
    nd_4316325008 = nd_4316325264.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=7.5711347092, oid="Node4316325008")
    nd_4316325008.edge.oid = "Edge4316325584"
    nd_4316325776 = nd_4316325648.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=24.0183292126, oid="Node4316325776")
    nd_4316325776.edge.oid = "Edge4316325840"
    nd_4316325520 = nd_4316325648.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=24.0183292126, oid="Node4316325520")
    nd_4316325520.edge.oid = "Edge4316325968"
    nd_4316326160 = nd_4316326032.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=22.1472140662, oid="Node4316326160")
    nd_4316326160.edge.oid = "Edge4316302544"
    nd_4316325904 = nd_4316326032.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=22.1472140662, oid="Node4316325904")
    nd_4316325904.edge.oid = "Edge4316326224"
    nd_4316326416 = nd_4316326288.new_child(label="Node4316326416", taxon=None, edge_length=12.2960670728, oid="Node4316326416")
    nd_4316326416.edge.oid = "Edge4316326480"
    nd_4316344784 = nd_4316326288.new_child(label="Node4316344784", taxon=None, edge_length=4.17834235089, oid="Node4316344784")
    nd_4316344784.edge.oid = "Edge4316344848"
    nd_4316326544 = nd_4316326416.new_child(label="Node4316326544", taxon=None, edge_length=1.62908011033, oid="Node4316326544")
    nd_4316326544.edge.oid = "Edge4316326608"
    nd_4316327184 = nd_4316326416.new_child(label="Node4316327184", taxon=None, edge_length=6.20414166284, oid="Node4316327184")
    nd_4316327184.edge.oid = "Edge4316327248"
    nd_4316326672 = nd_4316326544.new_child(label="Node4316326672", taxon=None, edge_length=11.7610480778, oid="Node4316326672")
    nd_4316326672.edge.oid = "Edge4316326736"
    nd_4316327056 = nd_4316326544.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=21.1022363729, oid="Node4316327056")
    nd_4316327056.edge.oid = "Edge4316327120"
    nd_4316326800 = nd_4316326672.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=9.34118829512, oid="Node4316326800")
    nd_4316326800.edge.oid = "Edge4316326864"
    nd_4316302480 = nd_4316326672.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=9.34118829512, oid="Node4316302480")
    nd_4316302480.edge.oid = "Edge4316326992"
    nd_4316327312 = nd_4316327184.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=16.5271748204, oid="Node4316327312")
    nd_4316327312.edge.oid = "Edge4316327376"
    nd_4316326928 = nd_4316327184.new_child(label="Node4316326928", taxon=None, edge_length=4.27339647744, oid="Node4316326928")
    nd_4316326928.edge.oid = "Edge4316327504"
    nd_4316327568 = nd_4316326928.new_child(label="Node4316327568", taxon=None, edge_length=7.91289243511, oid="Node4316327568")
    nd_4316327568.edge.oid = "Edge4316327632"
    nd_4316344656 = nd_4316326928.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=12.253778343, oid="Node4316344656")
    nd_4316344656.edge.oid = "Edge4316344720"
    nd_4316327696 = nd_4316327568.new_child(label="Node4316327696", taxon=None, edge_length=1.89524872622, oid="Node4316327696")
    nd_4316327696.edge.oid = "Edge4316327760"
    nd_4316344528 = nd_4316327568.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=4.34088590785, oid="Node4316344528")
    nd_4316344528.edge.oid = "Edge4316344592"
    nd_4316327824 = nd_4316327696.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=2.44563718163, oid="Node4316327824")
    nd_4316327824.edge.oid = "Edge4316327888"
    nd_4316327440 = nd_4316327696.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.44563718163, oid="Node4316327440")
    nd_4316327440.edge.oid = "Edge4316344464"
    nd_4316344912 = nd_4316344784.new_child(label="Node4316344912", taxon=None, edge_length=5.03710806168, oid="Node4316344912")
    nd_4316344912.edge.oid = "Edge4316344976"
    nd_4316345424 = nd_4316344784.new_child(label="Node4316345424", taxon=None, edge_length=4.09757269601, oid="Node4316345424")
    nd_4316345424.edge.oid = "Edge4316345616"
    nd_4316345040 = nd_4316344912.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=25.8119331435, oid="Node4316345040")
    nd_4316345040.edge.oid = "Edge4316345104"
    nd_4316344400 = nd_4316344912.new_child(label="Node4316344400", taxon=None, edge_length=11.4337686931, oid="Node4316344400")
    nd_4316344400.edge.oid = "Edge4316345232"
    nd_4316345296 = nd_4316344400.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=14.3781644505, oid="Node4316345296")
    nd_4316345296.edge.oid = "Edge4316345360"
    nd_4316345168 = nd_4316344400.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=14.3781644505, oid="Node4316345168")
    nd_4316345168.edge.oid = "Edge4316345488"
    nd_4316345680 = nd_4316345424.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=26.7514685092, oid="Node4316345680")
    nd_4316345680.edge.oid = "Edge4316345744"
    nd_4316345552 = nd_4316345424.new_child(label="Node4316345552", taxon=None, edge_length=2.33638652753, oid="Node4316345552")
    nd_4316345552.edge.oid = "Edge4316345872"
    nd_4316345936 = nd_4316345552.new_child(label="Node4316345936", taxon=None, edge_length=17.4436719435, oid="Node4316345936")
    nd_4316345936.edge.oid = "Edge4316346000"
    nd_4316346320 = nd_4316345552.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=24.4150819817, oid="Node4316346320")
    nd_4316346320.edge.oid = "Edge4316346384"
    nd_4316346064 = nd_4316345936.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=6.97141003814, oid="Node4316346064")
    nd_4316346064.edge.oid = "Edge4316346128"
    nd_4316345808 = nd_4316345936.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=6.97141003814, oid="Node4316345808")
    nd_4316345808.edge.oid = "Edge4316346256"
    nd_4316346576 = nd_4316346448.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=22.9852293488, oid="Node4316346576")
    nd_4316346576.edge.oid = "Edge4316346640"
    nd_4316346192 = nd_4316346448.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=22.9852293488, oid="Node4316346192")
    nd_4316346192.edge.oid = "Edge4316346768"
    nd_4316346960 = nd_4316346704.new_child(label="Node4316346960", taxon=None, edge_length=6.75664224871, oid="Node4316346960")
    nd_4316346960.edge.oid = "Edge4316347024"
    nd_4316347472 = nd_4316346704.new_child(label="Python regius", taxon=tax_4313760144, edge_length=35.9720033585, oid="Node4316347472")
    nd_4316347472.edge.oid = "Edge4316347664"
    nd_4316347088 = nd_4316346960.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=29.2153611098, oid="Node4316347088")
    nd_4316347088.edge.oid = "Edge4316347152"
    nd_4316346832 = nd_4316346960.new_child(label="Node4316346832", taxon=None, edge_length=8.15978945225, oid="Node4316346832")
    nd_4316346832.edge.oid = "Edge4316347280"
    nd_4316347344 = nd_4316346832.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=21.0555716576, oid="Node4316347344")
    nd_4316347344.edge.oid = "Edge4316347408"
    nd_4316347216 = nd_4316346832.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=21.0555716576, oid="Node4316347216")
    nd_4316347216.edge.oid = "Edge4316347536"
    tree_4316347600 = dendropy.Tree(label="Tree07", taxon_set=tree_list.taxon_set, oid="Tree4316347600")
    tree_list.append(tree_4316347600, reindex_taxa=False)
    tree_4316347600.seed_node.oid = 'Node4316347984'
    nd_4316348112 = tree_4316347600.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=124.564186516, oid="Node4316348112")
    nd_4316348112.edge.oid = "Edge4316348176"
    nd_4316348304 = tree_4316347600.seed_node.new_child(label="Node4316348304", taxon=None, edge_length=36.3676780441, oid="Node4316348304")
    nd_4316348304.edge.oid = "Edge4316348368"
    nd_4316368976 = nd_4316348304.new_child(label="Node4316368976", taxon=None, edge_length=11.1789504571, oid="Node4316368976")
    nd_4316368976.edge.oid = "Edge4316369040"
    nd_4316397200 = nd_4316348304.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=88.196508472, oid="Node4316397200")
    nd_4316397200.edge.oid = "Edge4316397264"
    nd_4316369104 = nd_4316368976.new_child(label="Node4316369104", taxon=None, edge_length=20.3346663059, oid="Node4316369104")
    nd_4316369104.edge.oid = "Edge4316369168"
    nd_4316396944 = nd_4316368976.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=77.0175580149, oid="Node4316396944")
    nd_4316396944.edge.oid = "Edge4316397136"
    nd_4316369232 = nd_4316369104.new_child(label="Node4316369232", taxon=None, edge_length=10.7040090023, oid="Node4316369232")
    nd_4316369232.edge.oid = "Edge4316369296"
    nd_4316396048 = nd_4316369104.new_child(label="Node4316396048", taxon=None, edge_length=19.1477140595, oid="Node4316396048")
    nd_4316396048.edge.oid = "Edge4316396240"
    nd_4316369360 = nd_4316369232.new_child(label="Node4316369360", taxon=None, edge_length=7.28944403695, oid="Node4316369360")
    nd_4316369360.edge.oid = "Edge4316369424"
    nd_4316395664 = nd_4316369232.new_child(label="Node4316395664", taxon=None, edge_length=20.035415736, oid="Node4316395664")
    nd_4316395664.edge.oid = "Edge4316395856"
    nd_4316369488 = nd_4316369360.new_child(label="Node4316369488", taxon=None, edge_length=2.76745367097, oid="Node4316369488")
    nd_4316369488.edge.oid = "Edge4316369552"
    nd_4316372752 = nd_4316369360.new_child(label="Node4316372752", taxon=None, edge_length=1.09335779327, oid="Node4316372752")
    nd_4316372752.edge.oid = "Edge4316393680"
    nd_4316369616 = nd_4316369488.new_child(label="Node4316369616", taxon=None, edge_length=12.5842551962, oid="Node4316369616")
    nd_4316369616.edge.oid = "Edge4316369680"
    nd_4316371280 = nd_4316369488.new_child(label="Node4316371280", taxon=None, edge_length=7.66239308607, oid="Node4316371280")
    nd_4316371280.edge.oid = "Edge4316371472"
    nd_4316369744 = nd_4316369616.new_child(label="Node4316369744", taxon=None, edge_length=12.0191169658, oid="Node4316369744")
    nd_4316369744.edge.oid = "Edge4316369808"
    nd_4316370128 = nd_4316369616.new_child(label="Node4316370128", taxon=None, edge_length=2.25737379321, oid="Node4316370128")
    nd_4316370128.edge.oid = "Edge4316370192"
    nd_4316369872 = nd_4316369744.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=11.3186128368, oid="Node4316369872")
    nd_4316369872.edge.oid = "Edge4316369936"
    nd_4316348240 = nd_4316369744.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=11.3186128368, oid="Node4316348240")
    nd_4316348240.edge.oid = "Edge4316370064"
    nd_4316370256 = nd_4316370128.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=21.0803560094, oid="Node4316370256")
    nd_4316370256.edge.oid = "Edge4316370320"
    nd_4316370000 = nd_4316370128.new_child(label="Node4316370000", taxon=None, edge_length=6.93019779087, oid="Node4316370000")
    nd_4316370000.edge.oid = "Edge4316370448"
    nd_4316370512 = nd_4316370000.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=14.1501582185, oid="Node4316370512")
    nd_4316370512.edge.oid = "Edge4316370576"
    nd_4316370384 = nd_4316370000.new_child(label="Node4316370384", taxon=None, edge_length=3.68709033509, oid="Node4316370384")
    nd_4316370384.edge.oid = "Edge4316347920"
    nd_4316347856 = nd_4316370384.new_child(label="Node4316347856", taxon=None, edge_length=4.34576715349, oid="Node4316347856")
    nd_4316347856.edge.oid = "Edge4316370704"
    nd_4316371152 = nd_4316370384.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=10.4630678834, oid="Node4316371152")
    nd_4316371152.edge.oid = "Edge4316371344"
    nd_4316370768 = nd_4316347856.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=6.11730072991, oid="Node4316370768")
    nd_4316370768.edge.oid = "Edge4316370832"
    nd_4316370640 = nd_4316347856.new_child(label="Node4316370640", taxon=None, edge_length=2.07374196069, oid="Node4316370640")
    nd_4316370640.edge.oid = "Edge4316370960"
    nd_4316371024 = nd_4316370640.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=4.04355876922, oid="Node4316371024")
    nd_4316371024.edge.oid = "Edge4316371088"
    nd_4316370896 = nd_4316370640.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=4.04355876922, oid="Node4316370896")
    nd_4316370896.edge.oid = "Edge4316371216"
    nd_4316371536 = nd_4316371280.new_child(label="Node4316371536", taxon=None, edge_length=5.79725116231, oid="Node4316371536")
    nd_4316371536.edge.oid = "Edge4316371600"
    nd_4316372048 = nd_4316371280.new_child(label="Node4316372048", taxon=None, edge_length=3.71811515781, oid="Node4316372048")
    nd_4316372048.edge.oid = "Edge4316372240"
    nd_4316371664 = nd_4316371536.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=22.4623407504, oid="Node4316371664")
    nd_4316371664.edge.oid = "Edge4316371728"
    nd_4316371408 = nd_4316371536.new_child(label="Node4316371408", taxon=None, edge_length=5.91246822929, oid="Node4316371408")
    nd_4316371408.edge.oid = "Edge4316371856"
    nd_4316371920 = nd_4316371408.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=16.5498725211, oid="Node4316371920")
    nd_4316371920.edge.oid = "Edge4316371984"
    nd_4316371792 = nd_4316371408.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=16.5498725211, oid="Node4316371792")
    nd_4316371792.edge.oid = "Edge4316372112"
    nd_4316372304 = nd_4316372048.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=24.5414767549, oid="Node4316372304")
    nd_4316372304.edge.oid = "Edge4316372368"
    nd_4316372176 = nd_4316372048.new_child(label="Node4316372176", taxon=None, edge_length=5.68085904285, oid="Node4316372176")
    nd_4316372176.edge.oid = "Edge4316372432"
    nd_4316372496 = nd_4316372176.new_child(label="Node4316372496", taxon=None, edge_length=11.551508813, oid="Node4316372496")
    nd_4316372496.edge.oid = "Edge4316372560"
    nd_4316372944 = nd_4316372176.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=18.8606177121, oid="Node4316372944")
    nd_4316372944.edge.oid = "Edge4316393552"
    nd_4316372624 = nd_4316372496.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=7.30910889905, oid="Node4316372624")
    nd_4316372624.edge.oid = "Edge4316372688"
    nd_4316372816 = nd_4316372496.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=7.30910889905, oid="Node4316372816")
    nd_4316372816.edge.oid = "Edge4316372880"
    nd_4316393744 = nd_4316372752.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=37.5960808765, oid="Node4316393744")
    nd_4316393744.edge.oid = "Edge4316393808"
    nd_4316393616 = nd_4316372752.new_child(label="Node4316393616", taxon=None, edge_length=2.73294846613, oid="Node4316393616")
    nd_4316393616.edge.oid = "Edge4316393872"
    nd_4316393936 = nd_4316393616.new_child(label="Node4316393936", taxon=None, edge_length=2.22916081797, oid="Node4316393936")
    nd_4316393936.edge.oid = "Edge4316394000"
    nd_4316395408 = nd_4316393616.new_child(label="Node4316395408", taxon=None, edge_length=13.2144705479, oid="Node4316395408")
    nd_4316395408.edge.oid = "Edge4316395472"
    nd_4316394064 = nd_4316393936.new_child(label="Node4316394064", taxon=None, edge_length=19.1660901297, oid="Node4316394064")
    nd_4316394064.edge.oid = "Edge4316394128"
    nd_4316394512 = nd_4316393936.new_child(label="Node4316394512", taxon=None, edge_length=7.55440984409, oid="Node4316394512")
    nd_4316394512.edge.oid = "Edge4316394576"
    nd_4316394192 = nd_4316394064.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=13.4678814627, oid="Node4316394192")
    nd_4316394192.edge.oid = "Edge4316394256"
    nd_4316394384 = nd_4316394064.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=13.4678814627, oid="Node4316394384")
    nd_4316394384.edge.oid = "Edge4316394448"
    nd_4316394640 = nd_4316394512.new_child(label="Node4316394640", taxon=None, edge_length=18.9933631628, oid="Node4316394640")
    nd_4316394640.edge.oid = "Edge4316394704"
    nd_4316395024 = nd_4316394512.new_child(label="Node4316395024", taxon=None, edge_length=5.77339164759, oid="Node4316395024")
    nd_4316395024.edge.oid = "Edge4316395088"
    nd_4316394768 = nd_4316394640.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=6.08619858554, oid="Node4316394768")
    nd_4316394768.edge.oid = "Edge4316394832"
    nd_4316394320 = nd_4316394640.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=6.08619858554, oid="Node4316394320")
    nd_4316394320.edge.oid = "Edge4316394960"
    nd_4316395152 = nd_4316395024.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=19.3061701007, oid="Node4316395152")
    nd_4316395152.edge.oid = "Edge4316395216"
    nd_4316394896 = nd_4316395024.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=19.3061701007, oid="Node4316394896")
    nd_4316394896.edge.oid = "Edge4316395344"
    nd_4316395536 = nd_4316395408.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=21.6486618625, oid="Node4316395536")
    nd_4316395536.edge.oid = "Edge4316395600"
    nd_4316395280 = nd_4316395408.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=21.6486618625, oid="Node4316395280")
    nd_4316395280.edge.oid = "Edge4316395728"
    nd_4316395920 = nd_4316395664.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=25.9434669707, oid="Node4316395920")
    nd_4316395920.edge.oid = "Edge4316395984"
    nd_4316395792 = nd_4316395664.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=25.9434669707, oid="Node4316395792")
    nd_4316395792.edge.oid = "Edge4316396112"
    nd_4316396304 = nd_4316396048.new_child(label="Node4316396304", taxon=None, edge_length=7.07244422644, oid="Node4316396304")
    nd_4316396304.edge.oid = "Edge4316396368"
    nd_4316396816 = nd_4316396048.new_child(label="Python regius", taxon=tax_4313760144, edge_length=37.5351776494, oid="Node4316396816")
    nd_4316396816.edge.oid = "Edge4316397008"
    nd_4316396432 = nd_4316396304.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=30.462733423, oid="Node4316396432")
    nd_4316396432.edge.oid = "Edge4316396496"
    nd_4316396176 = nd_4316396304.new_child(label="Node4316396176", taxon=None, edge_length=14.4432592042, oid="Node4316396176")
    nd_4316396176.edge.oid = "Edge4316396624"
    nd_4316396688 = nd_4316396176.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=16.0194742188, oid="Node4316396688")
    nd_4316396688.edge.oid = "Edge4316396752"
    nd_4316396560 = nd_4316396176.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=16.0194742188, oid="Node4316396560")
    nd_4316396560.edge.oid = "Edge4316396880"
    tree_4316397072 = dendropy.Tree(label="Tree08", taxon_set=tree_list.taxon_set, oid="Tree4316397072")
    tree_list.append(tree_4316397072, reindex_taxa=False)
    tree_4316397072.seed_node.oid = 'Node4316397456'
    nd_4316418128 = tree_4316397072.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=95.8502441646, oid="Node4316418128")
    nd_4316418128.edge.oid = "Edge4316418192"
    nd_4316418320 = tree_4316397072.seed_node.new_child(label="Node4316418320", taxon=None, edge_length=21.8741644934, oid="Node4316418320")
    nd_4316418320.edge.oid = "Edge4316418384"
    nd_4316418448 = nd_4316418320.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=73.9760796713, oid="Node4316418448")
    nd_4316418448.edge.oid = "Edge4316418512"
    nd_4316418256 = nd_4316418320.new_child(label="Node4316418256", taxon=None, edge_length=9.52951598189, oid="Node4316418256")
    nd_4316418256.edge.oid = "Edge4316418640"
    nd_4316418704 = nd_4316418256.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=64.4465636894, oid="Node4316418704")
    nd_4316418704.edge.oid = "Edge4316418768"
    nd_4316418576 = nd_4316418256.new_child(label="Node4316418576", taxon=None, edge_length=22.9882151659, oid="Node4316418576")
    nd_4316418576.edge.oid = "Edge4316418896"
    nd_4316418960 = nd_4316418576.new_child(label="Node4316418960", taxon=None, edge_length=6.43697452247, oid="Node4316418960")
    nd_4316418960.edge.oid = "Edge4316419024"
    nd_4316445776 = nd_4316418576.new_child(label="Node4316445776", taxon=None, edge_length=13.8002436509, oid="Node4316445776")
    nd_4316445776.edge.oid = "Edge4316445968"
    nd_4316419088 = nd_4316418960.new_child(label="Node4316419088", taxon=None, edge_length=6.7231540226, oid="Node4316419088")
    nd_4316419088.edge.oid = "Edge4316419152"
    nd_4316445392 = nd_4316418960.new_child(label="Node4316445392", taxon=None, edge_length=15.1871703268, oid="Node4316445392")
    nd_4316445392.edge.oid = "Edge4316445584"
    nd_4316397392 = nd_4316419088.new_child(label="Node4316397392", taxon=None, edge_length=1.03257897787, oid="Node4316397392")
    nd_4316397392.edge.oid = "Edge4316397328"
    nd_4316443088 = nd_4316419088.new_child(label="Node4316443088", taxon=None, edge_length=1.21568464049, oid="Node4316443088")
    nd_4316443088.edge.oid = "Edge4316443408"
    nd_4316419216 = nd_4316397392.new_child(label="Node4316419216", taxon=None, edge_length=7.52518318022, oid="Node4316419216")
    nd_4316419216.edge.oid = "Edge4316419280"
    nd_4316421136 = nd_4316397392.new_child(label="Node4316421136", taxon=None, edge_length=2.39608223465, oid="Node4316421136")
    nd_4316421136.edge.oid = "Edge4316421200"
    nd_4316419344 = nd_4316419216.new_child(label="Node4316419344", taxon=None, edge_length=2.42060556338, oid="Node4316419344")
    nd_4316419344.edge.oid = "Edge4316419408"
    nd_4316421008 = nd_4316419216.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=19.7404578203, oid="Node4316421008")
    nd_4316421008.edge.oid = "Edge4316421072"
    nd_4316419472 = nd_4316419344.new_child(label="Node4316419472", taxon=None, edge_length=8.29264517113, oid="Node4316419472")
    nd_4316419472.edge.oid = "Edge4316419536"
    nd_4316419856 = nd_4316419344.new_child(label="Node4316419856", taxon=None, edge_length=6.70163613113, oid="Node4316419856")
    nd_4316419856.edge.oid = "Edge4316419920"
    nd_4316419600 = nd_4316419472.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=9.02720708579, oid="Node4316419600")
    nd_4316419600.edge.oid = "Edge4316419664"
    nd_4316418832 = nd_4316419472.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=9.02720708579, oid="Node4316418832")
    nd_4316418832.edge.oid = "Edge4316419792"
    nd_4316419984 = nd_4316419856.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=10.6182161258, oid="Node4316419984")
    nd_4316419984.edge.oid = "Edge4316420048"
    nd_4316419728 = nd_4316419856.new_child(label="Node4316419728", taxon=None, edge_length=3.47880840545, oid="Node4316419728")
    nd_4316419728.edge.oid = "Edge4316420176"
    nd_4316420240 = nd_4316419728.new_child(label="Node4316420240", taxon=None, edge_length=4.27223311967, oid="Node4316420240")
    nd_4316420240.edge.oid = "Edge4316420304"
    nd_4316420752 = nd_4316419728.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=7.13940772034, oid="Node4316420752")
    nd_4316420752.edge.oid = "Edge4316420944"
    nd_4316420368 = nd_4316420240.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=2.86717460067, oid="Node4316420368")
    nd_4316420368.edge.oid = "Edge4316420432"
    nd_4316420112 = nd_4316420240.new_child(label="Node4316420112", taxon=None, edge_length=0.371215520464, oid="Node4316420112")
    nd_4316420112.edge.oid = "Edge4316420560"
    nd_4316420624 = nd_4316420112.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=2.49595908021, oid="Node4316420624")
    nd_4316420624.edge.oid = "Edge4316420688"
    nd_4316420496 = nd_4316420112.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.49595908021, oid="Node4316420496")
    nd_4316420496.edge.oid = "Edge4316420816"
    nd_4316421264 = nd_4316421136.new_child(label="Node4316421264", taxon=None, edge_length=3.24355281662, oid="Node4316421264")
    nd_4316421264.edge.oid = "Edge4316421328"
    nd_4316421776 = nd_4316421136.new_child(label="Node4316421776", taxon=None, edge_length=3.2027013644, oid="Node4316421776")
    nd_4316421776.edge.oid = "Edge4316421968"
    nd_4316421392 = nd_4316421264.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=21.6260059492, oid="Node4316421392")
    nd_4316421392.edge.oid = "Edge4316421456"
    nd_4316420880 = nd_4316421264.new_child(label="Node4316420880", taxon=None, edge_length=9.1533399099, oid="Node4316420880")
    nd_4316420880.edge.oid = "Edge4316421584"
    nd_4316421648 = nd_4316420880.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=12.4726660393, oid="Node4316421648")
    nd_4316421648.edge.oid = "Edge4316421712"
    nd_4316421520 = nd_4316420880.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=12.4726660393, oid="Node4316421520")
    nd_4316421520.edge.oid = "Edge4316421840"
    nd_4316422032 = nd_4316421776.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=21.6668574015, oid="Node4316422032")
    nd_4316422032.edge.oid = "Edge4316422096"
    nd_4316421904 = nd_4316421776.new_child(label="Node4316421904", taxon=None, edge_length=4.33602073619, oid="Node4316421904")
    nd_4316421904.edge.oid = "Edge4316442768"
    nd_4316442832 = nd_4316421904.new_child(label="Node4316442832", taxon=None, edge_length=11.5233229214, oid="Node4316442832")
    nd_4316442832.edge.oid = "Edge4316442896"
    nd_4316443216 = nd_4316421904.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=17.3308366653, oid="Node4316443216")
    nd_4316443216.edge.oid = "Edge4316443280"
    nd_4316442960 = nd_4316442832.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=5.8075137439, oid="Node4316442960")
    nd_4316442960.edge.oid = "Edge4316443024"
    nd_4316442704 = nd_4316442832.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=5.8075137439, oid="Node4316442704")
    nd_4316442704.edge.oid = "Edge4316443152"
    nd_4316443472 = nd_4316443088.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=27.0825353379, oid="Node4316443472")
    nd_4316443472.edge.oid = "Edge4316443536"
    nd_4316443344 = nd_4316443088.new_child(label="Node4316443344", taxon=None, edge_length=1.37229916019, oid="Node4316443344")
    nd_4316443344.edge.oid = "Edge4316443664"
    nd_4316443728 = nd_4316443344.new_child(label="Node4316443728", taxon=None, edge_length=2.64946637554, oid="Node4316443728")
    nd_4316443728.edge.oid = "Edge4316443792"
    nd_4316445136 = nd_4316443344.new_child(label="Node4316445136", taxon=None, edge_length=11.4050202795, oid="Node4316445136")
    nd_4316445136.edge.oid = "Edge4316445200"
    nd_4316443856 = nd_4316443728.new_child(label="Node4316443856", taxon=None, edge_length=13.5545767859, oid="Node4316443856")
    nd_4316443856.edge.oid = "Edge4316443920"
    nd_4316444240 = nd_4316443728.new_child(label="Node4316444240", taxon=None, edge_length=4.67390676307, oid="Node4316444240")
    nd_4316444240.edge.oid = "Edge4316444304"
    nd_4316443984 = nd_4316443856.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=9.50619301624, oid="Node4316443984")
    nd_4316443984.edge.oid = "Edge4316444048"
    nd_4316443600 = nd_4316443856.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=9.50619301624, oid="Node4316443600")
    nd_4316443600.edge.oid = "Edge4316444176"
    nd_4316444368 = nd_4316444240.new_child(label="Node4316444368", taxon=None, edge_length=12.8995814401, oid="Node4316444368")
    nd_4316444368.edge.oid = "Edge4316444432"
    nd_4316444752 = nd_4316444240.new_child(label="Node4316444752", taxon=None, edge_length=1.38849394051, oid="Node4316444752")
    nd_4316444752.edge.oid = "Edge4316444816"
    nd_4316444496 = nd_4316444368.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=5.48728159904, oid="Node4316444496")
    nd_4316444496.edge.oid = "Edge4316444560"
    nd_4316444112 = nd_4316444368.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=5.48728159904, oid="Node4316444112")
    nd_4316444112.edge.oid = "Edge4316444688"
    nd_4316444880 = nd_4316444752.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=16.9983690986, oid="Node4316444880")
    nd_4316444880.edge.oid = "Edge4316444944"
    nd_4316444624 = nd_4316444752.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=16.9983690986, oid="Node4316444624")
    nd_4316444624.edge.oid = "Edge4316445072"
    nd_4316445264 = nd_4316445136.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=14.3052158982, oid="Node4316445264")
    nd_4316445264.edge.oid = "Edge4316445328"
    nd_4316445008 = nd_4316445136.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=14.3052158982, oid="Node4316445008")
    nd_4316445008.edge.oid = "Edge4316445456"
    nd_4316445648 = nd_4316445392.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=19.8342036742, oid="Node4316445648")
    nd_4316445648.edge.oid = "Edge4316445712"
    nd_4316445520 = nd_4316445392.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=19.8342036742, oid="Node4316445520")
    nd_4316445520.edge.oid = "Edge4316445840"
    nd_4316446032 = nd_4316445776.new_child(label="Node4316446032", taxon=None, edge_length=3.93234771773, oid="Node4316446032")
    nd_4316446032.edge.oid = "Edge4316446096"
    nd_4316446544 = nd_4316445776.new_child(label="Python regius", taxon=tax_4313760144, edge_length=27.6581048725, oid="Node4316446544")
    nd_4316446544.edge.oid = "Edge4316467280"
    nd_4316446160 = nd_4316446032.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=23.7257571548, oid="Node4316446160")
    nd_4316446160.edge.oid = "Edge4316446224"
    nd_4316445904 = nd_4316446032.new_child(label="Node4316445904", taxon=None, edge_length=7.63118798191, oid="Node4316445904")
    nd_4316445904.edge.oid = "Edge4316446352"
    nd_4316446416 = nd_4316445904.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=16.0945691729, oid="Node4316446416")
    nd_4316446416.edge.oid = "Edge4316446480"
    nd_4316446288 = nd_4316445904.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=16.0945691729, oid="Node4316446288")
    nd_4316446288.edge.oid = "Edge4316446608"
    tree_4316446672 = dendropy.Tree(label="Tree09", taxon_set=tree_list.taxon_set, oid="Tree4316446672")
    tree_list.append(tree_4316446672, reindex_taxa=False)
    tree_4316446672.seed_node.oid = 'Node4316467472'
    nd_4316467600 = tree_4316446672.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=126.147943419, oid="Node4316467600")
    nd_4316467600.edge.oid = "Edge4316467664"
    nd_4316467344 = tree_4316446672.seed_node.new_child(label="Node4316467344", taxon=None, edge_length=40.5157695372, oid="Node4316467344")
    nd_4316467344.edge.oid = "Edge4316467728"
    nd_4316467792 = nd_4316467344.new_child(label="Node4316467792", taxon=None, edge_length=13.6608326978, oid="Node4316467792")
    nd_4316467792.edge.oid = "Edge4316467856"
    nd_4316512592 = nd_4316467344.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=85.6321738818, oid="Node4316512592")
    nd_4316512592.edge.oid = "Edge4316512656"
    nd_4316467920 = nd_4316467792.new_child(label="Node4316467920", taxon=None, edge_length=17.1918011574, oid="Node4316467920")
    nd_4316467920.edge.oid = "Edge4316467984"
    nd_4316512336 = nd_4316467792.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=71.971341184, oid="Node4316512336")
    nd_4316512336.edge.oid = "Edge4316512528"
    nd_4316468048 = nd_4316467920.new_child(label="Node4316468048", taxon=None, edge_length=5.77299961102, oid="Node4316468048")
    nd_4316468048.edge.oid = "Edge4316468112"
    nd_4316490896 = nd_4316467920.new_child(label="Node4316490896", taxon=None, edge_length=21.166651853, oid="Node4316490896")
    nd_4316490896.edge.oid = "Edge4316491088"
    nd_4316468176 = nd_4316468048.new_child(label="Node4316468176", taxon=None, edge_length=13.7010200713, oid="Node4316468176")
    nd_4316468176.edge.oid = "Edge4316468240"
    nd_4316490576 = nd_4316468048.new_child(label="Node4316490576", taxon=None, edge_length=24.6978541753, oid="Node4316490576")
    nd_4316490576.edge.oid = "Edge4316490704"
    nd_4316468304 = nd_4316468176.new_child(label="Node4316468304", taxon=None, edge_length=0.617658990864, oid="Node4316468304")
    nd_4316468304.edge.oid = "Edge4316468368"
    nd_4316470480 = nd_4316468176.new_child(label="Node4316470480", taxon=None, edge_length=1.35670146604, oid="Node4316470480")
    nd_4316470480.edge.oid = "Edge4316470544"
    nd_4316468432 = nd_4316468304.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=34.6878613533, oid="Node4316468432")
    nd_4316468432.edge.oid = "Edge4316468496"
    nd_4316467408 = nd_4316468304.new_child(label="Node4316467408", taxon=None, edge_length=2.23641376685, oid="Node4316467408")
    nd_4316467408.edge.oid = "Edge4316468624"
    nd_4316468688 = nd_4316467408.new_child(label="Node4316468688", taxon=None, edge_length=10.5685372807, oid="Node4316468688")
    nd_4316468688.edge.oid = "Edge4316468752"
    nd_4316469072 = nd_4316467408.new_child(label="Node4316469072", taxon=None, edge_length=2.11685229318, oid="Node4316469072")
    nd_4316469072.edge.oid = "Edge4316469136"
    nd_4316468816 = nd_4316468688.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=21.8829103058, oid="Node4316468816")
    nd_4316468816.edge.oid = "Edge4316468880"
    nd_4316468560 = nd_4316468688.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=21.8829103058, oid="Node4316468560")
    nd_4316468560.edge.oid = "Edge4316469008"
    nd_4316469200 = nd_4316469072.new_child(label="Node4316469200", taxon=None, edge_length=17.1757724208, oid="Node4316469200")
    nd_4316469200.edge.oid = "Edge4316469264"
    nd_4316469584 = nd_4316469072.new_child(label="Node4316469584", taxon=None, edge_length=3.96372676423, oid="Node4316469584")
    nd_4316469584.edge.oid = "Edge4316469648"
    nd_4316469328 = nd_4316469200.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=13.1588228725, oid="Node4316469328")
    nd_4316469328.edge.oid = "Edge4316469392"
    nd_4316468944 = nd_4316469200.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=13.1588228725, oid="Node4316468944")
    nd_4316468944.edge.oid = "Edge4316469520"
    nd_4316469712 = nd_4316469584.new_child(label="Node4316469712", taxon=None, edge_length=19.5683902852, oid="Node4316469712")
    nd_4316469712.edge.oid = "Edge4316469776"
    nd_4316470096 = nd_4316469584.new_child(label="Node4316470096", taxon=None, edge_length=4.82785669688, oid="Node4316470096")
    nd_4316470096.edge.oid = "Edge4316470160"
    nd_4316469840 = nd_4316469712.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=6.80247824391, oid="Node4316469840")
    nd_4316469840.edge.oid = "Edge4316469904"
    nd_4316469456 = nd_4316469712.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=6.80247824391, oid="Node4316469456")
    nd_4316469456.edge.oid = "Edge4316470032"
    nd_4316470224 = nd_4316470096.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=21.5430118322, oid="Node4316470224")
    nd_4316470224.edge.oid = "Edge4316470288"
    nd_4316469968 = nd_4316470096.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=21.5430118322, oid="Node4316469968")
    nd_4316469968.edge.oid = "Edge4316470416"
    nd_4316470608 = nd_4316470480.new_child(label="Node4316470608", taxon=None, edge_length=10.3401216686, oid="Node4316470608")
    nd_4316470608.edge.oid = "Edge4316470672"
    nd_4316488848 = nd_4316470480.new_child(label="Node4316488848", taxon=None, edge_length=7.03488295134, oid="Node4316488848")
    nd_4316488848.edge.oid = "Edge4316489040"
    nd_4316470736 = nd_4316470608.new_child(label="Node4316470736", taxon=None, edge_length=2.31259127041, oid="Node4316470736")
    nd_4316470736.edge.oid = "Edge4316470800"
    nd_4316488592 = nd_4316470608.new_child(label="Node4316488592", taxon=None, edge_length=13.9826784484, oid="Node4316488592")
    nd_4316488592.edge.oid = "Edge4316488656"
    nd_4316470864 = nd_4316470736.new_child(label="Node4316470864", taxon=None, edge_length=7.79402846351, oid="Node4316470864")
    nd_4316470864.edge.oid = "Edge4316470928"
    nd_4316488336 = nd_4316470736.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=21.2961059391, oid="Node4316488336")
    nd_4316488336.edge.oid = "Edge4316488528"
    nd_4316470992 = nd_4316470864.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=13.5020774756, oid="Node4316470992")
    nd_4316470992.edge.oid = "Edge4316471056"
    nd_4316470352 = nd_4316470864.new_child(label="Node4316470352", taxon=None, edge_length=3.5072599877, oid="Node4316470352")
    nd_4316470352.edge.oid = "Edge4316471184"
    nd_4316471248 = nd_4316470352.new_child(label="Node4316471248", taxon=None, edge_length=6.20487512451, oid="Node4316471248")
    nd_4316471248.edge.oid = "Edge4316487760"
    nd_4316488208 = nd_4316470352.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=9.99481748791, oid="Node4316488208")
    nd_4316488208.edge.oid = "Edge4316488400"
    nd_4316487824 = nd_4316471248.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=3.78994236341, oid="Node4316487824")
    nd_4316487824.edge.oid = "Edge4316487888"
    nd_4316471120 = nd_4316471248.new_child(label="Node4316471120", taxon=None, edge_length=1.25204312348, oid="Node4316471120")
    nd_4316471120.edge.oid = "Edge4316488016"
    nd_4316488080 = nd_4316471120.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=2.53789923993, oid="Node4316488080")
    nd_4316488080.edge.oid = "Edge4316488144"
    nd_4316487952 = nd_4316471120.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.53789923993, oid="Node4316487952")
    nd_4316487952.edge.oid = "Edge4316488272"
    nd_4316488720 = nd_4316488592.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=9.62601876119, oid="Node4316488720")
    nd_4316488720.edge.oid = "Edge4316488784"
    nd_4316488464 = nd_4316488592.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=9.62601876119, oid="Node4316488464")
    nd_4316488464.edge.oid = "Edge4316488912"
    nd_4316489104 = nd_4316488848.new_child(label="Node4316489104", taxon=None, edge_length=4.39672345568, oid="Node4316489104")
    nd_4316489104.edge.oid = "Edge4316489168"
    nd_4316489616 = nd_4316488848.new_child(label="Node4316489616", taxon=None, edge_length=3.94778865925, oid="Node4316489616")
    nd_4316489616.edge.oid = "Edge4316489808"
    nd_4316489232 = nd_4316489104.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=22.5172124711, oid="Node4316489232")
    nd_4316489232.edge.oid = "Edge4316489296"
    nd_4316488976 = nd_4316489104.new_child(label="Node4316488976", taxon=None, edge_length=11.9061835258, oid="Node4316488976")
    nd_4316488976.edge.oid = "Edge4316489424"
    nd_4316489488 = nd_4316488976.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=10.6110289453, oid="Node4316489488")
    nd_4316489488.edge.oid = "Edge4316489552"
    nd_4316489360 = nd_4316488976.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=10.6110289453, oid="Node4316489360")
    nd_4316489360.edge.oid = "Edge4316489680"
    nd_4316489872 = nd_4316489616.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=22.9661472676, oid="Node4316489872")
    nd_4316489872.edge.oid = "Edge4316489936"
    nd_4316489744 = nd_4316489616.new_child(label="Node4316489744", taxon=None, edge_length=3.7713704432, oid="Node4316489744")
    nd_4316489744.edge.oid = "Edge4316490064"
    nd_4316490128 = nd_4316489744.new_child(label="Node4316490128", taxon=None, edge_length=11.9409692012, oid="Node4316490128")
    nd_4316490128.edge.oid = "Edge4316490192"
    nd_4316490512 = nd_4316489744.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=19.1947768244, oid="Node4316490512")
    nd_4316490512.edge.oid = "Edge4316490384"
    nd_4316490256 = nd_4316490128.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=7.25380762314, oid="Node4316490256")
    nd_4316490256.edge.oid = "Edge4316490320"
    nd_4316490000 = nd_4316490128.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=7.25380762314, oid="Node4316490000")
    nd_4316490000.edge.oid = "Edge4316490448"
    nd_4316490768 = nd_4316490576.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=24.3086862402, oid="Node4316490768")
    nd_4316490768.edge.oid = "Edge4316490832"
    nd_4316490640 = nd_4316490576.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=24.3086862402, oid="Node4316490640")
    nd_4316490640.edge.oid = "Edge4316490960"
    nd_4316491152 = nd_4316490896.new_child(label="Node4316491152", taxon=None, edge_length=5.80983351578, oid="Node4316491152")
    nd_4316491152.edge.oid = "Edge4316491216"
    nd_4316491664 = nd_4316490896.new_child(label="Python regius", taxon=tax_4313760144, edge_length=33.6128881735, oid="Node4316491664")
    nd_4316491664.edge.oid = "Edge4316512400"
    nd_4316491280 = nd_4316491152.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=27.8030546577, oid="Node4316491280")
    nd_4316491280.edge.oid = "Edge4316491344"
    nd_4316491024 = nd_4316491152.new_child(label="Node4316491024", taxon=None, edge_length=10.4395768579, oid="Node4316491024")
    nd_4316491024.edge.oid = "Edge4316491472"
    nd_4316491536 = nd_4316491024.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=17.3634777999, oid="Node4316491536")
    nd_4316491536.edge.oid = "Edge4316491600"
    nd_4316491408 = nd_4316491024.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=17.3634777999, oid="Node4316491408")
    nd_4316491408.edge.oid = "Edge4316491728"
    tree_4316512464 = dendropy.Tree(label="Tree10", taxon_set=tree_list.taxon_set, oid="Tree4316512464")
    tree_list.append(tree_4316512464, reindex_taxa=False)
    tree_4316512464.seed_node.oid = 'Node4316512848'
    nd_4316512976 = tree_4316512464.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=146.770054852, oid="Node4316512976")
    nd_4316512976.edge.oid = "Edge4316513040"
    nd_4316513168 = tree_4316512464.seed_node.new_child(label="Node4316513168", taxon=None, edge_length=49.9930471528, oid="Node4316513168")
    nd_4316513168.edge.oid = "Edge4316513232"
    nd_4316513296 = nd_4316513168.new_child(label="Node4316513296", taxon=None, edge_length=13.4634525107, oid="Node4316513296")
    nd_4316513296.edge.oid = "Edge4316513360"
    nd_4316558096 = nd_4316513168.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=96.7770076997, oid="Node4316558096")
    nd_4316558096.edge.oid = "Edge4316558160"
    nd_4316513424 = nd_4316513296.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=83.313555189, oid="Node4316513424")
    nd_4316513424.edge.oid = "Edge4316513488"
    nd_4316513104 = nd_4316513296.new_child(label="Node4316513104", taxon=None, edge_length=30.4776798161, oid="Node4316513104")
    nd_4316513104.edge.oid = "Edge4316513616"
    nd_4316513680 = nd_4316513104.new_child(label="Node4316513680", taxon=None, edge_length=8.15039409252, oid="Node4316513680")
    nd_4316513680.edge.oid = "Edge4316513744"
    nd_4316536464 = nd_4316513104.new_child(label="Node4316536464", taxon=None, edge_length=14.4274535052, oid="Node4316536464")
    nd_4316536464.edge.oid = "Edge4316536592"
    nd_4316513808 = nd_4316513680.new_child(label="Node4316513808", taxon=None, edge_length=9.51023675819, oid="Node4316513808")
    nd_4316513808.edge.oid = "Edge4316513872"
    nd_4316512784 = nd_4316513680.new_child(label="Node4316512784", taxon=None, edge_length=25.3895116055, oid="Node4316512784")
    nd_4316512784.edge.oid = "Edge4316536272"
    nd_4316513936 = nd_4316513808.new_child(label="Node4316513936", taxon=None, edge_length=1.94845077373, oid="Node4316513936")
    nd_4316513936.edge.oid = "Edge4316514000"
    nd_4316534096 = nd_4316513808.new_child(label="Node4316534096", taxon=None, edge_length=1.20705242999, oid="Node4316534096")
    nd_4316534096.edge.oid = "Edge4316534160"
    nd_4316514064 = nd_4316513936.new_child(label="Node4316514064", taxon=None, edge_length=3.25905094181, oid="Node4316514064")
    nd_4316514064.edge.oid = "Edge4316514128"
    nd_4316515728 = nd_4316513936.new_child(label="Node4316515728", taxon=None, edge_length=10.0348009179, oid="Node4316515728")
    nd_4316515728.edge.oid = "Edge4316515792"
    nd_4316514192 = nd_4316514064.new_child(label="Node4316514192", taxon=None, edge_length=4.82232736143, oid="Node4316514192")
    nd_4316514192.edge.oid = "Edge4316514256"
    nd_4316514704 = nd_4316514064.new_child(label="Node4316514704", taxon=None, edge_length=4.58004825968, oid="Node4316514704")
    nd_4316514704.edge.oid = "Edge4316514896"
    nd_4316514320 = nd_4316514192.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=25.1454154452, oid="Node4316514320")
    nd_4316514320.edge.oid = "Edge4316514384"
    nd_4316513552 = nd_4316514192.new_child(label="Node4316513552", taxon=None, edge_length=10.5788836702, oid="Node4316513552")
    nd_4316513552.edge.oid = "Edge4316514512"
    nd_4316514576 = nd_4316513552.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=14.566531775, oid="Node4316514576")
    nd_4316514576.edge.oid = "Edge4316514640"
    nd_4316514448 = nd_4316513552.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=14.566531775, oid="Node4316514448")
    nd_4316514448.edge.oid = "Edge4316514768"
    nd_4316514960 = nd_4316514704.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=25.387694547, oid="Node4316514960")
    nd_4316514960.edge.oid = "Edge4316515024"
    nd_4316514832 = nd_4316514704.new_child(label="Node4316514832", taxon=None, edge_length=5.75635440071, oid="Node4316514832")
    nd_4316514832.edge.oid = "Edge4316515152"
    nd_4316515216 = nd_4316514832.new_child(label="Node4316515216", taxon=None, edge_length=12.7115868187, oid="Node4316515216")
    nd_4316515216.edge.oid = "Edge4316515280"
    nd_4316515600 = nd_4316514832.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=19.6313401463, oid="Node4316515600")
    nd_4316515600.edge.oid = "Edge4316515664"
    nd_4316515344 = nd_4316515216.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=6.91975332756, oid="Node4316515344")
    nd_4316515344.edge.oid = "Edge4316515408"
    nd_4316515088 = nd_4316515216.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=6.91975332756, oid="Node4316515088")
    nd_4316515088.edge.oid = "Edge4316515536"
    nd_4316515856 = nd_4316515728.new_child(label="Node4316515856", taxon=None, edge_length=14.4811675935, oid="Node4316515856")
    nd_4316515856.edge.oid = "Edge4316515920"
    nd_4316516240 = nd_4316515728.new_child(label="Node4316516240", taxon=None, edge_length=0.795030531054, oid="Node4316516240")
    nd_4316516240.edge.oid = "Edge4316516304"
    nd_4316515984 = nd_4316515856.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=8.71082523711, oid="Node4316515984")
    nd_4316515984.edge.oid = "Edge4316516048"
    nd_4316515472 = nd_4316515856.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=8.71082523711, oid="Node4316515472")
    nd_4316515472.edge.oid = "Edge4316516176"
    nd_4316532816 = nd_4316516240.new_child(label="Node4316532816", taxon=None, edge_length=7.95009719313, oid="Node4316532816")
    nd_4316532816.edge.oid = "Edge4316532880"
    nd_4316533840 = nd_4316516240.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=22.3969622995, oid="Node4316533840")
    nd_4316533840.edge.oid = "Edge4316534032"
    nd_4316532944 = nd_4316532816.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=14.4468651064, oid="Node4316532944")
    nd_4316532944.edge.oid = "Edge4316533008"
    nd_4316516112 = nd_4316532816.new_child(label="Node4316516112", taxon=None, edge_length=2.27822479328, oid="Node4316516112")
    nd_4316516112.edge.oid = "Edge4316533136"
    nd_4316533200 = nd_4316516112.new_child(label="Node4316533200", taxon=None, edge_length=6.98040063458, oid="Node4316533200")
    nd_4316533200.edge.oid = "Edge4316533264"
    nd_4316533648 = nd_4316516112.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=12.1686403131, oid="Node4316533648")
    nd_4316533648.edge.oid = "Edge4316533904"
    nd_4316533328 = nd_4316533200.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=5.18823967855, oid="Node4316533328")
    nd_4316533328.edge.oid = "Edge4316533392"
    nd_4316533072 = nd_4316533200.new_child(label="Node4316533072", taxon=None, edge_length=1.71064757594, oid="Node4316533072")
    nd_4316533072.edge.oid = "Edge4316533456"
    nd_4316533520 = nd_4316533072.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=3.47759210261, oid="Node4316533520")
    nd_4316533520.edge.oid = "Edge4316533584"
    nd_4316533712 = nd_4316533072.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=3.47759210261, oid="Node4316533712")
    nd_4316533712.edge.oid = "Edge4316533776"
    nd_4316534224 = nd_4316534096.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=33.9681920922, oid="Node4316534224")
    nd_4316534224.edge.oid = "Edge4316534288"
    nd_4316533968 = nd_4316534096.new_child(label="Node4316533968", taxon=None, edge_length=1.28486469944, oid="Node4316533968")
    nd_4316533968.edge.oid = "Edge4316534416"
    nd_4316534480 = nd_4316533968.new_child(label="Node4316534480", taxon=None, edge_length=12.4520799939, oid="Node4316534480")
    nd_4316534480.edge.oid = "Edge4316534544"
    nd_4316534864 = nd_4316533968.new_child(label="Node4316534864", taxon=None, edge_length=1.68023264943, oid="Node4316534864")
    nd_4316534864.edge.oid = "Edge4316534928"
    nd_4316534608 = nd_4316534480.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=20.2312473989, oid="Node4316534608")
    nd_4316534608.edge.oid = "Edge4316534672"
    nd_4316534352 = nd_4316534480.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=20.2312473989, oid="Node4316534352")
    nd_4316534352.edge.oid = "Edge4316534800"
    nd_4316534992 = nd_4316534864.new_child(label="Node4316534992", taxon=None, edge_length=19.0383478987, oid="Node4316534992")
    nd_4316534992.edge.oid = "Edge4316535056"
    nd_4316535376 = nd_4316534864.new_child(label="Node4316535376", taxon=None, edge_length=4.75943584051, oid="Node4316535376")
    nd_4316535376.edge.oid = "Edge4316535440"
    nd_4316535120 = nd_4316534992.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=11.9647468446, oid="Node4316535120")
    nd_4316535120.edge.oid = "Edge4316535184"
    nd_4316534736 = nd_4316534992.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=11.9647468446, oid="Node4316534736")
    nd_4316534736.edge.oid = "Edge4316535312"
    nd_4316535504 = nd_4316535376.new_child(label="Node4316535504", taxon=None, edge_length=17.1180008393, oid="Node4316535504")
    nd_4316535504.edge.oid = "Edge4316535568"
    nd_4316535888 = nd_4316535376.new_child(label="Node4316535888", taxon=None, edge_length=2.6817531508, oid="Node4316535888")
    nd_4316535888.edge.oid = "Edge4316535952"
    nd_4316535632 = nd_4316535504.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=9.12565806357, oid="Node4316535632")
    nd_4316535632.edge.oid = "Edge4316535696"
    nd_4316535248 = nd_4316535504.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=9.12565806357, oid="Node4316535248")
    nd_4316535248.edge.oid = "Edge4316535824"
    nd_4316536016 = nd_4316535888.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=23.561905752, oid="Node4316536016")
    nd_4316536016.edge.oid = "Edge4316536080"
    nd_4316512720 = nd_4316535888.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=23.561905752, oid="Node4316512720")
    nd_4316512720.edge.oid = "Edge4316535760"
    nd_4316536336 = nd_4316512784.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=19.2959696748, oid="Node4316536336")
    nd_4316536336.edge.oid = "Edge4316536400"
    nd_4316536208 = nd_4316512784.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=19.2959696748, oid="Node4316536208")
    nd_4316536208.edge.oid = "Edge4316536528"
    nd_4316536656 = nd_4316536464.new_child(label="Node4316536656", taxon=None, edge_length=6.97033769, oid="Node4316536656")
    nd_4316536656.edge.oid = "Edge4316536720"
    nd_4316557776 = nd_4316536464.new_child(label="Python regius", taxon=tax_4313760144, edge_length=38.4084218677, oid="Node4316557776")
    nd_4316557776.edge.oid = "Edge4316557968"
    nd_4316536784 = nd_4316536656.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=31.4380841777, oid="Node4316536784")
    nd_4316536784.edge.oid = "Edge4316557392"
    nd_4316557520 = nd_4316536656.new_child(label="Node4316557520", taxon=None, edge_length=13.5557299793, oid="Node4316557520")
    nd_4316557520.edge.oid = "Edge4316557584"
    nd_4316557648 = nd_4316557520.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=17.8823541984, oid="Node4316557648")
    nd_4316557648.edge.oid = "Edge4316557712"
    nd_4316557456 = nd_4316557520.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=17.8823541984, oid="Node4316557456")
    nd_4316557456.edge.oid = "Edge4316557840"
    tree_4316557904 = dendropy.Tree(label="Tree11", taxon_set=tree_list.taxon_set, oid="Tree4316557904")
    tree_list.append(tree_4316557904, reindex_taxa=False)
    tree_4316557904.seed_node.oid = 'Node4316558352'
    nd_4316558480 = tree_4316557904.seed_node.new_child(label="Node4316558480", taxon=None, edge_length=71.3865451194, oid="Node4316558480")
    nd_4316558480.edge.oid = "Edge4316558544"
    nd_4316607440 = tree_4316557904.seed_node.new_child(label="Candola aspera", taxon=tax_4313742032, edge_length=211.46817309, oid="Node4316607440")
    nd_4316607440.edge.oid = "Edge4316607632"
    nd_4316558608 = nd_4316558480.new_child(label="Xenopeltis unicolor", taxon=tax_4313760592, edge_length=140.081627971, oid="Node4316558608")
    nd_4316558608.edge.oid = "Edge4316558672"
    nd_4316558800 = nd_4316558480.new_child(label="Node4316558800", taxon=None, edge_length=26.0234610563, oid="Node4316558800")
    nd_4316558800.edge.oid = "Edge4316558864"
    nd_4316558928 = nd_4316558800.new_child(label="Loxocemus bicolor", taxon=tax_4313742480, edge_length=114.058166915, oid="Node4316558928")
    nd_4316558928.edge.oid = "Edge4316558992"
    nd_4316558736 = nd_4316558800.new_child(label="Node4316558736", taxon=None, edge_length=40.1699176918, oid="Node4316558736")
    nd_4316558736.edge.oid = "Edge4316559120"
    nd_4316559184 = nd_4316558736.new_child(label="Node4316559184", taxon=None, edge_length=12.7558559915, oid="Node4316559184")
    nd_4316559184.edge.oid = "Edge4316559248"
    nd_4316606544 = nd_4316558736.new_child(label="Node4316606544", taxon=None, edge_length=17.7722054357, oid="Node4316606544")
    nd_4316606544.edge.oid = "Edge4316606736"
    nd_4316559312 = nd_4316559184.new_child(label="Node4316559312", taxon=None, edge_length=12.1691499629, oid="Node4316559312")
    nd_4316559312.edge.oid = "Edge4316559376"
    nd_4316585744 = nd_4316559184.new_child(label="Node4316585744", taxon=None, edge_length=28.055060848, oid="Node4316585744")
    nd_4316585744.edge.oid = "Edge4316585808"
    nd_4316559440 = nd_4316559312.new_child(label="Node4316559440", taxon=None, edge_length=4.4227398576, oid="Node4316559440")
    nd_4316559440.edge.oid = "Edge4316559504"
    nd_4316583696 = nd_4316559312.new_child(label="Node4316583696", taxon=None, edge_length=3.26334313874, oid="Node4316583696")
    nd_4316583696.edge.oid = "Edge4316583760"
    nd_4316559568 = nd_4316559440.new_child(label="Node4316559568", taxon=None, edge_length=5.98387013923, oid="Node4316559568")
    nd_4316559568.edge.oid = "Edge4316559632"
    nd_4316561232 = nd_4316559440.new_child(label="Node4316561232", taxon=None, edge_length=17.9619739892, oid="Node4316561232")
    nd_4316561232.edge.oid = "Edge4316561296"
    nd_4316559696 = nd_4316559568.new_child(label="Node4316559696", taxon=None, edge_length=7.53311752135, oid="Node4316559696")
    nd_4316559696.edge.oid = "Edge4316559760"
    nd_4316560208 = nd_4316559568.new_child(label="Node4316560208", taxon=None, edge_length=4.18557870369, oid="Node4316560208")
    nd_4316560208.edge.oid = "Edge4316560400"
    nd_4316559824 = nd_4316559696.new_child(label="Morelia carinata", taxon=tax_4313742800, edge_length=31.0235157502, oid="Node4316559824")
    nd_4316559824.edge.oid = "Edge4316559888"
    nd_4316559056 = nd_4316559696.new_child(label="Node4316559056", taxon=None, edge_length=12.7654788618, oid="Node4316559056")
    nd_4316559056.edge.oid = "Edge4316560016"
    nd_4316560080 = nd_4316559056.new_child(label="Morelia viridisS", taxon=tax_4313759952, edge_length=18.2580368884, oid="Node4316560080")
    nd_4316560080.edge.oid = "Edge4316560144"
    nd_4316559952 = nd_4316559056.new_child(label="Morelia viridisN", taxon=tax_4313759888, edge_length=18.2580368884, oid="Node4316559952")
    nd_4316559952.edge.oid = "Edge4316560272"
    nd_4316560464 = nd_4316560208.new_child(label="Antaresia maculosa", taxon=tax_4313741328, edge_length=34.3710545678, oid="Node4316560464")
    nd_4316560464.edge.oid = "Edge4316560528"
    nd_4316560336 = nd_4316560208.new_child(label="Node4316560336", taxon=None, edge_length=7.81168349566, oid="Node4316560336")
    nd_4316560336.edge.oid = "Edge4316560656"
    nd_4316560720 = nd_4316560336.new_child(label="Node4316560720", taxon=None, edge_length=15.350690271, oid="Node4316560720")
    nd_4316560720.edge.oid = "Edge4316560784"
    nd_4316561104 = nd_4316560336.new_child(label="Antaresia perthensis", taxon=tax_4313741584, edge_length=26.5593710722, oid="Node4316561104")
    nd_4316561104.edge.oid = "Edge4316561168"
    nd_4316560848 = nd_4316560720.new_child(label="Antaresia childreni", taxon=tax_4313741136, edge_length=11.2086808012, oid="Node4316560848")
    nd_4316560848.edge.oid = "Edge4316560912"
    nd_4316560592 = nd_4316560720.new_child(label="Antaresia stimsoni", taxon=tax_4313741840, edge_length=11.2086808012, oid="Node4316560592")
    nd_4316560592.edge.oid = "Edge4316561040"
    nd_4316561360 = nd_4316561232.new_child(label="Node4316561360", taxon=None, edge_length=1.39777518086, oid="Node4316561360")
    nd_4316561360.edge.oid = "Edge4316581968"
    nd_4316583312 = nd_4316561232.new_child(label="Node4316583312", taxon=None, edge_length=13.6008570384, oid="Node4316583312")
    nd_4316583312.edge.oid = "Edge4316583376"
    nd_4316582032 = nd_4316561360.new_child(label="Morelia oenpelliensis", taxon=tax_4313743248, edge_length=25.1807542407, oid="Node4316582032")
    nd_4316582032.edge.oid = "Edge4316582096"
    nd_4316560976 = nd_4316561360.new_child(label="Node4316560976", taxon=None, edge_length=7.6242060025, oid="Node4316560976")
    nd_4316560976.edge.oid = "Edge4316582224"
    nd_4316582288 = nd_4316560976.new_child(label="Morelia tracyae", taxon=tax_4313759824, edge_length=17.5565482382, oid="Node4316582288")
    nd_4316582288.edge.oid = "Edge4316582352"
    nd_4316582160 = nd_4316560976.new_child(label="Node4316582160", taxon=None, edge_length=3.73213849687, oid="Node4316582160")
    nd_4316582160.edge.oid = "Edge4316582480"
    nd_4316582544 = nd_4316582160.new_child(label="Node4316582544", taxon=None, edge_length=8.62088071739, oid="Node4316582544")
    nd_4316582544.edge.oid = "Edge4316582608"
    nd_4316583056 = nd_4316582160.new_child(label="Morelia amethistina", taxon=tax_4313742608, edge_length=13.8244097414, oid="Node4316583056")
    nd_4316583056.edge.oid = "Edge4316583248"
    nd_4316582672 = nd_4316582544.new_child(label="Morelia clastolepis", taxon=tax_4313742928, edge_length=5.20352902397, oid="Node4316582672")
    nd_4316582672.edge.oid = "Edge4316582736"
    nd_4316582416 = nd_4316582544.new_child(label="Node4316582416", taxon=None, edge_length=2.83199057731, oid="Node4316582416")
    nd_4316582416.edge.oid = "Edge4316582864"
    nd_4316582928 = nd_4316582416.new_child(label="Morelia kinghorni", taxon=tax_4313743056, edge_length=2.37153844665, oid="Node4316582928")
    nd_4316582928.edge.oid = "Edge4316582992"
    nd_4316582800 = nd_4316582416.new_child(label="Morelia nauta", taxon=tax_4313743120, edge_length=2.37153844665, oid="Node4316582800")
    nd_4316582800.edge.oid = "Edge4316583120"
    nd_4316583440 = nd_4316583312.new_child(label="Morelia spilota", taxon=tax_4313743312, edge_length=12.9776723832, oid="Node4316583440")
    nd_4316583440.edge.oid = "Edge4316583504"
    nd_4316583184 = nd_4316583312.new_child(label="Morelia bredli", taxon=tax_4313742736, edge_length=12.9776723832, oid="Node4316583184")
    nd_4316583184.edge.oid = "Edge4316583632"
    nd_4316583824 = nd_4316583696.new_child(label="Node4316583824", taxon=None, edge_length=3.5026210108, oid="Node4316583824")
    nd_4316583824.edge.oid = "Edge4316583888"
    nd_4316585104 = nd_4316583696.new_child(label="Node4316585104", taxon=None, edge_length=0.85376044557, oid="Node4316585104")
    nd_4316585104.edge.oid = "Edge4316585168"
    nd_4316583952 = nd_4316583824.new_child(label="Node4316583952", taxon=None, edge_length=23.3327851176, oid="Node4316583952")
    nd_4316583952.edge.oid = "Edge4316584016"
    nd_4316584336 = nd_4316583824.new_child(label="Node4316584336", taxon=None, edge_length=12.8697548268, oid="Node4316584336")
    nd_4316584336.edge.oid = "Edge4316584400"
    nd_4316584080 = nd_4316583952.new_child(label="Antaresia melanocephalus", taxon=tax_4313741456, edge_length=18.8644940012, oid="Node4316584080")
    nd_4316584080.edge.oid = "Edge4316584144"
    nd_4316583568 = nd_4316583952.new_child(label="Antaresia ramsayi", taxon=tax_4313741712, edge_length=18.8644940012, oid="Node4316583568")
    nd_4316583568.edge.oid = "Edge4316584272"
    nd_4316584464 = nd_4316584336.new_child(label="Node4316584464", taxon=None, edge_length=22.0000981488, oid="Node4316584464")
    nd_4316584464.edge.oid = "Edge4316584528"
    nd_4316558288 = nd_4316584336.new_child(label="Node4316558288", taxon=None, edge_length=5.24037108992, oid="Node4316558288")
    nd_4316558288.edge.oid = "Edge4316584720"
    nd_4316584592 = nd_4316584464.new_child(label="Liasis fuscus", taxon=tax_4313742224, edge_length=7.32742614319, oid="Node4316584592")
    nd_4316584592.edge.oid = "Edge4316584656"
    nd_4316584208 = nd_4316584464.new_child(label="Liasis mackloti", taxon=tax_4313742288, edge_length=7.32742614319, oid="Node4316584208")
    nd_4316584208.edge.oid = "Edge4316584784"
    nd_4316584848 = nd_4316558288.new_child(label="Apodora papuana", taxon=tax_4313741904, edge_length=24.0871532021, oid="Node4316584848")
    nd_4316584848.edge.oid = "Edge4316584912"
    nd_4316558224 = nd_4316558288.new_child(label="Liasis olivaceus", taxon=tax_4313742352, edge_length=24.0871532021, oid="Node4316558224")
    nd_4316558224.edge.oid = "Edge4316585040"
    nd_4316585232 = nd_4316585104.new_child(label="Node4316585232", taxon=None, edge_length=19.3585494438, oid="Node4316585232")
    nd_4316585232.edge.oid = "Edge4316585296"
    nd_4316585616 = nd_4316585104.new_child(label="Morelia boeleni", taxon=tax_4313742672, edge_length=44.846139684, oid="Node4316585616")
    nd_4316585616.edge.oid = "Edge4316585680"
    nd_4316585360 = nd_4316585232.new_child(label="Bothrochilus boa", taxon=tax_4313741968, edge_length=25.4875902402, oid="Node4316585360")
    nd_4316585360.edge.oid = "Edge4316585424"
    nd_4316584976 = nd_4316585232.new_child(label="Liasis albertisii", taxon=tax_4313742160, edge_length=25.4875902402, oid="Node4316584976")
    nd_4316584976.edge.oid = "Edge4316585552"
    nd_4316585872 = nd_4316585744.new_child(label="Python timoriensis", taxon=tax_4313760464, edge_length=33.0773323832, oid="Node4316585872")
    nd_4316585872.edge.oid = "Edge4316585936"
    nd_4316585488 = nd_4316585744.new_child(label="Python reticulatus", taxon=tax_4313760272, edge_length=33.0773323832, oid="Node4316585488")
    nd_4316585488.edge.oid = "Edge4316606608"
    nd_4316606800 = nd_4316606544.new_child(label="Node4316606800", taxon=None, edge_length=13.6364666177, oid="Node4316606800")
    nd_4316606800.edge.oid = "Edge4316606864"
    nd_4316607312 = nd_4316606544.new_child(label="Python regius", taxon=tax_4313760144, edge_length=56.1160437871, oid="Node4316607312")
    nd_4316607312.edge.oid = "Edge4316607504"
    nd_4316606928 = nd_4316606800.new_child(label="Python curtus", taxon=tax_4313760016, edge_length=42.4795771694, oid="Node4316606928")
    nd_4316606928.edge.oid = "Edge4316606992"
    nd_4316606672 = nd_4316606800.new_child(label="Node4316606672", taxon=None, edge_length=16.6495052056, oid="Node4316606672")
    nd_4316606672.edge.oid = "Edge4316607120"
    nd_4316607184 = nd_4316606672.new_child(label="Python sebae", taxon=tax_4313760336, edge_length=25.8300719638, oid="Node4316607184")
    nd_4316607184.edge.oid = "Edge4316607248"
    nd_4316607056 = nd_4316606672.new_child(label="Python molurus", taxon=tax_4313760080, edge_length=25.8300719638, oid="Node4316607056")
    nd_4316607056.edge.oid = "Edge4316607376"

    # set labels of nodes with taxa to taxon label, else oid (for consistent
    # identification in debugging)
    for t in tree_list:
        t.assign_node_labels_from_taxon_or_oid()

    return tree_list


def reference_continuous_dict():
    val_dict = {}
    val_dict["Python molurus"] = [-0.08956193469026609, -0.28853328529693084, -0.026318238908102854, -0.1616168944925791, -0.08745045274904639, 0.21900696711880685, -0.4250888576431998, 0.12489737960863644, 0.23656542058445135, -0.18529445143751505, 0.05007982819570443, -0.04256927071540935, -0.7758625303514285, -0.13295463608223285, 0.22750687284929147, 0.11171990225139619, -0.05939327657280469, 0.4143105353728031, -0.29176115271231534, 0.17547718867426998, 0.1303324212218117, -0.11761014406847911, 0.2893407037611951, -0.9558654990063262, 0.2979897591694161, 0.10966108573965055, 0.5545675185446975, -0.5483632765139086, -0.06919318399096928, -0.2933349574456793, 0.011627517097510578, 0.25139736074648145, -0.49236401122642415, 0.4076908676571471, -0.024508799819400035, -0.12194022047098237, -0.13811559640164686, -0.017444561311853887, -0.1502547314031303, -0.1890324840476297, -0.083063689907993, -0.05269469924105295, 0.22006852486694065, 0.4575694069092294, 0.17970748702031067, -0.00401967756799814, -0.25165681730083783, -0.11393639250878206, 0.08433650666546753, 0.24898509997522997, -0.04766142562142892, 0.05540276033780106, 0.2861957754020195, -0.035790991480086125, 0.1577369485338326, 0.1693985828586302, -0.14255198505209465, 0.18567524069453895, 0.035504492026016424, 0.06735363520833434, -0.3514024644717678, -0.023985055897920198, 0.025022452646233508, -0.03916426209375296, 0.05058529302542099, 0.08232627724996472, -0.07060656429668419, -0.08008459485144284, -0.3224507234506318, -0.040343736134774315, -0.012089796980072402, -0.21605050117195093, 0.031453681103264075, 0.2153382331572617, 0.016182647236383813, 0.34937989760611166, 0.012435781315143807, 0.4545557650519242, -0.14206233728717885, -0.05575142869407321, -0.1864135496173435, 0.22240109593355464, 0.5849846786652463, 0.12966864574755094, -0.3057330194548949, 0.1053100317983198, 0.30403433551196746, -0.4263643753105065, 0.2670091929818761, -0.1572411016491901, -0.6098105918532315, -0.006806433026180957, 0.28917920313659745, 0.23229130959627, -0.0412871298810328, 0.4227013818414445, -0.386175284797848, -0.4926632378619112, -0.24777502213421387, 0.19540017504788812],
    val_dict["Python sebae"] = [0.02079349205340751, 0.5905776220171024, 0.35299882804702665, -0.3619374731335605, -0.2609679093527202, -0.02975274323873578, -0.12252115708503611, 0.06099880918276404, -0.15038528296069226, 0.177726922137512, 0.9103667327297166, 0.09586032954160995, -0.13405314905050766, 0.13602550376368228, -0.47571319781521115, 0.3876430354058918, 0.41727800763777867, 0.08467580385887138, -0.07887272313429783, 0.24929393763269375, 0.6073190855818995, -0.4527792111677741, -0.3774028894642513, -0.18304539133427167, 0.18144735168566664, 0.205856248552433, 0.15308629353821845, -0.42304135377192226, 0.15241269937561885, 0.3577616531301235, 0.4426614554711241, -0.8176512177201586, 0.677652869156132, -0.22512436653902138, 0.2013035454614075, -0.21054101354588137, 0.31889463667004814, -0.1864096748428744, -0.5415562468129854, 0.3072855887039624, -0.4981605974794945, -0.40277481366445383, -0.008143591342755163, 0.129008728579238, 0.0511511711133876, -0.05823362486678049, -0.3426475419894246, -0.1266351371509198, 0.16635883436407092, 0.020869347715237102, 0.23449718667102876, -0.4693451543449362, -0.474234003652534, 0.31269810193173553, -0.44421810760942204, 0.3525926425526843, 0.21944810313886137, 0.22450563310835783, 0.1116063450112496, -0.4193005097115841, 0.16699259479742382, 0.205053730732426, -0.5441580429606051, 0.24116461188153296, -0.09996659179120725, -0.46553549089788737, -0.12341996829518986, -0.6459892625746568, -0.1822523524917187, 0.3462574560798293, -0.3722266041470155, 0.1323121303616743, -0.16034748222068082, 0.19675600722855663, -0.20160669403381773, -0.0035643987966246515, -0.2801979759533486, 0.64013893821795, -0.3996982478978147, -0.3138136840634379, 0.16008855537122554, -0.06921753135488129, 0.35679133716655315, -0.9170235309309878, 0.36487649737219435, -0.37415177602986494, 0.024101343199300088, -0.5357829670058332, -0.2724601152182902, -0.0384052820751759, -0.04045319882358911, 0.02358487007868727, -0.5962246048506735, -0.3433469339202841, 0.6377260115949419, 0.6269490246004985, -0.6389125443696233, -0.8156618406019004, 0.07310870887063739, -0.12365256327553362],
    val_dict["Morelia tracyae"] = [0.3299153327872225, 0.03154613266271997, -0.4605622025648676, -0.27894497662396994, -0.5221366182308549, -0.5581762613691341, 0.23009026763990842, 0.38436487731312374, -0.43548292065464134, -0.3641472047872083, 0.4017620957701974, 0.05396197966112505, 0.2889347253019653, -0.09889740462902737, -0.48661411203904104, -0.2833855393612847, 0.47033869570956577, -0.15722891885978832, -0.2741107554079465, -0.03323951295986177, 0.4965629585319527, -0.21146994268260713, 0.2499822681046352, -0.37620214118001566, 0.9244422160295285, -0.18926877898388583, 0.26998850487353804, -0.03489424968850896, 0.5649837532487195, -0.15872594430464051, -0.029134262251640758, -0.35356665432389583, 0.25842660450126265, -0.19916932280173802, -0.07265890108900613, 0.4915378818900439, 0.04191425210993566, -0.45602811826758116, 0.054696437852233465, -0.5195962958143968, -0.8153651341448813, -0.07259550665308546, 0.3128404553240919, 0.44773426503860797, -0.9790417136954974, -0.8052175762336544, -0.6567227183778757, 0.6098604012013978, -9.846029976655069E-4, 0.539311848125118, -0.30603845802732443, 0.15553461954785922, -0.3106349409175177, 0.5717096497526701, 0.029278292857297966, 0.18506562743253774, 0.4783382451244484, 0.5519044216800233, 0.28607712105766675, -0.14367931271811155, 0.5690481613000609, 0.5730600698644394, 0.9205817217574831, 0.7483806281461188, -0.16524882731071433, -0.4777781938795892, -0.006042118376374448, 0.1079919934627685, -0.5297231932570604, 0.022983654649559628, 0.6121104433358946, 0.0485785556942728, 0.4951105545020358, 0.251892457496696, 0.2898524991416508, 0.6539199399138708, 0.6781614514007518, -0.22313155053626926, 0.026250006662555, 0.3381677103557192, -0.04100230631751034, 0.11824599589843, 0.0692325951343957, -1.1611330101760586, 0.4758348519438633, -0.09595266449668896, 0.002951541266410257, -0.04057618370669977, 0.03691190324598015, 0.11992458480646678, 0.539919352516118, 0.36389833099862723, 0.13720790580277384, -0.531612520584613, 0.3856972419201484, -0.25918008787677627, -0.3654176530232903, 0.0222635292188525, -0.3516284620466121, 0.3160044386023662],
    val_dict["Morelia amethistina"] = [0.4451272188943145, 0.45776466182930675, -0.15787064619749938, -0.4336488509933112, -0.9825019225568242, -0.5588417770710905, -0.15971014930800365, -0.06856502153564449, -0.7646821927858499, -0.002986595653355642, 0.49430673766638017, -0.10024829537737812, 0.47252526601043077, -0.018437070682932544, -0.8873355801370083, -0.4321424897481603, 0.49353480535199123, -0.13687496803884802, -0.3881252035956667, 0.047391542969368525, -0.002556025912033584, -0.4773497015850155, 0.6055442129525632, -0.7767512665470901, -0.0013349642505478196, -0.05022645510824442, 0.1308124300171305, -0.11627019076788381, 0.6250733771980926, 0.29193670748270495, 0.09070720706526832, -0.2600624305333218, 0.16741290769249872, -0.3418726157009966, 0.2287968951002014, 0.5817554048782303, 0.016401527281726058, -0.189050068779187, 0.25884225809107997, 0.0257874794687206, -0.5768821460248578, -0.05756178219568807, 0.3953447429625403, -0.2744707604063407, -0.5159863650302206, -0.5893289091467807, -0.31889388432344246, 0.3319342962418501, -0.3836512396099837, 0.5905078610758661, -0.20949598768843133, 0.6076879018658099, -0.2446463695411369, 0.17779134797920626, -0.008089354058430642, -0.10298896382928278, -0.26176281671076795, 0.9651238508483337, 0.603363597873782, -0.2553329927771758, -0.05170879984948916, 0.2923769408200074, 0.5786696899123153, 0.2853310716172128, -0.1701548471233669, -0.8159787975083026, -0.080759451548423, 0.2639301325559496, -0.6898144192422688, 0.19264118155970644, 0.5659900621787252, 0.04881051074627937, 0.5824020675895885, 0.36159025368471037, 0.43909587188864435, 0.14138523856087373, 0.24013838582264846, 0.032843036238126386, 0.36888663459951165, -0.07855661428613231, -0.46433691419747725, -0.16270583082160242, -0.3310065300221858, -0.8356864141942214, -0.0318697143856542, -0.2991133410960931, 0.18502469848042885, 0.2200494926738405, -0.1285210112016325, 0.1444684206750289, 0.5500618257076536, 0.5482891855256047, -0.036576512744920825, -0.6965176386978081, -0.10210158903890787, -0.10385474655935036, -0.2183014634355553, 0.20663437469213175, -0.5835154749433265, 0.21005564286189132],
    val_dict["Morelia nauta"] = [0.35531045339679185, 0.451472484049604, -0.155546931095427, -0.6606018229983057, -0.7577085360552097, -0.3201080213321279, 0.35516083831884065, -0.09322790432469552, -0.5310005016008591, 0.05056088699748315, 0.36049251958102746, -0.27387019669602963, 0.5316176510687503, -0.015237105233307395, -0.4804000425668195, 0.26367260470452486, 0.3536920759912896, 0.2687555496404111, -0.021726616721537972, -0.1709651220922535, 0.08509509058046962, -0.13040404752032375, 0.5919969533288145, -0.48710642842280905, 0.2136492378863657, -0.1005901292760602, 0.32653014291868465, -0.26524712684684476, 0.42020509463787103, -0.08151128165916742, 0.18029479817353522, -0.5511366146417007, 0.5224433145548532, 0.029934109692894832, 0.0953687525951887, 0.11802164858225286, 0.15055722656562617, -0.006462988651472393, 0.19573552478736722, -0.3303632195373853, -0.616421548222383, -0.3512244100687701, 0.18952872451991892, 0.1285257649856241, -0.3499523306235097, -0.7966041933860343, -0.5022259207059379, 0.6231772561193545, -0.28466047859483157, 0.49038221340531807, 0.004160387450733827, 0.12008963778435719, -0.5746521653340917, 0.3463086037893004, -0.2199593657872554, -0.023574179205180828, 0.328671373270493, 0.5283864599263094, 0.07456248250025774, -0.2089134620819154, 0.43749620815207024, 0.24872587858320047, 0.7886818622296559, 0.42465748191677466, -0.26767899978305576, -0.3317228029887408, -0.027873525527065218, 0.5912356871903081, -0.7472544623223, -0.17211760634176076, 0.49585188721161155, 0.39695768137449006, 0.10136646683608978, 0.33068465327339197, 0.39217303224930145, 0.5252474220349376, 0.5981740147371875, 0.11353335408834435, -0.01707607034253361, 0.08690145875422961, -0.4020559935873691, 0.41587627887800677, -0.16309385445203312, -0.5916485164961282, -0.024002574851809844, -0.11990659000533224, -0.061396310758528146, -0.0943582309620423, -0.18965009641348088, 0.35233407331229954, 0.6918684078332206, 0.16523301219200764, -0.20202980049203373, -0.7225350922532985, -0.07373262781088534, -0.2209825049673137, -0.46061794115146887, -0.17305822308584878, -0.2987028056535943, 0.1789548531720234],
    val_dict["Morelia kinghorni"] = [0.3562459592391911, 0.2336460210059938, -0.23493665161288713, -0.5143423578843724, -0.7248847246985228, -0.419034768643285, -0.04921566684487383, 0.06584401269668494, -0.5471876792370683, 0.1138912058041421, 0.4192230811740859, -0.3106102992525145, 0.4312815888237066, 0.028107095293351994, -0.5568209684733971, 0.32683353643332563, 0.09709980621086647, 0.4426699153901462, -0.2206351420867941, -0.086339107948403, -0.14351080424636486, -0.17363964578353722, 0.28440165343282087, -0.47511609118731835, 0.3446474373722793, -0.1190089181397072, 0.29194530892806325, -0.12588252996298668, 0.4281547401010116, 0.00824236415296295, 0.27441514672162903, -0.6554912520268008, 0.6927739751740899, 0.16855796402711035, 0.2858787450035942, -0.04671274103832947, 0.16061170380977524, -0.13082274991231493, 0.210254147963036, -0.2439781221423772, -0.4727091752899368, -0.4923811949175467, 0.24788141078981923, -0.03414004251386266, -0.46873692735068867, -0.8418336033075132, -0.3048639012975817, 0.6110732981539415, -0.3719649474166067, 0.7377708844797686, 7.485408902767332E-5, 0.3269403193191086, -0.554493503312531, 0.5507715742133752, -0.07607418715694259, 0.17602298828127028, 0.32221952331877773, 0.39933325542145076, 0.3366539695903092, -0.40509055165883523, 0.07152926051718292, 0.21124767313832993, 0.6926414335769467, 0.34833423110699935, -0.30713141307242137, -0.32836791950463207, 0.048042938388430184, 0.4247890039246299, -0.6262702556652694, 0.030597240296335815, 0.39807449463095046, 0.42270408433241974, 0.256995810145309, 0.3760351077843726, 0.6146571082233101, 0.6059954033896839, 0.5118304913647002, -0.03721295933388047, -0.08418721477554357, 0.0964484023445107, -0.2918178905000687, 0.32500679652886333, -0.2945655160906069, -0.7125931274046791, -0.19509464562450157, 0.06063262144806111, 0.09173244992998658, -0.13011891211230925, -0.38601262967288724, 0.44619956979642206, 0.4964435179813399, 0.2990038018619171, -0.260949869192768, -0.48654158626423766, 0.003147205481219304, -0.1669558676569498, -0.5169018586666162, -0.07020029457050682, -0.3522090220838987, 0.29257295638418324],
    val_dict["Morelia clastolepis"] = [0.3282418853628491, 0.4446671059708568, -0.250195164318769, -0.4233338026411278, -0.8212680657255227, -0.45136479232484555, 0.2206348633900354, -0.17837258857198143, -0.4994617650895665, -0.05624623143755157, 0.4453145890140769, -0.03520350723657912, 0.5025086593795628, -0.07319523916760932, -0.6754949835959225, 0.1435988836234029, 0.29670023890972164, 0.31343477895786387, -0.04034603572530862, -0.03953944291412391, 0.055307374573465594, -0.14705818268446152, 0.3707219871383655, -0.7099159331828173, 0.24295888141497385, -0.20443908895945972, 0.3148321559357898, -0.10717628696009056, 0.5428720085441888, -0.2825433355310089, 0.15368978468028135, -0.4309984099295538, 0.5605430854756396, -0.04967245120231356, 0.06983996826894803, 0.11140726600164252, 0.04666918610402328, -0.3298876535406888, 0.10355515798309362, -0.25997959315414465, -0.5924882171929291, -0.43594706841878916, 0.2892114289204776, 0.06609036175564158, -0.3915260163707659, -0.7791225064207091, -0.31328658995549985, 0.5742624012817908, -0.45532701296397105, 0.6000142287011858, -0.1554743484421422, 0.2999422380447742, -0.28836915148514286, 0.17178328551442096, -0.11143037561059976, -0.14150009612512718, 0.28154781144446916, 0.38589138192844247, 0.18496388851357276, -0.2608522950153542, 0.2382851019183672, 0.12633549569097113, 0.710078393999058, 0.42386475133097035, -0.2741496344049054, -0.3989463544518551, 0.0011175758170762044, 0.5203539555275268, -0.6483508876337016, 0.11943171291244797, 0.5494679439164927, 0.3752145218878577, 0.34881263931198814, 0.2554856233375256, 0.4566022948126858, 0.5791369297129998, 0.46084035808066764, 0.06377933244033995, -0.18524755155546616, 0.05457536773559496, -0.49613007102840634, 0.18732956953447133, -0.2663839955244267, -0.7248133519088399, -0.2783031344327298, 0.10011982369351817, -0.10432585944646805, -0.08070902101078052, -0.18662700087682094, 0.3314628405333469, 0.6952886121478667, 0.32749400427417863, -0.2168627584712678, -0.5934914037590646, 0.10277176160146692, -0.18874704570532314, -0.3922250345322343, -0.1073848577676714, -0.32646162079303676, 0.15598331776335198],
    val_dict["Morelia boeleni"] = [-0.3951167574323824, 0.13930789461054718, -0.6200937065088064, -0.34927848023744335, -0.32819861178373716, -0.22020424044829784, 0.05298734209931326, 0.6878233696534698, -0.5048922994132848, 0.03739635983636666, 0.4543699160991788, -0.0012664221060070119, 0.28559265930562566, 0.5791305666146839, -0.5008617094003486, 0.14488442047701863, 0.9577579255665494, 0.08903708579668335, 0.3348388043882265, 0.5150247313212906, -0.38988281011685816, 0.06255282549120139, 0.0934117567463047, -0.10097328806287273, 0.5368764896319407, -0.07932377827132314, 0.4268993903121583, -0.18667385171982437, 0.9379743795187836, 0.4901821310968843, -0.4537578619680156, 0.49325981579886125, -0.10219719553572038, -0.355834071862981, -0.4057518004805367, 0.13125107104005587, 0.6456078188266888, -0.6819307459668251, 0.561535918351116, 0.039599983297034194, -0.08000761443624388, -0.0816479070519518, -0.13031242856494513, -0.243160146814355, -0.21560087925981594, -0.9040079240924793, -0.4750396672414062, 0.08213333325683897, -0.0906945038638057, 0.3945067446234555, -0.5003064269809014, 0.5358895899057794, -0.15733367312523527, 0.5741189426757162, -0.28941757947282415, -0.13913375592329075, 0.206415277477198, 0.6135365833360711, 0.1839217177655649, -0.4157957571545637, 0.7450932495457159, 0.2994356417624822, 0.4030825855972029, -0.05994268127039963, 0.15410938117980413, -0.5442643115397771, -0.2899600457263084, 0.10444645624070203, -0.8036186824582134, 0.02446271805079378, 0.7287637244359124, 0.5153724693087773, 0.18931018480021888, 0.24603844242714282, 0.049337312365471064, 0.148178099531118, 0.4904179011855051, 0.23211182635543057, -0.017319563126749364, 0.32555416794495107, -0.7311540650147033, 0.6315587382055559, -0.37003287864226303, -0.5413935436242416, 0.11133922187630363, -0.0011349446247987394, 0.48639303415574503, 0.2762388716342138, -0.7960962129437918, 0.14385630900066207, 0.9470740272080493, 0.3894496081418093, 0.3310193819550525, -0.760374263617934, 0.12433620344919039, -0.39342771284707223, -0.7751460548828208, -0.7651070036387009, -0.07867055167486986, -0.4144106833849736],
    val_dict["Python reticulatus"] = [0.820514542218693, 0.07923602218824337, -0.4683757552600455, -0.3956918048406352, -0.3415521652426199, -0.30680725960698757, -0.9802039180396239, 0.2839025196888503, -0.7012637234222271, -0.5639068617084875, -0.14803985952207277, -0.02581189630344638, 0.5863968285869265, 0.0498817479121265, -1.3215533424312331, -0.6713599797351483, 0.6460889631948993, 0.669132689303595, -0.32841004807487123, 0.2728988371676344, -0.24953429174108294, -0.11196162593613729, -0.47074072505378695, 0.15709590798400172, 0.3030206680809097, 0.0531470417030801, 0.19820158709920085, 0.6058230046605428, 0.6787622750091784, 0.8222950740992301, 0.8729263788375141, 0.007418716714204354, 0.5316176674306201, -0.21320302930515286, -0.30678867758900463, -0.08555217062199427, 0.3542445262494628, -0.013273076660186595, -0.17933435934887734, 0.19857113991305364, -1.0291086547104653, -0.12527421766304325, -0.1395270216025994, 0.025725079445258286, 0.32067826965926494, -0.3893096718374978, -0.6675832547084605, -0.9296363637092755, 0.24581122013814188, 0.26331448418904796, -0.39604998410575765, 0.3912110386468404, -0.3595742850959207, 0.10345480848622304, -0.041243577370919715, 0.0516544816859284, 0.2547509602393403, 0.1931190610096974, -0.2735268987052523, -0.3645967993883974, -0.27012071292680956, 0.33556908360506793, 0.7899268317402146, -0.5819469477589266, -0.5933960891676351, 0.17724483093764012, -0.44617383030766966, -0.14515767896265652, 0.17619455771544196, 0.8769934916429751, 0.4931067080887478, 0.2227654341510159, 0.22369488793630468, -0.08276302228371311, 0.7037772828318114, 0.9733863841875947, 0.05320456286757204, 0.49624871876478777, -0.42884261025928083, 0.9257442495638244, -0.1692357518039948, 0.8991698812556987, 0.08987539511160078, -0.20976484429916165, -0.11468472362818956, 0.1508653369307145, -0.472283578776552, -0.5599977955779174, 0.03118235904849559, 0.1957532450241422, 1.0959940014277934, 0.5833423840653876, 0.2433633313006159, -0.09395655167210078, -0.689653555646675, 0.1419823149540385, -0.16114813981603546, 0.11591781878436885, 0.16276284481228942, 0.6949973404224826],
    val_dict["Python timoriensis"] = [-0.2072629136962617, -0.37862026790897063, -0.5224645133769704, -0.6071937790872736, -0.48248752394057237, -0.7021078096156035, 0.04731300728579585, 0.08828807213096374, -0.5906746995765046, 0.4841446452642575, -0.5387540307697416, 0.43412579044383626, 0.7741265018989738, -0.2296141443811607, -1.1851538418056904, -0.8002333503701269, 0.6099485236275158, 0.6872580696989983, -0.2999413633704826, 0.571188406014784, 0.032103618008636356, 0.0465185518208671, 0.3418828300187857, 0.11995955078701342, 0.27166920881347545, -0.3086948500486003, 0.11244794151115853, 0.43895672608469305, 0.43417658474429127, 0.8547197020889677, 1.1155650919194575, 0.15995008643969386, 0.4256297618504405, -0.007674751574534733, 0.0634826562082601, 0.570974990377628, 0.5398495273237125, -0.09126337813220121, 0.19487119445896697, -0.33306687856952955, -0.709027817817303, -0.036577643207165195, 0.06427153415079484, 0.12003968661470885, -0.32200578498601207, -0.1850410433040751, -0.4795468102879332, -0.5181845611408591, 0.32974248985564514, 0.08918921681013836, -0.3441275496024393, 0.33667816943926937, -0.62058584875621, -0.07294866899583846, -0.49583059553221764, -0.04482645220195264, 0.40689134499135426, -0.13866739037211778, 0.02152797006660999, 0.09271632923631962, -0.13203195128416767, 0.47561102832254, 0.4786106488000729, 0.1492781673647793, -0.2542425924849953, -0.637207903400157, -0.21696112123253344, 0.35713919337000566, 0.1520251511068032, 0.2029156484833991, -0.08272801241371563, 1.0667534189008914, 0.1275764084695354, -0.22545894642117614, 0.7676124735056156, 0.3340189178345957, -0.36642338472881497, -0.29378746145395596, -0.2623936485112499, 1.104751540282142, -0.5080109080467697, 0.8319490733853876, -0.3701032474072984, -1.3812291196778679, -0.4655246268034487, -0.16634875605783866, -0.2238823529582014, -0.4690865444640213, -0.10310357572273698, 0.39802845261494463, 1.0727933965699175, 0.5485219909109751, 0.0455067825890567, 9.874132613632569E-5, -0.3381112083104875, 0.15690966660215763, -0.11331346733984804, 0.05743983490477411, -0.0026095795713535763, 0.8004625278729537],
    val_dict["Morelia oenpelliensis"] = [0.37068195831250994, 0.831218476334904, 0.3889442083394443, -0.3579988404968421, -0.6848056358544848, -0.2850395580148955, 0.08421826212246841, 0.25989640369349787, -0.7779069242965754, -0.2580624084716267, 0.48205719296360566, -0.13539596231296616, 0.7470721536631983, 0.4794370214938752, -0.5505084903739808, -0.34003037020842963, 0.22008399318330302, 0.11509939122438446, -0.6258594795446504, 0.5506595396282344, 0.09220196987847505, -0.34938065598950585, -0.1798981364245853, 0.3493987109417555, -0.15148265585921386, 0.34758589131794343, 0.6309940778180112, -0.06189846636782953, 0.2240734251643709, 0.18004907237866177, 0.6970550471249288, -0.19788204326589368, 0.26970460664199536, -0.012377081101556608, -0.44329636900538527, -0.1472410666614234, 0.27359907422784613, -0.17782172640500638, 0.18321306518474814, 0.01226435629216957, -0.4915967536367546, -0.24087780912813264, 0.05707590677539569, -0.0879041576808651, -0.23308150490181684, -0.565252402082032, -0.5931287769433852, -0.2433640048481378, 0.2212496633247068, 0.39433054872099815, -0.4700907068188156, 1.0965291074832937, -0.23936111379962322, -0.12940655342987128, -0.37921826933391484, 0.21342761637911314, 0.5719154169889299, 0.35968083669104545, 0.48095810234487635, -0.002548820731696627, -0.19297520258114687, -0.15751038821935617, 0.5295979402964204, 0.17194954247055041, -0.47312413966827266, -0.21054186318406043, -0.17636185481016967, 0.4163652484241288, -0.5804390786178926, -0.6388500128216074, 0.7221405674788139, -0.24387422059263944, -0.04961431406801298, 0.3777585629474475, 0.7303873503038565, 0.4413052250414945, -0.5211567970588192, 0.24184108454337702, 0.6701139610490874, 0.315534442463492, -0.14939246272943368, -0.2713880304526892, -0.5135423843026448, -0.6431268291238781, 0.68298528161508, -0.001597808914985854, -0.25017275601972355, -0.2910093843387781, 0.2536590703657681, -0.09204564365659343, 1.0483353205972588, 0.24944066228154005, 0.016360336695804764, -0.8024996425904287, -0.5317902980132994, 0.43130012405904145, 0.07805922539208776, 1.0110265291645576E-4, 0.4107934755364231, 0.5674938751585901],
    val_dict["Morelia viridis"] = [-0.6011884467134684, 0.38481082693488455, 0.29003184665399206, -0.5347410332437568, -0.7043331112499467, -0.9969001038110967, -0.4483916152305475, -0.07646085925825924, -0.21015790342010493, -0.14395326326829105, 0.54161419651437, 0.7919312719673424, 0.5029092026413968, -0.30248144746932265, 0.23209436219395585, -0.9248536250095305, 0.16464376752525373, 0.5496878906732034, -0.2392819717035675, -0.08328176625780126, 0.09506174524266767, -0.2858757106568467, -0.19855354781767165, 0.3782282278529056, 0.24441315386182105, 0.2664312153832061, -0.03821460998768896, -0.47956085391647163, 0.5788223204306115, 0.01993891194985281, 0.8741352800983941, 0.6481195294940938, -0.0125103681784123, -0.8196771850234815, 0.09640069556985792, 0.7355689360782615, 0.38196303113180297, 0.3217663629469017, 0.4643980911588187, -0.05943844299696308, -0.40542818015745513, 0.08023380139961873, 0.6227178852360619, 0.502051861123307, -0.8210956088979137, -1.1106435587124355, -0.5273242615945315, 0.6161319744089666, -0.23907255985717973, 0.8942634641337405, -0.08286756256720401, 0.4725537992465066, 0.33544510362064667, 0.41682965512618514, 0.1101023485269296, 0.3245001315905657, 0.3110457640859807, 0.40095109652521127, 0.226346308367625, -0.42083621709739005, 0.045170763871919345, 0.4048650578403752, 0.8521280151187343, 0.16266334741384392, 0.945202372423849, -0.2670849223897458, -0.002402177960328253, 0.011999587997459205, -0.8368437936773168, -0.03478429208726101, 0.36316289906596105, -0.07783151180101339, -0.2376396422226325, 0.5490770919945687, 0.8759097037425082, -0.12116789679558458, -0.17740614412711558, -0.11860307195289502, 0.3823636703456419, 0.7919780961832283, -0.2701915581705872, -0.13505897411430767, 0.5949139201883813, 0.052158620817918, 0.2654106816416792, 0.45351721114327787, 0.39751398509873925, 0.13099623953483056, 0.02567723435160224, -0.11455639185402272, 0.8428777385524379, 0.7577979254299139, -0.42442084695118065, -0.8058483875049753, 0.2552973952894207, -0.7072952368724428, 0.23552563546295355, -0.19890265551460776, 0.11195456761766401, -0.07658268750439973],
    val_dict["Morelia carinata"] = [-0.3649083410547622, -0.1224514734129345, 0.23099883433455753, -0.6709299647819214, -0.42275097337573353, -0.1218499385752726, 0.09900838167970373, 0.22388317832937943, -0.5263888350721391, 0.2582041784218332, 0.30424289509521374, -0.31863034625662956, 0.554318285213588, -0.1901817822059116, 0.146374965214425, -0.5616158615907442, 0.5864050701527983, 0.43381227586572124, -0.5368627138082377, -0.16985674365339173, -0.9271361016923181, -0.21318859169964444, 0.07921227784642465, -0.1908264327088383, 0.62337793907538, -0.26159288880728426, -0.0911035615087726, -0.6744655153528043, 0.4266618534702292, -0.267742644224292, 0.5122560011510309, -0.09474238258371721, 0.16399573212652807, -0.3475564607170596, -0.4059163145055923, 0.41016485093892713, -0.08080096617087582, -0.18321546400481326, 0.15569016052002002, -0.11684131557792235, -0.6310114761877221, 0.40199615384949516, 0.43115404630206283, -0.13494196492935037, -0.23016972603118016, -0.9316823615687617, -0.1974019010569864, -0.2757677041772189, 0.19313576710914096, 0.1749025733677383, -0.9132536689015347, 0.6370965076806904, -0.3942117444699949, -0.13513784015643104, 0.05731132018482023, 0.06328579314671118, 0.10812868406909587, 0.09025995111425554, 0.5505179168977603, -0.15677415985166443, -0.11422106736427853, 0.3640643119967963, 0.7784642249081726, 0.023460618191408697, 0.1998627981398275, -0.16626082738953638, -0.17563270707926104, -0.07112714582850523, -0.7368214544637348, -0.5131586560594237, 0.7165465310785271, -0.0507749073223871, -0.7200163523584532, 0.8232558928390263, 0.40887094502939214, 0.5508650726777168, -0.07524414867343518, -0.20906440535874696, -0.281854282509165, 0.32279009086782334, 0.31778638097471407, 0.0469823621087604, 0.1807664163729606, -0.4521617448739612, 0.06511217087741925, 0.13419355561837673, -0.5241663842955001, -0.012546813935759044, 0.10680070633183916, 0.5575385308632705, 0.934804235340205, 0.6760471906963296, 0.508953417785927, -0.47327784672903345, -0.09908610385573716, 0.06695245438886371, 0.057132128079106156, 0.08805307370120766, -0.02309422448656709, -0.3757025133066144],
    val_dict["Morelia spilota"] = [0.048068694592795023, 0.163862619209744, 0.3168845382739038, -0.71162227106853, -0.32168575648835995, -0.5933762754801118, -0.04489711893781856, 0.13841229938045302, -0.7048335250940103, 0.13737497954157576, 0.7055957931225807, 0.25964864209093685, 0.4012451138326582, 0.03838591486465118, -0.065816771939216, -0.3978389182225909, 0.4141982413098622, 0.5668273601478439, -0.5474122315708587, 0.2243111502486576, -0.3449846185485714, -0.12352006670097553, -0.28247549581734444, -0.1880409334976428, -0.5178512713622349, 0.4583874759348794, 0.3159978309955572, -0.47494903248021514, 0.6223123451673267, -0.3394833997264299, 0.21700097161602586, -0.3027438711605978, -0.1257383952071134, -0.4412477797295088, -0.3930052653132274, 0.7415109867475027, -0.336900378299817, -0.06924963947106877, 0.10163262256054766, -0.15585928253400627, -0.5615950118806344, 0.23365858953504356, 0.16935508744549702, 0.10696644771030883, -0.5648881910332691, -0.35011994245446487, -0.12759210506475002, 0.4055275754419289, 0.654084630410117, 0.5395250564029039, -0.2388094776221183, 0.33243462224762965, -0.3321631646883555, 0.43254914262082167, 0.2260820264428874, 0.7725768642194972, 0.6337833779469129, 0.42005259036190623, 0.2628688713906745, -0.22045479108949678, -0.6269997235516872, 0.16706206836161608, 0.8638175094017796, -0.3265621628441466, 0.49176294995298614, -0.7215658468188233, -0.08284041661586956, -0.37210385622768366, -0.4489966813938753, 0.04345497779418356, 0.41877364802783856, -0.010254981750056291, 0.19973498660908112, 0.39960229046533313, 8.153160638926377E-4, 0.19513463382757967, 0.014002145117072165, 0.23692553508720499, 0.3030770053459718, 0.5286471078035369, -0.44855145196760937, -0.10090488967379903, 0.037814195584843585, -0.832213776394642, -0.05495540220310029, 0.13043370758097486, 0.0443705733919706, -0.02705098479065121, -0.12231030185401268, -0.12834090791393338, 0.7462477836322738, 0.0148560720963678, 0.3661012325274382, -0.6846024144604729, 0.02660903835409083, 0.2060483838317743, 0.36546868906824626, 0.0791783846005098, -0.01243245687489497, 0.2217730164825047],
    val_dict["Morelia bredli"] = [0.07731487683535951, -0.0034099431477830067, 0.604093513230699, -0.8956849505396182, -0.7174619789634562, -0.48477461381455506, 0.23256350558709227, 0.019682521408267012, -0.7670977631086829, -0.052576503888626815, 0.6257447594882667, 0.8948282430727843, 0.23733452884761555, -0.13646304164353212, 0.033020179856460774, -0.3827888641718163, 0.4899032144284424, 0.2618927204841119, -0.39971646077247064, 0.20016667506915792, -0.2597822738075245, 0.21861936668897802, -0.1252677724367359, -0.4849960031134408, -0.4068585579602465, 0.1415360421068353, 0.38796082862391157, -0.48564664115242706, 0.4115066750116939, -0.3663106961266804, 0.3203746578307555, -0.17929794731079982, 0.21205250765200628, -0.6065945270586861, -0.4928461644651927, 0.9762754318584193, -0.42805414341705855, 0.0655696521190545, -0.1382099593968963, 0.02183513257714262, -0.8102585261004323, 0.3923447052049781, 0.5107296552779939, -0.1503152912536372, -0.45954711285832817, -0.4112273980958088, -0.2611748178922172, 0.4564915358026633, 0.8537047364993811, 0.3129655681662914, -0.34752605784324364, 0.19679627320238452, -0.032077935465145854, 0.6675909114668519, 0.28344975922537563, 0.7528329711622544, 0.4662344739980309, 0.09146082290772486, 0.20918731243646627, -0.12808123509592056, -0.05308923206610222, -0.010263337967043235, 0.7606776280044109, -0.2567965703907963, 0.5362707563927606, -0.4718137713627058, -0.13348268102307886, -0.6505276932823317, -0.27738489546222656, 0.40180185516425304, 0.7015274709193166, 0.2552750016826402, 0.2807446531865529, 0.2548716454047487, -0.13605239640814024, 0.11044270980630819, 0.12337889913057769, 0.05596540666393658, 0.055771877413591076, 0.3314586598441821, -0.3191475320364527, -0.04771010539957657, 0.12998490676235064, -0.7578019914188924, -0.35052789165635195, -0.11329800535390738, 0.20439310663456348, -0.009287010014500663, -0.1380381579133581, 0.36658264785416894, 0.9418556657892142, 0.2157682928189018, 0.09546176217452707, -0.7475382244588441, -0.16819099407276927, 0.23465227003709613, 0.04122773342813807, -0.21527275599190693, 0.20255893420329873, -0.1009102032868531],
    val_dict["Antaresia maculosa"] = [-0.8321016199069953, 0.5740331851509757, 0.6058966882425709, -0.5218260492646062, -0.03978798012603624, -0.49227603967232286, 0.6868972621570989, 0.15701952218848958, -0.18635353485372702, -0.06019917112563311, -0.21432127976588067, 0.3279592272900459, 0.6474602504571034, 0.0038889691079857514, 0.6266080786172936, -1.2434084690209983, 0.9636624130685144, 0.7102610811539656, 0.2749500128468406, 0.08482744453115523, -0.17044491594591366, 0.43739993402382493, -0.035138137385717305, 0.2405021754261976, -0.09172430128143451, -0.15373588063589066, -0.43180570204817564, -0.22917004592799878, 0.3109843259362097, -0.9826468792144798, 0.9548926359125478, 0.41981862148358584, 0.2884539636576724, -0.1877884894721562, -0.4950866927612764, -0.0920884986492656, -0.8049141431246529, 0.6840738739633678, -0.031055280632988957, -0.6605235711954747, -0.2313889776031815, -0.2852524663331515, 0.5712084255609694, 0.46106530429637294, -0.6884747682278657, -0.29835052727881883, -0.22577097408963673, 0.8755318851109384, 0.4933743817040066, 0.848776384869709, -1.0138000115540515, 0.3121510723381837, -0.13464008354699486, -0.12171089740795932, 0.422494903161884, 0.1819223256135022, 0.1390720115245953, 0.01621184046429064, 0.26134112392805103, -0.7415405145243157, 0.1189473477641974, 0.5066795975977587, 0.9237879859548757, -0.18261283039690393, -0.07127700561447547, 0.04169655292355873, 0.2631487812756023, 0.6112761231910822, -0.7852539373485745, 0.21440117134471137, 1.1498287101725302, 0.5642096479095324, -0.27367852955551475, 0.449883535562697, 0.6942484814220015, -0.40199833025840404, 0.05750853287360771, 0.30570663375845814, 0.08255715031767426, 0.931216090930187, -0.4631309133007102, 0.5164156403153031, -0.22789839924946698, -0.6562804823439052, 0.37283848747842363, 0.392013114784017, -0.2801782974087877, 0.48285598018793346, -0.5956267795549056, -0.4907829482936285, 0.8717410321649701, 0.2732819425352813, 0.2917736667290534, -1.0075979152140784, -0.20141315931541554, 0.35767347490406093, 0.022320774383275507, 0.23618123680832892, 0.477531351128489, -0.051389092641904915],
    val_dict["Antaresia childreni"] = [-0.5796433524200388, 0.2613131132913179, -0.3879655076510392, -0.6208893533143592, -0.19215789132071182, -0.5629053868771987, -0.607542254360111, 0.7464003169352552, -0.37113140900681, 0.4615573497823538, 0.3602658939335724, -0.05845710802489425, 0.8922184271728829, -0.08087841415772676, 0.37796584127773697, -1.4675774455255683, 0.7739147142500287, 0.2475800335631762, -0.4262031744883737, 0.628568398666972, -0.1116352776077077, -0.3389915631716427, -0.3575641749292304, -0.08789852212608519, -0.28004616047822595, -0.06899567952562002, -0.5152202112877408, -0.5288762006403549, -0.02571628956762312, -1.1707257776666171, 0.5581212206700906, -0.6765022793654935, -0.042546671170473775, -0.9989109973307242, -0.13844240964063304, 0.6672884783118522, -0.37964190427230304, 0.487656659992695, -0.4582522860761169, 0.46782685261758156, 0.298636326586846, -0.32344903392588753, 0.47032356580687085, 0.4232879750689992, -1.0474160149947103, -0.8382933164391374, -0.7278393438795127, 0.4856417418218345, 0.44230958080955013, 0.7381807513414604, -0.13293747263580266, 0.061422099363614036, 0.27320793230461526, -0.09288658461763483, 0.12542668126346282, -0.4770099287527104, 0.0794504160782413, 0.6099590055884831, -0.3917261956596817, -0.2575689910811042, -0.2580082085340216, 0.4643882612109915, 0.600483654028489, -0.6279504968813056, 0.27665487036445374, -0.5698600848461338, 0.3200401028113309, -0.2066204273384749, -0.7815897332853465, -0.3218082287927039, 0.3446098569505919, -0.0015589058881497753, -0.20470335546622365, 0.6809658769093474, 0.4567301800539355, -0.561274064041086, -0.04875321546814443, -0.615538826895355, -0.35307655110375596, 0.5513825300722629, 0.12190349932982518, 0.052821984093169905, 0.5399775906721933, -1.144263089733151, -0.25427597714893946, 0.05891926811459126, -0.32836221770771873, 1.015080767262193, -0.058720222228657554, 0.06305740963584835, 0.9194972210198715, 0.09986240484775816, -0.5249762683214563, -0.5994068191680474, 0.04632790739234764, -0.23944579592651138, 0.202841819490592, -0.13398733716366973, -0.13154229548279756, -0.17160258796435268],
    val_dict["Antaresia stimsoni"] = [-0.7311268866891303, 0.05189783738036455, -0.12487634615016888, -0.45139784401088456, -0.17009222127622675, -0.7333145744213531, -0.6498285778250642, 0.592062949588059, -0.07599674363717501, 0.4773311788942205, -0.002157271819130002, -0.24075924863987508, 0.8405917879389952, -0.4133524584368803, 0.40248864670715895, -0.9780034523898945, 0.623713954130417, 0.07368582310080962, -0.4048613821685768, 0.4205202961585218, -0.7250505147035659, -0.49659062722199376, -0.23869739444247148, -0.17722560214625527, -0.39006345817681337, 0.031493378761064145, -0.3891582362608474, -0.5331493400210899, 0.3868519620199843, -1.0391964115936245, 0.47002350283574307, -0.24780996525124072, 0.06368264234074945, -1.2511059007808425, -0.04832169514730271, 0.7190093319704692, 0.02888862937167813, 0.2205097300415494, -0.08192987762117418, 0.5335827736462765, 0.10974799367959462, -0.4336764933783772, 0.6929496023837209, 0.5079502760313376, -0.955368172912586, -1.0347426326110036, -0.614276719705225, 0.3348376767872352, 0.2019490692436288, 0.5945731706788362, -0.46477154185464137, 0.06429848481650663, 0.4535613704592157, 0.13780136519221936, 0.2939548579726361, -0.28393470158445233, -0.001841577408981901, 0.5957179558498242, -0.16573317133837862, -0.2550953803457052, 0.15078598926293968, 0.6786778731926959, 0.8036892747098151, -0.6838316571552291, 0.18717404677514837, -0.2295382499108392, 0.5686866817348227, 0.26399434409464595, -0.8425736301836144, -0.38356436931432913, 0.5950424156705437, 0.262668145040298, 0.08798716926615552, 0.45992252616008933, 0.6032561105821319, -0.6004090592403856, -0.18609659442391696, -0.3795284990903761, -0.18127141439300318, 0.4027596848092526, 0.24214750877321115, -0.43908294745349, 0.3484081466142734, -0.8456403723822224, -0.4880079026641227, 0.2210123388065596, -0.3719734211835389, 0.7040403650480351, -0.22964648977900032, 0.05588769976549096, 0.6764938754257839, 0.3996755714023172, 0.03848110455728379, -0.6805078122566885, 0.12198293023617295, 0.44253059451026794, -0.07225601992709559, -0.4169626425448369, 0.06166919110899547, -0.38277034251332875],
    val_dict["Antaresia perthensis"] = [-0.5411032038555792, 0.03540184707851246, 0.4691701265202398, -0.4643035801801999, 0.1477093834719152, -0.74616612706555, -0.33902150886287097, 0.1288332201048099, -0.5205672473822794, 0.39771708141366247, 0.3201564901470565, -0.3253539969659185, 0.6537804775104386, 0.1859841828392057, -0.05710515184834056, -0.3374315611992733, 0.8538323544909576, 0.6452459055390699, -0.6468012599328304, 0.3205341729436997, -0.37690624142723794, -0.2896867308902026, -0.44939174141015314, 0.05842612153034593, 0.24763220777416722, 0.15588865904492782, -0.47244600687781796, -0.7736443212513664, 0.3388805274566974, -0.5929420543523722, 0.23573823343531208, 0.12605716610402073, 0.24319536623294571, -1.4575317357490776, 0.02063730420048443, 0.8552682140691218, -0.09448296137132917, -0.1303503216757613, 0.12592838640066067, 0.2148823908043028, -0.5790735095353414, -1.3201957888422222, -0.0019054753708105387, 0.36718063015136304, -1.0470953210312437, -0.4850298888593329, -0.6876935202462651, 0.102138212221419, 0.6467713464051448, 0.6666137474170334, -0.01846646710819927, 0.5947645456664419, 0.41988390914953355, 0.16131214137107133, 0.4992477775865012, 0.3476518866470112, 1.2807896184068892, 0.4012363746925885, 0.30304336686182626, -0.5022684878652104, -0.1988361337633505, 0.5034629859902598, 1.4637916408893399, -0.4169422759886574, 0.4062093638750266, -0.44835647922145355, 0.9128558638557871, -0.31434511189256914, -0.7525391301413176, 0.28712381510250845, 0.3222199429854229, 0.04178238480234164, -0.05120580834071486, 0.4928539274514989, -0.03964564790841374, -0.008970647535434417, 0.1225847041206296, 0.3043769874012895, -0.8756455274098662, 0.5167888289005171, 0.06216209948535234, -0.5087361908078001, -0.3253280278962525, -1.2781722195748142, -0.38988980606717516, -0.5461141296679947, -0.21956451221702655, 0.5062176355382699, -0.7278346540114577, 0.13409320005956837, 1.116544119077754, 0.6200338417152992, -0.01767667259027562, -0.36789921566571393, -0.39150457298384284, -0.28824617453538237, -0.06167465469608291, 0.06876492487713325, -0.21405868620418148, 0.4064154825472954],
    val_dict["Leiopython albertisii"] = [-0.5998231274204151, 0.6462665450656897, -0.5632813845565258, -0.08523249292307783, -0.4583170997572317, -0.5474828976454414, 0.10697826565637127, 0.6170961858884487, -0.1887331741603404, 0.08984453775728285, 0.7103085736713134, 0.3446912225082618, 0.07952032995553714, -0.13048607721878577, -0.4528821182188478, -0.038147277686385705, 0.6859336477169007, -0.025888649391319685, -0.5263771567164868, 0.05908160161527004, 0.08778287287861974, -1.1007492965807177, 0.05711295535723179, -0.08089570470748608, -0.940315943441733, -0.14585241872592564, 0.5522199369804102, 0.05202469509822183, 0.48240287889635913, 0.18889710874655624, 0.07034197762103417, -0.11895338702769613, -0.23924072452096873, -0.09389150820485692, -0.4680050670285495, 0.927075335755692, -0.33700547091973376, -0.3995191165928913, 0.5274180050879999, 0.10387071488488023, -0.5719915466841509, -0.339777188978158, -0.45105216619652083, -0.34140905874601674, -0.5387201371067307, -1.2923891553002245, -0.09857100778227329, 0.20926570855018978, 0.5283346735686492, 0.420576225871362, -0.2444744577076644, 0.18912429810356968, -0.38330023445001943, 0.12311978654305872, -0.33926755843858214, -0.29123771245674684, -0.5168458836772161, -0.009698111726158934, 0.16475858092078932, 0.24389587381213865, 0.014928310163265818, -0.02537269609035568, 0.4799725631032216, 0.08645074906472404, -0.13288214287169697, -0.4527808283704991, -0.24712976198487085, 0.11128714337516668, -0.48454107661145207, -0.4027937388825462, 0.6727405400153922, -0.13978236357506058, 0.0694493252496344, 0.11311815477293374, 0.005864818463599097, -0.3135858323624224, -0.3518295088567638, 0.7011737039629751, 0.3310916741636899, 0.32011393337481336, -0.5006171826485688, -0.246568536829375, 0.05890491373333037, -0.31593321132584273, -0.07272585092509526, -0.1997344834848739, -0.19537818942690827, -0.08661363942727132, -0.19294510059638986, 0.5818708275244329, 1.275372484269948, -0.3352712532927868, 0.26042341739872926, 0.0822572823645788, 0.3119152016300118, 0.01698004408436009, 0.018461084764971616, 0.2829723662583192, -0.20894539455540767, -0.01655154450363762],
    val_dict["Bothrochilus boa"] = [-0.6535558060210913, 0.4436631514567151, 0.1000625949047606, -0.5762768374173881, 0.12630600542466724, -0.19302434848920758, -0.03731705179523631, 0.2719209388549012, 0.30103011930493007, -0.11619935790132083, 0.4678604897668539, -0.21129305405268206, -0.1132076279773952, 0.25424243391409374, -0.46544109657279475, 0.08717752043595266, -0.1438665118908321, 1.0390536482788981, 0.017573724427164095, 0.03598383526454038, 0.1922932994126063, -0.6631062927029672, -0.10584500502151134, 0.16572325445975306, -0.1619588312199356, 0.42265710545206725, 0.2887066241554257, -0.20539720915285936, 1.0134673582631666, 0.45247328579618884, 0.5318514829946733, -0.5096928004531419, -0.2961371653475703, -0.6558417321936678, 0.3075009363038467, 0.2699733354304047, -0.024282843534304238, 0.4490027229895941, 0.7210380696828821, -0.1587813051034898, -0.6431826922776538, -0.5148674202376615, -0.4372558472042586, -0.6076362357552729, -0.5051995139887193, -0.4163202979685773, -0.019456157664108895, 0.24648031330218334, 0.19841845205523634, 0.7794447280714473, -0.4103730517379157, 0.7433461388047989, -0.23708415167124067, 0.46540848231901477, -0.4521117920563468, -0.07943699988234537, 0.3353563274313148, -0.19892070250972432, -0.23663958550737613, -0.19731409777240005, -0.46536376332802193, -0.3465383867575613, 0.3454284695113518, -0.0156619166471771, -0.32024335193274556, -0.8402474831496303, -0.23075926584579734, -0.5154890345539422, -0.980408789136324, -0.27361885741085185, 0.4696194141600995, 0.08536549533367487, 0.3560090243928915, 0.36219511605216004, -0.05843425523065704, -0.12046673900404087, -0.3480732669105509, 0.569825496171955, 0.2711090817781534, 0.40470400369436516, 0.20318706209663145, -0.0036033313552955096, 0.6369615760442426, -1.0685110493392047, 0.06230235592455257, -0.023811766587732208, 0.42025654420281777, -0.0214084321545287, -0.4049815483911212, 0.6311212537067125, 0.821577473807087, 0.28017554778445675, -0.08377959306788209, -0.4598298825925746, -0.08383192657451564, -0.00987802113863132, -9.673427017853134E-4, 0.27757638038935345, 0.5499996240854729, -0.17151882093050902],
    val_dict["Liasis olivaceus"] = [0.28615756297845407, 0.14909330295801365, -0.4146212989080808, -0.3700505175875797, -0.23776619844956015, -1.017293751704718, -0.06614490554113801, -0.05080352413090844, -0.7556573867595833, -0.06849251482474128, 0.3482462535881704, 0.5273112323130497, 0.33015645338808486, 0.0762818658042201, -0.08418480025693964, -0.5084065450644819, 0.12869736992135478, 0.8423788266647291, 0.07644204274215716, -0.0377116125717617, 0.2207747690040015, -0.2963486097611801, -0.07087260182535145, 0.26691539175888573, 0.07950871249605218, -0.3336583453518465, -0.24313304107415062, -0.6929737143626966, 0.31838005670130715, 0.019313416453130533, 0.36969860157432843, 0.027542300775651063, -0.03440667565903921, -0.32042319734498304, -0.31249222715762515, 0.8146354652586358, -0.37100484063495354, -0.009442042443446957, 0.24861520145504967, -0.07577174515183581, -0.0969735864919965, 0.06506687263709293, -0.16371097950487987, 0.07964321125562951, -0.20570970109508813, -0.7435662209888859, -0.5176204594788005, 0.40965600352038545, 0.020800612725101947, 0.17495629502929194, -0.5091315436719555, 0.10375956549801224, -0.32534303879154347, -0.017621719546406622, -0.2164687463965116, 0.050691247360273084, -0.076329289009754, 0.15121712937228274, 0.23771151207286578, 0.022036905496738743, 0.45737387654736883, -0.10522518713636247, 0.32685810402659465, -0.3002423262469815, -0.07139219014762802, -0.7441142157923473, 0.4235454960480864, 0.05853314994062192, -0.37168030715366507, -0.10340588373436499, 0.7246535807447999, 0.04317389585744236, 0.04967288081788218, -0.13710249552328138, 0.5022193653210615, 0.07638194040519469, -0.10519929537450345, 0.7047956509826943, 0.2521064617029982, 0.4818956599018263, 0.036263474376034804, -0.11229381273483913, -0.09369146000117612, -0.47042169420907526, 0.2146074110770963, -0.7577412266664161, 0.3558870020567058, -0.50922042957762, -0.10228791408930928, 0.28189040067915877, 0.3954909737980713, 0.37880939178777917, -0.06873838171578497, -0.729697174283147, 0.057419051672888285, -0.12857308222824304, 0.14692507540666425, -0.053853938995614234, -0.11841480075001518, -0.40252575091360143],
    val_dict["Liasis mackloti"] = [0.6216849199292218, -0.06942225843775048, 0.16043100666880755, 0.030427950882496627, -0.13852894843762595, -0.35395097140470677, -0.31847899474632746, 0.89223709488102, -0.8759588198647401, 0.024264939971556293, 0.40372079515312587, 0.20589592759960038, 0.45666133832248434, 0.37127298926356234, -0.30449506187467723, -0.6184406409636344, 0.015499117295670328, 0.5921257649980645, -0.18044869149723258, -0.42304389906340467, 0.06885244067989549, -0.050594951737415156, 0.08731823451222845, 0.3556830065625072, -0.28599900399079414, -0.3863681631020801, -0.23496702555882057, -0.7828100932639931, 0.020207247396875072, -0.1352627527877486, 0.3108825023545576, 0.37564316230086403, 0.2959797403394069, -0.7213214546504176, 0.36160521797716993, -0.3837707935563273, -0.18447102998265624, 0.021129203805988525, 0.30885005289658524, 0.3511520861071172, -0.4146105238148439, -0.13098944289732167, 0.07297106260105113, -0.20494132535476295, -0.4903658262916001, -0.5266659752893073, -0.6141604059858126, -0.08511920197879108, 0.29451609676526874, 0.35467546188301613, -0.150519841136015, 0.01114011448311808, -0.45467134774711004, -0.3130013914204719, 0.028965113045880314, 0.32027517482458084, 0.2151485074288313, 0.13532842260899208, -0.016778538834564843, 0.6304130978190696, -0.11771360731725646, -0.4271945051895951, 0.20762463995250877, -0.16264099617702657, -0.11825041570205795, -0.29749495676073345, -0.14110856331158206, 0.7469656007605284, -0.015209686810588319, -0.2859935918135428, 0.6630408083334247, 0.053819165848313644, 0.14426016711314857, 0.9129704378498704, 0.2371756892156998, 0.24689261036492602, -0.20004586064303637, 0.37155165224940845, 0.8257229299856756, 0.31212071249512624, -0.10653534500523808, 0.3251142451000767, 0.44958915761888585, -0.8033013718407247, -0.04817857814941597, -0.35905980917868474, -0.41046731202727454, 0.12788185370673447, 0.220194445779646, 0.2806458272197776, -0.3752963445596108, 0.21069929929534298, 0.1487508968310407, -0.5891096149990097, 0.3294196588913677, 0.504523695598132, -0.06878317039299903, -0.2362226374584391, 0.4172635655127714, -0.007868116194360869],
    val_dict["Liasis fuscus"] = [0.2597922970274713, -0.07721646027275493, -0.1627078040321705, 0.064350228718202, -0.19882392912845287, -0.513920636622111, -0.10793497919082677, 0.9765654182064859, -0.6074627484561713, 0.0695435438417189, 0.6521767837769922, 0.3498949823075982, 0.21717417715848664, 0.22031432010083232, -0.45125051798809945, -0.41269990474895146, 0.4598602215101102, 0.8121712368444279, 0.0035900910161019983, -0.2208743495140616, 0.14748729780738948, 0.06918186096969918, 0.26592379662203436, -0.03801704091314664, -0.2801296371309002, -0.2553426874121216, 0.2470492247993837, -0.7670995999587267, 0.10547993472100992, -0.39571043452882393, 0.034740095143012545, 0.16064889341711158, 0.22326375404290139, -0.636157033999569, 0.15197580763930907, -0.2527352222918598, -0.1808762110535197, -0.2166913068854985, 0.33254459705410033, -0.06573810936316847, -0.33835481283711877, 0.02325161759595243, 0.08014308740816833, -0.6344907564868036, -0.3741413914620149, -0.37785070198807996, -0.47018529657791625, -0.10245251860000862, -0.004186471256085816, 0.5507948646996104, -0.5084368554786524, 0.11354842083885242, -0.40434991885926613, -0.5150009871150228, 0.03632064311782314, 0.35811488635681565, -0.2502615057911402, 0.2583958034752183, -0.039230379210517224, 0.15666846314630883, -0.026361816081769674, -0.35119293042743593, 0.4352423694136638, -0.12022241791301855, -0.2796800023424668, -0.45341070221782737, -0.1602193472019895, 0.8278525095459475, -0.4280788888199612, -0.480789991370767, 0.7413900998124094, 0.30693733421494634, 0.2232978816109952, 0.5349971906909186, 0.25733503156321513, 0.3040326040383715, 0.07240907405448849, 0.5886530675811036, 0.6596911316577585, 0.6227743896810503, -0.08374910338340458, 0.23157883908476778, 0.23019443667103326, -0.5990754159266707, -0.10174760612362382, -0.4911723482270741, -0.33771564714350805, 0.17526009789600164, -0.12240177140111944, 0.5533570441714588, -0.4137169282377443, 0.20141622719774788, 0.05666981675419289, -0.5995239256832884, 0.0531314622281886, 0.6056623441193894, -0.17161072404667613, -0.20252683708982772, 0.5709546977137563, 0.011246583063304977],
    val_dict["Apodora papuana"] = [0.19374963287839972, -0.25747732131930734, -0.2787648989049549, -0.6307526431769149, -0.548627315995751, -0.853724844497691, 0.02240508190010676, 0.3530508495963939, -0.7473756692737074, 0.1828743031411024, 0.6924024776571949, 0.4476850637021738, 0.3725283810608097, -0.17024658880973398, -0.5986082885838072, -0.25604743602138336, 0.34464955046934004, 0.6612521955052824, 0.16757525704512508, -0.14434187168712842, 0.13493455422081269, -0.2883333225146908, -0.19067379572358475, -0.16704200016074666, -0.4434918759247087, -0.2552430262246452, 0.19520300034342783, -0.2885432199080579, 0.18734353123809389, 0.22316846668633422, 0.40465283066825847, 0.39216439349592164, 0.1002224504393337, -0.2668649491509302, -0.3946362253583445, -0.3680926285489763, 0.04570605789391311, -0.02627563460046288, 0.5049883500292311, -0.1513326441424141, -0.5656918525853603, -0.12359714174641959, 0.01828625767586329, 0.16180662342005975, -0.4552827701555713, -0.5305779654643459, -0.3242873316202202, -0.06817236048315711, 0.28598417222146383, -0.33550352541267175, -0.210603339545965, 0.2582801073969163, -0.3537369225150701, 0.6257700482342353, -0.03415769980692686, 0.22315167938310404, 0.0812269359105608, 0.043559205138205034, 0.2754146499354442, 0.28956273632264784, -0.18014130447356472, -0.3632485743056848, 0.8592359804190471, -0.3636623078181991, -0.17469984141561418, -0.2466051551537058, 0.028275997182511767, -0.2952793360912651, -0.6008471693194296, -0.31312132094585826, 0.7456583258752595, -0.07968876061588273, -0.10494181931732634, 0.1945339450853687, 0.4660580306339115, 0.0031241688215104257, -0.12227334099603673, 0.5884563371682967, -0.2426465700349213, 0.43438490389770895, 0.0680589124027326, -0.5909874017999222, 0.053694438903062955, -0.8245090602129982, 0.8410309362532014, 0.4597081411016156, -0.17003060501104308, 0.1591051075834413, -0.34398420590262424, 0.018180552257624993, 0.22703380561195216, 0.23182007984459482, 0.035183410803370424, -0.5219163564878726, -0.04172508003758471, 0.16012569192113424, -0.27120896819263807, 0.26782660055523205, -0.3157178893725667, -0.04240838583910131],
    val_dict["Aspidites ramsayi"] = [0.4030260025710066, -0.10286588506345738, -0.16321328168044422, 0.03748972228927583, -0.724440219218074, -0.5639981063291993, 0.8023904475363992, 0.15345156359735584, -0.42070757593958097, -0.37828868972382423, 0.5041890658669012, 1.020123057746476, 0.29906384867871544, 0.5269083918838595, -0.0033060491380123708, -0.08381890757712629, 0.12631472451432413, -0.12356798063260735, -0.19501439596519768, 0.5017767111329378, 0.5201847694455622, -0.32853640457820915, 0.3461297697902012, 0.03470118282150897, 0.24448982830294733, 0.10546525963337965, 0.439182947956873, -0.4713190364281659, 0.15921162307534778, 0.11747699621099732, -0.025084668428475676, -0.6766300937955984, 0.08853888344346834, -0.5116684177043874, -0.12683948913816462, 0.35537347282405907, 0.6330366074819131, -0.018017870966062027, 0.2810962292179431, 0.03636032330516391, 0.05588801656601369, -0.24840035916256573, 0.12029593968499086, -0.06848420581275448, 0.00858093651483348, -0.9720799334417111, -0.4113235337604839, 0.055318039332803465, 0.18472213190023684, 0.26564142512483907, -0.368392939933609, 0.25989153314650365, -0.5844785892609106, 0.1928110313789009, 2.993978827068311E-4, 0.06503393294097462, 0.15778793106991845, -0.094499277761285, 0.2859275220232446, -0.6090208465294878, 0.10735463293458378, -0.07086985574189861, 0.5550458360115251, 0.03403833104999435, -0.04691138234262605, -0.04947697993755168, -0.6079430839995335, 0.624085394314729, -0.0997208939560337, 0.32354836570906254, 0.27475245902651396, -0.16235315795541155, -0.13893069802896824, -0.14342272868808786, 0.2937495315529048, -0.43958932681042073, -0.1123933574161803, 0.08978750055188657, 0.430396544782953, -0.3961127571497057, -0.016428320030104304, 0.09793884042983438, 0.45667166746088617, -0.8553640494198054, 0.20134186981454233, 0.24247702461720455, -0.28845491642759574, -0.09557595989605763, -0.13882635558681639, -0.3453978910234963, 0.5677945456339374, 0.5349339632002921, 0.27712086476445, -0.2472472668569199, -0.32194256133639637, 0.3255899923304224, -0.47327729036570787, -0.32754597899135557, -0.1353495815326594, 0.2323470235104575],
    val_dict["Aspidites melanocephalus"] = [0.09085174249950884, -0.2544256976392304, 0.23087115438984218, -0.4867110320227567, -0.5964646909148928, -0.256928507633295, 0.19990799373057083, 0.5938483438825557, -0.8110456396313994, -0.13558643994766534, 0.21885622339198108, 0.9333683671324637, 0.05048549952457433, 0.433083306496001, -0.38518757965522205, -0.18838334562882705, 0.21527333041467872, 0.20036357297483887, -0.573645182207695, 0.2322435662986294, 0.40946950124161907, 0.25750350925500265, 0.2102608901511162, 0.025703445511069692, 0.1442587663163839, 0.0916744935644474, 0.2276555944579371, -0.4125795078916638, 0.5539930870076318, -0.375671728456889, 0.5512312886694494, -0.20174009738702214, -0.20330319467377733, -0.637738821536075, -0.008171077470225108, 0.2018954839110449, 0.35163489764253386, -0.07913697066852624, 0.5715095638602992, -0.14579876566232836, 0.06506510435759755, -0.1398754464240233, -0.4321172167536613, 0.21663800212984255, -0.3711115524054085, -1.224398763569613, -0.5426101533846094, 0.44332531471704983, -0.2540370623146857, 6.629274604203483E-4, -0.4459856198648705, 1.0349320080631808, -0.2411303657587437, -0.18217596279311987, -0.54712096642909, -0.22935241243673654, 0.021131944333535796, 0.22097013415444614, 0.07149508996271824, -0.6328313243900353, 0.22211247404743675, -0.10572742285863088, 0.6486220627734074, 0.39166103060987656, -0.5373985491055496, -0.49835950693075803, -0.43011466835057577, 0.6975635537985513, -0.1194293490116965, 0.41254169966121174, 0.30960202218772953, -0.3244249851390915, -0.3944857746370066, 0.12424400717542776, 1.0251921537931201, -0.34853041268731977, 0.3741381556138912, 0.3769829647964014, 0.024774408989706703, -1.1728992336464867, 0.006410385104130695, -0.30168334200048214, 0.30535033982233684, -0.6780296652320257, 0.12468069581757622, 0.051530481331514305, -0.1240437474918601, -0.059346265344719956, -0.24679489394436627, -0.04437341858147278, 0.3709632106702173, -0.3749699026894046, 0.2273598499579883, -0.7369971441929207, -0.02764979131514897, 0.4304842859795367, -0.1596566781664847, -0.11565277544265701, -0.09409359031747405, 0.608857176570636],
    val_dict["Python brongersmai"] = [-0.2039340729508378, -0.11392713678355426, 0.590206255653644, -0.128117127041372, 0.3089397218920988, 0.20512923893243679, 0.27286181644699026, 0.05540252192650534, -0.15995466055149227, 0.02089286044294634, 0.2784241616613561, 0.1803862862853101, 0.5975087631744883, 0.4746932578021958, 0.21394615758383095, 0.5510791450596301, -0.19898943406421657, 0.10612694164404289, -0.16989871052518096, 0.23119147909872126, 0.0098319686348932, 0.09878497681391352, 0.04146304123227589, -0.7797049666075705, 0.30091691775049045, -0.05774306161080924, 0.20106484727842136, 0.5433289672635966, -0.2318164096999961, -0.038144379790516214, -0.03169176281197516, -0.35032050904377704, -0.31249387396377337, 0.23275824558798025, 0.641631844046534, 0.21874458450310835, 0.23390233320809117, 0.8256336286144799, 0.7138653376235308, -0.048301841289484934, 0.33007100677246265, -0.09725459852426666, 0.07987020604953812, -0.46050185399049953, -0.22949313271687397, 0.5509825485977724, 0.6674270448261223, 0.19740468729110516, 0.21286078610885506, -0.23876852119739583, 0.17434840675154697, 0.3756155010789448, -0.25832593892759187, -0.1451444428342736, 0.7127208390299653, 0.006678601129493939, -0.17535124017388734, 0.07789157299263583, 0.21320320991859154, -0.2250511713608291, 0.3143000821904198, 0.34340126329001464, 0.4936557170162722, 0.3720625318179708, -0.4707675394013388, 0.2933879159561611, 0.07978596111464459, -0.3001929562664045, 0.013414612331979412, 0.7295413867111992, -0.03597007379206921, 0.39216428187690877, 0.20621480740192016, -0.036378933386593586, 0.16003393803321023, -0.4241719298312064, 0.026352332446898382, -0.3996908179755169, 0.3143506733933547, -0.08635207493312852, -0.4197623962804345, -0.16728987843264878, 0.08570559962847475, 0.19110951425995173, 0.27307357227512713, -0.11367701184263194, 0.5642482831166731, 0.2170871127272347, 0.7665512664295213, -0.037315405065321765, -0.16170122250344737, -0.149319469994399, -0.0959909268473001, 0.2974764413773902, 0.01573826230862506, 0.26998595194888747, -0.3675535658575683, -0.05464759909805548, 0.948317922392519, -0.07059719818866884],
    val_dict["Python regius"] = [-0.0326060993540331, -0.6902389365500093, 0.4669134590478878, 0.25595646512200876, 0.6221737480944566, -0.3160452687535183, 0.8331955083563803, -1.1331146213677332, 0.13039264827471111, 0.10749084270522627, -0.36634367500368087, 0.12467970025617078, -0.28994239695837754, -0.30618813003314194, 0.4276292573594561, -0.5645592982163212, -0.052157242772052206, -0.16697784566453192, 0.6380392991309564, 0.4727576629312811, -0.15245369689155025, 0.5094913158518717, -0.3657990580974933, -0.018869881358612072, 0.566644366920235, -0.018663075086019866, 0.04382015017242124, 0.22705014446511754, -0.1990821841267509, -0.06842233943714464, -0.15058997131113125, 0.6033962209179695, 0.33039707784434774, 0.6975762748175097, -0.1255522914738201, 0.34033822559832094, -0.44136275200866365, 0.00541106498033185, -0.3564383455904651, 0.14396689275349314, 0.8097497722336567, 0.28692800506849764, 0.40989425401249197, -0.17677889957037965, 0.4117575307997894, 0.08097479841551322, -0.22359579470168092, 0.8786899858705725, 0.1439204291065142, -0.13872650065095593, 0.2799714940593139, -0.08885876295968488, 0.2663000671771633, 0.7757123268058487, 0.05541467077806049, -0.14066117526569805, -0.44402003888541236, 0.216166897052318, -0.2617845856728444, 0.8966164048121148, -0.49380638392947884, -0.2840389035140871, -0.25990806271855305, -0.2679348408252006, 0.30493220343470806, 0.2775534079916471, -0.5281818263389475, -0.7833757330529801, -0.27478339907762495, -0.3820446619785443, 0.3396928314874575, -0.8855711380040459, 0.4099121918327466, -0.251659912667403, 0.032118865092848045, -0.004384921846509141, -0.17061396656347141, -0.2423370430893241, -0.5055651298783994, -0.8875870976295854, 0.2916328080573296, 0.05131876964149176, 0.3019906074433673, 0.05840754354089001, 0.2367392047226879, -0.27247145272713136, -0.10568398918504983, -1.1505648213962472, -0.3644001645770412, 0.0456021291807535, 0.20389648548472605, -0.05113286489878706, -0.350760083919202, 0.5771471855996521, -0.206285286316378, 0.5354423381946355, 0.11255639279243163, -0.45008565896311775, 0.8776851356896243, -0.3450359163373094],
    return val_dict

def reference_dna_dict():
    dna_dict = {}
    dna_dict['Xenopeltis unicolor'] = '???????????????????????????????????????????????????????AAACT--AAAATTCCCATTT-CCCACATATAT----GGATATT--ACGAAG--AAAAAA---GGACTAAAAAAAG---TCCTCTCTCGACCCCCCCCCTACCCCCCCCC-ACAGTT---AGTACGGGT------TTTCC-ATATATGTAACTCTTATAGATTTGCCTATCAAGGC-ATACTATGTATAATCATACATTAATGGCTTGCCCCATGAATATTAAACAGGAATTTCCCTTTAAATATTTTAGCCTAAAAAAGCCTTCGTACAGAACTTTAATA----CCACATTTCT-CAGTCGTTCAATGAAGCACGGAT-ATAGTA--TTGTT-GATAACCATGACTATCC--ACATCCAACTTGTCTTACAGGATCTTGCTA-TTCACGTGAAATCCTCTATCCTTTCATAGCAGGCATACCATTCGACTTCTCACGTCCATAACTTGTTTTAACCTCCCTAATGGCTTTTTCCAAGGCCGCTGGTTACACCTTCAAGATCATCTCGACGGTCCGGAACCATCCCTCTATACTAGCTTTTTCCAAGGCCTTTGGTCGCACCCTTTATAGTGGTACATTTCTCCTCATGATCTGATCACCTATGCCAGTTCAACCACTGGTTTCCTTTTTTT-CTCTGTACCTTTCACCTGACACCCATATATGCACACA----GTCAATGCCCCACCACTACATTCTAATAATCTTCGGCCTACTACCAGTAGCAACAAACATCTCAACATGATGGAACTTCGGCTCAATATTACTAACTTGTCTCGCCATACAAGTACTAACCGGATTCTTTCTAGCCGTCCACTATACAGCCAACATCAATCTAGCATTCTCATCAATTGTACATATCACACGAGACGTCCCCTACGGCTGAATAATACAAAACATTCATGCAATCGGAGCATCCCTATTCTTCATCTGCATTTACATTCACATCGCACGTGGACTATACTATGGATCCTACCTTAACAAAGAAACCTGAATATCAGGCATCACACTTCTCATTACATTAATAGCCACTGCCTTCTTCGGATATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCCACAGTCATCACAAACCTACTCACCGCAATACCATACCTAGGCACATCACTAACAACATGACTTTGAGGAGGATTTGCCATCAACGATCCAACATTAACCCGATTCTTTGCCCTACACTTCATTCTCCCATTCTTAATCATCTCCTTATCCTCACTTCATGTTATCCTTCTTCACGAAGAAGGATCCAGCAACCCTCTAGGAACAAACCCAGACATCGACAAAATCCCATTCCACCCATATCACACACACAAAGATCTCCTCCTACTAACACTTATAATTACCAGCCTATTTATCGTAGTTTCATTCTTCCCAGACATCTTCAACGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTCACACCACAACACATCAAACCAGAGTGATATTTTCTATTTGCCTACGGCATTCTCCGATCAATCCCAAATAAACTAGGAGGAGCACTAGCCCTAGTAATAGCGATCATAATCTTATTTACCCTACCATTCACACACACAGCCCACCTACGACCCATAACCTTCCGACCATTATCACAACTAATATTTTGAACATTAATTTCAACATTCGCAACCATCACATGAGCAGCCACAAAACCCGTAGAACCCCCATTTACCGCAATTAGCCAGGCAGCCTCAATCACCTACTTCACATATTTTGCCTCAAACCCAATTCTCGGCTGAACCGAAAACAAAAT????????????????????????????????????????ACTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACACTCAACCTAGAGGAGCCTGTCTACTAACCGATAATCCACGATTAACCCGACCACTTCTAGCCC-CCAGTCTATATACCGCCGTCGCCAGCCCACCTTGTGAAAGAAACAAAATGAGCTAAACAGTCATA-CACTAACACGACAGGTCGAGGTGTAACTAATGAAGTGGACCAAGATGGGCTACATTTTCTAACACAGACTATACGAATAAAGTCATGAAACTTACCTTAGAAGGTGGATTTAGCAGTAAGATGAGAACACAACACTCAGCTGAAACCATTGCAATATTAAAGGCAACGCCTGCCCAGTGAGA-TCTCTTTAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCTACATGAGAGTTGAACTGTCTCTTGTAATAAATCAATAAAACTGATCTCCTAGTCCAAAAGCTAGAATACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTGCATTAAACCAATTAACAAGTACTTTCAGTTGGGGCGACTTTGGAACAAAACAAAACTTCCAAACACCATGAGCCTTACC-TCATACCT---AGGCCGACAAGCCATTATAAGACCCAGTCTTACTGATTAATGAACCAAGTTACTCCAGGGATAACAGCGCTATCTTCTTCAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCTAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACG????????????'
    dna_dict['Loxocemus bicolor'] = '????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????ATGCCCCACCACTATGTTTTAACATTATTCGGCCTACTGCCAGTAGCTACTAATATCTCTACCTGATGAAACTTCGGATCCATACTACTAACATGCCTAACATTACAAGTAGCGACTGGCTTTTTTCTAGCTGTCCACTACACAGCAAACATTAACCTAGCATTCTCATCCATCGTACATATTACCCGAGACGTCCCATATGGATG-AATGTACAAAACCTACATGCTATCGGAGCATCCATATTCTTTATCTGCATTTACATTCACATCGCACGCGGACTATACTATGGCTCCTACATGAACAAAGAAACATGAATATCAGGAATTACCCTCCTAATTACATTAATAGCAACCGCCTTCTTTGGCTATGTTCTTCCATGGGGACAAATATCATTCTGAGCCGCAACCGTTATTACAAACCTACTCACTGCTGTACCATACCTAGGCACATCACTCACTACCTGATTATGGGGGGGATTTGCAATTAATGACCCTACTCTAACACGATTCTTTGCACTACACTTTATTCTACCATTCCTAATTATTTCACTGTCATCACTACACATTATTTTACTTCATGAAGAAGGGTCAAGCAATCCACTAGGAACCAACCCAGACATCGATAAAATCCCATTCCATCCGTATCACTCCTACAAAGACTTTCTTATACTAACACTAATAATTATAATACTATTAATTATTATATCCTTCTTACCAAACATTTTTAATGACCCAGACAATTTCTCAAAAGCCAACCCACTAGTAACACCCCAACATATTAAACCAGAGTGATATTTCCTATTCGCCTACGGAATTCTACGATCTATCCCCAACAAACTCGGAGGCGCACTAGCACTAACAATATCCATCATAATTCTACTTACAACCCCATTTACACATACAGCCCTATTCCGACCAATAACATTCCGCCCTCTAGCACAACTAACATTCTGAACACTAATTTCAACATTTACAACCATTACATGAGTAGCCACAAAACCAGTAGAACCCCCATACATCGCCATCAGCCAAGTAACCACAGCCCTCTACTTCACATTCTTCATCTCTATCCCATTAACCGGATGAATAGAAAATAAAATAATAAAAATCTAA???????????????????????????GCTC---GCCAAATTACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACACTCAACCTAGAGGAGCCTGTCTACTAACCGATAATCCACGATTAACCCAGCCACTTCTAGCCA-TCAGTCTATATACCGCCGTCGCCAGCCTACCTTGCAAAAGAAACAAAGCGGGCATAATAGTGTCAGCACTAACACGACAGGTCGAGGTGTAACTAATGAAGTGGACCAAGATGGGCTACATTTTCTAACATAGACTATACGAACAAAGATATGAAATCATCTTTCAAAGGCGGATTTAGCAGTAAGCTGAGAACATAATACTCAACTGAAATCAATGCATTATTAAAGGCAACGCCTGCCCAGTGAGAATCTCTTTAACGGCCGCGGTATCCTAACCGTGCAAAGGTAGCGTAATCATTCGTCTATTAATTGTAGACCCGTATGAAAGGCTATATGAGAGTCGAACTGTCTCTTGTAATAAGTCAATTAAACTGATCTTCTAGTACAAAAGCTAGAATAATGATATAAGACCAGAAGACCCTGTGAAGCTTAAACTACTCTATTT-AAACCTATTAATAGATACTTTCAGTTGGGGCGACTTTGGAACAAAAAGAAACTTCCAAATACTACAAGCTATAAC-TCGTATTAA--AGGTTAACAAACCACTATAAGACCCAGTATTACTGATTATTGAACCAAGTTACTCCAGGGATAACAGCGCTATCTTCTTCAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACGCCCAAATGGTGCAGACGCTATTAAAGGTTCGTTTGTTCAACG????????????'
    dna_dict['Morelia spilota'] = '?GCCACACCCCTCACTTCCTCC-------------------CAACCATAGTCTGTAA-TTTACAGACTATGGT--CCATGCCTTAATATA-AAGCCAAAAATCCATATAATTTACCACAAAATAAAG-----CTCTCTC-TCGGCCCCCCCCCTACCCCCCCCC---AARAA-CATTGGGGAR------ACCGGCACACAAAACCA--TTARAAAACTCTTAACAAACCT--CTCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTTCCCTTTTATTATTTTAGTCTAAAATGGCCTTTGTACAAAATATTCTG----TCCTCATTCTCTTGGTCGTTCTATGCAGCACGAGTT--AACTA-ATCTT-ATTAATCATGGATATTC-TCAAC-CTAAGGGTGTCTCTTAGTCTAGCG-CTTCCCGTGAAATCCTCTATCCTTCCATAGAATGCTAACCATTCGACTTCTCACGTCCATATCATGMTA-ATCCTCCCTACTAGCTCTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCCTTATA-TGGTACATATCACCTCATGTTCTGATCATCTATGTCTAT-CCACCATTGGTAGCTCTCTTTTT-TCTGTACCTTTCATCTGACCACCATATATGCACACACACAGTTAATGCCCCACCACTACATCCTAACCTTATTTGGCCTTCTCCCCGTAGCAACCAATATCTCAACATGATGAAACTTCGGCTCAATACTATTAACATGCCTAGCCCTACAAGTTCTAACCGGCTTCTTCTTAGCTGTCCACTACACAGCAAACATTAACCTGGCATTCTCATCCATCATTCACATTACCCGAGACGTCCCATACGGCTGGATAATACAAAACCTACACGCCATCGGAGCATCTATATTCTTCATTTGCATCTACATCCATATTGCACGTGGATTATACTACGGATCCTATCTCAACAAAGAAACCTGAATATCCGGCATTACCCTACTCATCACACTAATAGCAACCGCCTTCTTCGGTTACGTCCTCCCATGAGGACAAATGTCATTCT?AGCCGCAACTGTAATTACAAACCTACTCACCGCCGTACCCTACCTAGGCACATCTCTAACAACCTGACTATGAGGCGGGTTCGCAATCAATGACCCCACCTTAACACGATTCTTCGCATTACACTTCATCCTACCATTCGCAATTATCTCTCTCTCCTCACTACACATTATTTTACTTCACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCTTACCACACCCACAAAGACCTACTCCTACTAACACTTATAATCCTGTTCTTATTCATTGTCGTCTCATTCCTCCCAGACATCTTTAATGACCCAGACAACTTCTCAAAAGCTAACCCTCTAGTAACACCACAGCACATTAAACCAGAGTGGTACTTCCTATTCGCCTATGGCATTCTACGATCCATCCCCAATAAATTAGGAGGCGCACTAGCCCTAGTAATATCAATCCTAATTCTATTCACAATCCCATTCATACACACAGCCTATCTCCGCCCCATAACCTTCCGCCCCCTGTCACAACTCATATTTTGAACACTAATCTCAACATTCGCCACCATTACATGAGCTGCCACAAAGCCAGTAGAACCCCCATTCATCATTATCAGCCAAGTAACCTCAACACTATACTTCACATTCTTCCTATCCATCCCTATTCTAGGGTGAATAGAAAACAAAATAATAAACATCTCCTCCATAACACAGCAGTACAACC--ACTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGTTTAATCCAACCACTTCTGGCCATCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTATT-ACACTAACACGACAGGTCGAGGTGTAACTCATGAAGTGGACAAAGATGGGCTACATTTTCTAACACAGACCATACAAACAAAGACATGAAATTGCCTTTCGAAGGTGGATTTAGCAGTAAGCTAAGGACACAATACCTAGCTGAAATCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACCTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGACTACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCCTATAATAATTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGTCATACC-TCATATCGC--AGGCCAACAAGCCACTATAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTAMGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAAGSATTAACAGTC??'
    dna_dict['Morelia bredli'] = '?GCCAC-CCCCTCACTTCCT---------------------TAACCATAGTCTGTAA-TTTACAGACTATGGT--CCATGCCTTAATATA-AAACCAAAAATCCATATAATTTACCACAAAATAAAG-----YTYTYTY-TYGGCCCCCCCCCTACMCCCCCCC--AAAGAA-CATTGGGAAA------ACCGGCACACAAAACTA--TTAGAAAACTCTTAACAAACCC--CTCTATGTATAATCTTACATTAATGGTTTGCCTCATGAATATTAAGCAGGAATTTCCCTTTTATTATTTTAGTCTAAAATGGCCTTTGTACAAAACATTCCG----TCCTCATTCTCCTGGTCGTTCTATGCAGCATGAGTT--AACCA-ATCTT-ATTGATCATGGATATTC-TTAAC-CTAAGGGTGTCTCTTATCCTAGCA-CTTCCCGTGAAATCCTCTATCCTTCCATAGAATGCTAACCATTCGACTTCTCACGTCCATATTATGCTA-ATCCTCCCTACTAGCCCTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATATCACCTCATGTTCTGATCACTTATGTCTAT-CCACCATTGGTTGCTCTCTTTTT-TCTGTACCTTTCATCTGACCACCATATATGCACACACACAGTTA????????????????????????????????????????????????????????????????????????TTCGGCTCAATACTATTAACATGCCTAGCCCTGCAAATCCTAACCGGCTTCTTTTTAGCGGTCCACTACACAGCAAACATCAACCTAGCATTCTCATCCATCATCCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAASCTACACGCCATCGGAGCATCCCTATTCTTCATCTGCATCTACATCCATATCGCACGTGGGTTATACTACGGATCCTATCTCAACAAAGAAACCTGAATATCCGGTATTACCCTACTCATCACACTAATAGCAACTGCCTTCTTCGGTTATGTCCTTCCATGAGGACAAATATCATTCTRRGCCGCAACTGTAATTACAAATCTACTCACCGCCGTACCATACCTGGGCACATCACTAACAACCTGACTATGAGGCGGATTCGCAATCAATGACCCCACCTTAACACGATTCTTCGCGCTACACTTCATCCTACCATTCGCAATCATCTCTCTCTCTTCACTACACATTATCTTACTTCACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACACCCACAAAGACCTGCTCCTACTAACGCTCATAATCCTGTTCCTATTCATTATCGTCTCATTCCTCCCAGATATCTTCAATGACCCAGACAACTTCTCAAAAGCTAACCCCTTGGTAACACCACAACACATTAAACCAGAGTGGTACTTCCTATTTGCCTATGGCATTCTACGATCCATCCCCAATAAACTAGGAGGCGCACTAGCCCTAATAATATCGATCCTAATTCTATTCACGATCCCATTCATACACACAGCCTATCTCCGCCCTATAACCTTCCGCCCCCTGTCACAACTTATATTTTGAACACTAATCTCAACATTCGCCACCATTACATGAGCTGCCACAAAACCAGTAGAACCCCCATTCATCATCATTAGCCAAGTAACCTCAACACTATACTTCACATTCTTCCTAACCATCCCAATTCTAGGGTGAATAGAAAACAAAATAATAAACATCTCCTCCATAACACGGCAATCAATC---ATTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCACCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCTGACAGTACT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAAGTGGACTAAGATGGGCTACATTTTCTAACCCAGACCACACGAACAAAGACATGAAACTACCTTTCGAAGGTGGATTTAGCAGTAAGCCAGGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACCTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTTATAATAATTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAATCATACC-TCATACCAC--AGGCCAACAAGCCACTATAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Morelia carinata'] = '?GCCACAACCCTCACTTCCT--------------------AAAACCATAGTCTGTAAA--TACAGACTATGGT--TCTTACCTCAATATA-AAGCCAAAAACCCATATAAAACGC-ACACAATAAAACG---CTCTC-C-TCGGCCCCCCCCCTACCCCCCCC--ATAATAAACATAGGAGAA------ATCAGCACACAAAACTA--CTGAAGATACCCCCTCATCTCT--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTTCCCTTCAAATATTTTAGTCTAAATAAGCCTTCGTACAGAATATTTAG----TCCTCATTTTC-TGGTCGTTCAATGCAACACGGATT--AATGG-ATCTT-ACTAACCATGGCTATCC-TTGAT-CAAGKGGKGTCTYTTAATCTAGTA-CTTCCCGTGAAACCCTCTATCCTTCCATAGAATGCTAACCATTCGACTTYTCACGTCCATATACTGCCA-ATCCTCCCTTATAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATCTTGCCTCATGTTCTGATCACGTATGTTAGT-CCACCACTGGTTTCCCTTTTTTTCTCTGTACCTTTCATCTGACTACCATATATGCATACACACAGTTA???????????????????????????????????????????????????????????????????????CTTCGGCTCGATACTATTAACATGTTTAGCCCTACAAGTATTAACCGGCTTCTTCTTAGCTGTTCACTACACAGCAAACATTAACCTAGCATTCTCATCCATCATTCACATCACCCGAGACGTCCCATACGGCTGAATAATACAAAATCTGCACGCCATCGGAGCATCCATATTCTTCATCTGCATTTACATTCATATTGCACGAGGACTATACTATGGGTCTTACCTCAACAAAGAAACCTGAATATCTGGTATCACCCTGCTAATTATCCTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTACCCTACTTAGGCACATCACTAACAACCTGGCTTTGAGGCGGATTCGCAATCAATGACCCAACTCTAACACGATTCTTCGCACTACACTTCATCCTACCATTCGCTATTATCTCCCTATCCTCACTACACATTATCTTGCTTCACGAAGAAGGTTCTAGCAACCCCTTAGGAACCAACCCGGACATCGACAAAATCCCATTCCACCCCTATCACACCTACAAAGATCTTCTTCTACTAACAGTAATAATCCTATTTTTATTCATTATCGTTTCATTCTTCCCAGACATTTTTAATGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGATACTTTCTATTTGCCTACGGAATTCTACGATCCATCCCCAATAAACTGGGAGGAGCATTAGCCCTAGTAATATCAATTATAATCCTATTCACCATTCCATTTATACACACAGCCCATCTTCGCCCAATAACCTTCCGCCCACTATCMCAACTAATATTTTGAACACTAATCTCAACATTTATCACTATCACATGAGCTGCCACAAAACCAGTAGAACCCCCATTCATTACCATCAGCCAAGCAACTTCAGCGCTATACTTCACGTTCTTCCTAACCACCCCAATTCTAGGATGAGTAGAAAATAAAATARTAAACATTCCCTCCATAACACAGCAGCACAACC--ACTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATCAATCCAACCACTTCTGGCCACCCAGTCTATATACCGCCGTCGCCAGCCCACCCTGCAAAAGAAACAAAGTGAGCCAAATAGTACT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAAGTGGACCAAGATGGGCTACATTTTCTAACACAGACCATACGAACAAAGAAATGAAACCGTCTTTCAAAGGTGGATTTAGCAGTAAGCTAGGGACACAACACCTAACTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCTATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACTTATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCCTATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAATAACATGAGCTATTTC-TCATAATAC--AGGCCAACAAGCCACTAAAAGACCCAGTAATACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCTATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Morelia amethistina'] = '?ACCACACCCCTCACTTCCTC--------------------CAACCATAGTCTGTAA-TTTACAGACTATGGT--CCATGCCTTAATATA-AAGCCAAAAATCCATATAATTTACCACAAAATAAAG-----CTCTCTC-TCGGCCCCCCCCCTACCCCCCCCC--AAAAAAACATTGGGGAA------ACCGGCACACAAAACCA--TTAAAAAACTCTTAACAAACCT--CTCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTTCCCTTTTATTATTTTAGTCTAAAATGGCCTTTGTACAAAATATTCTG----TCCTCATTTTCTTGGTCGTTCTATGCAGCATGAGCT--AACTA-ATCTT-ATTAATCATGGATATTC-TTAAC-CTARGGGTGTCTCTTAGTCTAGCG-CTTCCCGTGAAATCCTCTATCCTTCCATAGAATGCTAACCATTCGACTTCTCACGTCCATATCATGCTA-ATCCTCCCTACTAGCCCTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCCTTATA-TGGTACATATCACCTCATGTTCTGATCATCTATGTCTAC-CCACCATTGGTTGCTCTCTTTTT-TCTGTACCTTTCATCTGACCACCATATATGCACACACACAGTTAATGCCCCACCACTACATCTTAACCTTATTTGGCCTCCTACCGGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGCCTTGCACTACAAGTACTAACTGGCTTCTTCCTAGCCGTACACTACACAGCGAACATTAACCTAGCATTCTCATCCATCATTCACATCACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTGCACGCTATTGGAGCATCCATATTCTTCATCTGCATCTACATTCACATCGCACGAGGACTATACTACGGATCATACCTTAACAAAGAAACCTGAATATCCGGCATTACCCTGCTCATCACACTAATAGCAACCGCTTTCTTCGGATACGTCCTCCCATGGGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTTACCGCCATCCCGTACCTAGGTACATCCCTAACAACCTGACTTTGGGGCGGATTTGCAATCAACGACCCCACCCTAACACGCTTTTTCGCATTACACTTCATTCTACCATTCGCAATCATCTCCTTATCCTCACTACACATCATCCTACTACACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGATATTGACAAAATCCCATTCCACCCATACCACTCCTACAAAGACCTCCTCCTACTAACACTAATAATCTTATTCCTATTCATCATCGTTTCATTCTTCCCTGATATCTTCAACGACCCAGACAACTTCTCAAAAGCCAATCCTCTAGTAACACCACAGCACATTAAACCAGAATGATACTTCCTATTCGCCTACGGCATCCTACGATCCATTCCCAACAAACTCGGGGGCGCACTAGCCCTAGTGATATCAATCATAATCCTATTTACCATCCCATTCATACACACAGCCCACCTCCGCCCCATAACCTTCCGCCCACTTTCACAATTAATATTTTGAACACTAGTCTCAACATTTATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCACCATACATTATAATCAGCCAAATAACCGCAACATTATACTTCACATTCTTCTTATCCATCCCCATCCTAGGATGAATAGAAAACAAAATAATAAATATCCCCTCCCAGCCATAACATGGCAAC--ATCTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAGTAACCGATAACCCACGATTAATCCAACCACTTCTGGCCWCCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAMGTGAGCCAAACAGTACC-ACACTAACACGACAGGTCGAGGTGTAACTTATGAAGTGGATCAAGATGGGCTACATTTTCTAACCCAGACAATACGAACAAAGACWTGAAACCATCTTTCGAAGGTGGATTTAGCAGTAAGCCAGGAACACAATACCCAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATGCCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACCAACCTATT--AAACCCTATAATAGCTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAA-ATGAGCTTTCCC-TCACACCTC--AGGCCAACAAGCCACTATAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTCAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAGTGGCGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Morelia oenpelliensis'] = '?GCCAC-MCCCTCACTTCCT---------------------TAACCATAGTCTGTAA-TTTACAGACTATGGT--CCATGCCTTAATATA-AAACCAGAAATCCATATAATTTACCACCAAATAAAG-----YTYTYTY-TYGGCCCCCCCCCTACCCCCCCCC--AAAGAA-CATTGGGAGA------ACCGGCACACAAAATTA--TTAGAAGACTTTTAACATACCC--CTCTATGGATAATTTTACATTAATGGTTTGCCTCATGRATATTAAGCAGGGAWTTCCCTTTTATTATTTTAGTCTAAAACGGCCCTTGTACCAGACATTCCG----TCCTCATTCTCCTGGTCGTTCTATGCAGCATGAGTT--AACCA-ATCTT-ATTGATCATGGATATTCCTTGAC-CTAAGGGTGTSTCTTATCCTAGCA-CTTCCCGTGAAATCCT-TATCCTTCCATAGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCTA-ATCCTCCCTACTAGCCCTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATATCACCTCATGTTCTGATCACTTATGTCTAT-CCACCATTGGTTGCTCTCTTTTT-TCTGTACCTTTCATCTGACCACCATATATGCACACACACAGTTA????????????????????????????????????????????????????????????????????????TTCGGCTCAATACTATTAACATGCCTAGCCCTACAAGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACAGCAAATATCAACCTAGCATTTTCATCCATTATCCACATCACCCGTGACGTCCCATACGGTTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATTTATATTCACATCGCTCGAGGACTATACTATGGGTCATACCTTAACAAAGAAACCTGAATATCCGGTATCACCCTACTCATTACACTAATAGCAACCGCCTTCTTCGGATATGTTCTTCCATGAGGACAGATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCCGTACCATACCTAGGCACATCACTAACAACATGACTATGAGGCGGATTCGCAATCAATGACCCAACCCTAACCCGATTCTTTGCATTACACTTCATCCTACCATTTGCAATTATCTCTTTATCCTCACTACATATCATCCTACTCCATGAAGAAGGTTCCAGCAATCCATTAGGAACCAATCCTGACATTGACAAAATCCCATTCCACCCATACCACTCCCACAAAGACCTACTCCTATTAACACTAATAACCCTACTCCTATTCATTATTGTCTCATTCTTCCCAGATATTTTCAACGACCCAGACAACTTCTCAAAAGCCAATCCCATAGTAACACCACAACACATCAAACCAGAATGGTACTTCCTATTTGCCTACGGCATCCTACGATCCATCCCAAACAAACTAGGAGGCGCATTAGCCCTAGTAATATCAATTATAATTCTATTCACCGCCCCATTCACACATACAGCCTACCTACGCCCTATAACCTTTCGCCCACTTTCACAACTAATATTCTGAGCACTAGTATCAACATTCATCACCATTACATGAGCCGCCACAAAACCAGTAGAACCACCATTTATCATCATCAGCCAAACAACTGCAACACTATACTTCACATTCTTCTTATCCATCCCCATCACAGGATGAATTGAAAACAAAATAATAAACACCCACTCC-TAACATGGCAACACACCG---TTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCATCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTATT-ACACTAAAACGACAGGTCGAGGTGTAACTTATGAAGTGGACCAAGATGGGCTACATTTTCTAACACAGACCATACGAAAAAAGACATGAAAACGTCTCTCAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGACCTCCCAGTACAAAGGCTGGGATACCCCTATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCCTATAATAGACACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCCATCCC-TCATACCAC--AGGCCAACAAGCCATTACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTCAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGSTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Morelia boeleni'] = '??AAACAACCCTCACTTCCTT--------------------CAACCATAGTCTGGA---TTCCAGACTATGGT--TGTTACCTAAAAAACTAAAGAAAAAATCCATATAAAC----------TAAAAA----CTCTCTC-TCGGCCCCCCCCCTACMCCCCCC---GGGTCAGCACAAAAAC-------ATCAC---------------CCAAAAATCCCCCTTTTT-CC--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATAATAAGCAGGAATTTCCCTTTAAATATTTTAGTCTAAAATAGCCTTT-TACATAAAATTATG----TCCTCATTTCT-TGGTCGTTCAATGCAGCACGGATT--AATAT-ATCTT-ATTGATCATGGATATCC-TTGGT-CTAATGGTGTCTCTTAGTCTAACA-CTTCCCGTGAAATCCTCTATCCTTCCATAGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCCA-ATCCTCCCTATCAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAARATCATCTCAATGGTCCGGAACCACCCCTCCATCCTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATCTCACCTCATGTTCTGATCATCTATGTCAAC-CCACCACTGGTAGCTCTTTTTTTCTCTGTACCTTTCATCTGACTACCATATATGCATACACACAGTTA????????????????????????????????????????????????????????????????????????TTCGGCTCTATACTATTAACATGCTTAGGCTTACAAGTAATAACCGGCTTCTTCCTAGCCGTACACTACACAGCAAACATCAACTTAGCATTCTCATCCATCATCCACATCACCCGAGACGTCCCATACGGCTGAATAATACAAAACTTGCACGCTATCGGAGCATCTATATTCTTCATCTGCATTTACATCCACATCGCACGAGGGTTGTACTACGGATCATACCTTAACAAAGAAACCTGAATATCTGGCATTACCCTACTTATCACATTAATAGCAACTGCCTTCTTTGGATACGTTCTCCCATGAGGACAAATATCATTCTGAGSSGCWRCMGTWATCACAAACCTACTCACTGCCATCCCTTATCTAGGCACATCACTAACAACTTGACTATGAGGCGGATTCGCAATCAATGATCCTACACTAACACGATTTTTCGCACTACACTTCATCCTTCCATTCGCAATCATTTCCTTATCCTCACTACACATCATCCTACTCCACGAAGAAGGTTCCAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCATACCACTCCTACAAAGATCTTCTCCTACTAACATTAATAACCCTGTTCTTATTTATCATCGTCTCATTCTTCCCAGATATTTTTAACGACCCAGACAACTTTTCAAAAGCTAATCCCCTAGTAACACCACAACACATCAAACCTGAGTGATACTTTCTATTCGCCTATGGCATCCTACGATCCATCCCCAACAAACTAGGAGGTGCATTAGCCCTAGTAATATCAATCATGATCCTGTTTACCATCCCGTTTACACATACAGCCCACCTCCGTCCTATAACCTTCCGTCCACTCTCACAACTAATATTCTGAATATTAGTATCAACATTCATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCACCATTCATTACCATCAGTCAAGTAACCTCAACACTTTACTTCACATTCTTCTTATCCATCCCCATCCTAGGATGGATAGAAAACAAAATAATAGACATTCCATCCGTAACACGGCAGTACAATC-----GCTGCCCGCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCATCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTATT-ACACTAACACGACAGGTCGAGGTGTAACTCATGAAGTGGACAAAGATGGGCTACATTTTCTAACACAGATCATACAAACAAAGACATGAAATTGCCTTTCGAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCAYTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACCTATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCACATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCCATTCC-TCATACCAC--AGGCCAACAAGCCACTATAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Morelia viridisS'] = '?A-ATCAACCCTCACTTCCTCC-------------------TAGCCATAGTCTGTAAG-TTACAGACTATGGCT--CATGCCTTAATATATAAACCAAAAACCCATATAAT-CACTGAACAATAAAA-----CTYTYTYCTCGGCCCCCCCCCTACCCCCCCC---GGAAAACCATAAAA---------ATCAGCACATAAATAAA--CCTACTAATCCCATTGCTTCCT---CCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTTCCCTTCAAATATTTTAGCCTAAATTAGCTTCCGTACAAAATATCTAG----CCCTCATTTTC-TGGTCGTTCAATGCAATCGGGGTT--AATAA-ATCTT-ACTAACCATGGATATCC-TTGAT-CAGGTGGTGTCTCTTAATTTAGTA-CTTCCCGTGAAATCCTCTATCCTTCCATAGAATGCTAACCATTCGACTTCTCACGTCCATATTATGTCA-ATCCTCCCTTCTAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATCCTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATCTTGTCTCATGTTCTGATCACCTATGCTAGT-CCACCACTGGTTTCCCTTTTTTTCTCGGTACCTTTCATCTGACTACCATATATGCACACACACAGTAA????????????????????????????????????????????????????????????????????????TTCGGYTCAATACTATTAACATGCCTAGCCCTACAAGTATTAACCGGCTTCTTCCTAGCCGTTCACTACACAGSAAACATTAATCTAGCATTCTCATCCATCATCCACATCTCCCGAGATGTTCCATACGGTTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATCTACATCCATATTGCACGAGGATTATACTATGGATCCTACCTCAACAAAGAAACCTGAATATCCGGTATTACCCTACTCATCACACTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTTACCGCTGTACCCTACCTGGGTACATCACTAACAACCTGATTATGAGGCGGATTCGCAATCAACGACCCAACCCTAACACGATTTTTTGCACTGCACTTCATCCTACCATTCGCAATCATCTCTCTATCCTCACTTCATGTTATCCTACTCCACGAAGAAGGCTCCAGCAATCCACTAGGAACTAATCCAGATATCGATAAAATCCCATTTCATCCATACCACTCCTACAAAGACCTACTCCTACTAACACTAATAATCCTATTCTTATTCATCATCGTTTCATTCTTCCCAGACATTTTTAACGATCCGGACAACTTCTCAAAAGCTAACCCATTAGTAACACCACAACACATCAAACCAGAATGGTATTTCCTATTCGCCTACGGCATTCTACGATCCATCCCCAACAAACTAGGAGGCGCATTAGCCTTAGTAATATCAATCATAATCCTATTTACCATCCCATTTACACACACAGCCTACCTCCGCCCCATAACCTTTCGTCCACTATCACAATTAATATTTTGAACATTGGTTGCAACATTCGCCACTATTACATGGGCTGCCACAAAACCAGTAGAACCCCCATTTATCCTCATTAGCCAAGTAACTTCAACACTATATTTCACATTCTTCTTATCCATCCCAATTTTAGGATGAATGGAAAATAAAATAATAAACATCCCCTCCATAACATGGCAGTACAATC--ACTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAGCCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCACCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAATAGTACT-ACACTAAAACGACAGGTCGAGGTGTAACTTATGAAGTGGACCAAGATGGGCTACATTTTCTAACACAGACTATACGAACAAGGACYTGAAATAGTCTCTCAAAGGTGGATTTAGCAGTAAGCTAAGAACATAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGTATATTCATATAAGACCAGAAGACCCTGTGAAGCTTAAATTAACCTATT--AAATCCCATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAATAACATGAGCTATTCC-TCATTACAC--AGGCCAACAAGCCACTACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Morelia viridisN'] = '?A-CTCAACCCTCACTTCCTTC-------------------CAGCCATAGTCTGTAAA-TTACAGGCTATGGCT--CATACCTTGATATATAAACCAAAAACCCATATAATTCACCACACAACAAAA-----CTCTCTCCTCGGCCCCCCCC-TACCCCCCCCC--GGAAAAACATAGAAGAA------GTCAGCACAATTAAACT--TACTGATAACCCCTTGCTTCCT--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTTCCCTTCAAATATTTTAGCCTAAAATAGCCTTTGTACAAAATACCTTG----TCCTCATTTTC-TGGTCGTYCAATGCAATCGGGTCT--AACAA-ATCTT-ACTAACCATGGATATCC-TTGAT-CAAGTTGTGTCTCTTAATCTAGTAACTTCCCGTGAAATCCTCTATCCT-CCATAGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCCA-ATCCTCCCTTTAAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATCCTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATCTTGTCTCATGTTCTGATCACCTATGTTAGT-CCGCCACTGGTTTCCCTTTTTTTCTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAA????????????????????????????????????????????????????????????????????????TTCGGCTCAATATTATTAACATGTTTAGCCCTACAAGTACTAACCGGCTTCTTCTTAGCCGTCCACTACACAGCAAACATCAACCTAGCATTCTCATCCATTATCCATATCACTCGAGATGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATTTACATCCACATCGCACGAGGACTATACTACGGATCCTACCTCAACAAAGAGACTTGAATATCCGGTATCACCCTACTCATCACATTAATGGCAACCGCCTTCTTTGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCTGCAACCGTAATCACAAACCTACTCACCGCCGTACCCTACCTGGGTACATCATTAACAACCTGGCTATGAGGCGGATTCGCAATCAACGACCCAACCCTCACACGATTCTTTGCACTACATTTCATCTTACCATTCGCAATCATCTCTCTATCCTCACTACACGTTATCCTTCTTCACGAAGAAGGCTCCAGCAACCCACTAGGAACTAACCCAGACATCGATAAAATTCCATTTCACCCATACCACTCCTACAAAGACCTGCTACTATTAACACTAATAATTTTATTCTTACTCATTATCGTTTCATTCTTCCCAGACATCTTTAACGACCCGGATAACTTCTCAAAAGCCAACCCCCTAGTCACACCACAACACATCAAACCAGAGTGGTATTTTCTATTCGCCTATGGCATCCTACGATCCATCCCCAACAAACTAGGAGGCGCATTAGCCCTAGTAATATCAATTATAATCCTATTTTCCACCCCATTCACACACACAGCCTACCTACGTCCAATAACTTTCCGTCCACTATCACAACTAATATTTTGAACATTAATCGCAACATTCGCCACCATTACATGAGCCGCCACAAAACCAGTAGAACCCCCATTCATCATCATTAGCCAAGTAACTTCAACACTATATTTCACATTCTTCTTATCCATCCCAATTCTAGGGTGGATAGAAAACAAAATAATAAACATCT???CCATAACATGGCAGTACAACC--ACTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCATTCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCTAAACAGTACT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAGGTGGACCAAGATGGGCTACATTTTCTAACACAGAMTATACGAACAAAGACATGAAACAGTCTCTCAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGCATATCTATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAATCTTATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAATAACATGAGCCATTCC-TCATCACAC--AGGCCAACAAGCCACTACAAGACCCAGTAAAACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Liasis olivaceus'] = '?GCCACAACCCTCACTTCCCC-----------------ACCTAACCATAGTCTGTAAA-TTACAGACTATGGT--TGATACCTTAATACA-AAGCCGAAACCCCATATAAACAGCACCACAACAAAA----CTCTACTC-TCGGCCCCCCCCCTACMCCCCCCC--ACAAAAACATAGGARAA------ATCAGCACAAACAATC---MCCTAAAATCCCCCCTTAACCC--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTCCCCTTTAAATATTTTAGTCTGAATTAGCCCTTGTACAAAAAATCTTG----TCCTCATTTTC-TGGTCGTTCAATGCAGCACGGATT--AATAG-ATCTT-ATTAACCATGGCTATCC-TTGAT-CTAGTGGTGTCCCATGATCTAGTA-CTTCCCGTGAAATCCTCTATCCTTCCATTGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCCA-ATCCTCCCTAACAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATTTCACCTCATGTTCTGATCAGCTATGTCAAT-CCACCACTGGTAGCTCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAAATGCCCCACCACTACATTCTAACCCTGTTCGGCCTCTTACCAGTAGCAACCAACATTTCAACATGATGAAACTTCGGTTCAATACTACTAACATGCCTAGTCCTACAAGTATTAACCGGTTTCTTCCTAGCTGTCCACTACACAGCAAACATCAATCTAGCATTCTCATCCATCGTTCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAACCTACACGCTATCGGAGCATCTATATTCTTCATTTGCATCTACATCCATATCGCACGAGGTCTATACTACGGATCATACCTTAACAAAGAAACCTGAATATCTGGTATCACCCTACTCATCACACTAATAGCAACCGCTTTCTTCGGATATGTCCTTCCATGGGGACAAATATCATTCTGGGCCGCAACCGTAATCACAAACCTACTCACTGCCGTACCCTATCTAGGCACATCACTAACAACCTGACTTTGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGATTCTTCGCACTACACTTCATCCTACCATTCGCAATCATCTCTTTATCCTCACTACACATCATCCTACTCCACGAAGAAGGATCTAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCCTATCACTCACACAAAGACCTTCTCCTATTAACACTAATAATAATATTCCTATTCATTATCGTATCATTCTTCCCAGATATTTTCAACGACCCAGATAACTTCTCAAAAGCCAACCCCTTAGTAACACCACAACACATTAAACCAGAATGATACTTCCTATTCGCCTACGGCATTCTACGATCTATCCCCAACAAACTTGGAGGCGCATTAGCTCTAGTAGCATCAATCATAATCCTATTCACCACCCCATTCACACACACAGCCAACCTCCGCCCTATAACCTTCCGCCCCCTGTCACAACTAATATTCTGAACATTAGTCTCAACATTCATCACTATTACATGAGCCGCCACAAAACCAGTAGAACCCCCATTCATCGCTATCAGCCAAGTAACCTCAATACTCTATTTCACATTTTTCTTGTCCATCCCCATCCTAGGATGAATAGAAAACAAAATAATAAACACCCCCTCCATAACACGGCAACACAACC--GTTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACCCCTAGCCACTCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTCAT-ACACTAACACGACAGGTCGAGGTGTAACTTATGGGGTGGACCAAGATGGGCTACATTTTCTAACATAGACCATACGAACAAAGACATGAAACCGTCTTTCAAAGGTGGATTTAGCAGTAAGCTAAGAACATAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGTATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTTTCAGTACAAAAGCTGAAATACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTATTAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCTATTCC-TCATAACAC--AGGCCAACAAGCCACTACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Liasis mackloti'] = '?ACCACAACCCTCACTTCCTT--------------------CAGCCATAGTCTGTAA-TTTACAGGCTATGGC--TGATACCTTAATATA-AAACCAAAATCCCATATAAATACCACCACAACAAAG-----CTCTCTC-TCGGCCCCCCCCCTACMCCCCCCC--ACCAAAACATAGAARAA------ATCAGCACAATAAATA---CTARAAGTATTTGCTTCCTTCC--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTTAGCAGGAATTTCCCTTTAAATATTTTAGCCTAAAATAGCCTTTGTACACAAAACTATG----TCCTCATTTCT-TGGTCGTTCAATGCAGCACGGATT--AATAG-ATTTT-AATAACCATGACTATCC-TTGAT-CTAGTGGTGTCCCATGATTTAGTA-CTTCCCGTGAAATCCTCTATCCTTCCATTGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCTA-ATCCTCCCTAACAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATTTCACCTCATGTTCTGATCAGCTATGTTACT-CCACCACTGGTATCCCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTTAATGCCCCACCACTACGTTCTAACCCTGTTTGGTCTCTTACCAGTAGCAACCAATATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTTTAGCCCTACAAGTACTAACCGGATTCTTCCTGGCTGTCCACTACACAGCAAATATTAACCTGGCATTCTCATCCATCGTTCACATCACTCGAGATGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCTATATTCTTCATTTGTATCTACATCCACATCGCCCGAGGCCTATACTACGGATCATACCTTAACAAAGAAACCTGAATATCCGGTATTACCCTGCTTATCACACTAATAGCAACCGCCTTCTTCGGATACGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTTCTCACCGCCGTACCCTACCTAGGCACATCCTTGACAACCTGGCTATGAGGGGGGTTCGCAATCAACGACCCCACCCTAACACGATTCTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCTTTATCCTCACTACACGTCATCCTCCTTCATGAAGAGGGGTCTAGCAACCCACTAGGAACTAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACTCCCACAAAGACCTTCTCCTACTAACACTAATAATAATATTCCTATTCATTATTGTTTCCTTTTTCCCAGACATCTTCAACGACCCAGACAATTTCTCAAAAGCCAACCCTCTAGTAACACCACAACACATTAAACCAGAGTGGTACTTCCTATTCGCCTACGGCATCCTACGATCTATCCCCAACAAACTTGGAGGAGCATTAGCCCTAGTAATATCAATCATAATCTTATTTTGTACCCCATTCACACACACAGCCCACCTCCGCCCCATAACTTTCCGCCCACTATCCCAACTAATATTCTGAACACTAGTCTCAACATTCATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTCATCACCATTAGCCAAGTAACCTCAATTCTATACTTCACATTTTTTCTATCCATCCCTATTCTAGGATGAGTAGAGAACAAAATTATAAACGCCCCCTCCATAACACGGCAACACARC--CATTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACCCCTGGCCATTCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTCAT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAGGTGGACCAAGATGGGCTACATTTTCTAACACAGACCACACGAACAAAGACATGAAACTGTCCTTCAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGGATACCCCTATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTTATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAAAATGAGCCATCCC-TCATACCAA--AGGCCAACAAGCCACCACAAGACCCAGTAAGACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCTATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCCGATGGTGCAGCCGCTATCAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Liasis fuscus'] = '?ACCACAACCCTCACTTCCTC--------------------CAGCCATAGTCTGTAA-TTTACAGGCTATGGC--TGATACCTTAATATA-AAACCAAAATCCCATATAAATACCACCACAACAAAG-----CTCTCTY-TCGGCCCCCCCCCTACCCCCCCCC--ACCAAAACATAGAAGAA------ATCAGCACA-AAATAACA-CTAGAAGTATTACTTCCTTGCC--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTTAGCAGGAATTTCCCTTCAAATATTTTAGCCTAAAATAGCCTTCATACATAAAATTATG----TCCTCATTTCT-TGGTCGTTCAATGCAGCACGGATT--AATGG-ATTTT-AATAACCATGACTATCC-TTGAT-CTAGTGGTGTCCCATGATTTAGTA-CTTCCCGTGAAATCCTCTATCCTTCCATTGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCTA-ATCCTCCCTAACAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATTTCACCTCATGTTCTGATCAGCTATGTTATT-CCACCACTGGTATCCCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAA????????????????????????????????????????????????ACCAATATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTTTAGCCCTACAAGTATTAACCGGATTCTTCCTGGCTGTCCACTATACAGCAAATATTGACCTGGCATTCTCATCCATCATCCACATCACTCGAGACGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCAATATTCTTCATTTGTATCTACATCCACATCGCCCGAGGCCTATACTACGGATCATACCTCAACAAAGAAACCTGAATATCCGGCATCACCCTACTTATCACACTAATAGCAACCGCCTTCTTCGGGTACGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTTCTTACCGCCGTACCCTACCTAGGCACATCCTTGACAACCTGATTATGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGATTCTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCTTTATCCTCACTACACGTCATCCTCCTCCACGAGGAGGGGTCTAGCAACCCACTAGGGACTAACCCAGACATCGACAAAATCCCATTCCACCCTTACCACTCCCACAAAGACCTTCTCCTACTAACACTAATAATAATATCCCTACTCATTATTGTTTCCTTTTTCCCAGACATCTTCAACGACCCAGACAACTTCTCAAAAGCCAACCCCTTGGTAACACCACAACACATTAAACCAGAATGATACTTCCTGTTCGCCTACGGCATCCTACGATCTATTCCCAACAAACTTGGAGGAGCATTAGCTCTAGTAATATCAATCATAATCTTATTTTCTACCCCATTCACACACACAGCCCACCTCCGCCCTATAACTTTCCGCCCACTATCCCAACTAATATTCTGAACACTAATCTCAACATTCATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTCATCATCATCAGCCAAGTAACCTCAATACTATACTTCACATTTTTCCTATCCATCCCCATTCTAGGATGGGTAGAGAACAAAATTATAAACACCCCCTCTATAACACGGCAATACAAC--CATTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACCCCTGGCCACCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTCAT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAGGTGGACTAAGATGGGCTACATTTTCTAACACAGACCACACGAACAAAGACATGAAACTGTCCTTCAAAGGCGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACCCCTATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTTATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAAAATGAGCCATCCC-TCATACCAA--AGGCCAACAAGCCACCACAAGACCCAGTAAAACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCTATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCCGATGGTGCAACCGCTATCAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Liasis albertisii'] = '????GCTCCTCTCACTTCCTC--------------------AGACCACAGTCTGCAA---TGCAGACTGTGGTTTTGTGCCCAGAATATA--AACCAAAAAACCATATAAACAACACCRCGACAAAAAAGA--TCTCTC-TCGGCCCCCCCCMTACMCCCCCCC--AAAAAAACATARAGGAA------ATCAG--------------------TTCATARACT--------CTCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTCCCCTTTAAATATTTTAATCTAAATTAGCCTTCGTACACAAAATTCAG----TCCTCATTCTC-TGGTCGTTCAATGCAGCACGGATT--AATCA-GTCTT-ACTAACCATGGATATCC-TTGAT-CTAGTCGTCTCTCTTAGTCTAACA-CTTCCCGTGAAACCCTCTATCCTTCCACTGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCTA-ATCCTCTCTACTGGCTCTTTCCAAGGSCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACACTTCACCTCATGTTCTGATCAGCTATGTCAAT-CCGCCTTTGGTTTCTCTTTTTTT-CCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTT?ATGCCCCACCACTACATTTTAACCCTATTTGGCCTCCTACCCGTAGCAACCAACATCTCAACATGATGAAACTTTGGTTCAATACTATTAACATGCTTAGCTCTACAGGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACAGCAAACATCAACCTAGCATTTTCATCCATCATCCACATTACCCGAGATGTCCCATTCGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATCTACATTCACATCGCACGGGGGCTCTACTACGGATCATACCTAAACAAAGAAACCTGAATATCCGGCATTACCCTACTCATCACACTGATAGCAACCGCCTTCTTCGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCTGCAACCGTAATCACAAACCTACTAACCGCCGTACCCTACCTAGGCACATCACTAACAACCTGATTATGGGGCGGATTTGCAATCAACGACCCTACCCTAACACGATTCTTCGCACTACACTTCATCCTACCTTTCGCAATCATCTCCTTATCTTCACTACACATTATCCTTCTTCACGAAGAAGGCTCTAGCAACCCACTAGGAACCAATCCAGACATCGACAAAATCCCATTCCACCCCTATCACTCCCACAAAGATCTTCTCCTACTGACACTAATAATACTATCTCTGCTCATCATCGTCTCATTCTTCCCAGACATCTTTAATGACCCAGATAACTTCTCTAAAGCCAACCCATTAGTAACGCCACAACACATCAAACCAGAATGATACTTCCTATTCGCCTATGGTATCCTACGATCCATCCCCAATAAACTAGGAGGCGCACTAGCCCTAGTAATATCAATCATAATTTTATTCACCATCCCGTTCACACACACAGCCCATCTTCGCCCCATAACCTTCCGCCCATTCTCACAACTAATATTTTGAACACTAGTCTCAACATTTATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTTATCGTCATCAGCCAAGTAACTTCAACACTATACTTTACATTCTTCTTACTCATCCCCATTTTGGGCTGAACAGAAAATAAAATAATAAACACCCTCTTCCTAACATGGCAACACAAT--AACTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCATTCAGTCTATATACCGCCGTCGCAAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAATAGTCCT-ACACTAACACGACAGGTCGAGGTGTAACTCATGAAGTGGACCAGGATGGGCTACATTTTCTAACACAGACCACACGAACAAAGATATGAAACTATCTTCCAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATGCCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTYTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATAAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACCAACATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCCCATAATAGCCACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCCATCCC-TCATACCTC--AGGCCAACAAGCCACTATAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Apodora papuana'] = '?GCCACAACCCTCAMTTCCTT--------------------CAGCCACAGTYTGTAA-TTTACAGACTGTGGC--CCATGCCTCAATATA-AAGCCGAAAATCCATATAAATAACACCAAAACAAAG-----CTYTCCC-TYGGCCCCCCCCCTACCCCCCCCC--AAAAAAATATAGAGAA------CTATAGAACAAATAACCA---CCAAGAAGTTCACTATCCCCC--TCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCCGGAATTTCCCTTTAAATATTTTAGTCTAAATATGCCCTTGTACACAAAATTCAG----TCCTCATTTCT-TGGTCGTTCAATGCAGCCAGGAAT--AATCA-ATCTT-ATTAACCATGGATATCC-TTGAT-CTAGTGGTGTCTCTTGGTCTAGTA-CTTCCCGTGAAATCCTCTATCCTTCCATTGAATGCTAACCATTCGACTTCTCACGTCCATGTATTGCCA-ATCCTCCCTAACAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATCCTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACACTTCACCTCATGTTCTGATCAGCTATGTCAAT-CCACCACTGGTAGCTCTCTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAAATGCCCCACCATTACATCCTAACCCTGTTCGGCCTCCTACCAGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGCCTAGCCCTACAAGTATTAACTGGCTTCTTCCTGGCCGTACACTACACAGCAAACATCAACCTAGCATTCTCATCCATCATTCACATCACCCGAGATGTCCCATACGGCTGAATAATACAAAACTTACACGCCATCGGAGCATCCATATTCTTCATCTGTATCTACATCCATATTGCACGGGGCCTATACTACGGATCGTACCTAAATAAAGAAACCTGAATATCTGGCATCACCCTACTCATCACACTAATAGCAACCGCCTTCTTCGGATATGTCCTTCCATGAGGACAAATGTCATTCTGAGCCGCAACTGTAATCACAAATCTGCTCACTGCAGTACCCTACCTGGGTACATCACTAACAACTTGATTATGAGGCGGCTTCGCAATCAATGACCCCACCCTGACACGATTCTTCGCACTACACTTCATCCTACCATTCGAAATCATCTCTTTATCCTCACTACACATCATCCTACTTCATGAAGAAGGCTCTAGTAACCCATTGGGAACTAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACTCTCACAAAGACCTTCTCCTATTAACACTAATAATTCTACTCCTATTCATCACCATATCATTCTTCCCAGATATTTTCAACGACCCAGACAACTTCTCAAAAGCCAACCCACTAGTAACACCACAACACATCAAACCAGAATGATACTTCCTATTCGCCTATGGGATCCTACGATCCATCCCCAACAAACTAGGGGGCGCATTAGCCCTAGTAACATCGATCATAATCCTATTCACCATCCCATTTACACACACAGCTCACCTCCGACCTATAACCTTCCGCCCCCTATCACAACTGATATTCTGAACATTAGTATCAACATTTATCACCATTACATGGGCCGCCACAAAACCAGTAGAACCCCCATTCATCATCATTAGCCAAGCAACCTCATTATTATACTTCACATTCTTCTTATCCTTCCCTATCCTAGGGTGGACAGAAAACAAAATAATAAACACCCCCTCCATAACATGGCAATACAACC-CATTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACCCCTAGCCTCTCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAATAGTTAT-ACACTAACACGACAGGTCGAGGTGTAGCTTATGGGGTGGCCCAAGATGGGCTACATTTTCTAACATAGACCATACGAACAAAGACATGAAACCGTCTTTAAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCGATGCCTGCCCAGTGAGAATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTCATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAATAACATGAGCCATCCC-TCATACCAA--AGGCCAACAAGCCACCACAAGACCCAGTAAGACTGATAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCCGATGGTGCAGCCGSTATCAAAGGTTCGTCTGTTCAACGATTAACASTCCT'
    dna_dict['Bothrochilus boa'] = '??GCACCGCCCTCACTTCCTC--------------------CGACCGCAGTCTGCC----AGCAGGCTGCGGTC-GCATGCCCAAAAACACAAACCAAAAAACCATATAAACAACGCCGCAACAAAAGG----YCYCYC-TCGGCCCCCCCCCTAC-CCCCCCC-ACAAAAAACATAGAGAAA------ATCAG--------------------TTTTCACAC---------CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTTCCCTTCAAATATTTTAATCTAAAATAGCCTTTGTACATAAAATTTGC----CCCTCATTTCT-TGGTCGTTCAATGCAGCATGGATT--AATCA-GTCTT-ATTAACCATGGATATTC-TCAGT-CTAGTTGTGTCTCTTAGCCTAACA-CTTCCCGTGAAATCCTCTATCCTTCCACTGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCTA-ATCCTCCCTACT-GGTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATTTTGTCTCATGTTCTGATCATCTATGCCTAC-CCTCCATTGGTTTGTCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTA?????????????????????????????????????????????????????????????????????????TTTGGCTCAATATTATTAACATGCCTGGCCCTACAAGTACTAACCGGCTTCTTCCTGGCCGTCCACTACACAGCAAACATCAACCTGGCATTCTCATCCATTATTCACATCACCCGAGATGTCCCATATGGCTGAATAATACAAAACCTGCACGCCATCGGAGCATCCATATTCTTCATTTGCGTATACATTCACATCGCACGAGGACTATACTACGGGTCATACCTAAACAAAGAAACCTGAATATCTGGCATTACCCTGCTCATCACACTAATAGCGACCGCCTTCTTTGGATATGTCCTCCCGTGAGGACAGATATCATTCTGAGCCGCAACCGTAATTACAAACCTGTTAACAGCAGTACCCTACCTGGGCACATCACTAACAACCTGGTTGTGAGGCGGGTTCGCAATCAATGACCCCACCCTAACACGATTCTTCGCACTACACTTCATCCTACCATTCGCAATCATCTCCCTATCCTCACTACACATCATCCTACTTCACGAGGAGGGCTCCAGCAACCCACTAGGAACCAACCCAGACATCGACAAAATCCCTTTCCACCCCTACCACTCCCACAAAGACTTTCTTCTTCTAACACTAATAACCCTATCCTTACTCATCATCGTCTTATTCTTCCCAGACATCTTTAACGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGATACTTCCTATTCGCCTACGGCATCCTACGTTCAATCCCCAATAAACTAGGAGGCGCACTAGCCTTAGTAATATCAATCATAATCCTATTCACCATCCCATTCACACACACAGCCCACCTTCGCCCCATAACCTTCCGTCCACTATCACAACTAATATTCTGAACATTAGTGTCAACATTTATCACTATCACATGGGCCGCCACAAAACCAGTAGAACCACCATTTATCACTATCAGCCAAACAACCTCAACACTATATTTTACATTCTTTTTACTTACCCCCATCCTAGGCTGAATAGAAAACAAAATAATAAAAACTCCCTCCCTAACACGGCAACACTATC--ATTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCATCCAGTCTATATACCGCCGTCGAAAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAATAGTTCT-ACACTAACACGACAGGTCGAGGTGTAACTCATGAAGTGGACCAGGATGGGCTACATTTTCTAACACAGACCATACGAACAAAGACATGAAACCGTCTTCCAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATGCTTAGCTGAAAC--TTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATAAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGGATACCAACATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCCCATAATAGCTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCCATTCC-TCATATCCT--AGGCCAACAAGCCACTACAAGACCCAGTAACACTGATAACTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTCAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTRTTCAACGATTAACAGTCCT'
    dna_dict['Antaresia maculosa'] = '?ACCACAACCCTCACTTCCTCC--------------------AGCCATAGTCTGTAAATTTACAGACTATGGC--TGATACCTCAACATA-CAGCCAAAATTCCATATAATAT-CCCCACAACAA-----CTYTYTYTCYTYGGCCCCCCCCCTACCCCCCCC---ATCCAAATATATAAGAA------ATCAGCACAATAAACCT--ACTAGGAATTGCCAATAACTCC--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTCCCCTTCAAATATTTTAGTCTAAAATAGCCTTTGTACAGAATATTTAG----TCCTCATTTCT-TGGTCGTTCAATGCAACACGGATT---ATCAGTTCTT-ACTAACCATGGATATCC-TTGAT-CTAGTGGTGTCTCTTAATCTAGTA-CTTCCCGTGAAATCCTCTATCCTTCCATAGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCCA-ATCCTCCCTATCAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCTTACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTTCATTTTGTCTCATGTTCTGATCACCTATGCTTAT-CCACCACTGGTAGCTCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAA????????????????????????????????????????????????????????????????????????TTCGGCTCAATATTATTAACATGSCTGGSCCTGCAAGTATTAACCGGATTCTTCTTGGCCGTCCATTACACAGCAAATATCAACCTAGCATTTTCATCCATTATTCACATCACCCGAGATGTCCCATACGGCTGAATAATACAAAACCTACACGCCATYGGAGCCTCCATATTCTTCATTTGCATTTATATCCACATTKYACGAGGACTATACTACGGATCTTACCTCAATAAAGAAACCTGAATATCTGGTATCACCCTTCTTATCACACTAATAGCAACAGCCTTCTTCGGTTACGTTCTCCCATGAGGACARATATCATTCTGAGCCGCAACCGTAATCACAAACTTACTTACCGCCGTCCCATATCTAGGCACATCACTAACAACATGATTATGAGGGGGCTTTGCAATCAATGATCCCACCCTGACACGATTCTTCGCACTACACTTTATTCTACCATTCGCAATTATCTCCTTATCCTCACTACACATTATTTTACTTCACGAAGAGGGCTCCAGCAACCCACTAGGAACCAACCCAGACATTGACAAAATTCCATTTCACCCCTACCACTCTCATAAAGACCTGCTCCTACTCACACTAATAATTCTACTCTTACTCACTATCGTCTCATTTCTCCCAGACATCTTCAATGACCCAGACAACTTCTCAAAAGCTAATCCCCTAGTAACACCACAACACATCAAACCAGAGTGATATTTCCTATTCGCCTATGGCATTTTACGATCCATCCCTAATAAACTGGGAGGCGCACTAGCCCTAGTAATATCAATCATAATCCTATTTACCATTCCATTCACACACACAGCCCGCCTACGCCCCATAACCTTCCGCCCACTATCACAACTAATATTTTGAACATTAGTATCAACATTCGCCACCATTACATGAGCCGCCACAAAACCAGTAGAACCCCCATTTATCATCATCAGCCAAGCAACTTCAATTCTATACTTCACATTCTTCCTATCTACCCCAATTCTAGGGTGGGTGGAAAACAAAATAATAAATATCTCCTCCCTAACACGGCAGTATAAAC--ACTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCATTCAGTCTATATACCACCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCTAAATAGTACT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAAGTGGACCAAGATGGGCTACATTTTCTAACACAGACTATACGAACAAAGACGTGAAACCGCCTTTCAAAGGTGGATTTAGCAGTAAGYTAAGAACATAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAAAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACTTATT--AAACCATATAATAGCTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACGTAAGCTATTCC-TTACATCACC-AGGCCAACAAGCCACTACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGCGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Antaresia stimsoni'] = '?ACCACAACCCTCAMTTCCTTTCAGCCATAGTCTGTAAATACAGCCATAGTCTGTAAA--TACAGACTATGGC--TGATACCGCCATATA-GAGCCGAAAACCCATATAATATGCCACACAATAAA------CTYTYTCCTYGGCCCCCCCCCTACCCCCCCCC--ATTAAAACATATGGGAA------AACAGCACAAATACATA--TTAAAGAATGTCCAATTAATCC--TCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTCCCCTTTAAATATTTTAGCCTAAAATGTCCTTCGTACAGAATATTAAG----TCCTCATTTTC-TGGTCGTTCAATGCAATCAGGATT--AATCA-TTCTT-ACTAACCATGGCTATCC-TTGAT-CTAGTGGTGTCCCTTAATTTAGTA-CTTCCCGTGAAATCCTCTATCCTTCCATAGAATGCTAACCATTCGACTTCTCACGTCCATATAGTGTTG-ATCCTCTT-AATAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTTCATTTAGTCTCATGTTCTGATCAGCTATGTCAAT-CCACCACTGGTAGCTCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAA????????????????????????????????????????????????????????????????????????TTCGGCTCAATACTATTAACATGTCTAGCCCTACAAGTATTAACCGGCTTTTTCCTAGCCGTTCATTATACAGCAAACATTAACCTAGCATTTTCATCCATCGTTCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTTATTTGTATTTATATTCACATCGCACGCGGACTATACTATGGATCCTACCTCAACAAAGAAACCTGAATATCCGGTATCACCCTGCTCATCACACTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTACCATACCTAGGCACATCGCTAACAACCTGACTGTGAGGGGGGTTCGCAATCAATGACCCCACCCTAACACGATTCTTCGCACTACACTTTATTCTACCATTCGCAATCATCTCCTTATCCTCCCTACACATTATCTTACTACACGAAGAAGGCTCAAGCAACCCACTAGGAACTAACCCAGACATCGACAAAATCCCATTCCACCCATACCACACCCACAAAGATCTACTCTTATTAACACTAATAGTTTTACTCCTATTCATTATCATTTCATTTTCCCCAGATATCTTTAATGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGGTACTTCCTATTTGCCTACGGTATCCTACGATCCATTCCCAACAAACTAGGAGGCGCACTAGCCTTAGTAATATCTATCATAATCCTATTCACCATCCCATTCACACACACAGCCCACCTACGCCCAATAACCTTCCGCCCACTCTCACAACTAACATTCTGAACACTGGTCTCAACATTCGCCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTTATCACCATCAGTCAAGTAACCTCAACACTATACTTCGCATTCTTCCTATCTATCCCAATTCTCGGATGAGTAGAAAACAAAATAATAAACATTTCATCCATAACACGGCAGTATTAAC---CTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGACAACCCACGATTAACCCAACCACTTCTAGCCACTCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTACT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAAGTGGACCAAGATGGGCTACATTTTCTAACACAGACTATACGAACAAAGACATGAAACTGTCTTTCAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACACCTATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTAATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCCTTCCC-TCATACCAA--AGGCCAACAAGCCATCACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGSTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Antaresia childreni'] = '?GCCACAACCCTCACTTCCTTCCAGCCATAGTCTGTAAATACAGCCATAGTCTGTAAA--TACAGACTATGGC--TGATGCCGCCATATA-GAGCCGAAAAACCATATAATATACCACACAATAAA------CTCTCTCCTCGGCCCCCCCCCTACCCCCCCCC--ATTAAAACATATGGGAA------AGCAGCACAAATACATA--TTAAAGAATGTCCAATTAATCC--TCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCAGGAATTCCCCTTTAAATATTTTAGTCTAAAATGTCCTTTGTACAGAATATTTAG----CCCTCATTCTC-TGGTCGTTCAATGCAATCGGGATT--AATCA-TTCTT-AATAACCATGACTATCC-TTGAT-CTAGTGGTGTCCCTTGATTTAGTA-CTTCCCGTGAAATCCTCTATCCTTTCATAGAATGCTAACCATTCGACTTCTCACGTCCATATAGTGTTG-ATCCTCTT-AATAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTTCATTTAGTCTCATGTTCTGATCAGCTATGTCAAT-CCACCACTGGTAGCTCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAAATGCCCCACCACTACATTCTAACCCTATTCGGCCTTCTGCCTGTAGCAACCAACATCTCAACATGATGAAACTTCGGCTCAATACTATTAACATGTCTAGCCCTACAAGTATTAACCGGTTTTTTCTTAGCTGTTCACTATACAGCAAACATTAACCTAGCATTTTCATCCATCGTTCACATTACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTTATTTGTATTTATATTCACATCGCACGCGGACTATACTATGGATCCTACCTCAACAAAGAAACCTGAATATCCGGTATCACCCTGCTCATCACACTAATAGCAACCGCCTTCTTCGGCTATGTTCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTACCATACCTGGGCACATCACTAACAACCTGACTGTGAGGAGGATTCGCAATCAATGACCCCACCCTAACGCGGTTCTTCGCGCTACACTTTATTCTACCATTCGCAATCATCTCCTTATCCTCCCTACACATTATTTTACTACACGAGGAAGGCTCCAGCAACCCACTAGGGACTAACCCAGACATCGATAAAATCCCATTCCACCCGTACCACACCCACAAAGACCTACTCCTACTAACACTAATAATTTTACTTCTATTAATTATCGTTTCATTTTCCCCGGACATCTTTCATGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATTAAACCAGAATTATACTTCCTATTTGCCTACGGCATCCTACGATCTATCCCCAACAAACTAGGAGGCGCACTGGCCTTAGTAATATCCATCATAATCCTATTCACCATCCCATTCACACACACAGCCCACCTACGGCCAATAACCTTCCGCCCACTCTCACAACTAATATTTTGAGCACTAATTTCAACATTCGTCACCATCACATGGGCCGCCACAAAGCCAGTAGAACCCCCATTTATCATCATCAGTCAAGTAACCTCAACACTATACTTCACATTCTTCCTATCTATCCCAATTCTCGGATGAGTAGAAAACAAAATAATAAACATTTCATCCATAACACGGCAGTATTAAC---CTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTAGCCACTCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTACT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAAGTGGGCCAAGATGGGCTACATTTTCTAACACAGACTATACGAACAAAGACATGAAACTGTCTTTCAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAACATTCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTCATAATAACTACTTTCAGTTGGGGCGASTTTGGAACAAAACAGAACTTCCAAACAACATGAGCCATCCCCTCATACAAC--AGGCCAACAAGCCACTACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCTTATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAACGGTACGTTTGTTCAACGATTACCAGTCCT'
    dna_dict['Antaresia perthensis'] = '?ATCA-AACCC------------------------AAAACTAAGCCACAGCCTGTTT--AAACAGGCTGTGGC--TGATGCCGCCATACA-AAGCCGAAATTCCATATAACACACCACAATATAAA------CTYTYTCCTYGGCCCCCCCCCTACCCCCCCCC--AACCAAACATATAAGAA------AACAGAACAGTGAACAA--TTAGAGATTCTCCAATTAACTC--TCCTATGTATAATCTTACATTAATGGTTTGCCCCATGGATATTAAGCAGGAATTTCCCTTCAAATATTTTAGTCTAAAATAGCCTTTGTACAGTCTATTTAG----TCCTCATTTTC-TGGTCGTTCAATGCAGCATGGATT--AATCA-TTCTT-ACCGATCATGACTATCC-TTGAT-CTAGTGGTGTCTCTTAATTTAGTA-CTTCCCGTGAAATCCTCTATCCTTCCCTAGAATGCTAACCATTCGACTTCTCACGTCCATATAATGCTA-ATCCTCCCTAAAAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTTCATTCTTTCTCATGTTCTGATCAGCTATGTCAAT-CCACCATTGGTTGCCCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAA????????????????????????????????????????????????????????????????????????TTCGGMTCAATACTACTAACATGTTTAGCCTTACAAGTACTAACCGGCTTTTTTTTGGCCGTCCACTACACAGCAAACATCAACCTGGCATTTTCATCCATCATTCACATTACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTGCACGCCATCGGAGCATCCATATTCTTCATTTGCATCTACATTCACATTGCACGCGGACTCTACTACGGATCCTACCTCAACAAAGAAACCTGGATATCGGGAATTACCCTCCTCATCACACTGATAGCTACCGCCTTCTTCGGCTACGTCCTCCCATGAGGACAGATATCATTCTGAGCCGCAACAGTAATCACCAACCTACTCACCGCTGTACCCTACCTAGGCACATCACTAACAACCTGACTATGAGGGGGGTTCGCAATCAACGATCCCACCCTGACGCGATTCTTCGCACTACACTTCATCCTACCATTCGCAATCATTTCTCTATCATCCTTACACATTATCTTACTACACGAAGAGGGCTCCAGCAACCCACTAGGAACCAACCCAGACATTGACAAAATCCCATTCCACCCATACCACTCCCACAAAGACCTGCTCCTACTAACACTTATAATTTTATTCTTATTCATTATCCTATCGTTCTCCCCGGACATCTTCAACGACCCAGATAACTTCTCAAAAGCTAACCCCCTAGTAACACCACAACACATTAAACCAGAATGGTACTTCCTGTTCGCCTACGGAATTCTACGATCCATTCCAAATAAATTAGGAGGCGCACTAGCCCTTGTGATATCAATCATAATCCTATTTACCATCCCATTCATACACACCGCCCACCTACGCCCAATAACCTTCCGCCCACTTTCACAACTTATATTCTGAACACTAGTCTCAACATTTGCCACCATTACATGAGCCGCCACAAAACCGGTAGAGCCCCCATACATCCTCATTAGCCAAGTGACCGCAACACTATACTTCACATTCTTTCTATCCATCCCAATCCTAGGATGAATAGAAAACAAAATAATAAACACCTCCTCCATAACACGGCASTACAACC--ACTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCTAATAACCGATAACCCACGATTAATCCAACCACTTCTAGCCGCTCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTAAT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAAGTGGACCAAGATGGGCTACATTTTCTAACACAGACCATACGAACAAAGGCATGAAACCGCCTTTCGAAGGTGGATTTAGCAGTAAGCTAAGAACACAACACCTAGCTGAAACCRTTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACCAACCTATT--AAACCCCATAATAGCTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAA-ATGAGCTCTCCC-TCACACCTC--AGGCCAACAAGCCACTATAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTCAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGCGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Antaresia melanocephalus'] = '?ACCACA----------CCTTCC-----------------CCAACCATAGTCTGTAACC--ACAGACTATGGT--CGATGTCTCAATATA-AAGCCAAAAATCTATATAAATAAA-ACACAATAAAG-----CTCTCTCCTCGGCCCCCCCC-TACMCCCCCC--ACAAGAAATATAGAAGAA------ACCAGCACATAAGACTA--TAAGGATTCCCCCCTTCTTTCC--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCCGGAATTCCCATTTAAATATTTTAATCTAAATTTGCCTTTGTACTTAAAATTCAG----TCCTCATTTCT-TGGTCGTTCAATGCAGCACGGATT--AATAG-ATCTT-ATTAACCATGGCTATCC-TTGAT-CTAGTGGTGTCCCATGATCTAGCT-CTTCCCGTGAAATCCTCTATCCTTCCTCTGAATGCTAACCATTCGACTTCTCACGTCCATATTATGCCA-ATCCTCCCTAACAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATTTTGCCTCATGTTCTGATCAGCTATGTTATT-CCACCACTGGTATCCCTTTTTT-CTCTGTACCTTTCATCTGACTACCATATATGCACACACACAGTAAATGCCCCACCACTACATCCTAACCCTATTTGGCCTTCTGCCTGTAGCAACTAACATCTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTTTAGGCCTACAAGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACCGCAAACATTAACCTGGCATTCTCATCTATCGTTCACATCTCCCGAGATGTCCCATACGGCTGAATAATACAAAACCTACATGCAATCGGAGCATCCATATTCTTCATCTGTATCTACATTCACATTGCACGAGGATTATACTACGGATCCTACCTTAACAAAGAAACCTGAATATCAGGCATCACACTACTCATCACACTAATAGCGACCGCTTTCTTCGGATATGTGCTTACATGAGGACAAATATCATTATGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTGCCCTACCTAGGCACATCTCTAACAACCTGACTATGAGGAGGATTCGCAATCAATGATCCTACCCTAACACGATTCTTCGCACTCCACTTCATTCTGCCATTCGCAATCATCTCCTTATCCTCACTACACATCATCCTACTTCACGAAGAAGGCTCCAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCATTCCACCCATACCACTCCCACAAAGATCTCCTCCTCCTAACACTAATAATCATAAGCCTATTCATCATCAGCTCGTTCTTCCCAGATATCTTCAACGACCCCGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAGCACATTAAACCAGAATGATACTTCCTATTCGCCTACGGAATTCTACGATCCATTCCCAACAAACTAGGAGGCGCATTAGCATTAGTAATATCAATTATAATCCTATTTATTATCCCATTTACACACACAGCCCGCCTTCGCCCTATAACCTTCCGCCCTCTATCACAACTAATATTTTGAACACTAGTATCAACATTCATCACCATTACATGGGCCGCCACAAAACCAGTAGAACCACCATTCATTGTAATCAGCCAAGTAACCTCATCACTATACTTCACATTCTTCTTATCAACACCCATTCTGGGATGAGCAGAAAACAAAATAATAAACATTTCATCCCAGCCATAACATGGCAACACACTTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCACCCAGTCTATATACCGCCGTCACCAGCCCACCTTGCAAAAGAAACATAGTGAGCCAAATAGTCCT-ACACTAACACGACAGGTCGAGGTGTAACTTATGAAGTGGATCAAGATGGGCTACATTTTCTAACACAGACCATACGAACAAAGACATGAAAACATCTTTCAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAGACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCCGGGATACCCATATAAGACCAGAAGACCCTGTGAAGCTTGAACTAACCTATT--AAACCACATAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAAAATGAGCCATCCC-TCATACCATTGAGGCCAACAAGCCACCACAAGACCCAGTAAAACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAGATGGTGCAGCCGCTATCAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Antaresia ramsayi'] = '?ACCACG----------CCTTCC-----------------CCA-CCATAGTCTGTAAA-TTACAGACTATGGT--CGTTGCCTCAACATA-AAGCCAAAAACCCATATAAACAAAAC--ATATAAA----CTCTCTCTCCTCGACCCCCCCC-TACCCCCCCC--ACAAGAAATATAGAAGAA------ACCAGCACATAAGACTA--TAAGGATTTCCCCCTCCTTTCC--CCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAGCCGGAATTCCCATTTAAATATTTTAATCTAAATT-GCCTTCGTACCTAAAATTCAG----TCCTCATTTCT-TGGTCGTTCAATGCAGCACGGATT--AATAG-ATCTT-ATTAACCATGGCTATCC-TTGAT-CTAGTGGTGTCCCATGATCTAGCT-CTTCCCGTGAAATCCTCTATCCTTCCTCTGAATGCTAACCATTCGACTTCTCACGTCCACGTTATGCCA-ATCCTCCCTAACAGCTTTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCAATGGTCCGGAACCATCCCTCCATACTAGCCTTTTCCAACACCTTTGGTTGCACCCTTTATA-TGGTACACTCTGCCTCATGTTCTGATCAGCTATGTTATT-CCACCACTGGTATCCCTTTTTT-CTCTG-ACCTTTCATCTGANTACCATATATGCACACACACAGTAA????????????????????????????????????????????????????????????????????????TTCGGCTCAATACTACTAACATGCTTAGGTYTACAAGTACTAACCGGCTTYTTCCTAGCCGTCCACTACACCGCAAACATTAACCTGGCATTCTCATCTATCGTTCACATCACCCGAGATGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGTATCTACATTCACATTGCACGAGGATTATACTACGGATCCTACCTTAACAAAGAAACCTGAATATCGGGTATTACATTACTCATTACACTAATAGCAACCGCCTTCTTCGGATATGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTCACCGCCGTACCATACCTAGGTACATCTCTAACAACCTGACTGTGAGGAGGATTCGCAATCAATGATCCCACCCTAACACGATTCTTTGCGCTACACTTCATCCTACCATTCGCAATCATCTCCTTGTCCTCACTACACATTATCTTACTTCACGAAGAAGGCTCCAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCGTACCACTCCCACAAAGACCTCCTCCTCCTAACGCTAATAATTATATCCTTATTTATTATCACCTCGTTCTTCCCAGACATCTTTAACGACCCCGACAATTTCTCAAAAGCCAACCCCCTAGTAACACCACAGCACATTAAACCAGAGTGGTACTTCCTATTTGCCTACGGCATCCTACGATCCATCCCCAACAAACTAGGAGGCGCACTAGCCCTAGTAATATCAATCATAATCCTATTTACCATTCCATTCACACACACAGCCTACCTTCGCCCTATAACCTTCCGCCCCCTATCACAACTAATATTCTGAACACTAGTCTCAACATTCATCACCATTACATGAGCCGCCACAAAACCAGTAGAACCACCATTCATTATAATTAGCCAAGTAACCTCAACACTATACTTCATATTCTTCTTATCAACACCCATCCTAGGATGAATAGAAAATAAAATAATAAACATTTCATCCATAACATGGCAACGCAA----TTTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAGCCTAGAGGAGCCTGTCCAATAACCGATAACCCACGATTAATCCAACCACTTTTGGCCACCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAATAGTTTC-ACACTAACACGACAGGTCGAGGTGTAACTCATGAAGTGGATCAAGATGGGCTACATTTTCTAACACAGACCATACGAACAAAGACATGAAAACATCTTTTAAAGGTGGATTTAGCAGTAAGCTAAGAACACAATACCTAGCTGAAACCGTTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATATACATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTTATAATAGGTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCCATCCC-TCATACCA---AGGCCAACAAGCCATTACAAGACCCAGTAACACTGATAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTCAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Python reticulatus'] = '??CCATCACCCTCACTTCCTCC--------------------AACCATAGCCAAATA-TTT---GGCTATGGTT-TCATGCCAAAATATATCAACCAAAAACCCATATTAATATAATGCTATAAAATGG-------TCCCTCGACCCCCCCCCTACCCCCCCC--AAAAAA--CATAAGGAAA------GTCCG-CACATCATAAACCTCGTACTTTTCCCTATTTTTT-GCTCCTATGTATAATCTTACATTAATGGCTTGCCCCATGGATAATAAGCAGGAATTTCCCTTTTAATATTTTAGTCTAAATTAGCCTTCGTACAGGTAATTCAGT----CCTCATTTTC-TGGTCGTTCAATGCAGCATGGATT--AATAA-TTGTT-GATAACCATGGATATCC-TTGAT-CTAGTTGTGTCCCTTGATTTAACA-CTTCCCGTGAAATCCTCTATCCTTCCGCGTAATGCTAACCATTCGACTTCTCACGTCCATTAATTGCTA-CTCCTCTTTACT-GGTTTTTCCAAGGCCGCTGGTTACACCTTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCCTTATA-TGGTACATATCACCTCATGTTCTGATCACCTATGCTAGA-CCACCACTGGTAGCTCTTTTTT-CTCTCTCCCTTTCACCTGACTACCATATATGCACACACACAGTAAATGCCCCACCATTATATCCTAACCTTATTTGGCCTTCTACCAGTAGCAACCAACATCTCAACCTGATGAAACTTCGGCTCAATATTACTAACATGTCTAGCCTTACAAGTACTAACCGGCTTTTTCCTAGCCGTCCATTACACAGCAAACATTAACCTAGCATTTTCATCCATCATCCACATCACCCGAGACGTCCCATACGGCTGAATAATACAAAACCTTCACGCTATCGGAGCATCCATATTCTTCATCTGCATCTACATCCACATCGCACGAGGCCTATACTACGGATCATACCTCAACAAAGAAACCTGAATATCAGGCATCACCCTACTCATCACACTAATAGCCACCGCTTTTTTCGGTTACGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTCACTGCCGTACCATACCTAGGTACATCACTAACAACCTGGCTTTGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGATTCTTTGCACTACATTTTATTTTACCATTCGCGATTATCTCATTATCCTCATTACACGTTATCTTACTCCACGAAGAAGGTTCTAGCAACCCCCTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCATACCACTCTTACAAAGACTTCCTCTTACTAACATTAATAGTCCTATCTCTATTCATTATCGTCTCATTCTTCCCAGATATCTTCAACGACCCAGACAATTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATTAAACCAGAATGATACTTCCTATTCGCTTACGGTATTCTACGATCCATCCCCAACCAACTTGGAGGAGCATTAGCCTTAGTAATATCTATTATAATCTTATTCACTATCCCATTCACACACACAGCTAATCTTCGTCCCATAACCTTCCGACCCCTCTATCAACTCATGTTCTGAACACTAGTCTCCACTTTTATTACTATCACATGAGCCGCTACAAAACCCGTAGAGCCCCCTTTTATCACTATTAGTCAAGTAACTTCAACACTTTATTTCACATTCTTCATCTCCATCCCATTCCTAGGATGAATAGAAAACAAAATAATACATCTCAACTCCATAACACGGCAATATAACCCA-TTGC?C---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCTAATAACCGATAACCCACGATTAACCCAACCACTTCTGGCCATCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAAAGAAACAAAGTGAGCCAAACAGTACT-ACACTAATACGACAGGTCGAGGTGTAACTAATGAAGTGGACCAAGATGGGCTACATTTTCTAACACAGAACACACGAACAAAGACATGAAACTGCCTTCCAAAGGTGGATTTAGCAGTAAGCTAAGGACAAAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAAAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTCCCAGTACAAAAGCTGGAATATCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACTTATT--AAACCCACTAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAAAACTTCCAAATAATATGAATAATCCC-TCATACCAT--AGGCCAACAAGCCACCACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTRTTCAACGATTAACAGTCCT'
    dna_dict['Python timoriensis'] = 'TA-CACCACCA------------------------------AGACCATAGTCGGTAAATC----GACTATGGTCTTTTTACGCCAAAAATACAACCAAAAATCCATATTAATATAGCAATATAAAATAG-------CCCCTCGACCCCCCCCCTACCCCCCCCC-ACAAAAA-TATAAAGAAA------ACCCG-TATGTCATAAACTCCGAATTTTTCCCTATTTTT--GCCCCTATGTATAATCATACATTATTGGCTTGCCCCATGGATAATAAGCAAGAATTCCCTTTTTAATATTTTAGTCTAAAATTGCCTTT-TACAAAAAACTCAGT----CCTCATTTCT-TGGTCGTTCAATGCAGCATGGGCT--AATAA-TTATT-AATAACCATGACTATCC-TTGAT-CTAGTTGTGTCTCTTAGTTTGGTA-CTTCCCGTGAAATCCTCTATCCTTCCGCGTAATGCTAACCATTCGACTTCTCACGTCCATTAATTGCTA-TTCCTCTTTACT-GGTTTTTCCAAGGCCGCTGGTTACACCTTCAAGATCATCTCAATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCCTTATA-TGGTACATTTTGCCTCATGTTCTGATCACCTATGCTAGA-CCACCTCTGGTAGCTCTTTTTT-CTCTCTCCCTTTCACCTGACTACCATATATGCACACACACAGTAA????????????????????????????????????????????????????????????????????????TTCGGCTCACTACTATTAACATGTCTAGCCCTACAAGTATTAACTGGTTTTTTCCTAGCCGTTCACTACACAGCAAACATTAACCTGGCATTTTCATCCATCATTCACATCACCCGAGACGTCCCATACGGTTGAATGATACAAAACCTCCACGCCATCGGAGCATCCATATTTTTCATTTGTATTTACATCCACATCGCACGAGGCCTATACTACGGATCATATYTTAACAAAGAAACTTGAATATCAGGCATCACCCTACTCATCACATTAATAGCTACTGCTTTCTTCGGATATGTTCTTCCATGAGGACAAATATCATTCTGRGCCGCAACTGTAATTACAAACCTACTTACAGCCGTACCATACCTGGGCACATCATTAACAACCTGACTCTGAGGCGGATTTGCAATCAACGACCCAACTCTAACACGATTCTTCGCACTACACTTTATCCTACCATTCGCAATCATCTCACTATCCTCACTACACATTATCTTACTCCATGAAGAAGGTTCTAGCAACCCCCTAGGAACTAACCCAGACATCGATAAAATCCCATTCCACCCCTATCATTCCCACAAAGACTTCCTCTTACTAATACTAATAATTCTATTTTTATTCATTATCGTTTCATTCTTCCCAGATATTTTCAACGACCCAGACAATTTCTCAAAAGCTAACCCACTAGTGACACCACAACACATTAAACCAGAATGATACTTCCTATTTGCCTATGGCATCCTACGATCTATTCCCAATAAACTAGGAGGAGCCCTAGCCCTAGTAATATCTATTATAATTCTATTCACCATCCCATTCACACACACAGCCTATCTTCGTCCAATAACCTTTCGCCCTTTCTCACAATTCATATTCTGAACACTAATCACTACATTCATCACCATCACATGAGCCGCTACAAAACCTGTAGAACCACCATTCATTATTATCAGCCAAGCGACATCAACACTATATTTCACCTTCTTCATTTCAATCCCTCTTCTAGGCTGAATAGAAAACAAAATAATACACCTCAATTCCATAACATGGCAATAAAACC--ACTGCTC---GCCAAACAACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCTAATAACCGATAACCCACGATTAACCCAACCATTTCTGGCCACCCAGTCTATATACCGCCGTCGCCAGCCCACCTTGCAAGAGAAACAAAGTGAGCCAAACAGTACT-ACACTAACACGACAGGTCGAGGTGTAACTAATGAAATGGACTAAGATGGGCTACATTTTCTAGCACAGAACACACGAACAAAGATATGAAACTATCTTCCAAAGGTGGATTTAGCAGTAAGCTAAGGACAAAATACCTAGCTGAAACCATTGCAATATTAAAGGCAATGCCTGCCCAGTGAGTATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTTTCAGTACAAAAGCTGAAATACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTATTAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCTATTCC-TCATAACAC--AGGCCAACAAGCCACTACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Python sebae'] = '?????????????????????????????????????????????????????????????????????????????CTTCCTCAGACAC-AAACTCA-ACCTCAAATAAAAATAAAAATAAT-----------CCTACCTCGGCCCCCCCCCTACCCCCCCC--ACTATTT-CATATGGAA-------TACAGGATATATAC-TTTGTTAGAAAAATCCATATTTTTTCTACCCTATGTATAATCTTACATTAATGGCTTGCCCCATGAATAATAAGCGGGAATTCCTAATAAAATATTTTAGCCTAAAATTGCCTTCGTACATAAAATT-AGC---TCCACATTTCTTTGGTCGTTCAATGCTGCANGGATTATAGTAC-TTCTT-AATACACATGACTATCC-TTGAT-CTAGTCGTCTCTCTTAACTTAACA-CTTCCCGTGAAATCCTCTATCCTTTCATA-CATGCTAACCATTCGACTTCTCACGTCCATAAGT-GCTA-CCCCTCTTCTCTTGCTCTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCGACGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATTTTGCCGCATGTTCTGATCACCTATGCTAGT-CAACCACTGGTATCCCTTTTTT-CTCT-TCCCTTTCACCTGACTACCATATATGCACACACACAGTTAATGCCACACCATTATATCTTAACCCTATTCGGACTCCTACCAGTAGCAACCAACATTTCAACATGATGAAATTTCGGCTCAATACTACTAACATGTTTAGCCTTACAAACGCTCACAGGCTTCTTCCTAGCTGTCCACTACACAGCAAACATTAACCTAGCATTCTCATCTATCATTCACATCATCCGTGACGTCCCACATGGCTGAATAATACAAAACCTGCACGCCATCGGCGCATCTATATTCTTTATTTGCATTTACATCCACATCGCACGAGGCCTATACTATGGATCCTATCTTAACAAAGAAACCTGAATATCAGGTATCACACTCCTCATCATCCTAATAGCAACCGCGTTCTTCGGCTACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACAGTAATCACCAACCTACTCACTGCTGTACCCTACCTAGGAACAACTCTAACAACCTGATTATGGGGAGGCTTCGCAATCAACGACCCAACCCTGACACGATTCTTCGCACTACACTTCATCTTACCATTCGCTATCATCTCTCTATCATCATTACACGTCATCCTACTACACGAAGAAGGCTCCAGCAACCCATTAGGTACCAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACTCATACAAAGACCTCCTTCTACTAACACTAATGATCCTTTCCCTACTAATCATCGTCTCATTCTTTCCAGATATTTTCAACGACCCAGACAACTTCTCAAAAGCCAATCCACTAGTCACACCCCAGCATATTAAACCAGAATGATACTTCCTATTCGCCTACGGAATCCTACGATCAATCCCAAACAAACTAGGAGGCGCCCTAGCCCTAGTAATATCAATCATAATTCTATTAACCATCCCATTCACACACACATCCACTATACGATCAATAACATTCCGACCATTATCACAACTAATATTCTGAACATTAGTATCAACATTCATCACTATCACATGAGCCGCCACAAAACCAGTAGAACCACCATTCATTATCATCAGCCAAGTAACCGCAACACTATACTTCACATTCTTCATTTCAACCCCCATCCTAGGATGACTAGAAAACAAAATAACAAATCACCCATCCATAACACGGCAACACTAT--AGTTGCTC---GCCAAACCACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCTAATAACCGATAACCCACGATTAATCCAACCACTTCTAGCCACCCAGTCTATATACCGCCGTCGCCAGCTCACCTTGCAAAAGAAACAAAGTGAACCAAATAGTCAT-ACACTAACACGACAGGTCGAGGTGTAACTAATGAAGTGGACTATGATGGGCTACATTTTCTAATACAGACCATACGAATAAAGACATGAAACTTTCTTTTGAAGGTGGATTTAGCAGTAAGATGAGAACACAATACTCAACTGAAACCATTGCAATATTAAAGGCAACGCCTGCCCAGTGAGAATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCACTTGTCTATTAATTGTAGACCCGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGAGATATACATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAAACTATT--AAACCCACTAATAGCTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCCATCCC-TCATATTAC--AGGCCAACAAGCCACCACAAGACCCAGTAACACTGATAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Python molurus'] = '?????????????????????????????????????????????????????????????????????????????CTTCCTCAGACAC-AAACTCA-ACCTCAAATAAAAATAAAAACAAT-----------CCTACCTCGGCCCCCCCCCTACCCCCCCCC-ACTATTT-CATATGGAA-------TACAGGATATATACATTTGTTAGAAAAATCCATATTTTTTCTACCCTATGTATAATCTTACATTAATGGCTTGCCCCATGAATAATAAGCGGGAATTCCTAATAAAATATTTTAGCCTAAAATTGCCTTCGTACATAAAATT-AGC---TCCACATTTCTTTGGTCGTTCAATGCTGCACGGATTATAGTAC-TTCTT-AATACACATGACTATCC-TTGAT-CTAGTCGTCTCTCTTAACTTAACA-CTTCCCGTGAAATCCTCTATCCTTTCATA-CATGCTAACCATTCGACTTCTCACGTCCATAAGT-GCTA-CCCCTCTTCTCTTGCTCTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCGACGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACATTTTGCCGCATGTTCTGATCACCTATGCTAGT-CAACCACTGGTATCCCTTTTTT-CTCTCTCCCTTTCACCTGACTACCATATATGCACACACACAGTTAATGCCCCACCACTATATCCTAACCTTATTTGGCCTCCTACCAGTAGCAACCAACATCTCAACATGATGAAACTTCGGCTCAATACTATTAGCATGCTTAGCCCTACAAGTATTAACCGGATTCTTCCTAGCCGTCCACTACACAGCAAACATCAACCTAGCATTCTCATCTATCATTCACATCACCCGCGATGTTCCATACGGCTGAATAATACAAAACCTACACGCTATCGGCGCATCCATATTCTTCATCTGCATCTACATTCACATCGCACGAGGACTATACTACGGCTCCTATCTAAATAAAGAAACCTGAATATCCGGAATTACACTACTCATCACACTTATGGCAACCGCCTTCTTCGGATATGTCCTCCCATGAGGGCAAATATCATTCTGAGCTGCAACCGTAATTACCAACCTATTAACCGCCGTACCATACTTAGGCACAACCCTAACAACCTGGTTATGAGGAGGATTCGCAATCAATGATCCCACCCTCACACGATTTTTTGCACTACATTTCATCCTACCATTCGCAATCATCTCCATATCATCACTACACATCATCCTACTCCACGAAGAAGGATCTAGCAACCCACTAGGAACAAACCCAGACATCGACAAAATCCCATTCCACCCATACCACTCATACAAAGACCTACTCTTCCTGACCCTAATAATCCTATTTATACTCATCATCGTCTCATTCTTCCCTGATATCTTCAACGACCCAGACAACTTCTCAAAAGCCAATCCACTAGTTACACCCCAACACATTAAACCAGAGTGGTACTTCCTATTCGCCTATGGGATCCTACGATCCATCCCAAACAAACTAGGTGGCGCATTAGCCCTAGTAATATCAATCATAATCCTATTTATTATCCCATTCACACATACAGCCCACTTCCGCCCAATAACTTTCCGCCCACTATCACAACTAATGTTCTGAACACTAGTATCAACATTCATCACTATCACATGAGCCGCCACAAAACCAGTAGAACCTCCATATATCATCATTAGCCAAGTAACAGCAACACTATACTTTATCTTCTTCATCTCTATACCCCTCCTAGGATGAATTGAAAACAAAATAACAAACACCCCCTCCATAACACGGCAACACTAT--AGTTGCTC---GCCAAACCA?TACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCTAATAACCGATAACCCACGATTAATCCAACCACTTCTAGCCACCCAGTCTATATACCGCCGTCGCCAGCTCACCTTGCAAAAGAAACAAAGTGAACCAAATAGTCAT-ACACTAACACGACAGGTCGAGGTGTAACTAATGAAGTGGACTATGATGGGCTACATTTTCTAATACAGACCATACGAATAAAGACATGAAACTTTCTTTTGAAGGTGGATTTAGCAGT?????????????????????????????????????TATTAAAGGCAATGCCTGCCCAGTGAGTATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTTTCAGTACAAAAGCTGAAATACCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCTATTAATAACTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCTATTCC-TCATAACAC--AGGCCAACAAGCCACTACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Python curtus'] = '?CCACAAAA-----------------------------------CCAT-----------------ATTAATYTT--CCCACCTATAAYTA-AACCCGAAATTCCCTATAAA--CACAACAAAAAATA-----CTCCTTCYTCGCCCCCCCCC-TACCCCCCCCCCAC-ATTT-AATATAAGAT------TCTGG--AATATACACACATCGTTAATTTCCATATTTTTT--ATGCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATAATAAGCCGGAATTCCATATTAAATATTTTAGCCTAAAATTGCCTTAGTACCTAAAACT-AGTCCTTCCTCATTTTC-TGGTCGTTCAATGCTGCATGGATT--AATCA-TTCTTTAACAGATATGTCTATCC-TTGAT-CTAGTCGTCTCTCTTAACCTGGCG-CTTCCCGTGAAATCCTCTATCCTTTCATA-CATGCTAACCATTCGACTTCTCACGTCCACATTT-GCTA-CCCCTCTTCTCTTGCTCTTTCCAAGGCCGCTGGTTACACTCTCAAGATCATCTCGACGGTCCGGAACCACCCCTCCATACTAGCCTTTTCCAAGACCTTTGGTCGCACCCCTTATA-TGGTACATTTTGCTTCATGTTCTGGTCACCTATGCTAGT-CAACCACTGGTTTCCCTTTTTTTCTCTGTACCTTTCATCTGACCACCATATATGCACACACACAGTTA????????????????????????????????????????????????????????????????????????TTCGGTTCAATATTACTCACTTGCCTAGTCCTACAAGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACAGCAAACATCAACCTAGCATTTTCCTCTATTATACACATCACCCGCGACGTCCCATACGGCTGAATAATACAAAACTTACACGCYATCGGCGCATCTATATTTTTCATCTGCATCTATATCCACATCGCACGAGGACTATATTACGGCTCCTATCTCAATAAAGAAACCTGAATGTCTGGCATTACACTCCTCATCACACTAATAGCAACCGCTTTTTTCGGATATGTCCTCCCATGAGGACAGATGTCATTCTGAGCCGCAACCGTAATCACCAATCTACTAACTGCTGTACCATACCTAGGCACAACCCTAACAACCTGATTATGAGGAGGGTTCGCAATCAACGACCCCACCCTCACACGATTCTTTGCACTACACTTCATCCTACCTTTCGCAATCATCTCTTTATCATCACTACACATTATTCTCCTTCATGAAGAAGGATCTAGCAATCCACTAGGAACCAACCCCGACATCGACAAAATCCCATTCCACCCATACCACTCTCACAAAGACTTCCTCCTACTCACACTATTAATCCTTTTCCTATTTATCATTGTCTCCTTCTTCCCAGACATTTTTAATGATCCAGATAACTTCTCAAAAGCTAACCCCCTTGTCACACCCCAACACATTAAACCAGAATGATACTTCTTATTCGCTTACGGAATCCTACGATCCATCCCAAACAAACTAGGTGGCGCATTAGCATTAGTAATATCAATCATAATCTTATTTACCATCCCATTCACACACACAGCCCATCTTCGCCCTATAACCTTCCGACCGTTCTCACAACTAATATTCTGAACACTAATCTCAACATTCATCACCATCACATGGGCCGCCACAAAACCAGTAGAACCCCCATATATTATCATCAGCCAAGCAACTGCAGCATTATACTTCACCTTCTTCATCTCTACACCCCTCCTGGGCTGAATAGAAAATAAAATAACAAACACTCCCT?CATAACACGGCAATACAA-----TTGCTC---GCCAAACCACTACGAGTGAAAACTTAAAACTTAAAGGACTTGACGGTACTTCACC--CAACCTAGAGGAGCCTGTCTAATAACCGATAACCCACGATTAATCCAACCACTTCTGGCCACCCAGTCTATATACCACCGTCGCCAGCTCACCTTGCAAAAGAAACAAAGTGAACAAAACAGTCATTACACTAACACGACAGGTCGAGGTGTAACTAATGAAGTGGACTAAGATGGGCTACATTTTCTAATACAGACCACACGAATAAAGACATGAAACTGTTCTCTGAAGGCGGATTTAGCAGTAAGCTAAGGACACAATACTCAACTGAAACCACTGCAATATTAAAGGCAATGCCTGCCCAGTGAGAACTCCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAAAGGCCACATGAGAGCCAAACTGTCTCTTGTAATCAATCAATTAAACTGATCTCCCAGTACAAAAGCTGAGATACTCATATAAGACCAGAAGACCCTGTGAAGCTTAAACCAACCTATT--AAACCTAATAATAGCTACTTTCAGTTGGGGCGACTTTGGAACAAAACAGAACTTCCAAACAACATGAGCAATTCC-TCATACCAC--AGGCCAACAAGCCACTACAAGACCCAGTAACACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTCAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGCGCAGCCGCTATTAAAGGKTCSTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Python regius'] = '?????????????????????????????????????????????????????????????????????????????TTACCTCAAT----AAACCCAAACCCACTATAAAAATATAA-----------------CCCCCTCGGCCCCCCCCCTTCCCCCCCCC-ACTTACA---TAGGAGGA------TTTAG-ATATATACACATATTAGGATTTTCCCTATCTTTTC-ACCCTATGTATAATCTTACATTAATGGTTTGCCCCATGAATATTAAACCAGAATTTCCAATTAAATATTTTAACCTAAAATTGCCTTCGTACACTACACC-AGT---CCCTCATTTCT-TGGTCGTTCAATGCTGCACGGATTATAGTAC-TTATT-AATGCTCATGTCTATCC-TTGGT-CTAGTGGTGTCTCTTAGTTTAACA-CTTCCCGTGAAATCCTCTATCCTTTCATA-CATGCTAACCATTCGACTTCTCACGTCCACAAGT-GCTA-CCCCTCCTCTCTTGCTCTTTCCAAGACCGCTGGTTACACTCTCAAGATCATCTCGATGGTCCGGAACCACCCCTCCATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATA-TGGTACAATTCACCTCATGTTCTGATCAGCTATGCCAGT-CCACCACTGGTATCCCTTTTTT-CTCTCTCCCTTTCACCTGACTACCATATATGCACACACACAGTTAATGCCCCACCACTATATCCTAACCCTCTTCGGCCTTCTACCAGTAGCAACCAACATCTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTCTAATGTTACAAGTACTTACCGGCTTCTTCCTAGCTGTCCACTATACAGCAAACATCAACCTAGCATTTTCATCCATTATCCATATTACCCGTGACGTCCCCTACGGCTGACTAATACAAAACCTACACGCCATCGGCGCATCGATATTCTTTATCTGCATCTACATTCACATCGCACGAGGACTATACTACGGCTCCCACCTCAATAAAGAAACCTGGATATCAGGTATTACACTTCTCATCACACTGATGGCAACCGCCTTCTTCGGGTATGTACTCCCATGAGGACAAATATCCTTCTGAGCCGCAACAGTAATTACCAACCTACTCACTGCTGTACCGTACCTAGGCGCAACCATAACCACCTGATTATGAGGAGGGTTCGCAATCAACGACCCCACCCTTACACGATTTTTCGCACTACACTTCATCCTACCATTCGAGATTATTTCCCTGTCATCATTACACATTATTTTACTTCACGAAGAAGGATCTAGCAACCCACTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCGTATCACTCATACAAAGACCTCCTACTACTAACACTAATACTACTAACACTTATAATCACCGTCTCCTTCTTCCCAGATATCTTCAACGACCCAGACAACTTCTCAAAAGCCAACCCATTAGTCACCCCCCAACACATTAAACCAGAGTGATACTTCCTATTCGCCTATGGCATCCTACGATCAATCCGAAATAAGCTTGGAGGGGCCTTAGCTCTAGTAATATCAATTATAATTATACTAACAGCCCCACTCACACACACAGCCCACCTCCGCCCAATAACCTTCCGACCACTTTCACAACTAATATTTTGAACCCTAATTTCAACATTCATTACCATTACATGAGCCGCCATAAAACCAGTAGAACCCCCATACATCATTATCAGCCAAACAACTTCAACACTATACTTCACCTTCTTTATCTCAACACCCATCCTAGGGTGAGTT?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????TATTAAAGGCAATGCCTGCCCAGTGAGAATTTCTTCAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTATTAATTGTAGACCAGTATGAA?GGCCACATGAAAGTCAAACTGTCTCTTGTAATTAATCAATTAAACTGATCTCCCAGTACAAAAGSTGGAATATCCATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACTTATT--AAACCCASTAATAACTACTTTCAGTTGGGGCRACTTTGGAACAAAACAAAACTTCCAAAYAAYATGAATAATCCC-TCATACCAT--AGGCCAACAAGCCACCACAAGACCCAGTAACACTGACANTTGANCCAAGTTACTCCAGGGATAACAGCGCCATCTTCTTTAARASCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCA?ATGG?????????????????????????????????????????????'
    dna_dict['Candola aspera'] = '?????????????????????????????????????????????????????????????????????????????????????????????????AAA----CTA-------------------------CT-CTCTGG-GACCCCCCCC-TACCCCCCCC--AGATAAACTATACTAAAATTTACCTGAGTACACTATGTAAATATTGTACATTAGTCTATATTTC--ATGCTATGTATAATCATACATTAATGATCTGCCCCATGGATAATAAGCAGGAATTTCCCTATTAATATTTCAGCCTATTAATGCCTTAGTACAGTCAGTGTGTC---ACCACATCAT--GGGTCGTTTTATGCAGCAAGGATTA-ACTA--TTATT-GGTAATCATGCCTATCC--TGATCCAAGTTGTC-CTCTTAATCTACCTA-CTCACGTGAAATCCTCTATCCTTCCAAGAATGGCTAACAGTCCTGCTTTTCACGTCCATATAT-GCTA-CCCCTCCTTTATTGCACTTTCCAAGGCCACTGGTTACTCTTTCAAGTCCATCTCAACGGTCCGGAACCATCCCTCTATACTAGCTTTTTCCAAGACCTTTGGTCGCACCCTTTATATTGGCGCAC-TGATCTCATGATCTGATCACCTATGCCAGTTCAACCACTGGTAACCTTTTTTT-CTCTCTACCTTTCACCTGACACCCATATAT-----------GTTGATGCCCCACCAACAAATACTAATTCTATTTGGACTACTACCAGTAGCTACTAACATTTCAACATGATGAAACTTCGGATCAATATTACTAACCTGTTCAGCATTACAAGTAATAACAGGATTCTCCTTATCAATACATTACACAGCAAATATTAACCTAGCATTCTGCTCCATTATCCATATTACACGAGACGTCCCACATGGATGAGTAATACAAAACCTACACGCCATCGGAGCATCGATATTTTTTATCTGTATCTACATATATATTGCACGAGGACTATACTATGGCTCGTACTTAAACAAAGAGACATGACTATCAGGCACCACACTACTAATCATACTTATAGCAACCGCATTCTTCGGATATGTACTACCATGAGGCCAAATATCATTCTGAGCAGCAACTGTTATCACTAATCTACTAACTGCAATTCCATACCTTGGAACTACAATGACAACTTGACTATGAGGCGGATTCGCAATCAACGATCCAACCCTGACCCGATTCTTTGCACTACACTTCATTCTACCATTCGGAATTATTTCAGTATCATCAATTCATATTATACTTCTGCACGAAGACGGATCAGGCAACCCATTAGGAACCAACTCAGACATTGATAAAATTCCATTCCATCCATACCACACATATAAAGATATTTTAATAATCTCAATTATAATTATTACACTACTACTAACCGTATCATTCTTCCCTGACATTATAAATGACCCAGAAAACTTCTCAAAAGCTAACCCTCTTGTCACACCACAACACATTAAACCAGAATGATATTTCTTATTTGCCTACGGAATTCTACGATCAATCCCAAACAAGCTCGGAGGGGCCCTAGCCCTAGTAATATCAATCATAATTCTATTCACCATGCCATTCACACACACATCACCTATACGATCACTAACATTCCGACCATTAATACAATTTATATTCTGAACATTAGTGGCAACATTCGTAATCATCACCTGAACCGCCACAAAACCTGTAGAACCACCATTCACTACCATCAGCCAAGTAGCATCAATTATCTACTTCACATTCTTCATATCTAACCCAATCCTAGGGTGACTTGAAAATAAAATTACAAAACACAACTACAAGGCAGATAAATACAA---AACTGCTC---GCCAAATTACTACGAGTGAAAACTTAAAACTTAAAAGACTTGACGGTACTTCACACTCAACCTAGAGGAGCCTGTCTAATAACCGATAATCCACGATTAACCCAACCACTCCTGGCCCCACAGCCTATATACCGCCGTCGCCAGCCCACCTTGTGAAAGAAACTAAGTGGACAAAACAG-CATCACACTAACACGACAGGTCGAGGTGTAGCATATGAAGTGGTCAAAGATGGGCTACATTTTCTAACTTAGACCAAACGAACAAGACAATGAAATAAACCTTTAAAGGCGGATTTAGTAGTAAGATGAGAACACAATACTCAACTGATACAAATGCAATATTAAAGGCAATGCCTGCCCAGTGAGAATATCTTTAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCGTAATCATTTGTCTACTAATTATAGACCAGTATGAAAGGCCACATGAGAGCCAAACTGTCTCTTATAGTTAATCAATTAAACTGATCTACTAGTACAAAAGCTAGAATAAACATATAAGACCAGAAGACCCTGTGAAGCTTAAACTAACCTATT--AAACCAACTAATAGCTACTTTAGGTTGGGGCGACCTTGGAACAAAACAGAACTTCCAAATAAAATGAGCTATAAC-TCATACAAC--AGACCAACAAGCCAATTTACGACCCAGTATAACTGACAATTGAACCAAGTTACTCCAGGGATAACAGCGCTATCTTCTTTAAGAGCCCATATCAAAAAGAAGGTTTACGACCTCGATGTTGGATCAGGACACCCAAATGGTGCAGAAGCTATTAAAGGTTCGTTTGTTCAACGATTAACAGTCCT'
    dna_dict['Morelia nauta'] = '????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????CCACTACATCTTAACCTTATTTGGCCTCCTACCGGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGCCTAGCGCTACAAGTACTAACCGGCTTCTTCCTAGCCGTACACTACACAGCGAACATTAACCTAGCATTTTCATCCATCATCCACATTACCCGAAACGTCCCATATGGCTGAATAATACAGAACCTACACGCTATCGGAGCATCCATATTCTTCATCTGCATTTACATCCACATCGCACGAGGACTATACTACGGATCATACCTAAACAAAGAAACCTGAATATCCGGCATCACCCTGCTCATCACACTAATAGCAACTGCCTTCTTCGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCTATTCCTTACCTAGGCACATCACTGACAACCTGACTTTGAGGCGGATTCGCAATCAACGACCCCACCTTAACACGTTTTTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCCTTATCCTCACTACACATCATCCTACTTCACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCATACCACACCTACAAAGACCTCCTCTTACTAACACTAATAATCCTGTTCCTATTCATCATCGTTTCATTCTTCCCTGATATC?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????'
    dna_dict['Morelia clastolepis'] = '????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????CCACTACATCTTAACCCTATTTGGCCTCCTACCAGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGCCTAGCACTACAAGTACTAACCGGCTTCTTCCTAGCCGTACACTACACAGCGAACATTAACCTAGCATTCTCATCCATCATCCACATTACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTACACGCTATCGGAGCATCCATATTCTTCATTTGCATTTACATTCACATCGCACGAGGACTATACTACGGATCTTACCTAAACAAAGAAACCTGAATATCCGGCATCACCCTGCTCATCACACTAATAGCAACTGCCTTCTTCGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCTATTCCTTACCTAGGCACATCACTGACAACCTGACTTTGAGGCGGATTCGCAATCAACGACCCCACCCTTACACGTTTTTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCCTTATCCTCACTACACGTCATCCTACTACACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCATACCACACCTACAAAGACCTCCTCTTACTAACACTAATAGTCCTGTTCCTATTCATCATCGTTTCATTCTTCCCTGATATC?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????'
    dna_dict['Morelia tracyae'] = '????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????CCACTACATCCTAACCCTATTTGGCCTCCTACCAGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTCTAGCTCTACAAGTACTAACCGGCTTCTTTCTAGCCGTACACTACACAGCAAACATTAACCTAGCATTTTCATCCATCATTCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAATCTACACGCTATCGGAGCATCCATATTCTTCATTTGCATTTACATCCACATCGCACGAGGACTATACTACGGATCTTACCTAAACAAAGAAACTTGAATATCAGGCATTACCCTACTCATCACACTAATAGCAACTGCCTTCTTTGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCCATCCCATACCTAGGCACATCTCTAACAACCTGACTTTGAGGCGGATTCGCAATTAACGACCCTACCCTAACACGCTTTTTCGCATTACACTTCATCCTACCATTCGCAATCATTTCCCTATCCTCACTACACATCATTCTACTCCACGAAGAAGGTTCTAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCATTCCATCCGTACCACTCCCACAAAGACCTCCTTTTACTAACACTAATAATCTTATTTCTATTCATCATCGTTTCATTCTTTCCTGATATT?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????'
    dna_dict['Morelia kinghorni'] = '????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????CCACTACATCTTAACCTTATTTGGCCTCCTACCAGTAGCAACCAACATTTCAACATGATGAAACTCCGGCTCAATACTACTAACATGTCTAGCACTACAAGTACTAACCGGCTTCTCCCTAGCTGTACACTACACAGCGAACATTAACCTAGCATTTTCATCCATCATCCACATTACCCGAGACGTCCCATATGGCTGAATAATACAGAACCTACACGCTATCGGAGCATCCATATTCTTCATCTGCATTTACATCCACATCGCACGAGGACTATACTACGGATCATACCTAAACAAAGAAACCTGAATATCCGGCATCACCCTGCTCATCACACTAATAGCAACTGCCTTCTTCGGATACGTTCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCTATTCCTTACCTAGGCACATCACTGACAACCTGACTTTGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGTTTTTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCCTTATCCTCACTACACGTCATCCTACTACACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCATACCACACCTACAAAGACCTCCTCTTACTAACACTAATAGTCCTGTTCCTATTCATCATCGTTTCATTCTTCCCTGATATT?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????'
    return dna_dict

def reference_taxon_set():
    taxon_set = dendropy.TaxonSet()
    taxon_set.new_taxon(label="Antaresia childreni", oid="Taxon4313741136")
    taxon_set.new_taxon(label="Antaresia maculosa", oid="Taxon4313741328")
    taxon_set.new_taxon(label="Antaresia melanocephalus", oid="Taxon4313741456")
    taxon_set.new_taxon(label="Antaresia perthensis", oid="Taxon4313741584")
    taxon_set.new_taxon(label="Antaresia ramsayi", oid="Taxon4313741712")
    taxon_set.new_taxon(label="Antaresia stimsoni", oid="Taxon4313741840")
    taxon_set.new_taxon(label="Apodora papuana", oid="Taxon4313741904")
    taxon_set.new_taxon(label="Bothrochilus boa", oid="Taxon4313741968")
    taxon_set.new_taxon(label="Candola aspera", oid="Taxon4313742032")
    taxon_set.new_taxon(label="Liasis albertisii", oid="Taxon4313742160")
    taxon_set.new_taxon(label="Liasis fuscus", oid="Taxon4313742224")
    taxon_set.new_taxon(label="Liasis mackloti", oid="Taxon4313742288")
    taxon_set.new_taxon(label="Liasis olivaceus", oid="Taxon4313742352")
    taxon_set.new_taxon(label="Loxocemus bicolor", oid="Taxon4313742480")
    taxon_set.new_taxon(label="Morelia amethistina", oid="Taxon4313742608")
    taxon_set.new_taxon(label="Morelia boeleni", oid="Taxon4313742672")
    taxon_set.new_taxon(label="Morelia bredli", oid="Taxon4313742736")
    taxon_set.new_taxon(label="Morelia carinata", oid="Taxon4313742800")
    taxon_set.new_taxon(label="Morelia clastolepis", oid="Taxon4313742928")
    taxon_set.new_taxon(label="Morelia kinghorni", oid="Taxon4313743056")
    taxon_set.new_taxon(label="Morelia nauta", oid="Taxon4313743120")
    taxon_set.new_taxon(label="Morelia oenpelliensis", oid="Taxon4313743248")
    taxon_set.new_taxon(label="Morelia spilota", oid="Taxon4313743312")
    taxon_set.new_taxon(label="Morelia tracyae", oid="Taxon4313759824")
    taxon_set.new_taxon(label="Morelia viridisN", oid="Taxon4313759888")
    taxon_set.new_taxon(label="Morelia viridisS", oid="Taxon4313759952")
    taxon_set.new_taxon(label="Python curtus", oid="Taxon4313760016")
    taxon_set.new_taxon(label="Python molurus", oid="Taxon4313760080")
    taxon_set.new_taxon(label="Python regius", oid="Taxon4313760144")
    taxon_set.new_taxon(label="Python reticulatus", oid="Taxon4313760272")
    taxon_set.new_taxon(label="Python sebae", oid="Taxon4313760336")
    taxon_set.new_taxon(label="Python timoriensis", oid="Taxon4313760464")
    taxon_set.new_taxon(label="Xenopeltis unicolor", oid="Taxon4313760592")
    return taxon_set

def reference_dna_matrix(taxon_set=None):
    if taxon_set is None:
        taxon_set = reference_taxon_set()
    dna = dendropy.DnaCharacterMatrix(taxon_set=taxon_set)
    assert len(dna.taxon_set)== 33
    sa = dna.default_state_alphabet
    dna_dict = reference_dna_dict()
    for t in dna.taxon_set:
        dna[t] = sa.get_states_as_vector(symbols=dna_dict[t.label])
    return dna

def reference_continuous_matrix(taxon_set=None):
    if taxon_set is None:
        taxon_set = reference_taxon_set()
    cvals = dendropy.ContinuousCharacterMatrix(taxon_set=taxon_set)
    assert len(cvals.taxon_set)== 33
    val_dict = reference_continuous_dict()
    for t in cvals.taxon_set:
        cvals[t] = dendropy.CharacterDataVector([dendropy.CharacterDataCell(value=v) for v in val_dict[t.label]])

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
    cells = [dendropy.CharacterDataCell(value=s, character_type=col_type) for s in col_type.state_alphabet.get_states(symbols=symbols)]
    return cells

def reference_standard_binary_matrix(taxon_set=None):
    if taxon_set is None:
        taxon_set = reference_taxon_set()
    ca1 = dendropy.StandardCharacterMatrix(taxon_set=taxon_set)
    assert len(ca1.taxon_set)== 33
    sa1 = _get_standard_state_alphabet("01")
    ca1.state_alphabets = [sa1]
    col_012 = dendropy.CharacterType(state_alphabet=sa1, label="COL_01")
    ca1.character_types = [col_012]
    for t in taxon_set:
        ca1[t] = dendropy.CharacterDataVector(_get_standard_cells(col_012, "0011000111-??0111001--?1101"))
    return ca1

def reference_standard_matrix(taxon_set=None):
    if taxon_set is None:
        taxon_set = reference_taxon_set()
    ca1 = dendropy.StandardCharacterMatrix(taxon_set=taxon_set)
    assert len(ca1.taxon_set)== 33
    sa1 = _get_standard_state_alphabet("012")
    sa2 = _get_standard_state_alphabet("XYZ")
    sa3 = _get_standard_state_alphabet("JKL")
    ca1.state_alphabets = [sa1, sa2, sa3]
    col_012 = dendropy.CharacterType(state_alphabet=sa1, label="COL_012")
    col_xyz = dendropy.CharacterType(state_alphabet=sa2, label="COL_XYZ")
    col_jkl = dendropy.CharacterType(state_alphabet=sa3, label="COL_JKL")
    ca1.character_types = [col_012, col_xyz, col_jkl]
    for t in taxon_set:
        ca1[t] = dendropy.CharacterDataVector(_get_standard_cells(col_012, "001122-??012")) \
               + dendropy.CharacterDataVector(_get_standard_cells(col_xyz, "XYZXYZ??-XXZ")) \
               + dendropy.CharacterDataVector(_get_standard_cells(col_jkl, "JKJLKL-??KJJ"))
        for c in ca1[t]:
            assert c.character_type is not None
            assert c.character_type in [col_012, col_xyz, col_jkl]
    return ca1

def reference_single_taxonset_dataset():
    taxon_set = reference_taxon_set()
    d = dendropy.DataSet(reference_tree_list(taxon_set=taxon_set),
                         reference_dna_matrix(taxon_set=taxon_set),
                         reference_standard_matrix(taxon_set=taxon_set),
                         taxon_set=taxon_set)
    return d

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
