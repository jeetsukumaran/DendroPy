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
Data and value generation
"""

from random import Random
from dendropy.test.support import pathmap
from dendropy.dataobject.tree import NodeRelationship
from dendropy.utility.messaging import get_logger
_LOG = get_logger(__name__)
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
    val_dict["Python regius"] = [-0.023008880157273925, -0.3273762612571307, -0.4836766440254265, 0.0868649474847456, -0.39427773481036776, -0.07527557074089436, -0.06524034584915558, -0.36041065830256214, 0.7356815654753754, -0.4794577538888527, 0.514275136235459, -0.10935592628052918, -0.08257435522293372, -0.7699934768255304, -0.301809357028146, 0.2569964973199855, 0.20845707860876755, -1.0039377269122811, 0.011427024404403255, -0.56110333944859, 0.08114980394297958, -0.3031120185981524, 0.025844253058697017, 0.14661192797306122, -0.40402607224351983, 0.028324650335011234, -0.1281376384530487, -0.158990133900242, 0.3476236449320067, -0.5372624194651937, -0.377377086233481, 0.2942868051603925, -0.31861217175662904, -0.7365595007629129, -0.10499728715227935, -0.49092061851032726, 0.11456767425541242, 0.16667869522628473, -0.9708471105782656, 0.9214351064712286, 0.20301748443453851, -0.8674284346005481, 0.03393927494767596, 0.18489927923783323, 0.5182827194940268, -0.3053455472446557, 0.591462049173616, 0.1403213687640485, 0.44177877110455166, -0.8729859965077695, -0.029392672059284405, -0.11207569247307034, -0.24187915988340464, -0.2481269348770918, 0.05982321355403042, -0.21703091240971573, 0.4376055185197054, -0.5937175739996388, -0.07532114765741295, 0.09974348711136494, -0.7018926675576702, 0.7416298629782743, -0.487022668518591, 0.027898594449524666, -0.3169829576009593, -0.3202825524267038, 0.40014335211730184, 0.7577264091654978, -0.8268201221277747, -0.3969504099222313, 0.26060658865237046, 0.1760548419179992, -0.3071403919071882, 0.36984903073645636, 0.12469930359403701, -0.08602172089702037, -0.6428260388458278, 0.3493662368519492, 0.1286609783240248, -0.7800889870960475, 0.6838919416665796, -0.08092064758534795, 0.012430953599999114, 0.3311213724636394, 0.15726063645851085, 0.05459789877660294, -0.2793240260289779, 0.7095734145354026, 0.03279814905206341, -0.18336299578541942, 0.2540440156687163, -0.14530463112989306, -0.32955372437709923, -0.45209483509633963, 0.04772568150951846, -0.5135170411593627, 0.022050264344369686, -0.14792364081446674, 0.6458504904989473, 0.07702042542341403]
    val_dict["Python sebae"] = [-0.31539438815801385, 0.2542472321993406, 0.4277957910429151, 1.1204917367330878, -0.10480849429563424, -0.047378202400136404, -0.05484326729951475, 0.15091974849759549, 0.6821437339086465, 0.17106263318820542, 0.08574313499160557, 0.10691436798512859, 0.6829924770468063, -0.2574790357667563, 0.07199740813502448, -0.21691967857652567, 0.43943899155221855, -0.24298286763192542, 0.560312742515378, -0.7732431356205428, 0.12255918564539699, -0.16732828684648998, -0.18824222279612723, -0.23866732748832048, 0.23509668633318276, 0.023790833540049558, -0.2193006328302687, 0.27188224215202444, 0.6488969006019866, -0.43464430339271465, -0.3031393355552258, -0.12607324430352584, -0.159314915022439, -0.4312131406594643, -0.2189961682067868, 0.08597968733585429, 0.13842254584386682, 0.4345512442553625, -0.6519164552750916, 0.5429560797650379, -0.5287023727865386, -0.41976550418965985, -0.7265493623072231, -0.010053744802354436, -0.03198441914101925, -0.5457191303743277, 0.8724204860322886, 0.1579480958332422, 0.03148277886212458, -0.2685972810205964, 0.23919220213541806, 0.12477931089598371, -0.4144066825532036, -0.00706844240416464, 0.5036314228115555, -0.15822396679476433, 0.5204105103821443, -0.6791936548933212, -0.07246240245464303, 0.457238680846623, -0.19722527835648, 0.025385565610340333, -0.6912120342742263, 0.437806713548006, -0.5608268129817706, -0.4323530490418329, 0.15501361513298056, 0.5864208477892863, 0.007905673355450171, -0.2741173991159548, 0.1842465461311708, 0.910320239826141, 0.21143069378304136, 0.511940915633296, 0.10524361666939416, 0.6127462390196827, -0.6529084655019325, 0.6687805027298241, 0.42945724596324636, -0.4769913917123535, 0.14854594103064334, -0.6157989070304246, 0.11100075234005997, 0.449244902296287, -0.4171928196186985, -0.2966852381179528, 0.052241095673798005, 0.26109514298953684, 0.25444482486477815, 0.034557658652939235, -0.0321916518173479, -0.05588489705559038, -0.05446081165685157, -0.5340517130076976, 0.4277679190538289, -0.2873237532207707, 0.16864517030323703, 0.01297026338748844, 0.19573370047869876, 0.8962510457938124]
    val_dict["Python molurus"] = [-0.18707621665457153, 0.009339909232670557, 0.48101828339609964, 0.6872941343985788, -0.12933553801885658, 0.3246601432624537, -0.2798314604805429, -0.05191684982259024, 1.1543756697339427, 0.11424755655053491, 0.3171405559795459, 0.4830167887955765, 0.022556511303798543, -0.6431787042845616, -0.2147367410954059, -0.5320113603868923, 0.019478928809939677, -0.5158401966917445, 0.0791908156819259, -0.510929188326607, -0.4967374584036555, 0.11747149214993245, 0.03775909586574568, -0.13329543976025357, 0.24625205652238358, -0.3743436778713754, -0.2789748251749228, -0.16790326189609162, 0.6650440743704755, -0.5408926759191599, 0.19606124956597887, 0.09808902634624866, -0.32042298484450576, -0.43661083598677514, 0.10356954385920919, 0.24365706521419306, -0.15807388756315877, 0.27840981735434484, -0.8208888556945719, 0.3952293503984462, -0.11059332411405412, -0.26947384175339195, -0.21841315034958905, 0.030991167496636163, 0.03515778826238011, -0.2550833554806984, -0.06457433911799226, -0.2967102358554325, 0.10855417567064216, 0.08197136099802156, 0.4082467733125769, -0.5013379023814333, -0.4088742197391203, -0.05191786238432959, 0.36935352798906623, 0.07672412305662885, 0.6208725090181777, -0.7481867711812662, -0.5486878469308446, 0.09199474912995081, -0.11745100642441111, 0.08715984862203398, -0.22059216131635662, 0.07230814229477961, -0.20542993738492352, -0.6069994629412333, -0.10857802735567153, 0.497918952130005, -0.3104345667385809, -0.21902056003533937, 0.11095515326758174, 1.1084774349193678, -0.17344101399717943, 0.7308810340948635, 0.17084927628927035, 0.3200035860403717, -0.5681245478267947, 0.16494112640962197, 0.3652197336801024, -0.6365091836150917, 0.039753165826301534, -0.05809751667733353, 0.3525315984406511, 0.04708225216629752, -0.5396504929432602, 0.1641927926267641, -0.14721710453657888, 0.41238555019116246, 0.36241275238929493, 0.05690343234286537, 0.3491127331445852, -0.18500338608598565, -0.21719746808264928, -0.44419850595907157, 0.27456020722027735, -0.2743770795596432, -0.3981173540956282, 0.20125746424159852, 0.041076267284013535, 0.3426372953471289]
    val_dict["Python curtus"] = [0.10976590468647396, 0.14431928383035347, 0.2746628252531358, 1.0220259801854428, -0.035337292515237914, 1.0085937636080395, -0.16662033511730787, -0.904859681052023, 1.1309704926383026, -0.0343109512236396, -0.4623018649848129, 0.44822809832079336, -0.01484409667445552, -0.1686894219652451, -0.06253363685225832, -0.16868911303621312, 0.27766080724353187, -0.888984029326404, -0.43592810444701713, -0.6063398862057153, -0.051699782732331556, 0.5387407543762467, 0.5627690529527994, -0.8603003182069475, 0.2913591012892804, -0.13456862508896744, -0.22009803185382573, -0.06707199355380193, 0.7066176550256982, -0.12942273941766158, -0.25380046982308546, 0.21774613727197412, -0.2524059180252121, -1.1485950891748256, 0.34470846632259183, -0.05614849900297224, 0.7116269439824894, -0.21874388413385554, -0.6617699325182923, 0.44202795028719516, -0.38968330108378824, -1.3366759869765446, -0.2837541596559585, -0.023947396246194325, 0.06770249542892497, -0.6591203670897263, 0.4951562792278729, -0.42539363397605956, -0.22074376028684325, 0.038657228679751354, 0.6433323965298784, 0.08486671950128488, -0.36199143653421306, 0.24035730589426132, 0.56539626675422, 0.2911016953418343, 0.6004804586108838, -0.3392155879350926, -0.5024848355477365, -0.986231417996906, -0.543732509413962, 0.11574954137613169, -0.13189636579027522, 0.0884139839151935, -0.3436828553504853, -0.3873512696551003, -0.018279335365908603, 0.4140095658694437, 0.44971246237092055, -0.32004418742453355, 0.06378113128313917, 0.39834485898730326, 0.19263161223211456, 0.18426666910186626, 0.40349744649809244, 0.138232930290323, 0.1279251796906782, 0.2819075502812161, 0.8388040097533934, -0.2663072092248807, -0.05650496431776583, 0.289481509631136, -0.1741697858750543, 0.09138703268922299, 0.18291047838665497, 0.2756235291883238, 0.5172462398743003, 0.3802548928370812, 0.5524308157228448, 0.8913015817204736, -0.5069497622357177, -0.30785495661469947, -0.1615423463460628, -0.4404592913259644, -0.43728523070518016, 0.8670933881461265, 2.76983225188987E-4, 0.05103074526163981, 0.7867307436088479, 0.1578040982109301]
    val_dict["Morelia bredli"] = [0.3587477886316148, 0.603433228555608, -0.25743663301448094, 0.8576350487121777, -0.8634438112836409, 0.5453676692004684, -0.03477988920819686, 0.25512395445519614, 0.08710780938236831, 0.8244038804607751, -0.6069177852531176, 0.19396854540716646, -0.044774284768091616, 0.4190914994603494, 0.784933323478815, -0.11566430995184585, 0.31957211700470367, -0.6869530619740116, -0.3139032841037713, -0.3122879521545662, -0.6221542037621389, 0.8042512977515924, -0.07096813810295684, 0.2956274116867744, 1.0542391554597481, -0.5853689396554248, -0.7684805525372509, 0.03953508939174083, 0.09271237637124541, 0.6604006350231054, -0.1377775499493728, 0.32251246909103104, 0.27144492812157106, 0.15420832638385545, 0.4347771367910068, -0.4674886996861548, 0.40603546524033846, 0.25098865086308747, -0.329657506834112, 0.7100730705945547, -0.3983613971604894, 0.001732724670675534, -0.31937172217847176, 0.0797072845187708, -0.3293833317454349, -0.07877139197883495, 0.14092665833506518, -0.3424852333781493, -0.059434709164277594, -0.4714225551527149, -0.59569331198294, 0.38472801493859765, 1.2397147124699979, 0.3000525560219655, -0.0416380475728704, 1.2716137839024793, -0.3493891785623035, -0.75652148488594, 0.3258123607275588, 0.01656467314779203, -0.4455230480171447, 0.38441617174148657, -0.6170811147193652, -0.4479049116839658, -0.04650368993803574, -0.47447774952634053, 0.5142166964293384, 0.8998182825087103, 0.5452867886245891, 0.7252866342867926, 0.2455111461520872, 4.628091506467258E-5, 0.01648494417402087, 0.27528068564820163, 0.07315242441885564, -0.4345011105560411, -0.2127845575024716, 0.5534662385920588, 0.24751793754246465, 0.22926636521490185, -0.06567271752872864, 0.18936429769730878, 0.03803490826996217, -0.11474220878589694, -0.03310168760787974, 0.14127310840168256, -0.4190810441798092, -0.021749266360020725, 0.41431250933502617, 0.1524044520437682, -0.47106175804657047, 0.323510414217908, -0.1851640445893023, 0.27053321222341975, 0.6470455677488843, 0.872025082800058, -0.631703749083969, -0.17404928738847897, -0.28517519190021146, 0.7245873335298507]
    val_dict["Morelia spilota"] = [0.4274232064675586, 0.4127166583223507, -0.17213747632024065, 0.7207213364972161, -1.2134152922877521, 0.48689867117513763, -0.14807705610777402, 0.2779523557279302, 0.45352011524967206, 1.0205507919981749, -0.5281556116020626, 0.27862214617334885, 0.04828078671152822, 0.6647755459959275, 0.29757133603327707, -0.06111118554550443, 0.27271211448778526, -0.8112922394329501, -0.06318261313265888, -0.4927260432987261, -0.17760549143502646, 0.7680093124884951, 0.01349067274594834, 0.5684842869138191, 0.89639700120945, -0.39012163028201485, -0.5904759685820444, 0.14825225757647617, 0.017425488011199152, 0.5776584812738281, -0.08103069763896442, -0.08722024730033584, 0.06358970899244071, 0.2082201000477834, 0.4444460409044422, -0.49569774166868763, 0.5314672412335587, 0.2001620192880268, -0.7267751548425443, 0.7144040252407743, -0.2911004414685827, 0.02632559169946222, -0.14300791326801696, -0.10418327241547506, -0.5935102241959289, -0.13695783742265882, -0.0023018929701700414, -0.2910660332892196, -0.06423742108894875, -0.7444370941795794, -0.9021865803074227, 0.5012581946792104, 0.8447057796612741, 0.17419913330600176, 0.1417323078815023, 1.068399782704533, -0.5117992088254663, -0.9773724638075312, 0.4415756908449959, 0.07110902602780028, -0.12942332283104782, 0.06422793823563994, -0.7988400847742853, -0.46858968900344605, -0.01969663259197009, -0.3651122195778427, 0.2864815379951528, 0.6394686577420778, 0.44142089627145065, 0.4024852931317912, 0.0745049723198878, -0.021783843014455392, -0.14789157321414953, 0.10484372720811824, -0.10587278386547641, -0.37480861295663154, -0.25800391593343924, 0.3421537454419167, 0.22437235740237577, 0.13159010038185287, -0.03768066717752068, 0.13469161405101585, -0.03452113805676958, -0.1871519854570888, -0.6776807203490841, -0.23163929155773472, -0.12514758704904955, -0.22710956513446195, 0.2230016995451608, 0.3250978371725558, -0.34178995436372833, 0.27739968208910765, -0.2479559594939194, 0.26978203102501397, 0.5955294101254911, 0.7185652390085704, -0.7541573753846796, -0.033938293489901114, 0.21508040410919785, 0.9051415456986822]
    val_dict["Morelia tracyae"] = [-0.1699644274726526, 0.0885724277268464, -0.312573568394239, 0.6682671713602445, -0.903327213579521, -0.05051272874694546, -0.2778686056552831, 0.5283218984921377, 0.3329129867539375, -0.049776500516163114, -0.5458217135149693, -0.0927924786512436, -0.07685830957481249, 0.546054096918058, 0.5484117927177553, -0.801109975680135, 0.003398839444569632, -1.0987650227383865, -0.3013508276454814, -1.0236757041283542, -0.37539391555130086, 0.4413412978008034, 0.13447930037681624, 0.20397658761752702, 0.5341261152974941, -0.27197123308740084, -0.9852716292639863, -0.5328247558615522, 0.08162085321256961, 0.4244008657584158, 0.10768938725062521, -0.47614196481158705, 0.0493877044237975, -0.26875599268291944, 0.598978441310301, -0.30675762240606963, 0.37843328505724627, -0.13831322463673698, -0.811483696647722, 0.3288735835128282, -0.1840593922363898, 0.38792213829243827, -0.21260299960165907, 0.043087278643165916, -0.3860305562630157, -0.06946912198409161, -0.2267874468363763, 0.06088025742237875, 0.17113855288355045, -0.33513988747336665, -0.5686471259700075, 0.3523903161407227, 0.7778972692846052, 0.06967286172286141, -0.20503831401730377, 1.0820182653438204, -0.36452523629858213, -0.5403786884686244, 0.17374278548001462, 0.9232062195071105, 0.06884857063983935, 0.01403483432378047, -0.5389457909847274, -0.09368640803593489, 0.22647564418759175, -0.606147240742557, 0.2826907960426558, 1.0197318988004629, 0.1523648370450871, 0.44709651914144155, -0.17571249999381358, 0.09861024155093844, -0.3119819442797643, 0.5038218577735465, -0.2400219862742705, -0.8078627723646374, -0.5879087636295605, 0.4672522645853993, -0.06402535539578105, 0.1653634507418068, -0.17953416056355515, -0.0520230385894202, -0.29455060776129266, -0.005712110722776859, -0.20696602811262807, 0.1415192691000448, -0.4330048684625141, -5.088974478654473E-4, -0.19781589971314717, -0.13379925588152577, -0.1512448435188552, 0.7719562359675525, 0.09458136596677363, 0.3267452576182409, 1.135111344247872, 1.004356835609907, -0.31750020610017704, 0.00948302529515066, 0.22060633226727172, 0.9243033431971135]
    val_dict["Morelia clastolepis"] = [-0.03707713735072906, 0.35794260531221317, -0.3436343101765865, 0.3087692847528298, -0.7144536349651675, 0.19222445394326743, 0.2604730829389691, 0.5103207929985335, 0.03338680665466398, 0.6123380576306876, -0.6988164886977764, 0.2953478009894051, 0.37434204855400144, 0.6198863127676062, 0.3133857965872694, -0.17563720788276177, -0.08501376204221198, -1.055950682838759, -0.45983095039902133, -0.7082099125568898, -0.4676323320893235, 0.7237338101284791, 0.5367973860983615, 0.1174717706107004, 1.0735824487324959, -0.3428712670835868, -0.7128332499010539, -0.15251690242456345, 0.0695013343477907, 0.49877105387146364, 0.5511200773012689, 0.08233334343698513, 0.2592690144186133, -0.5243775511447256, 0.33586932830147176, -0.2804284190056173, 0.4632650056274087, -0.18957190946172459, -0.74011549643159, 0.5037049339561083, -0.2868561617445218, 0.5094766476279525, -0.24115085411601772, -0.12479093914311916, -0.7238513716708153, 0.02233031552518923, -0.14588235404330147, 0.04989646160576294, 0.19301525805301167, -0.315107498078348, -0.3880724043729443, 0.44006913279055426, 1.076607277220065, 0.13631968546842693, -0.2539951566562218, 1.2289638982321138, -0.3134980898542859, -1.0635067052903178, 0.4248984453191168, 0.8425001504179472, -0.4166378462367024, -0.08345111840854302, -0.6503558793841607, -0.2825275532077805, 0.13244384775176657, -0.04528669913795777, 0.2473673680021425, 0.9662793523429363, 0.3002262425489777, 0.4754475627981618, 0.055176096957377366, 0.13463620798455384, -0.06491061689806411, 0.6629935591996923, -0.15269817834966262, -0.2099230615682772, -0.4256639794792506, 0.4652577401016138, -0.2453566621328856, 0.11236535957532579, -0.09809801089035208, 0.2877567619633156, 0.17069813422896074, 0.3051430209064465, 0.07795429142749007, 0.22788588041045788, -0.4057792733724809, -0.23338217525508148, -0.306197062198234, -0.15714812277484333, -0.5069065365271902, 0.6247229898061321, -0.3583715237531594, 0.4776799659951227, 0.9060492101005422, 0.9116115087352226, -0.6358295308561062, 0.20058080105077017, 0.03410299922098825, 0.6050071602021015]
    val_dict["Morelia kinghorni"] = [0.033501897358283114, 0.24843212652132496, -0.542815434295232, 0.27927182764493363, -0.6652476850699474, 0.23666683516589837, 0.36615109101539545, 0.4456386357787948, -0.10311723212733281, 0.6735569757940988, -0.7269477882992463, 0.28065118198756384, 0.3421706577011629, 0.7727904073150063, 0.31645178048987355, -0.13182025534879416, 0.038300862401558874, -0.9883648924886275, -0.5058502328006267, -0.5036388119519281, -0.37139643401381045, 0.7570575964296683, 0.32075885315862757, 0.3147312826281775, 1.122188180660774, -0.6492656107874067, -0.824549333780003, -0.24110326872513874, -0.17393070027826668, 0.4866745341114429, 0.3848530389165462, 6.091745261443959E-4, 0.29058854906385895, -0.5634418045843007, 0.24505460858739037, -0.3329521668648981, 0.4143543462905697, -0.03038344738928056, -0.596530119779163, 0.7485946168533237, -0.21169276361269856, 0.4639232555627056, -0.34271763936106436, -0.0770545677407272, -0.8459306564573983, 0.054532619529658446, -0.17448672232462364, 0.10140709849422955, 0.4519013304999042, -0.24288997165273385, -0.45157237224048435, 0.5745460495926141, 1.191688739648127, 0.1471308549697509, -0.4806257355383587, 1.2693497407772887, -0.12347377650439495, -1.0202889061554121, 0.3228606521631392, 0.8304855768131103, -0.31973888067432077, 0.04389790258002052, -0.6014066842501365, -0.23479427034274608, 0.10099055666118893, -0.08781129616386808, 0.47436339912329334, 0.925009956657815, 0.3036100340678936, 0.45094749197294015, 0.08973095040610148, 0.34924339660084897, -0.09493592022937528, 0.7665423477063743, -0.26345065159054715, -0.22803736688158738, -0.516405270017238, 0.3463257795820102, -0.1395101100138324, 0.063611244772595, -0.22375050408149422, 0.3582104596725286, 0.034611735821928835, 0.271681304894012, -0.06811663292395072, 0.204554709402806, -0.5539455457979836, -0.11574337806760931, -0.17428203136072123, -0.39033779018972303, -0.4604005808458809, 0.47746160965857964, -0.2037148842061201, 0.6987032772742676, 0.8383587755883207, 0.9770953251130824, -0.6945300320513321, 0.21276829379189968, 0.12736266966769405, 0.7310187102559746]
    val_dict["Morelia nauta"] = [0.005155782765228317, 0.38982553022833, -0.5239570598464205, 0.2721545211556563, -0.6655319239896669, 0.279306545237324, 0.35109968079983245, 0.40022016870507104, -0.04285374347944469, 0.46792383816542915, -0.7113070980447387, 0.17300053582030311, 0.33841838935905655, 0.6284475024629211, 0.3017839732644943, -0.12593050304814604, 0.027921894756392444, -1.0824364489331035, -0.7581768645040029, -0.6015931839878501, -0.3417306605710268, 0.8476781181882003, 0.4481664883488351, 0.11149606309712853, 1.0715873500293962, -0.5903360623005984, -0.9222291989674589, -0.07064310807172947, -0.018034139350114067, 0.5512526213387676, 0.7605706158522398, 0.13713674334596787, 0.3565614601559064, -0.6424439239407926, 0.3777757481389343, -0.16788655835247518, 0.5520184613146368, -0.2426522034917611, -0.6528738540117817, 0.7012243618697661, -0.2850147504336109, 0.4135177703378298, -0.3368673826115862, -0.19734234726408312, -0.7542579492107804, -0.1661597119719877, -0.11391864271760904, 0.010959643664941184, 0.46525627486572385, -0.2194758961620895, -0.5396915199454166, 0.43089005761198695, 1.0693058957587491, 0.15219534591756398, -0.5051834741304853, 1.3554694886111514, -0.2394943079838152, -1.0670281869028693, 0.45294754395722603, 0.8685508942587543, -0.2577050326796167, -0.048560419681001596, -0.5636898055631361, -0.3751642621366058, 0.07057383690624258, -0.10333545241650131, 0.2478929331753579, 1.010313685457738, 0.2633783579100697, 0.4225128730199012, 0.0876796869191772, 0.23063862894609277, -0.15827028871246318, 0.8841365481196617, -0.30479012748450834, -0.25692060573790937, -0.521523982949179, 0.3509684887674411, -0.19685047292078067, 0.08572239162255318, -0.09553660313920384, 0.038477480823738586, 0.033990037117208244, 0.26995373194252664, -0.008402980278248873, 0.11134311014843125, -0.4785168293030596, -0.1948972458343316, 0.10267289845004979, -0.46677571820278907, -0.5627160164626901, 0.5551308243518304, -0.18971813827214537, 0.5277690323838837, 0.9690298905403145, 0.9716352262160423, -0.6435916990324523, 0.2376971157863378, -0.08863507614595567, 0.7303977110115302]
    val_dict["Morelia amethistina"] = [-0.011836817897467992, 0.22523067795540225, -0.3473532299883031, 0.45647102776933834, -0.8486148445752227, 0.10260481586938078, 0.23130986646348495, 0.6134258889969036, 0.3778916907630002, 0.36172342032621435, -0.7256967594092795, 0.08830258451717658, 0.30997931667047746, 0.3379251673891973, 0.5034969822035824, -0.24023399082139263, 0.011694428104738992, -1.0109094111587489, -0.08534783938466067, -0.6150405933487324, -0.2086896966479364, 0.48422681891703145, 0.3043031168525088, -0.08206759560126176, 0.9176652817223657, -0.3379345063873965, -0.7891590608505402, -0.4026802854472461, 0.0022496589232980924, 0.3836891825982074, 0.4201000188120245, 0.08829207937100113, 0.13313937577724153, -0.06914072537053803, 0.2283159502534427, -0.2757351194040416, 0.4502406294342753, -0.27957518406336346, -0.5712577821416347, 0.46563904081075, -0.31088386776905014, 0.45611263816758385, -0.2929646252917377, 0.06972118428505142, -0.3412322551088887, -0.04070122328777086, -0.33426512327982266, 0.3138886611731808, -0.10051500673530507, -0.3162230398923152, -0.15275715135490253, 0.4205601253767368, 0.6920977048405492, -0.21609789660226722, -0.04123140009337322, 1.0342114887821126, -0.18224851830340114, -0.9144632290980939, 0.6180915060468283, 0.8035045964281375, -0.2813121833274621, -0.23716881272529816, -0.5730105696249338, -0.11286105918672076, 0.4120743605836106, -0.2540016273817333, 0.7045242032965231, 0.6848735986497184, 0.6337227823378955, 0.5350067627939491, -0.04276500748768968, 0.12439273579574528, 0.05638010807919622, 0.5059278534153366, -0.46655757890553623, -0.5650073946006489, -0.5080603958262591, 0.40600223521392625, 0.08794928313569492, 0.2670802890073899, -0.07677142770029025, 0.5148726415563641, -0.016773983983018183, 0.04881751825654437, -0.030897264833240234, 0.39099834435382613, -0.017738813955932603, -0.34705143698121227, -0.31559315604458105, -0.0456248612471494, -0.4386851967374133, 0.5170122050500962, -0.3493850640105582, 0.10108358198814377, 0.9392337946722048, 1.18721532364272, -0.7038974688390974, 0.23008084069941698, 0.3547542673313451, 0.538874667488326]
    val_dict["Morelia oenpelliensis"] = [0.185946506263553, 0.0064292744073738625, -0.12703326126907524, 0.5400976701764154, 0.03376477746817741, 0.8660623948496551, -0.10584911413931275, 0.4366998840173582, -0.08174493808267858, 0.973439429416308, -0.5965996928724259, 0.622830547558247, 0.18870609417119422, 0.6590370829342076, 0.3755720066561764, -0.2259051444694596, 0.38192207642933973, -0.8999768184879795, -0.5922300001235541, -1.0277939958014473, -0.46609045262306503, 1.0816098475768794, -0.4070986123824028, 0.11537324066845459, 1.5428593531794583, -0.526781266169535, -0.9493217996212178, -0.022312114788644477, 0.27740018583642384, 0.01403935416924429, 0.3900515959374839, -0.014207563164054138, -0.2201965552851654, -0.19216378252492514, 0.014805736095012345, -0.010357227855924772, 0.6921806317320275, -0.47560884106152285, -0.8604524363073407, 1.0006636993587696, -0.27622849880430117, 0.39162911841150344, 0.28992120626563495, 0.014719084273593419, -0.4758797129934135, 0.01872620105898684, -0.11595306296894908, 0.2932830308267763, 0.15017193730141332, -0.41750768319324616, -0.6784944967112507, 0.3939253906503808, 0.8301429509132872, 0.8959161154914758, 0.14758405946046305, 1.0899140656966286, -0.6353549426773507, -0.9451058407192237, 0.0884102702698478, 0.3423630071809738, -0.054176889419385144, -0.15733581764466192, -0.2921659639182025, -0.6969248771053735, -0.023774721356303125, -0.7187070159952755, 0.4171448288115374, 0.4834846910312626, 0.319895693762593, 0.909705766445309, -0.00278604177619432, -0.2153726521398257, -0.5322216458051067, 0.41999389447498564, 0.25452890032166353, 0.09446140702974498, 0.15103999773841953, 0.24318784603114946, -0.04275984319719496, 0.25608955396190847, 0.3957743488561476, 0.6272717824505529, 0.3583087067551982, 0.39544527680539526, -0.2629745146241981, -0.11477769813469824, 0.12439943817434158, 0.06711976992808129, 0.009102960630672119, 0.136393386703647, -0.1322402146728261, 0.32681026493465215, -0.2436537744428148, 0.12966969495815944, 0.49528965192762575, 0.3807859308181276, -0.8989841503677528, -0.10795811260801982, 0.17544716490811962, 1.233753578821423]
    val_dict["Antaresia maculosa"] = [0.002839626096003925, -0.032826359994880594, 0.09055758501757713, 0.9384643233064239, -0.5883900682747603, 0.3858370656553249, -0.11765316481096551, -0.3039699280804926, 0.3572216689914637, 0.7860704351758644, -0.270028557836763, 0.5178161759976889, 0.3912709564369733, 0.45812330860031825, -0.07765400633882441, -0.7245807672772353, 0.13471544478572722, -0.7282350752842146, -0.5662304124619115, -0.4965513830050981, -0.6907827956554651, 0.894558386874939, -0.2516373792026423, 0.8052591282902875, 0.33857274128364845, -0.7127307840296929, -0.9101082892647924, -0.06762352943822039, -0.5102120849994295, 0.23713946752398565, 0.6969843741359603, 0.2751055707887321, 0.28575480321622226, -0.4414332010519032, 0.7171139055843061, -0.267193325301383, 0.6191429248417899, -0.11924289408540861, -0.07709449355065146, 0.6187001022206062, -0.05106586857309061, 0.09696283633223227, -0.18533142256690685, -0.08383351700699028, -0.19774866251322154, -0.004423538848287417, 0.2939910684527002, 0.02690332989236871, 0.7024442762489947, 0.13697391074476037, -0.7301682397607685, -0.21386286067362187, 0.504416572107127, 0.6252096808901859, -0.7806178452916135, 0.7274914794573065, -0.11999009268983736, -0.47972256889892095, 0.5994834501157689, 0.4905435323002394, -0.7908614486146822, 0.38588197713324346, -0.7385650513534967, 0.5703829709677355, 0.4273344904037312, -0.25435538696337695, -0.04933284845024405, 0.766988296559713, 0.05578040454293465, 0.9217792682306938, -0.6723107434067742, 0.1582683248237647, 0.09960786653471304, -0.4229515489220638, 0.38824664595140523, -0.3573176601898326, -0.4863401233707075, 0.5728511091287488, 0.42144151819169984, 0.5421260475193026, -0.1457432863042542, -0.0418800032100174, 0.2702141394033819, 0.033742194555535995, 0.6287933663627078, 0.7136115142919387, -0.16488723262067245, 0.6397711154057586, 0.06857671348894434, 0.8806133633271775, -0.14158849190269163, -0.1639784483095842, -0.2943627190974247, 0.22809976972131574, 0.81881388972258, 0.7111540477793534, -0.36822111047538936, 0.016752503773740907, 0.21821643709680852, 0.7612166543991153]
    val_dict["Antaresia perthensis"] = [-0.08913600765363369, 0.4772211192998014, 0.2857961914820324, 0.7623465814870646, -0.6784935131993634, 0.06943574070695963, -0.5117161745759823, -0.1563693527538069, 1.3805260617283888, -0.253712007287484, -0.6537942254907322, -0.08251453081296925, -0.0623851897936985, 0.6097106462388975, -0.48506793150604677, -0.7439642152806942, 0.32695513798027287, -0.9057207073491097, -0.37760018430804676, -0.9735181146148673, -0.3779440489806595, 0.6795910805647434, -0.35670190452959977, 0.17745281765945997, 0.05053576539719351, -0.8318503783805611, -0.6463877540942058, -0.3687523217746055, 0.0842879784969879, 0.5552718042998698, 0.7331367472252088, 0.24249180919088814, -0.06742659364712823, -0.03674428903428581, 0.6351092072748779, -0.5680734059146336, 0.29557030633516357, 0.37020961188766305, -0.7101726995206713, 0.6612696126349908, -0.25174826071767487, -0.1874891805263852, -0.5297190463137529, -0.12762941387616558, 0.15773380649260782, -0.22453908869729264, 0.0960441200801529, 0.06528317962928522, 0.1720241464331843, -0.32058852035367313, -0.4683438442762015, 0.18177287434722883, 0.6354912206453549, 1.1654120416848315, -0.5177185457458535, 0.7702899916577255, 0.12378260256105779, -0.28645416766207527, 0.20123889439908155, 1.2537032652758366, -0.2848441371271925, 0.18035481714639917, -0.6296350623588352, 0.18057443677917684, 0.5804343599162849, -0.5430064711865462, 0.638368291720877, 0.7140777015474563, 0.27585536878564754, 0.3006734520949043, -0.04423864905804701, 0.290992833196964, -0.2814715267292712, -0.26824260010745526, -0.5519976406771887, -0.44599598122238926, -0.18439113414643643, 0.4296850907436057, 0.5495404939119948, 0.28042421197525935, 0.7047971259046433, -0.1656362993147985, 0.8360113347973778, 0.08593306964683345, -0.09817307644671948, 0.6364487996495719, -0.33604436900110024, 0.24437076478147035, -2.157232701678613E-4, 0.6632656914324483, -0.19522189057814704, -0.09390757982812091, -0.2921986029516171, 0.3788043568495052, 0.7347430000331268, 0.18911213680462297, -0.25832423563186474, -0.2698744439449491, 0.007999538669731632, 0.35168639456225903]
    val_dict["Antaresia stimsoni"] = [-0.40648719521134047, 0.23675535091710873, 0.08619387929948397, 1.2098848810878944, -0.5098204213748665, -0.11281491246366218, -0.3400339117918796, -0.3089824962890344, 0.1528222851533498, 0.5241039052265243, -0.8566694736282878, 0.29026441263038116, 0.5166558682432536, 0.7090467205316988, -0.05310529959167569, -0.6137533610792403, 0.05097927408284032, -0.7402565251368216, -0.555055588794788, -0.3731888724639882, -0.7501253361180468, 0.8198073694960062, 0.35655587741014616, 0.8813061905834532, 0.9445966138910141, -1.01153089916223, -0.842140817530806, -0.37886335427849316, -0.2782412105980887, 0.5049078457782598, 0.2561486778503227, 0.2386306052565317, 0.23936469036844427, -0.23619128795055205, 0.47231115890551706, -0.8650216204695769, 0.21029997921013546, -0.14191283205189995, -0.42879818722828084, 0.0369937879664933, -0.3854288754632681, -0.06399002093808477, 0.1714173275641168, 0.08569620806594597, -0.1272820392937999, -0.6178642172618601, 0.41116774522117194, 0.32449054063081856, 0.3129869367567312, 0.15239335660725473, -0.4369369682566825, -0.07457465785339185, 0.4311629002072114, 0.9232626020621153, -0.46909564071331855, 0.7876832587271773, -0.12110742332825422, -0.17121001117381962, 0.40262976415412205, 0.3916939334193109, -0.39670971740708344, 0.3978544673854532, -0.8508386807592044, 0.11210881449811538, 0.2259631116959701, -0.4189383888085422, 0.6354537741236648, 0.9234225430742603, 0.2574869253518442, 0.8990158865300201, -0.20430681668705017, 0.047429077687684354, -0.02768872059711415, -0.40374645385610713, 0.07021659000368052, -0.1729614480085439, -0.0464045757577297, 0.27457725942297656, -0.2990747757933044, 0.29561314671435046, 0.2998438121479311, 0.10983804705780681, 0.37355013362621114, -0.21725873926331885, 0.3257572925270171, 0.1888726676846511, -0.3343908801882921, 0.3066641463634656, -0.31516124799014433, 0.6758135953058167, -0.28549665373222594, -0.47996288415216853, -0.27328132080402034, 0.5930934404652679, 0.5254577149500171, 0.5534517650761928, -0.11439196673340137, -0.4189417069143072, 0.6785247374385972, 0.8684520366668568]
    val_dict["Antaresia childreni"] = [-0.31062540062709704, 0.084849629931467, 0.22787777149680238, 0.9429549766362486, -0.9232678443921927, -0.047332456307310505, -0.48709756554976896, -0.09983015495622685, 0.15252787349966374, 0.6235755428786839, -0.9402799908290855, 0.3289880628928645, 0.5586158314925544, 0.5729999314960169, 0.01914169643263971, -0.5502832733668974, -0.09064103467635883, -0.8812305105884695, -0.37590836104569325, -0.6106810059440985, -0.7416644527673136, 0.8272199155257772, 0.05299430770195698, 0.7369286383010298, 0.9403240349498985, -0.9025467154353102, -0.8445856520438539, -0.8103465534774779, 0.23377137374423318, 0.7993562025746643, 0.2951695535405939, 0.5751473080373566, 0.1206979587368674, -0.2845816826536531, 0.3528886224012758, -0.6861087668531184, 0.6060988120476374, 0.0062742853095476975, -0.5356094139690066, 0.27303078307672657, -0.20614297085831543, 0.1495786132785417, 0.018721209559665584, 0.1861768512392372, -0.38562441628571137, -0.4005698686143036, 0.24174890217810466, 0.07637061419935004, 0.03656920812470732, 0.20014816311577943, -0.3205314685941117, -0.011842105957329778, 0.29805074255209074, 0.8445982680331127, -0.1576369747927494, 0.7239069884224961, 0.07475258317584747, -0.01039591372245588, 0.3048224454071031, 0.6311290703431625, -0.48861202291372907, 0.5200315960335431, -0.984286833685339, 0.1948323189602291, 0.37449639245536825, -0.3438052378941151, 0.4526109455670057, 0.7625665762205454, 0.24735660412844895, 0.9375275980513265, -0.26472762393774246, -0.08839930400120807, -0.15693246445591152, -0.49214143465237437, -0.2759271688109757, -0.5053065808007744, 0.17586533958958867, 0.4957591699888322, -0.4562780949113041, 0.2940099956711146, 0.5397418798206065, -0.0993213213016534, 0.1588485896048032, -0.12440225649620897, 0.19841785799697728, 0.15018188004570113, 0.04249468564620407, 0.13807738620071297, -0.3375430715960652, 0.6213478497315497, -0.13024081191475054, -0.20043802126243362, -0.2119138869732462, 0.25918758745470216, 0.4550408974350064, 0.3567639932481735, 0.005022343269516338, -0.14774297082274862, 0.8005450643768083, 0.6084142679707094]
    val_dict["Morelia carinata"] = [0.3785900475395713, 0.21567477683554748, 0.13882961265359836, -0.05576574757196806, -0.41140034954383736, 0.33965605305768093, -0.5636074826452853, 0.19001897269295437, 0.9211415906636906, 0.29439499560056015, 0.03649384612567731, 0.2368684648875956, 0.08303756311682975, 0.694877244945812, -0.3606070317997774, -0.6889344732311767, -0.004855534113217284, -0.5493196062843791, -0.427947675030049, -0.9465353489866981, -0.39913064393284614, 0.8853217328497509, -0.5053791433840663, 0.6331241193952866, 0.63218804293164, -0.15399332170723307, -0.9863274750059097, 0.34982632519337886, -0.10609477246021237, 0.3591288213610088, 0.029320303034672035, -0.06780198754373218, -0.28563669558866, -0.5150793694419057, 0.12466130263066186, -0.43037241722854586, 0.45202561321413787, 0.49697909244826444, -1.3586850219859996, 0.8244416410232918, -0.04123579269102813, 0.006978708705367723, 0.09421352919728171, -0.09350202712133102, -0.05797816307102416, 0.10141909188517675, 0.3174319400361289, 0.28080049638250965, 0.060640735064414364, 0.01157678726507902, -0.3784656012107348, 0.10220973312702399, 0.7450117545295424, 0.6986369798187295, -0.4040509520934814, 0.33040268539078177, -0.14559118083996045, -0.3609347321861235, 1.1804553892682146, 0.6914164388944599, -0.14683879812102688, 0.7240933905576815, -1.0627028372209069, -0.5198763382162008, -0.3274969894290152, -0.7978308619292962, 0.6354680899586537, 0.06930703831156032, 0.1127750751770942, 0.4236386690252006, -0.3506717602208738, 0.5611349217230506, -0.07998413762303376, -0.5478825185097036, -0.21255587995317077, 0.04864939670070573, -0.47093702296031226, 1.0622134290741054, 0.6442150350957606, 0.5494547972424904, 0.5618475450720151, 0.013008139846752642, 0.3773495575610658, -0.16480917054104047, 0.3286552726701872, 0.5189511108009458, -0.6201569407240146, -0.5437392447562247, 0.18594799463115813, 0.8394115875729895, -0.7352525781703072, 0.2573816294072063, -0.3207632106612393, -0.01947227500237514, 0.7406986569381651, 0.07188249312273344, -0.20151300063626526, 0.07967751395234135, 0.09321659905110619, 0.28593129441591]
    val_dict["Morelia viridisN"] = [0.1111714344724616, -0.14480593458026614, 0.3497682367410207, 0.42329399493422437, -0.3100331837493862, 0.7315475054644202, -0.5995065777691273, 0.06748691980944516, 0.5549179460268524, 0.4968700511752563, -0.5369704035703842, 0.15510112793232042, 0.11521346720736216, 0.27737158510014026, -0.006097906462938799, -0.6362650768359486, 0.0916060079509427, -0.9812902115463947, -0.5780216809910884, -0.20097814572375516, -0.09195890670315499, 1.0556765048220804, -0.3668004559653733, 0.18863863559788896, 0.3615183407211907, -0.23880211647479244, -0.7963033625906867, 0.28541224100964324, 0.07602269275018286, 0.6365195731825851, 0.30641807453111336, 0.22770752379512682, -0.08463402468317308, -0.5027230402853575, 0.20014350912278014, -0.2672234452209207, 0.7184030101842863, -0.26535640235770597, -0.6986311368441767, 1.0005499298179026, -0.2934536796087359, 0.5322564001183132, 0.20699237603032444, -0.03654545814511247, -0.5128612077024721, 0.47985011716735926, 0.18883149882401243, 0.5157722736041314, -0.02481577119318429, 0.0074674876019556954, -0.22958029545310477, -0.06981375520340398, 0.1925812465697932, 0.6687848481507198, -0.27742134706498484, 0.6903773563654616, 0.023279699540754303, -0.6305595731730661, 0.7759311648499077, 0.7437188582678536, 0.08808450578124087, 0.741707374047365, -0.5961439811406901, -0.05535754160644714, 0.24385579578888716, -0.6211064740703236, 0.30763030606173, 0.04100025360151949, 0.47951361996525355, 0.27763012948168275, -0.20814854186318998, 0.2501099275815717, -0.24427877702518844, -0.6060337335963181, 0.4149237013887703, -0.5966320425536867, -0.6082462926547997, 0.36178054625571643, 0.3219710257762677, 0.20838815878602032, 0.21213104197191157, 0.2766663003276345, 0.16703219513792525, 0.20913159499206002, 0.3743392677302981, 0.5889468914685685, -0.11862392520583752, 0.275772357497963, -0.018333404560623323, 0.14793166813777178, -0.018216009766440666, -0.15758177156978068, -0.5373326950362689, 0.09541146454547254, 0.5444673226090256, 0.07025625281743571, -0.19832446729969122, 0.025663031151681892, 0.5953134501674955, 0.6946983855293218]
    val_dict["Morelia viridisS"] = [-0.19564184352811165, 0.18086062783775178, 0.41926106879093145, 0.7164136661357271, -0.3055499845222507, 0.6843729122992561, -0.6652716694777989, -0.21027372636673344, 0.7948884941362839, 0.5020506560562915, -0.45842990521870225, 0.22105149184089257, 0.2562372315816349, 0.25781480812778945, 0.29654908952617215, -0.7376415641824722, 0.45427652005681396, -1.0601721020389854, -1.206675705982461, -0.6478229320294691, 0.05680302514770896, 0.8138779607280914, -0.3365524577453585, 0.3130776851766619, 0.24769415615713966, -0.35286828678834725, -1.1265908458134217, 0.15591795923003782, 0.01075588988578148, 0.003926082551473531, 0.3212135724804471, -0.07605868318812646, -0.3466584975562842, -0.7532612535245857, 0.16144778803963286, -0.3085885734224757, 0.8142003958521301, 0.4062164678507304, -0.7227547676463887, 0.5718180564904107, -0.31130960483241427, 0.11552835591806357, 0.17323886357201881, 0.11774395052330544, -0.4394313162475426, -0.22654548653134618, 0.6233475746098402, 0.5377257467178648, -0.08980184772481109, -0.21177453974855506, -0.5289122826125073, -0.11450283755721047, 0.29451748196187777, 0.21442992000319405, -0.44416646884205, 0.7669859235337378, -0.37316717369519625, -0.9550589990503437, 0.7552007265031644, 0.9121045371761772, 0.15406330134320154, 0.704856090545985, -0.15271653951733882, -0.05620394530022581, 0.9470845523259186, -0.6608066796425952, 0.8183453292377049, 0.18222403844415674, 0.06126172583148301, 0.9294947731349995, -0.14316870845578225, 0.6641534493667429, -0.1653459869241536, -0.38330189675802456, -0.10509958335881726, -0.44042747871490073, -0.0015583904001526122, 0.9885255724938098, 0.1611592008974035, 0.5356655133351005, 0.3741846921996519, 0.3077915062367475, -0.12908196435050406, 0.0165436428407936, 0.09479797702590892, 0.9949600444711435, 0.0587626058205421, 0.3692965472891331, -0.4715454201704829, 0.09444285617109022, 0.14067039575487647, -0.2574262962397543, -0.2778023139719551, 0.18877809510208005, 0.4199474045982234, -0.21389748456914023, 0.22201484822266287, 0.14769941281638863, 0.37241996522110854, 0.8434744106176317]
    val_dict["Apodora papuana"] = [-0.333349760315575, 0.9379996812088196, 0.19041684194223066, 0.6338375265053364, 0.2903113251734366, 0.3417238122253093, 0.45252819042042225, 0.4259076534852059, -0.22218802345771782, 0.7234116255195612, -0.19811505924279887, 0.0486237497874519, 0.3929235242572352, 0.9537723434854293, 0.24142780106204398, -0.758469925674385, 0.01869659276528507, -0.08538282821731702, 0.25086010066942643, -0.05129166155828291, -0.33419583288593857, 0.6137143748908893, -0.933545403872212, 0.6513735215717615, 0.4097760396619297, -0.659606444172233, 0.26172193964023127, -0.09177325960626154, -0.33267483084959903, 0.41294065227937793, 0.30942996488585806, 0.5579073672034756, 0.14557678384312076, -0.3058497883140783, 0.45478076761572117, -0.41618780324600374, 0.07512295737119218, 0.1971571271388458, -0.9644923476273801, 0.41387610678703085, -0.46521485209780633, 0.2775529146707806, 0.320882648143625, -0.03757471492394843, -0.2534946846422394, -0.20062385335197044, -0.2618660558299431, 0.6190147405478831, 0.3110893428331795, -0.2085543201528925, -0.607873455511589, 0.4240419928596113, 0.6099351667380486, -0.15101861078846712, -0.2036166295879093, 0.8860147217679363, -0.4652492473118644, -0.7167175226151423, 0.37307955936896053, 0.7498341043273146, -0.6885330503294839, -0.04329215822805724, -0.6225934171689274, 0.08909625840833078, -0.13474707326011656, -0.501894840412574, -0.035746035762169476, -0.11371662578528899, 0.5572799418642219, 0.25598536895343316, -0.07919520026776655, 0.12945647125517185, -0.2092505918708472, 0.1922574785178165, -0.1497481808331125, -0.031038992315084557, -0.1492525280966372, 0.8686993607522779, -0.19851693379900714, 0.34238960607818086, 0.5250315859569173, -0.1737213123886101, 0.10644816478620384, 0.11293151897403667, -0.13664610496335278, 0.5210817011041057, -0.540620884261954, 0.6597626291216273, 0.15869520249020055, 0.27277691435423074, 0.17102993414582873, 0.10278798520296702, 0.2583660786056935, 0.33106638220832113, 0.12283401311470414, 0.183847466985188, -0.0011858007569514606, -0.37685816228180835, 0.3007088223897098, 0.1830570346654184]
    val_dict["Liasis olivaceus"] = [0.11615122835010966, 0.4578578136174916, 0.6747894923078218, 0.437907751253062, -0.08890103278036512, 0.5738523577811792, -0.3367883090818582, 0.5952867116227536, 0.4291871108141066, 0.5392351515140121, 0.028598507806426546, 0.5697474619073446, 0.2973236758733811, 0.056330768671995046, 0.686024235857726, -0.6694596075872083, -0.49947799024617423, -0.18494890387597374, -0.16640185430094784, -0.7977067194210208, -0.29772076721232704, 0.5553769441912235, -0.5782540904070304, 0.9612453370465484, 0.5406291678254335, -0.6145689659375533, -0.3809234319072589, 0.09605593568973045, 0.41030614142052235, 0.8622353573804502, -0.1196921232283682, -0.06286865689756418, -0.3638647473989844, -0.21527458528778928, -0.1002372060754051, -0.45664129803933384, 0.849049249443527, -0.38491961057424806, -0.8208642922185846, 0.7639279978875683, -0.6454034294831577, 0.2519392931638466, -0.11524348051533903, -0.19022761494069265, -0.5092389863404093, 0.18097251699606748, 0.3469055304953014, 0.23570105885418854, 0.41007196304578336, -0.1209625467542689, -0.5850642925670713, -0.3713858486537106, 0.34101214797649404, 0.7642091548410606, 0.04088257752637173, 0.4400902827630767, -0.31424679198249134, -0.5330852586725126, 0.5480594346581111, 0.1947290683036414, -0.5967973012600587, 0.3990346160253636, -0.3422228167532898, 0.06509855398010281, 0.2048072160927839, -0.03566436104937537, 0.41539808490327546, 0.11316966218390016, 0.6300335593010278, 0.3430739112145449, -0.1491009260885835, 0.11239848442514723, 0.26278220841547195, 0.1532415915542205, 0.0901231653196563, -0.2021330598478559, -0.3410565388033444, 0.19144561457710463, -0.16280299495595196, 0.6716868058677158, 0.23669471252951457, -0.20005209977272953, 0.36837680974294035, 0.11635221570709653, -0.08604098063965387, 0.5932960659658059, -0.344354977997956, 0.866641102459285, 0.1930961546386517, 0.10406470739921234, -0.2795602072850429, -0.03012397136283619, 0.6131583151931838, 0.23235696324096317, -0.1823250295865933, -0.1583613901864687, -0.46534176827615603, 0.08512159006217182, 0.6892997000805889, -0.06485173256509219]
    val_dict["Liasis fuscus"] = [0.2429816339889912, 0.047215713097623194, 0.9055511875510993, 0.3916999392235513, 0.18523888916598305, 0.3444424259927913, -0.2601595867061947, 0.6868266093786489, -0.02005727567847465, 0.36684562114975433, -0.6761265204652485, 0.10207496811019212, 0.1974615767383426, 0.20252017202962125, -0.03219898840480273, -0.24461669815664686, -0.7668768154523948, -0.6984658964738907, 0.22435685261074123, -0.619804571654976, -0.04609698024052085, 0.7315986510356908, 0.09064077792556077, -0.08180346568468141, 1.2640362513967471, -0.4320567810247121, 0.015574560957995254, -0.04211489515776491, 0.08380855772701613, 0.9946889541844568, -0.07849187436612562, -0.5978693379094346, -0.05901023693534459, -0.7803266627524379, 0.038882407875277206, -0.18220796999786676, 2.464685320449078E-4, 0.12681390476979676, -0.47575693654186485, 0.8402732192061915, 0.03171875511653924, 0.4334915347705033, -0.06384054480021137, 0.10597721248395278, -0.9277569372637124, -0.18899261391367952, 0.24632675731941683, 0.5957073151621006, 0.17551850015428574, 0.002261943385033681, -0.5054657062122668, -0.30797870442219155, 0.31009820985576164, 0.7180805435406501, -0.2990435869692172, 0.3465664441169066, -0.3631894152241846, -0.17474457633559984, -0.23271187522681233, 0.4919628749111339, 0.1239453765340133, -0.1173452025629842, -0.30820679651090954, 0.5455033241750702, 0.48752183379298136, -0.056979314313508334, 0.5266407648649886, 0.5544937522024675, 0.3454534836018408, -0.12874118486875435, 0.08621662660045762, -0.16728185730034378, -0.5077114759807885, -0.5161003791961218, -0.06673070054467724, -0.26629161148347213, -0.5277005519437552, -0.00864777709727544, 0.051304964785276796, 0.8024936631771782, 0.08708757643714327, 0.09657011291634049, -0.08762759513949456, -0.004200833322482361, 0.1334301787128342, 0.05802466346885538, -0.508365499990975, 0.18683762312621283, 0.00397934623736762, 0.5967963472054332, 0.25184666914517273, -0.5041233565410496, -0.1926783442826732, 0.27584878964854587, 0.03245285320437058, 0.05600373502385324, -0.7255926025637067, -0.10337218186663426, 0.4758651685913429, 0.41588341116647737]
    val_dict["Liasis mackloti"] = [0.24911728635205824, 0.08564515973522596, 0.4861426472723624, 0.4819535328186636, 0.06411710072232749, 0.34244661423205597, -0.2388899711903706, 0.4587576122088772, 0.09387481380779211, 0.34024393901768735, -0.38186850350680335, 0.08307133903373749, 0.14248402908364502, 0.2366902392607363, -0.19026956561948288, -0.4703880467588185, -0.8626033209599094, -0.6009742370125993, 0.13820825965494804, -0.6187043767449669, -0.03809728432353597, 0.8662073566061005, -0.24725455813905195, -0.054256979351333265, 1.0848565668755372, -0.9369467338601207, -0.43098528220575083, -0.07159746657226429, 0.20525832788499324, 0.7988899741836558, 0.12115049727400766, -0.4864719541210042, 0.040135694468528174, -0.4455379148861992, 0.06256760289463932, -0.1260467351685523, 0.34767249812222434, -0.2056092894664645, -0.379891871288721, 0.8805329252132924, -0.2250897185706266, 0.3881639984651543, -0.05125761445515661, -0.17494556691480223, -0.8478325167117698, -0.4777469317917513, 0.23735970596202163, 0.36247829808799714, 0.2419532930455089, -0.0881913816042755, -0.45979548373875634, -0.09497776894619958, 0.45012524178688657, 0.6626694419173393, -0.41493817579915054, 0.5805197267943354, -0.23842936072333698, -0.33374415555430026, -0.21520685457262057, 0.47856608847902465, -0.010949092979807881, 0.0046618750534558814, -0.5578745359581071, 0.13179137524780843, 0.4040569588663244, -0.10534847791991656, 0.40595517757339517, 0.39969927551697026, 0.6616764339442434, -0.045125657294826274, 0.014768829158566782, 0.0982904656210586, -0.3465356284641471, -0.7923968625931006, -0.04074220482779807, -0.3265451342975338, -0.5719654930266136, -0.04708298768059306, 0.04116364422845492, 0.6884493643261643, 0.2694580677427365, 0.01812526533188364, -0.07853279195336091, -0.07690601080221673, 0.22864838120437175, -0.10241509946412033, -0.3769923628391627, 0.16253838992359204, 0.26926381744451533, 0.3946795944383332, 0.18513975306625854, -0.34880042465953875, -0.3336552978510261, 0.4511912186473644, -0.19529239606682536, 0.23240909562235743, -0.5389184903590284, -0.13028072126133233, 0.4903143648130598, 0.44560725794387895]
    val_dict["Antaresia melanocephalus"] = [-0.1828763266295226, 0.39881764673546316, 0.21998570133560938, 0.7721832888909403, -1.0548556179491246, 0.5649552603212312, 0.037112215197546594, -0.09235646994301505, 0.8526679420108187, 0.3771586587252546, -0.16453343392752717, 0.632027327135727, 0.49665477373244793, -0.24798665400043762, -0.033913615594106994, -1.0699789962257487, -0.6411947147474546, -0.8620701005872653, -0.3492871533543157, -0.9119897714688808, 0.19686243669444753, 0.3741610719080052, -0.30190128765507596, 0.4141682799339195, 0.12473419148761744, -0.5285516384154648, -0.7212251393084541, -0.011283804775608627, -0.6633987794031978, 0.34949108495299885, 0.21133107487281133, -0.6411445278239039, 0.36109575297867413, -0.5706094725715397, 0.19548119854872817, -0.0779157949284488, 0.5373498753300834, -0.39312273815907134, -0.38539421289516085, 0.9598720630805664, -0.14205408880539366, 0.5341239693594197, -0.578028865550934, -0.44193179406120653, -0.7242813837018315, 0.0030601088587691727, 0.035798222407850175, 0.23048124578969392, 0.19540624375249285, -0.737851855137215, -0.4430454004938424, -0.09746626815378003, 0.942852251615075, 1.1790271192285526, -0.3063806585222091, 0.8755837043363024, 0.18565314553013665, -0.8627062344699071, 0.33293962770286595, 0.9074834357731019, -0.03736356289181386, 0.5981397953665046, -0.25696807687172174, 0.21628199798728837, 0.622344693886917, -0.3221852509859547, 0.6560443340547655, 0.06774475194646042, 0.7188110267000476, 0.7536063090179518, -0.3012149192234952, 0.10120887341493474, -0.03031843981298784, -0.018023554007912904, 0.06782587181152834, -0.16861646839619104, -0.40660289882431344, 0.5397715767565382, 0.3897044873265439, -0.3153164970601905, 0.2128151841530713, 0.7888299358755099, 0.30169916947324, 0.5632494474432732, 0.054163549475901546, 0.44177627117524065, -0.3590257255127771, 0.45707671478021716, -0.021752343358095644, 0.6489369089509813, -0.17444558903035864, -0.27685873630755276, 0.45031727005078126, 0.11275487956068997, 0.28806801548916294, 0.11970392972857646, -0.3547189713955756, -0.06475189310851587, 0.6382951103938268, 0.580535522591896]
    val_dict["Antaresia ramsayi"] = [0.0550136437894001, 0.34877214117083044, 0.4265010988914893, 0.8660687362836188, -0.6981509093419199, 0.5734553803890774, 0.0054702755840887185, 0.07543298937261317, 1.0003971652821537, 0.3516301943443956, -0.2286291611327847, 0.6347829113985894, 0.3803020118481135, -0.08981207028837489, -0.2709917210008783, -0.6182571092527998, -0.24690648656310987, -0.3168739742084861, -0.2045937444807841, -0.8983758456580146, 0.11048567978007713, 0.4432336003464372, -0.34463007653094696, 0.15019751004300907, 0.34828767120216864, -0.9840074502210592, -0.4199786243983529, 0.05841622225608528, -0.6922829985571581, 0.41230640939975893, 0.21536519638255888, -0.30444418105956184, 0.3986126839802461, -0.3859676200912196, 0.09636076023745846, -0.6353218161748273, 0.7260192504376082, -0.35967245966461714, -0.5179092834337204, 0.9920559419017999, -0.5951022815273046, 0.21664980147033097, -0.4609225690786949, -0.5598100841553817, -0.7380524343705033, 0.23339822550526618, 0.10635969533011985, -0.18433085980637193, 0.13546479326211963, -0.26367206630046636, -0.48017068783372097, -0.23978422697854446, 0.6180466757205332, 0.7022380629074493, -0.5126141794233905, 0.45914467478257115, -0.04732436529438497, -0.9715178045537647, 0.4707360940580251, 0.4248959605998251, 0.08696005864781456, -0.2280070003819314, -0.17484775251279622, -0.08847595688827967, -0.24458650749431984, -0.025297932066814008, 0.3347441931055716, 0.09265044303254766, 0.5781234828950566, 0.7340635138739218, -0.32437230221454955, -0.008426688593853055, -0.2538853704559949, 0.11607327349278454, 0.38411855842080594, 0.1397167588843119, -0.19063599169332276, 0.2556972055215963, 0.2578694668650857, -0.010880040826870205, -0.0021964920234608526, 0.39402450990967547, 0.4441385828736813, 0.20235903662018492, -0.20167313814206717, 0.17916927322294018, -0.3971742782759914, 0.242555853735351, 0.2089789563710438, 0.5565502608103653, -0.26156552970203706, -0.4575624465961444, 0.21436341340778184, -0.06619894525474732, 0.5033506634639126, -0.19278956639180203, 0.05421160222119592, -0.025508195892467098, 0.6502931155457493, 0.34158737022406405]
    val_dict["Liasis albertisii"] = [0.1338563720558809, 1.2909106895533127, 0.2225271954486219, 0.28414401847704673, -0.4821480223627333, 0.7879293558485955, 0.2583137413439375, 0.06263997371566238, 0.5483776094398648, 0.8093874315422102, -0.9294977479536474, 0.008454618480482692, 0.5011822305321163, 0.24728021578985204, 0.32747461152905, -0.5170517582196297, -0.12699067135975378, -0.8674333254961388, 0.23346564593136931, -0.31969612055820495, 0.28978870847775834, 1.058040269558188, -0.33416711347069916, 0.44253991506777046, 0.5674559086067751, -0.9298249698151855, -0.48466454371939993, -0.01164584172082922, -0.21013734721004992, 0.2637635050871464, 0.2213756113959211, -0.2605212754518249, -0.016410343739178324, -0.9249575610235341, -0.051512298520922595, -0.6567331399942367, 0.20168322158713997, 0.3369932328204918, -0.3998567868419581, 1.0029129825136676, -0.8859348886009459, 0.14571103362823035, -0.11605273862388604, -0.2730552390199467, -0.517956721144009, 0.07992722591989226, 0.18501493664318974, 0.4627454883223793, -0.2295526985289616, 0.17403520108670711, -0.6626941823080046, -0.2709983625780059, 0.5387235799818687, 0.7054860549973809, 0.0683098190675481, 0.4104795868221382, 0.5477087330831726, -0.7022813444361233, -0.11438465811233217, 0.3900819124736265, -0.4597479737144726, 0.49046898682969464, -0.4380186702799934, -0.045746912115448915, -0.206655947311818, 0.29374247278031723, 0.8375656575721964, 0.1574399411201867, 0.019442515310440306, 0.07588311230410694, -0.38654134488023834, 0.25538755329148516, -0.3610380080549943, -0.18356308957615725, -0.19315694257756957, -0.2540028879018247, -0.4471480500459226, 0.16042857685762738, 0.1112931140391353, 0.219558109796302, -0.13059986655144545, -0.28523735139330425, 0.2924993443868444, 0.4927264150256246, -0.36075194929241705, 0.5670193215831443, 0.10954171578704558, -0.16846002478390737, 0.24532955267406362, 1.303340352128946, 0.1615737409504714, -0.21058750403727083, 0.5505869299059208, -0.099837477163212, 0.22425483877359084, 0.10253510636711805, -0.48598525130759496, -0.28477687271518803, 0.3529792881941173, 0.346471263222687]
    val_dict["Bothrochilus boa"] = [-0.46584898422673415, 1.0211198644378225, 0.06229618453044433, 0.3589712911587379, -0.7832861224253536, 0.42019409069509006, 0.23948299893758151, -0.1058703460108053, 0.4239333065667278, -0.15159778929602763, -0.6855106397185754, -0.009168554180891589, 0.4040167113894405, -0.170083482456017, 0.004528876490866002, -0.6492248755245399, 0.32336706020405775, -0.595730400781251, -0.37344781628239176, -0.27295906891532634, -0.31621240706210246, 0.7587618715917581, -1.1576509094365166, 0.9074673315889867, 0.4101556905845832, -0.5990980265463374, -0.5269185944117846, -0.42951357262472073, 0.36860465333315406, 0.45353065649829105, -0.22934650011669827, -0.020808390010334937, 0.11725464051713491, -0.5617660252771319, 0.04792251004242884, -0.7496392596961657, 0.4734180530301766, -0.29930360765387587, -0.9809725099595799, 0.9802497312115896, -0.4792697856981302, 0.44507563778558434, -0.15003106114373022, -0.22793042689467397, -0.0919306150853951, 0.4353298766824597, 0.5373183403828572, 0.9120263367476178, 0.19678525379528536, -0.5757224975353232, -0.09322852640720664, -0.3370373272782041, 0.7706659582706259, 0.5522682419620584, 0.5636778148442024, 0.3964093158846792, 0.34406394855723577, -0.29980520243169717, 0.33722560108462746, 0.5119540722004793, -0.7531806522704954, 0.2824882730519933, 0.08625082896890884, -0.5647222068066741, -0.07612158570598404, -0.045931112360406944, 0.5186289077361061, 0.23290407438443056, 0.2745103026421708, 0.7371257286411994, -0.2999255943786685, -0.18471695984238107, -0.26983007950453713, 0.5192322942246012, -0.12935843320719173, 0.009146526710048036, -0.21540925208122744, 0.5519406463347722, -0.16086549489306523, 0.6427058123052194, 0.022749210123851577, -0.07127023059187784, 0.4068105229203746, 0.33562919656785106, -0.4515793139082209, 0.28962611697916585, -0.3412926259106599, 0.012940035813878994, -0.007286719210231052, 0.4686970037145938, -0.31319564763342617, 0.2353577115241414, -0.016775604047478302, -0.1428703340501273, -0.33067787871998666, 0.7883345989260337, -0.4315776043219617, -0.48998319107168375, 0.21576462250467737, 0.16696583739818582]
    val_dict["Morelia boeleni"] = [0.4495624904860971, 0.3815964914770747, 0.21637271179531564, 0.15444617556079288, -0.4513541292518066, 0.95450705435944, -0.6623743217543999, 0.15563575074580221, 0.7774339717234413, 0.8669297465566268, -0.3478442337847126, 0.6048550818351615, 0.14202231755017036, 0.32473472509658063, -0.3713878952421488, -1.314823957055675, 0.21900684203753037, -0.7096029724481906, -0.2764564311012366, -1.247472331117699, -0.9660461169610105, 0.1114976966772872, -0.5513304026928239, 0.3126361898969109, 0.5461011473434418, -0.288954333688125, -0.39887256523199033, 0.05249858012281258, -0.18949311992284265, 1.0847711993261222, 0.4193858657461046, 0.022641418713128317, -0.4196626685035914, -0.3006757839330384, 0.6015727644187413, -0.616492623954492, 0.619807509067704, 0.20894078552022377, -0.7368194690260614, 0.12388951898219913, -0.4883544856263653, 0.3246301626200404, -0.23503437117473497, -0.1384634060650008, -0.09533302273449495, 0.10866437624776112, 0.0042901153670109465, 0.12792214589343354, -0.674972589378737, -0.24807916515697312, -0.8396394442680646, -0.48185409942392793, -0.3946645867575628, 0.7610333692580763, -0.4353282458564297, 0.027420421240016757, -0.09510253190488882, -0.7405317470020447, 0.6328058342351901, 0.5453572496568836, 0.09575028816530912, 0.5590036889433486, -0.27045435550990615, 0.04354598216919278, 0.28614259704595735, 0.3282513290550281, 0.36614437989699666, 0.25370961033565875, 0.5590202951276149, 0.4291846318087165, -0.1657651546048891, 0.4423311233093086, -0.39897963681597537, -0.6882474322270227, -0.5834140239021541, -0.7168671214985056, 0.025102231036946865, 0.6852358705584138, -0.07874036134996823, -0.27374914317702026, -0.232868632104452, -0.2136966773272628, 0.3555238366134038, 0.486276045880913, 0.06061535070236322, 0.5919152267522342, -0.1805142566794374, 0.010212570270757615, 0.13724086475355174, 0.5132255486454709, -0.2885113908440649, -0.2212688529416054, 0.3625915352506578, 0.34931476595162325, -0.2051558592776081, 0.26085479159182434, -0.45355548062404044, 0.6402941190207132, 0.0841832467919536, -0.06904903716543082]
    val_dict["Python timoriensis"] = [-0.2993201952683183, -0.35830343878764426, 0.5116408222475933, 0.5305369642656759, -0.030973435176103327, 0.1583505614499111, -0.35510994794180156, -0.08351416075953814, 0.3929282717613925, 0.4676229075816749, -0.32382457007092336, 0.18897655795268706, -0.3919509550296981, 0.4123922893138603, -0.16591407229967253, -0.8383458163517846, 0.20974273200711782, -1.0396921143485527, 0.011387753387965271, 0.03649608189686693, -0.7300003829873433, 0.7847554100800485, 0.16915446479166707, -0.18801540776872444, -0.770962674757409, 0.03303062947128654, -0.7413563884254578, -0.2075716118720919, 0.3160057153489539, 0.19990817026899238, -0.1629169979857616, -0.19085665470747715, -0.48422044242204676, -0.49643034296798344, -0.39948557003269225, -0.929774667605282, 0.8468707615646343, 0.3495220590603937, -1.2318551919761225, 0.13880620507072822, -0.06128494520196698, 0.24891008593282435, -0.22721472307457558, 0.3945321442425566, -0.1485523673729286, 0.23497519875384926, -0.47939165716966403, -0.4289327751504128, -0.19428607990334373, 0.13714739209040167, -0.45707156087654544, -0.2625027700012786, 0.5211160781081379, 0.09563348719559811, 0.6542426208149114, 0.5145004262601014, 0.18789858702015727, 0.22138124690764355, 0.6321807115020192, 0.6246795557617644, 0.27467235891963293, -0.7205080099497054, 0.23695670878558925, -0.16856352763441282, -0.20086495866337756, -0.5680408511032506, 0.5072739187322157, 0.6455205606718868, 0.07869998971494935, 0.16250942189779233, 0.25553015927146505, 0.26007885078845794, 0.1774199536149157, -0.0141026324551689, -0.0636720476502875, 0.025433874396383338, 0.2351007196155972, -0.31121358255897447, 0.7316883075698194, 0.21071393298449845, 0.6817068458667898, 0.059829327046325265, -0.06169684369106204, 0.4101573168988367, -0.010314640292229837, 0.6772005665474999, -0.2817652690346034, 0.11507502908950215, 0.26520991398272, 0.39038396367824846, -0.1547378249648072, -0.10231875252422909, -0.27171914600486713, -0.11049465485041932, 0.04056466677452855, 0.4617333082080983, 0.33564222788548415, 0.3519082344623879, 0.5606699515351866, 0.27045261509119556]
    val_dict["Python reticulatus"] = [0.36053825556239366, -0.3996557593931853, 0.5276408995367758, 0.4453086141410851, -0.1731183379430418, 0.35320270279850835, 0.1410785492929132, -0.3665501035257141, 0.6667189150156314, 0.7663347174445609, -0.6452225877576345, 0.23093856287807768, 0.14675768625253505, 0.4307258436916486, -0.31556529553172585, -0.6671825116622025, -0.3785085766823774, -0.47567080534375566, -0.30866091063593415, -0.19399154534077712, -0.6696672571284988, 0.6750685559454414, 0.23916203149734178, 0.18379013803175548, -0.48431882379068014, -0.8325591807962615, -0.5065666358568384, -0.055834541235954625, 0.35797194179466413, 0.3056755576532644, 0.4057140492089184, -0.7076112117795632, -0.26081290370568827, -0.27645484674622567, -0.3839011164572983, -0.8901486169429839, 0.28965242008519687, 0.14613729493829913, -0.7030632906151676, -0.084815172652981, 0.14943757935329327, 0.05828422691780166, -0.4785418468809623, -0.41726637808122174, -0.020953026000723407, 3.932229070136017E-9, -0.7730810378032282, 0.24461036749267057, 0.4396523427384369, 0.015372341546792595, -0.7208757745620907, -0.043136623733328994, 0.3256626693570088, 0.2549516581302459, 0.43081579264701614, 0.5267961832941067, -0.17670091414319555, 0.19485016053765553, 0.3331389264815913, 0.4529241289673547, -0.28138400437847993, -0.3687803239995521, -0.05070183347080717, 0.021881503621266285, 0.3744105915367215, -0.8043729003868205, 0.03780046067975319, 0.6016587974796537, 0.3324560048737688, 0.2976302265434488, 0.43310703618410457, 0.9074247524397998, 0.02814284897020536, -0.06803528736753992, 0.3187821680515953, 0.010605168181892777, -0.23859828713068637, 0.15719449446518105, 1.077035381562189, 0.4659555238657302, 1.339808423127497, 0.008655343041064134, 0.25496833822438414, 0.20775177731602384, 0.3717570776454674, 0.4025177167181309, -0.10333271824769813, 0.2618446685182717, 0.08547681212459626, 0.19337595168090893, 0.07291742651610539, 0.2669641896054383, 0.13353020198857068, -0.17838663115838876, 0.2360598066976941, 0.0553815581514251, -0.2305352419380728, 0.891577915954662, 0.3131944491257824, -0.05460131032097429]
    val_dict["Xenopeltis unicolor"] = [-0.12664106005264758, 0.032903292119431715, -0.2975642099638368, -0.09485841088059496, 1.0140972412875509, -0.08397391347800681, -0.8612280000285408, 0.19499993244766595, 0.03607688720349113, -0.4052475437123401, 0.5568617290624406, 0.05781611014771472, 0.3107592868324674, 0.10167422165908602, 0.6082531040786614, 0.4743909708565066, 0.4366151098046645, 1.1930739210792687, -0.40636604781875707, 0.2124462464776434, -0.6691888810001302, 0.24028923183602796, 0.38801397954602335, -0.0031786761555590917, 0.21841132571698238, 0.3573166087190817, 0.00800800008631418, -0.23290590784922421, 0.3155144951139559, -1.449047945324792, 0.5727526333342419, -0.19544912823299632, -0.433945938379185, -0.06515545687842532, 0.5335380880239993, -0.16572333786595114, -0.5886997999448966, 0.12768068441040248, 0.2505951595864218, 0.2764464736291565, -0.08996854724011957, -0.3582471528972136, 0.05496233862667338, 0.4372575724209106, 0.46197266058285913, 0.3532827785625013, -0.24642184352183516, -0.34891776658359785, 0.9371830238921739, -0.6848750507493656, 0.42395621732374833, 0.5391840263119081, 0.11358122579520803, 0.6270619764631561, -0.5586491957079646, -0.02361167124070565, -0.2317031647583469, 0.07394342045562274, 0.26086209310746167, 0.3235067273469298, 0.197508295932402, 0.12493811956164128, 0.12008958851512316, -0.22074068326343893, 0.8138149952186611, 0.30233351672615166, -0.7759284377282983, -0.767660183032841, -0.69990821857353, 0.17351416794139846, 0.3617147249551747, 0.038912012090951735, -0.18127528850650432, -0.18475518113959072, 0.3284186629594207, -0.4803444004781864, -0.2572311859334843, 0.06517999213316311, 0.34192127399589484, 0.31868262143171633, -0.23065528757982273, 0.1952952202621504, -0.06200217797385384, 0.2189095779012651, -0.35807491522943774, -0.06544202112880163, 0.2447991237392162, 0.014678872832097886, -0.22142242012855284, -0.3087125895638537, -0.22106888711914494, 0.04884077769541521, 0.17918805308038308, -0.05900691528255396, -0.06796468480904705, -0.3069703030581442, 0.11318680151621358, -0.4280540562943652, -0.23412558086577967, -0.4330840337183335]
    val_dict["Candola aspera"] = [0.5415073439145329, -0.8261851752182238, 0.280374228361929, -1.0458119475720775, -0.5879476187175094, 0.5700143496596618, -0.8292641444852076, -0.6702694622550589, -0.9480901201022143, 0.33502021617202604, 0.6805009708608114, 0.15586856863635418, 0.7892489646811306, -1.1453557369820766, 0.36452062455213907, 1.16607256090316, 0.3953526215120202, 1.6684493806562049, 0.007996743558983055, -0.4819036225554018, 0.026559591460595944, 0.7352572544755298, 0.5906512380081679, 0.3307887896111144, 0.4778958810730767, -0.21822902282627166, -0.531956616578794, -0.12486253574358203, 0.10178018009526718, -1.0369031726485534, 1.0284204610799061, -1.2722975514473445, 0.8620084479071548, -0.11569624902603305, -0.7191078082719772, -0.899343756933466, 0.44130152278027224, -0.15572845260190246, -0.31908732876951, 0.3601623647259761, 0.20818887861984042, -0.17901884053495137, 0.19913207776418712, -0.28585811979075954, 0.7100003485577421, 0.15904064834949316, 0.44181645392713426, -1.3352679486697252, 0.11002396926434246, -0.08057061752566023, 0.7278491628895039, 0.4961433359529824, 0.26448506279377204, -0.3750879236157483, -0.17158111107259355, 0.7164493380433936, -0.04663475174269294, -1.0613670475874053, 0.5164015353746929, -0.3284568109962387, 0.32100832256571316, -0.23273112120667738, -0.8454346791189112, -0.11843661482725729, 0.25958614260850615, 0.3168475351394587, -0.5737954792276314, 0.7330343625512763, 0.20065983974677323, 0.44004347514617265, 0.13544629554463464, 0.0526790338053572, -0.6171058144117554, 0.9346687975574848, 0.8149920923948215, -0.10195544966143313, 0.10628599636940993, 0.513322825545871, -0.002869011867480231, -1.3800004773929637, 0.32168833317412837, 0.053923037760593, -0.2272100361040466, 0.04968323055576585, -0.046484592583790396, 1.29869175242832, 1.0113826728236226, -0.4819390754660011, 0.517698551084355, -0.6916519124931837, -1.6619561957545632, 0.7634336590013959, -0.023458656557591218, 1.6501522159310666, -1.042169009622521, 0.2102005354339974, -1.9142014377267351, -0.12223995948995559, 0.025274502076111066, 0.6752210375038439]
    val_dict["Loxocemus bicolor"] = [-0.48582867963460824, 0.5561677763841473, 0.17458074714015517, 0.385053440786788, -0.01136784021396468, 0.14234444644739344, 0.5596480217143147, -0.18704884299254537, -0.34727537883051995, 0.6280451374869649, -0.49004203146040265, -0.25601854985973743, 0.013125321179269893, -0.045253108716276136, -1.043071151609052, 0.23750821694679286, -1.461426784929922, -0.7363587867138378, -0.8956761285626337, -0.3186875130659576, -0.365186396948155, 0.4309709036972196, -0.7874881193590557, -0.3071991866432301, 0.2656581195072468, -0.652382901754644, -0.5050261697993766, 0.02892116680076118, 0.47897322157257805, -0.3866140935241261, 0.40375917395157845, -0.8412943755940032, -0.33746571813159637, -0.8861676606913741, 0.3854716598117076, -0.263540860897671, -0.22409483173995542, 0.5116153366018729, -0.29466962436532856, -0.06044962816591291, -1.3530962136410358, 0.6663578607799884, -0.01643874341404085, 0.34488420881733756, -0.2948416741073028, 0.30935481350796995, 0.2179794475089686, -0.21246231946534694, 0.04457495561461377, 0.24863419519216778, -0.16838005417564145, -0.30447206584727393, 1.064132048787358, -0.17146894748012473, 0.4615188958276311, 0.4524310656704919, 0.3921272551693488, 0.20959524529756857, 0.42416031944879284, 0.8040531560249723, -0.2489623887225627, -0.3465257340592563, -0.7119106289902093, 0.3502391152482162, -0.39896574569713894, 0.023136181148250293, -0.4333099378574654, -0.09307312305902067, 0.08029022826963969, -0.5395373285896503, 0.9740929814793209, -0.029423015088022232, -0.4266813690703312, -1.1969795093341153, -0.11692905615539284, 0.20846505568479645, 0.32076473938523875, -0.42146251375819205, 0.087655482278029, 0.316767110244248, 0.4605433380142203, -0.1723002176024076, -0.6502998758026927, 0.18181951081107403, -1.242457947293979, 0.8152265947468574, 0.45329419203745924, -0.5064672242877676, 0.39255582013311513, -0.20021681063959576, 0.1306008971866452, 0.41314982676479584, -0.13635375564364613, 1.1236931264065555, 0.34594618178840486, 0.25961282630301724, 0.530198507464853, 0.5900915701480034, 0.670247130924944, 0.08122391648367935]
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
        __init__ creates an instance with the optional argument `x` controls seeding,
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
        #_LOG.debug("RepeatedRandom returning %f" % r)
        return r
