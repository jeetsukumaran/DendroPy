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

def reference_trees_filename(format):
    if format == "nexus":
        return "pythonidae.reference-trees.nexus"
    elif format == "newick":
        return "pythonidae.reference-trees.newick"
    else:
        raise ValueError("Reference trees not available in '%s' format" % format)

def reference_tree_list_postorder_node_labels():
    return [
        ['Morelia oenpelliensis','Morelia spilota','Node36995248','Apodora papuana','Aspidites melanocephalus','Morelia carinata','Node36995920','Node36995600','Bothrochilus boa','Morelia amethistina','Node36995952','Node36995536','Morelia nauta','Node36995472','Morelia tracyae','Python regius','Node36996400','Node36995440','Python reticulatus','Python brongersmai','Liasis olivaceus','Node36996880','Node36996464','Python timoriensis','Node36996560','Morelia clastolepis','Antaresia maculosa','Node36997072','Liasis fuscus','Liasis mackloti','Node36997360','Node36997168','Node36996176','Morelia boeleni','Node36996528','Node36995344','Python molurus','Antaresia stimsoni','Python sebae','Node36997712','Morelia kinghorni','Node36997840','Node36997616','Aspidites ramsayi','Antaresia childreni','Leiopython albertisii','Node36998032','Antaresia perthensis','Node36998320','Node36997968','Morelia viridis','Morelia bredli','Node36998704','Node36998064','Node36997392','Node36962224','Node36995184'],
        ['Morelia boeleni','Leiopython albertisii','Bothrochilus boa','Node36999024','Morelia carinata','Morelia spilota','Morelia bredli','Node37061072','Antaresia maculosa','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37061584','Node37061456','Node37061392','Node37060944','Morelia viridis','Node37061040','Node37060720','Node36962160','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37062064','Node37061872','Apodora papuana','Node37061744','Node36995152','Morelia oenpelliensis','Morelia amethistina','Morelia clastolepis','Morelia kinghorni','Morelia nauta','Node37062896','Node37062768','Node37062544','Morelia tracyae','Node37062448','Python reticulatus','Python timoriensis','Node37063184','Node37062512','Node37062288','Python sebae','Python regius','Python brongersmai','Python molurus','Node37063664','Node37063536','Node37063248','Aspidites ramsayi','Aspidites melanocephalus','Node37063760','Node37063088','Node37062192','Node36999056','Node36995120'],
        ['Python sebae','Python molurus','Python brongersmai','Node37064432','Python regius','Node37064496','Node37064048','Aspidites ramsayi','Aspidites melanocephalus','Node37654640','Node37063696','Morelia carinata','Morelia spilota','Morelia bredli','Node37655120','Node37064656','Antaresia maculosa','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37655536','Node37655408','Node37655152','Node37654896','Morelia viridis','Node37654704','Node36998896','Leiopython albertisii','Bothrochilus boa','Node37655728','Node37064208','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37656272','Node37655632','Apodora papuana','Node37655984','Morelia boeleni','Morelia amethistina','Morelia nauta','Morelia kinghorni','Node37656912','Morelia clastolepis','Node37656976','Node37656304','Morelia tracyae','Node37656688','Python reticulatus','Python timoriensis','Node37657392','Node37656528','Morelia oenpelliensis','Node37656656','Node37656560','Node37656048','Node37064144'],
        ['Python sebae','Python regius','Python molurus','Python brongersmai','Node37658128','Node37658000','Node36962032','Aspidites ramsayi','Aspidites melanocephalus','Node37658224','Node37657776','Python reticulatus','Python timoriensis','Node37658160','Morelia amethistina','Morelia nauta','Morelia kinghorni','Morelia clastolepis','Node37712496','Node37712368','Node37711984','Morelia tracyae','Node37712016','Morelia boeleni','Node37658576','Node37658288','Morelia oenpelliensis','Node37657648','Morelia viridis','Antaresia maculosa','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37713488','Node37713360','Node37713072','Morelia spilota','Morelia bredli','Node37713648','Node37713200','Node37712400','Morelia carinata','Node37712816','Leiopython albertisii','Bothrochilus boa','Node37714000','Node37712848','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37714352','Node37714128','Apodora papuana','Node37714192','Node37712912','Node37658384','Node37657712'],
        ['Python reticulatus','Python timoriensis','Node37714768','Morelia viridis','Morelia carinata','Morelia spilota','Morelia bredli','Node37715504','Node37715184','Antaresia maculosa','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37715760','Node37715792','Node37715536','Node37715312','Node37714992','Leiopython albertisii','Bothrochilus boa','Node37715600','Node37714960','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37773840','Node37773616','Apodora papuana','Node37773680','Node37714512','Morelia amethistina','Morelia clastolepis','Morelia nauta','Morelia kinghorni','Node37774608','Node37774480','Node37774160','Morelia tracyae','Node37773904','Morelia boeleni','Node37774032','Morelia oenpelliensis','Node37774128','Python sebae','Python regius','Python brongersmai','Python molurus','Node37775312','Node37775184','Node37774928','Aspidites ramsayi','Aspidites melanocephalus','Node37775408','Node37775024','Node37773872','Node37715024','Node37714704'],
        ['Antaresia maculosa','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37776432','Node37776304','Node37776016','Morelia viridis','Morelia spilota','Morelia bredli','Node37776784','Node37776144','Node37775952','Morelia carinata','Node37775472','Leiopython albertisii','Bothrochilus boa','Node37777072','Node37775568','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37777360','Node37777232','Apodora papuana','Node37777264','Node37064112','Morelia amethistina','Morelia clastolepis','Morelia kinghorni','Morelia nauta','Node37831344','Node37831216','Node37830928','Morelia tracyae','Node37777040','Morelia boeleni','Node37830704','Morelia oenpelliensis','Node37830864','Node37775856','Python reticulatus','Python timoriensis','Node37831504','Python sebae','Python brongersmai','Python molurus','Node37832048','Python regius','Node37832176','Node37831888','Aspidites ramsayi','Aspidites melanocephalus','Node37832400','Node37832016','Node37831600','Node37775792'],
        ['Python sebae','Python molurus','Python regius','Python brongersmai','Node37833136','Node37833008','Node37832848','Python reticulatus','Python timoriensis','Node37833168','Morelia viridis','Morelia spilota','Morelia bredli','Node37833904','Morelia carinata','Node37834032','Node37833520','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37834512','Node37834288','Antaresia maculosa','Node37834256','Node37833680','Leiopython albertisii','Bothrochilus boa','Node37834736','Node37833456','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37101968','Node37101808','Apodora papuana','Node37101712','Node37833744','Node37833328','Morelia oenpelliensis','Morelia tracyae','Morelia amethistina','Morelia nauta','Morelia kinghorni','Node37102608','Morelia clastolepis','Node37102672','Node37102544','Node37102288','Morelia boeleni','Node37102352','Node37101936','Node37833360','Aspidites ramsayi','Aspidites melanocephalus','Node37103024','Node37833200','Node37832784'],
        ['Morelia tracyae','Morelia nauta','Morelia clastolepis','Node37103920','Morelia kinghorni','Node37103824','Morelia amethistina','Node37103888','Node37103568','Morelia boeleni','Node37103504','Apodora papuana','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37104720','Node37104592','Node37104336','Leiopython albertisii','Bothrochilus boa','Node37104752','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37105392','Node37105200','Antaresia maculosa','Node37104880','Morelia viridis','Morelia carinata','Morelia spilota','Morelia bredli','Node37163248','Node37105232','Node37105584','Node37104976','Node37104784','Node37104432','Node37103120','Python reticulatus','Python timoriensis','Node37163312','Morelia oenpelliensis','Node37103280','Node37102896','Aspidites ramsayi','Aspidites melanocephalus','Node37163504','Node37714672','Python molurus','Python brongersmai','Python regius','Node37164080','Node37163696','Node37103408','Python sebae','Node37103344'],
        ['Morelia viridis','Morelia spilota','Morelia bredli','Node37164880','Morelia carinata','Node37164720','Antaresia maculosa','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37165456','Node37165328','Node37165264','Node37164848','Node37164304','Python reticulatus','Python timoriensis','Node37165616','Apodora papuana','Leiopython albertisii','Bothrochilus boa','Node37166192','Node37165808','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37166512','Node37166384','Node37165840','Morelia tracyae','Morelia kinghorni','Morelia nauta','Node37166992','Morelia clastolepis','Node37166896','Morelia amethistina','Node37166960','Node37166704','Morelia boeleni','Node37166544','Morelia oenpelliensis','Node37166576','Node37165968','Node37165712','Node37775760','Aspidites ramsayi','Aspidites melanocephalus','Node37220432','Node37164368','Python regius','Python brongersmai','Node37220944','Python molurus','Node37220976','Node37164496','Python sebae','Node37164432'],
        ['Python sebae','Python regius','Python molurus','Python brongersmai','Node37221776','Node37221616','Morelia oenpelliensis','Morelia tracyae','Morelia kinghorni','Morelia nauta','Node37222448','Morelia clastolepis','Node37222288','Morelia amethistina','Node37222416','Node37222096','Morelia boeleni','Node37222224','Node37221872','Liasis olivaceus','Liasis fuscus','Liasis mackloti','Node37223184','Node37222864','Python reticulatus','Python timoriensis','Node37223376','Apodora papuana','Node37223248','Node37222960','Leiopython albertisii','Bothrochilus boa','Node37223408','Morelia carinata','Morelia spilota','Morelia bredli','Node37224208','Node37223888','Morelia viridis','Node37223728','Antaresia maculosa','Antaresia perthensis','Antaresia childreni','Antaresia stimsoni','Node37277904','Node37224368','Node37224304','Node37223760','Node37223632','Node37222736','Node37221968','Aspidites ramsayi','Aspidites melanocephalus','Node37224432','Node37221936','Node37221584','Node37221424'],
    ]

def reference_tree_list_newick_string():
    return """\
        (('Morelia oenpelliensis':5.50104114707,'Morelia spilota':5.50104114707)Node36995248:61.9915926669,(((((('Apodora papuana':0.707490222359,('Aspidites melanocephalus':0.236984491852,'Morelia carinata':0.236984491852)Node36995920:0.470505730507)Node36995600:2.09759626003,('Bothrochilus boa':1.06956986931,'Morelia amethistina':1.06956986931)Node36995952:1.73551661307)Node36995536:0.351619575146,'Morelia nauta':3.15670605753)Node36995472:1.09375666234,('Morelia tracyae':0.585328937986,'Python regius':0.585328937986)Node36996400:3.66513378189)Node36995440:33.0195611855,(((('Python reticulatus':1.32521964733,('Python brongersmai':0.626131270775,'Liasis olivaceus':0.626131270775)Node36996880:0.699088376552)Node36996464:3.58364386958,'Python timoriensis':4.9088635169)Node36996560:4.90381799328,(('Morelia clastolepis':0.520293181797,'Antaresia maculosa':0.520293181797)Node36997072:3.89664227731,('Liasis fuscus':0.193228819804,'Liasis mackloti':0.193228819804)Node36997360:4.2237066393)Node36997168:5.39574605108)Node36996176:3.49206820101,'Morelia boeleni':13.3047497112)Node36996528:23.9652741942)Node36995344:0.377327561902,(('Python molurus':12.8622487575,(('Antaresia stimsoni':1.20132174738,'Python sebae':1.20132174738)Node36997712:3.83058430268,'Morelia kinghorni':5.03190605005)Node36997840:7.83034270743)Node36997616:21.1685275644,(('Aspidites ramsayi':2.78732559217,(('Antaresia childreni':0.105428136944,'Leiopython albertisii':0.105428136944)Node36998032:1.30019225223,'Antaresia perthensis':1.40562038917)Node36998320:1.381705203)Node36997968:4.80697376221,('Morelia viridis':2.42117796541,'Morelia bredli':2.42117796541)Node36998704:5.17312138897)Node36998064:26.4364769675)Node36997392:3.61657514544)Node36962224:29.8452823467)Node36995184;
        ('Morelia boeleni':26.0384706379,(((('Leiopython albertisii':15.2611605145,'Bothrochilus boa':15.2611605145)Node36999024:2.95693808858,('Morelia carinata':16.0625616183,((('Morelia spilota':7.05059979983,'Morelia bredli':7.05059979983)Node37061072:7.4747693113,('Antaresia maculosa':13.0854581588,('Antaresia perthensis':10.2267471766,('Antaresia childreni':4.8707020453,'Antaresia stimsoni':4.8707020453)Node37061584:5.35604513133)Node37061456:2.85871098218)Node37061392:1.43991095233)Node37060944:0.166805015451,'Morelia viridis':14.6921741266)Node37061040:1.37038749174)Node37060720:2.15553698471)Node36962160:1.64883491058,(('Liasis olivaceus':10.850645624,('Liasis fuscus':2.83048300382,'Liasis mackloti':2.83048300382)Node37062064:8.02016262021)Node37061872:1.94094943657,'Apodora papuana':12.7915950606)Node37061744:7.07533845302)Node36995152:0.952265882909,(('Morelia oenpelliensis':19.051793087,((('Morelia amethistina':8.12042923083,('Morelia clastolepis':2.0657194471,('Morelia kinghorni':0.999651724272,'Morelia nauta':0.999651724272)Node37062896:1.06606772283)Node37062768:6.05470978372)Node37062544:1.9541248425,'Morelia tracyae':10.0745540733)Node37062448:7.19864380037,('Python reticulatus':6.6158514927,'Python timoriensis':6.6158514927)Node37063184:10.657346381)Node37062512:1.77859521326)Node37062288:1.42657874327,(('Python sebae':14.3734153538,('Python regius':11.7636662419,('Python brongersmai':10.5318200836,'Python molurus':10.5318200836)Node37063664:1.23184615824)Node37063536:2.6097491119)Node37063248:3.93324105164,('Aspidites ramsayi':9.72719847603,'Aspidites melanocephalus':9.72719847603)Node37063760:8.5794579294)Node37063088:2.1717154248)Node37062192:0.340827566292)Node36999056:5.21927124134)Node36995120;
        ((((('Python sebae':16.9224723285,(('Python molurus':13.9521024428,'Python brongersmai':13.9521024428)Node37064432:0.218628811139,'Python regius':14.170731254)Node37064496:2.75174107457)Node37064048:2.97395412534,('Aspidites ramsayi':7.82754382801,'Aspidites melanocephalus':7.82754382801)Node37654640:12.0688826259)Node37063696:3.261105684,((('Morelia carinata':11.7896562791,('Morelia spilota':4.05445265776,'Morelia bredli':4.05445265776)Node37655120:7.73520362135)Node37064656:2.97877570893,('Antaresia maculosa':13.5106259111,('Antaresia perthensis':11.2360265701,('Antaresia childreni':3.48847926746,'Antaresia stimsoni':3.48847926746)Node37655536:7.74754730265)Node37655408:2.27459934102)Node37655152:1.25780607691)Node37654896:2.9034791189,'Morelia viridis':17.6719111069)Node37654704:5.48562103095)Node36998896:1.06929723627,('Leiopython albertisii':13.3580917481,'Bothrochilus boa':13.3580917481)Node37655728:10.8687376261)Node37064208:0.0364728670888,((('Liasis olivaceus':6.21487951846,('Liasis fuscus':2.60536108558,'Liasis mackloti':2.60536108558)Node37656272:3.60951843288)Node37655632:1.91471562327,'Apodora papuana':8.12959514173)Node37655984:14.878020141,('Morelia boeleni':20.4730397427,(((('Morelia amethistina':4.43324641637,(('Morelia nauta':0.811519404534,'Morelia kinghorni':0.811519404534)Node37656912:0.622029416262,'Morelia clastolepis':1.4335488208)Node37656976:2.99969759558)Node37656304:3.28625146667,'Morelia tracyae':7.71949788304)Node37656688:9.88501788667,('Python reticulatus':5.78102399509,'Python timoriensis':5.78102399509)Node37657392:11.8234917746)Node37656528:0.296424924481,'Morelia oenpelliensis':17.9009406942)Node37656656:2.57209904849)Node37656560:2.53457554003)Node37656048:1.25568695854)Node37064144;
        ((('Python sebae':16.5921325784,('Python regius':14.3799352093,('Python molurus':13.3821288803,'Python brongersmai':13.3821288803)Node37658128:0.997806328959)Node37658000:2.21219736911)Node36962032:3.43039474792,('Aspidites ramsayi':11.703431691,'Aspidites melanocephalus':11.703431691)Node37658224:8.31909563528)Node37657776:3.09272413642,(((('Python reticulatus':6.46326841791,'Python timoriensis':6.46326841791)Node37658160:13.3467104377,((('Morelia amethistina':4.67092080519,('Morelia nauta':2.08880085866,('Morelia kinghorni':1.48123046362,'Morelia clastolepis':1.48123046362)Node37712496:0.607570395039)Node37712368:2.58211994653)Node37711984:1.71754130617,'Morelia tracyae':6.38846211136)Node37712016:3.03031595643,'Morelia boeleni':9.41877806779)Node37658576:10.3912007878)Node37658288:0.121667790832,'Morelia oenpelliensis':19.9316466464)Node37657648:2.71760756257,(((('Morelia viridis':15.9488284198,(('Antaresia maculosa':14.0002983201,('Antaresia perthensis':11.4477602235,('Antaresia childreni':5.72844294344,'Antaresia stimsoni':5.72844294344)Node37713488:5.71931728003)Node37713360:2.5525380966)Node37713072:1.55388791109,('Morelia spilota':4.49531637063,'Morelia bredli':4.49531637063)Node37713648:11.0588698605)Node37713200:0.394642188595)Node37712400:2.23472642535,'Morelia carinata':18.1835548451)Node37712816:1.84844364601,('Leiopython albertisii':12.7456406695,'Bothrochilus boa':12.7456406695)Node37714000:7.28635782162)Node37712848:1.89223670238,(('Liasis olivaceus':13.5304034382,('Liasis fuscus':1.74454580066,'Liasis mackloti':1.74454580066)Node37714352:11.7858576376)Node37714128:1.71076794287,'Apodora papuana':15.2411713811)Node37714192:6.6830638124)Node37712912:0.725019015485)Node37658384:0.465997253751)Node37657712;
        (('Python reticulatus':10.8300265774,'Python timoriensis':10.8300265774)Node37714768:17.6343232815,(((('Morelia viridis':14.0077690635,(('Morelia carinata':11.7785955814,('Morelia spilota':4.97579770159,'Morelia bredli':4.97579770159)Node37715504:6.80279787983)Node37715184:1.82500345118,('Antaresia maculosa':12.3274764924,('Antaresia perthensis':8.00332145017,('Antaresia childreni':4.01522437806,'Antaresia stimsoni':4.01522437806)Node37715760:3.98809707211)Node37715792:4.32415504222)Node37715536:1.27612254021)Node37715312:0.404170030911)Node37714992:4.88508122224,('Leiopython albertisii':11.1450426742,'Bothrochilus boa':11.1450426742)Node37715600:7.74780761159)Node37714960:1.31772152158,(('Liasis olivaceus':10.9742846954,('Liasis fuscus':3.16669963676,'Liasis mackloti':3.16669963676)Node37773840:7.80758505868)Node37773616:3.78712208743,'Apodora papuana':14.7614067829)Node37773680:5.44916502448)Node37714512:4.52906062503,((((('Morelia amethistina':5.50319615164,('Morelia clastolepis':1.57257658111,('Morelia nauta':1.4357285473,'Morelia kinghorni':1.4357285473)Node37774608:0.136848033805)Node37774480:3.93061957053)Node37774160:4.50938662462,'Morelia tracyae':10.0125827763)Node37773904:6.85886795982,'Morelia boeleni':16.8714507361)Node37774032:5.48583318771,'Morelia oenpelliensis':22.3572839238)Node37774128:1.8069504389,(('Python sebae':16.535465064,('Python regius':14.962419097,('Python brongersmai':13.1354320851,'Python molurus':13.1354320851)Node37775312:1.82698701194)Node37775184:1.573045967)Node37774928:4.52064732483,('Aspidites ramsayi':14.6640030997,'Aspidites melanocephalus':14.6640030997)Node37775408:6.39210928917)Node37775024:3.10812197386)Node37773872:0.57539806967)Node37715024:3.72471742653)Node37714704;
        ((((((('Antaresia maculosa':12.6196623606,('Antaresia perthensis':10.1068643474,('Antaresia childreni':4.39318906371,'Antaresia stimsoni':4.39318906371)Node37776432:5.71367528366)Node37776304:2.51279801319)Node37776016:1.83427060068,('Morelia viridis':13.41617982,('Morelia spilota':4.03982713639,'Morelia bredli':4.03982713639)Node37776784:9.37635268364)Node37776144:1.03775314121)Node37775952:0.835745557531,'Morelia carinata':15.2896785188)Node37775472:1.65141371467,('Leiopython albertisii':12.7629580607,'Bothrochilus boa':12.7629580607)Node37777072:4.17813417278)Node37775568:0.754225130388,(('Liasis olivaceus':10.6502403013,('Liasis fuscus':3.67165568794,'Liasis mackloti':3.67165568794)Node37777360:6.97858461332)Node37777232:4.06497756637,'Apodora papuana':14.7152178676)Node37777264:2.9800994962)Node37064112:4.09971291135,(((('Morelia amethistina':4.6447123226,('Morelia clastolepis':2.17174905761,('Morelia kinghorni':0.904095798097,'Morelia nauta':0.904095798097)Node37831344:1.26765325952)Node37831216:2.47296326498)Node37830928:4.27056899374,'Morelia tracyae':8.91528131634)Node37777040:3.1657791763,'Morelia boeleni':12.0810604926)Node37830704:6.58446507703,'Morelia oenpelliensis':18.6655255697)Node37830864:3.12950470552)Node37775856:1.01117536037,(('Python reticulatus':11.5585932535,'Python timoriensis':11.5585932535)Node37831504:11.1718400909,(('Python sebae':15.9676320675,(('Python brongersmai':13.8599706405,'Python molurus':13.8599706405)Node37832048:1.04487219988,'Python regius':14.9048428404)Node37832176:1.06278922707)Node37831888:5.31414730047,('Aspidites ramsayi':15.5786974573,'Aspidites melanocephalus':15.5786974573)Node37832400:5.70308191065)Node37832016:1.44865397646)Node37831600:0.075772291124)Node37775792;
        (('Python sebae':29.0711850036,('Python molurus':22.4124826989,('Python regius':21.9264827074,'Python brongersmai':21.9264827074)Node37833136:0.485999991488)Node37833008:6.65870230472)Node37832848:2.12261395029,(((('Python reticulatus':14.5492404256,'Python timoriensis':14.5492404256)Node37833168:5.89251749724,(((('Morelia viridis':14.834011977,(('Morelia spilota':6.2738399067,'Morelia bredli':6.2738399067)Node37833904:7.33173338162,'Morelia carinata':13.6055732883)Node37834032:1.22843868872)Node37833520:0.322649062335,(('Antaresia perthensis':10.1091415386,('Antaresia childreni':3.48320079048,'Antaresia stimsoni':3.48320079048)Node37834512:6.62594074811)Node37834288:4.22599220252,'Antaresia maculosa':14.3351337411)Node37834256:0.821527298272)Node37833680:2.93802844014,('Leiopython albertisii':9.28422587488,'Bothrochilus boa':9.28422587488)Node37834736:8.81046360463)Node37833456:0.826870765744,(('Liasis olivaceus':11.0151449224,('Liasis fuscus':4.0880111593,'Liasis mackloti':4.0880111593)Node37101968:6.92713376313)Node37101808:2.70301025103,'Apodora papuana':13.7181551735)Node37101712:5.2034050718)Node37833744:1.52019767757)Node37833328:1.10888028588,('Morelia oenpelliensis':17.6009445592,(('Morelia tracyae':5.80279597323,('Morelia amethistina':4.38107217442,(('Morelia nauta':1.8029155441,'Morelia kinghorni':1.8029155441)Node37102608:0.0774488112983,'Morelia clastolepis':1.8803643554)Node37102672:2.50070781902)Node37102544:1.42172379881)Node37102288:6.37126259081,'Morelia boeleni':12.174058564)Node37102352:5.4268859952)Node37101936:3.94969364946)Node37833360:1.24633634472,('Aspidites ramsayi':9.84704548855,'Aspidites melanocephalus':9.84704548855)Node37103024:12.9499290649)Node37833200:8.39682440052)Node37832784;
        ((((((('Morelia tracyae':8.31432086046,((('Morelia nauta':0.912676316658,'Morelia clastolepis':0.912676316658)Node37103920:0.71440797492,'Morelia kinghorni':1.62708429158)Node37103824:2.88132393295,'Morelia amethistina':4.50840822453)Node37103888:3.80591263593)Node37103568:4.13809660753,'Morelia boeleni':12.452417468)Node37103504:9.50779575096,(('Apodora papuana':17.2449081924,('Liasis olivaceus':12.2505582513,('Liasis fuscus':6.29219869082,'Liasis mackloti':6.29219869082)Node37104720:5.95835956047)Node37104592:4.9943499411)Node37104336:1.50938344234,(('Leiopython albertisii':13.0633246673,'Bothrochilus boa':13.0633246673)Node37104752:5.19337056498,((('Antaresia perthensis':10.3499331345,('Antaresia childreni':3.38573289378,'Antaresia stimsoni':3.38573289378)Node37105392:6.96420024077)Node37105200:3.43791705079,'Antaresia maculosa':13.7878501853)Node37104880:2.06371429021,('Morelia viridis':15.2096070352,('Morelia carinata':13.6526925966,('Morelia spilota':4.92163119076,'Morelia bredli':4.92163119076)Node37163248:8.73106140586)Node37105232:1.5569144386)Node37105584:0.641957440329)Node37104976:2.40513075673)Node37104784:0.497596402447)Node37104432:3.20592158422)Node37103120:0.983173088098,(('Python reticulatus':14.7952666012,'Python timoriensis':14.7952666012)Node37163312:6.61343105392,'Morelia oenpelliensis':21.4086976551)Node37103280:1.53468865192)Node37102896:4.2262164672,('Aspidites ramsayi':6.63522564694,'Aspidites melanocephalus':6.63522564694)Node37163504:20.5343771273)Node37714672:6.41982933271,('Python molurus':23.4003588995,('Python brongersmai':23.3469263821,'Python regius':23.3469263821)Node37164080:0.0534325174359)Node37163696:10.1890732074)Node37103408:8.43656365605,'Python sebae':42.025995763)Node37103344;
        ((((('Morelia viridis':14.5888315231,((('Morelia spilota':4.61249165082,'Morelia bredli':4.61249165082)Node37164880:5.06534441374,'Morelia carinata':9.67783606456)Node37164720:4.89600409139,('Antaresia maculosa':14.05750906,('Antaresia perthensis':7.60849369834,('Antaresia childreni':1.86358271521,'Antaresia stimsoni':1.86358271521)Node37165456:5.74491098314)Node37165328:6.44901536169)Node37165264:0.516331095919)Node37164848:0.0149913671254)Node37164304:2.82624159524,(('Python reticulatus':8.60012281197,'Python timoriensis':8.60012281197)Node37165616:8.7380438705,((('Apodora papuana':13.0669390597,('Leiopython albertisii':9.33840615844,'Bothrochilus boa':9.33840615844)Node37166192:3.7285329013)Node37165808:1.38724843406,('Liasis olivaceus':10.9703192915,('Liasis fuscus':7.13435131963,'Liasis mackloti':7.13435131963)Node37166512:3.83596797183)Node37166384:3.48386820234)Node37165840:1.98993837612,((('Morelia tracyae':8.32099940228,((('Morelia kinghorni':1.07398430658,'Morelia nauta':1.07398430658)Node37166992:0.774826953003,'Morelia clastolepis':1.84881125958)Node37166896:2.02523285213,'Morelia amethistina':3.87404411171)Node37166960:4.44695529056)Node37166704:1.7419442928,'Morelia boeleni':10.0629436951)Node37166544:3.52728125011,'Morelia oenpelliensis':13.5902249452)Node37166576:2.85390092474)Node37165968:0.894040812546)Node37165712:0.0769064358513)Node37775760:4.09158345085,('Aspidites ramsayi':9.1538829568,'Aspidites melanocephalus':9.1538829568)Node37220432:12.3527736124)Node37164368:3.29007229288,(('Python regius':22.6656294361,'Python brongersmai':22.6656294361)Node37220944:0.5196540118,'Python molurus':23.1852834479)Node37220976:1.61144541412)Node37164496:4.51223406477,'Python sebae':29.3089629268)Node37164432;
        ('Python sebae':41.037475788,(('Python regius':32.1860160849,('Python molurus':28.8078489217,'Python brongersmai':28.8078489217)Node37221776:3.37816716317)Node37221616:7.88042810387,((('Morelia oenpelliensis':20.2736124855,(('Morelia tracyae':8.92774508942,((('Morelia kinghorni':3.15335248923,'Morelia nauta':3.15335248923)Node37222448:0.318849288103,'Morelia clastolepis':3.47220177733)Node37222288:4.46813855011,'Morelia amethistina':7.94034032744)Node37222416:0.98740476198)Node37222096:3.53762335708,'Morelia boeleni':12.4653684465)Node37222224:7.80824403904)Node37221872:3.04566241261,((('Liasis olivaceus':14.4158757568,('Liasis fuscus':4.30247792095,'Liasis mackloti':4.30247792095)Node37223184:10.1133978359)Node37222864:6.09401457319,(('Python reticulatus':11.1261504133,'Python timoriensis':11.1261504133)Node37223376:9.21506648422,'Apodora papuana':20.3412168975)Node37223248:0.168673432469)Node37222960:1.17144499717,(('Leiopython albertisii':16.554642571,'Bothrochilus boa':16.554642571)Node37223408:4.80149248356,((('Morelia carinata':13.6491727777,('Morelia spilota':5.95375576867,'Morelia bredli':5.95375576867)Node37224208:7.69541700903)Node37223888:0.59310677926,'Morelia viridis':14.242279557)Node37223728:1.51569357225,('Antaresia maculosa':13.2829808571,('Antaresia perthensis':9.80844705042,('Antaresia childreni':5.81147051386,'Antaresia stimsoni':5.81147051386)Node37277904:3.99697653657)Node37224368:3.47453380669)Node37224304:2.47499227209)Node37223760:5.59816192539)Node37223632:0.325200272563)Node37222736:1.63793957099)Node37221968:5.61233091052,('Aspidites ramsayi':9.3156335341,'Aspidites melanocephalus':9.3156335341)Node37224432:19.6159722746)Node37221936:11.1348383801)Node37221584:0.971031599202)Node37221424;
    """

def reference_tree_list_node_relationships():
    treelist_node_references = [
        {
            'Node36995184' : NodeRelationship(parent_label=None, child_labels=['Node36995248','Node36962224'], edge_length=None, taxon_label=None),
            'Node36995248' : NodeRelationship(parent_label='Node36995184', child_labels=['Morelia oenpelliensis','Morelia spilota'], edge_length=61.9915926669, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node36995248', child_labels=[], edge_length=5.50104114707, taxon_label='Morelia oenpelliensis'),
            'Morelia spilota' : NodeRelationship(parent_label='Node36995248', child_labels=[], edge_length=5.50104114707, taxon_label='Morelia spilota'),
            'Node36962224' : NodeRelationship(parent_label='Node36995184', child_labels=['Node36995344','Node36997392'], edge_length=29.8452823467, taxon_label=None),
            'Node36995344' : NodeRelationship(parent_label='Node36962224', child_labels=['Node36995440','Node36996528'], edge_length=0.377327561902, taxon_label=None),
            'Node36995440' : NodeRelationship(parent_label='Node36995344', child_labels=['Node36995472','Node36996400'], edge_length=33.0195611855, taxon_label=None),
            'Node36995472' : NodeRelationship(parent_label='Node36995440', child_labels=['Node36995536','Morelia nauta'], edge_length=1.09375666234, taxon_label=None),
            'Node36995536' : NodeRelationship(parent_label='Node36995472', child_labels=['Node36995600','Node36995952'], edge_length=0.351619575146, taxon_label=None),
            'Node36995600' : NodeRelationship(parent_label='Node36995536', child_labels=['Apodora papuana','Node36995920'], edge_length=2.09759626003, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node36995600', child_labels=[], edge_length=0.707490222359, taxon_label='Apodora papuana'),
            'Node36995920' : NodeRelationship(parent_label='Node36995600', child_labels=['Aspidites melanocephalus','Morelia carinata'], edge_length=0.470505730507, taxon_label=None),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node36995920', child_labels=[], edge_length=0.236984491852, taxon_label='Aspidites melanocephalus'),
            'Morelia carinata' : NodeRelationship(parent_label='Node36995920', child_labels=[], edge_length=0.236984491852, taxon_label='Morelia carinata'),
            'Node36995952' : NodeRelationship(parent_label='Node36995536', child_labels=['Bothrochilus boa','Morelia amethistina'], edge_length=1.73551661307, taxon_label=None),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node36995952', child_labels=[], edge_length=1.06956986931, taxon_label='Bothrochilus boa'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node36995952', child_labels=[], edge_length=1.06956986931, taxon_label='Morelia amethistina'),
            'Morelia nauta' : NodeRelationship(parent_label='Node36995472', child_labels=[], edge_length=3.15670605753, taxon_label='Morelia nauta'),
            'Node36996400' : NodeRelationship(parent_label='Node36995440', child_labels=['Morelia tracyae','Python regius'], edge_length=3.66513378189, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node36996400', child_labels=[], edge_length=0.585328937986, taxon_label='Morelia tracyae'),
            'Python regius' : NodeRelationship(parent_label='Node36996400', child_labels=[], edge_length=0.585328937986, taxon_label='Python regius'),
            'Node36996528' : NodeRelationship(parent_label='Node36995344', child_labels=['Node36996176','Morelia boeleni'], edge_length=23.9652741942, taxon_label=None),
            'Node36996176' : NodeRelationship(parent_label='Node36996528', child_labels=['Node36996560','Node36997168'], edge_length=3.49206820101, taxon_label=None),
            'Node36996560' : NodeRelationship(parent_label='Node36996176', child_labels=['Node36996464','Python timoriensis'], edge_length=4.90381799328, taxon_label=None),
            'Node36996464' : NodeRelationship(parent_label='Node36996560', child_labels=['Python reticulatus','Node36996880'], edge_length=3.58364386958, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node36996464', child_labels=[], edge_length=1.32521964733, taxon_label='Python reticulatus'),
            'Node36996880' : NodeRelationship(parent_label='Node36996464', child_labels=['Python brongersmai','Liasis olivaceus'], edge_length=0.699088376552, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node36996880', child_labels=[], edge_length=0.626131270775, taxon_label='Python brongersmai'),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node36996880', child_labels=[], edge_length=0.626131270775, taxon_label='Liasis olivaceus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node36996560', child_labels=[], edge_length=4.9088635169, taxon_label='Python timoriensis'),
            'Node36997168' : NodeRelationship(parent_label='Node36996176', child_labels=['Node36997072','Node36997360'], edge_length=5.39574605108, taxon_label=None),
            'Node36997072' : NodeRelationship(parent_label='Node36997168', child_labels=['Morelia clastolepis','Antaresia maculosa'], edge_length=3.89664227731, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node36997072', child_labels=[], edge_length=0.520293181797, taxon_label='Morelia clastolepis'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node36997072', child_labels=[], edge_length=0.520293181797, taxon_label='Antaresia maculosa'),
            'Node36997360' : NodeRelationship(parent_label='Node36997168', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=4.2237066393, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node36997360', child_labels=[], edge_length=0.193228819804, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node36997360', child_labels=[], edge_length=0.193228819804, taxon_label='Liasis mackloti'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node36996528', child_labels=[], edge_length=13.3047497112, taxon_label='Morelia boeleni'),
            'Node36997392' : NodeRelationship(parent_label='Node36962224', child_labels=['Node36997616','Node36998064'], edge_length=3.61657514544, taxon_label=None),
            'Node36997616' : NodeRelationship(parent_label='Node36997392', child_labels=['Python molurus','Node36997840'], edge_length=21.1685275644, taxon_label=None),
            'Python molurus' : NodeRelationship(parent_label='Node36997616', child_labels=[], edge_length=12.8622487575, taxon_label='Python molurus'),
            'Node36997840' : NodeRelationship(parent_label='Node36997616', child_labels=['Node36997712','Morelia kinghorni'], edge_length=7.83034270743, taxon_label=None),
            'Node36997712' : NodeRelationship(parent_label='Node36997840', child_labels=['Antaresia stimsoni','Python sebae'], edge_length=3.83058430268, taxon_label=None),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node36997712', child_labels=[], edge_length=1.20132174738, taxon_label='Antaresia stimsoni'),
            'Python sebae' : NodeRelationship(parent_label='Node36997712', child_labels=[], edge_length=1.20132174738, taxon_label='Python sebae'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node36997840', child_labels=[], edge_length=5.03190605005, taxon_label='Morelia kinghorni'),
            'Node36998064' : NodeRelationship(parent_label='Node36997392', child_labels=['Node36997968','Node36998704'], edge_length=26.4364769675, taxon_label=None),
            'Node36997968' : NodeRelationship(parent_label='Node36998064', child_labels=['Aspidites ramsayi','Node36998320'], edge_length=4.80697376221, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node36997968', child_labels=[], edge_length=2.78732559217, taxon_label='Aspidites ramsayi'),
            'Node36998320' : NodeRelationship(parent_label='Node36997968', child_labels=['Node36998032','Antaresia perthensis'], edge_length=1.381705203, taxon_label=None),
            'Node36998032' : NodeRelationship(parent_label='Node36998320', child_labels=['Antaresia childreni','Leiopython albertisii'], edge_length=1.30019225223, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node36998032', child_labels=[], edge_length=0.105428136944, taxon_label='Antaresia childreni'),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node36998032', child_labels=[], edge_length=0.105428136944, taxon_label='Leiopython albertisii'),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node36998320', child_labels=[], edge_length=1.40562038917, taxon_label='Antaresia perthensis'),
            'Node36998704' : NodeRelationship(parent_label='Node36998064', child_labels=['Morelia viridis','Morelia bredli'], edge_length=5.17312138897, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node36998704', child_labels=[], edge_length=2.42117796541, taxon_label='Morelia viridis'),
            'Morelia bredli' : NodeRelationship(parent_label='Node36998704', child_labels=[], edge_length=2.42117796541, taxon_label='Morelia bredli'),
        },
        {
            'Node36995120' : NodeRelationship(parent_label=None, child_labels=['Morelia boeleni','Node36999056'], edge_length=None, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node36995120', child_labels=[], edge_length=26.0384706379, taxon_label='Morelia boeleni'),
            'Node36999056' : NodeRelationship(parent_label='Node36995120', child_labels=['Node36995152','Node37062192'], edge_length=5.21927124134, taxon_label=None),
            'Node36995152' : NodeRelationship(parent_label='Node36999056', child_labels=['Node36962160','Node37061744'], edge_length=0.952265882909, taxon_label=None),
            'Node36962160' : NodeRelationship(parent_label='Node36995152', child_labels=['Node36999024','Node37060720'], edge_length=1.64883491058, taxon_label=None),
            'Node36999024' : NodeRelationship(parent_label='Node36962160', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=2.95693808858, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node36999024', child_labels=[], edge_length=15.2611605145, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node36999024', child_labels=[], edge_length=15.2611605145, taxon_label='Bothrochilus boa'),
            'Node37060720' : NodeRelationship(parent_label='Node36962160', child_labels=['Morelia carinata','Node37061040'], edge_length=2.15553698471, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node37060720', child_labels=[], edge_length=16.0625616183, taxon_label='Morelia carinata'),
            'Node37061040' : NodeRelationship(parent_label='Node37060720', child_labels=['Node37060944','Morelia viridis'], edge_length=1.37038749174, taxon_label=None),
            'Node37060944' : NodeRelationship(parent_label='Node37061040', child_labels=['Node37061072','Node37061392'], edge_length=0.166805015451, taxon_label=None),
            'Node37061072' : NodeRelationship(parent_label='Node37060944', child_labels=['Morelia spilota','Morelia bredli'], edge_length=7.4747693113, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37061072', child_labels=[], edge_length=7.05059979983, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37061072', child_labels=[], edge_length=7.05059979983, taxon_label='Morelia bredli'),
            'Node37061392' : NodeRelationship(parent_label='Node37060944', child_labels=['Antaresia maculosa','Node37061456'], edge_length=1.43991095233, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37061392', child_labels=[], edge_length=13.0854581588, taxon_label='Antaresia maculosa'),
            'Node37061456' : NodeRelationship(parent_label='Node37061392', child_labels=['Antaresia perthensis','Node37061584'], edge_length=2.85871098218, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37061456', child_labels=[], edge_length=10.2267471766, taxon_label='Antaresia perthensis'),
            'Node37061584' : NodeRelationship(parent_label='Node37061456', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=5.35604513133, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37061584', child_labels=[], edge_length=4.8707020453, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37061584', child_labels=[], edge_length=4.8707020453, taxon_label='Antaresia stimsoni'),
            'Morelia viridis' : NodeRelationship(parent_label='Node37061040', child_labels=[], edge_length=14.6921741266, taxon_label='Morelia viridis'),
            'Node37061744' : NodeRelationship(parent_label='Node36995152', child_labels=['Node37061872','Apodora papuana'], edge_length=7.07533845302, taxon_label=None),
            'Node37061872' : NodeRelationship(parent_label='Node37061744', child_labels=['Liasis olivaceus','Node37062064'], edge_length=1.94094943657, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37061872', child_labels=[], edge_length=10.850645624, taxon_label='Liasis olivaceus'),
            'Node37062064' : NodeRelationship(parent_label='Node37061872', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=8.02016262021, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37062064', child_labels=[], edge_length=2.83048300382, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37062064', child_labels=[], edge_length=2.83048300382, taxon_label='Liasis mackloti'),
            'Apodora papuana' : NodeRelationship(parent_label='Node37061744', child_labels=[], edge_length=12.7915950606, taxon_label='Apodora papuana'),
            'Node37062192' : NodeRelationship(parent_label='Node36999056', child_labels=['Node37062288','Node37063088'], edge_length=0.340827566292, taxon_label=None),
            'Node37062288' : NodeRelationship(parent_label='Node37062192', child_labels=['Morelia oenpelliensis','Node37062512'], edge_length=1.42657874327, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37062288', child_labels=[], edge_length=19.051793087, taxon_label='Morelia oenpelliensis'),
            'Node37062512' : NodeRelationship(parent_label='Node37062288', child_labels=['Node37062448','Node37063184'], edge_length=1.77859521326, taxon_label=None),
            'Node37062448' : NodeRelationship(parent_label='Node37062512', child_labels=['Node37062544','Morelia tracyae'], edge_length=7.19864380037, taxon_label=None),
            'Node37062544' : NodeRelationship(parent_label='Node37062448', child_labels=['Morelia amethistina','Node37062768'], edge_length=1.9541248425, taxon_label=None),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37062544', child_labels=[], edge_length=8.12042923083, taxon_label='Morelia amethistina'),
            'Node37062768' : NodeRelationship(parent_label='Node37062544', child_labels=['Morelia clastolepis','Node37062896'], edge_length=6.05470978372, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37062768', child_labels=[], edge_length=2.0657194471, taxon_label='Morelia clastolepis'),
            'Node37062896' : NodeRelationship(parent_label='Node37062768', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=1.06606772283, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37062896', child_labels=[], edge_length=0.999651724272, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node37062896', child_labels=[], edge_length=0.999651724272, taxon_label='Morelia nauta'),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37062448', child_labels=[], edge_length=10.0745540733, taxon_label='Morelia tracyae'),
            'Node37063184' : NodeRelationship(parent_label='Node37062512', child_labels=['Python reticulatus','Python timoriensis'], edge_length=10.657346381, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37063184', child_labels=[], edge_length=6.6158514927, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37063184', child_labels=[], edge_length=6.6158514927, taxon_label='Python timoriensis'),
            'Node37063088' : NodeRelationship(parent_label='Node37062192', child_labels=['Node37063248','Node37063760'], edge_length=2.1717154248, taxon_label=None),
            'Node37063248' : NodeRelationship(parent_label='Node37063088', child_labels=['Python sebae','Node37063536'], edge_length=3.93324105164, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node37063248', child_labels=[], edge_length=14.3734153538, taxon_label='Python sebae'),
            'Node37063536' : NodeRelationship(parent_label='Node37063248', child_labels=['Python regius','Node37063664'], edge_length=2.6097491119, taxon_label=None),
            'Python regius' : NodeRelationship(parent_label='Node37063536', child_labels=[], edge_length=11.7636662419, taxon_label='Python regius'),
            'Node37063664' : NodeRelationship(parent_label='Node37063536', child_labels=['Python brongersmai','Python molurus'], edge_length=1.23184615824, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node37063664', child_labels=[], edge_length=10.5318200836, taxon_label='Python brongersmai'),
            'Python molurus' : NodeRelationship(parent_label='Node37063664', child_labels=[], edge_length=10.5318200836, taxon_label='Python molurus'),
            'Node37063760' : NodeRelationship(parent_label='Node37063088', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=8.5794579294, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37063760', child_labels=[], edge_length=9.72719847603, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37063760', child_labels=[], edge_length=9.72719847603, taxon_label='Aspidites melanocephalus'),
        },
        {
            'Node37064144' : NodeRelationship(parent_label=None, child_labels=['Node37064208','Node37656048'], edge_length=None, taxon_label=None),
            'Node37064208' : NodeRelationship(parent_label='Node37064144', child_labels=['Node36998896','Node37655728'], edge_length=0.0364728670888, taxon_label=None),
            'Node36998896' : NodeRelationship(parent_label='Node37064208', child_labels=['Node37063696','Node37654704'], edge_length=1.06929723627, taxon_label=None),
            'Node37063696' : NodeRelationship(parent_label='Node36998896', child_labels=['Node37064048','Node37654640'], edge_length=3.261105684, taxon_label=None),
            'Node37064048' : NodeRelationship(parent_label='Node37063696', child_labels=['Python sebae','Node37064496'], edge_length=2.97395412534, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node37064048', child_labels=[], edge_length=16.9224723285, taxon_label='Python sebae'),
            'Node37064496' : NodeRelationship(parent_label='Node37064048', child_labels=['Node37064432','Python regius'], edge_length=2.75174107457, taxon_label=None),
            'Node37064432' : NodeRelationship(parent_label='Node37064496', child_labels=['Python molurus','Python brongersmai'], edge_length=0.218628811139, taxon_label=None),
            'Python molurus' : NodeRelationship(parent_label='Node37064432', child_labels=[], edge_length=13.9521024428, taxon_label='Python molurus'),
            'Python brongersmai' : NodeRelationship(parent_label='Node37064432', child_labels=[], edge_length=13.9521024428, taxon_label='Python brongersmai'),
            'Python regius' : NodeRelationship(parent_label='Node37064496', child_labels=[], edge_length=14.170731254, taxon_label='Python regius'),
            'Node37654640' : NodeRelationship(parent_label='Node37063696', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=12.0688826259, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37654640', child_labels=[], edge_length=7.82754382801, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37654640', child_labels=[], edge_length=7.82754382801, taxon_label='Aspidites melanocephalus'),
            'Node37654704' : NodeRelationship(parent_label='Node36998896', child_labels=['Node37654896','Morelia viridis'], edge_length=5.48562103095, taxon_label=None),
            'Node37654896' : NodeRelationship(parent_label='Node37654704', child_labels=['Node37064656','Node37655152'], edge_length=2.9034791189, taxon_label=None),
            'Node37064656' : NodeRelationship(parent_label='Node37654896', child_labels=['Morelia carinata','Node37655120'], edge_length=2.97877570893, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node37064656', child_labels=[], edge_length=11.7896562791, taxon_label='Morelia carinata'),
            'Node37655120' : NodeRelationship(parent_label='Node37064656', child_labels=['Morelia spilota','Morelia bredli'], edge_length=7.73520362135, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37655120', child_labels=[], edge_length=4.05445265776, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37655120', child_labels=[], edge_length=4.05445265776, taxon_label='Morelia bredli'),
            'Node37655152' : NodeRelationship(parent_label='Node37654896', child_labels=['Antaresia maculosa','Node37655408'], edge_length=1.25780607691, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37655152', child_labels=[], edge_length=13.5106259111, taxon_label='Antaresia maculosa'),
            'Node37655408' : NodeRelationship(parent_label='Node37655152', child_labels=['Antaresia perthensis','Node37655536'], edge_length=2.27459934102, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37655408', child_labels=[], edge_length=11.2360265701, taxon_label='Antaresia perthensis'),
            'Node37655536' : NodeRelationship(parent_label='Node37655408', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=7.74754730265, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37655536', child_labels=[], edge_length=3.48847926746, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37655536', child_labels=[], edge_length=3.48847926746, taxon_label='Antaresia stimsoni'),
            'Morelia viridis' : NodeRelationship(parent_label='Node37654704', child_labels=[], edge_length=17.6719111069, taxon_label='Morelia viridis'),
            'Node37655728' : NodeRelationship(parent_label='Node37064208', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=10.8687376261, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node37655728', child_labels=[], edge_length=13.3580917481, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node37655728', child_labels=[], edge_length=13.3580917481, taxon_label='Bothrochilus boa'),
            'Node37656048' : NodeRelationship(parent_label='Node37064144', child_labels=['Node37655984','Node37656560'], edge_length=1.25568695854, taxon_label=None),
            'Node37655984' : NodeRelationship(parent_label='Node37656048', child_labels=['Node37655632','Apodora papuana'], edge_length=14.878020141, taxon_label=None),
            'Node37655632' : NodeRelationship(parent_label='Node37655984', child_labels=['Liasis olivaceus','Node37656272'], edge_length=1.91471562327, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37655632', child_labels=[], edge_length=6.21487951846, taxon_label='Liasis olivaceus'),
            'Node37656272' : NodeRelationship(parent_label='Node37655632', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=3.60951843288, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37656272', child_labels=[], edge_length=2.60536108558, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37656272', child_labels=[], edge_length=2.60536108558, taxon_label='Liasis mackloti'),
            'Apodora papuana' : NodeRelationship(parent_label='Node37655984', child_labels=[], edge_length=8.12959514173, taxon_label='Apodora papuana'),
            'Node37656560' : NodeRelationship(parent_label='Node37656048', child_labels=['Morelia boeleni','Node37656656'], edge_length=2.53457554003, taxon_label=None),
            'Morelia boeleni' : NodeRelationship(parent_label='Node37656560', child_labels=[], edge_length=20.4730397427, taxon_label='Morelia boeleni'),
            'Node37656656' : NodeRelationship(parent_label='Node37656560', child_labels=['Node37656528','Morelia oenpelliensis'], edge_length=2.57209904849, taxon_label=None),
            'Node37656528' : NodeRelationship(parent_label='Node37656656', child_labels=['Node37656688','Node37657392'], edge_length=0.296424924481, taxon_label=None),
            'Node37656688' : NodeRelationship(parent_label='Node37656528', child_labels=['Node37656304','Morelia tracyae'], edge_length=9.88501788667, taxon_label=None),
            'Node37656304' : NodeRelationship(parent_label='Node37656688', child_labels=['Morelia amethistina','Node37656976'], edge_length=3.28625146667, taxon_label=None),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37656304', child_labels=[], edge_length=4.43324641637, taxon_label='Morelia amethistina'),
            'Node37656976' : NodeRelationship(parent_label='Node37656304', child_labels=['Node37656912','Morelia clastolepis'], edge_length=2.99969759558, taxon_label=None),
            'Node37656912' : NodeRelationship(parent_label='Node37656976', child_labels=['Morelia nauta','Morelia kinghorni'], edge_length=0.622029416262, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node37656912', child_labels=[], edge_length=0.811519404534, taxon_label='Morelia nauta'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37656912', child_labels=[], edge_length=0.811519404534, taxon_label='Morelia kinghorni'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37656976', child_labels=[], edge_length=1.4335488208, taxon_label='Morelia clastolepis'),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37656688', child_labels=[], edge_length=7.71949788304, taxon_label='Morelia tracyae'),
            'Node37657392' : NodeRelationship(parent_label='Node37656528', child_labels=['Python reticulatus','Python timoriensis'], edge_length=11.8234917746, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37657392', child_labels=[], edge_length=5.78102399509, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37657392', child_labels=[], edge_length=5.78102399509, taxon_label='Python timoriensis'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37656656', child_labels=[], edge_length=17.9009406942, taxon_label='Morelia oenpelliensis'),
        },
        {
            'Node37657712' : NodeRelationship(parent_label=None, child_labels=['Node37657776','Node37658384'], edge_length=None, taxon_label=None),
            'Node37657776' : NodeRelationship(parent_label='Node37657712', child_labels=['Node36962032','Node37658224'], edge_length=3.09272413642, taxon_label=None),
            'Node36962032' : NodeRelationship(parent_label='Node37657776', child_labels=['Python sebae','Node37658000'], edge_length=3.43039474792, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node36962032', child_labels=[], edge_length=16.5921325784, taxon_label='Python sebae'),
            'Node37658000' : NodeRelationship(parent_label='Node36962032', child_labels=['Python regius','Node37658128'], edge_length=2.21219736911, taxon_label=None),
            'Python regius' : NodeRelationship(parent_label='Node37658000', child_labels=[], edge_length=14.3799352093, taxon_label='Python regius'),
            'Node37658128' : NodeRelationship(parent_label='Node37658000', child_labels=['Python molurus','Python brongersmai'], edge_length=0.997806328959, taxon_label=None),
            'Python molurus' : NodeRelationship(parent_label='Node37658128', child_labels=[], edge_length=13.3821288803, taxon_label='Python molurus'),
            'Python brongersmai' : NodeRelationship(parent_label='Node37658128', child_labels=[], edge_length=13.3821288803, taxon_label='Python brongersmai'),
            'Node37658224' : NodeRelationship(parent_label='Node37657776', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=8.31909563528, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37658224', child_labels=[], edge_length=11.703431691, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37658224', child_labels=[], edge_length=11.703431691, taxon_label='Aspidites melanocephalus'),
            'Node37658384' : NodeRelationship(parent_label='Node37657712', child_labels=['Node37657648','Node37712912'], edge_length=0.465997253751, taxon_label=None),
            'Node37657648' : NodeRelationship(parent_label='Node37658384', child_labels=['Node37658288','Morelia oenpelliensis'], edge_length=2.71760756257, taxon_label=None),
            'Node37658288' : NodeRelationship(parent_label='Node37657648', child_labels=['Node37658160','Node37658576'], edge_length=0.121667790832, taxon_label=None),
            'Node37658160' : NodeRelationship(parent_label='Node37658288', child_labels=['Python reticulatus','Python timoriensis'], edge_length=13.3467104377, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37658160', child_labels=[], edge_length=6.46326841791, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37658160', child_labels=[], edge_length=6.46326841791, taxon_label='Python timoriensis'),
            'Node37658576' : NodeRelationship(parent_label='Node37658288', child_labels=['Node37712016','Morelia boeleni'], edge_length=10.3912007878, taxon_label=None),
            'Node37712016' : NodeRelationship(parent_label='Node37658576', child_labels=['Node37711984','Morelia tracyae'], edge_length=3.03031595643, taxon_label=None),
            'Node37711984' : NodeRelationship(parent_label='Node37712016', child_labels=['Morelia amethistina','Node37712368'], edge_length=1.71754130617, taxon_label=None),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37711984', child_labels=[], edge_length=4.67092080519, taxon_label='Morelia amethistina'),
            'Node37712368' : NodeRelationship(parent_label='Node37711984', child_labels=['Morelia nauta','Node37712496'], edge_length=2.58211994653, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node37712368', child_labels=[], edge_length=2.08880085866, taxon_label='Morelia nauta'),
            'Node37712496' : NodeRelationship(parent_label='Node37712368', child_labels=['Morelia kinghorni','Morelia clastolepis'], edge_length=0.607570395039, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37712496', child_labels=[], edge_length=1.48123046362, taxon_label='Morelia kinghorni'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37712496', child_labels=[], edge_length=1.48123046362, taxon_label='Morelia clastolepis'),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37712016', child_labels=[], edge_length=6.38846211136, taxon_label='Morelia tracyae'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node37658576', child_labels=[], edge_length=9.41877806779, taxon_label='Morelia boeleni'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37657648', child_labels=[], edge_length=19.9316466464, taxon_label='Morelia oenpelliensis'),
            'Node37712912' : NodeRelationship(parent_label='Node37658384', child_labels=['Node37712848','Node37714192'], edge_length=0.725019015485, taxon_label=None),
            'Node37712848' : NodeRelationship(parent_label='Node37712912', child_labels=['Node37712816','Node37714000'], edge_length=1.89223670238, taxon_label=None),
            'Node37712816' : NodeRelationship(parent_label='Node37712848', child_labels=['Node37712400','Morelia carinata'], edge_length=1.84844364601, taxon_label=None),
            'Node37712400' : NodeRelationship(parent_label='Node37712816', child_labels=['Morelia viridis','Node37713200'], edge_length=2.23472642535, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node37712400', child_labels=[], edge_length=15.9488284198, taxon_label='Morelia viridis'),
            'Node37713200' : NodeRelationship(parent_label='Node37712400', child_labels=['Node37713072','Node37713648'], edge_length=0.394642188595, taxon_label=None),
            'Node37713072' : NodeRelationship(parent_label='Node37713200', child_labels=['Antaresia maculosa','Node37713360'], edge_length=1.55388791109, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37713072', child_labels=[], edge_length=14.0002983201, taxon_label='Antaresia maculosa'),
            'Node37713360' : NodeRelationship(parent_label='Node37713072', child_labels=['Antaresia perthensis','Node37713488'], edge_length=2.5525380966, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37713360', child_labels=[], edge_length=11.4477602235, taxon_label='Antaresia perthensis'),
            'Node37713488' : NodeRelationship(parent_label='Node37713360', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=5.71931728003, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37713488', child_labels=[], edge_length=5.72844294344, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37713488', child_labels=[], edge_length=5.72844294344, taxon_label='Antaresia stimsoni'),
            'Node37713648' : NodeRelationship(parent_label='Node37713200', child_labels=['Morelia spilota','Morelia bredli'], edge_length=11.0588698605, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37713648', child_labels=[], edge_length=4.49531637063, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37713648', child_labels=[], edge_length=4.49531637063, taxon_label='Morelia bredli'),
            'Morelia carinata' : NodeRelationship(parent_label='Node37712816', child_labels=[], edge_length=18.1835548451, taxon_label='Morelia carinata'),
            'Node37714000' : NodeRelationship(parent_label='Node37712848', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=7.28635782162, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node37714000', child_labels=[], edge_length=12.7456406695, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node37714000', child_labels=[], edge_length=12.7456406695, taxon_label='Bothrochilus boa'),
            'Node37714192' : NodeRelationship(parent_label='Node37712912', child_labels=['Node37714128','Apodora papuana'], edge_length=6.6830638124, taxon_label=None),
            'Node37714128' : NodeRelationship(parent_label='Node37714192', child_labels=['Liasis olivaceus','Node37714352'], edge_length=1.71076794287, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37714128', child_labels=[], edge_length=13.5304034382, taxon_label='Liasis olivaceus'),
            'Node37714352' : NodeRelationship(parent_label='Node37714128', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=11.7858576376, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37714352', child_labels=[], edge_length=1.74454580066, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37714352', child_labels=[], edge_length=1.74454580066, taxon_label='Liasis mackloti'),
            'Apodora papuana' : NodeRelationship(parent_label='Node37714192', child_labels=[], edge_length=15.2411713811, taxon_label='Apodora papuana'),
        },
        {
            'Node37714704' : NodeRelationship(parent_label=None, child_labels=['Node37714768','Node37715024'], edge_length=None, taxon_label=None),
            'Node37714768' : NodeRelationship(parent_label='Node37714704', child_labels=['Python reticulatus','Python timoriensis'], edge_length=17.6343232815, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37714768', child_labels=[], edge_length=10.8300265774, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37714768', child_labels=[], edge_length=10.8300265774, taxon_label='Python timoriensis'),
            'Node37715024' : NodeRelationship(parent_label='Node37714704', child_labels=['Node37714512','Node37773872'], edge_length=3.72471742653, taxon_label=None),
            'Node37714512' : NodeRelationship(parent_label='Node37715024', child_labels=['Node37714960','Node37773680'], edge_length=4.52906062503, taxon_label=None),
            'Node37714960' : NodeRelationship(parent_label='Node37714512', child_labels=['Node37714992','Node37715600'], edge_length=1.31772152158, taxon_label=None),
            'Node37714992' : NodeRelationship(parent_label='Node37714960', child_labels=['Morelia viridis','Node37715312'], edge_length=4.88508122224, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node37714992', child_labels=[], edge_length=14.0077690635, taxon_label='Morelia viridis'),
            'Node37715312' : NodeRelationship(parent_label='Node37714992', child_labels=['Node37715184','Node37715536'], edge_length=0.404170030911, taxon_label=None),
            'Node37715184' : NodeRelationship(parent_label='Node37715312', child_labels=['Morelia carinata','Node37715504'], edge_length=1.82500345118, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node37715184', child_labels=[], edge_length=11.7785955814, taxon_label='Morelia carinata'),
            'Node37715504' : NodeRelationship(parent_label='Node37715184', child_labels=['Morelia spilota','Morelia bredli'], edge_length=6.80279787983, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37715504', child_labels=[], edge_length=4.97579770159, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37715504', child_labels=[], edge_length=4.97579770159, taxon_label='Morelia bredli'),
            'Node37715536' : NodeRelationship(parent_label='Node37715312', child_labels=['Antaresia maculosa','Node37715792'], edge_length=1.27612254021, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37715536', child_labels=[], edge_length=12.3274764924, taxon_label='Antaresia maculosa'),
            'Node37715792' : NodeRelationship(parent_label='Node37715536', child_labels=['Antaresia perthensis','Node37715760'], edge_length=4.32415504222, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37715792', child_labels=[], edge_length=8.00332145017, taxon_label='Antaresia perthensis'),
            'Node37715760' : NodeRelationship(parent_label='Node37715792', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=3.98809707211, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37715760', child_labels=[], edge_length=4.01522437806, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37715760', child_labels=[], edge_length=4.01522437806, taxon_label='Antaresia stimsoni'),
            'Node37715600' : NodeRelationship(parent_label='Node37714960', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=7.74780761159, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node37715600', child_labels=[], edge_length=11.1450426742, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node37715600', child_labels=[], edge_length=11.1450426742, taxon_label='Bothrochilus boa'),
            'Node37773680' : NodeRelationship(parent_label='Node37714512', child_labels=['Node37773616','Apodora papuana'], edge_length=5.44916502448, taxon_label=None),
            'Node37773616' : NodeRelationship(parent_label='Node37773680', child_labels=['Liasis olivaceus','Node37773840'], edge_length=3.78712208743, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37773616', child_labels=[], edge_length=10.9742846954, taxon_label='Liasis olivaceus'),
            'Node37773840' : NodeRelationship(parent_label='Node37773616', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=7.80758505868, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37773840', child_labels=[], edge_length=3.16669963676, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37773840', child_labels=[], edge_length=3.16669963676, taxon_label='Liasis mackloti'),
            'Apodora papuana' : NodeRelationship(parent_label='Node37773680', child_labels=[], edge_length=14.7614067829, taxon_label='Apodora papuana'),
            'Node37773872' : NodeRelationship(parent_label='Node37715024', child_labels=['Node37774128','Node37775024'], edge_length=0.57539806967, taxon_label=None),
            'Node37774128' : NodeRelationship(parent_label='Node37773872', child_labels=['Node37774032','Morelia oenpelliensis'], edge_length=1.8069504389, taxon_label=None),
            'Node37774032' : NodeRelationship(parent_label='Node37774128', child_labels=['Node37773904','Morelia boeleni'], edge_length=5.48583318771, taxon_label=None),
            'Node37773904' : NodeRelationship(parent_label='Node37774032', child_labels=['Node37774160','Morelia tracyae'], edge_length=6.85886795982, taxon_label=None),
            'Node37774160' : NodeRelationship(parent_label='Node37773904', child_labels=['Morelia amethistina','Node37774480'], edge_length=4.50938662462, taxon_label=None),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37774160', child_labels=[], edge_length=5.50319615164, taxon_label='Morelia amethistina'),
            'Node37774480' : NodeRelationship(parent_label='Node37774160', child_labels=['Morelia clastolepis','Node37774608'], edge_length=3.93061957053, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37774480', child_labels=[], edge_length=1.57257658111, taxon_label='Morelia clastolepis'),
            'Node37774608' : NodeRelationship(parent_label='Node37774480', child_labels=['Morelia nauta','Morelia kinghorni'], edge_length=0.136848033805, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node37774608', child_labels=[], edge_length=1.4357285473, taxon_label='Morelia nauta'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37774608', child_labels=[], edge_length=1.4357285473, taxon_label='Morelia kinghorni'),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37773904', child_labels=[], edge_length=10.0125827763, taxon_label='Morelia tracyae'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node37774032', child_labels=[], edge_length=16.8714507361, taxon_label='Morelia boeleni'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37774128', child_labels=[], edge_length=22.3572839238, taxon_label='Morelia oenpelliensis'),
            'Node37775024' : NodeRelationship(parent_label='Node37773872', child_labels=['Node37774928','Node37775408'], edge_length=3.10812197386, taxon_label=None),
            'Node37774928' : NodeRelationship(parent_label='Node37775024', child_labels=['Python sebae','Node37775184'], edge_length=4.52064732483, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node37774928', child_labels=[], edge_length=16.535465064, taxon_label='Python sebae'),
            'Node37775184' : NodeRelationship(parent_label='Node37774928', child_labels=['Python regius','Node37775312'], edge_length=1.573045967, taxon_label=None),
            'Python regius' : NodeRelationship(parent_label='Node37775184', child_labels=[], edge_length=14.962419097, taxon_label='Python regius'),
            'Node37775312' : NodeRelationship(parent_label='Node37775184', child_labels=['Python brongersmai','Python molurus'], edge_length=1.82698701194, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node37775312', child_labels=[], edge_length=13.1354320851, taxon_label='Python brongersmai'),
            'Python molurus' : NodeRelationship(parent_label='Node37775312', child_labels=[], edge_length=13.1354320851, taxon_label='Python molurus'),
            'Node37775408' : NodeRelationship(parent_label='Node37775024', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=6.39210928917, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37775408', child_labels=[], edge_length=14.6640030997, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37775408', child_labels=[], edge_length=14.6640030997, taxon_label='Aspidites melanocephalus'),
        },
        {
            'Node37775792' : NodeRelationship(parent_label=None, child_labels=['Node37775856','Node37831600'], edge_length=None, taxon_label=None),
            'Node37775856' : NodeRelationship(parent_label='Node37775792', child_labels=['Node37064112','Node37830864'], edge_length=1.01117536037, taxon_label=None),
            'Node37064112' : NodeRelationship(parent_label='Node37775856', child_labels=['Node37775568','Node37777264'], edge_length=4.09971291135, taxon_label=None),
            'Node37775568' : NodeRelationship(parent_label='Node37064112', child_labels=['Node37775472','Node37777072'], edge_length=0.754225130388, taxon_label=None),
            'Node37775472' : NodeRelationship(parent_label='Node37775568', child_labels=['Node37775952','Morelia carinata'], edge_length=1.65141371467, taxon_label=None),
            'Node37775952' : NodeRelationship(parent_label='Node37775472', child_labels=['Node37776016','Node37776144'], edge_length=0.835745557531, taxon_label=None),
            'Node37776016' : NodeRelationship(parent_label='Node37775952', child_labels=['Antaresia maculosa','Node37776304'], edge_length=1.83427060068, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37776016', child_labels=[], edge_length=12.6196623606, taxon_label='Antaresia maculosa'),
            'Node37776304' : NodeRelationship(parent_label='Node37776016', child_labels=['Antaresia perthensis','Node37776432'], edge_length=2.51279801319, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37776304', child_labels=[], edge_length=10.1068643474, taxon_label='Antaresia perthensis'),
            'Node37776432' : NodeRelationship(parent_label='Node37776304', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=5.71367528366, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37776432', child_labels=[], edge_length=4.39318906371, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37776432', child_labels=[], edge_length=4.39318906371, taxon_label='Antaresia stimsoni'),
            'Node37776144' : NodeRelationship(parent_label='Node37775952', child_labels=['Morelia viridis','Node37776784'], edge_length=1.03775314121, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node37776144', child_labels=[], edge_length=13.41617982, taxon_label='Morelia viridis'),
            'Node37776784' : NodeRelationship(parent_label='Node37776144', child_labels=['Morelia spilota','Morelia bredli'], edge_length=9.37635268364, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37776784', child_labels=[], edge_length=4.03982713639, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37776784', child_labels=[], edge_length=4.03982713639, taxon_label='Morelia bredli'),
            'Morelia carinata' : NodeRelationship(parent_label='Node37775472', child_labels=[], edge_length=15.2896785188, taxon_label='Morelia carinata'),
            'Node37777072' : NodeRelationship(parent_label='Node37775568', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=4.17813417278, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node37777072', child_labels=[], edge_length=12.7629580607, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node37777072', child_labels=[], edge_length=12.7629580607, taxon_label='Bothrochilus boa'),
            'Node37777264' : NodeRelationship(parent_label='Node37064112', child_labels=['Node37777232','Apodora papuana'], edge_length=2.9800994962, taxon_label=None),
            'Node37777232' : NodeRelationship(parent_label='Node37777264', child_labels=['Liasis olivaceus','Node37777360'], edge_length=4.06497756637, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37777232', child_labels=[], edge_length=10.6502403013, taxon_label='Liasis olivaceus'),
            'Node37777360' : NodeRelationship(parent_label='Node37777232', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=6.97858461332, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37777360', child_labels=[], edge_length=3.67165568794, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37777360', child_labels=[], edge_length=3.67165568794, taxon_label='Liasis mackloti'),
            'Apodora papuana' : NodeRelationship(parent_label='Node37777264', child_labels=[], edge_length=14.7152178676, taxon_label='Apodora papuana'),
            'Node37830864' : NodeRelationship(parent_label='Node37775856', child_labels=['Node37830704','Morelia oenpelliensis'], edge_length=3.12950470552, taxon_label=None),
            'Node37830704' : NodeRelationship(parent_label='Node37830864', child_labels=['Node37777040','Morelia boeleni'], edge_length=6.58446507703, taxon_label=None),
            'Node37777040' : NodeRelationship(parent_label='Node37830704', child_labels=['Node37830928','Morelia tracyae'], edge_length=3.1657791763, taxon_label=None),
            'Node37830928' : NodeRelationship(parent_label='Node37777040', child_labels=['Morelia amethistina','Node37831216'], edge_length=4.27056899374, taxon_label=None),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37830928', child_labels=[], edge_length=4.6447123226, taxon_label='Morelia amethistina'),
            'Node37831216' : NodeRelationship(parent_label='Node37830928', child_labels=['Morelia clastolepis','Node37831344'], edge_length=2.47296326498, taxon_label=None),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37831216', child_labels=[], edge_length=2.17174905761, taxon_label='Morelia clastolepis'),
            'Node37831344' : NodeRelationship(parent_label='Node37831216', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=1.26765325952, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37831344', child_labels=[], edge_length=0.904095798097, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node37831344', child_labels=[], edge_length=0.904095798097, taxon_label='Morelia nauta'),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37777040', child_labels=[], edge_length=8.91528131634, taxon_label='Morelia tracyae'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node37830704', child_labels=[], edge_length=12.0810604926, taxon_label='Morelia boeleni'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37830864', child_labels=[], edge_length=18.6655255697, taxon_label='Morelia oenpelliensis'),
            'Node37831600' : NodeRelationship(parent_label='Node37775792', child_labels=['Node37831504','Node37832016'], edge_length=0.075772291124, taxon_label=None),
            'Node37831504' : NodeRelationship(parent_label='Node37831600', child_labels=['Python reticulatus','Python timoriensis'], edge_length=11.1718400909, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37831504', child_labels=[], edge_length=11.5585932535, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37831504', child_labels=[], edge_length=11.5585932535, taxon_label='Python timoriensis'),
            'Node37832016' : NodeRelationship(parent_label='Node37831600', child_labels=['Node37831888','Node37832400'], edge_length=1.44865397646, taxon_label=None),
            'Node37831888' : NodeRelationship(parent_label='Node37832016', child_labels=['Python sebae','Node37832176'], edge_length=5.31414730047, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node37831888', child_labels=[], edge_length=15.9676320675, taxon_label='Python sebae'),
            'Node37832176' : NodeRelationship(parent_label='Node37831888', child_labels=['Node37832048','Python regius'], edge_length=1.06278922707, taxon_label=None),
            'Node37832048' : NodeRelationship(parent_label='Node37832176', child_labels=['Python brongersmai','Python molurus'], edge_length=1.04487219988, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node37832048', child_labels=[], edge_length=13.8599706405, taxon_label='Python brongersmai'),
            'Python molurus' : NodeRelationship(parent_label='Node37832048', child_labels=[], edge_length=13.8599706405, taxon_label='Python molurus'),
            'Python regius' : NodeRelationship(parent_label='Node37832176', child_labels=[], edge_length=14.9048428404, taxon_label='Python regius'),
            'Node37832400' : NodeRelationship(parent_label='Node37832016', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=5.70308191065, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37832400', child_labels=[], edge_length=15.5786974573, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37832400', child_labels=[], edge_length=15.5786974573, taxon_label='Aspidites melanocephalus'),
        },
        {
            'Node37832784' : NodeRelationship(parent_label=None, child_labels=['Node37832848','Node37833200'], edge_length=None, taxon_label=None),
            'Node37832848' : NodeRelationship(parent_label='Node37832784', child_labels=['Python sebae','Node37833008'], edge_length=2.12261395029, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node37832848', child_labels=[], edge_length=29.0711850036, taxon_label='Python sebae'),
            'Node37833008' : NodeRelationship(parent_label='Node37832848', child_labels=['Python molurus','Node37833136'], edge_length=6.65870230472, taxon_label=None),
            'Python molurus' : NodeRelationship(parent_label='Node37833008', child_labels=[], edge_length=22.4124826989, taxon_label='Python molurus'),
            'Node37833136' : NodeRelationship(parent_label='Node37833008', child_labels=['Python regius','Python brongersmai'], edge_length=0.485999991488, taxon_label=None),
            'Python regius' : NodeRelationship(parent_label='Node37833136', child_labels=[], edge_length=21.9264827074, taxon_label='Python regius'),
            'Python brongersmai' : NodeRelationship(parent_label='Node37833136', child_labels=[], edge_length=21.9264827074, taxon_label='Python brongersmai'),
            'Node37833200' : NodeRelationship(parent_label='Node37832784', child_labels=['Node37833360','Node37103024'], edge_length=8.39682440052, taxon_label=None),
            'Node37833360' : NodeRelationship(parent_label='Node37833200', child_labels=['Node37833328','Node37101936'], edge_length=1.24633634472, taxon_label=None),
            'Node37833328' : NodeRelationship(parent_label='Node37833360', child_labels=['Node37833168','Node37833744'], edge_length=1.10888028588, taxon_label=None),
            'Node37833168' : NodeRelationship(parent_label='Node37833328', child_labels=['Python reticulatus','Python timoriensis'], edge_length=5.89251749724, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37833168', child_labels=[], edge_length=14.5492404256, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37833168', child_labels=[], edge_length=14.5492404256, taxon_label='Python timoriensis'),
            'Node37833744' : NodeRelationship(parent_label='Node37833328', child_labels=['Node37833456','Node37101712'], edge_length=1.52019767757, taxon_label=None),
            'Node37833456' : NodeRelationship(parent_label='Node37833744', child_labels=['Node37833680','Node37834736'], edge_length=0.826870765744, taxon_label=None),
            'Node37833680' : NodeRelationship(parent_label='Node37833456', child_labels=['Node37833520','Node37834256'], edge_length=2.93802844014, taxon_label=None),
            'Node37833520' : NodeRelationship(parent_label='Node37833680', child_labels=['Morelia viridis','Node37834032'], edge_length=0.322649062335, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node37833520', child_labels=[], edge_length=14.834011977, taxon_label='Morelia viridis'),
            'Node37834032' : NodeRelationship(parent_label='Node37833520', child_labels=['Node37833904','Morelia carinata'], edge_length=1.22843868872, taxon_label=None),
            'Node37833904' : NodeRelationship(parent_label='Node37834032', child_labels=['Morelia spilota','Morelia bredli'], edge_length=7.33173338162, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37833904', child_labels=[], edge_length=6.2738399067, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37833904', child_labels=[], edge_length=6.2738399067, taxon_label='Morelia bredli'),
            'Morelia carinata' : NodeRelationship(parent_label='Node37834032', child_labels=[], edge_length=13.6055732883, taxon_label='Morelia carinata'),
            'Node37834256' : NodeRelationship(parent_label='Node37833680', child_labels=['Node37834288','Antaresia maculosa'], edge_length=0.821527298272, taxon_label=None),
            'Node37834288' : NodeRelationship(parent_label='Node37834256', child_labels=['Antaresia perthensis','Node37834512'], edge_length=4.22599220252, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37834288', child_labels=[], edge_length=10.1091415386, taxon_label='Antaresia perthensis'),
            'Node37834512' : NodeRelationship(parent_label='Node37834288', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=6.62594074811, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37834512', child_labels=[], edge_length=3.48320079048, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37834512', child_labels=[], edge_length=3.48320079048, taxon_label='Antaresia stimsoni'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37834256', child_labels=[], edge_length=14.3351337411, taxon_label='Antaresia maculosa'),
            'Node37834736' : NodeRelationship(parent_label='Node37833456', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=8.81046360463, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node37834736', child_labels=[], edge_length=9.28422587488, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node37834736', child_labels=[], edge_length=9.28422587488, taxon_label='Bothrochilus boa'),
            'Node37101712' : NodeRelationship(parent_label='Node37833744', child_labels=['Node37101808','Apodora papuana'], edge_length=5.2034050718, taxon_label=None),
            'Node37101808' : NodeRelationship(parent_label='Node37101712', child_labels=['Liasis olivaceus','Node37101968'], edge_length=2.70301025103, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37101808', child_labels=[], edge_length=11.0151449224, taxon_label='Liasis olivaceus'),
            'Node37101968' : NodeRelationship(parent_label='Node37101808', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=6.92713376313, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37101968', child_labels=[], edge_length=4.0880111593, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37101968', child_labels=[], edge_length=4.0880111593, taxon_label='Liasis mackloti'),
            'Apodora papuana' : NodeRelationship(parent_label='Node37101712', child_labels=[], edge_length=13.7181551735, taxon_label='Apodora papuana'),
            'Node37101936' : NodeRelationship(parent_label='Node37833360', child_labels=['Morelia oenpelliensis','Node37102352'], edge_length=3.94969364946, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37101936', child_labels=[], edge_length=17.6009445592, taxon_label='Morelia oenpelliensis'),
            'Node37102352' : NodeRelationship(parent_label='Node37101936', child_labels=['Node37102288','Morelia boeleni'], edge_length=5.4268859952, taxon_label=None),
            'Node37102288' : NodeRelationship(parent_label='Node37102352', child_labels=['Morelia tracyae','Node37102544'], edge_length=6.37126259081, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37102288', child_labels=[], edge_length=5.80279597323, taxon_label='Morelia tracyae'),
            'Node37102544' : NodeRelationship(parent_label='Node37102288', child_labels=['Morelia amethistina','Node37102672'], edge_length=1.42172379881, taxon_label=None),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37102544', child_labels=[], edge_length=4.38107217442, taxon_label='Morelia amethistina'),
            'Node37102672' : NodeRelationship(parent_label='Node37102544', child_labels=['Node37102608','Morelia clastolepis'], edge_length=2.50070781902, taxon_label=None),
            'Node37102608' : NodeRelationship(parent_label='Node37102672', child_labels=['Morelia nauta','Morelia kinghorni'], edge_length=0.0774488112983, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node37102608', child_labels=[], edge_length=1.8029155441, taxon_label='Morelia nauta'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37102608', child_labels=[], edge_length=1.8029155441, taxon_label='Morelia kinghorni'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37102672', child_labels=[], edge_length=1.8803643554, taxon_label='Morelia clastolepis'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node37102352', child_labels=[], edge_length=12.174058564, taxon_label='Morelia boeleni'),
            'Node37103024' : NodeRelationship(parent_label='Node37833200', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=12.9499290649, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37103024', child_labels=[], edge_length=9.84704548855, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37103024', child_labels=[], edge_length=9.84704548855, taxon_label='Aspidites melanocephalus'),
        },
        {
            'Node37103344' : NodeRelationship(parent_label=None, child_labels=['Node37103408','Python sebae'], edge_length=None, taxon_label=None),
            'Node37103408' : NodeRelationship(parent_label='Node37103344', child_labels=['Node37714672','Node37163696'], edge_length=8.43656365605, taxon_label=None),
            'Node37714672' : NodeRelationship(parent_label='Node37103408', child_labels=['Node37102896','Node37163504'], edge_length=6.41982933271, taxon_label=None),
            'Node37102896' : NodeRelationship(parent_label='Node37714672', child_labels=['Node37103120','Node37103280'], edge_length=4.2262164672, taxon_label=None),
            'Node37103120' : NodeRelationship(parent_label='Node37102896', child_labels=['Node37103504','Node37104432'], edge_length=0.983173088098, taxon_label=None),
            'Node37103504' : NodeRelationship(parent_label='Node37103120', child_labels=['Node37103568','Morelia boeleni'], edge_length=9.50779575096, taxon_label=None),
            'Node37103568' : NodeRelationship(parent_label='Node37103504', child_labels=['Morelia tracyae','Node37103888'], edge_length=4.13809660753, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37103568', child_labels=[], edge_length=8.31432086046, taxon_label='Morelia tracyae'),
            'Node37103888' : NodeRelationship(parent_label='Node37103568', child_labels=['Node37103824','Morelia amethistina'], edge_length=3.80591263593, taxon_label=None),
            'Node37103824' : NodeRelationship(parent_label='Node37103888', child_labels=['Node37103920','Morelia kinghorni'], edge_length=2.88132393295, taxon_label=None),
            'Node37103920' : NodeRelationship(parent_label='Node37103824', child_labels=['Morelia nauta','Morelia clastolepis'], edge_length=0.71440797492, taxon_label=None),
            'Morelia nauta' : NodeRelationship(parent_label='Node37103920', child_labels=[], edge_length=0.912676316658, taxon_label='Morelia nauta'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37103920', child_labels=[], edge_length=0.912676316658, taxon_label='Morelia clastolepis'),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37103824', child_labels=[], edge_length=1.62708429158, taxon_label='Morelia kinghorni'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37103888', child_labels=[], edge_length=4.50840822453, taxon_label='Morelia amethistina'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node37103504', child_labels=[], edge_length=12.452417468, taxon_label='Morelia boeleni'),
            'Node37104432' : NodeRelationship(parent_label='Node37103120', child_labels=['Node37104336','Node37104784'], edge_length=3.20592158422, taxon_label=None),
            'Node37104336' : NodeRelationship(parent_label='Node37104432', child_labels=['Apodora papuana','Node37104592'], edge_length=1.50938344234, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node37104336', child_labels=[], edge_length=17.2449081924, taxon_label='Apodora papuana'),
            'Node37104592' : NodeRelationship(parent_label='Node37104336', child_labels=['Liasis olivaceus','Node37104720'], edge_length=4.9943499411, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37104592', child_labels=[], edge_length=12.2505582513, taxon_label='Liasis olivaceus'),
            'Node37104720' : NodeRelationship(parent_label='Node37104592', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=5.95835956047, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37104720', child_labels=[], edge_length=6.29219869082, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37104720', child_labels=[], edge_length=6.29219869082, taxon_label='Liasis mackloti'),
            'Node37104784' : NodeRelationship(parent_label='Node37104432', child_labels=['Node37104752','Node37104976'], edge_length=0.497596402447, taxon_label=None),
            'Node37104752' : NodeRelationship(parent_label='Node37104784', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=5.19337056498, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node37104752', child_labels=[], edge_length=13.0633246673, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node37104752', child_labels=[], edge_length=13.0633246673, taxon_label='Bothrochilus boa'),
            'Node37104976' : NodeRelationship(parent_label='Node37104784', child_labels=['Node37104880','Node37105584'], edge_length=2.40513075673, taxon_label=None),
            'Node37104880' : NodeRelationship(parent_label='Node37104976', child_labels=['Node37105200','Antaresia maculosa'], edge_length=2.06371429021, taxon_label=None),
            'Node37105200' : NodeRelationship(parent_label='Node37104880', child_labels=['Antaresia perthensis','Node37105392'], edge_length=3.43791705079, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37105200', child_labels=[], edge_length=10.3499331345, taxon_label='Antaresia perthensis'),
            'Node37105392' : NodeRelationship(parent_label='Node37105200', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=6.96420024077, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37105392', child_labels=[], edge_length=3.38573289378, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37105392', child_labels=[], edge_length=3.38573289378, taxon_label='Antaresia stimsoni'),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37104880', child_labels=[], edge_length=13.7878501853, taxon_label='Antaresia maculosa'),
            'Node37105584' : NodeRelationship(parent_label='Node37104976', child_labels=['Morelia viridis','Node37105232'], edge_length=0.641957440329, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node37105584', child_labels=[], edge_length=15.2096070352, taxon_label='Morelia viridis'),
            'Node37105232' : NodeRelationship(parent_label='Node37105584', child_labels=['Morelia carinata','Node37163248'], edge_length=1.5569144386, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node37105232', child_labels=[], edge_length=13.6526925966, taxon_label='Morelia carinata'),
            'Node37163248' : NodeRelationship(parent_label='Node37105232', child_labels=['Morelia spilota','Morelia bredli'], edge_length=8.73106140586, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37163248', child_labels=[], edge_length=4.92163119076, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37163248', child_labels=[], edge_length=4.92163119076, taxon_label='Morelia bredli'),
            'Node37103280' : NodeRelationship(parent_label='Node37102896', child_labels=['Node37163312','Morelia oenpelliensis'], edge_length=1.53468865192, taxon_label=None),
            'Node37163312' : NodeRelationship(parent_label='Node37103280', child_labels=['Python reticulatus','Python timoriensis'], edge_length=6.61343105392, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37163312', child_labels=[], edge_length=14.7952666012, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37163312', child_labels=[], edge_length=14.7952666012, taxon_label='Python timoriensis'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37103280', child_labels=[], edge_length=21.4086976551, taxon_label='Morelia oenpelliensis'),
            'Node37163504' : NodeRelationship(parent_label='Node37714672', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=20.5343771273, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37163504', child_labels=[], edge_length=6.63522564694, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37163504', child_labels=[], edge_length=6.63522564694, taxon_label='Aspidites melanocephalus'),
            'Node37163696' : NodeRelationship(parent_label='Node37103408', child_labels=['Python molurus','Node37164080'], edge_length=10.1890732074, taxon_label=None),
            'Python molurus' : NodeRelationship(parent_label='Node37163696', child_labels=[], edge_length=23.4003588995, taxon_label='Python molurus'),
            'Node37164080' : NodeRelationship(parent_label='Node37163696', child_labels=['Python brongersmai','Python regius'], edge_length=0.0534325174359, taxon_label=None),
            'Python brongersmai' : NodeRelationship(parent_label='Node37164080', child_labels=[], edge_length=23.3469263821, taxon_label='Python brongersmai'),
            'Python regius' : NodeRelationship(parent_label='Node37164080', child_labels=[], edge_length=23.3469263821, taxon_label='Python regius'),
            'Python sebae' : NodeRelationship(parent_label='Node37103344', child_labels=[], edge_length=42.025995763, taxon_label='Python sebae'),
        },
        {
            'Node37164432' : NodeRelationship(parent_label=None, child_labels=['Node37164496','Python sebae'], edge_length=None, taxon_label=None),
            'Node37164496' : NodeRelationship(parent_label='Node37164432', child_labels=['Node37164368','Node37220976'], edge_length=4.51223406477, taxon_label=None),
            'Node37164368' : NodeRelationship(parent_label='Node37164496', child_labels=['Node37775760','Node37220432'], edge_length=3.29007229288, taxon_label=None),
            'Node37775760' : NodeRelationship(parent_label='Node37164368', child_labels=['Node37164304','Node37165712'], edge_length=4.09158345085, taxon_label=None),
            'Node37164304' : NodeRelationship(parent_label='Node37775760', child_labels=['Morelia viridis','Node37164848'], edge_length=2.82624159524, taxon_label=None),
            'Morelia viridis' : NodeRelationship(parent_label='Node37164304', child_labels=[], edge_length=14.5888315231, taxon_label='Morelia viridis'),
            'Node37164848' : NodeRelationship(parent_label='Node37164304', child_labels=['Node37164720','Node37165264'], edge_length=0.0149913671254, taxon_label=None),
            'Node37164720' : NodeRelationship(parent_label='Node37164848', child_labels=['Node37164880','Morelia carinata'], edge_length=4.89600409139, taxon_label=None),
            'Node37164880' : NodeRelationship(parent_label='Node37164720', child_labels=['Morelia spilota','Morelia bredli'], edge_length=5.06534441374, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37164880', child_labels=[], edge_length=4.61249165082, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37164880', child_labels=[], edge_length=4.61249165082, taxon_label='Morelia bredli'),
            'Morelia carinata' : NodeRelationship(parent_label='Node37164720', child_labels=[], edge_length=9.67783606456, taxon_label='Morelia carinata'),
            'Node37165264' : NodeRelationship(parent_label='Node37164848', child_labels=['Antaresia maculosa','Node37165328'], edge_length=0.516331095919, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37165264', child_labels=[], edge_length=14.05750906, taxon_label='Antaresia maculosa'),
            'Node37165328' : NodeRelationship(parent_label='Node37165264', child_labels=['Antaresia perthensis','Node37165456'], edge_length=6.44901536169, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37165328', child_labels=[], edge_length=7.60849369834, taxon_label='Antaresia perthensis'),
            'Node37165456' : NodeRelationship(parent_label='Node37165328', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=5.74491098314, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37165456', child_labels=[], edge_length=1.86358271521, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37165456', child_labels=[], edge_length=1.86358271521, taxon_label='Antaresia stimsoni'),
            'Node37165712' : NodeRelationship(parent_label='Node37775760', child_labels=['Node37165616','Node37165968'], edge_length=0.0769064358513, taxon_label=None),
            'Node37165616' : NodeRelationship(parent_label='Node37165712', child_labels=['Python reticulatus','Python timoriensis'], edge_length=8.7380438705, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37165616', child_labels=[], edge_length=8.60012281197, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37165616', child_labels=[], edge_length=8.60012281197, taxon_label='Python timoriensis'),
            'Node37165968' : NodeRelationship(parent_label='Node37165712', child_labels=['Node37165840','Node37166576'], edge_length=0.894040812546, taxon_label=None),
            'Node37165840' : NodeRelationship(parent_label='Node37165968', child_labels=['Node37165808','Node37166384'], edge_length=1.98993837612, taxon_label=None),
            'Node37165808' : NodeRelationship(parent_label='Node37165840', child_labels=['Apodora papuana','Node37166192'], edge_length=1.38724843406, taxon_label=None),
            'Apodora papuana' : NodeRelationship(parent_label='Node37165808', child_labels=[], edge_length=13.0669390597, taxon_label='Apodora papuana'),
            'Node37166192' : NodeRelationship(parent_label='Node37165808', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=3.7285329013, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node37166192', child_labels=[], edge_length=9.33840615844, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node37166192', child_labels=[], edge_length=9.33840615844, taxon_label='Bothrochilus boa'),
            'Node37166384' : NodeRelationship(parent_label='Node37165840', child_labels=['Liasis olivaceus','Node37166512'], edge_length=3.48386820234, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37166384', child_labels=[], edge_length=10.9703192915, taxon_label='Liasis olivaceus'),
            'Node37166512' : NodeRelationship(parent_label='Node37166384', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=3.83596797183, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37166512', child_labels=[], edge_length=7.13435131963, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37166512', child_labels=[], edge_length=7.13435131963, taxon_label='Liasis mackloti'),
            'Node37166576' : NodeRelationship(parent_label='Node37165968', child_labels=['Node37166544','Morelia oenpelliensis'], edge_length=2.85390092474, taxon_label=None),
            'Node37166544' : NodeRelationship(parent_label='Node37166576', child_labels=['Node37166704','Morelia boeleni'], edge_length=3.52728125011, taxon_label=None),
            'Node37166704' : NodeRelationship(parent_label='Node37166544', child_labels=['Morelia tracyae','Node37166960'], edge_length=1.7419442928, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37166704', child_labels=[], edge_length=8.32099940228, taxon_label='Morelia tracyae'),
            'Node37166960' : NodeRelationship(parent_label='Node37166704', child_labels=['Node37166896','Morelia amethistina'], edge_length=4.44695529056, taxon_label=None),
            'Node37166896' : NodeRelationship(parent_label='Node37166960', child_labels=['Node37166992','Morelia clastolepis'], edge_length=2.02523285213, taxon_label=None),
            'Node37166992' : NodeRelationship(parent_label='Node37166896', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=0.774826953003, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37166992', child_labels=[], edge_length=1.07398430658, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node37166992', child_labels=[], edge_length=1.07398430658, taxon_label='Morelia nauta'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37166896', child_labels=[], edge_length=1.84881125958, taxon_label='Morelia clastolepis'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37166960', child_labels=[], edge_length=3.87404411171, taxon_label='Morelia amethistina'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node37166544', child_labels=[], edge_length=10.0629436951, taxon_label='Morelia boeleni'),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37166576', child_labels=[], edge_length=13.5902249452, taxon_label='Morelia oenpelliensis'),
            'Node37220432' : NodeRelationship(parent_label='Node37164368', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=12.3527736124, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37220432', child_labels=[], edge_length=9.1538829568, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37220432', child_labels=[], edge_length=9.1538829568, taxon_label='Aspidites melanocephalus'),
            'Node37220976' : NodeRelationship(parent_label='Node37164496', child_labels=['Node37220944','Python molurus'], edge_length=1.61144541412, taxon_label=None),
            'Node37220944' : NodeRelationship(parent_label='Node37220976', child_labels=['Python regius','Python brongersmai'], edge_length=0.5196540118, taxon_label=None),
            'Python regius' : NodeRelationship(parent_label='Node37220944', child_labels=[], edge_length=22.6656294361, taxon_label='Python regius'),
            'Python brongersmai' : NodeRelationship(parent_label='Node37220944', child_labels=[], edge_length=22.6656294361, taxon_label='Python brongersmai'),
            'Python molurus' : NodeRelationship(parent_label='Node37220976', child_labels=[], edge_length=23.1852834479, taxon_label='Python molurus'),
            'Python sebae' : NodeRelationship(parent_label='Node37164432', child_labels=[], edge_length=29.3089629268, taxon_label='Python sebae'),
        },
        {
            'Node37221424' : NodeRelationship(parent_label=None, child_labels=['Python sebae','Node37221584'], edge_length=None, taxon_label=None),
            'Python sebae' : NodeRelationship(parent_label='Node37221424', child_labels=[], edge_length=41.037475788, taxon_label='Python sebae'),
            'Node37221584' : NodeRelationship(parent_label='Node37221424', child_labels=['Node37221616','Node37221936'], edge_length=0.971031599202, taxon_label=None),
            'Node37221616' : NodeRelationship(parent_label='Node37221584', child_labels=['Python regius','Node37221776'], edge_length=7.88042810387, taxon_label=None),
            'Python regius' : NodeRelationship(parent_label='Node37221616', child_labels=[], edge_length=32.1860160849, taxon_label='Python regius'),
            'Node37221776' : NodeRelationship(parent_label='Node37221616', child_labels=['Python molurus','Python brongersmai'], edge_length=3.37816716317, taxon_label=None),
            'Python molurus' : NodeRelationship(parent_label='Node37221776', child_labels=[], edge_length=28.8078489217, taxon_label='Python molurus'),
            'Python brongersmai' : NodeRelationship(parent_label='Node37221776', child_labels=[], edge_length=28.8078489217, taxon_label='Python brongersmai'),
            'Node37221936' : NodeRelationship(parent_label='Node37221584', child_labels=['Node37221968','Node37224432'], edge_length=11.1348383801, taxon_label=None),
            'Node37221968' : NodeRelationship(parent_label='Node37221936', child_labels=['Node37221872','Node37222736'], edge_length=5.61233091052, taxon_label=None),
            'Node37221872' : NodeRelationship(parent_label='Node37221968', child_labels=['Morelia oenpelliensis','Node37222224'], edge_length=3.04566241261, taxon_label=None),
            'Morelia oenpelliensis' : NodeRelationship(parent_label='Node37221872', child_labels=[], edge_length=20.2736124855, taxon_label='Morelia oenpelliensis'),
            'Node37222224' : NodeRelationship(parent_label='Node37221872', child_labels=['Node37222096','Morelia boeleni'], edge_length=7.80824403904, taxon_label=None),
            'Node37222096' : NodeRelationship(parent_label='Node37222224', child_labels=['Morelia tracyae','Node37222416'], edge_length=3.53762335708, taxon_label=None),
            'Morelia tracyae' : NodeRelationship(parent_label='Node37222096', child_labels=[], edge_length=8.92774508942, taxon_label='Morelia tracyae'),
            'Node37222416' : NodeRelationship(parent_label='Node37222096', child_labels=['Node37222288','Morelia amethistina'], edge_length=0.98740476198, taxon_label=None),
            'Node37222288' : NodeRelationship(parent_label='Node37222416', child_labels=['Node37222448','Morelia clastolepis'], edge_length=4.46813855011, taxon_label=None),
            'Node37222448' : NodeRelationship(parent_label='Node37222288', child_labels=['Morelia kinghorni','Morelia nauta'], edge_length=0.318849288103, taxon_label=None),
            'Morelia kinghorni' : NodeRelationship(parent_label='Node37222448', child_labels=[], edge_length=3.15335248923, taxon_label='Morelia kinghorni'),
            'Morelia nauta' : NodeRelationship(parent_label='Node37222448', child_labels=[], edge_length=3.15335248923, taxon_label='Morelia nauta'),
            'Morelia clastolepis' : NodeRelationship(parent_label='Node37222288', child_labels=[], edge_length=3.47220177733, taxon_label='Morelia clastolepis'),
            'Morelia amethistina' : NodeRelationship(parent_label='Node37222416', child_labels=[], edge_length=7.94034032744, taxon_label='Morelia amethistina'),
            'Morelia boeleni' : NodeRelationship(parent_label='Node37222224', child_labels=[], edge_length=12.4653684465, taxon_label='Morelia boeleni'),
            'Node37222736' : NodeRelationship(parent_label='Node37221968', child_labels=['Node37222960','Node37223632'], edge_length=1.63793957099, taxon_label=None),
            'Node37222960' : NodeRelationship(parent_label='Node37222736', child_labels=['Node37222864','Node37223248'], edge_length=1.17144499717, taxon_label=None),
            'Node37222864' : NodeRelationship(parent_label='Node37222960', child_labels=['Liasis olivaceus','Node37223184'], edge_length=6.09401457319, taxon_label=None),
            'Liasis olivaceus' : NodeRelationship(parent_label='Node37222864', child_labels=[], edge_length=14.4158757568, taxon_label='Liasis olivaceus'),
            'Node37223184' : NodeRelationship(parent_label='Node37222864', child_labels=['Liasis fuscus','Liasis mackloti'], edge_length=10.1133978359, taxon_label=None),
            'Liasis fuscus' : NodeRelationship(parent_label='Node37223184', child_labels=[], edge_length=4.30247792095, taxon_label='Liasis fuscus'),
            'Liasis mackloti' : NodeRelationship(parent_label='Node37223184', child_labels=[], edge_length=4.30247792095, taxon_label='Liasis mackloti'),
            'Node37223248' : NodeRelationship(parent_label='Node37222960', child_labels=['Node37223376','Apodora papuana'], edge_length=0.168673432469, taxon_label=None),
            'Node37223376' : NodeRelationship(parent_label='Node37223248', child_labels=['Python reticulatus','Python timoriensis'], edge_length=9.21506648422, taxon_label=None),
            'Python reticulatus' : NodeRelationship(parent_label='Node37223376', child_labels=[], edge_length=11.1261504133, taxon_label='Python reticulatus'),
            'Python timoriensis' : NodeRelationship(parent_label='Node37223376', child_labels=[], edge_length=11.1261504133, taxon_label='Python timoriensis'),
            'Apodora papuana' : NodeRelationship(parent_label='Node37223248', child_labels=[], edge_length=20.3412168975, taxon_label='Apodora papuana'),
            'Node37223632' : NodeRelationship(parent_label='Node37222736', child_labels=['Node37223408','Node37223760'], edge_length=0.325200272563, taxon_label=None),
            'Node37223408' : NodeRelationship(parent_label='Node37223632', child_labels=['Leiopython albertisii','Bothrochilus boa'], edge_length=4.80149248356, taxon_label=None),
            'Leiopython albertisii' : NodeRelationship(parent_label='Node37223408', child_labels=[], edge_length=16.554642571, taxon_label='Leiopython albertisii'),
            'Bothrochilus boa' : NodeRelationship(parent_label='Node37223408', child_labels=[], edge_length=16.554642571, taxon_label='Bothrochilus boa'),
            'Node37223760' : NodeRelationship(parent_label='Node37223632', child_labels=['Node37223728','Node37224304'], edge_length=5.59816192539, taxon_label=None),
            'Node37223728' : NodeRelationship(parent_label='Node37223760', child_labels=['Node37223888','Morelia viridis'], edge_length=1.51569357225, taxon_label=None),
            'Node37223888' : NodeRelationship(parent_label='Node37223728', child_labels=['Morelia carinata','Node37224208'], edge_length=0.59310677926, taxon_label=None),
            'Morelia carinata' : NodeRelationship(parent_label='Node37223888', child_labels=[], edge_length=13.6491727777, taxon_label='Morelia carinata'),
            'Node37224208' : NodeRelationship(parent_label='Node37223888', child_labels=['Morelia spilota','Morelia bredli'], edge_length=7.69541700903, taxon_label=None),
            'Morelia spilota' : NodeRelationship(parent_label='Node37224208', child_labels=[], edge_length=5.95375576867, taxon_label='Morelia spilota'),
            'Morelia bredli' : NodeRelationship(parent_label='Node37224208', child_labels=[], edge_length=5.95375576867, taxon_label='Morelia bredli'),
            'Morelia viridis' : NodeRelationship(parent_label='Node37223728', child_labels=[], edge_length=14.242279557, taxon_label='Morelia viridis'),
            'Node37224304' : NodeRelationship(parent_label='Node37223760', child_labels=['Antaresia maculosa','Node37224368'], edge_length=2.47499227209, taxon_label=None),
            'Antaresia maculosa' : NodeRelationship(parent_label='Node37224304', child_labels=[], edge_length=13.2829808571, taxon_label='Antaresia maculosa'),
            'Node37224368' : NodeRelationship(parent_label='Node37224304', child_labels=['Antaresia perthensis','Node37277904'], edge_length=3.47453380669, taxon_label=None),
            'Antaresia perthensis' : NodeRelationship(parent_label='Node37224368', child_labels=[], edge_length=9.80844705042, taxon_label='Antaresia perthensis'),
            'Node37277904' : NodeRelationship(parent_label='Node37224368', child_labels=['Antaresia childreni','Antaresia stimsoni'], edge_length=3.99697653657, taxon_label=None),
            'Antaresia childreni' : NodeRelationship(parent_label='Node37277904', child_labels=[], edge_length=5.81147051386, taxon_label='Antaresia childreni'),
            'Antaresia stimsoni' : NodeRelationship(parent_label='Node37277904', child_labels=[], edge_length=5.81147051386, taxon_label='Antaresia stimsoni'),
            'Node37224432' : NodeRelationship(parent_label='Node37221936', child_labels=['Aspidites ramsayi','Aspidites melanocephalus'], edge_length=19.6159722746, taxon_label=None),
            'Aspidites ramsayi' : NodeRelationship(parent_label='Node37224432', child_labels=[], edge_length=9.3156335341, taxon_label='Aspidites ramsayi'),
            'Aspidites melanocephalus' : NodeRelationship(parent_label='Node37224432', child_labels=[], edge_length=9.3156335341, taxon_label='Aspidites melanocephalus'),
        },
    ]
    return treelist_node_references

def reference_tree_list(taxon_set=None):
    tree_list = dendropy.TreeList(label=None, oid="TreeList4523184", taxon_set=taxon_set)
    tax_36960624 = tree_list.taxon_set.require_taxon(label="Antaresia childreni", oid="Taxon36960624")
    tax_36960560 = tree_list.taxon_set.require_taxon(label="Antaresia maculosa", oid="Taxon36960560")
    tax_36960592 = tree_list.taxon_set.require_taxon(label="Antaresia perthensis", oid="Taxon36960592")
    tax_36960432 = tree_list.taxon_set.require_taxon(label="Antaresia stimsoni", oid="Taxon36960432")
    tax_36960656 = tree_list.taxon_set.require_taxon(label="Aspidites melanocephalus", oid="Taxon36960656")
    tax_36960528 = tree_list.taxon_set.require_taxon(label="Aspidites ramsayi", oid="Taxon36960528")
    tax_36960784 = tree_list.taxon_set.require_taxon(label="Bothrochilus boa", oid="Taxon36960784")
    tax_36960688 = tree_list.taxon_set.require_taxon(label="Leiopython albertisii", oid="Taxon36960688")
    tax_36960720 = tree_list.taxon_set.require_taxon(label="Liasis fuscus", oid="Taxon36960720")
    tax_36960912 = tree_list.taxon_set.require_taxon(label="Liasis mackloti", oid="Taxon36960912")
    tax_36960816 = tree_list.taxon_set.require_taxon(label="Liasis olivaceus", oid="Taxon36960816")
    tax_36960848 = tree_list.taxon_set.require_taxon(label="Apodora papuana", oid="Taxon36960848")
    tax_36960880 = tree_list.taxon_set.require_taxon(label="Morelia amethistina", oid="Taxon36960880")
    tax_36960944 = tree_list.taxon_set.require_taxon(label="Morelia boeleni", oid="Taxon36960944")
    tax_36961008 = tree_list.taxon_set.require_taxon(label="Morelia bredli", oid="Taxon36961008")
    tax_36961072 = tree_list.taxon_set.require_taxon(label="Morelia carinata", oid="Taxon36961072")
    tax_36961264 = tree_list.taxon_set.require_taxon(label="Morelia clastolepis", oid="Taxon36961264")
    tax_36960752 = tree_list.taxon_set.require_taxon(label="Morelia kinghorni", oid="Taxon36960752")
    tax_36961328 = tree_list.taxon_set.require_taxon(label="Morelia nauta", oid="Taxon36961328")
    tax_36961456 = tree_list.taxon_set.require_taxon(label="Morelia oenpelliensis", oid="Taxon36961456")
    tax_36961200 = tree_list.taxon_set.require_taxon(label="Morelia spilota", oid="Taxon36961200")
    tax_36961392 = tree_list.taxon_set.require_taxon(label="Morelia tracyae", oid="Taxon36961392")
    tax_36961136 = tree_list.taxon_set.require_taxon(label="Morelia viridis", oid="Taxon36961136")
    tax_36961584 = tree_list.taxon_set.require_taxon(label="Python brongersmai", oid="Taxon36961584")
    tax_36961648 = tree_list.taxon_set.require_taxon(label="Python molurus", oid="Taxon36961648")
    tax_36961712 = tree_list.taxon_set.require_taxon(label="Python regius", oid="Taxon36961712")
    tax_36961904 = tree_list.taxon_set.require_taxon(label="Python reticulatus", oid="Taxon36961904")
    tax_36961776 = tree_list.taxon_set.require_taxon(label="Python sebae", oid="Taxon36961776")
    tax_36961520 = tree_list.taxon_set.require_taxon(label="Python timoriensis", oid="Taxon36961520")
    tree_36962288 = dendropy.Tree(label="Tree01", taxon_set=tree_list.taxon_set, oid="Tree36962288")
    tree_list.append(tree_36962288, reindex_taxa=False)
    tree_36962288.seed_node.oid = 'Node36995184'
    nd_36995248 = tree_36962288.seed_node.new_child(label="Node36995248", taxon=None, edge_length=61.9915926669, oid="Node36995248")
    nd_36995248.edge.oid = "Edge36995280"
    nd_36962224 = tree_36962288.seed_node.new_child(label="Node36962224", taxon=None, edge_length=29.8452823467, oid="Node36962224")
    nd_36962224.edge.oid = "Edge36995376"
    nd_36961968 = nd_36995248.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=5.50104114707, oid="Node36961968")
    nd_36961968.edge.oid = "Edge36995312"
    nd_36962096 = nd_36995248.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=5.50104114707, oid="Node36962096")
    nd_36962096.edge.oid = "Edge36995408"
    nd_36995344 = nd_36962224.new_child(label="Node36995344", taxon=None, edge_length=0.377327561902, oid="Node36995344")
    nd_36995344.edge.oid = "Edge36995504"
    nd_36997392 = nd_36962224.new_child(label="Node36997392", taxon=None, edge_length=3.61657514544, oid="Node36997392")
    nd_36997392.edge.oid = "Edge36997680"
    nd_36995440 = nd_36995344.new_child(label="Node36995440", taxon=None, edge_length=33.0195611855, oid="Node36995440")
    nd_36995440.edge.oid = "Edge36995568"
    nd_36996528 = nd_36995344.new_child(label="Node36996528", taxon=None, edge_length=23.9652741942, oid="Node36996528")
    nd_36996528.edge.oid = "Edge36996592"
    nd_36995472 = nd_36995440.new_child(label="Node36995472", taxon=None, edge_length=1.09375666234, oid="Node36995472")
    nd_36995472.edge.oid = "Edge36995632"
    nd_36996400 = nd_36995440.new_child(label="Node36996400", taxon=None, edge_length=3.66513378189, oid="Node36996400")
    nd_36996400.edge.oid = "Edge36996048"
    nd_36995536 = nd_36995472.new_child(label="Node36995536", taxon=None, edge_length=0.351619575146, oid="Node36995536")
    nd_36995536.edge.oid = "Edge36995696"
    nd_36995984 = nd_36995472.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=3.15670605753, oid="Node36995984")
    nd_36995984.edge.oid = "Edge36996240"
    nd_36995600 = nd_36995536.new_child(label="Node36995600", taxon=None, edge_length=2.09759626003, oid="Node36995600")
    nd_36995600.edge.oid = "Edge36995760"
    nd_36995952 = nd_36995536.new_child(label="Node36995952", taxon=None, edge_length=1.73551661307, oid="Node36995952")
    nd_36995952.edge.oid = "Edge36995888"
    nd_36995664 = nd_36995600.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=0.707490222359, oid="Node36995664")
    nd_36995664.edge.oid = "Edge36995824"
    nd_36995920 = nd_36995600.new_child(label="Node36995920", taxon=None, edge_length=0.470505730507, oid="Node36995920")
    nd_36995920.edge.oid = "Edge36995728"
    nd_36995792 = nd_36995920.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=0.236984491852, oid="Node36995792")
    nd_36995792.edge.oid = "Edge36995856"
    nd_36996016 = nd_36995920.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=0.236984491852, oid="Node36996016")
    nd_36996016.edge.oid = "Edge36996080"
    nd_36996144 = nd_36995952.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=1.06956986931, oid="Node36996144")
    nd_36996144.edge.oid = "Edge36996112"
    nd_36996208 = nd_36995952.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=1.06956986931, oid="Node36996208")
    nd_36996208.edge.oid = "Edge36996272"
    nd_36996368 = nd_36996400.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=0.585328937986, oid="Node36996368")
    nd_36996368.edge.oid = "Edge36996336"
    nd_36996496 = nd_36996400.new_child(label="Python regius", taxon=tax_36961712, edge_length=0.585328937986, oid="Node36996496")
    nd_36996496.edge.oid = "Edge36996432"
    nd_36996176 = nd_36996528.new_child(label="Node36996176", taxon=None, edge_length=3.49206820101, oid="Node36996176")
    nd_36996176.edge.oid = "Edge36996304"
    nd_36997424 = nd_36996528.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=13.3047497112, oid="Node36997424")
    nd_36997424.edge.oid = "Edge36997520"
    nd_36996560 = nd_36996176.new_child(label="Node36996560", taxon=None, edge_length=4.90381799328, oid="Node36996560")
    nd_36996560.edge.oid = "Edge36996656"
    nd_36997168 = nd_36996176.new_child(label="Node36997168", taxon=None, edge_length=5.39574605108, oid="Node36997168")
    nd_36997168.edge.oid = "Edge36997104"
    nd_36996464 = nd_36996560.new_child(label="Node36996464", taxon=None, edge_length=3.58364386958, oid="Node36996464")
    nd_36996464.edge.oid = "Edge36996720"
    nd_36997040 = nd_36996560.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=4.9088635169, oid="Node36997040")
    nd_36997040.edge.oid = "Edge36996688"
    nd_36996624 = nd_36996464.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=1.32521964733, oid="Node36996624")
    nd_36996624.edge.oid = "Edge36996784"
    nd_36996880 = nd_36996464.new_child(label="Node36996880", taxon=None, edge_length=0.699088376552, oid="Node36996880")
    nd_36996880.edge.oid = "Edge36996816"
    nd_36996752 = nd_36996880.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=0.626131270775, oid="Node36996752")
    nd_36996752.edge.oid = "Edge36996848"
    nd_36997008 = nd_36996880.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=0.626131270775, oid="Node36997008")
    nd_36997008.edge.oid = "Edge36996944"
    nd_36997072 = nd_36997168.new_child(label="Node36997072", taxon=None, edge_length=3.89664227731, oid="Node36997072")
    nd_36997072.edge.oid = "Edge36996976"
    nd_36997360 = nd_36997168.new_child(label="Node36997360", taxon=None, edge_length=4.2237066393, oid="Node36997360")
    nd_36997360.edge.oid = "Edge36997264"
    nd_36996912 = nd_36997072.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=0.520293181797, oid="Node36996912")
    nd_36996912.edge.oid = "Edge36997232"
    nd_36997328 = nd_36997072.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=0.520293181797, oid="Node36997328")
    nd_36997328.edge.oid = "Edge36997136"
    nd_36997296 = nd_36997360.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=0.193228819804, oid="Node36997296")
    nd_36997296.edge.oid = "Edge36997200"
    nd_36997488 = nd_36997360.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=0.193228819804, oid="Node36997488")
    nd_36997488.edge.oid = "Edge36997552"
    nd_36997616 = nd_36997392.new_child(label="Node36997616", taxon=None, edge_length=21.1685275644, oid="Node36997616")
    nd_36997616.edge.oid = "Edge36997584"
    nd_36998064 = nd_36997392.new_child(label="Node36998064", taxon=None, edge_length=26.4364769675, oid="Node36998064")
    nd_36998064.edge.oid = "Edge36998192"
    nd_36997648 = nd_36997616.new_child(label="Python molurus", taxon=tax_36961648, edge_length=12.8622487575, oid="Node36997648")
    nd_36997648.edge.oid = "Edge36997744"
    nd_36997840 = nd_36997616.new_child(label="Node36997840", taxon=None, edge_length=7.83034270743, oid="Node36997840")
    nd_36997840.edge.oid = "Edge36997776"
    nd_36997712 = nd_36997840.new_child(label="Node36997712", taxon=None, edge_length=3.83058430268, oid="Node36997712")
    nd_36997712.edge.oid = "Edge36997808"
    nd_36998128 = nd_36997840.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=5.03190605005, oid="Node36998128")
    nd_36998128.edge.oid = "Edge36998096"
    nd_36997872 = nd_36997712.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=1.20132174738, oid="Node36997872")
    nd_36997872.edge.oid = "Edge36997936"
    nd_36998000 = nd_36997712.new_child(label="Python sebae", taxon=tax_36961776, edge_length=1.20132174738, oid="Node36998000")
    nd_36998000.edge.oid = "Edge36997456"
    nd_36997968 = nd_36998064.new_child(label="Node36997968", taxon=None, edge_length=4.80697376221, oid="Node36997968")
    nd_36997968.edge.oid = "Edge36997904"
    nd_36998704 = nd_36998064.new_child(label="Node36998704", taxon=None, edge_length=5.17312138897, oid="Node36998704")
    nd_36998704.edge.oid = "Edge36998608"
    nd_36998160 = nd_36997968.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=2.78732559217, oid="Node36998160")
    nd_36998160.edge.oid = "Edge36998256"
    nd_36998320 = nd_36997968.new_child(label="Node36998320", taxon=None, edge_length=1.381705203, oid="Node36998320")
    nd_36998320.edge.oid = "Edge36998384"
    nd_36998032 = nd_36998320.new_child(label="Node36998032", taxon=None, edge_length=1.30019225223, oid="Node36998032")
    nd_36998032.edge.oid = "Edge36998352"
    nd_36998416 = nd_36998320.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=1.40562038917, oid="Node36998416")
    nd_36998416.edge.oid = "Edge36998480"
    nd_36998224 = nd_36998032.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=0.105428136944, oid="Node36998224")
    nd_36998224.edge.oid = "Edge36998448"
    nd_36998512 = nd_36998032.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=0.105428136944, oid="Node36998512")
    nd_36998512.edge.oid = "Edge36998576"
    nd_36998288 = nd_36998704.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=2.42117796541, oid="Node36998288")
    nd_36998288.edge.oid = "Edge36998672"
    nd_36998800 = nd_36998704.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=2.42117796541, oid="Node36998800")
    nd_36998800.edge.oid = "Edge36998544"
    tree_36998864 = dendropy.Tree(label="Tree02", taxon_set=tree_list.taxon_set, oid="Tree36998864")
    tree_list.append(tree_36998864, reindex_taxa=False)
    tree_36998864.seed_node.oid = 'Node36995120'
    nd_36998960 = tree_36998864.seed_node.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=26.0384706379, oid="Node36998960")
    nd_36998960.edge.oid = "Edge36998992"
    nd_36999056 = tree_36998864.seed_node.new_child(label="Node36999056", taxon=None, edge_length=5.21927124134, oid="Node36999056")
    nd_36999056.edge.oid = "Edge36999088"
    nd_36995152 = nd_36999056.new_child(label="Node36995152", taxon=None, edge_length=0.952265882909, oid="Node36995152")
    nd_36995152.edge.oid = "Edge36998768"
    nd_37062192 = nd_36999056.new_child(label="Node37062192", taxon=None, edge_length=0.340827566292, oid="Node37062192")
    nd_37062192.edge.oid = "Edge37062320"
    nd_36962160 = nd_36995152.new_child(label="Node36962160", taxon=None, edge_length=1.64883491058, oid="Node36962160")
    nd_36962160.edge.oid = "Edge36999152"
    nd_37061744 = nd_36995152.new_child(label="Node37061744", taxon=None, edge_length=7.07533845302, oid="Node37061744")
    nd_37061744.edge.oid = "Edge37061360"
    nd_36999024 = nd_36962160.new_child(label="Node36999024", taxon=None, edge_length=2.95693808858, oid="Node36999024")
    nd_36999024.edge.oid = "Edge37060688"
    nd_37060720 = nd_36962160.new_child(label="Node37060720", taxon=None, edge_length=2.15553698471, oid="Node37060720")
    nd_37060720.edge.oid = "Edge37060848"
    nd_36999120 = nd_36999024.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=15.2611605145, oid="Node36999120")
    nd_36999120.edge.oid = "Edge37060752"
    nd_37060816 = nd_36999024.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=15.2611605145, oid="Node37060816")
    nd_37060816.edge.oid = "Edge37060880"
    nd_37060784 = nd_37060720.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=16.0625616183, oid="Node37060784")
    nd_37060784.edge.oid = "Edge37060656"
    nd_37061040 = nd_37060720.new_child(label="Node37061040", taxon=None, edge_length=1.37038749174, oid="Node37061040")
    nd_37061040.edge.oid = "Edge37060976"
    nd_37060944 = nd_37061040.new_child(label="Node37060944", taxon=None, edge_length=0.166805015451, oid="Node37060944")
    nd_37060944.edge.oid = "Edge37061008"
    nd_37061840 = nd_37061040.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=14.6921741266, oid="Node37061840")
    nd_37061840.edge.oid = "Edge37061808"
    nd_37061072 = nd_37060944.new_child(label="Node37061072", taxon=None, edge_length=7.4747693113, oid="Node37061072")
    nd_37061072.edge.oid = "Edge37061136"
    nd_37061392 = nd_37060944.new_child(label="Node37061392", taxon=None, edge_length=1.43991095233, oid="Node37061392")
    nd_37061392.edge.oid = "Edge37061264"
    nd_37060912 = nd_37061072.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=7.05059979983, oid="Node37060912")
    nd_37060912.edge.oid = "Edge37061200"
    nd_37061296 = nd_37061072.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=7.05059979983, oid="Node37061296")
    nd_37061296.edge.oid = "Edge37061168"
    nd_37061232 = nd_37061392.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=13.0854581588, oid="Node37061232")
    nd_37061232.edge.oid = "Edge37061328"
    nd_37061456 = nd_37061392.new_child(label="Node37061456", taxon=None, edge_length=2.85871098218, oid="Node37061456")
    nd_37061456.edge.oid = "Edge37061488"
    nd_37061520 = nd_37061456.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=10.2267471766, oid="Node37061520")
    nd_37061520.edge.oid = "Edge37061424"
    nd_37061584 = nd_37061456.new_child(label="Node37061584", taxon=None, edge_length=5.35604513133, oid="Node37061584")
    nd_37061584.edge.oid = "Edge37061616"
    nd_37061648 = nd_37061584.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=4.8707020453, oid="Node37061648")
    nd_37061648.edge.oid = "Edge37061552"
    nd_37061712 = nd_37061584.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=4.8707020453, oid="Node37061712")
    nd_37061712.edge.oid = "Edge37061776"
    nd_37061872 = nd_37061744.new_child(label="Node37061872", taxon=None, edge_length=1.94094943657, oid="Node37061872")
    nd_37061872.edge.oid = "Edge37061904"
    nd_37062096 = nd_37061744.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=12.7915950606, oid="Node37062096")
    nd_37062096.edge.oid = "Edge37062256"
    nd_37061680 = nd_37061872.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=10.850645624, oid="Node37061680")
    nd_37061680.edge.oid = "Edge37061968"
    nd_37062064 = nd_37061872.new_child(label="Node37062064", taxon=None, edge_length=8.02016262021, oid="Node37062064")
    nd_37062064.edge.oid = "Edge37061104"
    nd_37061936 = nd_37062064.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=2.83048300382, oid="Node37061936")
    nd_37061936.edge.oid = "Edge37062000"
    nd_37062160 = nd_37062064.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=2.83048300382, oid="Node37062160")
    nd_37062160.edge.oid = "Edge37062224"
    nd_37062288 = nd_37062192.new_child(label="Node37062288", taxon=None, edge_length=1.42657874327, oid="Node37062288")
    nd_37062288.edge.oid = "Edge37062128"
    nd_37063088 = nd_37062192.new_child(label="Node37063088", taxon=None, edge_length=2.1717154248, oid="Node37063088")
    nd_37063088.edge.oid = "Edge37063376"
    nd_37062032 = nd_37062288.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=19.051793087, oid="Node37062032")
    nd_37062032.edge.oid = "Edge37062416"
    nd_37062512 = nd_37062288.new_child(label="Node37062512", taxon=None, edge_length=1.77859521326, oid="Node37062512")
    nd_37062512.edge.oid = "Edge37062384"
    nd_37062448 = nd_37062512.new_child(label="Node37062448", taxon=None, edge_length=7.19864380037, oid="Node37062448")
    nd_37062448.edge.oid = "Edge37062480"
    nd_37063184 = nd_37062512.new_child(label="Node37063184", taxon=None, edge_length=10.657346381, oid="Node37063184")
    nd_37063184.edge.oid = "Edge37063056"
    nd_37062544 = nd_37062448.new_child(label="Node37062544", taxon=None, edge_length=1.9541248425, oid="Node37062544")
    nd_37062544.edge.oid = "Edge37062608"
    nd_37062928 = nd_37062448.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=10.0745540733, oid="Node37062928")
    nd_37062928.edge.oid = "Edge37062992"
    nd_37062352 = nd_37062544.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=8.12042923083, oid="Node37062352")
    nd_37062352.edge.oid = "Edge37062672"
    nd_37062768 = nd_37062544.new_child(label="Node37062768", taxon=None, edge_length=6.05470978372, oid="Node37062768")
    nd_37062768.edge.oid = "Edge37062640"
    nd_37062576 = nd_37062768.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=2.0657194471, oid="Node37062576")
    nd_37062576.edge.oid = "Edge37062704"
    nd_37062896 = nd_37062768.new_child(label="Node37062896", taxon=None, edge_length=1.06606772283, oid="Node37062896")
    nd_37062896.edge.oid = "Edge37062736"
    nd_37062832 = nd_37062896.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=0.999651724272, oid="Node37062832")
    nd_37062832.edge.oid = "Edge37062864"
    nd_37063024 = nd_37062896.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=0.999651724272, oid="Node37063024")
    nd_37063024.edge.oid = "Edge37062960"
    nd_37063152 = nd_37063184.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=6.6158514927, oid="Node37063152")
    nd_37063152.edge.oid = "Edge37062800"
    nd_37063280 = nd_37063184.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=6.6158514927, oid="Node37063280")
    nd_37063280.edge.oid = "Edge37063120"
    nd_37063248 = nd_37063088.new_child(label="Node37063248", taxon=None, edge_length=3.93324105164, oid="Node37063248")
    nd_37063248.edge.oid = "Edge37063344"
    nd_37063760 = nd_37063088.new_child(label="Node37063760", taxon=None, edge_length=8.5794579294, oid="Node37063760")
    nd_37063760.edge.oid = "Edge37063888"
    nd_37063312 = nd_37063248.new_child(label="Python sebae", taxon=tax_36961776, edge_length=14.3734153538, oid="Node37063312")
    nd_37063312.edge.oid = "Edge37063440"
    nd_37063536 = nd_37063248.new_child(label="Node37063536", taxon=None, edge_length=2.6097491119, oid="Node37063536")
    nd_37063536.edge.oid = "Edge37063472"
    nd_37063408 = nd_37063536.new_child(label="Python regius", taxon=tax_36961712, edge_length=11.7636662419, oid="Node37063408")
    nd_37063408.edge.oid = "Edge37063504"
    nd_37063664 = nd_37063536.new_child(label="Node37063664", taxon=None, edge_length=1.23184615824, oid="Node37063664")
    nd_37063664.edge.oid = "Edge37063600"
    nd_37063216 = nd_37063664.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=10.5318200836, oid="Node37063216")
    nd_37063216.edge.oid = "Edge37063632"
    nd_37063792 = nd_37063664.new_child(label="Python molurus", taxon=tax_36961648, edge_length=10.5318200836, oid="Node37063792")
    nd_37063792.edge.oid = "Edge37063728"
    nd_37063568 = nd_37063760.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=9.72719847603, oid="Node37063568")
    nd_37063568.edge.oid = "Edge37063856"
    nd_37063952 = nd_37063760.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=9.72719847603, oid="Node37063952")
    nd_37063952.edge.oid = "Edge37064016"
    tree_37063920 = dendropy.Tree(label="Tree03", taxon_set=tree_list.taxon_set, oid="Tree37063920")
    tree_list.append(tree_37063920, reindex_taxa=False)
    tree_37063920.seed_node.oid = 'Node37064144'
    nd_37064208 = tree_37063920.seed_node.new_child(label="Node37064208", taxon=None, edge_length=0.0364728670888, oid="Node37064208")
    nd_37064208.edge.oid = "Edge37064240"
    nd_37656048 = tree_37063920.seed_node.new_child(label="Node37656048", taxon=None, edge_length=1.25568695854, oid="Node37656048")
    nd_37656048.edge.oid = "Edge37656016"
    nd_36998896 = nd_37064208.new_child(label="Node36998896", taxon=None, edge_length=1.06929723627, oid="Node36998896")
    nd_36998896.edge.oid = "Edge37064272"
    nd_37655728 = nd_37064208.new_child(label="Node37655728", taxon=None, edge_length=10.8687376261, oid="Node37655728")
    nd_37655728.edge.oid = "Edge37655344"
    nd_37063696 = nd_36998896.new_child(label="Node37063696", taxon=None, edge_length=3.261105684, oid="Node37063696")
    nd_37063696.edge.oid = "Edge37064336"
    nd_37654704 = nd_36998896.new_child(label="Node37654704", taxon=None, edge_length=5.48562103095, oid="Node37654704")
    nd_37654704.edge.oid = "Edge37654736"
    nd_37064048 = nd_37063696.new_child(label="Node37064048", taxon=None, edge_length=2.97395412534, oid="Node37064048")
    nd_37064048.edge.oid = "Edge37064400"
    nd_37654640 = nd_37063696.new_child(label="Node37654640", taxon=None, edge_length=12.0688826259, oid="Node37654640")
    nd_37654640.edge.oid = "Edge37654608"
    nd_37064304 = nd_37064048.new_child(label="Python sebae", taxon=tax_36961776, edge_length=16.9224723285, oid="Node37064304")
    nd_37064304.edge.oid = "Edge37064464"
    nd_37064496 = nd_37064048.new_child(label="Node37064496", taxon=None, edge_length=2.75174107457, oid="Node37064496")
    nd_37064496.edge.oid = "Edge37063824"
    nd_37064432 = nd_37064496.new_child(label="Node37064432", taxon=None, edge_length=0.218628811139, oid="Node37064432")
    nd_37064432.edge.oid = "Edge37063984"
    nd_37064624 = nd_37064496.new_child(label="Python regius", taxon=tax_36961712, edge_length=14.170731254, oid="Node37064624")
    nd_37064624.edge.oid = "Edge37064368"
    nd_37064528 = nd_37064432.new_child(label="Python molurus", taxon=tax_36961648, edge_length=13.9521024428, oid="Node37064528")
    nd_37064528.edge.oid = "Edge37064592"
    nd_37064688 = nd_37064432.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=13.9521024428, oid="Node37064688")
    nd_37064688.edge.oid = "Edge37064560"
    nd_37654672 = nd_37654640.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=7.82754382801, oid="Node37654672")
    nd_37654672.edge.oid = "Edge37654576"
    nd_37654768 = nd_37654640.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=7.82754382801, oid="Node37654768")
    nd_37654768.edge.oid = "Edge37654832"
    nd_37654896 = nd_37654704.new_child(label="Node37654896", taxon=None, edge_length=2.9034791189, oid="Node37654896")
    nd_37654896.edge.oid = "Edge37654800"
    nd_37655184 = nd_37654704.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=17.6719111069, oid="Node37655184")
    nd_37655184.edge.oid = "Edge37655760"
    nd_37064656 = nd_37654896.new_child(label="Node37064656", taxon=None, edge_length=2.97877570893, oid="Node37064656")
    nd_37064656.edge.oid = "Edge37654960"
    nd_37655152 = nd_37654896.new_child(label="Node37655152", taxon=None, edge_length=1.25780607691, oid="Node37655152")
    nd_37655152.edge.oid = "Edge37655312"
    nd_37654864 = nd_37064656.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=11.7896562791, oid="Node37654864")
    nd_37654864.edge.oid = "Edge37655024"
    nd_37655120 = nd_37064656.new_child(label="Node37655120", taxon=None, edge_length=7.73520362135, oid="Node37655120")
    nd_37655120.edge.oid = "Edge37655056"
    nd_37654992 = nd_37655120.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=4.05445265776, oid="Node37654992")
    nd_37654992.edge.oid = "Edge37655088"
    nd_37655248 = nd_37655120.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=4.05445265776, oid="Node37655248")
    nd_37655248.edge.oid = "Edge37654928"
    nd_37655280 = nd_37655152.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=13.5106259111, oid="Node37655280")
    nd_37655280.edge.oid = "Edge37655216"
    nd_37655408 = nd_37655152.new_child(label="Node37655408", taxon=None, edge_length=2.27459934102, oid="Node37655408")
    nd_37655408.edge.oid = "Edge37655440"
    nd_37655472 = nd_37655408.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=11.2360265701, oid="Node37655472")
    nd_37655472.edge.oid = "Edge37655376"
    nd_37655536 = nd_37655408.new_child(label="Node37655536", taxon=None, edge_length=7.74754730265, oid="Node37655536")
    nd_37655536.edge.oid = "Edge37655568"
    nd_37655600 = nd_37655536.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=3.48847926746, oid="Node37655600")
    nd_37655600.edge.oid = "Edge37655504"
    nd_37655664 = nd_37655536.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=3.48847926746, oid="Node37655664")
    nd_37655664.edge.oid = "Edge37655696"
    nd_37655824 = nd_37655728.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=13.3580917481, oid="Node37655824")
    nd_37655824.edge.oid = "Edge37655792"
    nd_37655920 = nd_37655728.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=13.3580917481, oid="Node37655920")
    nd_37655920.edge.oid = "Edge37655952"
    nd_37655984 = nd_37656048.new_child(label="Node37655984", taxon=None, edge_length=14.878020141, oid="Node37655984")
    nd_37655984.edge.oid = "Edge37655856"
    nd_37656560 = nd_37656048.new_child(label="Node37656560", taxon=None, edge_length=2.53457554003, oid="Node37656560")
    nd_37656560.edge.oid = "Edge37656240"
    nd_37655632 = nd_37655984.new_child(label="Node37655632", taxon=None, edge_length=1.91471562327, oid="Node37655632")
    nd_37655632.edge.oid = "Edge37656112"
    nd_37656400 = nd_37655984.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=8.12959514173, oid="Node37656400")
    nd_37656400.edge.oid = "Edge37656464"
    nd_37655888 = nd_37655632.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=6.21487951846, oid="Node37655888")
    nd_37655888.edge.oid = "Edge37656176"
    nd_37656272 = nd_37655632.new_child(label="Node37656272", taxon=None, edge_length=3.60951843288, oid="Node37656272")
    nd_37656272.edge.oid = "Edge37656144"
    nd_37656080 = nd_37656272.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=2.60536108558, oid="Node37656080")
    nd_37656080.edge.oid = "Edge37656208"
    nd_37656368 = nd_37656272.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=2.60536108558, oid="Node37656368")
    nd_37656368.edge.oid = "Edge37656432"
    nd_37656336 = nd_37656560.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=20.4730397427, oid="Node37656336")
    nd_37656336.edge.oid = "Edge37656496"
    nd_37656656 = nd_37656560.new_child(label="Node37656656", taxon=None, edge_length=2.57209904849, oid="Node37656656")
    nd_37656656.edge.oid = "Edge37656592"
    nd_37656528 = nd_37656656.new_child(label="Node37656528", taxon=None, edge_length=0.296424924481, oid="Node37656528")
    nd_37656528.edge.oid = "Edge37656624"
    nd_37657136 = nd_37656656.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=17.9009406942, oid="Node37657136")
    nd_37657136.edge.oid = "Edge37657552"
    nd_37656688 = nd_37656528.new_child(label="Node37656688", taxon=None, edge_length=9.88501788667, oid="Node37656688")
    nd_37656688.edge.oid = "Edge37656752"
    nd_37657392 = nd_37656528.new_child(label="Node37657392", taxon=None, edge_length=11.8234917746, oid="Node37657392")
    nd_37657392.edge.oid = "Edge37657328"
    nd_37656304 = nd_37656688.new_child(label="Node37656304", taxon=None, edge_length=3.28625146667, oid="Node37656304")
    nd_37656304.edge.oid = "Edge37656816"
    nd_37657200 = nd_37656688.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=7.71949788304, oid="Node37657200")
    nd_37657200.edge.oid = "Edge37657040"
    nd_37656720 = nd_37656304.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=4.43324641637, oid="Node37656720")
    nd_37656720.edge.oid = "Edge37656880"
    nd_37656976 = nd_37656304.new_child(label="Node37656976", taxon=None, edge_length=2.99969759558, oid="Node37656976")
    nd_37656976.edge.oid = "Edge37656848"
    nd_37656912 = nd_37656976.new_child(label="Node37656912", taxon=None, edge_length=0.622029416262, oid="Node37656912")
    nd_37656912.edge.oid = "Edge37656944"
    nd_37657264 = nd_37656976.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=1.4335488208, oid="Node37657264")
    nd_37657264.edge.oid = "Edge37657232"
    nd_37657008 = nd_37656912.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=0.811519404534, oid="Node37657008")
    nd_37657008.edge.oid = "Edge37657072"
    nd_37657168 = nd_37656912.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=0.811519404534, oid="Node37657168")
    nd_37657168.edge.oid = "Edge37657104"
    nd_37657296 = nd_37657392.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=5.78102399509, oid="Node37657296")
    nd_37657296.edge.oid = "Edge37657360"
    nd_37657488 = nd_37657392.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=5.78102399509, oid="Node37657488")
    nd_37657488.edge.oid = "Edge37656784"
    tree_37657456 = dendropy.Tree(label="Tree04", taxon_set=tree_list.taxon_set, oid="Tree37657456")
    tree_list.append(tree_37657456, reindex_taxa=False)
    tree_37657456.seed_node.oid = 'Node37657712'
    nd_37657776 = tree_37657456.seed_node.new_child(label="Node37657776", taxon=None, edge_length=3.09272413642, oid="Node37657776")
    nd_37657776.edge.oid = "Edge37657808"
    nd_37658384 = tree_37657456.seed_node.new_child(label="Node37658384", taxon=None, edge_length=0.465997253751, oid="Node37658384")
    nd_37658384.edge.oid = "Edge37658448"
    nd_36962032 = nd_37657776.new_child(label="Node36962032", taxon=None, edge_length=3.43039474792, oid="Node36962032")
    nd_36962032.edge.oid = "Edge37657840"
    nd_37658224 = nd_37657776.new_child(label="Node37658224", taxon=None, edge_length=8.31909563528, oid="Node37658224")
    nd_37658224.edge.oid = "Edge37658352"
    nd_37657584 = nd_36962032.new_child(label="Python sebae", taxon=tax_36961776, edge_length=16.5921325784, oid="Node37657584")
    nd_37657584.edge.oid = "Edge37657904"
    nd_37658000 = nd_36962032.new_child(label="Node37658000", taxon=None, edge_length=2.21219736911, oid="Node37658000")
    nd_37658000.edge.oid = "Edge37657936"
    nd_37657872 = nd_37658000.new_child(label="Python regius", taxon=tax_36961712, edge_length=14.3799352093, oid="Node37657872")
    nd_37657872.edge.oid = "Edge37657968"
    nd_37658128 = nd_37658000.new_child(label="Node37658128", taxon=None, edge_length=0.997806328959, oid="Node37658128")
    nd_37658128.edge.oid = "Edge37658064"
    nd_37657616 = nd_37658128.new_child(label="Python molurus", taxon=tax_36961648, edge_length=13.3821288803, oid="Node37657616")
    nd_37657616.edge.oid = "Edge37658096"
    nd_37658256 = nd_37658128.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=13.3821288803, oid="Node37658256")
    nd_37658256.edge.oid = "Edge37658032"
    nd_37658192 = nd_37658224.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=11.703431691, oid="Node37658192")
    nd_37658192.edge.oid = "Edge37658320"
    nd_37658416 = nd_37658224.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=11.703431691, oid="Node37658416")
    nd_37658416.edge.oid = "Edge37658480"
    nd_37657648 = nd_37658384.new_child(label="Node37657648", taxon=None, edge_length=2.71760756257, oid="Node37657648")
    nd_37657648.edge.oid = "Edge37657520"
    nd_37712912 = nd_37658384.new_child(label="Node37712912", taxon=None, edge_length=0.725019015485, oid="Node37712912")
    nd_37712912.edge.oid = "Edge37712880"
    nd_37658288 = nd_37657648.new_child(label="Node37658288", taxon=None, edge_length=0.121667790832, oid="Node37658288")
    nd_37658288.edge.oid = "Edge37658544"
    nd_37712752 = nd_37657648.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=19.9316466464, oid="Node37712752")
    nd_37712752.edge.oid = "Edge37712688"
    nd_37658160 = nd_37658288.new_child(label="Node37658160", taxon=None, edge_length=13.3467104377, oid="Node37658160")
    nd_37658160.edge.oid = "Edge37658608"
    nd_37658576 = nd_37658288.new_child(label="Node37658576", taxon=None, edge_length=10.3912007878, oid="Node37658576")
    nd_37658576.edge.oid = "Edge37712144"
    nd_37658512 = nd_37658160.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=6.46326841791, oid="Node37658512")
    nd_37658512.edge.oid = "Edge37711952"
    nd_37712048 = nd_37658160.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=6.46326841791, oid="Node37712048")
    nd_37712048.edge.oid = "Edge37711920"
    nd_37712016 = nd_37658576.new_child(label="Node37712016", taxon=None, edge_length=3.03031595643, oid="Node37712016")
    nd_37712016.edge.oid = "Edge37712080"
    nd_37712784 = nd_37658576.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=9.41877806779, oid="Node37712784")
    nd_37712784.edge.oid = "Edge37712528"
    nd_37711984 = nd_37712016.new_child(label="Node37711984", taxon=None, edge_length=1.71754130617, oid="Node37711984")
    nd_37711984.edge.oid = "Edge37712208"
    nd_37712592 = nd_37712016.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=6.38846211136, oid="Node37712592")
    nd_37712592.edge.oid = "Edge37712720"
    nd_37712112 = nd_37711984.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=4.67092080519, oid="Node37712112")
    nd_37712112.edge.oid = "Edge37712272"
    nd_37712368 = nd_37711984.new_child(label="Node37712368", taxon=None, edge_length=2.58211994653, oid="Node37712368")
    nd_37712368.edge.oid = "Edge37712240"
    nd_37712304 = nd_37712368.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=2.08880085866, oid="Node37712304")
    nd_37712304.edge.oid = "Edge37712336"
    nd_37712496 = nd_37712368.new_child(label="Node37712496", taxon=None, edge_length=0.607570395039, oid="Node37712496")
    nd_37712496.edge.oid = "Edge37712432"
    nd_37712176 = nd_37712496.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=1.48123046362, oid="Node37712176")
    nd_37712176.edge.oid = "Edge37712464"
    nd_37712624 = nd_37712496.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=1.48123046362, oid="Node37712624")
    nd_37712624.edge.oid = "Edge37712560"
    nd_37712848 = nd_37712912.new_child(label="Node37712848", taxon=None, edge_length=1.89223670238, oid="Node37712848")
    nd_37712848.edge.oid = "Edge37712656"
    nd_37714192 = nd_37712912.new_child(label="Node37714192", taxon=None, edge_length=6.6830638124, oid="Node37714192")
    nd_37714192.edge.oid = "Edge37714160"
    nd_37712816 = nd_37712848.new_child(label="Node37712816", taxon=None, edge_length=1.84844364601, oid="Node37712816")
    nd_37712816.edge.oid = "Edge37712976"
    nd_37714000 = nd_37712848.new_child(label="Node37714000", taxon=None, edge_length=7.28635782162, oid="Node37714000")
    nd_37714000.edge.oid = "Edge37713904"
    nd_37712400 = nd_37712816.new_child(label="Node37712400", taxon=None, edge_length=2.23472642535, oid="Node37712400")
    nd_37712400.edge.oid = "Edge37713040"
    nd_37713808 = nd_37712816.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=18.1835548451, oid="Node37713808")
    nd_37713808.edge.oid = "Edge37713776"
    nd_37712944 = nd_37712400.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=15.9488284198, oid="Node37712944")
    nd_37712944.edge.oid = "Edge37713104"
    nd_37713200 = nd_37712400.new_child(label="Node37713200", taxon=None, edge_length=0.394642188595, oid="Node37713200")
    nd_37713200.edge.oid = "Edge37713136"
    nd_37713072 = nd_37713200.new_child(label="Node37713072", taxon=None, edge_length=1.55388791109, oid="Node37713072")
    nd_37713072.edge.oid = "Edge37713168"
    nd_37713648 = nd_37713200.new_child(label="Node37713648", taxon=None, edge_length=11.0588698605, oid="Node37713648")
    nd_37713648.edge.oid = "Edge37713712"
    nd_37713232 = nd_37713072.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=14.0002983201, oid="Node37713232")
    nd_37713232.edge.oid = "Edge37713296"
    nd_37713360 = nd_37713072.new_child(label="Node37713360", taxon=None, edge_length=2.5525380966, oid="Node37713360")
    nd_37713360.edge.oid = "Edge37713424"
    nd_37713392 = nd_37713360.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=11.4477602235, oid="Node37713392")
    nd_37713392.edge.oid = "Edge37713328"
    nd_37713488 = nd_37713360.new_child(label="Node37713488", taxon=None, edge_length=5.71931728003, oid="Node37713488")
    nd_37713488.edge.oid = "Edge37713520"
    nd_37713552 = nd_37713488.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=5.72844294344, oid="Node37713552")
    nd_37713552.edge.oid = "Edge37713456"
    nd_37713616 = nd_37713488.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=5.72844294344, oid="Node37713616")
    nd_37713616.edge.oid = "Edge37713008"
    nd_37713584 = nd_37713648.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=4.49531637063, oid="Node37713584")
    nd_37713584.edge.oid = "Edge37713744"
    nd_37713840 = nd_37713648.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=4.49531637063, oid="Node37713840")
    nd_37713840.edge.oid = "Edge37713264"
    nd_37713872 = nd_37714000.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=12.7456406695, oid="Node37713872")
    nd_37713872.edge.oid = "Edge37713680"
    nd_37714064 = nd_37714000.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=12.7456406695, oid="Node37714064")
    nd_37714064.edge.oid = "Edge37713936"
    nd_37714128 = nd_37714192.new_child(label="Node37714128", taxon=None, edge_length=1.71076794287, oid="Node37714128")
    nd_37714128.edge.oid = "Edge37714096"
    nd_37714320 = nd_37714192.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=15.2411713811, oid="Node37714320")
    nd_37714320.edge.oid = "Edge37714544"
    nd_37714032 = nd_37714128.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=13.5304034382, oid="Node37714032")
    nd_37714032.edge.oid = "Edge37714256"
    nd_37714352 = nd_37714128.new_child(label="Node37714352", taxon=None, edge_length=11.7858576376, oid="Node37714352")
    nd_37714352.edge.oid = "Edge37714224"
    nd_37713968 = nd_37714352.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=1.74454580066, oid="Node37713968")
    nd_37713968.edge.oid = "Edge37714288"
    nd_37714448 = nd_37714352.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=1.74454580066, oid="Node37714448")
    nd_37714448.edge.oid = "Edge37714480"
    tree_37714608 = dendropy.Tree(label="Tree05", taxon_set=tree_list.taxon_set, oid="Tree37714608")
    tree_list.append(tree_37714608, reindex_taxa=False)
    tree_37714608.seed_node.oid = 'Node37714704'
    nd_37714768 = tree_37714608.seed_node.new_child(label="Node37714768", taxon=None, edge_length=17.6343232815, oid="Node37714768")
    nd_37714768.edge.oid = "Edge37714800"
    nd_37715024 = tree_37714608.seed_node.new_child(label="Node37715024", taxon=None, edge_length=3.72471742653, oid="Node37715024")
    nd_37715024.edge.oid = "Edge37714896"
    nd_36961840 = nd_37714768.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=10.8300265774, oid="Node36961840")
    nd_36961840.edge.oid = "Edge37714832"
    nd_37714928 = nd_37714768.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=10.8300265774, oid="Node37714928")
    nd_37714928.edge.oid = "Edge37714864"
    nd_37714512 = nd_37715024.new_child(label="Node37714512", taxon=None, edge_length=4.52906062503, oid="Node37714512")
    nd_37714512.edge.oid = "Edge37714576"
    nd_37773872 = nd_37715024.new_child(label="Node37773872", taxon=None, edge_length=0.57539806967, oid="Node37773872")
    nd_37773872.edge.oid = "Edge37774064"
    nd_37714960 = nd_37714512.new_child(label="Node37714960", taxon=None, edge_length=1.31772152158, oid="Node37714960")
    nd_37714960.edge.oid = "Edge37715088"
    nd_37773680 = nd_37714512.new_child(label="Node37773680", taxon=None, edge_length=5.44916502448, oid="Node37773680")
    nd_37773680.edge.oid = "Edge37773648"
    nd_37714992 = nd_37714960.new_child(label="Node37714992", taxon=None, edge_length=4.88508122224, oid="Node37714992")
    nd_37714992.edge.oid = "Edge37715152"
    nd_37715600 = nd_37714960.new_child(label="Node37715600", taxon=None, edge_length=7.74780761159, oid="Node37715600")
    nd_37715600.edge.oid = "Edge37773360"
    nd_37715056 = nd_37714992.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=14.0077690635, oid="Node37715056")
    nd_37715056.edge.oid = "Edge37715216"
    nd_37715312 = nd_37714992.new_child(label="Node37715312", taxon=None, edge_length=0.404170030911, oid="Node37715312")
    nd_37715312.edge.oid = "Edge37715248"
    nd_37715184 = nd_37715312.new_child(label="Node37715184", taxon=None, edge_length=1.82500345118, oid="Node37715184")
    nd_37715184.edge.oid = "Edge37715280"
    nd_37715536 = nd_37715312.new_child(label="Node37715536", taxon=None, edge_length=1.27612254021, oid="Node37715536")
    nd_37715536.edge.oid = "Edge37715728"
    nd_37715344 = nd_37715184.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=11.7785955814, oid="Node37715344")
    nd_37715344.edge.oid = "Edge37715408"
    nd_37715504 = nd_37715184.new_child(label="Node37715504", taxon=None, edge_length=6.80279787983, oid="Node37715504")
    nd_37715504.edge.oid = "Edge37715440"
    nd_37715376 = nd_37715504.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=4.97579770159, oid="Node37715376")
    nd_37715376.edge.oid = "Edge37715472"
    nd_37715632 = nd_37715504.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=4.97579770159, oid="Node37715632")
    nd_37715632.edge.oid = "Edge37715120"
    nd_37715568 = nd_37715536.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=12.3274764924, oid="Node37715568")
    nd_37715568.edge.oid = "Edge37715664"
    nd_37715792 = nd_37715536.new_child(label="Node37715792", taxon=None, edge_length=4.32415504222, oid="Node37715792")
    nd_37715792.edge.oid = "Edge37715824"
    nd_37715856 = nd_37715792.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=8.00332145017, oid="Node37715856")
    nd_37715856.edge.oid = "Edge37714416"
    nd_37715760 = nd_37715792.new_child(label="Node37715760", taxon=None, edge_length=3.98809707211, oid="Node37715760")
    nd_37715760.edge.oid = "Edge37715920"
    nd_37715888 = nd_37715760.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=4.01522437806, oid="Node37715888")
    nd_37715888.edge.oid = "Edge37715696"
    nd_37714640 = nd_37715760.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=4.01522437806, oid="Node37714640")
    nd_37714640.edge.oid = "Edge37715952"
    nd_37773456 = nd_37715600.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=11.1450426742, oid="Node37773456")
    nd_37773456.edge.oid = "Edge37773392"
    nd_37773552 = nd_37715600.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=11.1450426742, oid="Node37773552")
    nd_37773552.edge.oid = "Edge37773584"
    nd_37773616 = nd_37773680.new_child(label="Node37773616", taxon=None, edge_length=3.78712208743, oid="Node37773616")
    nd_37773616.edge.oid = "Edge37773488"
    nd_37773968 = nd_37773680.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=14.7614067829, oid="Node37773968")
    nd_37773968.edge.oid = "Edge37773808"
    nd_37773520 = nd_37773616.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=10.9742846954, oid="Node37773520")
    nd_37773520.edge.oid = "Edge37773744"
    nd_37773840 = nd_37773616.new_child(label="Node37773840", taxon=None, edge_length=7.80758505868, oid="Node37773840")
    nd_37773840.edge.oid = "Edge37773424"
    nd_37773712 = nd_37773840.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=3.16669963676, oid="Node37773712")
    nd_37773712.edge.oid = "Edge37773776"
    nd_37773936 = nd_37773840.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=3.16669963676, oid="Node37773936")
    nd_37773936.edge.oid = "Edge37774000"
    nd_37774128 = nd_37773872.new_child(label="Node37774128", taxon=None, edge_length=1.8069504389, oid="Node37774128")
    nd_37774128.edge.oid = "Edge37774096"
    nd_37775024 = nd_37773872.new_child(label="Node37775024", taxon=None, edge_length=3.10812197386, oid="Node37775024")
    nd_37775024.edge.oid = "Edge37774992"
    nd_37774032 = nd_37774128.new_child(label="Node37774032", taxon=None, edge_length=5.48583318771, oid="Node37774032")
    nd_37774032.edge.oid = "Edge37774192"
    nd_37774960 = nd_37774128.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=22.3572839238, oid="Node37774960")
    nd_37774960.edge.oid = "Edge37774832"
    nd_37773904 = nd_37774032.new_child(label="Node37773904", taxon=None, edge_length=6.85886795982, oid="Node37773904")
    nd_37773904.edge.oid = "Edge37774256"
    nd_37774896 = nd_37774032.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=16.8714507361, oid="Node37774896")
    nd_37774896.edge.oid = "Edge37774640"
    nd_37774160 = nd_37773904.new_child(label="Node37774160", taxon=None, edge_length=4.50938662462, oid="Node37774160")
    nd_37774160.edge.oid = "Edge37774320"
    nd_37774800 = nd_37773904.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=10.0125827763, oid="Node37774800")
    nd_37774800.edge.oid = "Edge37774704"
    nd_37774224 = nd_37774160.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=5.50319615164, oid="Node37774224")
    nd_37774224.edge.oid = "Edge37774384"
    nd_37774480 = nd_37774160.new_child(label="Node37774480", taxon=None, edge_length=3.93061957053, oid="Node37774480")
    nd_37774480.edge.oid = "Edge37774352"
    nd_37774416 = nd_37774480.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=1.57257658111, oid="Node37774416")
    nd_37774416.edge.oid = "Edge37774448"
    nd_37774608 = nd_37774480.new_child(label="Node37774608", taxon=None, edge_length=0.136848033805, oid="Node37774608")
    nd_37774608.edge.oid = "Edge37774544"
    nd_37774288 = nd_37774608.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=1.4357285473, oid="Node37774288")
    nd_37774288.edge.oid = "Edge37774576"
    nd_37774736 = nd_37774608.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=1.4357285473, oid="Node37774736")
    nd_37774736.edge.oid = "Edge37774512"
    nd_37774928 = nd_37775024.new_child(label="Node37774928", taxon=None, edge_length=4.52064732483, oid="Node37774928")
    nd_37774928.edge.oid = "Edge37774672"
    nd_37775408 = nd_37775024.new_child(label="Node37775408", taxon=None, edge_length=6.39210928917, oid="Node37775408")
    nd_37775408.edge.oid = "Edge37775216"
    nd_37774768 = nd_37774928.new_child(label="Python sebae", taxon=tax_36961776, edge_length=16.535465064, oid="Node37774768")
    nd_37774768.edge.oid = "Edge37775088"
    nd_37775184 = nd_37774928.new_child(label="Node37775184", taxon=None, edge_length=1.573045967, oid="Node37775184")
    nd_37775184.edge.oid = "Edge37775120"
    nd_37775056 = nd_37775184.new_child(label="Python regius", taxon=tax_36961712, edge_length=14.962419097, oid="Node37775056")
    nd_37775056.edge.oid = "Edge37775152"
    nd_37775312 = nd_37775184.new_child(label="Node37775312", taxon=None, edge_length=1.82698701194, oid="Node37775312")
    nd_37775312.edge.oid = "Edge37775248"
    nd_37774864 = nd_37775312.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=13.1354320851, oid="Node37774864")
    nd_37774864.edge.oid = "Edge37775280"
    nd_37775440 = nd_37775312.new_child(label="Python molurus", taxon=tax_36961648, edge_length=13.1354320851, oid="Node37775440")
    nd_37775440.edge.oid = "Edge37775376"
    nd_37775536 = nd_37775408.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=14.6640030997, oid="Node37775536")
    nd_37775536.edge.oid = "Edge37775504"
    nd_37775600 = nd_37775408.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=14.6640030997, oid="Node37775600")
    nd_37775600.edge.oid = "Edge37775344"
    tree_37775728 = dendropy.Tree(label="Tree06", taxon_set=tree_list.taxon_set, oid="Tree37775728")
    tree_list.append(tree_37775728, reindex_taxa=False)
    tree_37775728.seed_node.oid = 'Node37775792'
    nd_37775856 = tree_37775728.seed_node.new_child(label="Node37775856", taxon=None, edge_length=1.01117536037, oid="Node37775856")
    nd_37775856.edge.oid = "Edge37775888"
    nd_37831600 = tree_37775728.seed_node.new_child(label="Node37831600", taxon=None, edge_length=0.075772291124, oid="Node37831600")
    nd_37831600.edge.oid = "Edge37831728"
    nd_37064112 = nd_37775856.new_child(label="Node37064112", taxon=None, edge_length=4.09971291135, oid="Node37064112")
    nd_37064112.edge.oid = "Edge37775920"
    nd_37830864 = nd_37775856.new_child(label="Node37830864", taxon=None, edge_length=3.12950470552, oid="Node37830864")
    nd_37830864.edge.oid = "Edge37830832"
    nd_37775568 = nd_37064112.new_child(label="Node37775568", taxon=None, edge_length=0.754225130388, oid="Node37775568")
    nd_37775568.edge.oid = "Edge37775984"
    nd_37777264 = nd_37064112.new_child(label="Node37777264", taxon=None, edge_length=2.9800994962, oid="Node37777264")
    nd_37777264.edge.oid = "Edge37777104"
    nd_37775472 = nd_37775568.new_child(label="Node37775472", taxon=None, edge_length=1.65141371467, oid="Node37775472")
    nd_37775472.edge.oid = "Edge37776048"
    nd_37777072 = nd_37775568.new_child(label="Node37777072", taxon=None, edge_length=4.17813417278, oid="Node37777072")
    nd_37777072.edge.oid = "Edge37776880"
    nd_37775952 = nd_37775472.new_child(label="Node37775952", taxon=None, edge_length=0.835745557531, oid="Node37775952")
    nd_37775952.edge.oid = "Edge37776112"
    nd_37776976 = nd_37775472.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=15.2896785188, oid="Node37776976")
    nd_37776976.edge.oid = "Edge37777008"
    nd_37776016 = nd_37775952.new_child(label="Node37776016", taxon=None, edge_length=1.83427060068, oid="Node37776016")
    nd_37776016.edge.oid = "Edge37776176"
    nd_37776144 = nd_37775952.new_child(label="Node37776144", taxon=None, edge_length=1.03775314121, oid="Node37776144")
    nd_37776144.edge.oid = "Edge37776208"
    nd_37776080 = nd_37776016.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=12.6196623606, oid="Node37776080")
    nd_37776080.edge.oid = "Edge37776240"
    nd_37776304 = nd_37776016.new_child(label="Node37776304", taxon=None, edge_length=2.51279801319, oid="Node37776304")
    nd_37776304.edge.oid = "Edge37776336"
    nd_37776368 = nd_37776304.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=10.1068643474, oid="Node37776368")
    nd_37776368.edge.oid = "Edge37776272"
    nd_37776432 = nd_37776304.new_child(label="Node37776432", taxon=None, edge_length=5.71367528366, oid="Node37776432")
    nd_37776432.edge.oid = "Edge37776464"
    nd_37776496 = nd_37776432.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=4.39318906371, oid="Node37776496")
    nd_37776496.edge.oid = "Edge37776400"
    nd_37776560 = nd_37776432.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=4.39318906371, oid="Node37776560")
    nd_37776560.edge.oid = "Edge37776624"
    nd_37776656 = nd_37776144.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=13.41617982, oid="Node37776656")
    nd_37776656.edge.oid = "Edge37776528"
    nd_37776784 = nd_37776144.new_child(label="Node37776784", taxon=None, edge_length=9.37635268364, oid="Node37776784")
    nd_37776784.edge.oid = "Edge37776720"
    nd_37776592 = nd_37776784.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=4.03982713639, oid="Node37776592")
    nd_37776592.edge.oid = "Edge37776752"
    nd_37776912 = nd_37776784.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=4.03982713639, oid="Node37776912")
    nd_37776912.edge.oid = "Edge37776688"
    nd_37776944 = nd_37777072.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=12.7629580607, oid="Node37776944")
    nd_37776944.edge.oid = "Edge37776816"
    nd_37777136 = nd_37777072.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=12.7629580607, oid="Node37777136")
    nd_37777136.edge.oid = "Edge37777168"
    nd_37777232 = nd_37777264.new_child(label="Node37777232", taxon=None, edge_length=4.06497756637, oid="Node37777232")
    nd_37777232.edge.oid = "Edge37776848"
    nd_37830768 = nd_37777264.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=14.7152178676, oid="Node37830768")
    nd_37830768.edge.oid = "Edge37830800"
    nd_37777200 = nd_37777232.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=10.6502403013, oid="Node37777200")
    nd_37777200.edge.oid = "Edge37775696"
    nd_37777360 = nd_37777232.new_child(label="Node37777360", taxon=None, edge_length=6.97858461332, oid="Node37777360")
    nd_37777360.edge.oid = "Edge37777296"
    nd_37775632 = nd_37777360.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=3.67165568794, oid="Node37775632")
    nd_37775632.edge.oid = "Edge37777328"
    nd_37777392 = nd_37777360.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=3.67165568794, oid="Node37777392")
    nd_37777392.edge.oid = "Edge37830736"
    nd_37830704 = nd_37830864.new_child(label="Node37830704", taxon=None, edge_length=6.58446507703, oid="Node37830704")
    nd_37830704.edge.oid = "Edge37830896"
    nd_37831696 = nd_37830864.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=18.6655255697, oid="Node37831696")
    nd_37831696.edge.oid = "Edge37831536"
    nd_37777040 = nd_37830704.new_child(label="Node37777040", taxon=None, edge_length=3.1657791763, oid="Node37777040")
    nd_37777040.edge.oid = "Edge37830992"
    nd_37831632 = nd_37830704.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=12.0810604926, oid="Node37831632")
    nd_37831632.edge.oid = "Edge37831376"
    nd_37830928 = nd_37777040.new_child(label="Node37830928", taxon=None, edge_length=4.27056899374, oid="Node37830928")
    nd_37830928.edge.oid = "Edge37831056"
    nd_37831440 = nd_37777040.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=8.91528131634, oid="Node37831440")
    nd_37831440.edge.oid = "Edge37831568"
    nd_37830960 = nd_37830928.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=4.6447123226, oid="Node37830960")
    nd_37830960.edge.oid = "Edge37831120"
    nd_37831216 = nd_37830928.new_child(label="Node37831216", taxon=None, edge_length=2.47296326498, oid="Node37831216")
    nd_37831216.edge.oid = "Edge37831088"
    nd_37831024 = nd_37831216.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=2.17174905761, oid="Node37831024")
    nd_37831024.edge.oid = "Edge37831152"
    nd_37831344 = nd_37831216.new_child(label="Node37831344", taxon=None, edge_length=1.26765325952, oid="Node37831344")
    nd_37831344.edge.oid = "Edge37831184"
    nd_37831280 = nd_37831344.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=0.904095798097, oid="Node37831280")
    nd_37831280.edge.oid = "Edge37831312"
    nd_37831472 = nd_37831344.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=0.904095798097, oid="Node37831472")
    nd_37831472.edge.oid = "Edge37831408"
    nd_37831504 = nd_37831600.new_child(label="Node37831504", taxon=None, edge_length=11.1718400909, oid="Node37831504")
    nd_37831504.edge.oid = "Edge37831248"
    nd_37832016 = nd_37831600.new_child(label="Node37832016", taxon=None, edge_length=1.44865397646, oid="Node37832016")
    nd_37832016.edge.oid = "Edge37831792"
    nd_37831664 = nd_37831504.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=11.5585932535, oid="Node37831664")
    nd_37831664.edge.oid = "Edge37831824"
    nd_37831920 = nd_37831504.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=11.5585932535, oid="Node37831920")
    nd_37831920.edge.oid = "Edge37831856"
    nd_37831888 = nd_37832016.new_child(label="Node37831888", taxon=None, edge_length=5.31414730047, oid="Node37831888")
    nd_37831888.edge.oid = "Edge37831760"
    nd_37832400 = nd_37832016.new_child(label="Node37832400", taxon=None, edge_length=5.70308191065, oid="Node37832400")
    nd_37832400.edge.oid = "Edge37832496"
    nd_37831952 = nd_37831888.new_child(label="Python sebae", taxon=tax_36961776, edge_length=15.9676320675, oid="Node37831952")
    nd_37831952.edge.oid = "Edge37832080"
    nd_37832176 = nd_37831888.new_child(label="Node37832176", taxon=None, edge_length=1.06278922707, oid="Node37832176")
    nd_37832176.edge.oid = "Edge37832144"
    nd_37832048 = nd_37832176.new_child(label="Node37832048", taxon=None, edge_length=1.04487219988, oid="Node37832048")
    nd_37832048.edge.oid = "Edge37831984"
    nd_37832464 = nd_37832176.new_child(label="Python regius", taxon=tax_36961712, edge_length=14.9048428404, oid="Node37832464")
    nd_37832464.edge.oid = "Edge37832240"
    nd_37832208 = nd_37832048.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=13.8599706405, oid="Node37832208")
    nd_37832208.edge.oid = "Edge37832272"
    nd_37832368 = nd_37832048.new_child(label="Python molurus", taxon=tax_36961648, edge_length=13.8599706405, oid="Node37832368")
    nd_37832368.edge.oid = "Edge37832304"
    nd_37832432 = nd_37832400.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=15.5786974573, oid="Node37832432")
    nd_37832432.edge.oid = "Edge37832112"
    nd_37832592 = nd_37832400.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=15.5786974573, oid="Node37832592")
    nd_37832592.edge.oid = "Edge37832336"
    tree_37832720 = dendropy.Tree(label="Tree07", taxon_set=tree_list.taxon_set, oid="Tree37832720")
    tree_list.append(tree_37832720, reindex_taxa=False)
    tree_37832720.seed_node.oid = 'Node37832784'
    nd_37832848 = tree_37832720.seed_node.new_child(label="Node37832848", taxon=None, edge_length=2.12261395029, oid="Node37832848")
    nd_37832848.edge.oid = "Edge37832880"
    nd_37833200 = tree_37832720.seed_node.new_child(label="Node37833200", taxon=None, edge_length=8.39682440052, oid="Node37833200")
    nd_37833200.edge.oid = "Edge37833232"
    nd_37657680 = nd_37832848.new_child(label="Python sebae", taxon=tax_36961776, edge_length=29.0711850036, oid="Node37657680")
    nd_37657680.edge.oid = "Edge37832912"
    nd_37833008 = nd_37832848.new_child(label="Node37833008", taxon=None, edge_length=6.65870230472, oid="Node37833008")
    nd_37833008.edge.oid = "Edge37832528"
    nd_37832944 = nd_37833008.new_child(label="Python molurus", taxon=tax_36961648, edge_length=22.4124826989, oid="Node37832944")
    nd_37832944.edge.oid = "Edge37832976"
    nd_37833136 = nd_37833008.new_child(label="Node37833136", taxon=None, edge_length=0.485999991488, oid="Node37833136")
    nd_37833136.edge.oid = "Edge37833072"
    nd_37832656 = nd_37833136.new_child(label="Python regius", taxon=tax_36961712, edge_length=21.9264827074, oid="Node37832656")
    nd_37832656.edge.oid = "Edge37833104"
    nd_37833264 = nd_37833136.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=21.9264827074, oid="Node37833264")
    nd_37833264.edge.oid = "Edge37833040"
    nd_37833360 = nd_37833200.new_child(label="Node37833360", taxon=None, edge_length=1.24633634472, oid="Node37833360")
    nd_37833360.edge.oid = "Edge37833296"
    nd_37103024 = nd_37833200.new_child(label="Node37103024", taxon=None, edge_length=12.9499290649, oid="Node37103024")
    nd_37103024.edge.oid = "Edge37103056"
    nd_37833328 = nd_37833360.new_child(label="Node37833328", taxon=None, edge_length=1.10888028588, oid="Node37833328")
    nd_37833328.edge.oid = "Edge37833424"
    nd_37101936 = nd_37833360.new_child(label="Node37101936", taxon=None, edge_length=3.94969364946, oid="Node37101936")
    nd_37101936.edge.oid = "Edge37102256"
    nd_37833168 = nd_37833328.new_child(label="Node37833168", taxon=None, edge_length=5.89251749724, oid="Node37833168")
    nd_37833168.edge.oid = "Edge37833488"
    nd_37833744 = nd_37833328.new_child(label="Node37833744", taxon=None, edge_length=1.52019767757, oid="Node37833744")
    nd_37833744.edge.oid = "Edge37833616"
    nd_37833392 = nd_37833168.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=14.5492404256, oid="Node37833392")
    nd_37833392.edge.oid = "Edge37833552"
    nd_37833648 = nd_37833168.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=14.5492404256, oid="Node37833648")
    nd_37833648.edge.oid = "Edge37833584"
    nd_37833456 = nd_37833744.new_child(label="Node37833456", taxon=None, edge_length=0.826870765744, oid="Node37833456")
    nd_37833456.edge.oid = "Edge37833712"
    nd_37101712 = nd_37833744.new_child(label="Node37101712", taxon=None, edge_length=5.2034050718, oid="Node37101712")
    nd_37101712.edge.oid = "Edge37101648"
    nd_37833680 = nd_37833456.new_child(label="Node37833680", taxon=None, edge_length=2.93802844014, oid="Node37833680")
    nd_37833680.edge.oid = "Edge37833808"
    nd_37834736 = nd_37833456.new_child(label="Node37834736", taxon=None, edge_length=8.81046360463, oid="Node37834736")
    nd_37834736.edge.oid = "Edge37832624"
    nd_37833520 = nd_37833680.new_child(label="Node37833520", taxon=None, edge_length=0.322649062335, oid="Node37833520")
    nd_37833520.edge.oid = "Edge37833872"
    nd_37834256 = nd_37833680.new_child(label="Node37834256", taxon=None, edge_length=0.821527298272, oid="Node37834256")
    nd_37834256.edge.oid = "Edge37834352"
    nd_37833776 = nd_37833520.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=14.834011977, oid="Node37833776")
    nd_37833776.edge.oid = "Edge37833936"
    nd_37834032 = nd_37833520.new_child(label="Node37834032", taxon=None, edge_length=1.22843868872, oid="Node37834032")
    nd_37834032.edge.oid = "Edge37833968"
    nd_37833904 = nd_37834032.new_child(label="Node37833904", taxon=None, edge_length=7.33173338162, oid="Node37833904")
    nd_37833904.edge.oid = "Edge37834000"
    nd_37834320 = nd_37834032.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=13.6055732883, oid="Node37834320")
    nd_37834320.edge.oid = "Edge37834096"
    nd_37834064 = nd_37833904.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=6.2738399067, oid="Node37834064")
    nd_37834064.edge.oid = "Edge37834128"
    nd_37834224 = nd_37833904.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=6.2738399067, oid="Node37834224")
    nd_37834224.edge.oid = "Edge37834160"
    nd_37834288 = nd_37834256.new_child(label="Node37834288", taxon=None, edge_length=4.22599220252, oid="Node37834288")
    nd_37834288.edge.oid = "Edge37833840"
    nd_37834608 = nd_37834256.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=14.3351337411, oid="Node37834608")
    nd_37834608.edge.oid = "Edge37832688"
    nd_37834192 = nd_37834288.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=10.1091415386, oid="Node37834192")
    nd_37834192.edge.oid = "Edge37834448"
    nd_37834512 = nd_37834288.new_child(label="Node37834512", taxon=None, edge_length=6.62594074811, oid="Node37834512")
    nd_37834512.edge.oid = "Edge37834544"
    nd_37834576 = nd_37834512.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=3.48320079048, oid="Node37834576")
    nd_37834576.edge.oid = "Edge37834480"
    nd_37834640 = nd_37834512.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=3.48320079048, oid="Node37834640")
    nd_37834640.edge.oid = "Edge37834672"
    nd_37834704 = nd_37834736.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=9.28422587488, oid="Node37834704")
    nd_37834704.edge.oid = "Edge37101616"
    nd_37834384 = nd_37834736.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=9.28422587488, oid="Node37834384")
    nd_37834384.edge.oid = "Edge37101680"
    nd_37101808 = nd_37101712.new_child(label="Node37101808", taxon=None, edge_length=2.70301025103, oid="Node37101808")
    nd_37101808.edge.oid = "Edge37101776"
    nd_37102000 = nd_37101712.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=13.7181551735, oid="Node37102000")
    nd_37102000.edge.oid = "Edge37102192"
    nd_37834416 = nd_37101808.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=11.0151449224, oid="Node37834416")
    nd_37834416.edge.oid = "Edge37101872"
    nd_37101968 = nd_37101808.new_child(label="Node37101968", taxon=None, edge_length=6.92713376313, oid="Node37101968")
    nd_37101968.edge.oid = "Edge37101744"
    nd_37101840 = nd_37101968.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=4.0880111593, oid="Node37101840")
    nd_37101840.edge.oid = "Edge37101904"
    nd_37102064 = nd_37101968.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=4.0880111593, oid="Node37102064")
    nd_37102064.edge.oid = "Edge37102128"
    nd_37102160 = nd_37101936.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=17.6009445592, oid="Node37102160")
    nd_37102160.edge.oid = "Edge37102224"
    nd_37102352 = nd_37101936.new_child(label="Node37102352", taxon=None, edge_length=5.4268859952, oid="Node37102352")
    nd_37102352.edge.oid = "Edge37102096"
    nd_37102288 = nd_37102352.new_child(label="Node37102288", taxon=None, edge_length=6.37126259081, oid="Node37102288")
    nd_37102288.edge.oid = "Edge37102320"
    nd_37102576 = nd_37102352.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=12.174058564, oid="Node37102576")
    nd_37102576.edge.oid = "Edge37102736"
    nd_37102384 = nd_37102288.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=5.80279597323, oid="Node37102384")
    nd_37102384.edge.oid = "Edge37102448"
    nd_37102544 = nd_37102288.new_child(label="Node37102544", taxon=None, edge_length=1.42172379881, oid="Node37102544")
    nd_37102544.edge.oid = "Edge37102480"
    nd_37102416 = nd_37102544.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=4.38107217442, oid="Node37102416")
    nd_37102416.edge.oid = "Edge37102512"
    nd_37102672 = nd_37102544.new_child(label="Node37102672", taxon=None, edge_length=2.50070781902, oid="Node37102672")
    nd_37102672.edge.oid = "Edge37102032"
    nd_37102608 = nd_37102672.new_child(label="Node37102608", taxon=None, edge_length=0.0774488112983, oid="Node37102608")
    nd_37102608.edge.oid = "Edge37102640"
    nd_37102960 = nd_37102672.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=1.8803643554, oid="Node37102960")
    nd_37102960.edge.oid = "Edge37102832"
    nd_37102704 = nd_37102608.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=1.8029155441, oid="Node37102704")
    nd_37102704.edge.oid = "Edge37102768"
    nd_37102864 = nd_37102608.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=1.8029155441, oid="Node37102864")
    nd_37102864.edge.oid = "Edge37102800"
    nd_37102992 = nd_37103024.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=9.84704548855, oid="Node37102992")
    nd_37102992.edge.oid = "Edge37102928"
    nd_37103152 = nd_37103024.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=9.84704548855, oid="Node37103152")
    nd_37103152.edge.oid = "Edge37103216"
    tree_37103088 = dendropy.Tree(label="Tree08", taxon_set=tree_list.taxon_set, oid="Tree37103088")
    tree_list.append(tree_37103088, reindex_taxa=False)
    tree_37103088.seed_node.oid = 'Node37103344'
    nd_37103408 = tree_37103088.seed_node.new_child(label="Node37103408", taxon=None, edge_length=8.43656365605, oid="Node37103408")
    nd_37103408.edge.oid = "Edge37103440"
    nd_37163888 = tree_37103088.seed_node.new_child(label="Python sebae", taxon=tax_36961776, edge_length=42.025995763, oid="Node37163888")
    nd_37163888.edge.oid = "Edge37164176"
    nd_37714672 = nd_37103408.new_child(label="Node37714672", taxon=None, edge_length=6.41982933271, oid="Node37714672")
    nd_37714672.edge.oid = "Edge37103472"
    nd_37163696 = nd_37103408.new_child(label="Node37163696", taxon=None, edge_length=10.1890732074, oid="Node37163696")
    nd_37163696.edge.oid = "Edge37163952"
    nd_37102896 = nd_37714672.new_child(label="Node37102896", taxon=None, edge_length=4.2262164672, oid="Node37102896")
    nd_37102896.edge.oid = "Edge37103536"
    nd_37163504 = nd_37714672.new_child(label="Node37163504", taxon=None, edge_length=20.5343771273, oid="Node37163504")
    nd_37163504.edge.oid = "Edge37163760"
    nd_37103120 = nd_37102896.new_child(label="Node37103120", taxon=None, edge_length=0.983173088098, oid="Node37103120")
    nd_37103120.edge.oid = "Edge37103600"
    nd_37103280 = nd_37102896.new_child(label="Node37103280", taxon=None, edge_length=1.53468865192, oid="Node37103280")
    nd_37103280.edge.oid = "Edge37163344"
    nd_37103504 = nd_37103120.new_child(label="Node37103504", taxon=None, edge_length=9.50779575096, oid="Node37103504")
    nd_37103504.edge.oid = "Edge37103664"
    nd_37104432 = nd_37103120.new_child(label="Node37104432", taxon=None, edge_length=3.20592158422, oid="Node37104432")
    nd_37104432.edge.oid = "Edge37104272"
    nd_37103568 = nd_37103504.new_child(label="Node37103568", taxon=None, edge_length=4.13809660753, oid="Node37103568")
    nd_37103568.edge.oid = "Edge37103728"
    nd_37104112 = nd_37103504.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=12.452417468, oid="Node37104112")
    nd_37104112.edge.oid = "Edge37104176"
    nd_37103632 = nd_37103568.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=8.31432086046, oid="Node37103632")
    nd_37103632.edge.oid = "Edge37103792"
    nd_37103888 = nd_37103568.new_child(label="Node37103888", taxon=None, edge_length=3.80591263593, oid="Node37103888")
    nd_37103888.edge.oid = "Edge37103760"
    nd_37103824 = nd_37103888.new_child(label="Node37103824", taxon=None, edge_length=2.88132393295, oid="Node37103824")
    nd_37103824.edge.oid = "Edge37103856"
    nd_37104304 = nd_37103888.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=4.50840822453, oid="Node37104304")
    nd_37104304.edge.oid = "Edge37104016"
    nd_37103920 = nd_37103824.new_child(label="Node37103920", taxon=None, edge_length=0.71440797492, oid="Node37103920")
    nd_37103920.edge.oid = "Edge37103984"
    nd_37104240 = nd_37103824.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=1.62708429158, oid="Node37104240")
    nd_37104240.edge.oid = "Edge37104208"
    nd_37103696 = nd_37103920.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=0.912676316658, oid="Node37103696")
    nd_37103696.edge.oid = "Edge37104048"
    nd_37104144 = nd_37103920.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=0.912676316658, oid="Node37104144")
    nd_37104144.edge.oid = "Edge37104080"
    nd_37104336 = nd_37104432.new_child(label="Node37104336", taxon=None, edge_length=1.50938344234, oid="Node37104336")
    nd_37104336.edge.oid = "Edge37104368"
    nd_37104784 = nd_37104432.new_child(label="Node37104784", taxon=None, edge_length=0.497596402447, oid="Node37104784")
    nd_37104784.edge.oid = "Edge37104624"
    nd_37103952 = nd_37104336.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=17.2449081924, oid="Node37103952")
    nd_37103952.edge.oid = "Edge37104496"
    nd_37104592 = nd_37104336.new_child(label="Node37104592", taxon=None, edge_length=4.9943499411, oid="Node37104592")
    nd_37104592.edge.oid = "Edge37104560"
    nd_37104464 = nd_37104592.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=12.2505582513, oid="Node37104464")
    nd_37104464.edge.oid = "Edge37104400"
    nd_37104720 = nd_37104592.new_child(label="Node37104720", taxon=None, edge_length=5.95835956047, oid="Node37104720")
    nd_37104720.edge.oid = "Edge37104656"
    nd_37104528 = nd_37104720.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=6.29219869082, oid="Node37104528")
    nd_37104528.edge.oid = "Edge37104688"
    nd_37104816 = nd_37104720.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=6.29219869082, oid="Node37104816")
    nd_37104816.edge.oid = "Edge37104848"
    nd_37104752 = nd_37104784.new_child(label="Node37104752", taxon=None, edge_length=5.19337056498, oid="Node37104752")
    nd_37104752.edge.oid = "Edge37104944"
    nd_37104976 = nd_37104784.new_child(label="Node37104976", taxon=None, edge_length=2.40513075673, oid="Node37104976")
    nd_37104976.edge.oid = "Edge37105136"
    nd_37104912 = nd_37104752.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=13.0633246673, oid="Node37104912")
    nd_37104912.edge.oid = "Edge37105008"
    nd_37105072 = nd_37104752.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=13.0633246673, oid="Node37105072")
    nd_37105072.edge.oid = "Edge37105104"
    nd_37104880 = nd_37104976.new_child(label="Node37104880", taxon=None, edge_length=2.06371429021, oid="Node37104880")
    nd_37104880.edge.oid = "Edge37105168"
    nd_37105584 = nd_37104976.new_child(label="Node37105584", taxon=None, edge_length=0.641957440329, oid="Node37105584")
    nd_37105584.edge.oid = "Edge37103184"
    nd_37105200 = nd_37104880.new_child(label="Node37105200", taxon=None, edge_length=3.43791705079, oid="Node37105200")
    nd_37105200.edge.oid = "Edge37105264"
    nd_37105648 = nd_37104880.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=13.7878501853, oid="Node37105648")
    nd_37105648.edge.oid = "Edge37105616"
    nd_37105040 = nd_37105200.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=10.3499331345, oid="Node37105040")
    nd_37105040.edge.oid = "Edge37105328"
    nd_37105392 = nd_37105200.new_child(label="Node37105392", taxon=None, edge_length=6.96420024077, oid="Node37105392")
    nd_37105392.edge.oid = "Edge37105424"
    nd_37105456 = nd_37105392.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=3.38573289378, oid="Node37105456")
    nd_37105456.edge.oid = "Edge37105360"
    nd_37105520 = nd_37105392.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=3.38573289378, oid="Node37105520")
    nd_37105520.edge.oid = "Edge37105552"
    nd_37105488 = nd_37105584.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=15.2096070352, oid="Node37105488")
    nd_37105488.edge.oid = "Edge37105296"
    nd_37105232 = nd_37105584.new_child(label="Node37105232", taxon=None, edge_length=1.5569144386, oid="Node37105232")
    nd_37105232.edge.oid = "Edge37163120"
    nd_37163056 = nd_37105232.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=13.6526925966, oid="Node37163056")
    nd_37163056.edge.oid = "Edge37163152"
    nd_37163248 = nd_37105232.new_child(label="Node37163248", taxon=None, edge_length=8.73106140586, oid="Node37163248")
    nd_37163248.edge.oid = "Edge37163088"
    nd_37163184 = nd_37163248.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=4.92163119076, oid="Node37163184")
    nd_37163184.edge.oid = "Edge37163216"
    nd_37163376 = nd_37163248.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=4.92163119076, oid="Node37163376")
    nd_37163376.edge.oid = "Edge37163280"
    nd_37163312 = nd_37103280.new_child(label="Node37163312", taxon=None, edge_length=6.61343105392, oid="Node37163312")
    nd_37163312.edge.oid = "Edge37163408"
    nd_37163728 = nd_37103280.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=21.4086976551, oid="Node37163728")
    nd_37163728.edge.oid = "Edge37163600"
    nd_37163440 = nd_37163312.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=14.7952666012, oid="Node37163440")
    nd_37163440.edge.oid = "Edge37163536"
    nd_37163632 = nd_37163312.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=14.7952666012, oid="Node37163632")
    nd_37163632.edge.oid = "Edge37163568"
    nd_37163664 = nd_37163504.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=6.63522564694, oid="Node37163664")
    nd_37163664.edge.oid = "Edge37163472"
    nd_37163856 = nd_37163504.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=6.63522564694, oid="Node37163856")
    nd_37163856.edge.oid = "Edge37163920"
    nd_37163792 = nd_37163696.new_child(label="Python molurus", taxon=tax_36961648, edge_length=23.4003588995, oid="Node37163792")
    nd_37163792.edge.oid = "Edge37163984"
    nd_37164080 = nd_37163696.new_child(label="Node37164080", taxon=None, edge_length=0.0534325174359, oid="Node37164080")
    nd_37164080.edge.oid = "Edge37164016"
    nd_37163824 = nd_37164080.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=23.3469263821, oid="Node37163824")
    nd_37163824.edge.oid = "Edge37164048"
    nd_37164208 = nd_37164080.new_child(label="Python regius", taxon=tax_36961712, edge_length=23.3469263821, oid="Node37164208")
    nd_37164208.edge.oid = "Edge37164144"
    tree_37164272 = dendropy.Tree(label="Tree09", taxon_set=tree_list.taxon_set, oid="Tree37164272")
    tree_list.append(tree_37164272, reindex_taxa=False)
    tree_37164272.seed_node.oid = 'Node37164432'
    nd_37164496 = tree_37164272.seed_node.new_child(label="Node37164496", taxon=None, edge_length=4.51223406477, oid="Node37164496")
    nd_37164496.edge.oid = "Edge37164528"
    nd_37221168 = tree_37164272.seed_node.new_child(label="Python sebae", taxon=tax_36961776, edge_length=29.3089629268, oid="Node37221168")
    nd_37221168.edge.oid = "Edge37220560"
    nd_37164368 = nd_37164496.new_child(label="Node37164368", taxon=None, edge_length=3.29007229288, oid="Node37164368")
    nd_37164368.edge.oid = "Edge37164560"
    nd_37220976 = nd_37164496.new_child(label="Node37220976", taxon=None, edge_length=1.61144541412, oid="Node37220976")
    nd_37220976.edge.oid = "Edge37220880"
    nd_37775760 = nd_37164368.new_child(label="Node37775760", taxon=None, edge_length=4.09158345085, oid="Node37775760")
    nd_37775760.edge.oid = "Edge37164624"
    nd_37220432 = nd_37164368.new_child(label="Node37220432", taxon=None, edge_length=12.3527736124, oid="Node37220432")
    nd_37220432.edge.oid = "Edge37220752"
    nd_37164304 = nd_37775760.new_child(label="Node37164304", taxon=None, edge_length=2.82624159524, oid="Node37164304")
    nd_37164304.edge.oid = "Edge37164688"
    nd_37165712 = nd_37775760.new_child(label="Node37165712", taxon=None, edge_length=0.0769064358513, oid="Node37165712")
    nd_37165712.edge.oid = "Edge37165648"
    nd_37164592 = nd_37164304.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=14.5888315231, oid="Node37164592")
    nd_37164592.edge.oid = "Edge37164752"
    nd_37164848 = nd_37164304.new_child(label="Node37164848", taxon=None, edge_length=0.0149913671254, oid="Node37164848")
    nd_37164848.edge.oid = "Edge37164656"
    nd_37164720 = nd_37164848.new_child(label="Node37164720", taxon=None, edge_length=4.89600409139, oid="Node37164720")
    nd_37164720.edge.oid = "Edge37164784"
    nd_37165264 = nd_37164848.new_child(label="Node37165264", taxon=None, edge_length=0.516331095919, oid="Node37165264")
    nd_37165264.edge.oid = "Edge37165136"
    nd_37164880 = nd_37164720.new_child(label="Node37164880", taxon=None, edge_length=5.06534441374, oid="Node37164880")
    nd_37164880.edge.oid = "Edge37164944"
    nd_37165200 = nd_37164720.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=9.67783606456, oid="Node37165200")
    nd_37165200.edge.oid = "Edge37165072"
    nd_37164816 = nd_37164880.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=4.61249165082, oid="Node37164816")
    nd_37164816.edge.oid = "Edge37165008"
    nd_37165104 = nd_37164880.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=4.61249165082, oid="Node37165104")
    nd_37165104.edge.oid = "Edge37164976"
    nd_37165232 = nd_37165264.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=14.05750906, oid="Node37165232")
    nd_37165232.edge.oid = "Edge37165040"
    nd_37165328 = nd_37165264.new_child(label="Node37165328", taxon=None, edge_length=6.44901536169, oid="Node37165328")
    nd_37165328.edge.oid = "Edge37165360"
    nd_37165392 = nd_37165328.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=7.60849369834, oid="Node37165392")
    nd_37165392.edge.oid = "Edge37165296"
    nd_37165456 = nd_37165328.new_child(label="Node37165456", taxon=None, edge_length=5.74491098314, oid="Node37165456")
    nd_37165456.edge.oid = "Edge37165488"
    nd_37165520 = nd_37165456.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=1.86358271521, oid="Node37165520")
    nd_37165520.edge.oid = "Edge37165424"
    nd_37165584 = nd_37165456.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=1.86358271521, oid="Node37165584")
    nd_37165584.edge.oid = "Edge37165552"
    nd_37165616 = nd_37165712.new_child(label="Node37165616", taxon=None, edge_length=8.7380438705, oid="Node37165616")
    nd_37165616.edge.oid = "Edge37165680"
    nd_37165968 = nd_37165712.new_child(label="Node37165968", taxon=None, edge_length=0.894040812546, oid="Node37165968")
    nd_37165968.edge.oid = "Edge37165168"
    nd_37164912 = nd_37165616.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=8.60012281197, oid="Node37164912")
    nd_37164912.edge.oid = "Edge37165776"
    nd_37165872 = nd_37165616.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=8.60012281197, oid="Node37165872")
    nd_37165872.edge.oid = "Edge37165744"
    nd_37165840 = nd_37165968.new_child(label="Node37165840", taxon=None, edge_length=1.98993837612, oid="Node37165840")
    nd_37165840.edge.oid = "Edge37165904"
    nd_37166576 = nd_37165968.new_child(label="Node37166576", taxon=None, edge_length=2.85390092474, oid="Node37166576")
    nd_37166576.edge.oid = "Edge37166736"
    nd_37165808 = nd_37165840.new_child(label="Node37165808", taxon=None, edge_length=1.38724843406, oid="Node37165808")
    nd_37165808.edge.oid = "Edge37166032"
    nd_37166384 = nd_37165840.new_child(label="Node37166384", taxon=None, edge_length=3.48386820234, oid="Node37166384")
    nd_37166384.edge.oid = "Edge37166256"
    nd_37165936 = nd_37165808.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=13.0669390597, oid="Node37165936")
    nd_37165936.edge.oid = "Edge37166096"
    nd_37166192 = nd_37165808.new_child(label="Node37166192", taxon=None, edge_length=3.7285329013, oid="Node37166192")
    nd_37166192.edge.oid = "Edge37166064"
    nd_37166128 = nd_37166192.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=9.33840615844, oid="Node37166128")
    nd_37166128.edge.oid = "Edge37166160"
    nd_37166288 = nd_37166192.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=9.33840615844, oid="Node37166288")
    nd_37166288.edge.oid = "Edge37166352"
    nd_37166416 = nd_37166384.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=10.9703192915, oid="Node37166416")
    nd_37166416.edge.oid = "Edge37166320"
    nd_37166512 = nd_37166384.new_child(label="Node37166512", taxon=None, edge_length=3.83596797183, oid="Node37166512")
    nd_37166512.edge.oid = "Edge37166448"
    nd_37166000 = nd_37166512.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=7.13435131963, oid="Node37166000")
    nd_37166000.edge.oid = "Edge37166480"
    nd_37166608 = nd_37166512.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=7.13435131963, oid="Node37166608")
    nd_37166608.edge.oid = "Edge37166640"
    nd_37166544 = nd_37166576.new_child(label="Node37166544", taxon=None, edge_length=3.52728125011, oid="Node37166544")
    nd_37166544.edge.oid = "Edge37166672"
    nd_37220720 = nd_37166576.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=13.5902249452, oid="Node37220720")
    nd_37220720.edge.oid = "Edge37220624"
    nd_37166704 = nd_37166544.new_child(label="Node37166704", taxon=None, edge_length=1.7419442928, oid="Node37166704")
    nd_37166704.edge.oid = "Edge37166800"
    nd_37220496 = nd_37166544.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=10.0629436951, oid="Node37220496")
    nd_37220496.edge.oid = "Edge37220400"
    nd_37166224 = nd_37166704.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=8.32099940228, oid="Node37166224")
    nd_37166224.edge.oid = "Edge37166864"
    nd_37166960 = nd_37166704.new_child(label="Node37166960", taxon=None, edge_length=4.44695529056, oid="Node37166960")
    nd_37166960.edge.oid = "Edge37166832"
    nd_37166896 = nd_37166960.new_child(label="Node37166896", taxon=None, edge_length=2.02523285213, oid="Node37166896")
    nd_37166896.edge.oid = "Edge37166928"
    nd_37220592 = nd_37166960.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=3.87404411171, oid="Node37220592")
    nd_37220592.edge.oid = "Edge37220464"
    nd_37166992 = nd_37166896.new_child(label="Node37166992", taxon=None, edge_length=0.774826953003, oid="Node37166992")
    nd_37166992.edge.oid = "Edge37167056"
    nd_37164112 = nd_37166896.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=1.84881125958, oid="Node37164112")
    nd_37164112.edge.oid = "Edge37220528"
    nd_37166768 = nd_37166992.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=1.07398430658, oid="Node37166768")
    nd_37166768.edge.oid = "Edge37164240"
    nd_37167024 = nd_37166992.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=1.07398430658, oid="Node37167024")
    nd_37167024.edge.oid = "Edge37167088"
    nd_37220688 = nd_37220432.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=9.1538829568, oid="Node37220688")
    nd_37220688.edge.oid = "Edge37220656"
    nd_37220848 = nd_37220432.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=9.1538829568, oid="Node37220848")
    nd_37220848.edge.oid = "Edge37220784"
    nd_37220944 = nd_37220976.new_child(label="Node37220944", taxon=None, edge_length=0.5196540118, oid="Node37220944")
    nd_37220944.edge.oid = "Edge37220912"
    nd_37221232 = nd_37220976.new_child(label="Python molurus", taxon=tax_36961648, edge_length=23.1852834479, oid="Node37221232")
    nd_37221232.edge.oid = "Edge37221008"
    nd_37220816 = nd_37220944.new_child(label="Python regius", taxon=tax_36961712, edge_length=22.6656294361, oid="Node37220816")
    nd_37220816.edge.oid = "Edge37221040"
    nd_37221136 = nd_37220944.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=22.6656294361, oid="Node37221136")
    nd_37221136.edge.oid = "Edge37221072"
    tree_37221104 = dendropy.Tree(label="Tree10", taxon_set=tree_list.taxon_set, oid="Tree37221104")
    tree_list.append(tree_37221104, reindex_taxa=False)
    tree_37221104.seed_node.oid = 'Node37221424'
    nd_37221488 = tree_37221104.seed_node.new_child(label="Python sebae", taxon=tax_36961776, edge_length=41.037475788, oid="Node37221488")
    nd_37221488.edge.oid = "Edge37221520"
    nd_37221584 = tree_37221104.seed_node.new_child(label="Node37221584", taxon=None, edge_length=0.971031599202, oid="Node37221584")
    nd_37221584.edge.oid = "Edge37221360"
    nd_37221616 = nd_37221584.new_child(label="Node37221616", taxon=None, edge_length=7.88042810387, oid="Node37221616")
    nd_37221616.edge.oid = "Edge37221200"
    nd_37221936 = nd_37221584.new_child(label="Node37221936", taxon=None, edge_length=11.1348383801, oid="Node37221936")
    nd_37221936.edge.oid = "Edge37221552"
    nd_37832752 = nd_37221616.new_child(label="Python regius", taxon=tax_36961712, edge_length=32.1860160849, oid="Node37832752")
    nd_37832752.edge.oid = "Edge37221680"
    nd_37221776 = nd_37221616.new_child(label="Node37221776", taxon=None, edge_length=3.37816716317, oid="Node37221776")
    nd_37221776.edge.oid = "Edge37221648"
    nd_37221712 = nd_37221776.new_child(label="Python molurus", taxon=tax_36961648, edge_length=28.8078489217, oid="Node37221712")
    nd_37221712.edge.oid = "Edge37221744"
    nd_37221904 = nd_37221776.new_child(label="Python brongersmai", taxon=tax_36961584, edge_length=28.8078489217, oid="Node37221904")
    nd_37221904.edge.oid = "Edge37221840"
    nd_37221968 = nd_37221936.new_child(label="Node37221968", taxon=None, edge_length=5.61233091052, oid="Node37221968")
    nd_37221968.edge.oid = "Edge37221808"
    nd_37224432 = nd_37221936.new_child(label="Node37224432", taxon=None, edge_length=19.6159722746, oid="Node37224432")
    nd_37224432.edge.oid = "Edge37278000"
    nd_37221872 = nd_37221968.new_child(label="Node37221872", taxon=None, edge_length=3.04566241261, oid="Node37221872")
    nd_37221872.edge.oid = "Edge37222064"
    nd_37222736 = nd_37221968.new_child(label="Node37222736", taxon=None, edge_length=1.63793957099, oid="Node37222736")
    nd_37222736.edge.oid = "Edge37222896"
    nd_37222000 = nd_37221872.new_child(label="Morelia oenpelliensis", taxon=tax_36961456, edge_length=20.2736124855, oid="Node37222000")
    nd_37222000.edge.oid = "Edge37222128"
    nd_37222224 = nd_37221872.new_child(label="Node37222224", taxon=None, edge_length=7.80824403904, oid="Node37222224")
    nd_37222224.edge.oid = "Edge37222160"
    nd_37222096 = nd_37222224.new_child(label="Node37222096", taxon=None, edge_length=3.53762335708, oid="Node37222096")
    nd_37222096.edge.oid = "Edge37222192"
    nd_37222800 = nd_37222224.new_child(label="Morelia boeleni", taxon=tax_36960944, edge_length=12.4653684465, oid="Node37222800")
    nd_37222800.edge.oid = "Edge37222704"
    nd_37222256 = nd_37222096.new_child(label="Morelia tracyae", taxon=tax_36961392, edge_length=8.92774508942, oid="Node37222256")
    nd_37222256.edge.oid = "Edge37222320"
    nd_37222416 = nd_37222096.new_child(label="Node37222416", taxon=None, edge_length=0.98740476198, oid="Node37222416")
    nd_37222416.edge.oid = "Edge37222032"
    nd_37222288 = nd_37222416.new_child(label="Node37222288", taxon=None, edge_length=4.46813855011, oid="Node37222288")
    nd_37222288.edge.oid = "Edge37222352"
    nd_37222832 = nd_37222416.new_child(label="Morelia amethistina", taxon=tax_36960880, edge_length=7.94034032744, oid="Node37222832")
    nd_37222832.edge.oid = "Edge37222480"
    nd_37222448 = nd_37222288.new_child(label="Node37222448", taxon=None, edge_length=0.318849288103, oid="Node37222448")
    nd_37222448.edge.oid = "Edge37222512"
    nd_37222768 = nd_37222288.new_child(label="Morelia clastolepis", taxon=tax_36961264, edge_length=3.47220177733, oid="Node37222768")
    nd_37222768.edge.oid = "Edge37222544"
    nd_37222384 = nd_37222448.new_child(label="Morelia kinghorni", taxon=tax_36960752, edge_length=3.15335248923, oid="Node37222384")
    nd_37222384.edge.oid = "Edge37222576"
    nd_37222672 = nd_37222448.new_child(label="Morelia nauta", taxon=tax_36961328, edge_length=3.15335248923, oid="Node37222672")
    nd_37222672.edge.oid = "Edge37222608"
    nd_37222960 = nd_37222736.new_child(label="Node37222960", taxon=None, edge_length=1.17144499717, oid="Node37222960")
    nd_37222960.edge.oid = "Edge37222928"
    nd_37223632 = nd_37222736.new_child(label="Node37223632", taxon=None, edge_length=0.325200272563, oid="Node37223632")
    nd_37223632.edge.oid = "Edge37223536"
    nd_37222864 = nd_37222960.new_child(label="Node37222864", taxon=None, edge_length=6.09401457319, oid="Node37222864")
    nd_37222864.edge.oid = "Edge37223024"
    nd_37223248 = nd_37222960.new_child(label="Node37223248", taxon=None, edge_length=0.168673432469, oid="Node37223248")
    nd_37223248.edge.oid = "Edge37222992"
    nd_37222640 = nd_37222864.new_child(label="Liasis olivaceus", taxon=tax_36960816, edge_length=14.4158757568, oid="Node37222640")
    nd_37222640.edge.oid = "Edge37223088"
    nd_37223184 = nd_37222864.new_child(label="Node37223184", taxon=None, edge_length=10.1133978359, oid="Node37223184")
    nd_37223184.edge.oid = "Edge37223120"
    nd_37223056 = nd_37223184.new_child(label="Liasis fuscus", taxon=tax_36960720, edge_length=4.30247792095, oid="Node37223056")
    nd_37223056.edge.oid = "Edge37223152"
    nd_37223280 = nd_37223184.new_child(label="Liasis mackloti", taxon=tax_36960912, edge_length=4.30247792095, oid="Node37223280")
    nd_37223280.edge.oid = "Edge37223344"
    nd_37223376 = nd_37223248.new_child(label="Node37223376", taxon=None, edge_length=9.21506648422, oid="Node37223376")
    nd_37223376.edge.oid = "Edge37223216"
    nd_37223664 = nd_37223248.new_child(label="Apodora papuana", taxon=tax_36960848, edge_length=20.3412168975, oid="Node37223664")
    nd_37223664.edge.oid = "Edge37223504"
    nd_37223312 = nd_37223376.new_child(label="Python reticulatus", taxon=tax_36961904, edge_length=11.1261504133, oid="Node37223312")
    nd_37223312.edge.oid = "Edge37223472"
    nd_37223568 = nd_37223376.new_child(label="Python timoriensis", taxon=tax_36961520, edge_length=11.1261504133, oid="Node37223568")
    nd_37223568.edge.oid = "Edge37223440"
    nd_37223408 = nd_37223632.new_child(label="Node37223408", taxon=None, edge_length=4.80149248356, oid="Node37223408")
    nd_37223408.edge.oid = "Edge37223600"
    nd_37223760 = nd_37223632.new_child(label="Node37223760", taxon=None, edge_length=5.59816192539, oid="Node37223760")
    nd_37223760.edge.oid = "Edge37223952"
    nd_37223696 = nd_37223408.new_child(label="Leiopython albertisii", taxon=tax_36960688, edge_length=16.554642571, oid="Node37223696")
    nd_37223696.edge.oid = "Edge37223792"
    nd_37223856 = nd_37223408.new_child(label="Bothrochilus boa", taxon=tax_36960784, edge_length=16.554642571, oid="Node37223856")
    nd_37223856.edge.oid = "Edge37223920"
    nd_37223728 = nd_37223760.new_child(label="Node37223728", taxon=None, edge_length=1.51569357225, oid="Node37223728")
    nd_37223728.edge.oid = "Edge37223984"
    nd_37224304 = nd_37223760.new_child(label="Node37224304", taxon=None, edge_length=2.47499227209, oid="Node37224304")
    nd_37224304.edge.oid = "Edge37224272"
    nd_37223888 = nd_37223728.new_child(label="Node37223888", taxon=None, edge_length=0.59310677926, oid="Node37223888")
    nd_37223888.edge.oid = "Edge37224048"
    nd_37224240 = nd_37223728.new_child(label="Morelia viridis", taxon=tax_36961136, edge_length=14.242279557, oid="Node37224240")
    nd_37224240.edge.oid = "Edge37224400"
    nd_37223824 = nd_37223888.new_child(label="Morelia carinata", taxon=tax_36961072, edge_length=13.6491727777, oid="Node37223824")
    nd_37223824.edge.oid = "Edge37224112"
    nd_37224208 = nd_37223888.new_child(label="Node37224208", taxon=None, edge_length=7.69541700903, oid="Node37224208")
    nd_37224208.edge.oid = "Edge37224144"
    nd_37224080 = nd_37224208.new_child(label="Morelia spilota", taxon=tax_36961200, edge_length=5.95375576867, oid="Node37224080")
    nd_37224080.edge.oid = "Edge37224176"
    nd_37224336 = nd_37224208.new_child(label="Morelia bredli", taxon=tax_36961008, edge_length=5.95375576867, oid="Node37224336")
    nd_37224336.edge.oid = "Edge37224016"
    nd_37277744 = nd_37224304.new_child(label="Antaresia maculosa", taxon=tax_36960560, edge_length=13.2829808571, oid="Node37277744")
    nd_37277744.edge.oid = "Edge37221264"
    nd_37224368 = nd_37224304.new_child(label="Node37224368", taxon=None, edge_length=3.47453380669, oid="Node37224368")
    nd_37224368.edge.oid = "Edge37221296"
    nd_37277840 = nd_37224368.new_child(label="Antaresia perthensis", taxon=tax_36960592, edge_length=9.80844705042, oid="Node37277840")
    nd_37277840.edge.oid = "Edge37277808"
    nd_37277904 = nd_37224368.new_child(label="Node37277904", taxon=None, edge_length=3.99697653657, oid="Node37277904")
    nd_37277904.edge.oid = "Edge37277968"
    nd_37277936 = nd_37277904.new_child(label="Antaresia childreni", taxon=tax_36960624, edge_length=5.81147051386, oid="Node37277936")
    nd_37277936.edge.oid = "Edge37277872"
    nd_37278032 = nd_37277904.new_child(label="Antaresia stimsoni", taxon=tax_36960432, edge_length=5.81147051386, oid="Node37278032")
    nd_37278032.edge.oid = "Edge37278096"
    nd_37277776 = nd_37224432.new_child(label="Aspidites ramsayi", taxon=tax_36960528, edge_length=9.3156335341, oid="Node37277776")
    nd_37277776.edge.oid = "Edge37278128"
    nd_37278224 = nd_37224432.new_child(label="Aspidites melanocephalus", taxon=tax_36960656, edge_length=9.3156335341, oid="Node37278224")
    nd_37278224.edge.oid = "Edge37278288"

    # set labels of nodes with taxa to taxon label, else oid (for consistent
    # identification in debugging)
    for t in tree_list:
        t.assign_node_labels_from_taxon_or_oid()

    return tree_list

def reference_dna_dict():
    dna_dict = {}
    dna_dict['Python regius'] = 'ATGCCCCACCACTATATCCTAACCCTCTTCGGCCTTCTACCAGTAGCAACCAACATCTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTCTAATGTTACAAGTACTTACCGGCTTCTTCCTAGCTGTCCACTATACAGCAAACATCAACCTAGCATTTTCATCCATTATCCATATTACCCGTGACGTCCCCTACGGCTGACTAATACAAAACCTACACGCCATCGGCGCATCGATATTCTTTATCTGCATCTACATTCACATCGCACGAGGACTATACTACGGCTCCCACCTCAATAAAGAAACCTGGATATCAGGTATTACACTTCTCATCACACTGATGGCAACCGCCTTCTTCGGGTATGTACTCCCATGAGGACAAATATCCTTCTGAGCCGCAACAGTAATTACCAACCTACTCACTGCTGTACCGTACCTAGGCGCAACCATAACCACCTGATTATGAGGAGGGTTCGCAATCAACGACCCCACCCTTACACGATTTTTCGCACTACACTTCATCCTACCATTCGAGATTATTTCCCTGTCATCATTACACATTATTTTACTTCACGAAGAAGGATCTAGCAACCCACTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCGTATCACTCATACAAAGACCTCCTACTACTAACACTAATACTACTAACACTTATAATCACCGTCTCCTTCTTCCCAGATATCTTCAACGACCCAGACAACTTCTCAAAAGCCAACCCATTAGTCACCCCCCAACACATTAAACCAGAGTGATACTTCCTATTCGCCTATGGCATCCTACGATCAATCCGAAATAAGCTTGGAGGGGCCTTAGCTCTAGTAATATCAATTATAATTATACTAACAGCCCCACTCACACACACAGCCCACCTCCGCCCAATAACCTTCCGACCACTTTCACAACTAATATTTTGAACCCTAATTTCAACATTCATTACCATTACATGAGCCGCCATAAAACCAGTAGAACCCCCATACATCATTATCAGCCAAACAACTTCAACACTATACTTCACCTTCTTTATCTCAACACCCATCCTAGGGTGAGTT-------------------------'
    dna_dict['Python sebae'] = 'ATGCCACACCATTATATCTTAACCCTATTCGGACTCCTACCAGTAGCAACCAACATTTCAACATGATGAAATTTCGGCTCAATACTACTAACATGTTTAGCCTTACAAACGCTCACAGGCTTCTTCCTAGCTGTCCACTACACAGCAAACATTAACCTAGCATTCTCATCTATCATTCACATCATCCGTGACGTCCCACATGGCTGAATAATACAAAACCTGCACGCCATCGGCGCATCTATATTCTTTATTTGCATTTACATCCACATCGCACGAGGCCTATACTATGGATCCTATCTTAACAAAGAAACCTGAATATCAGGTATCACACTCCTCATCATCCTAATAGCAACCGCGTTCTTCGGCTACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACAGTAATCACCAACCTACTCACTGCTGTACCCTACCTAGGAACAACTCTAACAACCTGATTATGGGGAGGCTTCGCAATCAACGACCCAACCCTGACACGATTCTTCGCACTACACTTCATCTTACCATTCGCTATCATCTCTCTATCATCATTACACGTCATCCTACTACACGAAGAAGGCTCCAGCAACCCATTAGGTACCAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACTCATACAAAGACCTCCTTCTACTAACACTAATGATCCTTTCCCTACTAATCATCGTCTCATTCTTTCCAGATATTTTCAACGACCCAGACAACTTCTCAAAAGCCAATCCACTAGTCACACCCCAGCATATTAAACCAGAATGATACTTCCTATTCGCCTACGGAATCCTACGATCAATCCCAAACAAACTAGGAGGCGCCCTAGCCCTAGTAATATCAATCATAATTCTATTAACCATCCCATTCACACACACATCCACTATACGATCAATAACATTCCGACCATTATCACAACTAATATTCTGAACATTAGTATCAACATTCATCACTATCACATGAGCCGCCACAAAACCAGTAGAACCACCATTCATTATCATCAGCCAAGTAACCGCAACACTATACTTCACATTCTTCATTTCAACCCCCATCCTAGGATGACTAGAAAACAAAATAACAAATCACCCAT'
    dna_dict['Python brongersmai'] = '------------------------------------------------------------------------TTCGGTTCAATATTACTCACTTGCCTAGTCCTACAAGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACAGCAAACATCAACCTAGCATTTTCCTCTATTATACACATCACCCGCGACGTCCCATACGGCTGAATAATACAAAACTTACACGCYATCGGCGCATCTATATTTTTCATCTGCATCTATATCCACATCGCACGAGGACTATATTACGGCTCCTATCTCAATAAAGAAACCTGAATGTCTGGCATTACACTCCTCATCACACTAATAGCAACCGCTTTTTTCGGATATGTCCTCCCATGAGGACAGATGTCATTCTGAGCCGCAACCGTAATCACCAATCTACTAACTGCTGTACCATACCTAGGCACAACCCTAACAACCTGATTATGAGGAGGGTTCGCAATCAACGACCCCACCCTCACACGATTCTTTGCACTACACTTCATCCTACCTTTCGCAATCATCTCTTTATCATCACTACACATTATTCTCCTTCATGAAGAAGGATCTAGCAATCCACTAGGAACCAACCCCGACATCGACAAAATCCCATTCCACCCATACCACTCTCACAAAGACTTCCTCCTACTCACACTATTAATCCTTTTCCTATTTATCATTGTCTCCTTCTTCCCAGACATTTTTAATGATCCAGATAACTTCTCAAAAGCTAACCCCCTTGTCACACCCCAACACATTAAACCAGAATGATACTTCTTATTCGCTTACGGAATCCTACGATCCATCCCAAACAAACTAGGTGGCGCATTAGCATTAGTAATATCAATCATAATCTTATTTACCATCCCATTCACACACACAGCCCATCTTCGCCCTATAACCTTCCGACCGTTCTCACAACTAATATTCTGAACACTAATCTCAACATTCATCACCATCACATGGGCCGCCACAAAACCAGTAGAACCCCCATATATTATCATCAGCCAAGCAACTGCAGCATTATACTTCACCTTCTTCATCTCTACACCCCTCCTGGGCTGAATAGAAAATAAAATAACAAACACTCCCT'
    dna_dict['Antaresia maculosa'] = 'ATGCCCCACCACTACATTCTAACCCTATTTGGTCTTCTACCTGTAGCAACAAATATCTCAACATGATGAAACTTCGGCTCAATATTACTAACATGTCTGGCCCTACAAGTATTGACCGGATTCTTCTTAGCCATCCACTACACAGCAAACATCAACTTAGCATTCTCATCTATTATTCACATCACCCGAGATGTCCCATATGGCTGAATAATACAAAACCTACACGCCATCGGAGCCTCCATATTCTTCATTTGCATTTACATTCACATTGCACGAGGACTATACTACGGATCCTACCTCAATAAAGAAACCTGAATGTCTGGCATCACCCTTCTTATCACACTAATAGCAACAGCCTTCTTCGGTTACGTTCTCCCATGGGGACAGATATCATTCTGAGCCGCAACCGTAATCACAAACTTACTTACCGCCGTCCCATACCTAGGCRTATCACTAACAACATGATTATGGGGGGGCTTTGCAATCAATGATCCTACACTGACACGATTCTTCGCACTACATTTCATCCTACCATTCGCAATCATCTCCTTATCCTCACTACACATTATTTTACTTCACGAAGAAGGCTCTAGCAACCCATTAGGAACCAATCCAGACATCGACAAAATTCCATTCCACCCCTACCACTCTCACAAAGACCTGCTCCTACTCACATTTATAATTTTACTCTTATTCACTATCGTCTCATTTCTCCCAGACATTTTCAATGACCCAGACAACTTCTCAAAAGCTAATCCCCTAGTAACACCGCAACACATTAAACCAGAATGATATTTCCTATTCGCCTATGGCATCCTACGATCTATCCCTAATAAACTGGGAGGCGCACTAGCTCTAGTAATATCGATCATAATCCTATTCATCATTCCATTCACACACACAGCCCGCTTACGCCCCATAACCTTCCGCCCACTATCACAACTAATATTCTGAACATTAGTATCAACATTCGCCACCATTACATGAGCCGCCACAAAACCAGTAGAGCCACCATTTATCATCATCAGTCAAACAACTTCAATACTTTACTTCACATTCTTCCTATCTACCCCAATTCTAGGATGAATAGAAAACAAAATAATAAACATCTCCT'
    dna_dict['Python timoriensis'] = '------------------------------------------------------------------------TTCGGCTCACTACTATTAACATGTCTAGCCCTACAAGTATTAACTGGTTTTTTCCTAGCCGTTCACTACACAGCAAACATTAACCTGGCATTTTCATCCATCATTCACATCACCCGAGACGTCCCATACGGTTGAATGATACAAAACCTCCACGCCATCGGAGCATCCATATTTTTCATTTGTATTTACATCCACATCGCACGAGGCCTATACTACGGATCATATYTTAACAAAGAAACTTGAATATCAGGCATCACCCTACTCATCACATTAATAGCTACTGCTTTCTTCGGATATGTTCTTCCATGAGGACAAATATCATTCTGRGCCGCAACTGTAATTACAAACCTACTTACAGCCGTACCATACCTGGGCACATCATTAACAACCTGACTCTGAGGCGGATTTGCAATCAACGACCCAACTCTAACACGATTCTTCGCACTACACTTTATCCTACCATTCGCAATCATCTCACTATCCTCACTACACATTATCTTACTCCATGAAGAAGGTTCTAGCAACCCCCTAGGAACTAACCCAGACATCGATAAAATCCCATTCCACCCCTATCATTCCCACAAAGACTTCCTCTTACTAATACTAATAATTCTATTTTTATTCATTATCGTTTCATTCTTCCCAGATATTTTCAACGACCCAGACAATTTCTCAAAAGCTAACCCACTAGTGACACCACAACACATTAAACCAGAATGATACTTCCTATTTGCCTATGGCATCCTACGATCTATTCCCAATAAACTAGGAGGAGCCCTAGCCCTAGTAATATCTATTATAATTCTATTCACCATCCCATTCACACACACAGCCTATCTTCGTCCAATAACCTTTCGCCCTTTCTCACAATTCATATTCTGAACACTAATCACTACATTCATCACCATCACATGAGCCGCTACAAAACCTGTAGAACCACCATTCATTATTATCAGCCAAGCGACATCAACACTATATTTCACCTTCTTCATTTCAATCCCTCTTCTAGGCTGAATAGAAAACAAAATAATACACCTCAATT'
    dna_dict['Python molurus'] = 'ATGCCCCACCACTATATCCTAACCTTATTTGGCCTCCTACCAGTAGCAACCAACATCTCAACATGATGAAACTTCGGCTCAATACTATTAGCATGCTTAGCCCTACAAGTATTAACCGGATTCTTCCTAGCCGTCCACTACACAGCAAACATCAACCTAGCATTCTCATCTATCATTCACATCACCCGCGATGTTCCATACGGCTGAATAATACAAAACCTACACGCTATCGGCGCATCCATATTCTTCATCTGCATCTACATTCACATCGCACGAGGACTATACTACGGCTCCTATCTAAATAAAGAAACCTGAATATCCGGAATTACACTACTCATCACACTTATGGCAACCGCCTTCTTCGGATATGTCCTCCCATGAGGGCAAATATCATTCTGAGCTGCAACCGTAATTACCAACCTATTAACCGCCGTACCATACTTAGGCACAACCCTAACAACCTGGTTATGAGGAGGATTCGCAATCAATGATCCCACCCTCACACGATTTTTTGCACTACATTTCATCCTACCATTCGCAATCATCTCCATATCATCACTACACATCATCCTACTCCACGAAGAAGGATCTAGCAACCCACTAGGAACAAACCCAGACATCGACAAAATCCCATTCCACCCATACCACTCATACAAAGACCTACTCTTCCTGACCCTAATAATCCTATTTATACTCATCATCGTCTCATTCTTCCCTGATATCTTCAACGACCCAGACAACTTCTCAAAAGCCAATCCACTAGTTACACCCCAACACATTAAACCAGAGTGGTACTTCCTATTCGCCTATGGGATCCTACGATCCATCCCAAACAAACTAGGTGGCGCATTAGCCCTAGTAATATCAATCATAATCCTATTTATTATCCCATTCACACATACAGCCCACTTCCGCCCAATAACTTTCCGCCCACTATCACAACTAATGTTCTGAACACTAGTATCAACATTCATCACTATCACATGAGCCGCCACAAAACCAGTAGAACCTCCATATATCATCATTAGCCAAGTAACAGCAACACTATACTTTATCTTCTTCATCTCTATACCCCTCCTAGGATGAATTGAAAACAAAATAACAAACACCCCCT'
    dna_dict['Morelia carinata'] = '-----------------------------------------------------------------------CTTCGGCTCGATACTATTAACATGTTTAGCCCTACAAGTATTAACCGGCTTCTTCTTAGCTGTTCACTACACAGCAAACATTAACCTAGCATTCTCATCCATCATTCACATCACCCGAGACGTCCCATACGGCTGAATAATACAAAATCTGCACGCCATCGGAGCATCCATATTCTTCATCTGCATTTACATTCATATTGCACGAGGACTATACTATGGGTCTTACCTCAACAAAGAAACCTGAATATCTGGTATCACCCTGCTAATTATCCTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTACCCTACTTAGGCACATCACTAACAACCTGGCTTTGAGGCGGATTCGCAATCAATGACCCAACTCTAACACGATTCTTCGCACTACACTTCATCCTACCATTCGCTATTATCTCCCTATCCTCACTACACATTATCTTGCTTCACGAAGAAGGTTCTAGCAACCCCTTAGGAACCAACCCGGACATCGACAAAATCCCATTCCACCCCTATCACACCTACAAAGATCTTCTTCTACTAACAGTAATAATCCTATTTTTATTCATTATCGTTTCATTCTTCCCAGACATTTTTAATGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGATACTTTCTATTTGCCTACGGAATTCTACGATCCATCCCCAATAAACTGGGAGGAGCATTAGCCCTAGTAATATCAATTATAATCCTATTCACCATTCCATTTATACACACAGCCCATCTTCGCCCAATAACCTTCCGCCCACTATCMCAACTAATATTTTGAACACTAATCTCAACATTTATCACTATCACATGAGCTGCCACAAAACCAGTAGAACCCCCATTCATTACCATCAGCCAAGCAACTTCAGCGCTATACTTCACGTTCTTCCTAACCACCCCAATTCTAGGATGAGTAGAAAATAAAATARTAAACATTCCCT'
    dna_dict['Morelia boeleni'] = '------------------------------------------------------------------------TTCGGCTCTATACTATTAACATGCTTAGGCTTACAAGTAATAACCGGCTTCTTCCTAGCCGTACACTACACAGCAAACATCAACTTAGCATTCTCATCCATCATCCACATCACCCGAGACGTCCCATACGGCTGAATAATACAAAACTTGCACGCTATCGGAGCATCTATATTCTTCATCTGCATTTACATCCACATCGCACGAGGGTTGTACTACGGATCATACCTTAACAAAGAAACCTGAATATCTGGCATTACCCTACTTATCACATTAATAGCAACTGCCTTCTTTGGATACGTTCTCCCATGAGGACAAATATCATTCTGAGSSGCWRCMGTWATCACAAACCTACTCACTGCCATCCCTTATCTAGGCACATCACTAACAACTTGACTATGAGGCGGATTCGCAATCAATGATCCTACACTAACACGATTTTTCGCACTACACTTCATCCTTCCATTCGCAATCATTTCCTTATCCTCACTACACATCATCCTACTCCACGAAGAAGGTTCCAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCATACCACTCCTACAAAGATCTTCTCCTACTAACATTAATAACCCTGTTCTTATTTATCATCGTCTCATTCTTCCCAGATATTTTTAACGACCCAGACAACTTTTCAAAAGCTAATCCCCTAGTAACACCACAACACATCAAACCTGAGTGATACTTTCTATTCGCCTATGGCATCCTACGATCCATCCCCAACAAACTAGGAGGTGCATTAGCCCTAGTAATATCAATCATGATCCTGTTTACCATCCCGTTTACACATACAGCCCACCTCCGTCCTATAACCTTCCGTCCACTCTCACAACTAATATTCTGAATATTAGTATCAACATTCATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCACCATTCATTACCATCAGTCAAGTAACCTCAACACTTTACTTCACATTCTTCTTATCCATCCCCATCCTAGGATGGATAGAAAACAAAATAATAGACATTCCAT'
    dna_dict['Antaresia perthensis'] = '------------------------------------------------------------------------TTCGGMTCAATACTACTAACATGTTTAGCCTTACAAGTACTAACCGGCTTTTTTTTGGCCGTCCACTACACAGCAAACATCAACCTGGCATTTTCATCCATCATTCACATTACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTGCACGCCATCGGAGCATCCATATTCTTCATTTGCATCTACATTCACATTGCACGCGGACTCTACTACGGATCCTACCTCAACAAAGAAACCTGGATATCGGGAATTACCCTCCTCATCACACTGATAGCTACCGCCTTCTTCGGCTACGTCCTCCCATGAGGACAGATATCATTCTGAGCCGCAACAGTAATCACCAACCTACTCACCGCTGTACCCTACCTAGGCACATCACTAACAACCTGACTATGAGGGGGGTTCGCAATCAACGATCCCACCCTGACGCGATTCTTCGCACTACACTTCATCCTACCATTCGCAATCATTTCTCTATCATCCTTACACATTATCTTACTACACGAAGAGGGCTCCAGCAACCCACTAGGAACCAACCCAGACATTGACAAAATCCCATTCCACCCATACCACTCCCACAAAGACCTGCTCCTACTAACACTTATAATTTTATTCTTATTCATTATCCTATCGTTCTCCCCGGACATCTTCAACGACCCAGATAACTTCTCAAAAGCTAACCCCCTAGTAACACCACAACACATTAAACCAGAATGGTACTTCCTGTTCGCCTACGGAATTCTACGATCCATTCCAAATAAATTAGGAGGCGCACTAGCCCTTGTGATATCAATCATAATCCTATTTACCATCCCATTCATACACACCGCCCACCTACGCCCAATAACCTTCCGCCCACTTTCACAACTTATATTCTGAACACTAGTCTCAACATTTGCCACCATTACATGAGCCGCCACAAAACCGGTAGAGCCCCCATACATCCTCATTAGCCAAGTGACCGCAACACTATACTTCACATTCTTTCTATCCATCCCAATCCTAGGATGAATAGAAAACAAAATAATAAACACCTCCT'
    dna_dict['Morelia viridis'] = '------------------------------------------------------------------------TTCGGYTCAATACTATTAACATGCCTAGCCCTACAAGTATTAACCGGCTTCTTCCTAGCCGTTCACTACACAGSAAACATTAATCTAGCATTCTCATCCATCATCCACATCTCCCGAGATGTTCCATACGGTTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATCTACATCCATATTGCACGAGGATTATACTATGGATCCTACCTCAACAAAGAAACCTGAATATCCGGTATTACCCTACTCATCACACTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTTACCGCTGTACCCTACCTGGGTACATCACTAACAACCTGATTATGAGGCGGATTCGCAATCAACGACCCAACCCTAACACGATTTTTTGCACTGCACTTCATCCTACCATTCGCAATCATCTCTCTATCCTCACTTCATGTTATCCTACTCCACGAAGAAGGCTCCAGCAATCCACTAGGAACTAATCCAGATATCGATAAAATCCCATTTCATCCATACCACTCCTACAAAGACCTACTCCTACTAACACTAATAATCCTATTCTTATTCATCATCGTTTCATTCTTCCCAGACATTTTTAACGATCCGGACAACTTCTCAAAAGCTAACCCATTAGTAACACCACAACACATCAAACCAGAATGGTATTTCCTATTCGCCTACGGCATTCTACGATCCATCCCCAACAAACTAGGAGGCGCATTAGCCTTAGTAATATCAATCATAATCCTATTTACCATCCCATTTACACACACAGCCTACCTCCGCCCCATAACCTTTCGTCCACTATCACAATTAATATTTTGAACATTGGTTGCAACATTCGCCACTATTACATGGGCTGCCACAAAACCAGTAGAACCCCCATTTATCCTCATTAGCCAAGTAACTTCAACACTATATTTCACATTCTTCTTATCCATCCCAATTTTAGGATGAATGGAAAATAAAATAATAAACATCCCCT'
    dna_dict['Aspidites ramsayi'] = '------------------------------------------------------------------------TTCGGCTCAATACTACTAACATGCTTAGGTYTACAAGTACTAACCGGCTTYTTCCTAGCCGTCCACTACACCGCAAACATTAACCTGGCATTCTCATCTATCGTTCACATCACCCGAGATGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGTATCTACATTCACATTGCACGAGGATTATACTACGGATCCTACCTTAACAAAGAAACCTGAATATCGGGTATTACATTACTCATTACACTAATAGCAACCGCCTTCTTCGGATATGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTCACCGCCGTACCATACCTAGGTACATCTCTAACAACCTGACTGTGAGGAGGATTCGCAATCAATGATCCCACCCTAACACGATTCTTTGCGCTACACTTCATCCTACCATTCGCAATCATCTCCTTGTCCTCACTACACATTATCTTACTTCACGAAGAAGGCTCCAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCGTACCACTCCCACAAAGACCTCCTCCTCCTAACGCTAATAATTATATCCTTATTTATTATCACCTCGTTCTTCCCAGACATCTTTAACGACCCCGACAATTTCTCAAAAGCCAACCCCCTAGTAACACCACAGCACATTAAACCAGAGTGGTACTTCCTATTTGCCTACGGCATCCTACGATCCATCCCCAACAAACTAGGAGGCGCACTAGCCCTAGTAATATCAATCATAATCCTATTTACCATTCCATTCACACACACAGCCTACCTTCGCCCTATAACCTTCCGCCCCCTATCACAACTAATATTCTGAACACTAGTCTCAACATTCATCACCATTACATGAGCCGCCACAAAACCAGTAGAACCACCATTCATTATAATTAGCCAAGTAACCTCAACACTATACTTCATATTCTTCTTATCAACACCCATCCTAGGATGAATAGAAAATAAAATAATAAACATTTCAT'
    dna_dict['Aspidites melanocephalus'] = 'ATGCCCCACCACTACATCCTAACCCTATTTGGCCTTCTGCCTGTAGCAACTAACATCTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTTTAGGCCTACAAGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACCGCAAACATTAACCTGGCATTCTCATCTATCGTTCACATCTCCCGAGATGTCCCATACGGCTGAATAATACAAAACCTACATGCAATCGGAGCATCCATATTCTTCATCTGTATCTACATTCACATTGCACGAGGATTATACTACGGATCCTACCTTAACAAAGAAACCTGAATATCAGGCATCACACTACTCATCACACTAATAGCGACCGCTTTCTTCGGATATGTGCTTACATGAGGACAAATATCATTATGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTGCCCTACCTAGGCACATCTCTAACAACCTGACTATGAGGAGGATTCGCAATCAATGATCCTACCCTAACACGATTCTTCGCACTCCACTTCATTCTGCCATTCGCAATCATCTCCTTATCCTCACTACACATCATCCTACTTCACGAAGAAGGCTCCAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCATTCCACCCATACCACTCCCACAAAGATCTCCTCCTCCTAACACTAATAATCATAAGCCTATTCATCATCAGCTCGTTCTTCCCAGATATCTTCAACGACCCCGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAGCACATTAAACCAGAATGATACTTCCTATTCGCCTACGGAATTCTACGATCCATTCCCAACAAACTAGGAGGCGCATTAGCATTAGTAATATCAATTATAATCCTATTTATTATCCCATTTACACACACAGCCCGCCTTCGCCCTATAACCTTCCGCCCTCTATCACAACTAATATTTTGAACACTAGTATCAACATTCATCACCATTACATGGGCCGCCACAAAACCAGTAGAACCACCATTCATTGTAATCAGCCAAGTAACCTCATCACTATACTTCACATTCTTCTTATCAACACCCATTCTGGGATGAGCAGAAAACAAAATAATAAACATTTCAT'
    dna_dict['Morelia oenpelliensis'] = '------------------------------------------------------------------------TTCGGCTCAATACTATTAACATGCCTAGCCCTACAAGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACAGCAAATATCAACCTAGCATTTTCATCCATTATCCACATCACCCGTGACGTCCCATACGGTTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATTTATATTCACATCGCTCGAGGACTATACTATGGGTCATACCTTAACAAAGAAACCTGAATATCCGGTATCACCCTACTCATTACACTAATAGCAACCGCCTTCTTCGGATATGTTCTTCCATGAGGACAGATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCCGTACCATACCTAGGCACATCACTAACAACATGACTATGAGGCGGATTCGCAATCAATGACCCAACCCTAACCCGATTCTTTGCATTACACTTCATCCTACCATTTGCAATTATCTCTTTATCCTCACTACATATCATCCTACTCCATGAAGAAGGTTCCAGCAATCCATTAGGAACCAATCCTGACATTGACAAAATCCCATTCCACCCATACCACTCCCACAAAGACCTACTCCTATTAACACTAATAACCCTACTCCTATTCATTATTGTCTCATTCTTCCCAGATATTTTCAACGACCCAGACAACTTCTCAAAAGCCAATCCCATAGTAACACCACAACACATCAAACCAGAATGGTACTTCCTATTTGCCTACGGCATCCTACGATCCATCCCAAACAAACTAGGAGGCGCATTAGCCCTAGTAATATCAATTATAATTCTATTCACCGCCCCATTCACACATACAGCCTACCTACGCCCTATAACCTTTCGCCCACTTTCACAACTAATATTCTGAGCACTAGTATCAACATTCATCACCATTACATGAGCCGCCACAAAACCAGTAGAACCACCATTTATCATCATCAGCCAAACAACTGCAACACTATACTTCACATTCTTCTTATCCATCCCCATCACAGGATGAATTGAAAACAAAATAATAAACACCCACT'
    dna_dict['Bothrochilus boa'] = '------------------------------------------------------------------------TTTGGCTCAATATTATTAACATGCCTGGCCCTACAAGTACTAACCGGCTTCTTCCTGGCCGTCCACTACACAGCAAACATCAACCTGGCATTCTCATCCATTATTCACATCACCCGAGATGTCCCATATGGCTGAATAATACAAAACCTGCACGCCATCGGAGCATCCATATTCTTCATTTGCGTATACATTCACATCGCACGAGGACTATACTACGGGTCATACCTAAACAAAGAAACCTGAATATCTGGCATTACCCTGCTCATCACACTAATAGCGACCGCCTTCTTTGGATATGTCCTCCCGTGAGGACAGATATCATTCTGAGCCGCAACCGTAATTACAAACCTGTTAACAGCAGTACCCTACCTGGGCACATCACTAACAACCTGGTTGTGAGGCGGGTTCGCAATCAATGACCCCACCCTAACACGATTCTTCGCACTACACTTCATCCTACCATTCGCAATCATCTCCCTATCCTCACTACACATCATCCTACTTCACGAGGAGGGCTCCAGCAACCCACTAGGAACCAACCCAGACATCGACAAAATCCCTTTCCACCCCTACCACTCCCACAAAGACTTTCTTCTTCTAACACTAATAACCCTATCCTTACTCATCATCGTCTTATTCTTCCCAGACATCTTTAACGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGATACTTCCTATTCGCCTACGGCATCCTACGTTCAATCCCCAATAAACTAGGAGGCGCACTAGCCTTAGTAATATCAATCATAATCCTATTCACCATCCCATTCACACACACAGCCCACCTTCGCCCCATAACCTTCCGTCCACTATCACAACTAATATTCTGAACATTAGTGTCAACATTTATCACTATCACATGGGCCGCCACAAAACCAGTAGAACCACCATTTATCACTATCAGCCAAACAACCTCAACACTATATTTTACATTCTTTTTACTTACCCCCATCCTAGGCTGAATAGAAAACAAAATAATAAAAACTCCCT'
    dna_dict['Morelia bredli'] = '------------------------------------------------------------------------TTCGGCTCAATACTATTAACATGCCTAGCCCTGCAAATCCTAACCGGCTTCTTTTTAGCGGTCCACTACACAGCAAACATCAACCTAGCATTCTCATCCATCATCCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAASCTACACGCCATCGGAGCATCCCTATTCTTCATCTGCATCTACATCCATATCGCACGTGGGTTATACTACGGATCCTATCTCAACAAAGAAACCTGAATATCCGGTATTACCCTACTCATCACACTAATAGCAACTGCCTTCTTCGGTTATGTCCTTCCATGAGGACAAATATCATTCTRRGCCGCAACTGTAATTACAAATCTACTCACCGCCGTACCATACCTGGGCACATCACTAACAACCTGACTATGAGGCGGATTCGCAATCAATGACCCCACCTTAACACGATTCTTCGCGCTACACTTCATCCTACCATTCGCAATCATCTCTCTCTCTTCACTACACATTATCTTACTTCACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACACCCACAAAGACCTGCTCCTACTAACGCTCATAATCCTGTTCCTATTCATTATCGTCTCATTCCTCCCAGATATCTTCAATGACCCAGACAACTTCTCAAAAGCTAACCCCTTGGTAACACCACAACACATTAAACCAGAGTGGTACTTCCTATTTGCCTATGGCATTCTACGATCCATCCCCAATAAACTAGGAGGCGCACTAGCCCTAATAATATCGATCCTAATTCTATTCACGATCCCATTCATACACACAGCCTATCTCCGCCCTATAACCTTCCGCCCCCTGTCACAACTTATATTTTGAACACTAATCTCAACATTCGCCACCATTACATGAGCTGCCACAAAACCAGTAGAACCCCCATTCATCATCATTAGCCAAGTAACCTCAACACTATACTTCACATTCTTCCTAACCATCCCAATTCTAGGGTGAATAGAAAACAAAATAATAAACATCTCCT'
    dna_dict['Morelia spilota'] = 'ATGCCCCACCACTACATCCTAACCTTATTTGGCCTTCTCCCCGTAGCAACCAATATCTCAACATGATGAAACTTCGGCTCAATACTATTAACATGCCTAGCCCTACAAGTTCTAACCGGCTTCTTCTTAGCTGTCCACTACACAGCAAACATTAACCTGGCATTCTCATCCATCATTCACATTACCCGAGACGTCCCATACGGCTGGATAATACAAAACCTACACGCCATCGGAGCATCTATATTCTTCATTTGCATCTACATCCATATTGCACGTGGATTATACTACGGATCCTATCTCAACAAAGAAACCTGAATATCCGGCATTACCCTACTCATCACACTAATAGCAACCGCCTTCTTCGGTTACGTCCTCCCATGAGGACAAATGTCATTCTAAGCCGCAACTGTAATTACAAACCTACTCACCGCCGTACCCTACCTAGGCACATCTCTAACAACCTGACTATGAGGCGGGTTCGCAATCAATGACCCCACCTTAACACGATTCTTCGCATTACACTTCATCCTACCATTCGCAATTATCTCTCTCTCCTCACTACACATTATTTTACTTCACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCTTACCACACCCACAAAGACCTACTCCTACTAACACTTATAATCCTGTTCTTATTCATTGTCGTCTCATTCCTCCCAGACATCTTTAATGACCCAGACAACTTCTCAAAAGCTAACCCTCTAGTAACACCACAGCACATTAAACCAGAGTGGTACTTCCTATTCGCCTATGGCATTCTACGATCCATCCCCAATAAATTAGGAGGCGCACTAGCCCTAGTAATATCAATCCTAATTCTATTCACAATCCCATTCATACACACAGCCTATCTCCGCCCCATAACCTTCCGCCCCCTGTCACAACTCATATTTTGAACACTAATCTCAACATTCGCCACCATTACATGAGCTGCCACAAAGCCAGTAGAACCCCCATTCATCATTATCAGCCAAGTAACCTCAACACTATACTTCACATTCTTCCTATCCATCCCTATTCTAGGGTGAATAGAAAACAAAATAATAAACATCTCCT'
    dna_dict['Antaresia stimsoni'] = '------------------------------------------------------------------------TTCGGCTCAATACTATTAACATGTCTAGCCCTACAAGTATTAACCGGCTTTTTCCTAGCCGTTCATTATACAGCAAACATTAACCTAGCATTTTCATCCATCGTTCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTTATTTGTATTTATATTCACATCGCACGCGGACTATACTATGGATCCTACCTCAACAAAGAAACCTGAATATCCGGTATCACCCTGCTCATCACACTAATAGCAACCGCCTTCTTCGGCTATGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTACCATACCTAGGCACATCGCTAACAACCTGACTGTGAGGGGGGTTCGCAATCAATGACCCCACCCTAACACGATTCTTCGCACTACACTTTATTCTACCATTCGCAATCATCTCCTTATCCTCCCTACACATTATCTTACTACACGAAGAAGGCTCAAGCAACCCACTAGGAACTAACCCAGACATCGACAAAATCCCATTCCACCCATACCACACCCACAAAGATCTACTCTTATTAACACTAATAGTTTTACTCCTATTCATTATCATTTCATTTTCCCCAGATATCTTTAATGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATCAAACCAGAATGGTACTTCCTATTTGCCTACGGTATCCTACGATCCATTCCCAACAAACTAGGAGGCGCACTAGCCTTAGTAATATCTATCATAATCCTATTCACCATCCCATTCACACACACAGCCCACCTACGCCCAATAACCTTCCGCCCACTCTCACAACTAACATTCTGAACACTGGTCTCAACATTCGCCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTTATCACCATCAGTCAAGTAACCTCAACACTATACTTCGCATTCTTCCTATCTATCCCAATTCTCGGATGAGTAGAAAACAAAATAATAAACATTTCAT'
    dna_dict['Antaresia childreni'] = 'ATGCCCCACCACTACATTCTAACCCTATTCGGCCTTCTGCCTGTAGCAACCAACATCTCAACATGATGAAACTTCGGCTCAATACTATTAACATGTCTAGCCCTACAAGTATTAACCGGTTTTTTCTTAGCTGTTCACTATACAGCAAACATTAACCTAGCATTTTCATCCATCGTTCACATTACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTTATTTGTATTTATATTCACATCGCACGCGGACTATACTATGGATCCTACCTCAACAAAGAAACCTGAATATCCGGTATCACCCTGCTCATCACACTAATAGCAACCGCCTTCTTCGGCTATGTTCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTACTCACCGCCGTACCATACCTGGGCACATCACTAACAACCTGACTGTGAGGAGGATTCGCAATCAATGACCCCACCCTAACGCGGTTCTTCGCGCTACACTTTATTCTACCATTCGCAATCATCTCCTTATCCTCCCTACACATTATTTTACTACACGAGGAAGGCTCCAGCAACCCACTAGGGACTAACCCAGACATCGATAAAATCCCATTCCACCCGTACCACACCCACAAAGACCTACTCCTACTAACACTAATAATTTTACTTCTATTAATTATCGTTTCATTTTCCCCGGACATCTTTCATGACCCAGACAACTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATTAAACCAGAATTATACTTCCTATTTGCCTACGGCATCCTACGATCTATCCCCAACAAACTAGGAGGCGCACTGGCCTTAGTAATATCCATCATAATCCTATTCACCATCCCATTCACACACACAGCCCACCTACGGCCAATAACCTTCCGCCCACTCTCACAACTAATATTTTGAGCACTAATTTCAACATTCGTCACCATCACATGGGCCGCCACAAAGCCAGTAGAACCCCCATTTATCATCATCAGTCAAGTAACCTCAACACTATACTTCACATTCTTCCTATCTATCCCAATTCTCGGATGAGTAGAAAACAAAATAATAAACATTTCAT'
    dna_dict['Leiopython albertisii'] = 'ATGCCCCACCACTACATTTTAACCCTATTTGGCCTCCTACCCGTAGCAACCAACATCTCAACATGATGAAACTTTGGTTCAATACTATTAACATGCTTAGCTCTACAGGTACTAACCGGCTTCTTCCTAGCCGTCCACTACACAGCAAACATCAACCTAGCATTTTCATCCATCATCCACATTACCCGAGATGTCCCATTCGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCCATATTCTTCATTTGCATCTACATTCACATCGCACGGGGGCTCTACTACGGATCATACCTAAACAAAGAAACCTGAATATCCGGCATTACCCTACTCATCACACTGATAGCAACCGCCTTCTTCGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCTGCAACCGTAATCACAAACCTACTAACCGCCGTACCCTACCTAGGCACATCACTAACAACCTGATTATGGGGCGGATTTGCAATCAACGACCCTACCCTAACACGATTCTTCGCACTACACTTCATCCTACCTTTCGCAATCATCTCCTTATCTTCACTACACATTATCCTTCTTCACGAAGAAGGCTCTAGCAACCCACTAGGAACCAATCCAGACATCGACAAAATCCCATTCCACCCCTATCACTCCCACAAAGATCTTCTCCTACTGACACTAATAATACTATCTCTGCTCATCATCGTCTCATTCTTCCCAGACATCTTTAATGACCCAGATAACTTCTCTAAAGCCAACCCATTAGTAACGCCACAACACATCAAACCAGAATGATACTTCCTATTCGCCTATGGTATCCTACGATCCATCCCCAATAAACTAGGAGGCGCACTAGCCCTAGTAATATCAATCATAATTTTATTCACCATCCCGTTCACACACACAGCCCATCTTCGCCCCATAACCTTCCGCCCATTCTCACAACTAATATTTTGAACACTAGTCTCAACATTTATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTTATCGTCATCAGCCAAGTAACTTCAACACTATACTTTACATTCTTCTTACTCATCCCCATTTTGGGCTGAACAGAAAATAAAATAATAAACACCCTCT'
    dna_dict['Python reticulatus'] = 'ATGCCCCACCATTATATCCTAACCTTATTTGGCCTTCTACCAGTAGCAACCAACATCTCAACCTGATGAAACTTCGGCTCAATATTACTAACATGTCTAGCCTTACAAGTACTAACCGGCTTTTTCCTAGCCGTCCATTACACAGCAAACATTAACCTAGCATTTTCATCCATCATCCACATCACCCGAGACGTCCCATACGGCTGAATAATACAAAACCTTCACGCTATCGGAGCATCCATATTCTTCATCTGCATCTACATCCACATCGCACGAGGCCTATACTACGGATCATACCTCAACAAAGAAACCTGAATATCAGGCATCACCCTACTCATCACACTAATAGCCACCGCTTTTTTTGGTTACGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTCACTGCCGTACCATACCTAGGTACATCACTAACAACCTGGCTTTGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGATTCTTTGCACTACATTTTATTTTACCATTCGCGATTATCTCATTATCCTCATTACACGTTATCTTACTCCACGAAGAAGGTTCTAGCAACCCCCTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCATACCACTCTTACAAAGACTTCCTCTTACTAACATTAATAGTCCTATCTCTATTCATTATCGTCTCATTCTTCCCAGATATCTTCAACGACCCAGACAATTTCTCAAAAGCCAACCCCCTAGTAACACCACAACACATTAAACCAGAATGATACTTCCTATTCGCTTACGGTATTCTACGATCCATCCCCAACCAACTTGGAGGAGCATTAGCCTTAGTAATATCTATTATAATCTTATTCACTATCCCATTCACACACACAGCTAATCTTCGTCCCATAACCTTCCGACCACTCTATCAACTCATGTTCTGAACACTAGTCTCCACTTTTATTACTATCACATGAGCCGCTACAAAACCCGTAGAGCCCCCTTTTATCACTATTAGTCAAGTAACTTCAACACTTTATTTCACATTCTTCATCTCCATCCCATTCCTAGGCTGAATAGAAAACAAAATAATACATCTCAACT'
    dna_dict['Morelia tracyae'] = '--------CCACTACATCCTAACCCTATTTGGCCTCCTACCAGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTCTAGCTCTACAAGTACTAACCGGCTTCTTTCTAGCCGTACACTACACAGCAAACATTAACCTAGCATTTTCATCCATCATTCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAATCTACACGCTATCGGAGCATCCATATTCTTCATTTGCATTTACATCCACATCGCACGAGGACTATACTACGGATCTTACCTAAACAAAGAAACTTGAATATCAGGCATTACCCTACTCATCACACTAATAGCAACTGCCTTCTTTGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCCATCCCATACCTAGGCACATCTCTAACAACCTGACTTTGAGGCGGATTCGCAATTAACGACCCTACCCTAACACGCTTTTTCGCATTACACTTCATCCTACCATTCGCAATCATTTCCCTATCCTCACTACACATCATTCTACTCCACGAAGAAGGTTCTAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCATTCCATCCGTACCACTCCCACAAAGACCTCCTTTTACTAACACTAATAATCTTATTTCTATTCATCATCGTTTCATTCTTTCCTGATATT-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    dna_dict['Morelia amethistina'] = '--------CCACTACATCCTAACCTTATTTGGCCTCCTACCGGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGCCTGGCACTACAAGTACTAACCGGCTTCTTCCTAGCCGTACACTACACAGCAAACATTAACCTAGCATTCTCATCCATCATCCACATCACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTGCACGCTATTGGAGCATCCATATTCTTCATCTGCATCTACATTCATATCGCACGAGGACTATACTACGGATCATACCTCAACAAAGAAACCTGAATATCCGGCATTACCCTGCTCATCACACTAATAGCAACCGCTTTCTTCGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCCATCCCATACCTAGGCACATCTCTAACAACCTGACTTTGAGGCGGATTTGCAATCAACGACCCCACCCTAACACGCTTTTTCGCATTACACTTCATTCTACCATTTGCAATCATTTCCTTATCCTCACTACATATCATCCTACTACACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGATATTGACAAAATCCCGTTCCACCCATACCACTCCTACAAAGACCTCCTCTTACTAACACTAATAATCCTATTCCTATTCACCATCGTTTCATTCTTCCCTGATATC-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    dna_dict['Morelia nauta'] = '--------CCACTACATCTTAACCTTATTTGGCCTCCTACCGGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGCCTAGCGCTACAAGTACTAACCGGCTTCTTCCTAGCCGTACACTACACAGCGAACATTAACCTAGCATTTTCATCCATCATCCACATTACCCGAGACGTCCCATATGGCTGAATAATACAGAACCTACACGCTATCGGAGCATCCATATTCTTCATCTGCATTTACATCCACATCGCACGAGGACTATACTACGGATCATACCTAAACAAAGAAACCTGAATATCCGGCATCACCCTGCTCATCACACTAATAGCAACTACCTTCTTCGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCTATTCCTTACCTAGGCACATCACTGACAACCTGACTTTGAGGCGGATTCGCAATCAACGACCCCACCTTAACACGTTTTTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCCTTATCCTCACTACACGTCATCCTACTTCACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCATACCACACCTACAAAGACCTCCTCTTACTAACACTAATAATCCTGTTCCTATTCATCATCGTTTCATTCTTCCCTGATATC-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    dna_dict['Morelia kinghorni'] = '--------CCACTACATCTTAACCTTATTTGGCCTCCTACCAGTAGCAACCAACATTTCAACATGATGAAACTCCGGCTCAATACTACTAACATGTCTAGCACTACAAGTACTAACCGGCTTCTCCCTAGCTGTACACTACACAGCGAACATTAACCTAGCATTTTCATCCATCATCCACATTACCCGAGACGTCCCATATGGCTGAATAATACAGAACCTACACGCTATCGGAGCATCCATATTCTTCATCTGCATTTACATCCACATCGCACGAGGACTATACTACGGATCATACCTAAACAAAGAAACCTGAATATCCGGCATCACCCTGCTCATCACACTAATAGCAACTGCCTTCTTCGGATACGTTCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCTATTCCTTACCTAGGCACATCACTGACAACCTGACTTTGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGTTTTTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCCTTATCCTCACTACACGTCATCCTACTACACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCATACCACACCTACAAAGACCTCCTCTTACTAACACTAATAGTCCTGTTCCTATTCATCATCGTTTCATTCTTCCCTGATATT-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    dna_dict['Morelia clastolepis'] = '--------CCACTACATCTTAACCCTATTTGGCCTCCTACCAGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGCCTAGCACTACAAGTACTAACCGGCTTCTTCCTAGCCGTACACTACACAGCGAACATTAACCTAGCATTCTCATCCATCATCCACATTACCCGAGACGTCCCATATGGCTGAATAATACAAAACCTACACGCTATCGGAGCATCCATATTCTTCATTTGCATTTACATTCACATCGCACGAGGACTATACTACGGATCTTACCTAAACAAAGAAACCTGAATATCCGGCATCACCCTGCTCATCACACTAATAGCAACTGCCTTCTTCGGATACGTCCTCCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTACTTACCGCTATTCCTTACCTAGGCACATCACTGACAACCTGACTTTGAGGCGGATTCGCAATCAACGACCCCACCCTTACACGTTTTTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCCTTATCCTCACTACACGTCATCCTACTACACGAAGAAGGCTCTAGCAACCCATTAGGAACCAACCCAGACATTGACAAAATCCCGTTCCACCCATACCACACCTACAAAGACCTCCTCTTACTAACACTAATAGTCCTGTTCCTATTCATCATCGTTTCATTCTTCCCTGATATC-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    dna_dict['Liasis fuscus'] = '------------------------------------------------ACCAATATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTTTAGCCCTACAAGTATTAACCGGATTCTTCCTGGCTGTCCACTATACAGCAAATATTGACCTGGCATTCTCATCCATCATCCACATCACTCGAGACGTCCCATACGGCTGAATAATACAAAACCTACACGCCATCGGAGCATCAATATTCTTCATTTGTATCTACATCCACATCGCCCGAGGCCTATACTACGGATCATACCTCAACAAAGAAACCTGAATATCCGGCATCACCCTACTTATCACACTAATAGCAACCGCCTTCTTCGGGTACGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATCACAAACCTTCTTACCGCCGTACCCTACCTAGGCACATCCTTGACAACCTGATTATGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGATTCTTCGCATTACACTTCATCCTACCATTCGCAATCATCTCTTTATCCTCACTACACGTCATCCTCCTCCACGAGGAGGGGTCTAGCAACCCACTAGGGACTAACCCAGACATCGACAAAATCCCATTCCACCCTTACCACTCCCACAAAGACCTTCTCCTACTAACACTAATAATAATATCCCTACTCATTATTGTTTCCTTTTTCCCAGACATCTTCAACGACCCAGACAACTTCTCAAAAGCCAACCCCTTGGTAACACCACAACACATTAAACCAGAATGATACTTCCTGTTCGCCTACGGCATCCTACGATCTATTCCCAACAAACTTGGAGGAGCATTAGCTCTAGTAATATCAATCATAATCTTATTTTCTACCCCATTCACACACACAGCCCACCTCCGCCCTATAACTTTCCGCCCACTATCCCAACTAATATTCTGAACACTAATCTCAACATTCATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTCATCATCATCAGCCAAGTAACCTCAATACTATACTTCACATTTTTCCTATCCATCCCCATTCTAGGATGGGTAGAGAACAAAATTATAAACACCCCCT'
    dna_dict['Liasis mackloti'] = 'ATGCCCCACCACTACGTTCTAACCCTATTTGGTCTCTTACCAGTAGCAACCAATATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGTTTAGCCCTACAAGTACTAACCGGATTCTTCCTGGCTGTCCACTACACAGCAAATATTAACCTGGCATTCTCATCCATCGTTCACATCACTCGAGATGTCCCATACGGCTGAATGATACAAAACCTACACGCCATCGGAGCATCTATATTCTTTATTTGTATCTACATCCACATCGCCCGAGGCCTATACTACGGATCATACCTTAACAAAGAAACCTGAATATCCGGTATTACCCTGCTTATCACACTAATAGCAACCGCCTTCTTCGGATACGTCCTTCCATGAGGACAAATATCATTCTGAGCCGCAACCGTAATTACAAACCTTCTCACCGCCGTACCCTACCTAGGCACATCCTTGACAACCTGGCTATGAGGGGGGTTCGCAATCAACGACCCCACCCTAACACGATTCTTAGCATTACACTTCATCCTACCATTCGCAATCATCTCTTTATCCTCCTCACACGTCATCCTCCTTCATGAAGAAGGGTCTAGCAACCCACTAGAAACTAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACTCCCACAAAGACCTTCTCCTACTAACACTAATAATAATATTCCTATTCATTATTGTTTCCTTTTTCCCAGACATCTTCAACGACCCAGACAATTTCTCAAAAGCTAACCCTCTAGTAACACCACAACACATTAAACCAGAGTGGTACTTCCTATTCGCCTACGGCATCCTACGATCTATCCCCAACAAACTTGGAGGAGCATTAGCCCTAGTAATATCAATCATAATCTTATTTTCTACCCCATTCACACACACAGCCCACCTCCGCCCCATAACTTTCCGCCCACTATCCCAACTAATATTCTGAACACTAGTCTCAACATTCATCACCATCACATGAGCCGCCACAAAACCAGTAGAACCCCCATTCATCACCATTAGCCAAGTAACCTCAATTCTATACTTCACATTTTTCCTATCCATCCCTATTCTAGGATGAGTAGAGAACAAAATTATAAACGCCCCCT'
    dna_dict['Liasis olivaceus'] = 'ATGCCCCACCACTACATTCTAACCCTGTTCGGCCTCTTACCAGTAGCAACCAACATTTCAACATGATGAAACTTCGGTTCAATACTACTAACATGCCTAGTCCTACAAGTATTAACCGGTTTCTTCCTAGCTGTCCACTACACAGCAAACATCAATCTAGCATTCTCATCCATCGTTCACATTACCCGAGACGTCCCATACGGCTGAATAATACAAAACCTACACGCTATCGGAGCATCTATATTCTTCATTTGCATCTACATCCATATCGCACGAGGTCTATACTACGGATCATACCTTAACAAAGAAACCTGAATATCTGGTATCACCCTACTCATCACACTAATAGCAACCGCTTTCTTCGGATATGTCCTTCCATGGGGACAAATATCATTCTGGGCCGCAACCGTAATCACAAACCTACTCACTGCCGTACCCTATCTAGGCACATCACTAACAACCTGACTTTGAGGCGGATTCGCAATCAACGACCCCACCCTAACACGATTCTTCGCACTACACTTCATCCTACCATTCGCAATCATCTCTTTATCCTCACTACACATCATCCTACTCCACGAAGAAGGATCTAGCAACCCATTAGGAACCAACCCAGACATCGACAAAATCCCATTCCACCCCTATCACTCACACAAAGACCTTCTCCTATTAACACTAATAATAATATTCCTATTCATTATCGTATCATTCTTCCCAGATATTTTCAACGACCCAGATAACTTCTCAAAAGCCAACCCCTTAGTAACACCACAACACATTAAACCAGAATGATACTTCCTATTCGCCTACGGCATTCTACGATCTATCCCCAACAAACTTGGAGGCGCATTAGCTCTAGTAGCATCAATCATAATCCTATTCACCACCCCATTCACACACACAGCCAACCTCCGCCCTATAACCTTCCGCCCCCTGTCACAACTAATATTCTGAACATTAGTCTCAACATTCATCACTATTACATGAGCCGCCACAAAACCAGTAGAACCCCCATTCATCGCTATCAGCCAAGTAACCTCAATACTCTATTTCACATTTTTCTTGTCCATCCCCATCCTAGGATGAATAGAAAACAAAATAATAAACACCCCCT'
    dna_dict['Apodora papuana'] = 'ATGCCCCACCATTACATCCTAACCCTGTTCGGCCTCCTACCAGTAGCAACCAACATTTCAACATGATGAAACTTCGGCTCAATACTACTAACATGCCTAGCCCTACAAGTATTAACTGGCTTCTTCCTGGCCGTACACTACACAGCAAACATCAACCTAGCATTCTCATCCATCATTCACATCACCCGAGATGTCCCATACGGCTGAATAATACAAAACTTACACGCCATCGGAGCATCCATATTCTTCATCTGTATCTACATCCATATTGCACGGGGCCTATACTACGGATCGTACCTAAATAAAGAAACCTGAATATCTGGCATCACCCTACTCATCACACTAATAGCAACCGCCTTCTTCGGATATGTCCTTCCATGAGGACAAATGTCATTCTGAGCCGCAACTGTAATCACAAATCTGCTCACTGCAGTACCCTACCTGGGTACATCACTAACAACTTGATTATGAGGCGGCTTCGCAATCAATGACCCCACCCTGACACGATTCTTCGCACTACACTTCATCCTACCATTCGAAATCATCTCTTTATCCTCACTACACATCATCCTACTTCATGAAGAAGGCTCTAGTAACCCATTGGGAACTAACCCAGACATCGACAAAATCCCATTCCACCCCTACCACTCTCACAAAGACCTTCTCCTATTAACACTAATAATTCTACTCCTATTCATCACCATATCATTCTTCCCAGATATTTTCAACGACCCAGACAACTTCTCAAAAGCCAACCCACTAGTAACACCACAACACATCAAACCAGAATGATACTTCCTATTCGCCTATGGGATCCTACGATCCATCCCCAACAAACTAGGGGGCGCATTAGCCCTAGTAACATCGATCATAATCCTATTCACCATCCCATTTACACACACAGCTCACCTCCGACCTATAACCTTCCGCCCCCTATCACAACTGATATTCTGAACATTAGTATCAACATTTATCACCATTACATGGGCCGCCACAAAACCAGTAGAACCCCCATTCATCATCATTAGCCAAGCAACCTCATTATTATACTTCACATTCTTCTTATCCTTCCCTATCCTAGGGTGGACAGAAAACAAAATAATAAACACCCCCT'
    return dna_dict

def reference_taxon_set():
    taxon_set = dendropy.TaxonSet()
    taxon_set.new_taxon(label="Antaresia childreni", oid="Taxon36960624")
    taxon_set.new_taxon(label="Antaresia maculosa", oid="Taxon36960560")
    taxon_set.new_taxon(label="Antaresia perthensis", oid="Taxon36960592")
    taxon_set.new_taxon(label="Antaresia stimsoni", oid="Taxon36960432")
    taxon_set.new_taxon(label="Aspidites melanocephalus", oid="Taxon36960656")
    taxon_set.new_taxon(label="Aspidites ramsayi", oid="Taxon36960528")
    taxon_set.new_taxon(label="Bothrochilus boa", oid="Taxon36960784")
    taxon_set.new_taxon(label="Leiopython albertisii", oid="Taxon36960688")
    taxon_set.new_taxon(label="Liasis fuscus", oid="Taxon36960720")
    taxon_set.new_taxon(label="Liasis mackloti", oid="Taxon36960912")
    taxon_set.new_taxon(label="Liasis olivaceus", oid="Taxon36960816")
    taxon_set.new_taxon(label="Apodora papuana", oid="Taxon36960848")
    taxon_set.new_taxon(label="Morelia amethistina", oid="Taxon36960880")
    taxon_set.new_taxon(label="Morelia boeleni", oid="Taxon36960944")
    taxon_set.new_taxon(label="Morelia bredli", oid="Taxon36961008")
    taxon_set.new_taxon(label="Morelia carinata", oid="Taxon36961072")
    taxon_set.new_taxon(label="Morelia clastolepis", oid="Taxon36961264")
    taxon_set.new_taxon(label="Morelia kinghorni", oid="Taxon36960752")
    taxon_set.new_taxon(label="Morelia nauta", oid="Taxon36961328")
    taxon_set.new_taxon(label="Morelia oenpelliensis", oid="Taxon36961456")
    taxon_set.new_taxon(label="Morelia spilota", oid="Taxon36961200")
    taxon_set.new_taxon(label="Morelia tracyae", oid="Taxon36961392")
    taxon_set.new_taxon(label="Morelia viridis", oid="Taxon36961136")
    taxon_set.new_taxon(label="Python brongersmai", oid="Taxon36961584")
    taxon_set.new_taxon(label="Python molurus", oid="Taxon36961648")
    taxon_set.new_taxon(label="Python regius", oid="Taxon36961712")
    taxon_set.new_taxon(label="Python reticulatus", oid="Taxon36961904")
    taxon_set.new_taxon(label="Python sebae", oid="Taxon36961776")
    taxon_set.new_taxon(label="Python timoriensis", oid="Taxon36961520")
    return taxon_set

def reference_dna_array(taxon_set=None):
    if taxon_set is None:
        taxon_set = reference_taxon_set()
    dna = dendropy.DnaCharacterArray(taxon_set=taxon_set)
    sa = dna.default_state_alphabet
    assert len(dna.taxon_set) == 29
    labels = labels=[ "Python regius",
                      "Python sebae",
                      "Python brongersmai",
                      "Antaresia maculosa",
                      "Python timoriensis",
                      "Python molurus",
                      "Morelia carinata",
                      "Morelia boeleni",
                      "Antaresia perthensis",
                      "Morelia viridis",
                      "Aspidites ramsayi",
                      "Aspidites melanocephalus",
                      "Morelia oenpelliensis",
                      "Bothrochilus boa",
                      "Morelia bredli",
                      "Morelia spilota",
                      "Antaresia stimsoni",
                      "Antaresia childreni",
                      "Leiopython albertisii",
                      "Python reticulatus",
                      "Morelia tracyae",
                      "Morelia amethistina",
                      "Morelia nauta",
                      "Morelia kinghorni",
                      "Morelia clastolepis",
                      "Liasis fuscus",
                      "Liasis mackloti",
                      "Liasis olivaceus",
                      "Apodora papuana"]
    assert dna.taxon_set.has_taxa(labels=labels)

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
    cells = [dendropy.CharacterDataCell(value=s, character_type=col_type) for s in col_type.state_alphabet.get_states(symbols=symbols)]
    return cells

def reference_standard_array(taxon_set=None):
    if taxon_set is None:
        taxon_set = reference_taxon_set()
    ca1 = dendropy.StandardCharacterArray(taxon_set=taxon_set)
    assert len(ca1.taxon_set) == 29
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
                         reference_dna_array(taxon_set=taxon_set),
                         reference_standard_array(taxon_set=taxon_set))
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
