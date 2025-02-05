#! /usr/bin/env python
# -*- coding: utf-8 -*-

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
Tests of protracted speciation model/process.
"""

import random
import itertools
import math
import unittest
import json
import dendropy
from dendropy.calculate import treecompare
from dendropy.model import protractedspeciation
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from support import pathmap

_ET_LINEAGE_ID = 0
_ET_PARENT_ID = 1
_ET_ORIGIN_T = 2
_ET_SPECIATION_T = 3
_ET_EXTINCTION_TIME = 4
_ET_SPECIES_ID = 5
_ETX_NODE = 6

def _load_json(filename):
    with open(pathmap.other_source_path(filename)) as src:
        j = json.load(src)
    return j

def _get_stock_psm():
    return protractedspeciation.ProtractedSpeciationProcess(
            speciation_initiation_from_orthospecies_rate=0.01,
            speciation_initiation_from_incipient_species_rate=0.01,
            speciation_completion_rate=0.01,
            orthospecies_extinction_rate=0.01,
            incipient_species_extinction_rate=0.01,
            lineage_label_format_template="S{species_id}-{species_id}-{lineage_id}", # for compatibility with PBD generated data
            species_label_format_template="S{species_id}-{species_id}-{lineage_id}", # for compatibility with PBD generated data
            )

class ProtractedSpeciationCalcs(unittest.TestCase):

    def test_expected_duration_of_speciation(self):
        refs = [
            ((0.01, 0.01, 0.0), 69.31471805599453),
            ((0.01, 0.01, 0.005), 60.96406423356387),
            ((0.01, 0.01, 0.01), 52.344553381571906),
            ((0.01, 0.05, 0.0), 18.232155679395458),
            ((0.01, 0.05, 0.005), 17.03243269197848),
            ((0.01, 0.05, 0.01), 15.945243534199998),
            ((0.01, 0.1, 0.0), 9.531017980432493),
            ((0.01, 0.1, 0.005), 9.151173562881866),
            ((0.01, 0.1, 0.01), 8.796646443907925),
            ((0.01, 0.5, 0.0), 1.980262729617973),
            ((0.01, 0.5, 0.005), 1.9614061471032902),
            ((0.01, 0.5, 0.01), 1.9428949825831419),
            ((0.01, 1.0, 0.0), 0.9950330853168092),
            ((0.01, 1.0, 0.005), 0.9901794650170307),
            ((0.01, 1.0, 0.01), 0.9853722709707085),
            ((0.05, 0.01, 0.0), 35.8351893845611),
            ((0.05, 0.01, 0.025), 34.610696533124994),
            ((0.05, 0.01, 0.05), 27.67607317551539),
            ((0.05, 0.05, 0.0), 13.862943611198906),
            ((0.05, 0.05, 0.025), 12.192812846712773),
            ((0.05, 0.05, 0.05), 10.46891067631438),
            ((0.05, 0.1, 0.0), 8.109302162163289),
            ((0.05, 0.1, 0.025), 7.261356156801292),
            ((0.05, 0.1, 0.05), 6.4859100952301745),
            ((0.05, 0.5, 0.0), 1.9062035960864987),
            ((0.05, 0.5, 0.025), 1.830234712576371),
            ((0.05, 0.5, 0.05), 1.7593292887815883),
            ((0.05, 1.0, 0.0), 0.9758032833886409),
            ((0.05, 1.0, 0.025), 0.9541399331844385),
            ((0.05, 1.0, 0.05), 0.9333523826096167),
            ((0.1, 0.01, 0.0), 23.978952727983707),
            ((0.1, 0.01, 0.05), 25.39055478519575),
            ((0.1, 0.01, 0.1), 20.28571191298721),
            ((0.1, 0.05, 0.0), 10.986122886681098),
            ((0.1, 0.05, 0.05), 9.802581434685468),
            ((0.1, 0.05, 0.1), 8.109302162163289),
            ((0.1, 0.1, 0.0), 6.931471805599453),
            ((0.1, 0.1, 0.05), 6.096406423356386),
            ((0.1, 0.1, 0.1), 5.23445533815719),
            ((0.1, 0.5, 0.0), 1.8232155679395459),
            ((0.1, 0.5, 0.05), 1.7032432691978485),
            ((0.1, 0.5, 0.1), 1.5945243534199998),
            ((0.1, 1.0, 0.0), 0.9531017980432493),
            ((0.1, 1.0, 0.05), 0.9151173562881855),
            ((0.1, 1.0, 0.1), 0.8796646443907942),
            ((0.5, 0.01, 0.0), 7.8636512654486515),
            ((0.5, 0.01, 0.25), 10.29910408256293),
            ((0.5, 0.01, 0.5), 9.485738357710689),
            ((0.5, 0.05, 0.0), 4.795790545596741),
            ((0.5, 0.05, 0.25), 5.07811095703915),
            ((0.5, 0.05, 0.5), 4.057142382597441),
            ((0.5, 0.1, 0.0), 3.58351893845611),
            ((0.5, 0.1, 0.25), 3.4610696533125007),
            ((0.5, 0.1, 0.5), 2.7676073175515397),
            ((0.5, 0.5, 0.0), 1.3862943611198906),
            ((0.5, 0.5, 0.25), 1.219281284671277),
            ((0.5, 0.5, 0.5), 1.046891067631438),
            ((0.5, 1.0, 0.0), 0.8109302162163288),
            ((0.5, 1.0, 0.25), 0.7261356156801294),
            ((0.5, 1.0, 0.5), 0.6485910095230177),
            ((1.0, 0.01, 0.0), 4.61512051684126),
            ((1.0, 0.01, 0.5), 6.465141234074228),
            ((1.0, 0.01, 1.0), 6.774422698356164),
            ((1.0, 0.05, 0.0), 3.044522437723423),
            ((1.0, 0.05, 0.5), 3.556910942852685),
            ((1.0, 0.05, 1.0), 2.9389333245105935),
            ((1.0, 0.1, 0.0), 2.3978952727983707),
            ((1.0, 0.1, 0.5), 2.539055478519575),
            ((1.0, 0.1, 1.0), 2.0285711912987203),
            ((1.0, 0.5, 0.0), 1.0986122886681098),
            ((1.0, 0.5, 0.5), 0.9802581434685471),
            ((1.0, 0.5, 1.0), 0.8109302162163288),
            ((1.0, 1.0, 0.0), 0.6931471805599453),
            ((1.0, 1.0, 0.5), 0.6096406423356385),
            ((1.0, 1.0, 1.0), 0.523445533815719),
        ]
        for params, exp_result in refs:
            obs_result = protractedspeciation.expected_duration_of_speciation(
                    speciation_initiation_rate=params[0],
                    speciation_completion_rate=params[1],
                    incipient_species_extinction_rate=params[2],
                    )
            self.assertAlmostEqual(obs_result, exp_result)

    def test_probability_of_duration_of_speciation(self):
        refs = [
            ((0.01, 0.01, 0.0, 0.8795348619957456), 0.009999226458320016),
            ((0.01, 0.01, 0.0, 0.14519027181737457), 0.009999978919814594),
            ((0.01, 0.01, 0.0, 34.684466460885076), 0.008887282087181123),
            ((0.01, 0.01, 0.0, 63.20061216231189), 0.0068703579303570775),
            ((0.01, 0.01, 0.005, 0.08536631008065432), 0.012802289144553825),
            ((0.01, 0.01, 0.005, 0.011137362557058923), 0.01280705070148593),
            ((0.01, 0.01, 0.005, 76.52294782192368), 0.005405202225298228),
            ((0.01, 0.01, 0.005, 72.84339288579332), 0.0057297611990511305),
            ((0.01, 0.01, 0.01, 0.33558726166937963), 0.016125950330407143),
            ((0.01, 0.01, 0.01, 0.20457203115105108), 0.016147205741964963),
            ((0.01, 0.01, 0.01, 19.079554492974783), 0.012922829672422163),
            ((0.01, 0.01, 0.01, 24.68247362845256), 0.011956431293373187),
            ((0.01, 0.1, 0.0, 0.5935352032024056), 0.09476554969639099),
            ((0.01, 0.1, 0.0, 0.7003470486083184), 0.09384632919173516),
            ((0.01, 0.1, 0.0, 81.78775148685204), 1.498228340113515e-05),
            ((0.01, 0.1, 0.0, 23.412278765465388), 0.009072720079842342),
            ((0.01, 0.1, 0.005, 0.7792322482865833), 0.09704508143428221),
            ((0.01, 0.1, 0.005, 0.963289936603788), 0.09533361321027416),
            ((0.01, 0.1, 0.005, 22.853616783069036), 0.00905378799214888),
            ((0.01, 0.1, 0.005, 65.4548715141025), 7.09734365226258e-05),
            ((0.01, 0.1, 0.01, 0.7004511480083776), 0.10172748819395837),
            ((0.01, 0.1, 0.01, 0.3859718489988695), 0.10501232743022515),
            ((0.01, 0.1, 0.01, 18.979988804020195), 0.013337248459689992),
            ((0.01, 0.1, 0.01, 22.757722164422592), 0.0085842630673195),
            ((0.01, 1.0, 0.0, 0.3130295270193749), 0.7328713233872153),
            ((0.01, 1.0, 0.0, 0.31154063743324156), 0.7339582420670433),
            ((0.01, 1.0, 0.0, 88.07424235690094), 2.3767511877766183e-39),
            ((0.01, 1.0, 0.0, 50.37034301380811), 8.209658743393956e-23),
            ((0.01, 1.0, 0.005, 0.9003232000282951), 0.4077827698121026),
            ((0.01, 1.0, 0.005, 0.12159170357864164), 0.8903084715946127),
            ((0.01, 1.0, 0.005, 81.55525905558959), 1.1585339262652681e-36),
            ((0.01, 1.0, 0.005, 72.67637856869439), 9.493325334690782e-33),
            ((0.01, 1.0, 0.01, 0.6890551804405882), 0.5050883931592011),
            ((0.01, 1.0, 0.01, 0.9501351410838736), 0.3878984254146435),
            ((0.01, 1.0, 0.01, 20.649073041437447), 7.368492351424008e-10),
            ((0.01, 1.0, 0.01, 40.24647761617133), 1.540927018047582e-18),
            ((0.1, 0.01, 0.0, 0.8113389029985384), 0.01075027226017838),
            ((0.1, 0.01, 0.0, 0.846887182506354), 0.010784061294478495),
            ((0.1, 0.01, 0.0, 59.62256474016563), 0.0016681779606101286),
            ((0.1, 0.01, 0.0, 91.65272435084043), 5.057688899425374e-05),
            ((0.1, 0.01, 0.05, 0.8760733456577112), 0.018023721654607777),
            ((0.1, 0.01, 0.05, 0.8409206587355896), 0.01799950113041887),
            ((0.1, 0.01, 0.05, 51.06343770607376), 0.006125364048638298),
            ((0.1, 0.01, 0.05, 21.687465805666022), 0.023289784397034694),
            ((0.1, 0.01, 0.1, 0.4792155173312089), 0.03683021581716226),
            ((0.1, 0.01, 0.1, 0.7263126839962858), 0.036728414105188914),
            ((0.1, 0.01, 0.1, 55.892356838488816), 0.002968979103551312),
            ((0.1, 0.01, 0.1, 81.82838924315757), 0.0005828087031488993),
            ((0.1, 0.1, 0.0, 0.8258322196685367), 0.09932109002004828),
            ((0.1, 0.1, 0.0, 0.2824909284607273), 0.09992024131094354),
            ((0.1, 0.1, 0.0, 53.32752907252014), 9.334034463585915e-06),
            ((0.1, 0.1, 0.0, 67.44160225894768), 5.548245175549618e-07),
            ((0.1, 0.1, 0.05, 0.5636100846520368), 0.12412769734346397),
            ((0.1, 0.1, 0.05, 0.9787149640606224), 0.12081870252730309),
            ((0.1, 0.1, 0.05, 55.67386202995424), 3.438050746152094e-06),
            ((0.1, 0.1, 0.05, 43.70942242490415), 4.049995603202447e-05),
            ((0.1, 0.1, 0.1, 0.3999889157738908), 0.1552139893311176),
            ((0.1, 0.1, 0.1, 0.6345726544794333), 0.15125759800964952),
            ((0.1, 0.1, 0.1, 97.62739362247933), 1.0215724557144834e-10),
            ((0.1, 0.1, 0.1, 24.52730240810091), 0.0012785880221897646),
            ((0.1, 1.0, 0.0, 0.8115065427563289), 0.4573450980393226),
            ((0.1, 1.0, 0.0, 0.4516132153221019), 0.6542319875705164),
            ((0.1, 1.0, 0.0, 40.872089582281994), 3.6075312673974335e-20),
            ((0.1, 1.0, 0.0, 24.255563803233162), 3.1283987494891642e-12),
            ((0.1, 1.0, 0.05, 0.3457237938404357), 0.7448748845936354),
            ((0.1, 1.0, 0.05, 0.576803770523732), 0.5879219371779756),
            ((0.1, 1.0, 0.05, 83.20031921894818), 7.200169398604266e-42),
            ((0.1, 1.0, 0.05, 91.43820228760569), 5.9464003697653356e-46),
            ((0.1, 1.0, 0.1, 0.031986829095719926), 1.057136415593703),
            ((0.1, 1.0, 0.1, 0.5519036055923441), 0.6128084744263474),
            ((0.1, 1.0, 0.1, 77.67899885507961), 1.554561011597814e-40),
            ((0.1, 1.0, 0.1, 20.09264881920249), 6.069727719207903e-11),
            ((1.0, 0.01, 0.0, 0.4446353921382625), 0.015494312486713001),
            ((1.0, 0.01, 0.0, 0.06884594296267747), 0.01070481974865263),
            ((1.0, 0.01, 0.0, 45.0622683427891), 1.7484182033029857e-18),
            ((1.0, 0.01, 0.0, 23.93500942694461), 3.2348325257220874e-09),
            ((1.0, 0.01, 0.5, 0.07678647673830719), 0.020373520291909725),
            ((1.0, 0.01, 0.5, 0.3885919433134037), 0.02369988607726227),
            ((1.0, 0.01, 0.5, 45.11804274928912), 6.08508209809671e-10),
            ((1.0, 0.01, 0.5, 20.568063200782436), 0.0002670604419545594),
            ((1.0, 0.01, 1.0, 0.224072444141194), 0.10483702314125201),
            ((1.0, 0.01, 1.0, 0.5226436349919803), 0.1042922825095739),
            ((1.0, 0.01, 1.0, 33.02193432605858), 0.0005111539346059228),
            ((1.0, 0.01, 1.0, 40.74024255784781), 0.00010917936840225393),
            ((1.0, 0.1, 0.0, 0.8362413301387851), 0.194015289449898),
            ((1.0, 0.1, 0.0, 0.6926320472771594), 0.17581989725130096),
            ((1.0, 0.1, 0.0, 92.62429483753988), 6.822494879649729e-44),
            ((1.0, 0.1, 0.0, 19.625562941290408), 5.095470692858589e-09),
            ((1.0, 0.1, 0.5, 0.5302925473060789), 0.20894273300712232),
            ((1.0, 0.1, 0.5, 0.041980900523581485), 0.17708370000687484),
            ((1.0, 0.1, 0.5, 52.81106500217504), 2.2071314894480564e-17),
            ((1.0, 0.1, 0.5, 25.460992608253292), 1.7080423876427843e-08),
            ((1.0, 0.1, 1.0, 0.9296668365468552), 0.31052945846276236),
            ((1.0, 0.1, 1.0, 0.12145631266764487), 0.36515087675508107),
            ((1.0, 0.1, 1.0, 50.374060615431524), 1.0868531521448435e-14),
            ((1.0, 0.1, 1.0, 11.579762731930778), 0.0006666049628160538),
            ((1.0, 1.0, 0.0, 0.09683251662922314), 0.9906817667603305),
            ((1.0, 1.0, 0.0, 0.10193348600169332), 0.989681117037441),
            ((1.0, 1.0, 0.0, 39.66911063337399), 1.3992918064758832e-34),
            ((1.0, 1.0, 0.0, 13.915328179066181), 3.2761177836735036e-12),
            ((1.0, 1.0, 0.5, 0.8635504500484475), 0.46002864196564114),
            ((1.0, 1.0, 0.5, 0.7909190840591199), 0.518603031188729),
            ((1.0, 1.0, 0.5, 86.89921906556235), 5.226596384633828e-78),
            ((1.0, 1.0, 0.5, 21.085927006318414), 4.3877525911023785e-19),
            ((1.0, 1.0, 1.0, 0.6306393480885447), 0.6311458539888442),
            ((1.0, 1.0, 1.0, 0.041694658805305165), 1.5492998515367982),
            ((1.0, 1.0, 1.0, 54.275592049617586), 6.057654470628204e-53),
            ((1.0, 1.0, 1.0, 79.41886144450676), 2.3193926494535927e-77),
        ]
        for params, exp_result in refs:
            obs_result = protractedspeciation.probability_of_duration_of_speciation(
                    tau=params[3],
                    speciation_initiation_rate=params[0],
                    speciation_completion_rate=params[1],
                    incipient_species_extinction_rate=params[2],
                    )
            self.assertAlmostEqual(obs_result, exp_result)
            obs_result = protractedspeciation.log_probability_of_duration_of_speciation(
                    tau=params[3],
                    speciation_initiation_rate=params[0],
                    speciation_completion_rate=params[1],
                    incipient_species_extinction_rate=params[2],
                    )
            self.assertAlmostEqual(obs_result, math.log(exp_result))

    def test_maximum_probability_duration_of_speciation(self):
        refs = [
            ((0.01, 0.0001, 0.0), 455.9574441572374),
            ((0.01, 0.0001, 0.005), 615.4051920172764),
            ((0.01, 0.0001, 0.01), 0.0),
            ((0.01, 0.001, 0.0), 209.32591754491338),
            ((0.01, 0.001, 0.005), 159.4094017279303),
            ((0.01, 0.001, 0.01), 0.0),
            ((0.01, 0.01, 0.0), 0.0),
            ((0.01, 0.01, 0.005), 0.0),
            ((0.01, 0.01, 0.01), 0.0),
            ((0.01, 100.0, 0.0), 0.0),
            ((0.01, 100.0, 0.005), 0.0),
            ((0.01, 100.0, 0.01), 0.0),
            ((0.5, 0.005, 0.0), 9.119148883144735),
            ((0.5, 0.005, 0.25), 12.308103840345511),
            ((0.5, 0.005, 0.5), 0.0),
            ((0.5, 0.05, 0.0), 4.186518350898265),
            ((0.5, 0.05, 0.25), 3.1881880345586024),
            ((0.5, 0.05, 0.5), 0.0),
            ((0.5, 0.5, 0.0), 0.0),
            ((0.5, 0.5, 0.25), 0.0),
            ((0.5, 0.5, 0.5), 0.0),
            ((0.5, 5000.0, 0.0), 0.0),
            ((0.5, 5000.0, 0.25), 0.0),
            ((0.5, 5000.0, 0.5), 0.0),
            ((1.0, 0.01, 0.0), 4.559574441572368),
            ((1.0, 0.01, 0.5), 6.154051920172756),
            ((1.0, 0.01, 1.0), 0.0),
            ((1.0, 0.1, 0.0), 2.0932591754491323),
            ((1.0, 0.1, 0.5), 1.5940940172793012),
            ((1.0, 0.1, 1.0), 0.0),
            ((1.0, 1.0, 0.0), 0.0),
            ((1.0, 1.0, 0.5), 0.0),
            ((1.0, 1.0, 1.0), 0.0),
            ((1.0, 10000.0, 0.0), 0.0),
            ((1.0, 10000.0, 0.5), 0.0),
            ((1.0, 10000.0, 1.0), 0.0),
            ((100, 1.0, 0.0), 0.045595744415723685),
            ((100, 1.0, 50.0), 0.06154051920172753),
            ((100, 1.0, 100), 0.0),
            ((100, 10.0, 0.0), 0.020932591754491324),
            ((100, 10.0, 50.0), 0.01594094017279301),
            ((100, 10.0, 100), 0.0),
            ((100, 100.0, 0.0), 0.0),
            ((100, 100.0, 50.0), 0.0),
            ((100, 100.0, 100), 0.0),
            ((100, 1000000, 0.0), 0.0),
            ((100, 1000000, 50.0), 0.0),
            ((100, 1000000, 100), 0.0),
            ((1000, 10.0, 0.0), 0.004559574441572368),
            ((1000, 10.0, 500.0), 0.006154051920172753),
            ((1000, 10.0, 1000), 0.0),
            ((1000, 100.0, 0.0), 0.0020932591754491327),
            ((1000, 100.0, 500.0), 0.0015940940172793014),
            ((1000, 100.0, 1000), 0.0),
            ((1000, 1000.0, 0.0), 0.0),
            ((1000, 1000.0, 500.0), 0.0),
            ((1000, 1000.0, 1000), 0.0),
            ((1000, 10000000, 0.0), 0.0),
            ((1000, 10000000, 500.0), 0.0),
            ((1000, 10000000, 1000), 0.0),
        ]
        for params, exp_result in refs:
            obs_result = protractedspeciation.maximum_probability_duration_of_speciation(
                    speciation_initiation_rate=params[0],
                    speciation_completion_rate=params[1],
                    incipient_species_extinction_rate=params[2],
                    )
            self.assertAlmostEqual(obs_result, exp_result)

class LineageQueueTestCase(unittest.TestCase):

    def setUp(self):
        self.pbd_tables = (
                    (
                        [1  , 0  , -1e-10             , 0                , -1               , 1] ,
                        [2  , 1  , 6.67615606179445   , -1               , -1               , 1] ,
                        [3  , 1  , 0                  , -1               , -1               , 1] ,
                        [4  , -3 , 5.03236379137566   , -1               , 12.4555382563192 , 1] ,
                        [5  , -4 , 5.82067825984527   , 8.46863636938638 , 10.0921769421767 , 2] ,
                        [6  , -4 , 7.43692788970062   , -1               , 9.49403765538272 , 1] ,
                        [7  , -3 , 9.3964900131076    , -1               , -1               , 1] ,
                        [8  , -4 , 12.0831255638222   , -1               , -1               , 1] ,
                        [9  , -7 , 13.2396812008651   , -1               , -1               , 1] ,
                        [10 , -8 , 13.437755189092    , 13.7529521422405 , -1               , 3]
                    ),
                    (
                        (1  , 0   , -1e-10            , 0                , -1 , 1) ,
                        (2  , 1   , 2.85912504426404  , 3.34543532862682 , -1 , 2) ,
                        (3  , 1   , 0                 , -1               , -1 , 1) ,
                        (4  , -3  , 0.772044112888897 , 9.58420144285785 , -1 , 4) ,
                        (5  , -3  , 5.17020330746366  , -1               , -1 , 1) ,
                        (6  , -3  , 6.91843597758094  , -1               , -1 , 1) ,
                        (7  , -6  , 8.0848153631966   , 9.09087035090051 , -1 , 3) ,
                        (8  , -5  , 8.17426394667072  , -1               , -1 , 1) ,
                        (9  , 4   , 9.67145552103779  , 11.527621978474  , -1 , 5) ,
                        (10 , -8  , 10.1553007690956  , -1               , -1 , 1) ,
                        (11 , -6  , 10.9176536767969  , -1               , -1 , 1) ,
                        (12 , -8  , 12.3660726549949  , 13.6104437800854 , -1 , 6) ,
                        (13 , -3  , 12.5586831168778  , -1               , -1 , 1) ,
                        (14 , -8  , 12.6877684466107  , -1               , -1 , 1) ,
                        (15 , -12 , 13.1641235715743  , -1               , -1 , 1) ,
                        (16 , -10 , 14.5891386243993  , -1               , -1 , 1) ,
                        (17 , -6  , 14.8032830504263  , -1               , -1 , 1)
                    ),
                )

    def test_build_from_collection(self):
        psm = _get_stock_psm()
        for tree_type in ("lineage", "species"):
            for is_drop_extinct in (False, True):
                if tree_type == "lineage" and not is_drop_extinct:
                    continue
                if tree_type == "lineage":
                    label_template_attr = "lineage_label_format_template"
                    node_attr = "lineage_node"
                elif tree_type == "species":
                    label_template_attr = "species_label_format_template"
                    node_attr = "species_node"
                for pbd_table in self.pbd_tables:
                    lineage_collection = [protractedspeciation.ProtractedSpeciationProcess._Lineage.from_pbd_entry(pbd_lineage_entry)
                            for pbd_lineage_entry in pbd_table]
                    lineage_queue = psm._build_lineage_queue(
                            lineage_collection,
                            max_time=100,
                            is_drop_extinct=is_drop_extinct,
                            node_attr="lineage_node",
                            label_template = getattr(psm, label_template_attr),
                            )
                    if is_drop_extinct:
                        pbd_table = [entry for entry in pbd_table if entry[_ET_EXTINCTION_TIME] < 0]
                    expected = sorted(pbd_table, key=lambda x: x[2], reverse=True)
                    self.assertEqual(len(lineage_queue), len(pbd_table))
                    for ei in expected:
                        lineage_entry = lineage_queue.pop_youngest_lineage()
                        self.assertEqual(ei[0], lineage_entry.lineage_id)
                        self.assertEqual(abs(ei[1]), lineage_entry.parent_lineage_id)
                        if ei[1] < 0:
                            self.assertFalse(lineage_entry.is_parent_orthospecies)
                        else:
                            self.assertTrue(lineage_entry.is_parent_orthospecies)
                        self.assertEqual(ei[2], lineage_entry.origin_time)
                        self.assertEqual(ei[3], lineage_entry.speciation_completion_time)
                        if ei[4] < 0:
                            self.assertIs(None, lineage_entry.extinction_time)
                            self.assertFalse(lineage_entry.is_extinct)
                        else:
                            self.assertEqual(ei[4], lineage_entry.extinction_time)
                            self.assertTrue(lineage_entry.is_extinct)
                        self.assertEqual(ei[5], lineage_entry.species_id)

class ProtractedSpeciationLowLevelTreeCompilationFromEventsTestCase(unittest.TestCase):
    test_reference = _load_json("protracted_speciation_process.json")

    def setUp(self):
        ## dummy PSM; use when not simulating process but checking other aspects
        self.psm = _get_stock_psm()

    def pbd_table_to_lineage_collection(self, pbd_table):
        lineage_collection = [protractedspeciation.ProtractedSpeciationProcess._Lineage.from_pbd_entry(pbd_lineage_entry)
                    for pbd_lineage_entry in pbd_table]
        return lineage_collection

    def assert_equal_trees(self, t0, t1):
        self.assertEqual(treecompare.unweighted_robinson_foulds_distance(t0, t1), 0)
        self.assertAlmostEqual(treecompare.weighted_robinson_foulds_distance(t0, t1), 0, 8)

    def iter_test_references(self):
        for tidx, test_ref in enumerate(self.test_reference):
            r = {}
            r["data"] = test_ref
            r["psm"] = protractedspeciation.ProtractedSpeciationProcess(
                    speciation_initiation_from_orthospecies_rate=test_ref["lineage_origination_rate"],
                    speciation_initiation_from_incipient_species_rate=test_ref["lineage_origination_rate"],
                    speciation_completion_rate=test_ref["speciation_completion_rate"],
                    orthospecies_extinction_rate=test_ref["extinction_rate"],
                    incipient_species_extinction_rate=test_ref["extinction_rate"],
                    )
            r["lineage_tree"] = dendropy.Tree.get(
                    data=test_ref["lineage_tree"],
                    schema="newick",
                    rooting="force-rooted",
                    )
            r["lineage_tree_incl_extinct"] = dendropy.Tree.get(
                    data=test_ref["lineage_tree_incl_extinct"],
                    schema="newick",
                    rooting="force-rooted",
                    )
            r["species_tree_oldest_samples"] = dendropy.Tree.get(
                    data=test_ref["species_tree_oldest_samples"],
                    schema="newick",
                    rooting="force-rooted",
                    )
            r["species_tree_youngest_samples"] = dendropy.Tree.get(
                    data=test_ref["species_tree_youngest_samples"],
                    schema="newick",
                    rooting="force-rooted",
                    )
            r["lineage_collection"] = self.pbd_table_to_lineage_collection(test_ref["lineage_table"])
            yield r

    def test_lineage_tree_compilation(self):
        for test_idx, test_ref in enumerate(self.iter_test_references()):
            lineage_collection = test_ref["lineage_collection"]
            t0 = test_ref["lineage_tree"]
            t1 = self.psm._compile_lineage_tree(
                    lineage_collection=lineage_collection,
                    max_time=test_ref["data"]["age"],
                    is_drop_extinct=True)
            self.psm._build_taxa(tree=t1, taxon_namespace=t0.taxon_namespace)
            self.assert_equal_trees(t0, t1)

    def test_lineage_tree_incl_extinct_compilation(self):
        for test_ref in self.iter_test_references():
            lineage_collection = test_ref["lineage_collection"]
            t0 = test_ref["lineage_tree_incl_extinct"]
            t1 = self.psm._compile_lineage_tree(
                    lineage_collection=lineage_collection,
                    max_time=test_ref["data"]["age"],
                    is_drop_extinct=False,
                    )
            self.psm._build_taxa(tree=t1, taxon_namespace=t0.taxon_namespace)
            self.assert_equal_trees(t0, t1)

    def test_species_tree_of_oldest_lineages(self):
        for test_idx, test_ref in enumerate(self.iter_test_references()):
            lineage_collection = test_ref["lineage_collection"]
            t0 = test_ref["species_tree_oldest_samples"]
            self.psm.species_lineage_sampling_scheme = "oldest"
            t1 = self.psm._compile_species_tree(
                    lineage_collection=lineage_collection,
                    max_time=test_ref["data"]["age"],
                    )
            self.psm._build_taxa(tree=t1, taxon_namespace=t0.taxon_namespace)
            self.assert_equal_trees(t0, t1)

    def test_species_tree_of_youngest_lineages(self):
        for test_idx, test_ref in enumerate(self.iter_test_references()):
            lineage_collection = test_ref["lineage_collection"]
            t0 = test_ref["species_tree_youngest_samples"]
            self.psm.species_lineage_sampling_scheme = "youngest"
            t1 = self.psm._compile_species_tree(
                    lineage_collection=lineage_collection,
                    max_time=test_ref["data"]["age"],
                    )
            self.psm._build_taxa(tree=t1, taxon_namespace=t0.taxon_namespace)
            self.assert_equal_trees(t0, t1)

    def test_multi_tree_compilation(self):
        for test_idx, test_ref in enumerate(self.iter_test_references()):
            lineage_collection = test_ref["lineage_collection"]
            t1_0 = test_ref["lineage_tree"]
            t1_1 = self.psm._compile_lineage_tree(
                    lineage_collection=lineage_collection,
                    max_time=test_ref["data"]["age"],
                    is_drop_extinct=True,
                    )
            self.psm._build_taxa(tree=t1_1, taxon_namespace=t1_0.taxon_namespace)
            self.assert_equal_trees(t1_0, t1_1)
            t2_0 = test_ref["lineage_tree_incl_extinct"]
            t2_1 = self.psm._compile_lineage_tree(
                    lineage_collection=lineage_collection,
                    max_time=test_ref["data"]["age"],
                    is_drop_extinct=False,
                    )
            self.psm._build_taxa(tree=t2_1, taxon_namespace=t2_0.taxon_namespace)
            self.assert_equal_trees(t2_0, t2_1)
            self.psm.species_lineage_sampling_scheme = "oldest"
            t3_0 = test_ref["species_tree_oldest_samples"]
            t3_1 = self.psm._compile_species_tree(
                    lineage_collection=lineage_collection,
                    max_time=test_ref["data"]["age"],
                    )
            self.psm._build_taxa(tree=t3_1, taxon_namespace=t3_0.taxon_namespace)
            self.assert_equal_trees(t3_0, t3_1)
            self.psm.species_lineage_sampling_scheme = "youngest"
            t4_0 = test_ref["species_tree_youngest_samples"]
            t4_1 = self.psm._compile_species_tree(
                    lineage_collection=lineage_collection,
                    max_time=test_ref["data"]["age"])
            self.psm._build_taxa(tree=t4_1, taxon_namespace=t4_0.taxon_namespace)
            self.assert_equal_trees(t4_0, t4_1)

class ProtractedSpeciationProcessGeneration(unittest.TestCase):

    def iter_psm_models(self, **kwargs):
        for splitting_rate in (0.1,):
            for extinction_rate_factor in (0.5, 0.0):
                extinction_rate = splitting_rate * extinction_rate_factor
                for speciation_completion_rate in (splitting_rate * 0.5,):
                    psm = protractedspeciation.ProtractedSpeciationProcess(
                            speciation_initiation_from_orthospecies_rate=splitting_rate,
                            orthospecies_extinction_rate=extinction_rate,
                            speciation_initiation_from_incipient_species_rate=splitting_rate,
                            speciation_completion_rate=speciation_completion_rate,
                            incipient_species_extinction_rate=extinction_rate,
                            **kwargs
                            )
                    yield psm

    def iter_samples(self, psm, additional_kwargs=None):
        if additional_kwargs is None:
            additional_kwargs = {}
        for kwargs in (
                {"max_time": 20},
                {"num_extant_orthospecies": 10},
                {"num_extant_lineages": 20},
                ):
            for is_initial_lineage_orthospecies in (True, False):
                psm.is_initial_lineage_orthospecies = is_initial_lineage_orthospecies
                kwargs.update(additional_kwargs)
                yield psm.generate_sample(**kwargs)

    def check(self, tree):
        tree.calc_node_ages()
        leaf_root_distances = tree.calc_node_root_distances(return_leaf_distances_only=True)
        for dist in leaf_root_distances:
            self.assertAlmostEqual(dist, tree.seed_node.age, 8)
        seen_taxa = set()
        seen_taxon_labels = set()
        expected_taxa = set([taxon for taxon in tree.taxon_namespace])
        num_leaves = 0
        for leaf in tree.leaf_node_iter():
            self.assertIn(leaf.taxon, expected_taxa)
            expected_taxa.remove(leaf.taxon)
            self.assertNotIn(leaf.taxon, seen_taxa)
            seen_taxa.add(leaf.taxon)
            self.assertNotIn(leaf.taxon.label, seen_taxon_labels)
            seen_taxon_labels.add(leaf.taxon.label)
            num_leaves += 1
        if len(expected_taxa) != 0:
            print(tree.as_string("newick"))
            print("Remaining: {}".format(", ".join(t.label for t in expected_taxa)))
        self.assertEqual(len(expected_taxa), 0)
        self.assertEqual(num_leaves, len(tree.taxon_namespace))
        for nd in tree.internal_nodes():
            self.assertEqual(len(nd.child_nodes()), 2)

    def test_by_num_time(self):
        for psm in self.iter_psm_models():
            for max_time in (10, 15, 20):
                lineage_tree, orthospecies_tree = psm.generate_sample(max_time=max_time)
                for tree_idx, tree in enumerate((lineage_tree, orthospecies_tree,)):
                    for nd in tree.leaf_node_iter():
                        self.assertAlmostEqual(nd._time, max_time, 8)
                    self.check(tree)

    def test_by_num_lineages(self):
        for psm in self.iter_psm_models():
            for num in (5, 25, 50):
                lineage_tree, orthospecies_tree = psm.generate_sample(num_extant_lineages=num)
                # self.assertEqual(len(lineage_tree.taxon_namespace), num)
                self.assertEqual(len(lineage_tree.taxon_namespace), num)
                for tree_idx, tree in enumerate((lineage_tree, orthospecies_tree,)):
                    self.check(tree)

    def test_by_num_orthospecies(self):
        for psm in self.iter_psm_models():
            for num in (5, 10, 20):
                # lineage_tree, orthospecies_tree = psm.generate_sample(num_extant_orthospecies=num)
                lineage_tree, orthospecies_tree = psm.generate_sample(num_extant_orthospecies=num)
                self.assertEqual(len(orthospecies_tree.taxon_namespace), num)
                for tree_idx, tree in enumerate((lineage_tree, orthospecies_tree,)):
                    self.check(tree)

    def test_taxon_assignment_and_namespace(self):
        for seed in itertools.chain((559, 631, 230, 212, 907, 237,), (random.randint(0, 1000) for i in range(10))):
            rng = random.Random(seed)
            for psm in self.iter_psm_models(rng=rng):
                for kwargs in (
                        {"max_time": 20},
                        {"num_extant_orthospecies": 10},
                        {"num_extant_lineages": 20},
                        ):
                    lineage_taxon_namespace = dendropy.TaxonNamespace()
                    species_taxon_namespace = dendropy.TaxonNamespace()
                    kwargs["lineage_taxon_namespace"] = lineage_taxon_namespace
                    kwargs["species_taxon_namespace"] = species_taxon_namespace
                    lineage_tree, orthospecies_tree = psm.generate_sample(**kwargs)
                    self.assertIs(lineage_tree.taxon_namespace, lineage_taxon_namespace)
                    self.assertIs(orthospecies_tree.taxon_namespace, species_taxon_namespace)
                    for tree in (lineage_tree, orthospecies_tree):
                        self.check(tree)

    def test_lineage_species_tip_correlation(self):
        psm = protractedspeciation.ProtractedSpeciationProcess(
                speciation_initiation_from_orthospecies_rate=0.1,
                speciation_initiation_from_incipient_species_rate=0.1,
                speciation_completion_rate=0.05,
                orthospecies_extinction_rate=0.0,
                incipient_species_extinction_rate=0.00,
                )
        lineage_tree, orthospecies_tree = psm.generate_sample(num_extant_orthospecies=5)
        # for seed in itertools.chain((559, 631, 230, 212, 907, 237,), (random.randint(0, 1000) for i in range(10))):
        for seed in itertools.chain((559, 631, 230, 212, 907, 237,), (random.randint(0, 1000) for i in range(20))):
            rng = random.Random(seed)
            for psm in self.iter_psm_models(rng=rng):
                for test_idx, (lineage_tree, orthospecies_tree) in enumerate(self.iter_samples(psm)):
                    for tree in (lineage_tree, orthospecies_tree):
                        self.check(tree)
                    lineage_tree_species_labels = set([taxon.label.split(".")[0] for taxon in lineage_tree.taxon_namespace])
                    species_tree_taxon_labels = set([taxon.label for taxon in orthospecies_tree.taxon_namespace])
                    self.assertEqual(lineage_tree_species_labels, species_tree_taxon_labels)
                    # print("\nLineages: {}\nSpecies: {}".format(",".join(lineage_tree_species_labels), ", ".join(species_tree_taxon_labels)))
                    check_species_node_lineage_nodes_map = {}
                    for lineage_node in lineage_tree.leaf_node_iter():
                        lineage_node_species_node = getattr(lineage_node, psm.lineage_tree_to_species_tree_node_attr)
                        species_node_lineage_nodes = getattr(lineage_node_species_node, psm.species_tree_to_lineage_tree_node_attr)
                        self.assertIn(lineage_node, species_node_lineage_nodes)
                        self.assertEqual(lineage_node.taxon.label.split(".")[0], lineage_node_species_node.taxon.label)
                        try:
                            check_species_node_lineage_nodes_map[lineage_node_species_node].add(lineage_node)
                        except KeyError:
                            check_species_node_lineage_nodes_map[lineage_node_species_node] = set([lineage_node])
                    for species_node in orthospecies_tree.leaf_node_iter():
                        species_node_lineage_nodes = getattr(species_node, psm.species_tree_to_lineage_tree_node_attr)
                        self.assertEqual(check_species_node_lineage_nodes_map[species_node], species_node_lineage_nodes)

    def test_(self):
        for psm in self.iter_psm_models():
            for test_idx, (lineage_tree, orthospecies_tree) in enumerate(self.iter_samples(psm)):
                pass

if __name__ == "__main__":
    unittest.main()
