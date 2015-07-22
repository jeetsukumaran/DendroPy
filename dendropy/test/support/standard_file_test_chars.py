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

import collections
from dendropy.datamodel import charstatemodel
from dendropy.datamodel import charmatrixmodel

def general_verify_taxa(
        test_case,
        char_matrix,
        checker_reference_class,
        taxon_namespace,
        check_annotations=True):
    test_case.assertEqual(len(taxon_namespace), len(checker_reference_class.labels))
    for taxon, label in zip(taxon_namespace, checker_reference_class.labels):
        test_case.assertEqual(taxon.label, label)

def general_verify_state_alphabet_symbols(
        test_case,
        state_alphabet,
        checker_reference_class):
    fundamental_symbols = list(state_alphabet.fundamental_symbol_iter())
    test_case.assertEqual(len(fundamental_symbols), len(checker_reference_class.state_alphabet_fundamental_symbols))
    for s1, s2 in zip(fundamental_symbols, checker_reference_class.state_alphabet_fundamental_symbols):
        test_case.assertEqual(s1, s2)

def general_verify_character_cell_states(
        test_case,
        char_matrix,
        checker_reference_class,
        c1, c2,
        check_annotations=True):
    if char_matrix.data_type != "continuous":
        test_case.assertIs(c1, c2)
    else:
        test_case.assertEqual(c1, c2)
    if check_annotations:
        pass # not yet implemented

def general_verify_sequences(
        test_case,
        char_matrix,
        checker_reference_class,
        s1, s2,
        check_sequence_annotations=True,
        check_cell_annotations=True):
    # x1 = "".join(str(c) for c in s1)
    # x2 = "".join(str(c) for c in s2)
    # print("{}\n{}".format(x1, x2))
    test_case.assertEqual(len(s1), len(s2))
    idx = 0
    for c1, c2 in zip(s1, s2):
        idx += 1
        general_verify_character_cell_states(
                test_case=test_case,
                char_matrix=char_matrix,
                checker_reference_class=checker_reference_class,
                c1=c1,
                c2=c2,
                check_annotations=check_cell_annotations)

def general_char_matrix_checker(
        test_case,
        char_matrix,
        checker_reference_class,
        check_taxon_annotations=True,
        check_matrix_annotations=True,
        check_sequence_annotations=True,
        check_column_annotations=True,
        check_cell_annotations=True):
    test_case.assertEqual(len(char_matrix), len(checker_reference_class.labels))
    general_verify_taxa(
        test_case=test_case,
        char_matrix=char_matrix,
        checker_reference_class=checker_reference_class,
        taxon_namespace=char_matrix.taxon_namespace,
        check_annotations=check_taxon_annotations)
    for taxon, label in zip(char_matrix, checker_reference_class.labels):
        test_case.assertEqual(taxon.label, label)
        s1 = char_matrix[taxon]
        s2 = checker_reference_class.label_sequence_map[label]
        general_verify_sequences(
            test_case=test_case,
            char_matrix=char_matrix,
            checker_reference_class=checker_reference_class,
            s1=s1,
            s2=s2,
            check_sequence_annotations=check_sequence_annotations,
            check_cell_annotations=check_cell_annotations)

class CharacterTestChecker(object):

    @staticmethod
    def create_class_fixtures(cls, matrix_type, states_lists, labels=None):
        cls.matrix_type = matrix_type
        cls.create_class_fixtures_labels(cls, labels=labels)
        cls.create_class_fixtures_label_sequence_map(cls, states_lists=states_lists)

    @staticmethod
    def create_class_fixtures_label_sequence_map(cls, states_lists):
        assert len(cls.labels) == len(states_lists)
        cls.label_sequence_map = collections.OrderedDict()
        for label, ss in zip(cls.labels, states_lists):
            cls.label_sequence_map[label] = ss

    @staticmethod
    def create_class_fixtures_labels(cls, labels=None):
        if labels is None:
            cls.labels = (
                    "a",
                    "b",
                    "c",
                    "e",
                    "f",
                    "g",
                    "h",
                    "i",
                    "j",
                    "k",
                    "l",
                    "m",
                    "n",
                    "o",
                    "p",
            )
        else:
            cls.labels = tuple(labels)

    @classmethod
    def get_char_matrix_from_class_data(cls, taxon_namespace=None):
        c = cls.matrix_type.from_dict(source_dict=cls.label_sequence_map,
                taxon_namespace=taxon_namespace)
        return c

    def verify_char_matrix(self,
            char_matrix,
            check_taxon_annotations=True,
            check_matrix_annotations=True,
            check_sequence_annotations=True,
            check_column_annotations=True,
            check_cell_annotations=True):
        general_char_matrix_checker(self,
                char_matrix,
                self.__class__,
                check_taxon_annotations=check_taxon_annotations,
                check_matrix_annotations=check_matrix_annotations,
                check_sequence_annotations=check_sequence_annotations,
                check_column_annotations=check_column_annotations,
                check_cell_annotations=check_cell_annotations,)

    def verify_get_from(self,
            matrix_type,
            src_filepath,
            schema,
            factory_kwargs,
            check_taxon_annotations=True,
            check_matrix_annotations=True,
            check_sequence_annotations=True,
            check_column_annotations=True,
            check_cell_annotations=True):
        char_matrix = matrix_type.get(
                path=src_filepath,
                schema=schema,
                **factory_kwargs)
        self.verify_char_matrix(char_matrix,
            check_taxon_annotations=check_taxon_annotations,
            check_matrix_annotations=check_matrix_annotations,
            check_sequence_annotations=check_sequence_annotations,
            check_column_annotations=check_column_annotations,
            check_cell_annotations=check_cell_annotations)

class GenericDiscreteCharacterTestChecker(CharacterTestChecker):

    @staticmethod
    def create_class_fixtures(cls,
            matrix_type,
            state_alphabet_fundamental_symbols,
            seq_symbols,
            labels=None):
        cls.matrix_type = matrix_type
        cls.state_alphabet_fundamental_symbols = list(state_alphabet_fundamental_symbols)
        CharacterTestChecker.create_class_fixtures_labels(cls, labels=labels)
        cls.states_symbols_lists = []
        for ss in seq_symbols:
            cls.states_symbols_lists.append(list(ss))

    @staticmethod
    def create_class_fixtures_label_sequence_map_based_on_state_alphabet(cls, state_alphabet):
        states = []
        for ss in cls.states_symbols_lists:
            states.append(state_alphabet.get_states_for_symbols(ss))
        CharacterTestChecker.create_class_fixtures_label_sequence_map(cls, states_lists=states)

    def verify_char_matrix(self,
            char_matrix,
            check_taxon_annotations=True,
            check_matrix_annotations=True,
            check_sequence_annotations=True,
            check_column_annotations=True,
            check_cell_annotations=True):
        general_verify_state_alphabet_symbols(
                self,
                char_matrix.default_state_alphabet,
                self.__class__)
        GenericDiscreteCharacterTestChecker.create_class_fixtures_label_sequence_map_based_on_state_alphabet(self.__class__, char_matrix.default_state_alphabet)
        general_char_matrix_checker(self,
                char_matrix,
                self.__class__,
                check_taxon_annotations=check_taxon_annotations,
                check_matrix_annotations=check_matrix_annotations,
                check_sequence_annotations=check_sequence_annotations,
                check_column_annotations=check_column_annotations,
                check_cell_annotations=check_cell_annotations,)

class Standard01234TestChecker(GenericDiscreteCharacterTestChecker):

    @classmethod
    def build(cls, labels=None, state_alphabet_fundamental_symbols="01234-"):
        seq_symbols = (
                "43?2423132-330-0??21334--33?330?32--??-40?2442343??20--???0314020222401?-04?14-4434--20--3343142240124243-3300?-43-1?203-130414232?3143-1143",
                "03334?3?102220-44?300??0304-3?-04?-0?1-43?22424433?002-0??33-4024--010010402-2?4441-?22203-?4300?1013--33--303342?3040213-?101043313-42011-3",
                "-30-2423323?143011?23341-330-24?204---?424?4410330--00-3-0043?3-32321-1--40-3421-13-2132??23024?01-0342201041?1-4??320??-11-314?221-112-?-43",
                "?3??142?01???0?1--02?13?--0-014-?232??431-2200?01?121--21243111233?-4001-??-2??4413-?414033-3342--040-043-2??13212-2320--33?1-?0221-1213?1-1",
                "4042-0330003-00--2433320-4??33043133-1-42-2212141??-231?1?2-34??0422401-4?--10032210-2??2402012--40320?244443?202332??434-?011?4224-122?2?--",
                "?-101213-31230-0-0020-44?40?-???321212--012-4034201300?-3213?140-34?4013?04024311?4?02344?-131002?213444-2?-04--43?1?-02303-3-4-0?-3421-1-?3",
                "1-03-130-0-34?-02331234011130322-?-4-?024?0142?334?2--0??-13242101-2144404--?1-?33?44300---3?241220-?4-43333303?43?1113411?04?0?1300-0213--3",
                "0341323-0?3-33-00040--22?314301-42?-?4?431?43?24??-301-0-??22--2414-331?-34?-?-440?--03---431303140122023-4-10?0-44412--?3?-?14110313-1412-3",
                "?1?24231?3210-34?1-?0410122??-2?32-?12-43320213441441004230??20210-2021?001??43413--40-2--003121122--0011-3-2-3-4303?2?13131014233?1--3???30",
                "330023-24-032-0131-102412310431203-4104-3--2-?43130?102?2143?3020?124441224?241413--430-?344?4144144214-4000?03--3-0020221144?23?2?1?4032-03",
                "4244000?042??0-214?2-2-4?-31?00424-04-1?400412230?11204?-104?-421422?-42?44121?-203442213031?31223-3043214?104333?04-4?3-13311-40?3?003311-?",
                "?432?-312?2433-03?40?-40-2?0--0100--?-?413?201244?422?3-?133-?0422311000-0---00?43--12---3?41-423000241004-343?443?044??4?124042--?01313122?",
                "43?332?00?-4?24024-0331-210?233302213?-402?3?20414000--4?12-34324441433?304?11?4431?-24-341104332-2-1-44313320?213302?023142403--440223014?0",
                "4243331132-?204-?22-3-?4230110313-4-?2?333-24114330402024?-3?200-?1-0?3??04104-303---03-032431042-3223-13-4013?3-3-40441-4?0244-32?33130103?",
                "43202-14024240?3???13-4?210404?3-3-4-41-01-340?422333--32123?40?02034303?0?1123?01--?0?-03-42??42011042?333324?-31?1140323?2333?30??13302221",
                )
        GenericDiscreteCharacterTestChecker.create_class_fixtures(
                cls,
                matrix_type=charmatrixmodel.StandardCharacterMatrix,
                state_alphabet_fundamental_symbols=state_alphabet_fundamental_symbols,
                seq_symbols=seq_symbols,
                labels=labels)

class FixedStateAlphabetCharacterTestChecker(CharacterTestChecker):

    @staticmethod
    def create_class_fixtures(cls,
            matrix_type,
            state_alphabet,
            seq_symbols,
            labels=None):
        cls.state_alphabet = state_alphabet
        states_lists = []
        for ss in seq_symbols:
            seq_states = tuple(cls.state_alphabet.get_states_for_symbols(ss))
            states_lists.append(seq_states)
        CharacterTestChecker.create_class_fixtures(
                cls,
                matrix_type=matrix_type,
                states_lists=states_lists,
                labels=labels)

class DnaTestChecker(FixedStateAlphabetCharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        seq_symbols = (
                "CaY?dbdNYWkKN?kHMXtKWXRYYwATrRbD-?hTacgSRNbkdgNSRymmNWwnKVkYGAkmRAXXhNtgwnmGramYH-ms?vNb-Ca?DNAsdwRDc-KYd-HTCvNyRsXVwmDbBYRaXd?mcRdyKwgrbdswmSGAYndNYVbbtAKCgK-KykswasSYmmNGYGNNKR-dwgMGGTgBdhSdsWvAGKBnMK--TVc-KMWkHdXgBAVTtwVntrcRNwCRVdVamYXYvhbh-GnswKRYtscr",
                "kvYYWbTKYCMKNHkdsXbHsaCTYHtXN-TDSk?AaggSrAhrX?YbBWNmWVtBGGSXBARmYAVMcwtgadmsraWtwXasTYX-s?asDmXDAbRDb?RYWNwRdYXyMGdB?R-ABKCRbMaacBcSyw-adddBvGKAYddNVcAHWAwtVAbKycXgVKdYAkNWDGXAvrvsHgBtGGsBKh-RnBRmr??RTCwAmwDBbGVarMCkBADChh-cbAc?NGCdtdwCcmHYvmVtRTccSScYtYcV",
                "rNarTVWSTXTKDhNSGbGKhXwwTsGvwDMDN?agmcXBNKbVdRN-K?mtwWBnWDks-hkBwaS-kvcDRVmarVmRX-s?raSAhDaGDmmgRwRScnWYwvmTnSNyRsXVtMYbBdRKSAywhRdKwyVBAAAwDNaySrvNyTYYVAWCXymKSssNVsmGvdNkYrywYRBvwTNBGDgwRhMBBW?TrKVgVSG-?BcvKdanrvDgSAVvMrvCwybwNhMG?WVagBXwNwrcytWmgBgSt-kk",
                "sgC?dsaBYhkrV?kYwmDKaXRYYwvT-?YD-mXhVmnKGMTkdSKsKkrmcdskKRkRAATSwCYXvNWkbSNGrSXdtvC??kcgMGaWMyAWdBHHvMYad-HVyGNCHVXDWWwrrYHtanYNcM?yKSVrdWsVrKSnT?tXvgbDrMaCgKarDrKwskHmgyXGYhmRK-asm-Xnt-WbG-gGTKvcNKgntNrnNHBWMrMkHGg?d?VVGvcTt-cgbwdHndsAmXYYAVdt-GnkgvRYdRcH",
                "kaAcdGnNwmgcGkXXMKNyAXwRV-KTsRyTWw?TMwgSRvRVtyycGYAtHWmhXVWYkAcMtCybwNtmWw-sKRXrydmdrvX?NkY?cscYmKtdc-KaKr-NTkNydKGRthSbTyRmkdbv?RvrKwBDS--wBdDARndba-S-tKRygKcKWraSsmXvbgmdmbMmcshAwnbGSTsCXrSGsTkWnKasRyhMTcgWYBTSsdmwSAVM?tmhWScRn?kRBcHcmSkYccdhSGwr-?SRC-Bn",
                "CgH?dcCwY?VTWvgCACShvNRXY?dHs?BTSYHdaMgcRSCdYgcSt???VWTXhDvdyAgmAsWDnASNRn-GWttSKkGN?yND-SXYYbAcdmR-mwmaCsHnmkbyBmHDKXDkWYGtYDcTvCvcSvGwbkbHymAkBKhNYayhckYhcy-GXhRnGBSKVANdAgDSnDmdH?MGNbnBXhThKHvaMNnGDKNKdAMSKMSTgHwcvyYcwYdrtDVVtDAHVdkyWyCsvAGWtGkVVMya-Gbr",
                "DRcwhbdwsDs-yMKsHThKdwKYsvAArcvwT?MTAMgmRntHdTNknwSCNvnrWKkXMCKtvrXYknSwMWvarRnAsVHsbvNbKHcXKbHMdGDTBW?w-dGcsvHWXsKVdKVnBvnNcYCWks-kyyWwaCsrdGXNvsNHyVDaHAtCbnXThAnwaBSXHhVskXNtygVdTdMcsNWnaMkwckaXwt-hWKWBwG?-rMb-HAvyDArbtRVNBRkggmNHHd-DXYasWASd-rmVkKwVthsa",
                "raYsbGsYBaWKkgkHSWCKGkmCWwwTvRbD?RhwvcNkR-GGbYASVDRmg?vn-?CCGHtmWgXghCtmhdHGWNtHKdmkchDbSSNgDMAsnrRDYdkhDC?HXCgRdmXKSyKrWdYVGYBmBaGNkmgrchWdmXGdNGGnkThStWmNkV-mVKasCnBYRmN-YkyGHCW-wwM?XGVaMkSvgWAAnKwMkWM?DHyaKR-GAYXkmRCVtRbywCvVwGWT-C?aHnhtWHkHmCysXXRbtMcr",
                "cHgKdtkgKKXygnsCtTGmVdRHXD?TAkmDrHsTdhVnwstkdXCGwy?ms?wnsVkYkYkYccXXTNtbwr?GM?yXn-msbvBwWRB?XNSSmcRHB--RbHywBaNnrMgrDD-n-wrKXy?mcN-MrnDrnKgs?StSBbVGtNaBvvYCw-VrM?sGasSCm?WdYXsWmABCwgMYBmhXMBtwmANyNGydGcsgyVXVrDRwXwX?BGXKtrkCtWKm?wCrhDgvmwSWsCCk-GmsKWVYtaca",
                "ryTBakDTwVwd-?kHMArVrXGbrGbCAVbaH?nTrygYMCbRhsaSby?rgBBngVCNmAk?KRdXb?VgGnNGgatYGBmsmvNS-Mw-wHACSvAKcgvWWNHtbvaaTbXXwcDYktgaXc-KXRMCacgybR--k-KNvAdNYVrTCGDv-d-KwwKKayRGDmn-VyDHnNWdnSkGGhSwAkNdhWAGCKwnhBw-TYXBKkW-HdXhdhnmdCGnTTcRwrXs?mmnXwVNrgmhkgMRgKantsct",
                "gvYXAkgCVryKdaYRHTdKW?crbwaTYnbDSbgcmRgXaNw-YkRS?dsBNWWwgmsNVGkwHHyBvbhXwH-Dr-SGDSwN-hvbMVawDTSYgwRWV-WYvTHTBvKTrdXVvbgbckTWhN?R?ndygwg-mdsVghV-YcKNsVcbTAR--mhTyAcwMbnCsAyKYMw-bVAdwmBXGKMyahSAAWSVTvB-MAGYmSs-KXykHt-gkKVvSdrnKHAaGwTRXnVrdwmHvamcvHtawYNTksSc",
                "GgYXWRTtbGhw??r?MyGHWcRwM?mRrkbm-mhTWvgMgNvgygMScyAwNYknnKaGmVGYRAX?htMT?kb?GCkbMVHs?NTysCD?hvHgwVTDkRWnmyHSyvDDRHRVXnRY-HRaCd?XMRYYTS?dHHVmmYbvaRRNDKXGtbhWsakrrakYDsKYmCGvtnNybhAhyNVgkydsctDkTnvAHtyAmXWtaCmGKGBNH-SMMKVTsNtrVbYRyccSWdhKCMV-a?sY-RnSWydYWWSr",
                "CrCAsCndkdmcSmvyTtdTWXdTAtAHTABXWhhGKcgaMcY?N??SGMbsncwvKVDWMXhmDMN-hNvXKHtNWGmYStm??aVX-GSYmS?sdnyDAMmvtXnnCrKynT?WwkydXYDaBcCTGrCvVm-RMwsTmbsyy-sWYwBmtDKhVKrX?a?masbYTmdYkGNADbXdKSygMKmrdMgdb?SGAKBnHK--TVckKSksgDvgAA-wkwR?mHXgswCbgkVY?YHAAhNkvygnd?RYdsht",
                "tWB?d?akmWw?ARCHVTKDCvarAvr-YabdbmrTawYDVNwkHXNNdhWrkW?nkBBTkSgcHHyXctwgwDSnGhVSH-yCgghbhgyMgmCtXwRdrRYVdMvYbVcBgRwyHmbvSdXAb-?GwsYtHwWAbDtkWrGw-KVvYhabM-WCWbBKbBsTaySGTsKGvGYtHHbdGYybW-gBAdWdhCKCKYanswvbsgsvSnhHRKAgb?yKtwKmtYMBHYHWMB?RmYSYnmWhCmVswdGSVsHM",
                "hWYKRNdvRWHKHWwHywkdXXKYYmbByWbbmThTvcySRnrk-sDYGV?sMHMkwKYdtcDVdVrTtSnCRnHyR?YYbNGGNvWbvagvDHTgMS?cc-KYv-HTCbDyRCwVaHCTYKsarAvKcygycVyhNKDyraGgH?mNMTDbD?NTgYtbtckwSbWCaaWtYDNCKRTdwKKwSnhyBASHsWbWRKKrW-bBdbS-ahnkHmXvdAVnmwVgBTGR?VgXCAHVmYX-VcNhVhgHwrgBkwtd",
                )
        FixedStateAlphabetCharacterTestChecker.create_class_fixtures(
                cls,
                matrix_type=charmatrixmodel.DnaCharacterMatrix,
                state_alphabet=charstatemodel.DNA_STATE_ALPHABET,
                seq_symbols=seq_symbols,
                labels=labels)

class RnaTestChecker(FixedStateAlphabetCharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        seq_symbols = (
                "Dg-aRmgGyUAVdrahhAbXbUnmvURSUuwyDWUwRSWykcvcAMXnmMhhYcvKAGbuDvXRDwnYSVGUudSNnBsUNAVSXAGYNGDVvDhdSBMygCmgUv?cdB-DRwdbcwvKG?XvYAWBavVWSdcDuXmuuvyCWnBCAw-CcksNKHYHUaV-UXdbdcmbYXGGmnyNVXKvukvbRCuWvKCVdBwvHKaSSNcDwD?VbYWHbkCrSDhnAwVc-S-NGhurbhrsKusWWyMHYRuuYWRV",
                "DB-nRnBGruAVnWMhaGbXhVumAMD-UXmXKXwDCuaysBsGAGnnduGhYcbbMHBuUcHcgCnhSVsNumWwSYADcwsShMrGDYDXwNAYNBrmg-vbbUNndu?yyXW-?ArKYDXcWAWByvK-Suy?cXWbHHCWgnMScVNsUwNNWmvUMaVNAXvCBaaXUGDannDNVYuC-MmMHSugRBwVdCBMswNXRKnGADhVXcWSbbhrghSnAnSuhSYNbngdbUH-VuuCWvSndM-uMWcM",
                "MWCasmyngUBbyrWUCHbmUbnnyWWmUnbSDbYmUgWBkSRcVhXnBMGhshWRMRAmmuBkAUnHaVHMDnRaYWsUnRuaVAGYNgm?vNhdCs?VgkXgyVUCd?kDmRMV?Uvkw?rGdAWgUmgh-WgUgAmdyHGskBRCUCv-CrsSdvmYUw?GHvwYdwmhYDbsmDMrcKUsukvRuMVWXSCUUuwWkKvK?rcDXXMNdYXcBkCUSDVrYHhcscVNVhdrbXnRUkHWXhrVYRMbYWwX",
                "dRswabUBmUnndsahDBunbhc-vhRSUSUHDXDSbaVykr-mVgcnBMChvcYBMmhunCDRgNGYCaHUgcmrBn-UDMUaXMgKBGDrvDWXSWmRsdNgwvrNdBBDWgSwwkDUsXCWCcWBABynSHbR-yRYmWmgwvaUdddncC?RygNnyyUkyDdbsyBAcUNubckHVAUNuMwayCuuRKbVdubvcVaDVXsDNr?KAUAabRCXknscAyVrsSNXY-YKbNYsUWUScbASrRcUYWa?",
                "BRdGb-yvaanKbgbhdcCXbXnVYvRSUuMvcWRw?GWuWcD?AcXvGBhhUDvM?GbmDMXygYmmunWCsd?CyBskGCVHXYvCC?RAvYVsSAKDgywwkcacwBhDuGBkWwNkCnDBWbNA-bKWkWRDCHrKXVmCCwHcARNCsDAKVHuYgkV-sXdN-cmWkvXkuKBSb?mgCuXnSUCWNnCkrgw-vV?vANRnVdhkHrWBvUBwR-hn-wkUCGANGNRrCBrsuusbNWvXcDkudgRn",
                "DgurCDgCyYgWGuDhhADXMUSmvCRGSuraDNnhwsCN?WvcdNWADBhW-bUgnSbXuRvkS?KhSYyMyKSNAhMsCAsggAVYNrrAYwhgAhMGwdmyUyRddSSUdVsHDSYKYKknSDGSMRdXNRdgKGmC?Ay-CgvBnXkMdauwKHY?waGYAXNDdkRMhXkcmWyryBsvnmrbdCURky-V?gBs-K?HHC?nG-yw-aKAACGvwGM-?wVkcDXg?dyDBCCbKNaGvyYUaCWYaaRc",
                "RHdYR-gskuabW?GhHHXDhUKdrUgcUu?yD?gsuSayccvSWAVHCXAhgMh-VgRy-BnRd-XH-AGKudcCmnUDvKVVnAGSDHDSkbNwSSkSgbmUNHhad-GcRwVhywAKGwasWbyBycyWadcXNWvHsMyHVnBDsnbVVDGnhrBHDHA-nHdcrAHXWBGGKnkgdCKXYYYrgCSucKCXmmgBnBs-YNXVScdrbYWuWsMc?DhXnwVgaa-dCvb-bgdDmusWkWCUGRuuWSHH",
                "DsYgDWKArwA-dAKgGAkHWyAsrkCRWaCWDgUwnUdDkMUcAyKysMhhYHAbAWsuwYSRXKUGmVvCDhGmvDVagAYSnGGYrwDWBr?ySBcyrCXWkDrXknDDbKNHKwbgN?svkAABrWrwSrURuXRWRvyCcCSCMGgKVGVNKmwdVasnrYdbnMmbYRGDYnc-KXHD-yW?SgXWkUbVKvwGduaMgvcDGM?VbdYgnGy?XsaKUYVmvA-dgMCANaXshusWBHYSunCuYvRr",
                "rWagAhsnyhrVHrRRRmbbWBDhdURHdgnummUNHVNhHavcKSbKKguhkAYDAGbWVvdNDwrvhubhRKnNDaB-NmBHkDG-nUrgvDhDHCSvyAMvKnDSdMwcNgDRcs-kG?YWvyWryvmWcAcDWwyGuvR-dGdK-rDCuBWrNYgHMVWvc?mdbMmbhXGVny-NcCCkV?HbSWbWUYKMrBcmbcHSSVdDas?kSAWMbDvWCDnnAbNVYRADADNyrCyYHCBWAyDGcvVubKcw",
                "vhDAwmY?knSgGNyarAAVHgsXVU?uRbyuwsUkD?AkkcvYBMsvmCb?VBkrAGv?HUXvWBWYaHGvvAuKnDsUXXVcMMHYbBDMMYXH-mMUXWMBaCHSYm-DRwBbUbVKkVsYVCmrkvVUMdHgUWhunmyMVVnC?BXmSXsAMHnngUV--XNns?vbg-XDnXCNMsBWGAgdVhrwsdNVNcbBCMarsRcsaVSVUYN-au-GW?dMgCvAkDKUSBVkWBhgvuSCyuvNYCsvCmkV",
                "mmBybMhGmvGVGNVWhhnWNdUSsVvDcywyXyRXrNayGkGUAMBsKMCWM?AKbWnuarghDannUVSBurSdHBXsNCmRYNGdYwnhakXdSYMYyChgBs?mdb-ARX-USSDnRkaYwcWbaySUBdAhCs-huvgwWDBCVr-rckvmKMYnUdCkYghwrvyDMKdkrmrNSY?HDkGRHCVWXyCbrKG-cSmXXdANwnCrsuraRcGBSDMUWWKXAb?kCnnrbhmmKuKGhgMvGRVhYVDy",
                "-baX-rKwGbAHVaaC?HKXbcSNHUBBKwwHDuUrCGWyvKgNVvabYMRA?GvVasKGDARdkwncSVbmAVyvnBHCkMWSXrVGYUCcvnDdYBsVwsYcCncchUYYWwuVcDRrS?sHBGwHBvhMAbdDuXucmXycubMCXgSkKg--BhRBUaUGMRNvaHNcVDyBbnkHVCKHXkNhGdAGvKCbYbsKVKMkXVmVuUVVbYaCWbnDhaRHhyrcaNVg-hcRh?rdKhsWWCUYXXcddCYV",
                "?dSa-ugKUUHkAUAGbAnnbCnmc?GSscwVDWWBKs?HkcvaM-CmaKca-HvNdGbuHNsRGAdYSWBVacSRUS-NKrvSXNCYNGnYCShHca-ygYHSyAcHCw-mYVNKcHvcsuuvmUbBHUNrrRHDhXRKGA?duWKvWycCwaUMKvnhRah-UXamdYmGKYnGNnggVXDv?WCnCsbHvnhC-uwRHKDWNmKWKDrVmXSsbk?vasawUwgXr?bNGhcrchXUmKmwWKMaRRCSYNVs",
                "WnkUsdhGyR?VHUyhD?SXbKUmRUrMvCwAdbUVDKnHdcVwHBXKKRWmvDmKkacSg?YHHsXwSbmXaGm?nhc?Nd?rX?GkmrvmNDCnNC?ygruCMvHUABSHRvYYc?bWARXRVBamAgVdMdUDuDW-aWCdXmkHd?KCcssu-ruDUGgw-XSCdawUW?yANdrRYHvcVCnVvSbgkwduKBwGHRKRXvDyMa?yUgWDBcCyvcRvVshU-S-GVHugbWbUrrYBnSBrURwu-s-r",
                "DgkVwmrSdhABWvaMWbdXuCvHKRHXCuwXrW?mvWuuCAvgAMnwrSRhGdWkrbbHmhrw-wwdvBanuCSNMkuVsMU?BgrYnBhWvaumSYXydnmHGGNy?M?CHWdbb-vrBYSMHRWBavXXgBcua?mMhsysdngyBgVScwHAcHsDHCVHvBASYcHbRmXg-byRnRMhAmRwD?w?cKCwdVYhC-hgwghWND-nrgSsWBmUXsWbnWVVVYdyGYucYh?vUHrNW?YRUSrGyWun",
                )
        FixedStateAlphabetCharacterTestChecker.create_class_fixtures(
                cls,
                matrix_type=charmatrixmodel.RnaCharacterMatrix,
                state_alphabet=charstatemodel.RNA_STATE_ALPHABET,
                seq_symbols=seq_symbols,
                labels=labels)

class ProteinTestChecker(FixedStateAlphabetCharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        seq_symbols = (
                "PGSD*nCdNZ-Pa?kDyNxtPdiXh-zNnmYyAgEqmZyXrkYXzhAWMWSwyRSTRGDkItcmXNhyKAnfYlvTRqwQXrapfmLkvibbXQGleGksliC-EmNewqgcARSlZDWnAzXwEaSGK*Vz?yVYDvMQABZNrCCrqEmCFmiShDrpWZBqCSQyblBEVSlCHhl?FITpyXklyxihtA*mcbCLnMq*pZprTSGRDregeEC-SIGNwZb?tHRxqCibZeAEAGnYeCk-TVSDgNDegPcXfwKCMziIPbnAltrFzwxFksB?yChwaZ-XpLEyDZdPB*qFip-dGtWHeyLGm?iPPvdkRnbEaht?yIVCrXhQswSwgCykeYChEgVKDLLKzmgxqmhkQ*KAbfKiFVXkt?HveNZbDExKWdctQDQvpVHIYHtXzNWNE?lFtenNrMAnX*XkrtVRvpptYQeKwDKxgcnXKtzB-TDHyHd-MiqaDvAFqatQ*??KlnVyaDAHkhkDkFKcrTiAq?a?mDige-aV-iBh",
                "ReMQIDvVAvwAavnawMZhglSXYnkYnfVyNgpiviyyNkffrhGf*cSwVHzELYhkDGDmeZmzLKBHPDvztkDaHa*NAmLBv-bNfEcTM*LszffdFHNfwWFIntSAaqAzDKcfEAHMaiQlfqpwlKdYglZlCCtPbEEXEWfShLVqCZfqCSKyygBMBaDQFBlaFYYCKCDXyfLHHqINemCvXaqaxXprwLeBhpezzEEtPcWqtWMgBf-NKPPtZWvErFlqewkKTTrDwcsflScXBN?AMtZfHYPtanrEDwTiqHGWlaW-an-ZDrEazQdiBqBgtEXlrwWteIgnS-KPPfXmzzXrWiCnywVC*EhlrxAPgCyZRYCTxNWKpVWKwmaz*KM*Qt-LAi-SbIziTxzAtacFD-xKWmcHYDexieHcAIwHTGNNRTAqwyrGrwHdXTSvKtVwBHpVmQSCmNkKl*nbKpLKITWWyi?acQsmLvzNqQRkb??KldaPzkhCtxQEZFzcqIqkdKM?YlXnTqmVhPHH",
                "PdpZ*BxHsZA*arThBiwtkd*XhWSmnyY?SVaqmpWeREG?*hDWVASbyXWLnGWXWnePREzyxpnfHivmvrCQKZrkcgpGvVbln*PylgfsrLZhRtsizYyVIFYfPDvnMqKnamSAK*VrWyE?-GzysNrcRNCzyCTKFmfahsrayrzbCqQyC*BCyIAZNWR?Ta*TpXAlGxT-iVBmcbNLwXQvpcIieqMPXlygWd*-NnqQbvXVahRxqgHsY?baAGnICskRLdiDaymegkhVQyix*zcw?YqFYwMCqwrFLxQKgCWHasii*QixDfMPVLqQlpRvrWZgenfLRCiyPvpfNBZeSvse?XrICDcwZqSwvCrtlY*KEgVCAfzDImQMInhfTAqRpLaZvVXkXgPiwMqqELQKdDrQgElKzVLaywtEGegsa?CCtAcQpPAGRcGlkHQevsntbmWiiXKxgcGkIIbm-LRHyndSM-DwNhXbqPzmnrFSH*s*zDgH-XkIkFLcWqNLqWaGmb?Qd-PayiwC",
                "PwSDvnadzFfEby-HyeptwmLPh*cNxmYCkgEqkcCdkzYgeiTgbQsTiRSTBnDEnL?iehcXmA?fBdCwbDwQFpadedBkikzbVQrRmzkGxCCPlpRTCqHdhRDLZQHyXyXvPe-Gx*aVlBavM-ZQW-xNhi*MrZbmwmQRhaqNkZal?vhDwpExyrlCGPlstsTW?ZLNXFSGBH*QnYbwPMIFpzxri*PlnrwlyEviRImww?XsydAxxNdbvzAAGAQBvxTEHh?ygWDLgQarfxrLMLnBXsnIlVbGztTITsRXysfYa*isaRgderEdy*GadxbZydFVbbN-v?qMnvcLKhNWdBLhHIgrrvAQZwSyzTyWpBecACsnBxLdhmggqKhcQiKnIIKlNzE-y*lA-qEblETKFdsdwXcATDaYYdtdWAWvWzaW-PHksD?WWQVPiHVINnwYtQmFeFixeWQdAbBV-TiQfHaYriIDbvnAncCcm?FKlz?yDDlksvnpkdHnQlvvi?a?eqyVqAyRYImc",
                "DGAQnnZdPkEkapki*SxmF-qXmZpaEmyvPTPTeVMqPkztVnbTQQddyRLMRC?kmSTMyTgaBKsfzSHxwLtPn?NpswBrAQeSXQGqIGMYCPmrGiCbwhcCqmSrZqWlMYvwhavtMBMCaFVyaqFNqgWNrzCAImVbbhySvDxq?rWvbSfhbmEEadrCwzGxPmppNgPl?dVYqeWvSTqBIYqPsqbrrSWqHrmgeIdwEVhXp-XttHAnWtbfFkglYGKnKtMtCMHDtG*eNYCKbnPHM?rygDWAFt*??wbxfNAHyAhAeZXayLEyMStPPVF*BiwsbDntfqZcqsiyETwaSNIreikWXkDIreFEsZihaIPxenvpEgMRlkssxhxPRQhksXvRLXqcSavdC?HvRQBElMVkx-ztKbMvpsGrHzmezeLNE?aFP?nHDvpv*xd*S*meaiXmYcatSDkxmbnbdBQa-pnHaGEsPiAEEviAScrbyBDvlnXnVmBGDEzeYXhfBc?AqLI?m*i?ewfVkiBh",
                "sWzDCnCiNyRPesGQMNgGcsHXa?FXnhXIAlEvphYAAmYGDaXwRPBRyRSIQEDrIMQlqXIyexIGaAxTkxwdn*zcTmgmLZ*bElnBBXgsziYLiDNFPBpfAvSlhtSvg?Nwh?YvCkIy?yvnYPmmpygLtVQighxqFzASwxrp-ZDdsRPbtlKEPGAgqHlyYIndLXTltFRmwnhDcyLLgTw*p*NpDfPpszaRfsCvSIGWEZSeyktHbAibRlwLdtDxDkVNdYiPvQW*Qfib*AEemyikEbSbzZtVNYItEsBGXCrtTD*Xpp-IDzMYrXIEidntMQvHryTEsgIIkShqdKBEQsYhVTQFM?zQTwarSL*iftErEg-BDqkPDDgxQmTrQGXAqfwDV*xkhEKBLeSvcstYNdttgpRAbGrCYmFnmyikk?fqtsnbrMDnW*WVpnaKtpiwIQbSwXxxycIZKtz*-hQpQHW-?iRIDEILqevYFtGKEqIySDMM?YpDdIiWreTRMXa?nkiTeVTknGnb",
                "RfSkGggaVrATryLDANbsidmkSP-LnYMyAgEWahZ*FkYtvhh?-BnarWGK*KakIIvTkNYyFan*YbvTl?VBXPaxfDFkmyCeGkGBsklsxiCZKmNXtqaczrrZWD?a-V?AslNwEKyzweaeDSMKGSXvracMFagnWAiThDhdWegtdP-skl?ZVRSccnlQtG*--XGFldixQRr-LwCRMEZrpZhGnYtkMEkgBqC-SGVNzsb?y*RWQciLkvvEfwExqziwdFfsCyvepivLEPACMPR*HlxnxprdZtbv-DQ?lCkrAZCXyGH-bZZSPrqkmtvMGPW*EmL-NBTtIEHZeNhHMcP?yGvpakhELfGwlGbsCwdhx*lzSyxPylLvqfFkQXLnGAKNKXcMd?BdmzxSbgxRPEPMiIZlYtsmrHtXyVlSbaBctdcNcwicz-XZeAEkdbepzEYEG-KxptnWYQnLNTpghVItIyqailAF?nEr**-?M-gLaWAXShdDkRecALLIY?e-GhZFs-aFRiBb",
                "yEws*nlbes?Px?ISyWBBAyt-drGgmmLyaZEmNZyEcZLXphAlsWQPHLSLGGHqCIyDwQIMvAGfYzvZtVTfpWyhPgwvxPzGXBmI*NyslVZSaTNeZDgmyELMxlcEvzdVEaSAK*Tz?yVYmRMeALEYGsryDVQFFkHMXkapDfYkWXQqQcAEzSStDhnQyLHXyQMwyviAZHRiNWILF*W*XkCrSwiTDx?CQiCpvTfNDvqbHx*pdCiHrbWEAiSQG?BHTfXgcwDeMPYXfw*WwiixTlYAbtTwpwXLXsl?WxlygwliptmyzZdnYmqtiGbVYQlyYhnRdEdvPGBzknZ*xmpzyeCLXXXQLwSlSByk?YNhEKprmskFCmQbmghNk*KesAQcmWskb?LagwsDLPCPWqawHyKvaTBI?HkwMPqlTcP-GXKzkMEFG*-*MBDRvpH?tQpNmkFwctWbiHtg-vMGYHN-fiCYrvXGlheH*iMbha?mpKLYXcaPpMZbrQZhqETiDziQN-mcIQMF",
                "ixAgxnctNkzIVpkDkXDVBADeZWfNVEYwVXbyIZqBPkrAQvEKMMqg-wETqDDhqAaLdYZgKLnbxlvcRxwQXMaRfWHkexVsSeGlwGBy-iCdErgFpwGckeSDVFZMLZX*laxggInX?kfYcGq?qlbNrIkMhEmiYirVwXrPHidvBPVLbQbEQClbabs?FITpqP-YyvyaYsZmNCmLIMdhLZfHTS-RIcGbwkpHSEkmtnr?bnmwBdRbveZKyGQhYC?fT*SIxQDBZmdZpHKDYgZXHbdAmtcztw?FvQTcDCzFeZcgNFyNVpqahMPSfz-dGNdIpKLKG?dHPwdMRMXKaGtXyFmCrXZQSBAH*XhMbbKXEymKntLYCqcaqANCWQYPVipvnVZbDMHDeEiMekLKWkcPLMdnhPDIhydCzTVNTP?zreFLcXpMeWpTytaEgpptYpMKwbKHgWZXBt*z-LGvCI*thyweqvfYq-dQdLdfInqPaXQwktSiyADwezqNQziY-DHE*-AhhNBh",
                "FcfIDnLdQZ-hQw*DypxNNNiYhxpNVyYgCsHPnEyXrTcXDDQyaaVYI*sSRdYfICgly?hwhAgNiFHyPTzYQMYpFGhxqihYXmDbeKLbNirXkENeNEgWFQLZImZXVSlwEaSarK?DdyQY*wM*VhgNlnCSIE-l*BicSryZFYbtpVKxrlDIZFKCtykQXIFpt?zlyHiKwmw?KZNSwMS*Zl?TZkGqErpgVEvqSBTNnKhhtYp*TrDbCg-EAbpYWeM-HsSsikDQyLQMBExCihvIeGEkeeaLYhSbHVsTqwSzZrAgpmySDNaP*hqezMmzeRfDeFZaRRqPtdlkRvTEHWaNF*fWlgfxbNnwhqKkslbWVFV?IYtyzcGxqrYMR*bX*MWiFdLiFyHBfewsMgxNWdgtCFkwZHHexB*nINanE?AFAQ?drSG*xdkkEZsLvpMbYyaTwcsGmKwyMTtwQLgcNf?-GSmXvLbHmTpgBTPXTiVwwDgHkHgztsKkpxMAd*TYImIWeLeVyifL",
                "xyVIPHSMEt-V*rCv-rsthbAYEWVNvvWFAtELZZItrPybcFxIHfd?DRHwp*fkGZxiaBDxFAINYl*YlLwPWfbpfcvgLvBfdrT*e-ENYil*mVkvmFMNApfKix-bAebwwRSGXaVF?ReGRrMFARdvlIHcqEmSXE*khKalYgnqIhCyfDNqnZsEsNlEEIPIrZklSCNK-?*mmDDZnfQZpZQCrSLrpP?y-yG*WrNWwRw?-HzxzlYpZRlbfaefeeGYYYkxqvNp*zlFbVqLtzYbZYnABQrQz?kQwh*WGbXDN*sXECEygEBwLrylwsrdGiDHNRdGS-nprNXrLnbLCXmcyiKxQLBQsbSFewyEeKlmNghKVRrZFyFb*lRmXsKA-f*PIpBLzebitNZYLl*KKdQtgHfvpRWI?vtqtNFRmyvtQXnFwbAnIAXRm*ADvVzIfwCKweyCK*YXyeKDLTvXKyzVEiNmK?W-nvLs*?ldnnRylAsfghvDIFK*yTsNAfadSMefeAymBe?h",
                "hb-t*LxhiqptnHX-eQayxaVsMzyPnmsTxhBMgBNXhIfRzYAwcWQSqHdMWHdkvKYNbLhvvAlfmlsTnzIQ?XipsYhIEcrelCqlhSksMiC-TVftPGgTC*cTAmv*LkRGEaMGbBVAFyfQDsrBA*aqFXPwwYtCCwhSLgCVnppHCWeFsaH-IhlBFhl?lQ--MiCr?V*fkCQmDRCaHKffdZtr-yGRi?etfSmIZI?DgRVlYtIxBbzbZeHYEHaiaCiPTSyfNiDewPp?GcWCRCiGtvABllrgPDbZBDizMKZnsFX-KB?wiZqa-btyTATYBegHVmLVAVPYPvdkLWTQTXtkxIVPGXHFfAIWpCEkFPfhxHyMwDxtTmbxmihGq*WHmfyCniX?EtYSGNfpnLxiPerxcDQhxYDIcCQLbilIyAThPFqCCXeGytXkgnVwIkptddIlw*nYITwWCtMXybHEslNLk?fK**AFqBTQeyiKCfYBACAAiGZwrxxirEBfCgP?ZsmgDFGEwXBB",
                "lQqDMRndBArwrnyDMFhtDlLM?crYtkYmAHHqaXT-RkbNQxhhMtmXyYSxVGLhreqdvwvSVAGznlNKIHcLX-atBZLkhabDXAhxeAMPyGA-ChmqzRgXAAD?YkXZD*WRyaMcwBWPvQVrDiFQcEZgrsPlYyfCFqXsSVPhWxVqpmDyTYBgFXGvZhl-RnTXmXFlQkw?yzVdcvc*KzsZEC*Wf?MYHceXMEXttIGNeZb*vFafqDqLVIVBA-nKe?arTVSCBNIHCPatziIqkzGsIWzEllzvDlxSHfQFvWMNitKXeytsDddPRkq-dw-dteFNexLLggiPaIYCRSbwFht?ZI?cc*?QHISweCykfyLhAyfxsLPKLmYIgbckB*KqcfcDX?tzkXHeIDabBEdFbkyMQcZvyhkrFflGIdMNbcnFBLyXdXqNiFpLmxBsvmxVVNBKzfsxagKCnMWQVprHVmhFCtqaRWCBQ*tpft??lyWfaaAwXl-DxFSN*YTAfME?Bwi*gbX-ChZh",
                "HhPiTngnKhBItWtDfezSLdePWpzaVFyvqEftLMtq?-*AZRwTWCew*RyHRGgAIDCskXgskpnvYyCmlmgbAyPRw*IklvbdXsvdZdasexwBEFnELqNcQEQlZDtMtScSyaSipGqn-RXQD*srtW-NEtnIGEwCBCxakWrIkxbfIIIcKh*sVSMZThwcHnTdyrWyPdlqghdDnbwkkKQ*NZVVyCYI?DEdeEEpeCL-nLdCVW*MWCxRtEAWBXvBcCxsTVllmi-mTIbpxIbDqIix??XPRfbyzGavkHt?NClHibXSSLtylNGMsrVqiAXqQHWr*-CGerqqhvY-mnBgEct?ErRQbSlSq*QihCIKfPCT?gBNVIkbzaIG*mXDpKKAhfcyFXCsg?WDLSgyDfWyWR*yeRhYXxVVYVZZkWWMb?IFtgnNrmFiFGv?PWhevqVpQQDKZkzAMsnyMtL*iTnpyfMqpGq-hHAkeaLlnlEVcnLQKD*pENXDKraV*pKqFMVtQMigTkgvYhIK",
                "PkSKIPpeNRheawRmfkRtsBxXhkHXGKYdgFEtmkyX?GeXzzacbYaFytKwngtkvFymkwBpnRycWAtx-Iz?SNPfdwrDDiRDlCfXgmHyH-Y-tKZiTGKTzARlxdbH*XQ*P*fciQZM?mPKKNLLiBZHwKmQMBTIYvWmhIr-LZ?AaQAXMlgfPElDWizAfpypvakeQCbEtRgRcKMWNPqHNeTdTiwl?ZTwflCEFelxTVZLtZRxZF*bKXa?Ag?Zk-tATbyHFenAgwcN*wKbNniCP-n-StrFzgkmQhyPEhrFly-XQiCycZlPFilriILeISChLRHqHFiPx*hGqnAEHmKWMVRCMshQGaYwgnybeEC?embKD-sgXazxgyYpQzppQeKvFlrktW-vTwlNDrYSbpMydlDyp?AVYmPSzDtPigi*tFnfTmAQbI-kSty?MqmtBTeBADKPpkbNwdhyYTaHSHL-LiIScqAK*CmQctsnhyhcawfekvvElYBQ*ADLq?KCcXLglxckeQwl",
                )
        FixedStateAlphabetCharacterTestChecker.create_class_fixtures(
                cls,
                matrix_type=charmatrixmodel.ProteinCharacterMatrix,
                state_alphabet=charstatemodel.PROTEIN_STATE_ALPHABET,
                seq_symbols=seq_symbols,
                labels=labels)

class ContinuousTestChecker(CharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        states_lists = (
                (-231.6391 ,  +972.4189 ,  +626.6717 ,  -328.6811 ,  -213.5738 ,  +464.3897 ,  -91.3483  ,  +349.8176 ,  +333.4800 ,  +521.4970 ,  -371.4108 ,  -821.4290 ,  -86.9872  ,  -804.4891 ,  +275.3547),
                (+104.4199 ,  +669.7402 ,  -68.6082  ,  +975.4302 ,  -874.4510 ,  -191.3305 ,  -179.8437 ,  +655.5611 ,  -657.4532 ,  -563.7863 ,  +39.0321  ,  +317.0017 ,  +887.7048 ,  +342.7651 ,  -184.7631),
                (-613.2947 ,  -600.7053 ,  -700.5140 ,  +438.6092 ,  +615.5268 ,  +640.7933 ,  +503.8948 ,  -159.7922 ,  +866.8036 ,  +274.0275 ,  +462.5738 ,  -506.4329 ,  -445.4251 ,  -343.7987 ,  -285.2830),
                (-654.7695 ,  +103.3806 ,  -971.8866 ,  +853.9164 ,  +653.5797 ,  +823.6672 ,  -476.6859 ,  +325.9331 ,  +456.0902 ,  -399.7095 ,  -930.6770 ,  +762.7456 ,  +851.4525 ,  +66.4253  ,  +914.3272),
                (-762.4904 ,  +808.3665 ,  +522.5775 ,  +250.6523 ,  -287.9786 ,  -995.4612 ,  +571.9263 ,  -793.3975 ,  -42.7027  ,  +186.8869 ,  -1.5874   ,  -758.0643 ,  -69.9948  ,  -395.9015 ,  -109.9725),
                (+399.7389 ,  +31.6152  ,  +372.6323 ,  -573.5724 ,  -505.0045 ,  -375.2316 ,  +454.9046 ,  -217.4422 ,  -434.9173 ,  -454.7752 ,  -597.9571 ,  -47.6864  ,  +326.4957 ,  +545.6246 ,  +437.5032),
                (+988.9931 ,  -654.4159 ,  -767.1182 ,  -91.8658  ,  +588.7146 ,  +184.9196 ,  +115.9319 ,  -52.5935  ,  -418.1644 ,  +633.3638 ,  +736.8064 ,  -967.5157 ,  -107.9049 ,  +352.9680 ,  -70.6195),
                (-962.4392 ,  -453.7332 ,  -451.0608 ,  +341.5584 ,  +394.6056 ,  -923.0757 ,  -746.9843 ,  -965.4329 ,  -947.0617 ,  -773.1573 ,  +730.2412 ,  +375.6009 ,  +915.0743 ,  +359.0937 ,  -399.3825),
                (-296.7620 ,  +410.3270 ,  -350.8748 ,  -780.7118 ,  -175.8489 ,  +309.5648 ,  +423.9918 ,  +969.7167 ,  -244.8730 ,  -373.5084 ,  -604.9683 ,  -897.4527 ,  -534.2310 ,  -281.9905 ,  -869.7215),
                (+911.3894 ,  -989.9314 ,  +749.8109 ,  +137.4197 ,  -586.6819 ,  -153.6740 ,  -380.2316 ,  +916.4085 ,  -999.7195 ,  -893.5339 ,  -89.0466  ,  -35.9522  ,  -90.6951  ,  +776.4369 ,  -537.6624),
                (+557.6844 ,  +582.3857 ,  -284.0375 ,  -93.6579  ,  +832.2125 ,  +417.9069 ,  +597.3379 ,  +617.9041 ,  +256.8092 ,  -174.4282 ,  +489.6727 ,  +657.5248 ,  +216.8542 ,  -751.1718 ,  -398.9243),
                (-520.3224 ,  +76.4691  ,  +440.5944 ,  -55.3567  ,  -53.8396  ,  -936.0467 ,  -609.4023 ,  -120.8625 ,  +679.0885 ,  -647.1835 ,  -864.9598 ,  -470.3704 ,  +153.9019 ,  +844.2103 ,  +567.2804),
                (-332.0693 ,  +332.9506 ,  -955.3792 ,  -988.5815 ,  -603.7497 ,  +691.1748 ,  -20.4197  ,  +436.6508 ,  +31.5305  ,  -791.6890 ,  -796.3053 ,  +309.7812 ,  -138.6499 ,  +866.3018 ,  -657.6447),
                (-501.6849 ,  +249.9836 ,  +389.5901 ,  -239.3963 ,  -701.4142 ,  -25.8607  ,  -56.0275  ,  +531.3697 ,  -133.7667 ,  -973.9617 ,  +480.1729 ,  +776.4012 ,  +413.0529 ,  -456.9224 ,  -772.2399),
                (-708.9703 ,  +374.4682 ,  +688.2557 ,  -818.8122 ,  +111.4564 ,  -770.8261 ,  -838.9334 ,  -483.0598 ,  +335.7136 ,  +650.4290 ,  -957.3401 ,  -773.5307 ,  +539.6006 ,  +321.6839 ,  +366.9738),
                )
        CharacterTestChecker.create_class_fixtures(
                cls,
                matrix_type=charmatrixmodel.ContinuousCharacterMatrix,
                states_lists=states_lists,
                labels=labels)


