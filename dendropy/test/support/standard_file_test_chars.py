#! /usr/bin/env python

import collections
from dendropy.datamodel import charstatemodel

def general_verify_taxa(# {{{
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
    if char_matrix.data_type_name != "continuous":
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
    test_case.assertEqual(len(s1), len(s2))
    for c1, c2 in zip(s1, s2):
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
        s2 = test_case.label_sequence_map[label]
        general_verify_sequences(
            test_case=test_case,
            char_matrix=char_matrix,
            checker_reference_class=checker_reference_class,
            s1=s1,
            s2=s2,
            check_sequence_annotations=check_sequence_annotations,
            check_cell_annotations=check_cell_annotations)

class CharacterTestChecker(object):

    @classmethod
    def build(cls, states_lists, labels=None):
        cls.build_labels(labels=labels)
        cls.build_label_sequence_map(states_lists=states_lists)

    @classmethod
    def build_label_sequence_map(cls, states_lists):
        assert len(cls.labels) == len(states_lists)
        cls.label_sequence_map = collections.OrderedDict()
        for label, ss in zip(cls.labels, states_lists):
            cls.label_sequence_map[label] = ss

    @classmethod
    def build_labels(cls, labels=None):
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
        char_matrix = matrix_type.get_from_path(
                src_filepath,
                schema,
                **factory_kwargs)
        self.verify_char_matrix(char_matrix,
            check_taxon_annotations=check_taxon_annotations,
            check_matrix_annotations=check_matrix_annotations,
            check_sequence_annotations=check_sequence_annotations,
            check_column_annotations=check_column_annotations,
            check_cell_annotations=check_cell_annotations)

class GenericDiscreteCharacterTestChecker(CharacterTestChecker):

    @classmethod
    def build(cls,
            state_alphabet_fundamental_symbols,
            seq_symbols,
            labels=None):
        cls.state_alphabet_fundamental_symbols = list(state_alphabet_fundamental_symbols)
        CharacterTestChecker.build_labels(labels=labels)
        cls.states_symbols_lists = []
        for ss in seq_symbols:
            cls.states_symbols_lists.append(list(ss))

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
        states = []
        for ss in self.__class__.states_symbols_lists:
            states.append(char_matrix.default_state_alphabet.get_states_for_symbols(ss))
        CharacterTestChecker.build_label_sequence_map(states_lists=states)
        general_char_matrix_checker(self,
                char_matrix,
                self.__class__,
                check_taxon_annotations=check_taxon_annotations,
                check_matrix_annotations=check_matrix_annotations,
                check_sequence_annotations=check_sequence_annotations,
                check_column_annotations=check_column_annotations,
                check_cell_annotations=check_cell_annotations,)# }}}

class Standard01234TestChecker(GenericDiscreteCharacterTestChecker):# {{{

    @classmethod
    def build(cls, labels=None):
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
        GenericDiscreteCharacterTestChecker.build(
                state_alphabet_fundamental_symbols="01234-",
                seq_symbols=seq_symbols,
                labels=labels)# }}}

class FixedStateAlphabetCharacterTestChecker(CharacterTestChecker):# {{{

    @classmethod
    def build(cls,
            state_alphabet,
            seq_symbols,
            labels=None):
        cls.state_alphabet = state_alphabet
        states_lists = []
        for ss in seq_symbols:
            seq_states = tuple(cls.state_alphabet.get_states_for_symbols(ss))
            states_lists.append(seq_states)
        CharacterTestChecker.build(
                states_lists=states_lists,
                labels=labels)# }}}

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
        FixedStateAlphabetCharacterTestChecker.build(
                state_alphabet=charstatemodel.DNA_STATE_ALPHABET,
                seq_symbols=seq_symbols,
                labels=labels)

class RnaTestChecker(FixedStateAlphabetCharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        seq_symbols = (
               "ACGU-NRYMWSKVHDB?ACGU-NRYMWSKVHDB?AXACGU-NRWSKVHDA",
               "ACGU-NRYMWSKU-NVHDB?XXACRYMWSKGUVHAB?ACG-NRWSKVHDC",
               "U-NRYMWSKVHDCGUXACGU-NRW-NRYMXWSB?AKVHDAKVHDB?ACGG",
               "NKVHDB?ACGU-NXXACGU-NRYMACGU-NWSKVADB?DB?XNXACGU-U",
               "RYMWSRYMWSKVHDB?ACGU-NRWSKVHDRYMWSAVHCGUA-RWSKVHD-",
                )
        FixedStateAlphabetCharacterTestChecker.build(
                state_alphabet=charstatemodel.RNA_STATE_ALPHABET,
                seq_symbols=seq_symbols,
                labels=labels)

class ProteinTestChecker(FixedStateAlphabetCharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        seq_symbols = (
               "Y?LTDCQWD-RIISDM*WZ?PTYBEWAT*LW?NYFILESB?VR?CAWYP-CXG?-FHVHZXTEKDW?R-FGES*KA-FZCGNX*LMVWZGFZMHEQPSD???P*RNLTPRSAINMGVQKBNMYIXBXKHACWWWHBKQQVV",
               "Y?LTDCQWD-QKBNMYIXBR-FGES*KA-FZCGNX*LMVWZGFZMHEQPSDRIISDM*WZ?PTRNLTPRSAINMGV???P*VXKHACWWWHBKQQVYBEWAT*LW?NYFILESB?VR?CAWYP-CXG?-FHVHZXTEKDW?",
               "Y?LTDCQWD-HVHZXTEKDW?R-FGES*KA-FW?NYFILESB?VR?CAAINRIISDM*WZ?PTYBEWYP-CXG?-FMGVQKBNMYIXBXKHACWWWHBKQQVZCGNX*LMVWZGFZMHEQPSD???P*RNLTPRSWAT*LV",
               "KHACWWWHBKAT*LW?NYFILESB?VR?CAWYP-CXG?-FHVHZXTEKDW?QQVDCQWD-RIISDM*WZ?PTYBEWR-FGES*KA-FZCGNX*LMVWZGY?LTFZMHEQPSD???P*RNLTPRSAINMGVQKBNMYIXBVV",
               "YKDW?R-FGES*KA-FZCGNX*LMVWZGFZMHEQPSD???P*RNLTPRSAINMGVQKBNMYIXBXKHACWWWHBKQQV?LTDCQWD-RIISDM*WZ?PTYBEWAT*LW?NYFILESB?VR?CAWYP-CXG?-FHVHZXTVV",
                )
        FixedStateAlphabetCharacterTestChecker.build(
                state_alphabet=charstatemodel.PROTEIN_STATE_ALPHABET,
                seq_symbols=seq_symbols,
                labels=labels)

class ContinuousTestChecker(CharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        states_lists = (
                    (-0.0230 , -0.3273 , -13.4836 , +0.0868 , -0.3942 , -0.0752 , -1020.0652),
                    (-1.3153 , +1.2542 , +14.4277 , +1.1204 , -0.1048 , -0.0473 , -2210.0548),
                    (-2.1870 , +2.0093 , +15.4810 , +0.6872 , -0.1293 , +0.3246 , +10301.2798),
                    (+3.1097 , +3.1443 , +16.2746 , +1.0220 , -0.0353 , +1.0085 , +20.1666),
                    (+4.4274 , +4.4127 , -17.1721 , +0.7207 , -1.2134 , +0.4868 , -0.1480),
                )
        CharacterTestChecker.build(
                states_lists=states_lists,
                labels=labels)


