#! /usr/bin/env python

import collections
from dendropy.datamodel import charstatemodel

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
                    "Red",
                    "Blue",
                    "Green",
                    "White",
                    "Black")
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
                check_cell_annotations=check_cell_annotations,)

class Standard01234TestChecker(GenericDiscreteCharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        seq_symbols = (
                "00-10000000000-00000110-112031101000000000100000012100001000?00101010010011000011?000101100410000010000010010000111033",
                "03000100?00001100011111001110001111011000001111000011130201002011111001100001101111111011011101111?????010?0???????0?1",
                "01000000?440100110001111100100011111?10000011010?0010010101012110111001100001100011111011011100011??????????0??0111??1",
                "030000000000011000111110011100011110?10000011110?0010011111032010111001100301101011111011011100111??????10???????????1",
                "01000000?000011000011110011103011110?10000011110?0010010111302110111011??????????1111101101?10?111???????????????????1",
                )
        GenericDiscreteCharacterTestChecker.build(
                state_alphabet_fundamental_symbols="01234-",
                seq_symbols=seq_symbols,
                labels=labels)

class FixedStateAlphabetCharacterTestChecker(CharacterTestChecker):

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
                labels=labels)

class DnaTestChecker(FixedStateAlphabetCharacterTestChecker):

    @classmethod
    def build(cls, labels=None):
        seq_symbols = (
            #  0         1         2         3         4         5
            #  012345678901234567890123456789012345678901234567890
               "ACGT-NRYMWSKVHDB?ACGT-NRYMWSKVHDB?AXACGT-NRWSKVHDA",
               "ACGT-NRYMWSKT-NVHDB?XXACRYMWSKGTVHAB?ACG-NRWSKVHDC",
               "T-NRYMWSKVHDCGTXACGT-NRW-NRYMXWSB?AKVHDAKVHDB?ACGG",
               "NKVHDB?ACGT-NXXACGT-NRYMACGT-NWSKVADB?DB?XNXACGT-T",
               "RYMWSRYMWSKVHDB?ACGT-NRWSKVHDRYMWSAVHCGTA-RWSKVHD-",
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


