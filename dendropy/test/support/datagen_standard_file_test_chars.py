#! /usr/bin/env python

import collections
from dendropy.datamodel import charstatemodel

class CharacterTestChecker(object):

    def verify_taxa(self,
            taxon_namespace,
            check_annotations=True):
        cls = self.__class__
        self.assertEqual(len(taxon_namespace), len(cls.labels))
        for taxon, label in zip(taxon_namespace, cls.labels):
            self.assertEqual(taxon.label, label)

    def verify_character_cell_states(self,
            c1, c2,
            check_annotations=True):
        raise NotImplementedError
        if check_annotations:
            pass # not yet implemented

    def verify_sequences(self,
            s1, s2,
            check_sequence_annotations=True,
            check_cell_annotations=True):
        self.assertEqual(len(s1), len(s2))
        for c1, c2 in zip(s1, s2):
            self.verify_character_cell_states(
                    c1, c2,
                    check_cell_annotations)

    def verify_char_matrix(self,
            char_matrix,
            check_taxon_annotations=True,
            check_matrix_annotations=True,
            check_sequence_annotations=True,
            check_column_annotations=True,
            check_cell_annotations=True):
        cls = self.__class__
        self.assertEqual(len(char_matrix), len(cls.labels))
        self.verify_taxa(char_matrix.taxon_namespace,
            check_annotations=check_taxon_annotations)
        for taxon, label in zip(char_matrix, cls.labels):
            self.assertEqual(taxon.label, label)
            s1 = char_matrix[taxon]
            s2 = self.label_sequence_map[label]
            self.verify_sequences(
                s1, s2,
                check_sequence_annotations=check_sequence_annotations,
                check_cell_annotations=check_cell_annotations)

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

class FixedStateAlphabetCharacterTestChecker(CharacterTestChecker):

    @classmethod
    def build(cls,
            state_alphabet,
            seq_symbols,
            labels=None):
        if labels is None:
            cls.labels = (
                    "Red",
                    "Blue",
                    "Green",
                    "White",
                    "Black",
                    )
        else:
            cls.labels = tuple(labels)
        cls.state_alphabet = state_alphabet
        cls.seq_symbols = seq_symbols
        assert len(cls.seq_symbols) == len(cls.labels)
        cls.seq_states = []
        cls.label_sequence_map = collections.OrderedDict()
        for label, ss in zip(cls.labels, cls.seq_symbols):
            seq_states = tuple(cls.state_alphabet.get_states_for_symbols(ss))
            cls.seq_states.append(seq_states)
            cls.label_sequence_map[label] = seq_states

    def verify_character_cell_states(self,
            c1, c2,
            check_annotations=True):
        self.assertIs(c1, c2)

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
            #  0         1         2         3         4         5
            #  012345678901234567890123456789012345678901234567890
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




