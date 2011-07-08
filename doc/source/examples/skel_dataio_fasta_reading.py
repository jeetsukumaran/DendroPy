
# When explicitly reading into
# a CharacterMatrix of a particular
# type, the `data_type` argument
# is not needed.
d = dendropy.DnaCharacterMatrix.get_from_path(
        "data.fas",
        "fasta",
        taxon_set=None,
        row_type='rich')

# Otherwise ...
d = dendropy.DataSet.get_from_path(
        "data.fas",
        "fasta",
        taxon_set=None,
        data_type="dna",
        row_type='rich')
# Or ..
d = dendropy.DataSet.get_from_path(
        "data.fas",
        "fasta",
        taxon_set=None,
        char_matrix_type=dendropy.DnaCharacterMatrix,
        row_type='rich')

