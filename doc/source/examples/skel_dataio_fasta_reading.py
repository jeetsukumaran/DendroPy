
# When explicitly reading into
# a CharacterMatrix of a particular
# type, the `data_type` argument
# is not needed.
d = dendropy.DnaCharacterMatrix.get_from_path(
        "data.fas",
        "fasta",
        taxon_namespace=None,
        row_type='rich')

# Otherwise ...
d = dendropy.DataSet.get_from_path(
        "data.fas",
        "fasta",
        taxon_namespace=None,
        data_type="dna",
        row_type='rich')
# Or ..
d = dendropy.DataSet.get_from_path(
        "data.fas",
        "fasta",
        taxon_namespace=None,
        char_matrix_type=dendropy.DnaCharacterMatrix,
        row_type='rich')

