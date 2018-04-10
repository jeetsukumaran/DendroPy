d = dendropy.StandardCharacterMatrix.get(
        path="data.fas",
        schema="fasta",
        label=None,
        taxon_namespace=None,
        matrix_offset=None,
        ignore_unrecognized_keyword_arguments=False,
        default_state_alphabet=None,
        )

d = dendropy.StandardCharacterMatrix.get(
        path="data.fas",
        schema="fasta",
        label=None,
        taxon_namespace=None,
        matrix_offset=None,
        ignore_unrecognized_keyword_arguments=False,
        default_state_alphabet=dendropy.new_standard_state_alphabet("0123456789"),
        )
