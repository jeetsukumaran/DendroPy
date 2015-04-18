d = dendropy.StandardCharacterMatrix.get(
        path="data.phylip",
        schema="phylip",
        label=None,
        taxon_namespace=None,
        matrix_offset=None,
        strict=False,
        interleaved=False,
        multispace_delimiter=False,
        underscore_to_spaces=False,
        ignore_invalid_chars=False,
        ignore_unrecognized_keyword_arguments=False,
        default_state_alphabet=None,
        )

d = dendropy.StandardCharacterMatrix.get(
        path="data.phylip",
        schema="phylip",
        label=None,
        taxon_namespace=None,
        matrix_offset=None,
        strict=False,
        interleaved=False,
        multispace_delimiter=False,
        underscore_to_spaces=False,
        ignore_invalid_chars=False,
        ignore_unrecognized_keyword_arguments=False,
        default_state_alphabet=dendropy.new_standard_state_alphabet("0123456789"),
        )
