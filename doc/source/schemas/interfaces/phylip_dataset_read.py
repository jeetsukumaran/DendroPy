d = dendropy.DataSet()
d.read(
        path="data.phylip",
        schema="phylip",
        label=None,
        taxon_namespace=None,
        strict=False,
        interleaved=False,
        multispace_delimiter=False,
        underscore_to_spaces=False,
        ignore_invalid_chars=False,
        ignore_unrecognized_keyword_arguments=False,
        data_type="dna",
        )

d = dendropy.DataSet()
d.read(
        path="data.phylip",
        schema="phylip",
        label=None,
        taxon_namespace=None,
        strict=False,
        interleaved=False,
        multispace_delimiter=False,
        underscore_to_spaces=False,
        ignore_invalid_chars=False,
        ignore_unrecognized_keyword_arguments=False,
        data_type="protein",
        )

d = dendropy.DataSet()
d.read(
        path="data.phylip",
        schema="phylip",
        label=None,
        taxon_namespace=None,
        strict=False,
        interleaved=False,
        multispace_delimiter=False,
        underscore_to_spaces=False,
        ignore_invalid_chars=False,
        ignore_unrecognized_keyword_arguments=False,
        data_type="standard",
        )

d = dendropy.DataSet()
d.read(
        path="data.phylip",
        schema="phylip",
        label=None,
        taxon_namespace=None,
        strict=False,
        interleaved=False,
        multispace_delimiter=False,
        underscore_to_spaces=False,
        ignore_invalid_chars=False,
        ignore_unrecognized_keyword_arguments=False,
        data_type="standard",
        default_state_alphabet=dendropy.new_standard_state_alphabet("abc"),
        )
