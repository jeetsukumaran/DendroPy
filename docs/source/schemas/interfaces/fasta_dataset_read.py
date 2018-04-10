d = dendropy.DataSet()
d.read(
        path="data.fas",
        schema="fasta",
        label=None,
        taxon_namespace=None,
        ignore_unrecognized_keyword_arguments=False,
        data_type="dna",
        )

d = dendropy.DataSet()
d.read(
        path="data.fas",
        schema="fasta",
        label=None,
        taxon_namespace=None,
        ignore_unrecognized_keyword_arguments=False,
        data_type="protein",
        )

d = dendropy.DataSet()
d.read(
        path="data.fas",
        schema="fasta",
        label=None,
        taxon_namespace=None,
        ignore_unrecognized_keyword_arguments=False,
        data_type="standard",
        )

d = dendropy.DataSet()
d.read(
        path="data.fas",
        schema="fasta",
        label=None,
        taxon_namespace=None,
        ignore_unrecognized_keyword_arguments=False,
        data_type="standard",
        default_state_alphabet=dendropy.new_standard_state_alphabet("abc"),
        )
