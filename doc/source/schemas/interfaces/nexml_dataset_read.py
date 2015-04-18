data_set = dendropy.DataSet()
data_set.read(
    path="path/to/file",
    schema="nexml",
    exclude_chars=False,
    exclude_trees=False,
    default_namespace="http://www.nexml.org/2009",
    case_sensitive_taxon_labels=False,
    ignore_unrecognized_keyword_arguments=False,
    )



