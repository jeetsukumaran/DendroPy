tree_list = dendropy.TreeList()
tree_list.read(
    path="path/to/file",
    schema="nexml",
    collection_offset=None,
    tree_offset=None,
    default_namespace="http://www.nexml.org/2009",
    case_sensitive_taxon_labels=False,
    ignore_unrecognized_keyword_arguments=False,
    )
