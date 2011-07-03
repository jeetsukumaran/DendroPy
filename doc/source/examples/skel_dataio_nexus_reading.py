trees = dendropy.TreeList.get_from_path(
        "data.nex",
        "nexus",
        taxon_set=None,
        exclude_trees=False,
        exclude_chars=False,
        as_rooted=None,
        as_unrooted=None,
        default_as_rooted=None,
        default_as_unrooted=None,
        edge_len_type=float,
        extract_comment_metadata=False,
        store_tree_weights=False,
        case_sensitive_taxon_labels=False,
        preserve_underscores=False,
        suppress_internal_node_taxa=False,
        allow_duplicate_taxon_labels=False,
        hyphens_as_tokens=False
        )

