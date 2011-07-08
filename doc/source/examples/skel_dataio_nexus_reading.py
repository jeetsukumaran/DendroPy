# Reading a mixed (trees and characters) NEXUS data source, and you want all
# the trees *and* all the characters.
data = dendropy.DataSet.get_from_path(
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
        hyphens_as_tokens=False)

# Reading a tree-only NEXUS or NEWICK data source, or you want just the trees
# from a mixed NEXUS data source. For NEWICK format data, replace "nexus" with
# "newick" below. You can use ``collection_offset`` to specify a particular
# tree block to read (integer; first block offset = 0), and ``tree_offset`` to
# skip to a particular tree (integer; offset of 0 skips no trees, offset of 1
# skips the first tree, etc.)
trees = dendropy.TreeList.get_from_path(
        "data.nex",
        "nexus",
        taxon_set=None,
        exclude_trees=False,
        exclude_chars=True,
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
        hyphens_as_tokens=False,
        collection_offset=0,
        tree_offset=0)

# Reading character-only NEXUS data source, or you want just characters from a
# mixed data source. You can use ``matrix_offset` to specify a particular
# character block to read (integer; first matrix offset = 0).
dna = dendropy.DnaCharacterSet.get_from_path(
        "data.nex",
        "nexus",
        taxon_set=None,
        exclude_trees=True,
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
        hyphens_as_tokens=False,
        matrix_offset=0)
