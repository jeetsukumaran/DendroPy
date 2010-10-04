#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.txt" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

import sys
import re
import dendropy
from dendropy.test.support import pathmap
from dendropy.test.support import datagen

def main():
    tlist = dendropy.TreeList(stream=pathmap.tree_source_stream(datagen.reference_trees_filename("nexus")), schema="nexus")
    for idx, t in enumerate(tlist):
        t.label = "Tree%02d" % (idx+1)
        t.assign_node_labels_from_taxon_or_oid()

    result = []

    result.append("def reference_tree_list_postorder_node_labels():")
    result.append("    return [")
    for t in tlist:
        nodes = [("'" + nd.label + "'") for nd in t.postorder_node_iter()]
        result.append("        [%s]," % (",".join(nodes)))
    result.append("    ]")
    result.append("")

    result.append("def reference_tree_list_newick_string(taxon_set=None):")
    result.append('    return """\\')
    for t in tlist:
        result.append('        %s;' % t.as_newick_string(include_internal_labels=True))
    result.append('    """')
    result.append("")

    result.append("def reference_tree_list_node_relationships():")
    result.append("    treelist_node_references = [")
    for t in tlist:
        result.append("        {")
        for nd in t:
            if nd.parent_node is not None:
                pn_label = "'" + nd.parent_node.label + "'"
            else:
                pn_label = 'None'
            if nd.taxon is not None:
                t_label = "'" + nd.taxon.label + "'"
            else:
                t_label = 'None'
            result.append("            '%s' : NodeRelationship(parent_label=%s, child_labels=[%s], edge_length=%s, taxon_label=%s)," % \
                (nd.label,
                 pn_label,
                 ",".join(["'"+c.label+"'" for c in nd.child_nodes()]),
                 nd.edge.length,
                 t_label))
        result.append("        },")
    result.append("    ]")
    result.append("    return treelist_node_references")
    result.append("")

    tree_list_name = 'tree_list'
    src_lines = tlist.as_python_source(tree_list_name=tree_list_name, oids=True).split("\n")
    result.append("def reference_tree_list(taxon_set=None):")
    src_lines[0] = re.sub("(dendropy\.TreeList\(label\=.*)\)", "\g<1>, taxon_set=taxon_set)", src_lines[0])
    for s in src_lines:
        result.append("    %s" % s)
    result.append("""\

    # set labels of nodes with taxa to taxon label, else oid (for consistent
    # identification in debugging)
    for t in %s:
        t.assign_node_labels_from_taxon_or_oid()

    return %s
    """ % (tree_list_name, tree_list_name))
    result.append("")

    sys.stdout.write("\n".join(result))

if __name__ == "__main__":
    main()
