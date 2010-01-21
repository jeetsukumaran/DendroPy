#! /usr/bin/env python

###############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2009 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

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
            result.append("            '%s' : NodeRelationship(parent_label=%s, child_labels=[%s], edge_length=%s, taxon_label=%s)," % \
                (nd.label,
                 ("'"+nd.parent_node.label+"'") if nd.parent_node is not None else 'None',
                 ",".join(["'"+c.label+"'" for c in nd.child_nodes()]),
                 nd.edge.length,
                 ("'"+nd.taxon.label+"'") if nd.taxon is not None else 'None'))
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
