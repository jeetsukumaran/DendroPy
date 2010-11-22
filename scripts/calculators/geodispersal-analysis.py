#!/usr/bin/env python
#!/usr/bin/env python
import sys
from dendropy.utility.messaging import get_logger
from dendropy.treecalc import fitch_down_pass, fitch_up_pass
from dendropy import DataSet
from dendropy.utility.error import DataParseError
_DEBUGGING = True
_LOG = get_logger('geodispersal')
verbose = False
gStateNames = None
col_width = 17
def write_as_nexus(stream, patterns, label):
    stream.write("\n[!%s ]\n" % label)
    p = patterns[0]
    num_chars = len(patterns)
    num_areas = len(p)
    stream.write("""Begin Data;
    Dimensions ntax = %d nchar = %d;
    Format datatype=standard symbols="012" ;
    Matrix \n""" % (num_areas, num_chars))
    for area_ind in range(num_areas):
        if gStateNames:
            name = gStateNames[area_ind]
        else:
            name = "area%d" % area_ind
        padding = ' ' * (col_width - len(name))
        stream.write("%s%s" % (name, padding))
        for p in patterns:
            stream.write(" %d" % p[area_ind])
        stream.write("\n")
    if num_areas > 13:
        search = "HSearch"
    else:
        search = "BAndB"
    stream.write(";\nEnd;\n\nbegin paup;\n   typeset * o = ord : 1-. ;\n    %s ;\n  SaveTrees from =1 to=100 file = %s.tre;\nend; \n" % (search, label))


def vicariance_patterns(node_list, num_areas):
    patterns = []
    root = node_list[0]
    assert(root.parent_node is None)
    curr_pat = [0]*num_areas
    for i in root.state_sets[0]:
        curr_pat[i] = 1
    patterns.append(curr_pat)

    count = 0
    for node in node_list[1:]:
        curr_pat = [0]*num_areas
        p = node.parent_node
        node.biogeo_number = count
        assert(p)
        par_area = p.state_sets[0]
        child_area = node.state_sets[0]

        diff = par_area - child_area
        _LOG.debug("count = %d Par = %s       Des = %s    diff = %s" %(count, str(par_area), str(child_area), str(diff)))
        if diff:
            # range contraction...
            # parent existed in areas that the descendant did not
            twos = set.intersection(child_area, par_area)
            ones = set.symmetric_difference(child_area, par_area)
            for i in twos:
                curr_pat[i] = 2
        else:
            ones = set.union(child_area, par_area)
        for i in ones:
            curr_pat[i] = 1
        _LOG.debug("curr_pat = %s" %(str(curr_pat)))
        patterns.append(curr_pat)
        count += 1

    return patterns

def dispersal_patterns(node_list, num_areas):
    patterns = []
    root = node_list[0]
    assert(root.parent_node is None)
    curr_pat = [0]*num_areas
    for i in root.state_sets[0]:
        curr_pat[i] = 1
    patterns.append(curr_pat)

    for node in node_list[1:]:
        p = node.parent_node
        assert(p)
        par_area = p.state_sets[0]
        child_area = node.state_sets[0]

        twos = child_area - par_area
        curr_pat = [0]*num_areas
        for i in twos:
            curr_pat[i] = 2
        for i in set.intersection(child_area, par_area):
            curr_pat[i] = 1
        _LOG.debug("dispersal: count = %d Par = %s, Des = %s, twos = %s, pattern= %s" % 
                    (node.biogeo_number, str(par_area), str(child_area), str(twos), str(curr_pat)))
        patterns.append(curr_pat)
    return patterns


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    (options, args) = parser.parse_args()
    if len(args) == 0:
        sys.exit("Expecting a filename as an argument")
    if len(args) == 2:
        labels_filename = args[1]
        gStateNames = [i.strip() for i in open(labels_filename, "rU").readlines()]

    char_index = 0
    char_mat_index = 0
    taxon_set_index = 0
    tree_list_index = 0
    tree_index = 0

    try:
        for f in args[:1]:
            fo = open(f, "rU")
            dataset = DataSet()
            try:
                dataset.read(stream=fo, schema="NEXUS")
            except DataParseError as dfe:
                raise ValueError(str(dfe))
            if len(dataset.taxon_sets) != 1:
                raise ValueError("Expecting one set of taxa in %s" % f)
            if len(dataset.tree_lists) != 1:
                raise ValueError("Expecting one tree in %s" % f)
            if len(dataset.char_matrices) != 1:
                raise ValueError("Expecting one character matrix in %s" % f)
            char_mat = dataset.char_matrices[char_mat_index]
            taxon_set = dataset.taxon_sets[taxon_set_index]
            tree = dataset.tree_lists[tree_list_index][tree_index]
            state_alphabet = char_mat.state_alphabets[char_index]

            taxon_to_state_indices = char_mat.create_taxon_to_state_set_map(char_indices=[char_index])

            if not tree.is_rooted:
                raise ValueError("Tree must be rooted")

            root = tree.seed_node
            root_children = root.child_nodes()
            if len(root_children) != 2:
                raise ValueError("Expecting a binary rooted tree.  Root has more than 2 children!")

            if len(root_children[0].child_nodes()) > 0:
                if len(root_children[1].child_nodes()) > 0:
                    raise ValueError("Expecting one of the children of the root to be a single taxon (leaf) outgroup")
                outgroup = root_children[1]
            else:
                if len(root_children[1].child_nodes()) == 0:
                    raise ValueError("A tree of more than 2 leaves is required")
                outgroup = root_children[0]

            # we get the preorder and reverse twice, rather than the post-order 
            #   and reverse it once, so that we traverse the edges left to right
            node_list = [i for i in tree.preorder_node_iter()]
            for nd in node_list:
                c = nd.child_nodes()
                if c:
                    if len(c) != 2:
                        raise ValueError("Tree must be fully resolved")
            
            node_list.reverse()
            fitch_down_pass(node_list, taxa_to_state_set_map=taxon_to_state_indices)
            node_list.reverse()
            root.state_sets = list(outgroup.state_sets)
            fitch_up_pass(node_list, taxa_to_state_set_map=taxon_to_state_indices)
            node_list.remove(outgroup)
            num_areas = len(state_alphabet.fundamental_states())
            vp = vicariance_patterns(node_list, num_areas)
            dp = dispersal_patterns(node_list, num_areas)

            sys.stdout.write("#NEXUS\n")
            write_as_nexus(sys.stdout, vp, "Vicariance")
            sys.stdout.write("\n\n\n\n")
            write_as_nexus(sys.stdout, dp, "Dispersal")
    except Exception as x:
        if _DEBUGGING:
            raise
        sys.exit(str(x))

