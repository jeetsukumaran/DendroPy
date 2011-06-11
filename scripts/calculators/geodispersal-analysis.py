#!/usr/bin/env python
#!/usr/bin/env python
import sys
import os
import subprocess
from dendropy.utility.messaging import get_logger
from dendropy.treecalc import fitch_down_pass, fitch_up_pass
from dendropy import DataSet
from dendropy.utility.error import DataParseError
from dendropy.utility.textutils import escape_nexus_token
_DEBUGGING = True
_LOG = get_logger('geodispersal')
verbose = False
AREA_NAME_LIST = []
col_width = 17

def warn(msg):
    _LOG.warn(msg)
LAST_COMMAND = ''
def write_as_nexus(stream, patterns, label):
    global LAST_COMMAND
    stream.write("\n[!%s ]\n" % label)
    p = patterns[0]
    num_chars = len(patterns)
    num_areas = len(p)
    if num_areas < len(AREA_NAME_LIST):
        warn('%d labels were found in the labels file, but only %d areas were read in the input NEXUS files' % (
                    len(AREA_NAME_LIST),
                    num_areas))
    elif num_areas > len(AREA_NAME_LIST):
        warn('Only %d labels were found in the labels file, but %d areas were read in the input NEXUS files' % (
                    len(AREA_NAME_LIST),
                    num_areas))
    stream.write("""Begin Data;
    Dimensions ntax = %d nchar = %d;
    Format datatype=standard symbols="012" ;
    Matrix \n""" % (num_areas, num_chars))
    for area_ind in range(num_areas):
        if AREA_NAME_LIST:
            try:
                name = AREA_NAME_LIST[area_ind]
            except:
                name = 'unlabelled area ' + str(1 + area_ind)
        else:
            name = "area%d" % area_ind
        name = escape_nexus_token(name)
        padding = ' ' * (col_width - len(name))
        stream.write("%s%s" % (name, padding))
        for p in patterns:
            stream.write(" %s" % str(p[area_ind]))
        stream.write("\n")
    if num_areas > 11:
        search = "HSearch"
    else:
        search = "BAndB"
    stream.write(""";
End;
Begin Paup;
    typeset * o = ord : 1-. ;
    %s ;
    SaveTrees from =1 to=100 file = %s.tre;
    %s
End;
""" % (search, label, LAST_COMMAND))
    stream.flush()


def vicariance_patterns(node_list, num_areas):
    '''Returns a transposed (characters as rows) matrix of 0, 1, 2 values for the
    vicariance matrix'''

    patterns = []
    root = node_list[0]
    assert(root.parent_node is None)
    curr_pat = [0]*num_areas
    for i in root.state_sets[0]:
        curr_pat[i] = 1
    patterns.append(curr_pat)
    areas_seen = set()
    count = 0
    for node in node_list[1:]:
        curr_pat = [0]*num_areas
        p = node.parent_node
        node.biogeo_number = count
        assert(p)
        par_area = p.state_sets[0]
        child_area = node.state_sets[0]
        areas_seen.update(par_area)
        areas_seen.update(child_area)

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

    if ABSENT_FOR_ALLTAXA_CODE != 0:
        replace_unseen_areas(patterns, num_areas, areas_seen)
    return patterns

def replace_unseen_areas(patterns, num_areas, areas_seen):
    for area_ind in range(num_areas):
        if area_ind not in areas_seen:
            for pattern in patterns:
                pattern[area_ind] = ABSENT_FOR_ALLTAXA_CODE

def dispersal_patterns(node_list, num_areas):
    '''Returns a transpose (characters as rows) matrix of 0, 1, 2 values for the
    dispersal matrix'''
    patterns = []
    root = node_list[0]
    assert(root.parent_node is None)
    curr_pat = [0]*num_areas
    for i in root.state_sets[0]:
        curr_pat[i] = 1
    patterns.append(curr_pat)

    areas_seen = set()

    for node in node_list[1:]:
        p = node.parent_node
        assert(p)
        par_area = p.state_sets[0]
        child_area = node.state_sets[0]

        areas_seen.update(par_area)
        areas_seen.update(child_area)

        twos = child_area - par_area
        curr_pat = [0]*num_areas
        for i in twos:
            curr_pat[i] = 2
        for i in par_area: #BRUCE's EMAIL mod...  par_area instead of: set.intersection(child_area, par_area):        for i in set.intersection(child_area, par_area):
            curr_pat[i] = 1
        _LOG.debug("dispersal: count = %d Par = %s, Des = %s, twos = %s, pattern= %s" % 
                    (node.biogeo_number, str(par_area), str(child_area), str(twos), str(curr_pat)))
        patterns.append(curr_pat)

    if ABSENT_FOR_ALLTAXA_CODE != 0:
        replace_unseen_areas(patterns, num_areas, areas_seen)
    return patterns

ABSENT_FOR_ALLTAXA_CODE = 0

def add_to_full(mat, full_mat, mat_el_len_list):
    prev_n_char = len(full_mat)
    prev_n_areas = prev_n_char > 0 and len(full_mat[0]) or 0
    n_additional_char = len(mat)
    assert(n_additional_char > 0)
    n_areas_mat = len(mat[0])
    if n_areas_mat > prev_n_areas:  
        for el in full_mat:
            el.extend([ABSENT_FOR_ALLTAXA_CODE]*(n_areas_mat - prev_n_areas))
    elif n_areas_mat < prev_n_areas:  
        for el in mat:
            el.extend([ABSENT_FOR_ALLTAXA_CODE]*(prev_n_areas - n_areas_mat))

    mat_el_len_list.append(n_additional_char)
    full_mat.extend(mat)

def check_file_overwrite(fn):
    if os.path.exists(fn):
        c = raw_input(fn  + ' exists. Are you sure that you would like to over-write it? (y/n)\n')
        if c.lower() != 'y':
            sys.exit('Did not get "y" in response to prompt to overwrite file. Exiting...')

if __name__ == '__main__':
    from optparse import OptionParser
    usage = '''geodispersal-analysis.py is a script to assist in analyses that use the
modified Brooks Parsimony Analysis of Lieberman and Eldredge (1996).

The input is a NEXUS file with a data matrix and a single tree. The matrix 
should be coded with in the standard datatype, and should countain a single 
character. The different states of the character represent the different 
geographical regions in the analysis. The {}-braces are used to indicate the 
occurrence of a taxon in mulitple regions.  The tree depicts the evolutionary
relationships for the taxa.

For instance:

##############################################
#NEXUS
BEGIN DATA;
	DIMENSIONS ntax = 4 nchar = 1;
	FORMAT datatype = standard symbols = "01234";
MATRIX
H_sapiens  {01234}
P_bonobo   0
G_gorilla  0
P_pygmaeus 2
;
END;
BEGIN TREES;
    TREE one = [&R] (P_pygmaeus,(G_gorilla,(P_bonobo,H_sapiens)));
END;
##############################################

is a valid (but uninteresting) input file.


'''
    parser = OptionParser(usage=usage)
    parser.add_option("--vicariance",
                      dest="vicariance",
                      default=None,
                      help="The name of the output file for the vicariance matrix")
    parser.add_option("--dispersal", 
                      dest="dispersal",
                      default=None,
                      help="The name of the output file for the dispersal matrix")
    parser.add_option("--labels", 
                      dest="labels",
                      default=None,
                      help="Name of an (optional) file with labels for the areas.  The format for this file is simply one label per line.  Note that the labels should correspond to that the state codes are listed in the SYMBOLS option of the FORMAT command input files.")
    parser.add_option('-p',
                      '--paup', 
                      dest="paup", 
                      default=False, 
                      action="store_true",
                      help="If specified, then PAUP* will be invoked to produce Vicariance.tre and Dispersal.tre Note that both the --vicariance and --dispersal options must also be used if you use the --paup option")
    parser.add_option('--absent-as-missing', 
                      dest="absent_is_missing", 
                      default=False, 
                      action="store_true",
                      help="If specified, and an area is absent for all taxa in one of the input analyses then the area will be coded as missing (?).  If this option is not used, then the area will be assigend a 0 for all characters generated from that analysis.")

    
    (options, args) = parser.parse_args()
    if len(args) == 0:
        sys.exit("Expecting a NEXUS filename as an argument")
    if options.dispersal:
        if os.path.exists(options.dispersal):
            check_file_overwrite(options.dispersal)
        disp_stream = open(options.dispersal, 'w')        
    else:
        disp_stream = sys.stdout
    if options.vicariance:
        if options.vicariance == options.dispersal:
            vic_stream = disp_stream
        else:
            check_file_overwrite(options.vicariance)
            vic_stream = open(options.vicariance, 'w')        
    else:
        vic_stream = sys.stdout
    if options.absent_is_missing:
        ABSENT_FOR_ALLTAXA_CODE = '?'
    if options.paup:
        if vic_stream == sys.stdout or disp_stream == sys.stdout:
            sys.exit("The --paup option can only be used if filenames are supplied using --vicariance and --dispersal")
        if options.vicariance == options.dispersal:
            sys.exit("The filenames are supplied in --vicariance and --dispersal must be distinct when using the --paup option")
        check_file_overwrite('Vicariance.tre')
        check_file_overwrite('Dispersal.tre')
        LAST_COMMAND = ' Quit; '

    if options.labels:
        labels_filename = options.labels
        try:
            raw_labels = [i.strip() for i in open(labels_filename, "rU").readlines()]
        except:
            sys.exit('Error reading the area labels file "%s"\n' % labels_filename)
        blank_found = False
        for label in raw_labels:
            if not label:
                blank_found = True
            else:
                if blank_found:
                    sys.exit('Blank line (or line with only whitespace) found in "%s"\nThis is not supported' % labels_filename)
                AREA_NAME_LIST.append(label)

    char_index = 0
    char_mat_index = 0
    taxon_set_index = 0
    tree_list_index = 0
    tree_index = 0

    full_vp = []
    full_dp = []
    vp_el_len = []
    dp_el_len = []
    full_fund_states = []
    try:
        for f in args:
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
            ffs = state_alphabet.fundamental_states()
            ffs_label = [str(s) for s in ffs]
            num_areas = len(ffs)
            if num_areas > len(full_fund_states):
                if ffs_label[:len(full_fund_states)] != full_fund_states:
                    sys.exit('States in file "%s" are not in the same order as states in the previous files. This will lead to incorrect coding.\nPlease edit the files so that the ordering of states is consistent in the symbols list.\n(%s != %s)' % (f, str(ffs_label), str(full_fund_states)))
                full_fund_states = ffs_label
            elif num_areas < len(full_fund_states):
                if full_fund_states[:len(ffs_label)] != ffs_label:
                    sys.exit('States in file "%s" are not in the same order as states in the previous files. This will lead to incorrect coding.\nPlease edit the files so that the ordering of states is consistent in the symbols list.\n(%s != %s)'  % (f, str(ffs_label), str(full_fund_states)))
                
            vp = vicariance_patterns(node_list, num_areas)
            add_to_full(vp, full_vp, vp_el_len)
            dp = dispersal_patterns(node_list, num_areas)
            add_to_full(dp, full_dp, dp_el_len)
    except Exception as x:
        if _DEBUGGING:
            raise
        sys.exit(str(x))


    vic_stream.write("#NEXUS\n")
    write_as_nexus(vic_stream, full_vp, "Vicariance")
    
    if vic_stream == disp_stream:
        disp_stream.write("\n\n\n\n")
    else:
        disp_stream.write("#NEXUS\n")
    
    write_as_nexus(disp_stream, full_dp, "Dispersal")
    
    if options.paup:
        try:
            vic_stream.close()
            disp_stream.close()
            rc = subprocess.call(['paup', '-n', options.vicariance])
            rc = subprocess.call(['paup', '-n', options.dispersal])
        except:
            sys.exit("Error running paup! (A python exception in subprocess).\nYou must have PAUP* in the PATH variable in your shell's environment to use the --paup option.\n") 
        if rc != 0:
            sys.exit("Error running paup! (PAUP exited with code = %d)\n" % rc)
        
        
        


