#! /usr/bin/env python

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Wrapper around calls to RAxML.
"""

import sys
import os
import subprocess
import tempfile
import dendropy
import random

from dendropy.utility.messaging import ConsoleMessenger
from dendropy.utility import processio

def get_messenger(verbosity=1):
    if verbosity == 0:
        messaging_level = ConsoleMessenger.ERROR_MESSAGING_LEVEL
    else:
        messaging_level = ConsoleMessenger.INFO_MESSAGING_LEVEL
    messenger = ConsoleMessenger(name="raxml-map-bipartitions",
                    messaging_level=messaging_level)
    return messenger

###############################################################################
## RAXML WRAPPER
# raxmlHPC[-SSE3|-PTHREADS|-PTHREADS-SSE3|-HYBRID|-HYBRID-SSE3]
#   -s sequenceFileName -n outputFileName -m substitutionModel
#   [-a weightFileName] [-A secondaryStructureSubstModel]
#   [-b bootstrapRandomNumberSeed] [-B wcCriterionThreshold]
#   [-c numberOfCategories] [-C] [-d] [-D]
#   [-e likelihoodEpsilon] [-E excludeFileName]
#   [-f a|A|b|c|d|e|E|F|g|h|i|I|j|J|m|n|o|p|r|s|S|t|u|v|w|x|y] [-F]
#   [-g groupingFileName] [-G placementThreshold] [-h]
#   [-i initialRearrangementSetting] [-I autoFC|autoMR|autoMRE|autoMRE_IGN]
#   [-j] [-J MR|MR_DROP|MRE|STRICT|STRICT_DROP] [-k] [-K] [-M]
#   [-o outGroupName1[,outGroupName2[,...]]]
#   [-p parsimonyRandomSeed] [-P proteinModel]
#   [-q multipleModelFileName] [-r binaryConstraintTree]
#   [-R binaryModelParamFile] [-S secondaryStructureFile] [-t userStartingTree]
#   [-T numberOfThreads] [-U] [-v] [-w outputDirectory] [-W slidingWindowSize]
#   [-x rapidBootstrapRandomNumberSeed] [-X] [-y]
#   [-z multipleTreesFile] [-#|-N numberOfRuns|autoFC|autoMR|autoMRE|autoMRE_IGN]
#   -a      Specify a column weight file name to assign individual weights to each column of
#           the alignment. Those weights must be integers separated by any type and number
#           of whitespaces whithin a separate file, see file "example_weights" for an example.
#   -A      Specify one of the secondary structure substitution models implemented in RAxML.
#           The same nomenclature as in the PHASE manual is used, available models:
#           S6A, S6B, S6C, S6D, S6E, S7A, S7B, S7C, S7D, S7E, S7F, S16, S16A, S16B
#           DEFAULT: 16-state GTR model (S16)
#   -b      Specify an integer number (random seed) and turn on bootstrapping
#           DEFAULT: OFF
#   -B      specify a floating point number between 0.0 and 1.0 that will be used as cutoff threshold
#           for the MR-based bootstopping criteria. The recommended setting is 0.03.
#           DEFAULT: 0.03 (recommended empirically determined setting)
#   -c      Specify number of distinct rate catgories for RAxML when modelOfEvolution
#           is set to GTRCAT or GTRMIX
#           Individual per-site rates are categorized into numberOfCategories rate
#           categories to accelerate computations.
#           DEFAULT: 25
#   -C      Conduct model parameter optimization on gappy, partitioned multi-gene alignments with per-partition
#           branch length estimates (-M enabled) using the fast method with pointer meshes described in:
#           Stamatakis and Ott: "Efficient computation of the phylogenetic likelihood function on multi-gene alignments and multi-core processors"
#           WARNING: We can not conduct useful tree searches using this method yet! Does not work with Pthreads version.
#   -d      start ML optimization from random starting tree
#           DEFAULT: OFF
#   -D      ML search convergence criterion. This will break off ML searches if the relative
#           Robinson-Foulds distance between the trees obtained from two consecutive lazy SPR cycles
#           is smaller or equal to 1%. Usage recommended for very large datasets in terms of taxa.
#           On trees with more than 500 taxa this will yield execution time improvements of approximately 50%
#           While yielding only slightly worse trees.
#           DEFAULT: OFF
#   -e      set model optimization precision in log likelihood units for final
#           optimization of tree topology under MIX/MIXI or GAMMA/GAMMAI
#           DEFAULT: 0.1   for models not using proportion of invariant sites estimate
#                    0.001 for models using proportion of invariant sites estimate
#   -E      specify an exclude file name, that contains a specification of alignment positions you wish to exclude.
#           Format is similar to Nexus, the file shall contain entries like "100-200 300-400", to exclude a
#           single column write, e.g., "100-100", if you use a mixed model, an appropriatly adapted model file
#           will be written.
#   -f      select algorithm:
#           "-f a": rapid Bootstrap analysis and search for best-scoring ML tree in one program run
#           "-f A": compute marginal ancestral states on a ROOTED reference tree provided with "t"
#           "-f b": draw bipartition information on a tree provided with "-t" based on multiple trees
#                   (e.g., from a bootstrap) in a file specifed by "-z"
#           "-f c": check if the alignment can be properly read by RAxML
#           "-f d": new rapid hill-climbing
#                   DEFAULT: ON
#           "-f e": optimize model+branch lengths for given input tree under GAMMA/GAMMAI only
#           "-f E": execute very fast experimental tree search, at present only for testing
#           "-f F": execute fast experimental tree search, at present only for testing
#           "-f g": compute per site log Likelihoods for one ore more trees passed via
#                   "-z" and write them to a file that can be read by CONSEL
#           "-f h": compute log likelihood test (SH-test) between best tree passed via "-t"
#                   and a bunch of other trees passed via "-z"
#           "-f i": EXPERIMENTAL do not use for real tree inferences: conducts a single cycle of fast lazy SPR moves
#                   on a given input tree, to be used in combination with -C and -M
#           "-f I": EXPERIMENTAL do not use for real tree inferences: conducts a single cycle of thorough lazy SPR moves
#                   on a given input tree, to be used in combination with -C and -M
#           "-f j": generate a bunch of bootstrapped alignment files from an original alignemnt file.
#                   You need to specify a seed with "-b" and the number of replicates with "-#"
#           "-f J": Compute SH-like support values on a given tree passed via "-t".
#           "-f m": compare bipartitions between two bunches of trees passed via "-t" and "-z"
#                   respectively. This will return the Pearson correlation between all bipartitions found
#                   in the two tree files. A file called RAxML_bipartitionFrequencies.outpuFileName
#                   will be printed that contains the pair-wise bipartition frequencies of the two sets
#           "-f n": compute the log likelihood score of all trees contained in a tree file provided by
#                   "-z" under GAMMA or GAMMA+P-Invar
#           "-f o": old and slower rapid hill-climbing without heuristic cutoff
#           "-f p": perform pure stepwise MP addition of new sequences to an incomplete starting tree and exit
#           "-f r": compute pairwise Robinson-Foulds (RF) distances between all pairs of trees in a tree file passed via "-z"
#                   if the trees have node labales represented as integer support values the program will also compute two flavors of
#                   the weighted Robinson-Foulds (WRF) distance
#           "-f s": split up a multi-gene partitioned alignment into the respective subalignments
#           "-f S": compute site-specific placement bias using a leave one out test inspired by the evolutionary placement algorithm
#           "-f t": do randomized tree searches on one fixed starting tree
#           "-f u": execute morphological weight calibration using maximum likelihood, this will return a weight vector.
#                   you need to provide a morphological alignment and a reference tree via "-t"
#           "-f v": classify a bunch of environmental sequences into a reference tree using thorough read insertions
#                   you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)
#           "-f w": compute ELW test on a bunch of trees passed via "-z"
#           "-f x": compute pair-wise ML distances, ML model parameters will be estimated on an MP
#                   starting tree or a user-defined tree passed via "-t", only allowed for GAMMA-based
#                   models of rate heterogeneity
#           "-f y": classify a bunch of environmental sequences into a reference tree using parsimony
#                   you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)
#           DEFAULT for "-f": new rapid hill climbing
#   -F      enable ML tree searches under CAT model for very large trees without switching to
#           GAMMA in the end (saves memory).
#           This option can also be used with the GAMMA models in order to avoid the thorough optimization
#           of the best-scoring ML tree in the end.
#           DEFAULT: OFF
#   -g      specify the file name of a multifurcating constraint tree
#           this tree does not need to be comprehensive, i.e. must not contain all taxa
#   -G      enable the ML-based evolutionary placement algorithm heuristics
#           by specifiyng a threshold value (fraction of insertion branches to be evaluated
#           using slow insertions under ML).
#   -h      Display this help message.
#   -i      Initial rearrangement setting for the subsequent application of topological
#           changes phase
#   -I      a posteriori bootstopping analysis. Use:
#          "-I autoFC" for the frequency-based criterion
#          "-I autoMR" for the majority-rule consensus tree criterion
#          "-I autoMRE" for the extended majority-rule consensus tree criterion
#          "-I autoMRE_IGN" for metrics similar to MRE, but include bipartitions under the threshold whether they are compatible
#                           or not. This emulates MRE but is faster to compute.
#           You also need to pass a tree file containg several bootstrap replicates via "-z"
#   -j      Specifies that intermediate tree files shall be written to file during the standard ML and BS tree searches.
#           DEFAULT: OFF
#   -J      Compute majority rule consensus tree with "-J MR" or extended majority rule consensus tree with "-J MRE"
#           or strict consensus tree with "-J STRICT".
#           Options "-J STRICT_DROP" and "-J MR_DROP" will execute an algorithm that identifies dropsets which contain
#           rogue taxa as proposed by Pattengale et al. in the paper "Uncovering hidden phylogenetic consensus".
#           You will also need to provide a tree file containing several UNROOTED trees via "-z"
#   -k      Specifies that bootstrapped trees should be printed with branch lengths.
#           The bootstraps will run a bit longer, because model parameters will be optimized
#           at the end of each run under GAMMA or GAMMA+P-Invar respectively.
#           DEFAULT: OFF
#   -K      Specify one of the multi-state substitution models (max 32 states) implemented in RAxML.
#           Available models are: ORDERED, MK, GTR
#           DEFAULT: GTR model
#   -m      Model of Binary (Morphological), Nucleotide, Multi-State, or Amino Acid Substitution:
#           BINARY:
#             "-m BINCAT"         : Optimization of site-specific
#                                   evolutionary rates which are categorized into numberOfCategories distinct
#                                   rate categories for greater computational efficiency. Final tree might be evaluated
#                                   automatically under BINGAMMA, depending on the tree search option
#             "-m BINCATI"        : Optimization of site-specific
#                                   evolutionary rates which are categorized into numberOfCategories distinct
#                                   rate categories for greater computational efficiency. Final tree might be evaluated
#                                   automatically under BINGAMMAI, depending on the tree search option
#             "-m BINGAMMA"       : GAMMA model of rate
#                                   heterogeneity (alpha parameter will be estimated)
#             "-m BINGAMMAI"      : Same as BINGAMMA, but with estimate of proportion of invariable sites
#           NUCLEOTIDES:
#             "-m GTRCAT"         : GTR + Optimization of substitution rates + Optimization of site-specific
#                                   evolutionary rates which are categorized into numberOfCategories distinct
#                                   rate categories for greater computational efficiency.  Final tree might be evaluated
#                                   under GTRGAMMA, depending on the tree search option
#             "-m GTRCATI"        : GTR + Optimization of substitution rates + Optimization of site-specific
#                                   evolutionary rates which are categorized into numberOfCategories distinct
#                                   rate categories for greater computational efficiency.  Final tree might be evaluated
#                                   under GTRGAMMAI, depending on the tree search option
#             "-m GTRGAMMA"       : GTR + Optimization of substitution rates + GAMMA model of rate
#                                   heterogeneity (alpha parameter will be estimated)
#             "-m GTRGAMMAI"      : Same as GTRGAMMA, but with estimate of proportion of invariable sites
#           MULTI-STATE:
#             "-m MULTICAT"         : Optimization of site-specific
#                                   evolutionary rates which are categorized into numberOfCategories distinct
#                                   rate categories for greater computational efficiency. Final tree might be evaluated
#                                   automatically under MULTIGAMMA, depending on the tree search option
#             "-m MULTICATI"        : Optimization of site-specific
#                                   evolutionary rates which are categorized into numberOfCategories distinct
#                                   rate categories for greater computational efficiency. Final tree might be evaluated
#                                   automatically under MULTIGAMMAI, depending on the tree search option
#             "-m MULTIGAMMA"       : GAMMA model of rate
#                                   heterogeneity (alpha parameter will be estimated)
#             "-m MULTIGAMMAI"      : Same as MULTIGAMMA, but with estimate of proportion of invariable sites
#             You can use up to 32 distinct character states to encode multi-state regions, they must be used in the following order:
#             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V
#             i.e., if you have 6 distinct character states you would use 0, 1, 2, 3, 4, 5 to encode these.
#             The substitution model for the multi-state regions can be selected via the "-K" option
#           AMINO ACIDS:
#             "-m PROTCATmatrixName[F]"         : specified AA matrix + Optimization of substitution rates + Optimization of site-specific
#                                                 evolutionary rates which are categorized into numberOfCategories distinct
#                                                 rate categories for greater computational efficiency.   Final tree might be evaluated
#                                                 automatically under PROTGAMMAmatrixName[f], depending on the tree search option
#             "-m PROTCATImatrixName[F]"        : specified AA matrix + Optimization of substitution rates + Optimization of site-specific
#                                                 evolutionary rates which are categorized into numberOfCategories distinct
#                                                 rate categories for greater computational efficiency.   Final tree might be evaluated
#                                                 automatically under PROTGAMMAImatrixName[f], depending on the tree search option
#             "-m PROTGAMMAmatrixName[F]"       : specified AA matrix + Optimization of substitution rates + GAMMA model of rate
#                                                 heterogeneity (alpha parameter will be estimated)
#             "-m PROTGAMMAImatrixName[F]"      : Same as PROTGAMMAmatrixName[F], but with estimate of proportion of invariable sites
#             Available AA substitution models: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, MTART, MTZOA, PMB, HIVB, HIVW, JTTDCMUT, FLU, GTR
#             With the optional "F" appendix you can specify if you want to use empirical base frequencies
#             Please note that for mixed models you can in addition specify the per-gene AA model in
#             the mixed model file (see manual for details). Also note that if you estimate AA GTR parameters on a partitioned
#             dataset, they will be linked (estimated jointly) across all partitions to avoid over-parametrization
#   -M      Switch on estimation of individual per-partition branch lengths. Only has effect when used in combination with "-q"
#           Branch lengths for individual partitions will be printed to separate files
#           A weighted average of the branch lengths is computed by using the respective partition lengths
#           DEFAULT: OFF
#   -n      Specifies the name of the output file.
#   -o      Specify the name of a single outgrpoup or a comma-separated list of outgroups, eg "-o Rat"
#           or "-o Rat,Mouse", in case that multiple outgroups are not monophyletic the first name
#           in the list will be selected as outgroup, don't leave spaces between taxon names!
#   -p      Specify a random number seed for the parsimony inferences. This allows you to reproduce your results
#           and will help me debug the program.
#   -P      Specify the file name of a user-defined AA (Protein) substitution model. This file must contain
#           420 entries, the first 400 being the AA substitution rates (this must be a symmetric matrix) and the
#           last 20 are the empirical base frequencies
#   -q      Specify the file name which contains the assignment of models to alignment
#           partitions for multiple models of substitution. For the syntax of this file
#           please consult the manual.
#   -r      Specify the file name of a binary constraint tree.
#           this tree does not need to be comprehensive, i.e. must not contain all taxa
#   -R      Specify the file name of a binary model parameter file that has previously been generated
#           with RAxML using the -f e tree evaluation option. The file name should be:
#           RAxML_binaryModelParameters.runID
#   -s      Specify the name of the alignment data file in PHYLIP format
#   -S      Specify the name of a secondary structure file. The file can contain "." for
#           alignment columns that do not form part of a stem and characters "()<>[]{}" to define
#           stem regions and pseudoknots
#   -t      Specify a user starting tree file name in Newick format
#   -T      PTHREADS VERSION ONLY! Specify the number of threads you want to run.
#           Make sure to set "-T" to at most the number of CPUs you have on your machine,
#           otherwise, there will be a huge performance decrease!
#   -U      Try to save memory by using SEV-based implementation for gap columns on large gappy alignments
#           WARNING: this will only work for DNA under GTRGAMMA and is still in an experimental state.
#   -v      Display version information
#   -w      FULL (!) path to the directory into which RAxML shall write its output files
#           DEFAULT: current directory
#   -W      Sliding window size for leave-one-out site-specific placement bias algorithm
#           only effective when used in combination with "-f S"
#           DEFAULT: 100 sites
#   -x      Specify an integer number (random seed) and turn on rapid bootstrapping
#           CAUTION: unlike in version 7.0.4 RAxML will conduct rapid BS replicates under
#           the model of rate heterogeneity you specified via "-m" and not by default under CAT
#   -X      EXPERIMENTAL OPTION: This option will do a per-site estimate of protein substitution models
#           by looping over all given, fixed models LG, WAG, JTT, etc and using their respective base frequencies to independently
#           assign a prot subst. model to each site via ML optimization
#           At present this option only works with the GTR+GAMMA model, unpartitioned datasets, and in the sequential
#           version only.
#           DEFAULT: OFF
#   -y      If you want to only compute a parsimony starting tree with RAxML specify "-y",
#           the program will exit after computation of the starting tree
#           DEFAULT: OFF
#   -z      Specify the file name of a file containing multiple trees e.g. from a bootstrap
#           that shall be used to draw bipartition values onto a tree provided with "-t",
#           It can also be used to compute per site log likelihoods in combination with "-f g"
#           and to read a bunch of trees for a couple of other options ("-f h", "-f m", "-f n").
#   -#|-N   Specify the number of alternative runs on distinct starting trees
#           In combination with the "-b" option, this will invoke a multiple boostrap analysis
#           Note that "-N" has been added as an alternative since "-#" sometimes caused problems
#           with certain MPI job submission systems, since "-#" is often used to start comments.
#           If you want to use the bootstopping criteria specify "-# autoMR" or "-# autoMRE" or "-# autoMRE_IGN"
#           for the majority-rule tree based criteria (see -I option) or "-# autoFC" for the frequency-based criterion.
#           Bootstopping will only work in combination with "-x" or "-b"
#           DEFAULT: 1 single analysis

class RaxmlRunner(object):

    def __init__(self,
            working_dir_path=None,
            replace=None,
            postclean=None,
            name=None,
            verbosity=1,
            raxml_path="raxmlHPC"):
        self.dirs_to_clean = []
        self.files_to_clean = []
        if working_dir_path is None:
            self.working_dir_path = None
            self.replace = False
            self.postclean = postclean if postclean is not None else True
        else:
            self.working_dir_path = working_dir_path
            self.replace = replace if replace is not None else True
            self.postclean = postclean if postclean is not None else False
        self._name = name
        self.verbosity = verbosity
        self.messenger = get_messenger(self.verbosity)
        self.input_format = "nexus"
        self.output_format = "nexus"
        self.raxml_path = raxml_path
        self.taxon_label_map = {}

    @property
    def name(self):
        if self._name is None:
            self._name = "dendropy_raxml"
        return self._name

    def _compose_fname(self, s):
        return "RAxML_%s.%s" % (s, self.name)

    @property
    def input_seq_fname(self):
        return "%s.seqs" % self.name

    @property
    def best_tree_fname(self):
        return self._compose_fname("bestTree")

    @property
    def bipartitions_fname(self):
        return self._compose_fname("bipartitions")

    @property
    def info_fname(self):
        return self._compose_fname("info")

    @property
    def log_fname(self):
        return self._compose_fname("log")

    @property
    def parsimony_tree_fname(self):
        return self._compose_fname("parsimonyTree")

    @property
    def result_fname(self):
        return self._compose_fname("result")

    def _raxml_output_filenames(self):
        return [self.input_seq_fname,
                self.best_tree_fname,
                self.bipartitions_fname,
                self.info_fname,
                self.log_fname,
                self.parsimony_tree_fname,
                self.result_fname]

    def _raxml_output_filepaths(self):
        return [os.path.join(self.working_dir_path, fp) for fp in self._raxml_output_filenames()]

    def _get_trees(self, tree_filepath, tree_list=None, **kwargs):
        if tree_list is None:
            tree_list = dendropy.TreeList()
        tree_list.read_from_path(tree_filepath,
                self.input_format,
                **kwargs)
        return tree_list

    def _expand_path(self, path):
        return os.path.expanduser(os.path.expandvars(path))

    def _check_overwrite(self, path):
        if os.path.exists(path) and not self.replace:
            ok = input("Overwrite existing file '{}'? (y/n/all [n])? ".format(path))
            if not ok:
                return False
            ok = ok[0].lower()
            if ok == "a":
                self.replace = True
                return True
            if ok == "y":
                return True
            return False
        else:
            return True

    # def _send_info(self, msg):
    #     self.messenger.send_info(msg, wrap=False)

    # def _send_warning(self, msg):
    #     self.messenger.send_warning(msg, wrap=False)

    # def _send_error(self, msg):
    #     self.messenger.send_info(msg, wrap=False)

    def _write_dummy_seqs(self, taxon_namespace, out):
        nchar = 19
        out.write("{} {}\n".format(len(taxon_namespace), nchar))
        bases = ["A", "C", "G", "T"]
        for idx, taxon in enumerate(taxon_namespace):
            base_seq = [random.choice(bases) for x in range(nchar)]
            out.write("{}    {}\n".format(taxon.label, "".join(base_seq)))

    def _remap_taxon_labels(self, taxa):
        for idx, taxon in enumerate(taxa):
            label = "T{}".format(idx)
            self.taxon_label_map[label] = taxon.label
            taxon.label = label

    def _create_working_dir(self):
        if self.working_dir_path is None:
            self.working_dir_path = tempfile.mkdtemp()
            self.dirs_to_clean.append(self.working_dir_path)
        if not os.path.exists(self.working_dir_path):
            # self._send_info("Creating work directory: {}".format(self.working_dir_path))
            os.makedirs(self.working_dir_path)

    def _clean_working_files(self):
        for fpath in self._raxml_output_filepaths():
            if not self._check_overwrite(fpath):
                sys.exit(0)
            if os.path.exists(fpath):
                os.remove(fpath)

    def _preclean_working_dir(self):
        self._clean_working_files()

    def _postclean_working_dir(self):
        if self.postclean:
            # self._send_info("Cleaning up run files")
            for fpath in self.files_to_clean:
                # self._send_info("Deleting file: {}".format(fpath))
                try:
                    os.remove(fpath)
                except OSError as e:
                    pass
            for dir_path in self.dirs_to_clean:
                # self._send_info("Deleting directory: {}".format(dir_path))
                try:
                    os.rmdir(dir_path)
                except OSError as e:
                    pass

    def estimate_tree(self,
            char_matrix,
            raxml_args=None):

        # set up taxa
        taxa = char_matrix.taxon_namespace

        # create working directory
        self._create_working_dir()

        # remap taxon labels
        self.taxon_label_map = {}
        self._remap_taxon_labels(taxa)

        # clean working directory of previous runs
        self._preclean_working_dir()

        # write input sequences
        raxml_seqs_filepath = os.path.join(self.working_dir_path, self.input_seq_fname)
        # self._send_info("Creating RAxML dummy sequences file: {}".format(raxml_seqs_filepath))
        # if not self._check_overwrite(raxml_seqs_filepath):
        #     sys.exit(0)
        raxml_seqs_filepath_out = open(raxml_seqs_filepath, "w")
        char_matrix.write_to_stream(raxml_seqs_filepath_out, "phylip")
        raxml_seqs_filepath_out.flush()
        raxml_seqs_filepath_out.close()
        self.files_to_clean.append(raxml_seqs_filepath)
        self.files_to_clean.append(raxml_seqs_filepath + ".reduced")

        # run RAxML
        if raxml_args is None:
            raxml_args = []
        cmd = [self.raxml_path,
                '-m',
                'GTRCAT',
                '-s', raxml_seqs_filepath,
                '-n', self.name,
                '-p', str(random.randint(0, sys.maxsize))] + raxml_args
        # self._send_info("Executing: {}".format(" ".join(cmd)))
        if self.verbosity >= 2:
            stdout_pipe = None
            stderr_pipe = None
        else:
            stdout_pipe = subprocess.PIPE
            stderr_pipe = subprocess.PIPE
        p = subprocess.Popen(cmd,
            stdout=stdout_pipe,
            stderr=stderr_pipe,
            cwd=self.working_dir_path)
        stdout, stderr = processio.communicate(p)
        if p.returncode != 0:
            sys.stderr.write("[RAxML run failed]:\n\n%s\n\n" % (" ".join(cmd)))
            sys.stdout.write(stdout)
            sys.stderr.write(stderr)
            sys.exit(p.returncode)

        # # read result
        raxml_best_tree_fpath = os.path.join(self.working_dir_path, self.best_tree_fname)
        if not os.path.exists(raxml_best_tree_fpath):
            self._send_error("RAxML result not found: {}".format(raxml_best_tree_fpath))
            sys.exit(1)
        best_tree = dendropy.Tree.get_from_path(raxml_best_tree_fpath,
                "newick",
                taxon_namespace=taxa)

        # remap labels
        for taxon in best_tree.taxon_namespace:
            taxon.label = self.taxon_label_map[taxon.label]

        # # write results
        # mapped_tree.write_to_stream(self.output_dest, self.output_format)

        # clean-up
        self._postclean_working_dir()

        # # return result
        return best_tree

    def map_bipartitions(self, target_tree_fpath, bootstrap_trees_fpaths):

        # set up taxa
        taxa = dendropy.TaxonNamespace()
        taxon_label_map = {}

        # read target tree
        target_tree_fpath = self._expand_path(target_tree_fpath)
        # self._send_info("Reading target tree file: {}".format(target_tree_fpath))
        target_tree = self._get_trees(target_tree_fpath, taxon_namespace=taxa)[0]

        # read boostrap trees
        boot_trees = dendropy.TreeList()
        for fpath in bootstrap_trees_fpaths:
            fpath = self._expand_path(fpath)
            # self._send_info("Reading bootstrap tree file: {}".format(fpath))
            self._get_trees(tree_filepath=fpath, tree_list=boot_trees, taxon_namespace=taxa)
        # self._send_info("Read: {} taxa, {} bootstrap trees".format(len(taxa), len(boot_trees)))

        # create working directory
        self._create_working_dir()

        # remap taxon labels
        self.taxon_label_map = {}
        self._remap_taxon_labels(taxa)

        # write input target tree
        raxml_target_tree_filepath = os.path.join(self.working_dir_path, "{}.target_tree".format(self.name))
        # self._send_info("Creating RAxML target tree file: {}".format(raxml_target_tree_filepath))
        if not self._check_overwrite(raxml_target_tree_filepath):
            sys.exit(0)
        target_tree.write_to_path(raxml_target_tree_filepath, "newick")
        self.files_to_clean.append(raxml_target_tree_filepath)

        # write input bootstrap trees
        raxml_bootstrap_trees_filepath = os.path.join(self.working_dir_path, "{}.boot_trees".format(self.name))
        # self._send_info("Creating RAxML bootstrap tree file: {}".format(raxml_bootstrap_trees_filepath))
        if not self._check_overwrite(raxml_bootstrap_trees_filepath):
            sys.exit(0)
        boot_trees.write_to_path(raxml_bootstrap_trees_filepath, "newick")
        self.files_to_clean.append(raxml_bootstrap_trees_filepath)

        # write input (dummy) sequences
        raxml_seqs_filepath = os.path.join(self.working_dir_path, "{}.seqs".format(self.name))
        # self._send_info("Creating RAxML dummy sequences file: {}".format(raxml_seqs_filepath))
        if not self._check_overwrite(raxml_seqs_filepath):
            sys.exit(0)
        raxml_seqs_filepath_out = open(raxml_seqs_filepath, "w")
        self._write_dummy_seqs(taxa, raxml_seqs_filepath_out)
        raxml_seqs_filepath_out.flush()
        raxml_seqs_filepath_out.close()
        self.files_to_clean.append(raxml_seqs_filepath)

        # clean working directory of previous runs
        self._preclean_working_dir()

        # run RAxML
        cmd = [self.raxml_path, '-f', 'b',
                '-t', os.path.basename(raxml_target_tree_filepath),
                '-z', os.path.basename(raxml_bootstrap_trees_filepath),
                '-s', os.path.basename(raxml_seqs_filepath),
                '-m', 'GTRCAT',
                '-n', self.name]
        # self._send_info("Executing: {}".format(" ".join(cmd)))
        if self.verbosity >= 2:
            stdout_pipe = None
            stderr_pipe = None
        else:
            stdout_pipe = subprocess.PIPE
            stderr_pipe = subprocess.PIPE
        p = subprocess.Popen(cmd,
            stdout=stdout_pipe,
            stderr=stderr_pipe,
            cwd=self.working_dir_path)
        stdout, stderr = processio.communicate(p)
        if p.returncode != 0:
            self._send_error("RAxML run failed")
            if self.verbosity < 2:
                sys.stdout.write(stdout)
                sys.stderr.write(stderr)
            sys.exit(p.returncode)

        # read result
        raxml_mapped_tree_fpath = os.path.join(self.working_dir_path, self.bipartitions_fname)
        if not os.path.exists(raxml_mapped_tree_fpath):
            self._send_error("RAxML result not found: {}".format(raxml_mapped_tree_fpath))
            sys.exit(1)
        mapped_tree = dendropy.Tree.get_from_path(raxml_mapped_tree_fpath, "newick")

        # remap labels
        for taxon in mapped_tree.taxon_namespace:
            taxon.label = taxon_label_map[taxon.label]

#         # write results
#         mapped_tree.write_to_stream(self.output_dest, self.output_format)

        # clean-up
        self.files_to_clean.append(raxml_mapped_tree_fpath)
        self.files_to_clean.append(self.info_fname)
        self._postclean_working_dir()

        # return result
        return mapped_tree

    # def estimate_tree(self,
    #         char_matrix,
    #         raxml_args=None):
    #     command = [self.raxml_path] + self.raxml_args + raxml_args
    #     for dup_cmd in ['-n', '-s']:
    #         if dup_cmd in command:
    #             raise Exception("Cannot specify '%s' option when running through wrapper"  % dup_cmd)

    #     # if '-w' not in command:
    #     #     work_dir = tempfile.mkdtemp()
    #     #     clean_up_dir = True
    #     #     command.a
    #     #     command.extend(['-w', work_dir])
    #     # else:
    #     #     arg_idx = command(command.index('-w')+1)
    #     #     work_dir = os.path.abspath(arg_idx)
    #     #     command[arg_idx] = work_dir
    #     #     clean_up_dir = False
    #     work_dir = os.path.abspath('.')

    #     input_data_path = os.path.join(work_dir, "input.phy")
    #     char_matrix.write_to_path(input_data_path, "phylip")
    #     command.extend(['-n', 'dendropy.raxml'])
    #     command.extend(['-s', input_data_path])
    #     run = subprocess.Popen(command,
    #             stdin=subprocess.PIPE,
    #             stdout=subprocess.PIPE,
    #             stderr=subprocess.PIPE)
    #     stdout, stderr = processio.communicate(run)
    #     results = stdout.split("\n")
    #     if run.returncode:
    #         sys.stderr.write("\n*** ERROR FROM RAxML:\n")
    #         sys.stderr.write(stdout)
    #         sys.stderr.write("\n\n*** COMMAND SENT TO RAxML:\n    ")
    #         sys.stderr.write(' '.join(command))
    #         sys.stderr.write("\n")
    #         sys.exit(1)
    #     return results

