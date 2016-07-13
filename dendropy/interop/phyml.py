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
Wrapper around calls to PhyML.
"""

import errno
import glob
import os
import re
import subprocess
import tempfile

import dendropy
from dendropy.utility import processio


def _get_path(seq_path, phyml_suffix):
    """
    Get the full path to a PhyML output file. Returns the path regardless
    of whether the filename has an extension (e.g. ".txt") or not.
    Return None if no file matches the constructed search pattern.
    """
    match_pattern = seq_path + phyml_suffix + '*'
    matches = glob.glob(match_pattern)
    if len(matches) == 1:
        path = matches[0]
    elif len(matches) == 0:
        path = None
    else:
        raise RuntimeError(
            'More than one file matches the pattern: ' +
            match_pattern + '\nFilenames: ' + matches)
    return path


def _read_phyml_file(seq_path, phyml_suffix, file_type="text"):
    """
    Read a PhyML output file.

    Parameters
    ----------
    seq_path : str
        Path to the sequence (input) file.
    phyml_suffix : str
        Filename suffix added by PhyML to the filename (e.g.
        "_phyml_tree"), but without potential filename extension.
    file_type : {"text", "tree", treelist"}

    Returns
    -------
    t : tuple with 2 elements
        The full path to the file and the output in a format
        determined by `file_type`.
    """
    path = _get_path(seq_path, phyml_suffix)
    if path:
        if file_type == "text":
            with open(path, "r") as f:
                output = f.read()
        elif file_type == "tree":
            output = dendropy.Tree.get_from_path(path, "newick")
        elif file_type == "treelist":
            output = dendropy.TreeList.get_from_path(path, "newick")
        else:
            raise ValueError(
                "File type should be one of text, tree or treelist")
    else:
        output = None
    t = (path, output)
    return t


def run_phyml(
        phyml_path, char_matrix, data_type=None, parsimony_starting_tree=False,
        bootstrap=None, subst_model=None, amino_acid_rates=None,
        state_freqs=None, ts_tv_ratio=None, prop_invar=None, gamma_cats=None,
        gamma_shape=None, gamma_cat_median=False, free_rates=False,
        integrated_branch_length=False, codon_position=None, search_move=None,
        starting_tree=None, optimization=None, random_starting_tree=False,
        num_random_starting_trees=None, random_seed=None,
        site_likelihoods=False, trace_search=False, run_id=None,
        alias_subpattern=False):
    """
    Wrapper for running PhyML via its command-line interface.

    A parameter value set to None will in most cases result in the PhyML
    default value. Check the output to verify that your analysis was set up
    properly. Consult the PhyML documentation for details on parameters and
    default values.

    Parameters
    ----------
    phyml_path : str
        Path to PhyML executable.
    char_matrix : |CharacterMatrix|
        Matrix with data to be analyzed.
    data_type : str
        PhyML data type: "nt" (default) for nucleotide, "aa" for
        amino-acid sequences, or "generic".
    parsimony_starting_tree : bool
        If True, a minimum parsimony starting tree is used. This
        option is taken into account when `starting_tree` is False
        and when tree topology modifications are to be done.
    bootstrap : int
        * > 0 : the number of bootstrap replicates to generate.
        *   0 : neither approximate likelihood ratio test nor bootstrap values
                will be computed.
        *  -1 : approximate likelihood ratio test returning aLRT statistics.
        *  -2 : approximate likelihood ratio test returning Chi2-based.
                parametric branch supports.
        *  -4 : SH-like branch supports alone.
        *  -5 : (default) approximate Bayes branch supports.
    subst_model : str
        Substitution model name.
        * Nucleotide-based models : "HKY85" (default), "JC69", "K80", "F81",
          "F84", "TN93", "GTR", or a custom GTR-family model, e.g. "00000".
        * Amino-acid based models : "LG" (default), "WAG" ,"JTT", "MtREV",
          "Dayhoff", "DCMut", "RtREV", "CpREV", "VT", "AB", "Blosum62",
          "MtMam", "MtArt", "HIVw", "HIVb", "custom".
    amino_acid_rates : str
        amino acid substitution rate matrix in PAML format. It is compulsory
        to use this option when analyzing amino acid sequences with the
        "custom" substitution model.
    state_freqs : str or list of floats
        * "e" : the character frequencies will be determined as follows :
            - Nucleotide sequences: (Empirical) the equilibrium base
              frequencies are estimated by counting the occurence of the
              different bases in the alignment.
            - Amino-acid sequences: (Empirical) the equilibrium amino-acid
              frequencies are estimated by counting the occurence of the
              different amino-acids in the alignment.
        * "m" : the character frequencies are determined as follows :
            - Nucleotide sequences: (ML) the equilibrium base frequencies are
              estimated using maximum likelihood.
            - Amino-acid sequences: (Model) the equilibrium amino-acid
              frequencies are estimated using the frequencies defined by
              the substitution model.
        * "fA,fC,fG,fT" : only valid for nucleotide-based models. fA, fC, fG
            and fT are floating numbers that correspond to the frequencies of
            A, C, G and T respectively (WARNING: do not use any blank space
            between your values of nucleotide frequencies, only commas!)
    ts_tv : float or str
        transition/transversion ratio. DNA sequences only. Can be a fixed
        positive value (ex: 4.0) or "e" to get the maximum likelihood
        estimate.
    prop_invar : float or str
        proportion of invariable sites. Can be a fixed value in the [0,1]
        range or "e" to get the maximum likelihood estimate.
    gamma_cats : int
        number of relative substitution rate categories. Must be a positive
        integer. Default value 4.
    gamma_shape : float or str
        distribution of the gamma distribution shape parameter. Can be a
        fixed positive value or "e" to get the maximum likelihood estimate.
    gamma_cat_median : bool
        If True, use median instead of mean as the middle of each substitution
        rate class in the discrete gamma distribution.
    free_rates : bool
        If True, the FreeRate model of substitution rate variation across
        sites will be used.
    integrated_branch_length : bool
        If True, the integrated length (IL) model will be used. The IL model
        can be considered as an approximation to the covarion model.
    codon_position : {1, 2, 3}
        When analyzing an alignment of coding sequences, use this option to
        consider only the first, second or the third coding position.
    search_move : {"NNI", "SPR", "BEST"}
        Tree topology search operation option. Can be either "NNI" (default,
        fast) or "SPR" (a bit slower than NNI) or "BEST" (best of NNI and SPR
        search).
    starting_tree : |Tree|
        User-provided starting tree.
    optimization : {"tlr", "tl", "lr", "l", "r", "n"}
        Specify which parameters to optmimize. Tree topology (t),
        branch lengths (l), rate parameters (r) and no parameter (n).
    random_starting_tree : bool
        If True, sets the initial tree to random. It is only valid if SPR
        searches are to be performed.
    num_random_starting_trees : int
        Number of initial random trees to be used. It is only valid if SPR
        searches are to be performed.
    random_seed : int
        Seed used to initiate the random number generator.
    site_likelihoods : bool
        If True, return likelood for each site.
    trace_search : bool
        If True, return each phylogeny explored during the tree search.
    run_id : str
        Append an ID-string to the PhyML output.
    alias_subpattern : bool
        If True, site aliasing is generalized at the subtree level.
        Sometimes lead to faster calculations. See Kosakovsky Pond SL,
        Muse SV, Sytematic Biology (2004) for an example.

    Returns
    -------
    result : :class:`~dendropy.interop.phyml.PhymlResult`
    """
    char_matrix_f = tempfile.NamedTemporaryFile()

    # Compose arguments
    args = []
    args.append(phyml_path)
    args.append("-i%s" % char_matrix_f.name)
    if data_type:
        args.append("-d%s" % str(data_type))
    args.append("-q")
    args.append("-n1")
    if parsimony_starting_tree:
        args.append("-p")
    if bootstrap:
        args.append("-b%s" % str(bootstrap))
    if subst_model:
        args.append("-m%s" % str(subst_model))
    if amino_acid_rates:
        args.extend(["--aa_rate_file", char_matrix_f.name + "_aa_rate"])
    if state_freqs:
        if isinstance(state_freqs, str):
            args.append("-f%s" % state_freqs)
        else:
            args.append("-f%s" % (",".join([str(s) for s in state_freqs])))
    if ts_tv_ratio:
        args.append("-t%s" % str(ts_tv_ratio))
    if prop_invar:
        args.append("-v%s" % str(prop_invar))
    if gamma_cats:
        args.append("-c%s" % str(gamma_cats))
    if gamma_shape:
        args.append("-a%s" % str(gamma_shape))
    if gamma_cat_median:
        args.append("--use_median")
    if free_rates:
        args.append("--freerates")
    if integrated_branch_length:
        args.append("--il")
    if codon_position:
        args.extend(["--codpos", str(codon_position)])
    if search_move:
        args.append("-s%s" % str(search_move))
    if starting_tree:
        args.append("-u%s" % str(char_matrix_f.name + "_starting_tree"))
    if optimization:
        args.append("-o%s" % str(optimization))
    if random_starting_tree:
        args.append("--rand_start")
    if num_random_starting_trees:
        args.extend(["--n_rand_starts", str(num_random_starting_trees)])
    if site_likelihoods:
        args.append("--print_site_lnl")
    if random_seed:
        args.extend(["--r_seed%s" % str(random_seed)])
    if trace_search:
        args.append("--print_trace")
    if run_id:
        args.extend(["--run_id", run_id])
    args.append("--quiet")
    args.append("--no_memory_check")
    if alias_subpattern:
        args.append("--alias_subpatt")
    command_line = ' '.join(args)
    try:
        # Write data to files
        char_matrix.write_to_path(
            char_matrix_f.name, "phylip", spaces_to_underscores=True)
        if starting_tree:
            starting_tree.write_to_path(
                char_matrix_f.name + "_starting_tree",
                "newick", preserve_spaces=False)
        if amino_acid_rates:
            with open(char_matrix_f.name + "_aa_rate", "w") as aa_rate_f:
                aa_rate_f.write(str(amino_acid_rates))
        # Call PhyML
        proc = subprocess.Popen(
            args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = processio.communicate(proc)
        # Check output
        if stderr or proc.returncode != 0:
            if stderr:
                raise RuntimeError(
                    "PhyML error: %s\n%s" % (command_line, stderr))
            else:
                raise RuntimeError("PhyML error: %s" % stdout)
        else:
            # Collect output
            result = PhymlResult()
            result.command_line = command_line
            result.stdout_text = stdout
            result.best_tree = dendropy.Tree.get_from_path(
                char_matrix_f.name + "_phyml_tree", "newick")
            with open(char_matrix_f.name + "_phyml_stats", "r") as stats_f:
                result.stats_text = stats_f.read()
            if os.path.isfile(char_matrix_f.name + "_phyml_boot_trees"):
                result.boot_trees = dendropy.TreeList.get_from_path(
                    char_matrix_f.name + "_phyml_boot_trees", "newick")
            if os.path.isfile(char_matrix_f.name + "_phyml_boot_stats"): 
                with open(
                    char_matrix_f.name + "_phyml_boot_stats", "r"
                ) as boot_stats_f:
                    result.boot_stats_text = boot_stats_f.read()
            if os.path.isfile(char_matrix_f.name + "_phyml_rand_trees"):  
                result.rand_trees = dendropy.TreeList.get_from_path(
                    char_matrix_f.name + "_phyml_rand_trees", "newick")
            if os.path.isfile(char_matrix_f.name + "_phyml_lk"):  
                with open(char_matrix_f.name + "_phyml_lk", "r") as lk_f:
                    result.site_likelihoods_text = lk_f.read()
            if os.path.isfile(char_matrix_f.name + "_phyml_trace"):  
                result.search_trace_trees = dendropy.TreeList.get_from_path(
                    char_matrix_f.name + "_phyml_trace", "newick")
    finally:
        # Clean up
        char_matrix_f.close()
        for phyml_file in glob.glob(char_matrix_f.name + "*"):
            os.remove(phyml_file)

    return result


class PhymlResult(object):
    """
    Container for PhyML output.

    Attributes
    ----------
    command_line : str
        The PhyML command that was run.
    stdout_text : str
        Text that was written to stdout (by PhyML).
    best_tree : |Tree|
        Maximum likelihood tree (file: "*._phyml_tree").
    stats_text : str
        Maximum likelihood model parameter estimates
        (file: "*._phyml_stats").
    boot_trees : |TreeList|
        Maximum likelihood trees from analyses of bootstrap
        replicates (file: "*_phyml_boot_trees").
    boot_stats : str
        Maximum likelihood parameter estimates from analyses
        of bootstrap replicates (file: "*_phyml_boot_stats").
    rand_trees : |Tree|
        Maximum likelihood trees from analyses with random
        starts (file: "*_phyml_rand_trees").
    site_likelihoods_text : str
        Likelihood for each site (file: "*_phyml_lk").
    search_trace_trees : |TreeList|
        Trees explored during the tree search (file: "*_phyml_trace").
    best_log_likelihood : float
        The best log-likelihood value.
    output_files : dict
        Paths to temporary files used to populate the result object.
    """

    def __init__(self):
        self.command_line = None
        self.stdout_text = None
        self.best_tree = None
        self.stats_text = None
        self.boot_trees = None
        self.boot_stats_text = None
        self.rand_trees = None
        self.site_likelihoods_text = None
        self.search_trace_trees = None
        self.output_files = None

    @property
    def best_log_likelihood(self):
        match_best = re.search(
            r"\. Log-likelihood:\s+(.+)\n", self.stats_text)
        if match_best is not None:
            return float(match_best.group(1))
