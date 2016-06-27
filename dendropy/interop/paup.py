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
Wrapper around calls to PAUP*, mainly for testing purposes rather than analysis.
"""

import os
import sys
import subprocess
import tempfile
import re
import csv

import dendropy
from dendropy.utility import textprocessing
from dendropy.utility import error
from dendropy.utility import metavar
from dendropy.utility import container
from dendropy.utility import messaging
from dendropy.utility import filesys
from dendropy.utility import processio
from dendropy.dataio import nexuswriter
_LOG = messaging.get_logger(__name__)

import dendropy

PAUP_PATH = os.environ.get(metavar.DENDROPY_PAUP_PATH_ENVAR, "NONE")
if PAUP_PATH == "NONE":
    DENDROPY_PAUP_INTEROPERABILITY = False
else:
    DENDROPY_PAUP_INTEROPERABILITY = True

STANDARD_PREAMBLE = "set warnreset=no increase=auto warnroot=no warnReset=no warnTree=no warnTSave=no warnBlkName=no errorStop=no errorBeep=no queryBeep=no"

class PaupService(object):

    @staticmethod
    def call(
            paup_commands,
            suppress_standard_preamble=False,
            ignore_error_returncode=False,
            ignore_nonempty_stderr=False,
            strip_extraneous_prompts_from_stdout=True,
            strip_extraneous_prompts_from_stderr=True,
            cwd=None,
            env=None,
            paup_path=PAUP_PATH,
            timeout=None,
            ):
        """
        Executes a sequence of commands in PAUP* and returns the results.

        Parameters
        ----------
        paup_commands : iterable of strings
            A list or some other iterable of strings representing PAUP
            commands.
        suppress_standard_preamble : bool
            If |True|, then the command sequence will not be prefaced by the
            standard preamble.
        ignore_error_returncode : bool
            If |True|, then a non-0 return code from the PAUP process will not
            result in an exception being raised.
        ignore_nonempty_stderr : bool
            If |True|, then the PAUP process writing to standard error will not
            result in an exception being raised.
        strip_extraneous_prompts_from_stdout : bool
            If |True|, then all occurrences of 'paup>' will be removed from the
            standard output contents.
        strip_extraneous_prompts_from_stderr : bool
            If |True|, then all occurrences of 'paup>' will be removed from the
            standard error contents.
        cwd : string
            Set the working directory of the PAUP* process to this directory.
        env : dictionary
            Environmental variables to set for the PAUP* process.
        paup_path : string
            Path to the PAUP* executable.

        Returns
        -------
        returncode : exit value of PAUP process.
        stdout : string
            Contents of the PAUP process standard output.
        stderr : string
            Contents of the PAUP process standard error.
        """
        if textprocessing.is_str_type(paup_commands):
            commands = [paup_commands]
        else:
            commands = list(paup_commands)
        if not suppress_standard_preamble:
            commands.insert(0, STANDARD_PREAMBLE)
        commands.append("quit")
        paup_block = ";\n".join(commands) + ";\n"
        invocation_command = [paup_path, "-n", "-u"]
        p = subprocess.Popen(
                invocation_command,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=cwd,
                env=env,
                )
        raw_stdout, raw_stderr = processio.communicate(p, paup_block, timeout=timeout)
        # try:
        #     raw_stdout, raw_stderr = processio.communicate(p, paup_block, timeout=timeout)
        # except TypeError as e:
        #     raise
        #     if str(e) == "communicate() got an unexpected keyword argument 'timeout'":
        #         raw_stdout, raw_stderr = processio.communicate(p, paup_block)
        #     else:
        #         raise
        stdout = raw_stdout
        stderr = raw_stderr
        if strip_extraneous_prompts_from_stdout:
            # weird dev/paup error ... lots or prompts spring up
            stdout = stdout.replace("paup>", "")
        if strip_extraneous_prompts_from_stderr:
            # weird dev/paup error ... lots or prompts spring up
            stderr = stderr.replace("paup>", "")
            chk_stderr = stderr
        else:
            chk_stderr = stderr.replace("paup>", "")
        if (p.returncode != 0 and not ignore_error_returncode) or (chk_stderr != "" and not ignore_nonempty_stderr):
            raise error.ExternalServiceError(
                    service_name="PAUP*",
                    invocation_command=invocation_command,
                    service_input=paup_block,
                    returncode = p.returncode,
                    stdout=raw_stdout,
                    stderr=raw_stderr)
        return p.returncode, stdout, stderr

    @staticmethod
    def bipartition_groups_to_split_bitmask(group_string, normalized=None):
        """
        This converts a PAUP* group representation (i.e. a string of askterisks
        and periods, where the asterisks denote the taxon index counting from
        left to right) to a mask representation:
            - a clade mask, where 1's represent descendents of the split/edge
                (with taxon index counting from right to left, i.e., first taxon
                is right-most bit)
            - a split mask, an unrooted normalized version of the above, where
                if the right most bit is not 1 the clade mask is complemented
                (and not changed otherwise).
        """
        group_string = group_string[::-1] # flip to get correct orientation
        split_bitmask = int(group_string.replace("*", "1").replace(".", "0"), 2)
        if normalized:
            mask=((2 ** len(group_string)) -1)
            return container.NormalizedBitmaskDict.normalize(split_bitmask, mask, 1)
        else:
            return split_bitmask

    def __init__(self,
            suppress_standard_preamble=False,
            ignore_error_returncode=False,
            strip_extraneous_prompts_from_stderr=True,
            strip_extraneous_prompts_from_stdout=True,
            cwd=None,
            env=None,
            paup_path=PAUP_PATH):
        self.suppress_standard_preamble = suppress_standard_preamble
        self.ignore_error_returncode = ignore_error_returncode
        self.strip_extraneous_prompts_from_stderr = strip_extraneous_prompts_from_stderr
        self.strip_extraneous_prompts_from_stdout = strip_extraneous_prompts_from_stdout
        self.cwd = cwd
        self.env = env
        self.paup_path = paup_path
        self._nexus_writer = nexuswriter.NexusWriter()
        self.commands = []

    def count_splits_from_files(self,
            tree_filepaths=None,
            is_rooted=None,
            use_tree_weights=None,
            burnin=None,
            taxa_definition_filepath=None,
            taxon_namespace=None):
        """
        Counts splits (bipartitions) in trees from files and returns the results.

        Parameters
        ----------
        tree_filepaths : iterable of strings
            A list or some other iterable of file paths containing trees in
            NEXUS format.
        is_rooted : bool
            If |True| then trees will be treated as rooted. If |False|, then
            rooting follows that specified in the tree statements, defaulting
            to unrooted if not specified.
        use_tree_weights : bool
            If |False| then tree weighting statements are disregarded.
            Otherwise, they will be regarded.
        burnin : integer
            Skip these many trees (from beginning of each source).
        taxa_definition_filepath : str
            Path of file containing TAXA block to execute. This is crucial to
            getting the taxon order (and hence, indexes, and hence, split
            bitmasks) correct. If not given, will use the first file
            given in ``tree_filepaths``.
        taxon_namespace : |TaxonNamespace|
            The |TaxonNamespace| object to populate.

        Returns
        -------
        d : dictionary
            A dictionary with the following keys and values:

                -   "bipartition_counts" : dictionary with split bitmasks as keys
                    and (weighted) counts of occurrences as values
                -   "bipartition_frequencies" : dictionary with split bitmasks as keys
                    and (weighted) proportional frequencies of occurrences as values
                -   "num_trees" : number of trees counted
                -   "taxon_namespace" : |TaxonNamespace| instance
                    corresponding to the taxa <=> split bitmask mapping
                -   "is_rooted" : indicates whether the trees were rooted or not
        """
        self.commands = []
        if taxa_definition_filepath is not None:
            self.stage_execute_file(
                    taxa_definition_filepath,
                    clear_trees=True)
        self.stage_load_trees(
            tree_filepaths=tree_filepaths,
            is_rooted=is_rooted,
            use_tree_weights=use_tree_weights,
            burnin=burnin,
            mode=7)
        self.stage_list_taxa()
        self.stage_tree_info()
        self.stage_count_splits(use_tree_weights=use_tree_weights)
        # print("\n".join(self.commands))
        returncode, stdout, stderr = self._execute_command_sequence()
        # print("\n".join(stdout))
        taxon_namespace = self.parse_taxon_namespace(stdout,
                taxon_namespace=taxon_namespace)
        is_rooted = self.parse_is_tree_rooted(stdout)
        tree_count, bipartition_counts, bipartition_freqs = self.parse_group_freqs(stdout, is_rooted=is_rooted)
        d = {
            "num_trees" : tree_count,
            "bipartition_counts" : bipartition_counts,
            "bipartition_freqs" : bipartition_freqs,
            "taxon_namespace" : taxon_namespace,
            "is_rooted" : is_rooted,
            }
        return d

    def get_split_distribution_from_files(self,
            tree_filepaths=None,
            is_rooted=None,
            use_tree_weights=None,
            burnin=None,
            taxa_definition_filepath=None,
            taxon_namespace=None,
            split_distribution=None):
        """
        Returns a SplitDistribution object based on splits given in
        tree files.

        tree_filepaths : iterable of strings
            A list or some other iterable of file paths containing trees in
            NEXUS format.
        is_rooted : bool
            If |True| then trees will be treated as rooted. If |False|, then
            rooting follows that specified in the tree statements, defaulting
            to unrooted if not specified.
        use_tree_weights : bool
            If |False| then tree weighting statements are disregarded.
            Otherwise, they will be regarded.
        burnin : integer
            Skip these many trees (from beginning of each source).
        taxa_definition_filepath : str
            Path of file containing TAXA block to execute. This is crucial to
            getting the taxon order (and hence, indexes, and hence, split
            bitmasks) correct. If not given, will use the first file
            given in ``tree_filepaths``.
        taxon_namespace : |TaxonNamespace|
            |TaxonNamespace| object to use.
        split_distribution : `SplitDistribution`
            `SplitDistribution object to use.
        """
        if split_distribution is None:
            split_distribution = dendropy.SplitDistribution(taxon_namespace=taxon_namespace)
            taxon_namespace = split_distribution.taxon_namespace
        else:
            if taxon_namespace is None:
                taxon_namespace = split_distribution.taxon_namespace
            else:
                assert split_distribution.taxon_namespace is taxon_namespace
        result = self.count_splits_from_files(
            tree_filepaths=tree_filepaths,
            is_rooted=is_rooted,
            use_tree_weights=use_tree_weights,
            burnin=burnin,
            taxa_definition_filepath=taxa_definition_filepath,
            taxon_namespace=taxon_namespace)
        for split in result["bipartition_counts"]:
            if not is_rooted:
                sd_split_key = split_distribution.normalize_bitmask(split)
            else:
                sd_split_key = split
            split_distribution.add_split_count(sd_split_key, result["bipartition_counts"][split])
        split_distribution.total_trees_counted = result["num_trees"]
        return split_distribution

    def stage_execute_file(self,
            filepath,
            clear_trees=False):
        """Executes file, optionally clearing trees from file if requested"""
        self.commands.append("execute {}".format(filepath))
        if clear_trees:
            self.commands.append("cleartrees")
        return commands

    def stage_load_trees(self,
            tree_filepaths,
            is_rooted=None,
            use_tree_weights=None,
            burnin=None,
            mode=7): # keep trees in memory, specify 3 to clear
        """
        Composes commands to load a set of trees into PAUP*, with the specified
        number of burnin dropped.
        """
        if textprocessing.is_str_type(tree_filepaths):
            raise Exception("expecting list of filepaths, not string")
        if is_rooted is None:
            rooting = ""
        elif is_rooted:
            rooting = "rooted=yes"
        else:
            rooting = "unrooted=yes"
        if use_tree_weights is None:
            treewts = ""
        elif use_tree_weights:
            treewts = "storetreewts=yes"
        else:
            treewts = "storetreewts=no"
        if burnin is None:
            burnin = 0
        gettree_template = "gett file= '{{tree_filepath}}' storebrlens=yes warntree=no {rooting} {treewts} from={burnin} mode={mode};".format(
                rooting=rooting,
                treewts=treewts,
                burnin=burnin+1,
                mode=mode)
        for tree_filepath in tree_filepaths:
            # self.commands.append(gettree_template.format(tree_filepath=tree_filepath))
            # using relpath because of a bug in PAUP* 4.0b10 with long paths passed to gettrees
            self.commands.append(gettree_template.format(tree_filepath=os.path.relpath(tree_filepath)))
        return self.commands

    def stage_list_taxa(self):
        """
        Given a data file in memory, this gets PAUP* to print a list of
        taxa that can be used to build a TaxaBlock later.
        """
        # self.commands.append("[!TAXON LIST BEGIN]\ntstatus / full;\n[!TAXON LIST END]\n")
        self.commands.append("[!TAXON LIST BEGIN]\ntstatus / full;\n[!TAXON LIST END]\n")
        return self.commands

    def stage_tree_info(self):
        self.commands.append("[!TREE INFO BEGIN]treeinfo;\n[!TREE INFO END]\n")
        return self.commands

    def stage_count_splits(self,
            use_tree_weights=None,
            majrule_filepath=None,
            majrule_freq=0.5):
        """
        Given trees in memory, this composes a command to count the split
        frequencies across the trees as well as a save the majority-rule
        consensus tree if a path is given.
        """
        percent = int(100 * majrule_freq)
        if majrule_filepath is None:
            treefile = ""
        else:
            treefile = " treefile={filepath} replace=yes "
        if use_tree_weights is None:
            treewts = ""
        elif use_tree_weights:
            treewts = "usetreewts=yes"
        else:
            treewts = "usetreewts=no"
        commands = []
        commands.append("[!SPLITS COUNT BEGIN]")
        commands.append("contree / strict=no {treefile} showtree=no grpfreq=yes majrule=yes percent={percent} {treewts}".format(
            treefile=treefile,
            percent=percent,
            treewts=treewts))
        commands.append("[!SPLITS COUNT END]")
        self.commands.extend(commands)
        return self.commands

    def stage_execute_file(self, filepath, clear_trees=False):
        """Executes file, optionally clearing trees from file if requested"""
        self.commands.append("execute '{}'".format(filepath))
        if clear_trees:
            self.commands.append("cleartrees")
        return self.commands

    ##############################################################################
    ## Processing of Output

    def parse_taxon_namespace(self, paup_output, taxon_namespace=None):
        """
        Given PAUP* output that includes a taxon listing as produced by
        ``stage_list_taxa``, this parses out and returns a taxon block.
        """
        taxlabels = []
        taxinfo_pattern = re.compile('\s*(\d+) (.*)\s+\-')
        idx = 0
        for line in paup_output:
            idx += 1
            if line == "TAXON LIST BEGIN":
                break
        for line in paup_output[idx:]:
            if line == "TAXON LIST END":
                break
            ti_match = taxinfo_pattern.match(line)
            if ti_match:
                label = ti_match.group(2).strip()
                taxlabels.append(label)
        if taxon_namespace is None:
            taxon_namespace = dendropy.TaxonNamespace()
        for taxlabel in taxlabels:
            taxon_namespace.require_taxon(label=taxlabel)
        return taxon_namespace

    def parse_is_tree_rooted(self, paup_output):
        """
        Given PAUP* output that includes a information produced by
        ``stage_tree_info``, this parses out and returns the rooting
        state of trees in memory
        """
        pattern = re.compile(r'\d+ (\w+) trees in memory')
        for line in paup_output:
            if line == "TREE INFO END":
                break
            match = pattern.match(line)
            if match:
                s = match.groups(1)[0]
                if s == "unrooted":
                    return False
                elif s == "rooted":
                    return True
                else:
                    return None
        raise Exception("Unable to find tree information")

    def parse_group_freqs(self, paup_output, is_rooted=None):
        """
        Given PAUP* output that includes a split counting procedure, this
        collects the splits and returns a dictionary of split bitmasks and their
        frequencies.
        """
        bipartitions = []
        bipartition_freqs = {}
        bipartition_counts = {}
        tree_count = None
        tree_count_pattern = re.compile('.*Majority-rule consensus of ([\d]*) tree.*', re.I)

        bipartition_section = re.compile('Bipartitions found in one or more trees and frequency of occurrence:')
        bp_full_row_with_perc_col = re.compile('([\.|\*]+)\s+([\d\.]+)\s+([\d\.]*)%')
        bp_full_row_with_no_perc_col = re.compile('([\.|\*]+)\s+([\d\.]+)')
        bp_row = re.compile('([\.|\*]+).*')

        # find tree count
        for idx, line in enumerate(paup_output):
            tp_match = tree_count_pattern.match(line)
            if tp_match:
                break
        if not tp_match:
            raise Exception("Failed to find tree count in PAUP* output")
        tree_count = int(tp_match.group(1))

        while not bp_row.match(paup_output[idx]):
            idx += 1

        split_idx = 0
        split_reps = {}
        for line in paup_output[idx:]:
            if line == "SPLITS COUNT END":
                    break
            bp_match = bp_full_row_with_perc_col.match(line)
            if not bp_match:
                bp_match = bp_full_row_with_no_perc_col.match(line)
            if bp_match:
                # full row, or end of partial rows
                if len(split_reps) == 0:
                    split_rep = bp_match.group(1)
                else:
                    split_rep = split_reps[split_idx] + bp_match.group(1)
                split_bitmask = PaupService.bipartition_groups_to_split_bitmask(split_rep, normalized=not is_rooted)
                bipartition_counts[split_bitmask] = float(bp_match.group(2))
                try:
                    bipartition_freqs[split_bitmask] = float(bp_match.group(3)) / 100
                except IndexError:
                    bipartition_freqs[split_bitmask] = bipartition_counts[split_bitmask] / 100
                split_idx += 1
            else:
                # either (1) partial row or (2) break between sections
                bp_match = bp_row.match(line)
                if not bp_match:
                    split_idx = 0
                else:
                    if split_idx in split_reps:
                        split_reps[split_idx] += bp_match.group(1)
                    else:
                        split_reps[split_idx] = bp_match.group(1)
                    split_idx += 1
        return tree_count, bipartition_counts, bipartition_freqs

    ##############################################################################
    ## Support

    def _execute_command_sequence(self):
        returncode, stdout, stderr = PaupService.call(self.commands)
        self.commands = []
        stdout = stdout.split("\n")
        stderr = stderr.split("\n")
        return returncode, stdout, stderr

##############################################################################
## Wrappers for PAUP* Services

def call(*args, **kwargs):
    return PaupService.call(*args, **kwargs)

def symmetric_difference(tree1, tree2):
    if tree1.taxon_namespace is not tree2.taxon_namespace:
        trees = dendropy.TreeList([dendropy.Tree(tree1), dendropy.Tree(tree2)])
    else:
        trees = dendropy.TreeList([tree1, tree2], taxon_namespace=tree1.taxon_namespace)
    tf = tempfile.NamedTemporaryFile("w", delete=True)
    trees.write_to_stream(tf, schema='nexus')
    tf.flush()
    assert tree1.is_rooted == tree2.is_rooted
    sd = get_split_distribution(
            tree_filepaths=[tf.name],
            taxa_filepath=tf.name,
            is_rooted=tree1.is_rooted,
            use_tree_weights=True,
            burnin=0)
    sf = sd.split_frequencies
    conflicts = 0
    for k, v in sf.items():
        if v < 1.0:
            conflicts += 1
    return conflicts

def pscore_trees(
        trees,
        char_matrix,
        pset_option_list=None,
        pscore_option_list=None,
        paup_path=PAUP_PATH):

    if pset_option_list is not None:
        pset = "pset " + " ".join(pset_option_list)
    else:
        pset = ""

    scorefile = tempfile.NamedTemporaryFile("w+", delete=True)
    pscore_command = "pscore / scorefile={}".format(scorefile.name)
    if pscore_option_list is not None:
        pscore_command = pscore_command + " ".join(pscore_option_list)
    else:
        pscore_command = pscore_command

    post_est_commands = """\
    set crit=parsimony;
    {pset}
    {pscore_command}
    """.format(pset=pset, pscore_command=pscore_command)

    paup_block = """\
    set warnreset=no;
    exe '{data_file}';
    gettrees file= '{intree_file}' warntree=no;
    {post_est_commands};
    """

    cf = tempfile.NamedTemporaryFile("w", delete=True)
    char_matrix.write_to_stream(cf, schema='nexus')
    cf.flush()
    input_tree_file_handle = tempfile.NamedTemporaryFile("w", delete=True)
    input_tree_filepath = input_tree_file_handle.name
    trees.write_to_stream(input_tree_file_handle, schema="nexus")
    input_tree_file_handle.flush()
    paup_args = {}
    paup_args["data_file"] = cf.name
    paup_args["intree_file"] = input_tree_filepath
    paup_args["post_est_commands"] = post_est_commands
    paup_block = paup_block.format(**paup_args)
    paup_run = subprocess.Popen(['%s -n' % paup_path],
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE)
    stdout, stderr = processio.communicate(paup_run, paup_block)
    if stderr:
        sys.stderr.write("\n*** ERROR FROM PAUP ***")
        sys.stderr.write(stderr)
        sys.exit(1)
    scores_str = open(scorefile.name, "r").read()
    score_rows = [r for r in scores_str.split("\n")[1:] if r != ""]
    assert len(score_rows) == len(trees)
    scores = [int(s.split()[1]) for s in score_rows]
    assert len(scores) == len(trees)
    cf.close()
    input_tree_file_handle.close()
    scorefile.close()
    return scores

def estimate_ultrametric_tree(
        char_matrix,
        topology_tree=None,
        paup_path=PAUP_PATH):
    post_est_commands = """\
    set crit=likelihood;
    root rootmethod=midpoint;
    lset userbr=no nst = 1 basefreq = eq rates = eq clock =yes;
    lscore;
    """
    if topology_tree is None:
        ultrametric_tree = estimate_tree(char_matrix,
                tree_est_criterion="nj",
                num_states=2,
                unequal_base_freqs=False,
                gamma_rates=False,
                prop_invar=False,
                extra_post_est_commands=post_est_commands)
        return ultrametric_tree
    else:
        paup_block = """\
        set warnreset=no;
        exe '%(data_file)s';
        gettrees file= '%(intree_file)s' warntree=no;
        %(post_est_commands)s;
        savetrees file=%(outtree_file)s format=nexus root=yes brlens=yes taxablk=yes maxdecimals=20;
        """
        cf = tempfile.NamedTemporaryFile("w", delete=True)
        char_matrix.write_to_stream(cf, schema='nexus')
        cf.flush()
        input_tree_file_handle = tempfile.NamedTemporaryFile("w", delete=True)
        input_tree_filepath = input_tree_file_handle.name
        topology_tree.write_to_stream(input_tree_file_handle, schema="nexus")
        input_tree_file_handle.flush()
        # output_tree_file_handle, output_tree_filepath = tempfile.mkstemp(text=True)
        output_tree_file_handle = tempfile.NamedTemporaryFile("w+", delete=True)
        output_tree_filepath = output_tree_file_handle.name
        paup_args = {}
        paup_args["data_file"] = cf.name
        paup_args["intree_file"] = input_tree_filepath
        paup_args["post_est_commands"] = post_est_commands
        paup_args["outtree_file"] = output_tree_filepath
        paup_block = paup_block % paup_args
        paup_run = subprocess.Popen(['%s -n' % paup_path],
                                    shell=True,
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE)
        stdout, stderr = processio.communicate(paup_run, paup_block)
        t = dendropy.Tree.get_from_path(output_tree_filepath, "nexus", taxon_namespace=char_matrix.taxon_namespace)
        cf.close()
        input_tree_file_handle.close()
        output_tree_file_handle.close()
        return t

def estimate_tree(char_matrix,
                    tree_est_criterion="likelihood",
                    num_states=6,
                    unequal_base_freqs=True,
                    gamma_rates=True,
                    prop_invar=True,
                    extra_pre_est_commands=None,
                    extra_post_est_commands=None,
                    paup_path='paup',
                    char_matrix_writing_kwargs=None,
                    timeout=None,
                    ):
    """
    Given a dataset, ``char_matrix``, estimates a tree using the given criterion.
    """
    paup_args = {
        'nst': num_states,
        'basefreq' : unequal_base_freqs and 'estimate' or 'equal',
        'rates' : gamma_rates and 'gamma' or 'equal',
        'pinvar' : prop_invar and 'estimate' or '0',
    }
    cf = tempfile.NamedTemporaryFile("w", delete=True)
    if not char_matrix_writing_kwargs:
        char_matrix_writing_kwargs = {}
    char_matrix.write_to_stream(cf, schema='nexus', **char_matrix_writing_kwargs)
    cf.flush()
    paup_args['datafile'] = cf.name
    # output_tree_file_handle, output_tree_filepath = tempfile.mkstemp(text=True)
    output_tree_file_handle = tempfile.NamedTemporaryFile("w+", delete=True)
    output_tree_filepath = output_tree_file_handle.name
    paup_args['est_tree_file'] = output_tree_filepath
    if extra_pre_est_commands:
        if textprocessing.is_str_type(extra_pre_est_commands):
            extra_pre_est_commands = [extra_pre_est_commands]
        paup_args["pre_est_commands"] = ";\n".join(extra_pre_est_commands)
    else:
        paup_args["pre_est_commands"] = ""
    if extra_post_est_commands:
        if textprocessing.is_str_type(extra_post_est_commands):
            extra_post_est_commands = [extra_post_est_commands]
        paup_args["post_est_commands"] = ";\n".join(extra_post_est_commands)
    else:
        paup_args["post_est_commands"] = ""
    paup_template = """\
    set warnreset=no;
    exe %(datafile)s;
    """
    if tree_est_criterion.startswith("like"):
        paup_template += """\
    lset tratio=estimate rmatrix=estimate nst=%(nst)s basefreq=%(basefreq)s rates=%(rates)s shape=estimate pinvar=%(pinvar)s userbrlens=yes;
    """
    if tree_est_criterion not in ["nj", "upgma"] :
        paup_template += """\
        set crit=%s;
        """ % tree_est_criterion
    paup_template += """\
    %(pre_est_commands)s;
    """

    if tree_est_criterion in ["nj", "upgma"] :
        paup_template += tree_est_criterion + ";"
    else:
        paup_template += "hsearch;"

    paup_template += """\
    %(post_est_commands)s;
    savetrees file=%(est_tree_file)s format=nexus root=yes brlens=yes taxablk=yes maxdecimals=20;
    """
    # paup_run = subprocess.Popen(['%s -n' % paup_path],
    #                             shell=True,
    #                             stdin=subprocess.PIPE,
    #                             stdout=subprocess.PIPE)
    # stdout, stderr = processio.communicate(paup_run, paup_template % paup_args)
    returncode, stdout, stderr = PaupService.call(
            paup_commands=paup_template % paup_args,
            paup_path=paup_path,
            timeout=timeout,
            )
    t = dendropy.Tree.get_from_path(output_tree_filepath, "nexus", taxon_namespace=char_matrix.taxon_namespace)
    cf.close()
    output_tree_file_handle.close()
    return t

def estimate_model(char_matrix,
                    tree_model=None,
                    num_states=6,
                    unequal_base_freqs=True,
                    gamma_rates=True,
                    prop_invar=True,
                    tree_est_criterion="likelihood",
                    tree_user_brlens=True,
                    paup_path='paup'):
    """
    Given a dataset, ``char_matrix``, uses client-supplied tree or estimates a
    tree, and character substitution model for the data.
    Returns a tuple, consisting of a trees block with the tree(s) used for the
    estimated character model, and a dictionary with estimates of rates, kappa,
    base_frequencies, alpha, prop_invar, etc. as well as likelihood.
    """
    paup_args = {
        'nst': num_states,
        'basefreq' : unequal_base_freqs and 'estimate' or 'equal',
        'rates' : gamma_rates and 'gamma' or 'equal',
        'pinvar' : prop_invar and 'estimate' or '0',
    }
    if tree_model is not None:
        assert tree_model.taxon_namespace is char_matrix.taxon_namespace
        tf = tempfile.NamedTemporaryFile("w", delete=True)
        tree_model.write_to_stream(tf, 'nexus')
        tf.flush()
        paup_args['tree'] = "gettrees file=%s storebrlens=yes;" % tf.name
    else:
        if tree_est_criterion in ["nj", "upgma"] :
            paup_args['tree'] = tree_est_criterion
        else:
            paup_args['tree'] = "set crit=%s; hsearch; set crit=like;" % tree_est_criterion
    if tree_user_brlens:
        paup_args['userbrlens'] = 'yes'
    else:
        paup_args['userbrlens'] = 'no'

    cf = tempfile.NamedTemporaryFile("w", delete=True)
    char_matrix.write_to_stream(cf, schema='nexus')
    cf.flush()
    paup_args['datafile'] = cf.name
    # output_tree_file_handle, output_tree_filepath = tempfile.mkstemp(text=True)
    output_tree_file_handle = tempfile.NamedTemporaryFile("w+", delete=True)
    output_tree_filepath = output_tree_file_handle.name
    paup_args['est_tree_file'] = output_tree_filepath
    paup_template = """\
    set warnreset=no;
    exe %(datafile)s;
    set crit=like;
    lset tratio=estimate rmatrix=estimate nst=%(nst)s basefreq=%(basefreq)s rates=%(rates)s shape=estimate pinvar=%(pinvar)s userbrlens=%(userbrlens)s;
    %(tree)s;
    lscore 1 / userbrlens=%(userbrlens)s;
    savetrees file=%(est_tree_file)s format=nexus root=yes brlens=yes taxablk=yes maxdecimals=20;
"""
    paup_run = subprocess.Popen(['%s -n' % paup_path],
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE)
    stdout, stderr = processio.communicate(paup_run, paup_template % paup_args)
    results = {}
    patterns = {
        'likelihood' : re.compile('-ln L\s+([\d\.]+)'),
        'rAC' : re.compile('  AC\s+([\d\.]+)'),
        'rAG' : re.compile('  AG\s+([\d\.]+)'),
        'rAT' : re.compile('  AT\s+([\d\.]+)'),
        'rCG' : re.compile('  CG\s+([\d\.]+)'),
        'rCT' : re.compile('  CT\s+([\d\.]+)'),
        'rGT' : re.compile('  GT\s+([\d\.]+)'),
        'kappa': re.compile('  kappa\s+([\d\.]+)'),
        'prop_invar' : re.compile('P_inv\s+([\d\.]+)'),
        'alpha' : re.compile('Shape\s+([\S]+)'),
        'pA' : re.compile('  A\s+([\d\.]+)'),
        'pC' : re.compile('  C\s+([\d\.]+)'),
        'pG' : re.compile('  G\s+([\d\.]+)'),
        'pT' : re.compile('  T\s+([\d\.]+)'),
    }
    for value_name in patterns:
        results[value_name] = None
    for line in stdout.split('\n'):
        for value_name in patterns:
            m = patterns[value_name].match(line)
            if m:
                results[value_name] = m.group(1)
    for value_name in results.keys():
        if value_name == 'likelihood':
            results[value_name] = -1 * float(results[value_name])
            results["log_likelihood"] = results[value_name]
        elif results[value_name] is not None:
            try:
                results[value_name] = float(results[value_name])
            except:
                pass
    t = dendropy.Tree.get_from_path(output_tree_filepath, "nexus", taxon_namespace=char_matrix.taxon_namespace)
    cf.close()
    output_tree_file_handle.close()
    return t, results

def prune_taxa_from_trees(trees, taxa, paup_path='paup'):
    """
    Drops Taxon objects given in container ``taxa`` from TreeList ``trees``
    """
    tf = tempfile.NamedTemporaryFile("w", delete=True)
    trees.write_to_stream(tf, schema='nexus')
    tf.flush()
    output_tree_file_handle = tempfile.NamedTemporaryFile("w+", delete=True)
    output_tree_filepath = output_tree_file_handle.name
    tax_idxs = [ str(trees.taxon_namespace.index(t)+1) for t in taxa ]
    tax_idxs = " ".join(tax_idxs)
    paup_template = """\
    set warnreset=no;
    exe %s;
    gett file=%s storebrlens=yes;
    delete %s / prune;
    savetrees file=%s format=nexus brlens=user taxablk=yes maxdecimals=20;
    """ % (tf.name,
            tf.name,
            tax_idxs,
            output_tree_filepath)
    paup_run = subprocess.Popen(['%s -n' % paup_path],
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE)
    stdout, stderr = processio.communicate(paup_run, paup_template)
    t = dendropy.TreeList.get_from_path(output_tree_filepath,
            "nexus",
            taxon_namespace=trees.taxon_namespace)
    output_tree_file_handle.close()
    return t

###############################################################################
## PAUP* WRAPPERS

class PaupSession(processio.Session):
    """
    Starts a PAUP* session, which remains active until explicitly closed.
    Various commands can get executed and results returned.
    """

    EOC_FLAG = "@@@END-OF-COMMAND@@@"
    FLAG_DETECT = re.compile(r'^\s*%s\s*$' % EOC_FLAG, re.MULTILINE)
    EOC_FLAG_STRIP = re.compile(r"^(paup>)*\s*(\[!)*" + EOC_FLAG + "(\])*\s*$", re.MULTILINE)
    # FLAG_DETECT = re.compile(r'[^\[]\s*%s\s*[^\]]' % EOC_FLAG, re.MULTILINE)

    def __init__(self, paup_path=None):
        processio.Session.__init__(self, join_err_to_out=False)
        if paup_path is None:
            self.paup_path = PAUP_PATH
        else:
            self.paup_path = paup_path
        self.start([self.paup_path])

    def __del__(self):
        self.stop()

    def stop(self):
        if self.process:
            try:
                self.process.terminate()
            except:
                pass
        self.process = None

    def send_command(self, command):
        command = command + ";\n"
        command = command + "[!" + self.EOC_FLAG + "]\n"
        self.process.stdin.write(command)
        self.process.stdin.flush()
        stdout_block = ""
        while True:
            stdout = self._stdout_reader.read()
            if stdout is not None:
                stdout_block = stdout_block + stdout
            if self.FLAG_DETECT.search(stdout_block):
                stdout_block = self.EOC_FLAG_STRIP.sub("", stdout_block)
                break
            # else:
            #     print stdout_block
        stderr_block = ""
        while True:
                stderr = self._stderr_reader.read()
                if stderr is not None:
                    stderr_block += stderr
                else:
                    break
        return stdout_block, stderr_block

    def execute_file(self, filepath):
        return self.send_command("set warnreset=no; execute %s;\n" % filepath)

    def read_data(self, data):
        """
        Writes ``data`` as NEXUS-formatted file and
        executes file within processio.
        """
        cf = tempfile.NamedTemporaryFile("w", delete=True)
        data.write_to_stream(cf, "nexus")
        cf.flush()
        stdout, stderr = self.execute_file(cf.name)
        return stdout, stderr

