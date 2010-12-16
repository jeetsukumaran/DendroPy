#! /usr/bin/env python

"""
Wraps RAxML 'map bipartition' functionality to map bootstrap support onto a target tree.
"""

import argparse
import tempfile
import random
import sys
import os
import subprocess

import dendropy
from dendropy.utility.messaging import ConsoleMessenger

def get_messenger(verbosity=1):
    if verbosity == 0:
        messaging_level = ConsoleMessenger.ERROR_MESSAGING_LEVEL
    else:
        messaging_level = ConsoleMessenger.INFO_MESSAGING_LEVEL
    messenger = ConsoleMessenger(name="raxml-map-bipartitions",
                    messaging_level=messaging_level)
    return messenger

class RAxML(object):

    def __init__(self,
            working_dir_path=None,
            replace=None,
            postclean=None,
            name=None,
            verbosity=1,
            raxml_path="raxmlHPC",
            output_dest=None):
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
        self.output_dest = output_dest if output_dest is not None else sys.stdout
        self.raxml_path = raxml_path

    @property
    def name(self):
        if self._name is None:
            self._name = "raxml_map_bipartitions"
        return self._name

    @property
    def info_fname(self):
        return "RAxML_info.{}".format(self.name)

    @property
    def bipartitions_fname(self):
        return "RAxML_bipartitions.{}".format(self.name)

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
            ok = raw_input("Overwrite existing file '{}'? (y/n/all [n])? ".format(path))
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

    def _send_info(self, msg):
        self.messenger.send_info(msg, wrap=False)

    def _send_warning(self, msg):
        self.messenger.send_warning(msg, wrap=False)

    def _send_error(self, msg):
        self.messenger.send_info(msg, wrap=False)

    def _write_dummy_seqs(self, taxon_set, out):
        nchar = 19
        out.write("{} {}\n".format(len(taxon_set), nchar))
        bases = ["A", "C", "G", "T"]
        for idx, taxon in enumerate(taxon_set):
            base_seq = [random.choice(bases) for x in range(nchar)]
            out.write("{}    {}\n".format(taxon.label, "".join(base_seq)))

    def map_bipartitions(self, target_tree_fpath, bootstrap_trees_fpaths):

        # set up taxa
        taxa = dendropy.TaxonSet()
        taxon_label_map = {}

        # read target tree
        target_tree_fpath = self._expand_path(target_tree_fpath)
        self._send_info("Reading target tree file: {}".format(target_tree_fpath))
        target_tree = self._get_trees(target_tree_fpath, taxon_set=taxa)[0]

        # read boostrap trees
        boot_trees = dendropy.TreeList()
        for fpath in bootstrap_trees_fpaths:
            fpath = self._expand_path(fpath)
            self._send_info("Reading bootstrap tree file: {}".format(fpath))
            self._get_trees(tree_filepath=fpath, tree_list=boot_trees, taxon_set=taxa)
        self._send_info("Read: {} taxa, {} bootstrap trees".format(len(taxa), len(boot_trees)))

        # remap taxon labels
        for idx, taxon in enumerate(taxa):
            label = "T{}".format(idx)
            taxon_label_map[label] = taxon.label
            taxon.label = label

        # create working directory
        if self.working_dir_path is None:
            self.working_dir_path = tempfile.mkdtemp()
            self.dirs_to_clean.append(self.working_dir_path)
        if not os.path.exists(self.working_dir_path):
            self._send_info("Creating work directory: {}".format(self.working_dir_path))
            os.makedirs(self.working_dir_path)

        # write input target tree
        raxml_target_tree_filepath = os.path.join(self.working_dir_path, "{}.target_tree".format(self.name))
        self._send_info("Creating RAxML target tree file: {}".format(raxml_target_tree_filepath))
        if not self._check_overwrite(raxml_target_tree_filepath):
            sys.exit(0)
        target_tree.write_to_path(raxml_target_tree_filepath, "newick")
        self.files_to_clean.append(raxml_target_tree_filepath)

        # write input bootstrap trees
        raxml_bootstrap_trees_filepath = os.path.join(self.working_dir_path, "{}.boot_trees".format(self.name))
        self._send_info("Creating RAxML bootstrap tree file: {}".format(raxml_bootstrap_trees_filepath))
        if not self._check_overwrite(raxml_bootstrap_trees_filepath):
            sys.exit(0)
        boot_trees.write_to_path(raxml_bootstrap_trees_filepath, "newick")
        self.files_to_clean.append(raxml_bootstrap_trees_filepath)

        # write input (dummy) sequences
        raxml_seqs_filepath = os.path.join(self.working_dir_path, "{}.seqs".format(self.name))
        self._send_info("Creating RAxML dummy sequences file: {}".format(raxml_seqs_filepath))
        if not self._check_overwrite(raxml_seqs_filepath):
            sys.exit(0)
        raxml_seqs_filepath_out = open(raxml_seqs_filepath, "w")
        self._write_dummy_seqs(taxa, raxml_seqs_filepath_out)
        raxml_seqs_filepath_out.flush()
        raxml_seqs_filepath_out.close()
        self.files_to_clean.append(raxml_seqs_filepath)

        # clean working directory of previous runs
        raxml_info_fpath = os.path.join(self.working_dir_path, self.info_fname)
        raxml_mapped_tree_fpath = os.path.join(self.working_dir_path, self.bipartitions_fname)
        if not self._check_overwrite(raxml_mapped_tree_fpath) or not self._check_overwrite(raxml_info_fpath):
            sys.exit(0)
        if os.path.exists(raxml_mapped_tree_fpath):
            os.remove(raxml_mapped_tree_fpath)
        if os.path.exists(raxml_info_fpath):
            os.remove(raxml_info_fpath)

        # run RAxML
        cmd = [self.raxml_path, '-f', 'b',
                '-t', os.path.basename(raxml_target_tree_filepath),
                '-z', os.path.basename(raxml_bootstrap_trees_filepath),
                '-s', os.path.basename(raxml_seqs_filepath),
                '-m', 'GTRCAT',
                '-n', self.name]
        self._send_info("Executing: {}".format(" ".join(cmd)))
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
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            self._send_error("RAxML run failed")
            if self.verbosity < 2:
                sys.stdout.write(stdout)
                sys.stderr.write(stderr)
            sys.exit(p.returncode)

        # read result
        if not os.path.exists(raxml_mapped_tree_fpath):
            self._send_error("RAxML result not found: {}".format(raxml_mapped_tree_fpath))
            sys.exit(1)
        mapped_tree = dendropy.Tree.get_from_path(raxml_mapped_tree_fpath, "newick")

        # remap labels
        for taxon in mapped_tree.taxon_set:
            taxon.label = taxon_label_map[taxon.label]

        # write results
        mapped_tree.write_to_stream(self.output_dest, self.output_format)

        # clean-up
        if self.postclean:
            self._send_info("Cleaning up run files")
            self.files_to_clean.append(raxml_mapped_tree_fpath)
            self.files_to_clean.append(raxml_info_fpath)
            for fpath in self.files_to_clean:
                self._send_info("Deleting file: {}".format(fpath))
                try:
                    os.remove(fpath)
                except OSError as e:
                    self._send_error("Failed: {}".format(str(e)))
            for dir_path in self.dirs_to_clean:
                self._send_info("Deleting directory: {}".format(dir_path))
                try:
                    os.rmdir(dir_path)
                except OSError as e:
                    self._send_error("Failed: {}".format(str(e)))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
            "-t", "--target-tree",
            dest="target_tree_file",
            required=True,
            type=str,
            default=None,
            help="path to target tree file")
    parser.add_argument(
            dest="bootstrap_files",
            nargs='+',
            metavar="FILE",
            type=str,
            help="path(s) to bootstrap tree files")
    parser.add_argument(
            "-w", "--working-directory",
            type=str,
            default=None,
            help="working directory for temporary files")
    parser.add_argument(
            "-n", "--name",
            type=str,
            default=None,
            help="name (suffix) for temporary output files")
    parser.add_argument(
            "-r", "--replace",
            action="store_true",
            default=False,
            help="replace/overwrite existing files with the same name as the temporary output files without prompting")
    parser.add_argument(
            "--no-clean",
            action="store_true",
            default=False,
            help="preserve (do not delete) temporary files")
    parser.add_argument(
            "-v", "--verbosity",
            type=int,
            default=1,
            help="progress message noise level (default = %(default)s)")
    parser.add_argument(
            "-a", "--app-path",
            dest="raxml_path",
            type=str,
            default="raxmlHPC",
            help="path to RAxML binary (default = '%(default)s')")

    args = parser.parse_args()

    rx = RAxML(working_dir_path=args.working_directory,
            replace=args.replace,
            postclean=not args.no_clean,
            name=args.name,
            verbosity=args.verbosity,
            raxml_path=args.raxml_path)

    rx.map_bipartitions(target_tree_fpath=args.target_tree_file,
            bootstrap_trees_fpaths=args.bootstrap_files)

if __name__ == "__main__":
    main()
