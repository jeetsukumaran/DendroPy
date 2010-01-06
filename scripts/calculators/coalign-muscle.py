#! /usr/bin/env python

############################################################################
##  coalign-muscle.py
##
##  Profile alignment of nucleotides.
##
##  Copyright 2009 Jeet Sukumaran.
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
############################################################################

import os
import sys
import subprocess
from optparse import OptionParser

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

from dendropy import datasets

def compose_default_output_filepath(infilepath, output_dir=None):
    if not output_dir:
        output_dir = os.path.dirname(infilepath)
    return os.path.join(output_dir, os.path.splitext(os.path.basename(infilepath))[0] + '.coaligned.fasta')

def compose_protein_filepath(infilepath, protein_dir):
    return os.path.join(protein_dir, os.path.basename(infilepath) + '.trans')

def compose_aligned_protein_filepath(infilepath, protein_dir):
    return os.path.join(protein_dir, os.path.splitext(os.path.basename(infilepath))[0] + '.trans.aligned')

def translate_nucleotide_file(infilepath, protfilepath, rna=False):
    infile = open(infilepath)
    ofile = file(protfilepath, 'w')
    d = datasets.Dataset()
    if rna:
        d.read(infile, "RNAFASTA")
    else:
        d.read(infile, "DNAFASTA")
    chars = d.char_blocks[0]
    for t in d.taxa_blocks[0]:
        s = chars[t]
        nucs = Seq(s.values_as_string(), generic_dna)
        prots = nucs.translate()
        ofile.write(">%s\n%s\n\n" % (t.label, prots.tostring()))

def load_fasta(fpath, schema="DNAFASTA"):
    d = datasets.Dataset()
    data = {}
    d.read(open(fpath, "rU"), schema)
    chars = d.char_blocks[0]
    for t in d.taxa_blocks[0]:
        data[t.label] = chars[t].values_as_string()
    return data

def muscle_align_codons(infilepath, outfilepath, protein_dir=None, rna=False, muscle_path="muscle", muscle_options=''):
    if not protein_dir:
        protein_dir = os.path.dirname(outfilepath)
    protfilepath = compose_protein_filepath(infilepath, protein_dir)
    translate_nucleotide_file(infilepath, protfilepath, rna)

    # align the proteins
    aligned_protfilepath = compose_aligned_protein_filepath(infilepath, protein_dir)
    command = '%s -in %s -out %s -stable %s' % (muscle_path, protfilepath, aligned_protfilepath, muscle_options)

    sys.stdout.write("%s\n" % command)
    p = subprocess.Popen(command, shell=True)
    p.wait()

    # create the nucleotide alignment, based on the protein alignment
    aligned_proteins = load_fasta(aligned_protfilepath, "PROTEINFASTA")
    if rna:
        nucleotides = load_fasta(infilepath, "RNAFASTA")
    else:
        nucleotides = load_fasta(infilepath, "DNAFASTA")

    ofile = file(outfilepath, 'w')

    for prot_title, prot_seq in aligned_proteins.items():
        nuc_index = 0
        aligned_nuc = []
        for prot_index in range(0, len(prot_seq)):
            if prot_seq[prot_index] == '-':
                aligned_nuc.append('---')
            else:
                aligned_nuc.append(nucleotides[prot_title][nuc_index])
                aligned_nuc.append(nucleotides[prot_title][nuc_index+1])
                aligned_nuc.append(nucleotides[prot_title][nuc_index+2])
                nuc_index = nuc_index + 3
        ofile.write(">%s\n%s\n\n" % (prot_title, ''.join(aligned_nuc)))

_prog_name = "COALIGN-Muscle"
_prog_version = "%prog version 1.0"

def main():

    description = "COALIGN-Muscle -- Codon Alignment using MUSCLE -- by Jeet Sukumaran"
    usage = "%prog [options] input [output]"

    parser = OptionParser(usage=usage, add_help_option=True, version = _prog_version, description=description)

    parser.add_option('-r','--rna',
                      action='store_true',
                      dest='rna',
                      default=False,
                      help="RNA nucleodtide data")

    parser.add_option('-a','--muscle',
                      dest='muscle_path',
                      default='muscle',
                      help='path to muscle application (default = "%default")')

    print 'COALIGN-Muscle -- Codon Alignment using MUSCLE'
    print 'by Jeet Sukumaran'

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print
        print 'Nothing to do: no input path specified.'
        sys.exit(1)
    elif len(args) == 1:
        input = args[0]
        output = compose_default_output_filepath(args[0])
    else:
        input = args[0]
        output = args[1]

    muscle_align_codons(infilepath=input, outfilepath=output, protein_dir=None, rna=options.rna, muscle_path=options.muscle_path, muscle_options='')

if __name__ == "__main__":
    main()
