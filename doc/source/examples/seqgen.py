#! /usr/bin/env python

import dendropy
from dendropy.interop import seqgen

trees = dendropy.TreeList.get(
        path="pythonidae.mcmc.nex",
        schema="nexus")
s = seqgen.SeqGen()

# generate one alignment per tree
# as substitution model is not specified, defaults to a JC model
# will result in a DataSet object with one DnaCharacterMatrix per input tree
d0 = s.generate(trees)
print(len(d0.char_matrices))
print(d0.char_matrices[0].as_string("nexus"))

# instruct Seq-Gen to scale branch lengths by factor of 0.1
# note that this does not modify the input trees
s.scale_branch_lens = 0.1

# more complex model
s.char_model = seqgen.SeqGen.GTR
s.state_freqs = [0.4, 0.4, 0.1, 0.1]
s.general_rates = [0.8, 0.4, 0.4, 0.2, 0.2, 0.1]
d1 = s.generate(trees)
print(len(d0.char_matrices))
print(d0.char_matrices[0].as_string("nexus"))
