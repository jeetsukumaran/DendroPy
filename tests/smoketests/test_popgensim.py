from dendropy.simulate.popgensim import FragmentedPopulations, pop_gen_tree
import dendropy
from . import marksmoke as pytestmark


def test_generate_pop_tree():
    populations = FragmentedPopulations(1)

    res = populations.generate_pop_tree("")
    assert isinstance(res, dendropy.Tree)

def test_generate_gene_tree():
    pop = FragmentedPopulations(1)

    res = pop.generate_gene_tree("")
    assert isinstance(res, dendropy.Tree)

def test_generate_sequences():
    populations = FragmentedPopulations(2, use_seq_gen=False)
    populations.generate_sequences("cheetah")

def test_pop_gen_tree():
    tree = dendropy.Tree()

    res = pop_gen_tree(tree=tree)
    assert isinstance(res, dendropy.Tree)
