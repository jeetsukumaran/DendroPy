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
    # TODO: https://github.com/jeetsukumaran/DendroPy/issues/179#issue-1965884280
    # populations = FragmentedPopulations(1)
    # populations.generate_sequences("2X2") 
    pass

def test_pop_gen_tree():
    tree = dendropy.Tree()
    
    res = pop_gen_tree(tree=tree)
    assert isinstance(res, dendropy.Tree)
