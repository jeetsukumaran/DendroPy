import pytest
from . import pytestmark

from dendropy.simulate.popgensim import FragmentedPopulations, pop_gen_tree

def test_popgenism():
    populations = FragmentedPopulations(1)
    # populations.generate_sequences("(A,(B,(C,D)));") # input to get_from_string() breaks; rest of file depends on this
