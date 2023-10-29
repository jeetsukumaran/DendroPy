import dendropy
from . import pytestmark

from dendropy.simulate.popgensim import FragmentedPopulations 

def test_popgenism():
    populations = FragmentedPopulations(1)
    # populations.generate_sequences("(A,(B,(C,D)));") # input to get_from_string() breaks; rest of file depends on this
