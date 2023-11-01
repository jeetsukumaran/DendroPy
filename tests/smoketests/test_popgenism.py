import dendropy
from dendropy.simulate.popgensim import FragmentedPopulations

from . import pytestmark


def test_popgenism():
    populations = FragmentedPopulations(1)
