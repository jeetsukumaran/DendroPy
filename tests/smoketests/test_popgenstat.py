import dendropy
from dendropy.calculate.popgenstat import (
    derived_state_matrix,
    unfolded_site_frequency_spectrum,
)

from . import pytestmark


def test_derived_state_matrix():
    m = dendropy.DnaCharacterMatrix.get_from_string(">s1", "fasta")

    res = derived_state_matrix(m)
    assert isinstance(res, dendropy.StandardCharacterMatrix)


def test_unfolded_site_frequency_spectrum():
    m = dendropy.DnaCharacterMatrix.get_from_string(">s1", "fasta")

    res = unfolded_site_frequency_spectrum(m)
    assert isinstance(res, dict)
