from dendropy.calculate.probability import (
    binomial_rv,
    chisq_pdf,
    exp_pdf,
    geometric_rv,
    hypergeometric_pmf,
    num_poisson_events,
    poisson_pmf,
    poisson_rv,
    sample_multinomial,
    z_pmf,
)

from . import marksmoke as pytestmark


def test_binomial_rv():
    res = binomial_rv(2, 2)
    assert isinstance(res, int)


def test_exp_pdf():
    res = exp_pdf(1, 1)
    assert isinstance(res, float)


def test_poisson_rv():
    res = poisson_rv(1)
    assert isinstance(res, int)


def test_num_poisson_events():
    res = num_poisson_events(1, 1)
    assert isinstance(res, int)


def test_poisson_pmf():
    res = poisson_pmf(1, 1)
    assert isinstance(res, float)


def test_sample_multinomial():
    res = sample_multinomial([])
    assert isinstance(res, int)


def test_chisq_pdf():
    res = chisq_pdf(1, 1)
    assert isinstance(res, float)


def test_z_pmf():
    res = z_pmf(1)
    assert isinstance(res, float)


def test_geometric_rv():
    res = geometric_rv(1)
    assert isinstance(res, int)


def test_hypergeometric_pmf():
    res = hypergeometric_pmf(1, 1, 1, 1)
    assert isinstance(res, float)
