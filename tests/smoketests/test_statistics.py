from dendropy.calculate.statistics import empirical_cdf, mode, quantile

from . import marksmoke as pytestmark


def test_mode():
    res = mode([1])
    assert isinstance(res, list)


def test_empirical_cdf():
    res = empirical_cdf([1], 0)
    assert isinstance(res, float)


def test_quantile():
    res1 = quantile([1.0], 1)
    res2 = quantile([1], 1)
    assert isinstance(res1, float)
    assert isinstance(res2, int)
