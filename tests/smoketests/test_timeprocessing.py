import datetime

from dendropy.utility.timeprocessing import (
    parse_timedelta,
    pretty_elapsed_datetime,
    pretty_timedelta,
    pretty_timestamp,
)

from . import marksmoke as pytestmark


def test_pretty_timestamp():
    res = pretty_timestamp()
    assert isinstance(res, str)


def test_pretty_elapsed_datetime():
    res = pretty_elapsed_datetime(datetime.datetime(1, 1, 1))
    assert isinstance(res, str)


def test_parse_timedelta():
    res = parse_timedelta(datetime.timedelta())
    assert isinstance(res, tuple)


def test_pretty_timedelta():
    res = pretty_timedelta(datetime.timedelta())
    assert isinstance(res, str)
