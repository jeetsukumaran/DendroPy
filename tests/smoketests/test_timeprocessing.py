from dendropy.utility.timeprocessing import pretty_timestamp, pretty_elapsed_datetime, parse_timedelta, pretty_timedelta
from . import marksmoke as pytestmark
import datetime

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