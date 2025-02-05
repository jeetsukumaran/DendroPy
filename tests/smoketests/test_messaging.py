import os

from dendropy.utility.messaging import ConsoleMessenger

from . import marksmoke as pytestmark


def test_error_leader():
    mess = ConsoleMessenger()

    res = mess.error_leader()
    isinstance(res, str)


def test_warning_leader():
    mess = ConsoleMessenger()

    res = mess.warning_leader()
    assert isinstance(res, str)


def test_info_leader():
    mess = ConsoleMessenger()

    res = mess.info_leader()
    assert isinstance(res, str)


def test_format_message():
    mess = ConsoleMessenger()

    res = mess.format_message("A", level=ConsoleMessenger.ERROR_MESSAGING_LEVEL)
    assert isinstance(res, str)


def test_log():
    mess = ConsoleMessenger()

    res = mess.log("A")
    assert res is None


def test_log_lines():
    mess = ConsoleMessenger()

    res = mess.log_lines("A", level=1)
    assert res is None


def test_error():
    with open(os.devnull, "w") as devnull:
        mess = ConsoleMessenger(dest=devnull)
        res = mess.error("A")
    assert res is None


def test_warning():
    with open(os.devnull, "w") as devnull:
        mess = ConsoleMessenger(dest=devnull)
        res = mess.warning("A")
    assert res is None


def test_info():
    with open(os.devnull, "w") as devnull:
        mess = ConsoleMessenger(dest=devnull)
        res = mess.info("A")
    assert res is None


def test_info_lines():
    with open(os.devnull, "w") as devnull:
        mess = ConsoleMessenger(dest=devnull)
        res = mess.info_lines("A")
    assert res is None


def test_info_raw():
    with open(os.devnull, "w") as devnull:
        mess = ConsoleMessenger(dest=devnull)
        res = mess.info_raw("A")
    assert res is None
