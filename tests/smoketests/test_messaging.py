from dendropy.utility.messaging import ConsoleMessenger
from . import marksmoke as pytestmark
import os

def test_error_leader():
    mess = ConsoleMessenger()
    mess.error_leader()

def test_warning_leader():
    mess = ConsoleMessenger()
    mess.warning_leader()

def test_info_leader():
    mess = ConsoleMessenger()
    mess.info_leader()

def test_format_message():
    mess = ConsoleMessenger()
    mess.format_message("A", level=ConsoleMessenger.ERROR_MESSAGING_LEVEL)

def test_log():
    mess = ConsoleMessenger()
    mess.log("A")

def test_log_lines():
    mess = ConsoleMessenger()
    mess.log_lines("A", level=1)

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