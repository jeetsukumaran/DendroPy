from dendropy.utility.error import ExternalServiceError

from . import marksmoke as pytestmark


def test_compose_message():
    err = ExternalServiceError("A", "A", "A", "A", "A", "A")

    res = err.compose_message()
    assert isinstance(res, str)
