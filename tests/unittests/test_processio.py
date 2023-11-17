import unittest
import dendropy
import pytest
import tempfile
from dendropy.utility import processio

class SessionTests(unittest.TestCase):
    @pytest.mark.asdf
    def test_session(self):
        with tempfile.NamedTemporaryFile() as file, open(file.name, mode="w") as filew, open(file.name, mode="r") as filer:
            reader = processio.SessionReader(filer)
            for i in range(3):  
                filew.write(f"asdf{i}\n")
                filew.flush()
                res = reader.read()
                assert res == f"asdf{i}\n"
        

if __name__ == "__main__":  
    unittest.main()