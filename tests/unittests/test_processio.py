import unittest
import dendropy
import pytest
import tempfile
from dendropy.utility import processio

class SessionTests(unittest.TestCase):
    def test_session(self):
        with tempfile.NamedTemporaryFile() as file, open(file.name, mode="w") as file_write, open(file.name, mode="r") as file_read:
            reader = processio.SessionReader(file_read)
            for i in range(3):  
                file_write.write(f"asdf{i}\n")
                file_write.flush()
                res = reader.read()
                assert res == f"asdf{i}\n"
        

if __name__ == "__main__":  
    unittest.main()