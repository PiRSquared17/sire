
from Sire.IO import *
from Sire.Mol import *

def test_exception():
    try:
        a = Amber()
        a.readCrdTop("test.crd", "test.top")
        assert(False)
    except:
        pass

if __name__ == "__main__":
    test_exception()
