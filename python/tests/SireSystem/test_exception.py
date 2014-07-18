
from Sire.System import *
from Sire.Mol import *

def test_exception():
    try:
        s = System()
        s[MolName("test")]
        assert(False)
    except:
        pass

if __name__ == "__main__":
    test_exception()
