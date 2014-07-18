
from Sire.MM import *
from Sire.Mol import *

def test_exception():
    try:
        cljff = InterCLJFF()
        cljff[MolName("test")]
        assert(False)
    except:
        pass

if __name__ == "__main__":
    test_exception()
