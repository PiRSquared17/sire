
from Sire.Move import *
from Sire.Mol import *

def test_exception():
    try:
        z = ZMatrix()
        z[ AtomIdx(0) ]
        assert(False)
    except:
        pass

if __name__ == "__main__":
    test_exception()
