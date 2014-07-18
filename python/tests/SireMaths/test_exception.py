
from Sire.Maths import *

def test_exception():
    try:
        v = NVector([1,2,3])
        v[4]
        assert(False)
    except:
        pass

if __name__ == "__main__":
    test_exception()
