
from Sire.Vol import *

def test_exception():
    try:
        space = Cartesian()
        space.volume()
        assert(False)
    except:
        pass

if __name__ == "__main__":
    test_exception()
