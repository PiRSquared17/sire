
from Sire.Base import *

def test_exception():
    try:
        p = Properties()
        p.property("cat")
        assert(False)
    except:
        pass

if __name__ == "__main__":
    test_exception()
