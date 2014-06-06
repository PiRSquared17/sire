
from Sire.Mol import *
from Sire.Vol import *
from Sire.IO import *

from nose.tools import assert_almost_equal

(mols, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

def _pvt_assert_maxval(map, maxval, verbose):

    i = -1

    for val in map.occupancy():
        i += 1

        if val > maxval:
            if verbose:
                print("Exceeded maxval %s at point %s on %s, %s" % \
                        (maxval, i, map.toString(), val))

            assert(False)

def _pvt_assert_contains(map, mol, expect_val, verbose):
    for atom in mol.atoms():
        coords = atom.property("coordinates")
        idx = map.gridInfo().pointToGridIndex(coords)
        val = map.occupancy()[ map.gridInfo().gridToArrayIndex(idx) ]
        
        if verbose:
            print("%s, %s: %s == %s (%s)" % (atom, coords, idx, val, expect_val))

        assert_almost_equal(val, expect_val, 2)

def test_volumemap(verbose = False):

    map = VolumeMap()

    for i in range(0,10):
        mol = mols[MolIdx(i)]
        map.add(mol)
        _pvt_assert_maxval(map, 0.1 + 1.0 / (i+1), verbose)

        for j in range(0,i+1):
            _pvt_assert_contains(map, mols[MolIdx(j)], 1.0 / (i+1), verbose)


if __name__ == "__main__":
    test_volumemap(True)

