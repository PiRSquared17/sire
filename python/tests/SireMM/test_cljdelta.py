
from Sire.IO import *
from Sire.MM import *
from Sire.Maths import *
from Sire.Mol import *

(mols,space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

cljboxes = CLJBoxes()
idxs = []

natoms = 0

for i in range(0,mols.nMolecules()):
    idxs.append( cljboxes.add( CLJAtoms(mols[MolIdx(i)]) ) )
    natoms += mols[MolIdx(i)].molecule().nAtoms()

def test_cljdelta(verbose = False):
    old_water = mols[MolIdx(0)].molecule()
    old_cljatoms = CLJAtoms(old_water)
    test_cljatoms = cljboxes.atoms(idxs[0])

    if verbose:
        print("\nMaking sure the box has the right number of atoms %s vs. %s",
               cljboxes.nAtoms(), natoms)

    assert( cljboxes.nAtoms() == natoms )

    if verbose:
        print("\nChecking I can get the old atoms back out again...")
        print("OLD:\n%s" % old_cljatoms)
        print("TEST:\n%s" % test_cljatoms)

    assert(old_cljatoms == test_cljatoms)

    new_water = old_water.move().translate( Vector(1) ).commit()
    new_cljatoms = CLJAtoms(new_water)

    cljdelta = CLJDelta(1, cljboxes.get(idxs[0]), CLJAtoms(new_water))

    test_cljatoms = cljdelta.newAtoms()

    if verbose:
        print("\nChecking changed atoms looks correct")
        print("CHANGED:\n%s" % cljdelta.changedAtoms())

        print("\nComparing new atoms are correctly in the delta")
        print("NEW:\n%s" % new_cljatoms)
        print("TEST:\n%s" % test_cljatoms)

    assert(new_cljatoms == test_cljatoms)

    test_cljatoms = cljdelta.oldAtoms()

    if verbose:
        print("\nComparing old atoms are correctly in the delta")
        print("OLD:\n%s" % old_cljatoms)
        print("TEST:\n%s" % test_cljatoms)

    assert(old_cljatoms == test_cljatoms)

if __name__ == "__main__":
    test_cljdelta(True)

