
from Sire.MM import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.IO import *
from Sire.Units import *

(mols, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

natoms = 0

coul_cutoff = 100 * angstrom
lj_cutoff = 100 * angstrom

for i in range(0,mols.nMolecules()):
    natoms += mols[MolIdx(i)].molecule().nAtoms()

def assert_almost_equal( x, y, delta = 0.5 ):
    if abs(x-y) > delta:
        print("ERROR: %s is not equal to %s within a delta of %s" % (x,y,delta))
        assert(False)

def assert_equal( x, y ):
    if not (x == y):
        print("ERROR: %s is not equal to %s" % (x, y))
        assert(False)

def test_adding(verbose = False):

    group = CLJGroup()

    assert(group.isEmpty())
    assert(group.molecules().isEmpty())
    assert(group.needsAccepting())
    assert(group.changedAtoms().isEmpty())
    assert(group.recalculatingFromScratch())

    group.add(mols)

    assert(group.molecules().nMolecules() == mols.nMolecules())
    assert(group.needsAccepting())

    group.accept()

    assert(not group.needsAccepting())
    assert(not group.recalculatingFromScratch())

    assert( group.cljBoxes().nOccupiedBoxes() > 0 )
    assert( group.cljBoxes().nAtoms() >= natoms )

    cljfunc = CLJShiftFunction(coul_cutoff,lj_cutoff)
    cljcalc = CLJCalculator()

    cljboxes = CLJBoxes( CLJAtoms(mols) )

    if verbose:
        print("nOccupiedBoxes: %s vs. %s" % (cljboxes.nOccupiedBoxes(),
                        group.cljBoxes().nOccupiedBoxes()) )

        print("nAtoms(): %s vs. %s" % (cljboxes.nAtoms(),
                        group.cljBoxes().nAtoms()) )

    assert( cljboxes.nOccupiedBoxes() == group.cljBoxes().nOccupiedBoxes() )
    assert( cljboxes.nAtoms() == group.cljBoxes().nAtoms() )

    (cnrg0,ljnrg0) = cljcalc.calculate(cljfunc, group.cljBoxes())
    (cnrg1,ljnrg1) = cljcalc.calculate(cljfunc, cljboxes )

    if verbose:
        print("\nTest total energy")
        print("cnrg0 = %s, ljnrg0 = %s" % (cnrg0,ljnrg0))
        print("cnrg1 = %s, ljnrg1 = %s" % (cnrg1,ljnrg1))

    assert_almost_equal(cnrg0, cnrg1)
    assert_almost_equal(ljnrg0, ljnrg1)

    changed_mol = mols[MolIdx(0)].molecule().move().translate( Vector(1,0,0) ).commit()

    mols.update(changed_mol)

    group.update(changed_mol)

    assert( not group.recalculatingFromScratch() )
    assert( group.needsAccepting() )

    changed_atoms = group.changedAtoms()

    assert(not changed_atoms.isEmpty())
    # twice, as old + new atoms
    assert_equal( changed_atoms.nAtoms(), 2 * changed_mol.nAtoms() )

    assert_equal( group.nChangedMolecules(), 1 )

    (c_atoms, old_atoms, new_atoms) = group.mergeChanges()

    assert_equal( c_atoms, changed_atoms )
    assert_equal( old_atoms, group.oldAtoms() )
    assert_equal( new_atoms, group.newAtoms() )

    test_mol = group.changedMolecules()[changed_mol.number()].molecule()
    assert_equal( test_mol.number(), changed_mol.number() )
    assert_equal( test_mol.version(), changed_mol.version() )
    assert_equal( test_mol, changed_mol )

    cljboxes = CLJBoxes( CLJAtoms(mols) )

    if verbose:
        print("\nTest single move")
        print("nOccupiedBoxes: %s vs. %s" % (cljboxes.nOccupiedBoxes(),
                        group.cljBoxes().nOccupiedBoxes()) )

        print("nAtoms(): %s vs. %s" % (cljboxes.nAtoms(),
                        group.cljBoxes().nAtoms()) )

    assert( cljboxes.nOccupiedBoxes() == group.cljBoxes().nOccupiedBoxes() )
    assert( cljboxes.nAtoms() == group.cljBoxes().nAtoms() )

    (dcnrg0,dljnrg0) = cljcalc.calculate(cljfunc, CLJBoxes(changed_atoms), group.cljBoxes())
    (cnrg1,ljnrg1) = cljcalc.calculate(cljfunc, cljboxes )

    cnrg0 += dcnrg0
    ljnrg0 += dljnrg0

    if verbose:
        print("\nSingle move energies")
        print("cnrg0 = %s, ljnrg0 = %s" % (cnrg0,ljnrg0))
        print("cnrg1 = %s, ljnrg1 = %s" % (cnrg1,ljnrg1))

    assert_almost_equal(cnrg0, cnrg1)
    assert_almost_equal(ljnrg0, ljnrg1)

    assert( group.needsAccepting() )
    assert( not group.recalculatingFromScratch() )

    group.accept()

    assert( not group.needsAccepting() )
    assert( not group.recalculatingFromScratch() )

    (cnrg0,ljnrg0) = cljcalc.calculate(cljfunc, group.cljBoxes())

    if verbose:
        print("\nCompare to full recalculation")
        print("cnrg0 = %s, ljnrg0 = %s" % (cnrg0,ljnrg0))
        print("cnrg1 = %s, ljnrg1 = %s" % (cnrg1,ljnrg1))

    assert_almost_equal(cnrg0, cnrg1)
    assert_almost_equal(ljnrg0, ljnrg1)

    changed_mol0 = mols[ MolIdx(100) ].molecule().move().translate( Vector(0.5,0,0) ).commit()
    changed_mol1 = mols[ MolIdx(101) ].molecule().move().translate( Vector(0.2,0,0) ).commit()

    mols.update(changed_mol0)
    mols.update(changed_mol1)

    group.update(changed_mol0)

    assert( group.needsAccepting() )
    assert( not group.recalculatingFromScratch() )
    assert_equal( group.nChangedMolecules(), 1 )

    group.update(changed_mol1)

    assert( group.needsAccepting() )
    assert( not group.recalculatingFromScratch() )
    assert_equal( group.nChangedMolecules(), 2 )

    (changed_atoms, old_atoms, new_atoms) = group.mergeChanges()

    assert_equal( changed_atoms, group.changedAtoms() )
    assert_equal( old_atoms, group.oldAtoms() )
    assert_equal( new_atoms, group.newAtoms() )

    assert_equal( changed_atoms.nAtoms(), 2*(changed_mol0.nAtoms()+
                                             changed_mol1.nAtoms()) )

    changed_mols = group.changedMolecules()

    assert_equal( changed_mols[changed_mol0.number()].molecule(), changed_mol0 )
    assert_equal( changed_mols[changed_mol1.number()].molecule(), changed_mol1 )

    cljboxes = CLJBoxes( CLJAtoms(mols) )

    assert_equal( group.cljBoxes().nAtoms(),
                  cljboxes.nAtoms() - changed_mol0.nAtoms() - changed_mol1.nAtoms() )

    (dcnrg0,dljnrg0) = cljcalc.calculate(cljfunc, CLJBoxes(changed_atoms), group.cljBoxes())
    (cnrg1,ljnrg1) = cljcalc.calculate(cljfunc, cljboxes )

    cnrg0 += dcnrg0
    ljnrg0 += dljnrg0

    #also need the change in energy of mol0 and mol1
    (old_cnrg,old_ljnrg) = cljfunc.calculate(old_atoms)
    (new_cnrg,new_ljnrg) = cljfunc.calculate(new_atoms)

    print("delta %s %s %s %s" % (old_cnrg,old_ljnrg,new_cnrg,new_ljnrg))

    print("delta %s %s" % (new_cnrg-old_cnrg,new_ljnrg-old_ljnrg))

    cnrg0 += new_cnrg - old_cnrg
    ljnrg0 += new_ljnrg - old_ljnrg

    assert( group.needsAccepting() )
    group.accept()
    assert( not group.needsAccepting() )

    (cnrg2,ljnrg2) = cljcalc.calculate(cljfunc, group.cljBoxes() )

    (cnrg3,ljnrg3) = cljcalc.calculate(cljfunc, CLJBoxes( CLJAtoms(group.molecules()) ) )

    if verbose:
        print("\nAfter double move")
        print("cnrg0 = %s, ljnrg0 = %s" % (cnrg0,ljnrg0))
        print("cnrg1 = %s, ljnrg1 = %s" % (cnrg1,ljnrg1))
        print("cnrg2 = %s, ljnrg2 = %s" % (cnrg2,ljnrg2))
        print("cnrg3 = %s, ljnrg3 = %s" % (cnrg3,ljnrg3))

    assert_almost_equal(cnrg0, cnrg1)
    assert_almost_equal(ljnrg0, ljnrg1)

    assert_almost_equal(cnrg0, cnrg2)
    assert_almost_equal(ljnrg0, ljnrg2)

    assert_almost_equal(cnrg0, cnrg3)
    assert_almost_equal(ljnrg0, ljnrg3)

if __name__ == "__main__":
    test_adding(True)

