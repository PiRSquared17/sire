

from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Move import *
from Sire.Maths import *
from Sire.Units import *
from Sire.Base import *
from Sire.Vol import *
from Sire.Qt import *

from Sire.Tools import AmberLoader

import Sire.Stream

import os

coul_cutoff = 15 * angstrom
lj_cutoff = 10 * angstrom

nmoves = 1000

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,
                                       lj_cutoff,lj_cutoff)

# First load the molecules with their default scale factors
amber = Amber()
print("Loading protein box with default 1-4 factors...")

if os.path.exists("../io/proteinbox.s3"):
    (mols, space) = Sire.Stream.load("../io/proteinbox.s3")
else:
    (mols,space) = amber.readCrdTop("../io/proteinbox.crd", "../io/proteinbox.top")
    Sire.Stream.save( (mols,space), "../io/proteinbox.s3" )

protein = mols[ MolWithResID("ALA") ].molecule()
ligand = mols[ MolWithResID("OSE") ].molecule()

flex = AmberLoader.generateFlexibility(ligand)
ligand = ligand.edit().setProperty("flexibility", flex).commit()
mols.update(ligand)

# create two pairs of forcefields - the old IntraCLJFF with
# an InternalFF that does not calculate 1-4 interactions...
old_intraclj = IntraCLJFF("intraclj")
old_intraff = InternalFF("intraff")

# ...and a new IntraFF with an InternalFF that does calculate
# 1-4 interactions. The sum of these two pairs should be the same
new_intraclj = IntraFF("intraclj")
new_intraff = InternalFF("intraff")

old_intraclj.setProperty("switchingFunction", switchfunc)
old_intraclj.setProperty("space", space)
old_intraclj.setShiftElectrostatics(True)
old_intraclj.setCombiningRules("arithmetic")

old_intraff.disable14Calculation()

new_intraclj.setCLJFunction( CLJIntraShiftFunction() )
new_intraclj.setProperty("space", space)
new_intraclj.setProperty("switchingFunction", switchfunc)
new_intraclj.setProperty("combiningRules", wrap("arithmetic"))

new_intraff.enable14Calculation()
new_intraff.setProperty("combiningRules", wrap("arithmetic"))

old_intraclj.add(protein)
old_intraff.add(protein)
new_intraclj.add(protein)
new_intraff.add(protein)
old_intraclj.add(ligand)
old_intraff.add(ligand)
new_intraclj.add(ligand)
new_intraff.add(ligand)

# now create systems to hold the forcefields
oldsys = System()
newsys = System()

oldsys.add(old_intraclj)
oldsys.add(old_intraff)

newsys.add(new_intraclj)
newsys.add(new_intraff)

# now create molecule groups to be used for the moves
ligmol = MoleculeGroup("ligand")
ligmol.add(ligand)

residues = MoleculeGroup("residues")
for residue in protein.residues():
    residues.add(residue)

oldsys.add(ligmol)
oldsys.add(residues)
newsys.add(ligmol)
newsys.add(residues)

print("Performing tests...\n")
t = QElapsedTimer()

def assert_almost_equal( x, y, delta = 0.5 ):
    if abs(x-y) > delta:
        print("ERROR: %s is not equal to %s within a delta of %s" % (x,y,delta))
        assert(False)

def assert_equal( x, y ):
    if not (x == y):
        print("ERROR: %s is not equal to %s" % (x, y))
        assert(False)

def _printComponents(nrgs):
    keys = list(nrgs.keys())
    keys.sort()

    for key in keys:
        print("%s: %s kcal mol-1" % (key, nrgs[key]))

def _check_energy(verbose):
    t.start()
    oldnrg = oldsys.energy().value()
    oldns = t.nsecsElapsed()

    t.start()
    newnrg = newsys.energy().value()
    newns = t.nsecsElapsed()

    if verbose:
        print("old: %s   new: %s" % (oldnrg,newnrg))
        print("old: %s ms vs. %s ms" % (0.000001*oldns,0.000001*newns))

    assert_almost_equal(oldnrg, newnrg)

def test_energy(verbose = False):
    print("\nInitial check")
    _check_energy(verbose)

    print("\nOLD ENERGY COMPONENTS")
    _printComponents(oldsys.energies())
    print("\nNEW ENERGY COMPONENTS")
    _printComponents(newsys.energies())

    print("\nRecalculate check")
    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()

    _check_energy(verbose)

def test_ligand_moves(verbose = False):

    moves = InternalMove(ligmol)
    moves.enableOptimisedMoves()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()

    # return the systems to their original conformation
    oldsys.update(mols)
    newsys.update(mols)

    if verbose:
        print("\nChecking total system energies...")

    _check_energy(verbose)

    if verbose:
        print("\nChecking system simulation...")

    moves.setGenerator( RanGenerator(42) )
    moves.clearStatistics()

    t.start()
    moves.move(oldsys, nmoves, False)
    oldns = t.nsecsElapsed()
    old_naccept = moves.nAccepted()
    old_nreject = moves.nRejected()

    moves.setGenerator( RanGenerator(42) )
    moves.clearStatistics()

    t.start()
    moves.move(newsys, nmoves, False)
    newns = t.nsecsElapsed()
    new_naccept = moves.nAccepted()
    new_nreject = moves.nRejected()

    moves.setGenerator( RanGenerator(42) )
    moves.clearStatistics()
       
    if verbose:
        print("OLD: accept = %s, reject = %s, time = %s ms" % \
                   (old_naccept, old_nreject, oldns*0.000001))
        print("NEW: accept = %s, reject = %s, time = %s ms" % \
                   (new_naccept, new_nreject, newns*0.000001))

        print("\nChecking final total energies...")

    _check_energy(verbose)

    if verbose:
        print("\nRechecking total energies...")

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()

    _check_energy(verbose)

    assert_equal( old_naccept, new_naccept )
    assert_equal( old_nreject, new_nreject )


def test_protein_moves(verbose = False):

    moves = RigidBodyMC(residues)
    moves.setMaximumTranslation(0.05 * angstrom)
    moves.setMaximumRotation(0.1 * degrees)
    moves.enableOptimisedMoves()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()

    # return the systems to their original conformation
    oldsys.update(mols)
    newsys.update(mols)

    if verbose:
        print("\nChecking total system energies...")

    _check_energy(verbose)

    if verbose:
        print("\nChecking system simulation...")

    moves.setGenerator( RanGenerator(42) )
    moves.clearStatistics()

    t.start()
    moves.move(oldsys, nmoves, False)
    oldns = t.nsecsElapsed()
    old_naccept = moves.nAccepted()
    old_nreject = moves.nRejected()

    moves.setGenerator( RanGenerator(42) )
    moves.clearStatistics()

    t.start()
    moves.move(newsys, nmoves, False)
    newns = t.nsecsElapsed()
    new_naccept = moves.nAccepted()
    new_nreject = moves.nRejected()
       
    if verbose:
        print("OLD: accept = %s, reject = %s, time = %s ms" % \
                   (old_naccept, old_nreject, oldns*0.000001))
        print("NEW: accept = %s, reject = %s, time = %s ms" % \
                   (new_naccept, new_nreject, newns*0.000001))

        print("\nChecking final total energies...")

    _check_energy(verbose)

    if verbose:
        print("\nRechecking total energies...")

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()

    _check_energy(verbose)

    assert_equal( old_naccept, new_naccept )
    assert_equal( old_nreject, new_nreject )

if __name__ == "__main__":
    test_energy(True)
    test_ligand_moves(True)
    test_protein_moves(True)

