
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Move import *
from Sire.Maths import *
from Sire.Units import *
from Sire.CAS import *
from Sire.Vol import *
from Sire.Qt import *

from Sire.Tools import AmberLoader

import Sire.Stream

import os

coul_cutoff = 15 * angstrom
lj_cutoff = 15 * angstrom

nmoves = 50

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,
                                       lj_cutoff,lj_cutoff)

# Need to load with 0 1-4 scale factors because
# the new IntraFF forcefield does not account for 1-4
# terms
amber = Amber()
amber.set14Factors(0,0)

print("Loading protein box with zero 1-4 factors...")

if os.path.exists("../io/proteinbox_zero14.s3"):
    (mols, space) = Sire.Stream.load("../io/proteinbox_zero14.s3")
else:
    (mols,space) = amber.readCrdTop("../io/proteinbox.crd", "../io/proteinbox.top")
    Sire.Stream.save( (mols,space), "../io/proteinbox_zero14.s3" )

#space = Cartesian()

protein = mols[ MolWithResID("ALA") ].molecule()
ligand = mols[ MolWithResID("OSE") ].molecule()

flex = AmberLoader.generateFlexibility(ligand)
ligand = ligand.edit().setProperty("flexibility", flex).commit()

mols.update(ligand)

all = protein.selection()
part0 = protein.selection()
part1 = protein.selection()
part0 = part0.selectNone()
part1 = part1.selectNone()

for i in range(0, int(protein.nResidues()/2)):
    part0 = part0.select( ResIdx(i) )

for i in range(int(protein.nResidues()/2), protein.nResidues()):
    part1 = part1.select( ResIdx(i) )

protein0 = PartialMolecule(protein, part0)
protein1 = PartialMolecule(protein, part1)

print("Protein: Part 0 contains %d atoms, part 1 contains %d atoms, all contains %d atoms" \
            % (part0.nSelected(), part1.nSelected(), all.nSelected()))

all = ligand.selection()
part0 = ligand.selection()
part0 = part0.selectNone()
part1 = ligand.selection()
part1 = part1.selectNone()

for i in range(0, int(ligand.nAtoms()/2)):
    part0 = part0.select( AtomIdx(i) )

for i in range(int(ligand.nAtoms()/2),ligand.nAtoms()):
    part1 = part1.select( AtomIdx(i) )

ligand0 = PartialMolecule(ligand, part0)
ligand1 = PartialMolecule(ligand, part1)

print("Ligand: Part 0 contains %d atoms, part 1 contains %d atoms, all contains %d atoms" \
            % (part0.nSelected(), part1.nSelected(), all.nSelected()))

test_against_old = False

if test_against_old:
    oldff = IntraCLJFF("oldff")
    oldff.setSwitchingFunction(switchfunc)
    oldff.setSpace(space)
    oldff.setShiftElectrostatics(True)

    old0ff = IntraCLJFF("old0ff")
    old0ff.setSwitchingFunction(switchfunc)
    old0ff.setSpace(space)
    old0ff.setShiftElectrostatics(True)

    old1ff = IntraCLJFF("old1ff")
    old1ff.setSwitchingFunction(switchfunc)
    old1ff.setSpace(space)
    old1ff.setShiftElectrostatics(True)
else:
    oldff = IntraFF("oldff")
    oldff.setProperty("cljFunction", CLJIntraShiftFunction())
    oldff.setProperty("switchingFunction",switchfunc)
    oldff.setProperty("space",space)
    oldff.disableParallelCalculation()

    old0ff = IntraFF("old0ff")
    old0ff.setProperty("cljFunction", CLJIntraShiftFunction())
    old0ff.setProperty("switchingFunction",switchfunc)
    old0ff.setProperty("space",space)
    old0ff.disableParallelCalculation()

    old1ff = IntraFF("old1ff")
    old1ff.setProperty("cljFunction", CLJIntraShiftFunction())
    old1ff.setProperty("switchingFunction",switchfunc)
    old1ff.setProperty("space",space)
    old1ff.disableParallelCalculation()

newff = IntraGroupFF("newff")
newff.setProperty("cljFunction", CLJIntraShiftFunction())
newff.setProperty("switchingFunction",switchfunc)
newff.setProperty("space",space)
newff.disableParallelCalculation()

parff = IntraGroupFF("parff")
parff.setProperty("cljFunction", CLJIntraShiftFunction())
parff.setProperty("switchingFunction",switchfunc)
parff.setProperty("space",space)
parff.enableParallelCalculation()

print("Adding molecules to forcefields...")
oldff.add(ligand)
oldff.add(protein)
old0ff.add(ligand0)
old0ff.add(protein0)
old1ff.add(ligand1)
old1ff.add(protein1)

newff.add(ligand0, MGIdx(0))
newff.add(ligand1, MGIdx(1))
newff.add(protein0, MGIdx(0))
newff.add(protein1, MGIdx(1))
parff.add(ligand0, MGIdx(0))
parff.add(ligand1, MGIdx(1))
parff.add(protein0, MGIdx(0))
parff.add(protein1, MGIdx(1))


ligmol = MoleculeGroup("ligand")
ligmol.add(ligand)

residues = MoleculeGroup("residues")

for residue in protein.residues():
    residues.add(residue)

# need to calculate required energy using differences as old IntraCLJFF is buggy!
oldsys = System()
oldsys.add(ligmol)
oldsys.add(residues)
oldsys.add(oldff)
oldsys.add(old0ff)
oldsys.add(old1ff)

etotal = Symbol("E_{total}")
ecoul = Symbol("E_{coul}")
elj = Symbol("E_{lj}")

oldsys.setComponent( etotal, oldff.components().total() - \
                            old0ff.components().total() - \
                            old1ff.components().total() )

if oldsys.totalComponent() != etotal:
    oldsys.setComponent( oldsys.totalComponent(), etotal )

oldsys.setComponent( ecoul, oldff.components().coulomb() - \
                            old0ff.components().coulomb() - \
                            old1ff.components().coulomb() )
oldsys.setComponent( elj, oldff.components().lj() - \
                          old0ff.components().lj() - \
                          old1ff.components().lj() )

newsys = System()
newsys.add(ligmol)
newsys.add(residues)
newsys.add(newff)
newsys.setComponent( ecoul, newff.components().coulomb() )
newsys.setComponent( elj, newff.components().lj() )

parsys = System()
parsys.add(ligmol)
parsys.add(residues)
parsys.add(parff)
parsys.setComponent( ecoul, parff.components().coulomb() )
parsys.setComponent( elj, parff.components().lj() )

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

def _checkEnergies(oldsys, newsys, parsys, verbose):
    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()
    
    oldcnrg = oldsys.energy( ecoul ).value()
    oldljnrg = oldsys.energy( elj ).value()
    
    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()
    
    newcnrg = newsys.energy( ecoul ).value()
    newljnrg = newsys.energy( elj ).value()
      
    t.start()
    nrgs = parsys.energies()
    parns = t.nsecsElapsed()
    
    parcnrg = parsys.energy( ecoul ).value()
    parljnrg = parsys.energy( elj ).value()

    if verbose:
        print("OLD SYS :  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                  0.000001*oldns))
        print("NEW SYS :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                  0.000001*newns))
        print("PAR SYS :  %s  %s  %s  : %s ms" % (parcnrg+parljnrg,parcnrg,parljnrg,
                                                  0.000001*parns))

    assert_almost_equal( oldcnrg, newcnrg, 0.5 )
    assert_almost_equal( oldljnrg, newljnrg, 0.5 )
    assert_almost_equal( parcnrg, newcnrg, 0.5 )
    assert_almost_equal( parljnrg, newljnrg, 0.5 )


def test_energy(verbose = False):
    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

    _checkEnergies(oldsys, newsys, parsys, verbose)


def test_ligand_moves(verbose = False):

    moves = InternalMove(ligmol)
    moves.enableOptimisedMoves()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

    # return the systems to their original conformation
    oldsys.update(mols)
    newsys.update(mols)
    parsys.update(mols)

    if verbose:
        print("\nChecking total system energies...")

    _checkEnergies(oldsys, newsys, parsys, verbose)

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

    t.start()
    moves.move(parsys, nmoves, False)
    parns = t.nsecsElapsed()
    par_naccept = moves.nAccepted()
    par_nreject = moves.nRejected()
       
    if verbose:
        print("OLD: accept = %s, reject = %s, time = %s ms" % \
                   (old_naccept, old_nreject, oldns*0.000001))
        print("NEW: accept = %s, reject = %s, time = %s ms" % \
                   (new_naccept, new_nreject, newns*0.000001))
        print("PAR: accept = %s, reject = %s, time = %s ms" % \
                   (par_naccept, par_nreject, parns*0.000001))

        print("\nChecking final total energies...")

    _checkEnergies(oldsys, newsys, parsys, verbose)

    if verbose:
        print("\nRechecking total energies...")

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

    _checkEnergies(oldsys, newsys, parsys, verbose)

    assert_equal( old_naccept, new_naccept )
    assert_equal( old_naccept, par_naccept )
    assert_equal( old_nreject, new_nreject )
    assert_equal( old_nreject, par_nreject )


def test_protein_moves(verbose = False):

    moves = RigidBodyMC(residues)
    moves.setMaximumTranslation(0.25 * angstrom)
    moves.setMaximumRotation(0.5 * degrees)
    moves.enableOptimisedMoves()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

    # return the systems to their original conformation
    oldsys.update(mols)
    newsys.update(mols)
    parsys.update(mols)

    if verbose:
        print("\nChecking total system energies...")

    _checkEnergies(oldsys, newsys, parsys, verbose)

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

    t.start()
    moves.move(parsys, nmoves, False)
    parns = t.nsecsElapsed()
    par_naccept = moves.nAccepted()
    par_nreject = moves.nRejected()
       
    if verbose:
        print("OLD: accept = %s, reject = %s, time = %s ms" % \
                   (old_naccept, old_nreject, oldns*0.000001))
        print("NEW: accept = %s, reject = %s, time = %s ms" % \
                   (new_naccept, new_nreject, newns*0.000001))
        print("PAR: accept = %s, reject = %s, time = %s ms" % \
                   (par_naccept, par_nreject, parns*0.000001))

        print("\nChecking final total energies...")

    _checkEnergies(oldsys, newsys, parsys, verbose)

    if verbose:
        print("\nRechecking total energies...")

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

    _checkEnergies(oldsys, newsys, parsys, verbose)

    assert_equal( old_naccept, new_naccept )
    assert_equal( old_naccept, par_naccept )
    assert_equal( old_nreject, new_nreject )
    assert_equal( old_nreject, par_nreject )

if __name__ == "__main__":
    test_energy(True)
    test_ligand_moves(True)
    test_protein_moves(True)
