
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.System import *
from Sire.Move import *
from Sire.Maths import *
from Sire.Units import *
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

oldff = IntraCLJFF("oldff")
oldff.setSwitchingFunction(switchfunc)
oldff.setSpace(space)
oldff.setShiftElectrostatics(True)

newff = IntraFF("newff")
newff.setProperty("cljFunction", CLJIntraShiftFunction())
newff.setProperty("switchingFunction",switchfunc)
newff.setProperty("space",space)
newff.disableParallelCalculation()

parff = IntraFF("parff")
parff.setProperty("cljFunction", CLJIntraShiftFunction())
parff.setProperty("switchingFunction",switchfunc)
parff.setProperty("space",space)
parff.enableParallelCalculation()

print("Adding molecules to forcefields...")
oldff.add(protein)
oldff.add(ligand)
newff.add(protein)
newff.add(ligand)
parff.add(protein)
parff.add(ligand)

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
    
    oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldsys.energy( oldff.components().lj() ).value()
    
    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()
    
    newcnrg = newsys.energy( newff.components().coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value()
      
    t.start()
    nrgs = parsys.energies()
    parns = t.nsecsElapsed()
    
    parcnrg = parsys.energy( parff.components().coulomb() ).value()
    parljnrg = parsys.energy( parff.components().lj() ).value()

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
    oldff.mustNowRecalculateFromScratch()
    newff.mustNowRecalculateFromScratch()
    parff.mustNowRecalculateFromScratch()

    t.start()
    nrgs = oldff.energies()
    oldns = t.nsecsElapsed()

    oldcnrg = oldff.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldff.energy( oldff.components().lj() ).value()

    t.start()
    nrgs = newff.energies()
    newns = t.nsecsElapsed()

    newcnrg = newff.energy( newff.components().coulomb() ).value()
    newljnrg = newff.energy( newff.components().lj() ).value()

    t.start()
    nrgs = parff.energies()
    parns = t.nsecsElapsed()

    parcnrg = parff.energy( parff.components().coulomb() ).value()
    parljnrg = parff.energy( parff.components().lj() ).value()

    if verbose:
        print("\nTotal energy")
        print("OLD FF :  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("PAR FF :  %s  %s  %s  : %s ms" % (parcnrg+parljnrg,parcnrg,parljnrg,
                                                 0.000001*parns))

    assert_almost_equal( oldcnrg, newcnrg, 0.5 )
    assert_almost_equal( oldljnrg, newljnrg, 0.5 )
    assert_almost_equal( parcnrg, newcnrg, 0.5 )
    assert_almost_equal( parljnrg, newljnrg, 0.5 )

def test_ligand_moves(verbose = False):
    mols = MoleculeGroup("mols")
    mols.add(ligand)

    oldsys = System()
    newsys = System()
    parsys = System()

    oldsys.add(oldff)
    oldsys.add(mols)

    newsys.add(newff)
    newsys.add(mols)

    parsys.add(parff)
    parsys.add(mols)

    moves = InternalMove(mols)
    moves.enableOptimisedMoves()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

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

def test_protein_moves(verbose = False):
    mols = MoleculeGroup("mols")
    
    for residue in protein.residues():
        mols.add(residue)

    oldsys = System()
    newsys = System()
    parsys = System()

    oldsys.add(oldff)
    oldsys.add(mols)

    newsys.add(newff)
    newsys.add(mols)

    parsys.add(parff)
    parsys.add(mols)

    moves = RigidBodyMC(mols)
    moves.setMaximumTranslation(0.25 * angstrom)
    moves.setMaximumRotation(1 * degrees)
    moves.enableOptimisedMoves()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

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


if __name__ == "__main__":
    test_energy(True)
    test_ligand_moves(True)
    test_protein_moves(True)
