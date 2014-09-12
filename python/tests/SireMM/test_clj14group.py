
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

# First load the molecules with their default scale factors
amber = Amber()
print("Loading protein box with default 1-4 factors...")

if os.path.exists("../io/proteinbox.s3"):
    (mols14, space14) = Sire.Stream.load("../io/proteinbox.s3")
else:
    (mols14,space14) = amber.readCrdTop("../io/proteinbox.crd", "../io/proteinbox.top")
    Sire.Stream.save( (mols,space), "../io/proteinbox.s3" )

# Need to load with 0 1-4 scale factors so that
# we can calculate the difference in energy
amber = Amber()
amber.set14Factors(0,0)

print("Loading protein box with zero 1-4 factors...")

if os.path.exists("../io/proteinbox_zero14.s3"):
    (mols, space) = Sire.Stream.load("../io/proteinbox_zero14.s3")
else:
    (mols,space) = amber.readCrdTop("../io/proteinbox.crd", "../io/proteinbox.top")
    Sire.Stream.save( (mols,space), "../io/proteinbox_zero14.s3" )

protein14 = mols14[ MolWithResID("ALA") ].molecule()
ligand14 = mols14[ MolWithResID("OSE") ].molecule()

flex = AmberLoader.generateFlexibility(ligand14)
ligand14 = ligand14.edit().setProperty("flexibility", flex).commit()

protein = mols[ MolWithResID("ALA") ].molecule()
ligand = mols[ MolWithResID("OSE") ].molecule()

flex = AmberLoader.generateFlexibility(ligand)
ligand = ligand.edit().setProperty("flexibility", flex).commit()

old_protein14ff = IntraCLJFF("old14ff")
old_protein14ff.setSwitchingFunction(switchfunc)
old_protein14ff.setSpace(space)
old_protein14ff.setShiftElectrostatics(True)
old_protein14ff.setCombiningRules("arithmetic")
old_protein14ff.add(protein14)

old_proteinff = IntraCLJFF("oldff")
old_proteinff.setSwitchingFunction(switchfunc)
old_proteinff.setSpace(space)
old_proteinff.setShiftElectrostatics(True)
old_proteinff.setCombiningRules("arithmetic")
old_proteinff.add(protein)

old_ligand14ff = IntraCLJFF("old14ff")
old_ligand14ff.setSwitchingFunction(switchfunc)
old_ligand14ff.setSpace(space)
old_ligand14ff.setShiftElectrostatics(True)
old_ligand14ff.setCombiningRules("arithmetic")
old_ligand14ff.add(ligand14)

old_ligandff = IntraCLJFF("oldff")
old_ligandff.setSwitchingFunction(switchfunc)
old_ligandff.setSpace(space)
old_ligandff.setShiftElectrostatics(True)
old_ligandff.setCombiningRules("arithmetic")
old_ligandff.add(ligand)

protein_clj14group = CLJ14Group(protein14)
ligand_clj14group = CLJ14Group(ligand14)

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

def _test_energy(group, ff14, ff, verbose):
    ff14_cnrg = ff14.energy( ff14.components().coulomb() ).value()
    ff14_ljnrg = ff14.energy( ff14.components().lj() ).value()
    ff_cnrg = ff.energy( ff.components().coulomb() ).value()
    ff_ljnrg = ff.energy( ff.components().lj() ).value()

    oldcnrg = ff14_cnrg - ff_cnrg
    oldljnrg = ff14_ljnrg - ff_ljnrg

    t.start()
    (newcnrg, newljnrg) = group.energy()
    ns = t.nsecsElapsed()

    if verbose:
        print("\nCOMPARE")
        print("COULOMB OLD = %s, NEW = %s" % (oldcnrg,newcnrg))
        print("LJ      OLD = %s, NEW = %s" % (oldljnrg,newljnrg))
        print("TIME %s ms" % (0.000001*ns))

    assert_almost_equal( oldcnrg, newcnrg )
    assert_almost_equal( oldljnrg, newljnrg )

    # Now recalculate from scratch
    ff14.mustNowRecalculateFromScratch()
    ff.mustNowRecalculateFromScratch()
    group.mustNowRecalculateFromScratch()

    ff14_cnrg = ff14.energy( ff14.components().coulomb() ).value()
    ff14_ljnrg = ff14.energy( ff14.components().lj() ).value()
    ff_cnrg = ff.energy( ff.components().coulomb() ).value()
    ff_ljnrg = ff.energy( ff.components().lj() ).value()

    oldcnrg = ff14_cnrg - ff_cnrg
    oldljnrg = ff14_ljnrg - ff_ljnrg

    t.start()
    (newcnrg, newljnrg) = group.energy()
    ns = t.nsecsElapsed()

    if verbose:
        print("\nRECALCULATE FROM SCRATCH")
        print("COULOMB OLD = %s, NEW = %s" % (oldcnrg,newcnrg))
        print("LJ      OLD = %s, NEW = %s" % (oldljnrg,newljnrg))
        print("TIME %s ms" % (0.000001*ns))

    assert_almost_equal( oldcnrg, newcnrg )
    assert_almost_equal( oldljnrg, newljnrg )


def test_ligand_energy(verbose = False):
    if verbose:
        print("\nLIGAND 14 ENERGY")

    _test_energy(ligand_clj14group, old_ligand14ff, old_ligandff, verbose)

def test_protein_energy(verbose = False):
    if verbose:
        print("\nPROTEIN 14 ENERGY")

    _test_energy(protein_clj14group, old_protein14ff, old_proteinff, verbose)

def test_ligand_move(verbose = False):
    if verbose:
        print("\nLIGAND 14 SIMULATION")

    dihedrals = ligand14.property("connectivity").getDihedrals()
    rangen = RanGenerator(42)

    for i in range(0,10):
        dihedral = dihedrals[rangen.randInt(0,len(dihedrals)-1)]
        new_ligand14 = ligand14.move().change( dihedral, 5*degrees ).commit()
        new_ligand = ligand.move().change( dihedral, 5*degrees ).commit()

        ligand_clj14group.update(new_ligand14)
        old_ligand14ff.update(new_ligand14)
        old_ligandff.update(new_ligand)

        if verbose:
            print("MOVE %s" % (i+1))

        _test_energy(ligand_clj14group, old_ligand14ff, old_ligandff, verbose)

def test_protein_move(verbose = False):
    if verbose:
        print("\nPROTEIN 14 SIMULATION")

    nresidues = protein14.nResidues()
    rangen = RanGenerator(42)

    for i in range(0,10):
        residx = ResIdx(rangen.randInt(0,nresidues-1))
        new_protein14 = protein14.residue(residx).move().translate( Vector(1.0) ).commit().molecule()
        new_protein = protein.residue(residx).move().translate( Vector(1.0) ).commit().molecule()

        protein_clj14group.update(new_protein14)
        old_protein14ff.update(new_protein14)
        old_proteinff.update(new_protein)

        if verbose:
            print("MOVE %s" % (i+1))

        _test_energy(protein_clj14group, old_protein14ff, old_proteinff, verbose)


if __name__ == "__main__":
    test_ligand_energy(True)
    test_protein_energy(True) 
    test_ligand_move(True)
    test_protein_move(True)