
import Sire.Stream

from Sire.MM import *
from Sire.Mol import *
from Sire.Vol import *
from Sire.Qt import *
from Sire.Maths import *
from Sire.Units import *
from Sire.IO import *

from nose.tools import assert_almost_equal

import sys
import os

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

group_space = space

protein = mols[ MolWithResID("ALA") ].molecule()

coul_cutoff = 15 * angstrom
lj_cutoff = 15 * angstrom

intraff = IntraCLJFF("old")
intraff.setSpace( Cartesian() )
intraff.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,lj_cutoff) )
intraff.setShiftElectrostatics( True )

intraff.add(protein)

intraff01 = IntraGroupCLJFF("old01")
intraff01.setSpace( Cartesian() )
intraff01.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,lj_cutoff) )
intraff01.setShiftElectrostatics( True )

nres = 1

select0 = protein.selection()
select0 = select0.selectNone()
select1 = protein.selection()
select1 = select1.selectNone()

for i in range(0,nres):
    select0 = select0.select( ResIdx(i) )

for i in range(nres, protein.nResidues()):
    select1 = select1.select( ResIdx(i) )

group0 = PartialMolecule( protein, select0 )
group1 = PartialMolecule( protein, select1 )

intraff01.add( group0, MGIdx(0) )
intraff01.add( group1, MGIdx(1) )

intrafunc = CLJIntraShiftFunction(coul_cutoff, lj_cutoff)
intrafunc.setConnectivity(protein)
intrafunc.setSpace( Cartesian() )

cljatoms = CLJAtoms(protein, CLJAtoms.USE_ATOMIDX)
cljboxes = CLJBoxes(cljatoms, 7.5*angstrom)

cljatoms0 = CLJAtoms(group0, CLJAtoms.USE_ATOMIDX)
cljatoms1 = CLJAtoms(group1, CLJAtoms.USE_ATOMIDX)

cljboxes0 = CLJBoxes(cljatoms0, 7.5*angstrom)
cljboxes1 = CLJBoxes(cljatoms1, 7.5*angstrom)

def pvt_compare_group(verbose):

    t = QElapsedTimer()
    t.start()
    intraff01.energies()
    old_ns = t.nsecsElapsed()

    old_cnrg = intraff01.energy( intraff01.components().coulomb() ).value()
    old_ljnrg = intraff01.energy( intraff01.components().lj() ).value()

    t.start()
    (new_cnrg, new_ljnrg) = intrafunc.calculate(cljatoms0, cljatoms1)
    new_ns = t.nsecsElapsed()

    cljcalculator = CLJCalculator()
    t.start()
    (box_cnrg, box_ljnrg) = cljcalculator.calculate(intrafunc, cljboxes0, cljboxes1)
    box_ns = t.nsecsElapsed()                

    if verbose:
        print("OLD: %s  %s  %s  %s ms" % (old_cnrg+old_ljnrg, old_cnrg, old_ljnrg, 0.000001*old_ns))
        print("NEW: %s  %s  %s  %s ms" % (new_cnrg+new_ljnrg, new_cnrg, new_ljnrg, 0.000001*new_ns))
        print("BOX: %s  %s  %s  %s ms" % (box_cnrg+box_ljnrg, box_cnrg, box_ljnrg, 0.000001*box_ns))

    assert_almost_equal( new_cnrg, old_cnrg, 2 )
    assert_almost_equal( new_ljnrg, old_ljnrg, 2 )
    assert_almost_equal( box_cnrg, old_cnrg, 2 )
    assert_almost_equal( box_ljnrg, old_ljnrg, 2 )
    

def pvt_compare(verbose):

    t = QElapsedTimer()
    t.start()
    intraff.energies()
    old_ns = t.nsecsElapsed()

    old_cnrg = intraff.energy( intraff.components().coulomb() ).value()
    old_ljnrg = intraff.energy( intraff.components().lj() ).value()

    t.start()
    (new_cnrg, new_ljnrg) = intrafunc.calculate(cljatoms)
    new_ns = t.nsecsElapsed()

    cljcalculator = CLJCalculator()
    t.start()
    (box_cnrg, box_ljnrg) = cljcalculator.calculate(intrafunc, cljboxes)
    box_ns = t.nsecsElapsed()

    if verbose:
        print("OLD: %s  %s  %s  %s ms" % (old_cnrg+old_ljnrg, old_cnrg, old_ljnrg, 0.000001*old_ns))
        print("NEW: %s  %s  %s  %s ms" % (new_cnrg+new_ljnrg, new_cnrg, new_ljnrg, 0.000001*new_ns))
        print("BOX: %s  %s  %s  %s ms" % (box_cnrg+box_ljnrg, box_cnrg, box_ljnrg, 0.000001*box_ns))
   
    assert_almost_equal( new_cnrg, old_cnrg, 2 )
    assert_almost_equal( new_ljnrg, old_ljnrg, 2 )
    assert_almost_equal( box_cnrg, old_cnrg, 2 )
    assert_almost_equal( box_ljnrg, old_ljnrg, 2 )


def test_compare_vacuum(verbose = False):
    intraff.setSpace( Cartesian() )
    intrafunc.setSpace( Cartesian() )

    intraff.setCombiningRules("arithmetic")
    intraff01.setCombiningRules("arithmetic")
    intrafunc.setCombiningRules( CLJFunction.ARITHMETIC )

    pvt_compare(verbose)

def test_compare_vacuum_geo(verbose = False):
    intraff.setSpace( Cartesian() )
    intrafunc.setSpace( Cartesian() )
    
    intraff.setCombiningRules("geometric")
    intraff01.setCombiningRules("geometric")
    intrafunc.setCombiningRules( CLJFunction.GEOMETRIC )
    
    pvt_compare(verbose)

def test_compare_box(verbose = False):
    intraff.setSpace( space )
    intrafunc.setSpace( space )

    intraff.setCombiningRules("arithmetic")
    intraff01.setCombiningRules("arithmetic")
    intrafunc.setCombiningRules( CLJFunction.ARITHMETIC )

    pvt_compare(verbose)

def test_compare_box_geo(verbose = False):
    intraff.setSpace( space )
    intrafunc.setSpace( space )

    intraff.setCombiningRules("geometric")
    intraff01.setCombiningRules("geometric")
    intrafunc.setCombiningRules( CLJFunction.GEOMETRIC )

    pvt_compare(verbose)

def test_compare_group(verbose = False):
    intraff01.setSpace( Cartesian() )
    intrafunc.setSpace( Cartesian() )

    intraff.setCombiningRules("arithmetic")
    intraff01.setCombiningRules("arithmetic")
    intrafunc.setCombiningRules( CLJFunction.ARITHMETIC )

    pvt_compare_group(verbose)

def test_compare_group_geo(verbose = False):
    intraff01.setSpace( Cartesian() )
    intrafunc.setSpace( Cartesian() )

    intraff.setCombiningRules("geometric")
    intraff01.setCombiningRules("geometric")
    intrafunc.setCombiningRules( CLJFunction.GEOMETRIC )

    pvt_compare_group(verbose)

def test_compare_group_box(verbose = False):
    intraff01.setSpace( group_space )
    intrafunc.setSpace( group_space )

    intraff.setCombiningRules("arithmetic")
    intraff01.setCombiningRules("arithmetic")
    intrafunc.setCombiningRules( CLJFunction.ARITHMETIC )

    pvt_compare_group(verbose)

def test_compare_group_box_geo(verbose = False):
    intraff01.setSpace( group_space )
    intrafunc.setSpace( group_space )
        
    intraff.setCombiningRules("geometric")
    intraff01.setCombiningRules("geometric")
    intrafunc.setCombiningRules( CLJFunction.GEOMETRIC )

    pvt_compare_group(verbose)

if __name__ == "__main__":

    print("\nGeometric combining rules")

    print("\nIntramolecular vacuum energy")
    test_compare_vacuum_geo(True)

    print("\nIntramolecular periodic boundaries energy")
    test_compare_box_geo(True)
        
    print("\nGroup intramolecular vacuum energy")
    test_compare_group_geo(True)
    
    print("\nGroup intramolecular periodic boundaries energy")
    test_compare_group_box_geo(True)

    print("\nArithmetic combining rules")

    print("\nIntramolecular vacuum energy")
    test_compare_vacuum(True)

    print("\nIntramolecular periodic boundaries energy")
    test_compare_box(True)

    print("\nGroup intramolecular vacuum energy")
    test_compare_group(True)

    print("\nGroup intramolecular periodic boundaries energy")
    test_compare_group_box(True)

