
from Sire.MM import *
from Sire.FF import *
from Sire.System import *
from Sire.Move import *
from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Vol import *
from Sire.Base import *
from Sire.Units import *
from Sire.Qt import *

(mols, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

use_only_two_waters = False

if use_only_two_waters:
    newmols = MoleculeGroup("close")
    newmols.add( mols[MolIdx(0)] )

    for i in range(1, mols.nMolecules()):
        if Vector.distance( mols[MolIdx(0)].evaluate().center(),
                            mols[MolIdx(i)].evaluate().center() ) < 5:
            newmols.add( mols[MolIdx(i)] )
            break

    mols = newmols

grid_spacing = 0.5 * angstrom
grid_buffer = 2.0 * angstrom

nmoves = 1000

#space = Cartesian()

reflect_sphere_center = Vector(10)
reflect_sphere_radius = 10 * angstrom

long_coul_cutoff = 25 * angstrom
coul_cutoff = 15 * angstrom
lj_cutoff = 10 * angstrom

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff)
grid_switchfunc = HarmonicSwitchingFunction(long_coul_cutoff,long_coul_cutoff,
                                            lj_cutoff,lj_cutoff)

cluster = MoleculeGroup("cluster")
waters = MoleculeGroup("waters")

for i in range(0, mols.nMolecules()):
    mol = mols[ MolIdx(i) ].molecule()

    if Vector.distance(mol.evaluate().center(),reflect_sphere_center) < reflect_sphere_radius.value() :
        cluster.add(mol)
    else:
        waters.add(mol)

newff = InterFF("newff")
newff.setCLJFunction( CLJShiftFunction() )
newff.setCLJFunction( "test", CLJShiftFunction() )
newff.setProperty("switchingFunction",switchfunc)
newff.setProperty("space",space)
newff.enableParallelCalculation()

new_clusterff = InterFF("new_clusterff")
new_clusterff.enableParallelCalculation()
new_clusterff.setCLJFunction( CLJShiftFunction() )
new_clusterff.setCLJFunction( "test", CLJShiftFunction() )
new_clusterff.setProperty("switchingFunction", grid_switchfunc )
new_clusterff.setProperty("space", Cartesian())
new_clusterff.setProperty("gridSpacing", LengthProperty(grid_spacing))
new_clusterff.setProperty("gridBuffer", LengthProperty(grid_buffer))

newff.add(mols)

new_clusterff.add(cluster)
new_clusterff.setFixedAtoms(waters.molecules())
new_clusterff.setUseGrid(True)

t = QElapsedTimer()

def assert_almost_equal( x, y, delta = 0.5 ):
    if abs(x-y) > delta:
        print("ERROR: %s is not equal to %s within a delta of %s" % (x,y,delta))
        assert(False)

def assert_equal( x, y ):
    if not (x == y):
        print("ERROR: %s is not equal to %s" % (x, y))
        assert(False)

def test_energy(verbose = False):
    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()
    newcnrg = newff.energy(newff.components().coulomb()).value()
    newljnrg = newff.energy(newff.components().lj()).value()
    tstcnrg = newff.energy(newff.components().coulomb("test")).value()
    tstljnrg = newff.energy(newff.components().lj("test")).value()

    if verbose:
        print("\nTotal energy")
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("TST FF :  %s  %s  %s  : %s ms" % (tstcnrg+tstljnrg,tstcnrg,tstljnrg,
                                                 0.000001*newns))

    assert_almost_equal( tstcnrg, newcnrg )
    assert_almost_equal( tstljnrg, newljnrg )

    if verbose:
        print("\nMoving a water molecule...\n")

    water = mols[ MolIdx(0) ].molecule()
    water = water.move().translate( Vector(1,0,0) ).commit()

    t.start()
    newff.update(water)
    new_change = t.nsecsElapsed()

    if verbose:
        print("ff.update(): %s ms" % (new_change*0.000001))

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()

    newcnrg = newff.energy( newff.components().coulomb() ).value()
    newljnrg = newff.energy( newff.components().lj() ).value()
    tstcnrg = newff.energy( newff.components().coulomb("test") ).value()
    tstljnrg = newff.energy( newff.components().lj("test") ).value()

    if verbose:
        print("\nChanged energy")
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("TST FF :  %s  %s  %s  : %s ms" % (tstcnrg+tstljnrg,tstcnrg,tstljnrg,
                                                 0.000001*newns))

    assert_almost_equal( tstcnrg, newcnrg )
    assert_almost_equal( tstljnrg, newljnrg  )

    newff.mustNowRecalculateFromScratch()

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()

    r_newcnrg = newff.energy( newff.components().coulomb() ).value()
    r_newljnrg = newff.energy( newff.components().lj() ).value()

    r_tstcnrg = newff.energy( newff.components().coulomb("test") ).value()
    r_tstljnrg = newff.energy( newff.components().lj("test") ).value()

    if verbose:
        print("\nRecalculated energy")
        print("NEW FF :  %s  %s  %s  : %s ms" % (r_newcnrg+r_newljnrg,r_newcnrg,r_newljnrg,
                                                 0.000001*newns))
        print("TST FF :  %s  %s  %s  : %s ms" % (r_tstcnrg+r_tstljnrg,r_tstcnrg,r_tstljnrg,
                                                 0.000001*newns))

    assert_almost_equal( tstcnrg, r_tstcnrg )
    assert_almost_equal( tstljnrg, r_tstljnrg )
    assert_almost_equal( newcnrg, r_newcnrg )
    assert_almost_equal( newljnrg, r_newljnrg )

def test_grid_energy(verbose = False):
    t.start()
    newnrgs = new_clusterff.energies()
    newns = t.nsecsElapsed()
    newcnrg = new_clusterff.energy(new_clusterff.components().coulomb()).value()
    newljnrg = new_clusterff.energy(new_clusterff.components().lj()).value()
    tstcnrg = new_clusterff.energy(new_clusterff.components().coulomb("test")).value()
    tstljnrg = new_clusterff.energy(new_clusterff.components().lj("test")).value()

    if verbose:
        print("\nTotal energy")
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("TST FF :  %s  %s  %s  : %s ms" % (tstcnrg+tstljnrg,tstcnrg,tstljnrg,
                                                 0.000001*newns))

    assert_almost_equal( tstcnrg, newcnrg )
    assert_almost_equal( tstljnrg, newljnrg )

    if verbose:
        print("\nMoving a water molecule...\n")

    water = mols[ MolIdx(0) ].molecule()
    water = water.move().translate( Vector(1,0,0) ).commit()

    t.start()
    new_clusterff.update(water)
    new_change = t.nsecsElapsed()

    if verbose:
        print("ff.update(): %s ms" % (new_change*0.000001))

    t.start()
    newnrgs = new_clusterff.energies()
    newns = t.nsecsElapsed()

    newcnrg = new_clusterff.energy( new_clusterff.components().coulomb() ).value()
    newljnrg = new_clusterff.energy( new_clusterff.components().lj() ).value()
    tstcnrg = new_clusterff.energy( new_clusterff.components().coulomb("test") ).value()
    tstljnrg = new_clusterff.energy( new_clusterff.components().lj("test") ).value()

    if verbose:
        print("\nChanged energy")
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("TST FF :  %s  %s  %s  : %s ms" % (tstcnrg+tstljnrg,tstcnrg,tstljnrg,
                                                 0.000001*newns))

    assert_almost_equal( tstcnrg, newcnrg )
    assert_almost_equal( tstljnrg, newljnrg  )

    new_clusterff.mustNowRecalculateFromScratch()

    t.start()
    newnrgs = new_clusterff.energies()
    newns = t.nsecsElapsed()

    r_newcnrg = new_clusterff.energy( new_clusterff.components().coulomb() ).value()
    r_newljnrg = new_clusterff.energy( new_clusterff.components().lj() ).value()

    r_tstcnrg = new_clusterff.energy( new_clusterff.components().coulomb("test") ).value()
    r_tstljnrg = new_clusterff.energy( new_clusterff.components().lj("test") ).value()

    if verbose:
        print("\nRecalculated energy")
        print("NEW FF :  %s  %s  %s  : %s ms" % (r_newcnrg+r_newljnrg,r_newcnrg,r_newljnrg,
                                                 0.000001*newns))
        print("TST FF :  %s  %s  %s  : %s ms" % (r_tstcnrg+r_tstljnrg,r_tstcnrg,r_tstljnrg,
                                                 0.000001*newns))

    assert_almost_equal( tstcnrg, r_tstcnrg )
    assert_almost_equal( tstljnrg, r_tstljnrg )
    assert_almost_equal( newcnrg, r_newcnrg )
    assert_almost_equal( newljnrg, r_newljnrg )

if __name__ == "__main__":
    test_energy(True)
    test_grid_energy(True)

