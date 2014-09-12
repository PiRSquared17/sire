
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

from nose.tools import assert_almost_equal

(mols, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

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
newff.setProperty("switchingFunction",switchfunc)
newff.setProperty("space",space)
newff.disableParallelCalculation()

parff = InterFF("parff")
parff.setProperty("switchingFunction",switchfunc)
parff.setProperty("space",space)
parff.enableParallelCalculation()

t = QElapsedTimer()
t.start()
cljmols = CLJBoxes( CLJAtoms(mols.molecules()) )
cljns = t.nsecsElapsed()

newff.add(mols)
parff.add(mols)

def test_sim(verbose = False):

    oldsys = System()
    newsys = System()
    parsys = System()

    oldsys.add(mols)
    newsys.add(mols)
    parsys.add(mols)

    oldsys.add(newff)
    newsys.add(newff)
    parsys.add(parff)

    t = QElapsedTimer()

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()

    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    t.start()
    nrgs = parsys.energies()
    parns = t.nsecsElapsed()

    oldcnrg = oldsys.energy( newff.components().coulomb() ).value()
    oldljnrg = oldsys.energy( newff.components().lj() ).value()

    newcnrg = newsys.energy( newff.components().coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value()

    parcnrg = parsys.energy( parff.components().coulomb() ).value()
    parljnrg = parsys.energy( parff.components().lj() ).value()

    if verbose:
        print("\nStarting energy")
        print("OLD SYS:  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("PAR SYS:  %s  %s  %s  : %s ms" % (parcnrg+parljnrg,parcnrg,parljnrg,
                                                 0.000001*parns))

    assert_almost_equal( oldcnrg, newcnrg )
    assert_almost_equal( oldljnrg, newljnrg )

    assert_almost_equal( parcnrg, newcnrg )
    assert_almost_equal( parljnrg, newljnrg )

    moves = RigidBodyMC(mols)
    moves.disableOptimisedMoves()
    moves.setGenerator( RanGenerator( 42 ) )

    optmoves = RigidBodyMC(mols)
    optmoves.enableOptimisedMoves()
    optmoves.setGenerator( RanGenerator( 42 ) )

    parmoves = RigidBodyMC(mols)
    parmoves.enableOptimisedMoves()
    parmoves.setGenerator( RanGenerator( 42 ) )

    t.start()
    moves.move(oldsys, nmoves, False)
    move_oldns = t.nsecsElapsed()

    old_naccepted = moves.nAccepted()
    old_nrejected = moves.nRejected()

    t.start()
    optmoves.move(newsys, nmoves, False)
    move_newns = t.nsecsElapsed()

    new_naccepted = optmoves.nAccepted()
    new_nrejected = optmoves.nRejected()

    t.start()
    parmoves.move(parsys, nmoves, False)
    move_parns = t.nsecsElapsed()

    par_naccepted = parmoves.nAccepted()
    par_nrejected = parmoves.nRejected()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()
    
    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()
    
    t.start()
    nrgs = parsys.energies()
    parns = t.nsecsElapsed()

    oldcnrg = oldsys.energy( newff.components().coulomb() ).value()
    oldljnrg = oldsys.energy( newff.components().lj() ).value()
    
    newcnrg = newsys.energy( newff.components().coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value()
    
    parcnrg = parsys.energy( parff.components().coulomb() ).value()
    parljnrg = parsys.energy( parff.components().lj() ).value()

    if verbose:
        print("\nMoves: %s ms vs. %s ms vs. %s ms" % (0.000001*move_oldns, 0.000001*move_newns, 0.000001*move_parns))
        print("OLD SYS:  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("nAccepted() = %s, nRejected() = %s" % (old_naccepted, old_nrejected))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("nAccepted() = %s, nRejected() = %s" % (new_naccepted, new_nrejected))
        print("PAR SYS:  %s  %s  %s  : %s ms" % (parcnrg+parljnrg,parcnrg,parljnrg,
                                                 0.000001*parns))
        print("nAccepted() = %s, nRejected() = %s" % (par_naccepted, par_nrejected))
    
    assert_almost_equal( old_naccepted, new_naccepted )
    assert_almost_equal( old_nrejected, new_nrejected )
    assert_almost_equal( oldcnrg, newcnrg )
    assert_almost_equal( oldljnrg, newljnrg )

    assert_almost_equal( par_naccepted, new_naccepted )
    assert_almost_equal( par_nrejected, new_nrejected )
    assert_almost_equal( parcnrg, newcnrg )
    assert_almost_equal( parljnrg, newljnrg )

    oldsys.mustNowRecalculateFromScratch()
    newsys.mustNowRecalculateFromScratch()
    parsys.mustNowRecalculateFromScratch()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()
    
    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    t.start()
    nrgs = parsys.energies()
    parns = t.nsecsElapsed()

    r_oldcnrg = oldsys.energy( newff.components().coulomb() ).value()
    r_oldljnrg = oldsys.energy( newff.components().lj() ).value()

    r_newcnrg = newsys.energy( newff.components().coulomb() ).value()
    r_newljnrg = newsys.energy( newff.components().lj() ).value()

    r_parcnrg = parsys.energy( parff.components().coulomb() ).value()
    r_parljnrg = parsys.energy( parff.components().lj() ).value()

    if verbose:
        print("\nRecalculated energy")
        print("OLD SYS:  %s  %s  %s  : %s ms" % (r_oldcnrg+r_oldljnrg,r_oldcnrg,r_oldljnrg,
                                                 0.000001*oldns))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (r_newcnrg+r_newljnrg,r_newcnrg,r_newljnrg,
                                                 0.000001*newns))
        print("PAR SYS:  %s  %s  %s  : %s ms" % (r_parcnrg+r_parljnrg,r_parcnrg,r_parljnrg,
                                                 0.000001*parns))

    assert_almost_equal( r_oldcnrg, r_newcnrg )
    assert_almost_equal( r_oldljnrg, r_newljnrg )
    assert_almost_equal( r_parljnrg, r_newljnrg )

if __name__ == "__main__":
    test_sim(True)
