
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

nmoves = 1000

#space = Cartesian()

coul_cutoff = 15 * angstrom
lj_cutoff = 10 * angstrom

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff)

waters0 = MoleculeGroup("waters0")
waters1 = MoleculeGroup("waters1")

for i in range(0, mols.nMolecules()):
    mol = mols[ MolIdx(i) ].molecule()

    if i % 2 == 0:
        waters0.add(mol)
    else:
        waters1.add(mol)

oldff = InterGroupCLJFF("oldff")
oldff.setSwitchingFunction(switchfunc)
oldff.setSpace(space)
oldff.setShiftElectrostatics(True)

newff = InterGroupFF("newff")
newff.setProperty("switchingFunction",switchfunc)
newff.setProperty("space",space)
newff.disableParallelCalculation()

parff = InterGroupFF("parff")
parff.setProperty("switchingFunction",switchfunc)
parff.setProperty("space",space)
parff.enableParallelCalculation()

t = QElapsedTimer()

t.start()
oldff.add( waters0, MGIdx(0) )
oldff.add( waters1, MGIdx(1) )
oldns = t.nsecsElapsed()

t.start()
newff.add( waters0, MGIdx(0) )
newff.add( waters1, MGIdx(1) )
newns = t.nsecsElapsed()

t.start()
parff.add( waters0, MGIdx(0) )
parff.add( waters1, MGIdx(1) )
parns = t.nsecsElapsed()

t.start()
cljmols0 = CLJBoxes( CLJAtoms(waters0.molecules()) )
cljmols1 = CLJBoxes( CLJAtoms(waters1.molecules()) )
cljns = t.nsecsElapsed()

print("Setup times: OLD %s ms, NEW %s ms, PAR %s ms (%s ms)" % \
         (0.000001*oldns,0.000001*newns,
          0.000001*parns,0.000001*cljns))

def assert_almost_equal( x, y, delta = 0.5 ):
    if abs(x-y) > delta:
        print("ERROR: %s is not equal to %s within a delta of %s" % (x,y,delta))
        assert(False)

def assert_equal( x, y ):
    if not (x == y):
        print("ERROR: %s is not equal to %s" % (x, y))
        assert(False)

def test_energy(verbose = False):
    cljfunc = CLJShiftFunction(coul_cutoff, lj_cutoff)
    cljfunc.setSpace(space)
    cljcalc = CLJCalculator()
    t.start()
    (cnrg,ljnrg) = cljcalc.calculate(cljfunc,cljmols0,cljmols1)
    ns = t.nsecsElapsed()

    t.start()
    oldnrgs = oldff.energies()
    oldns = t.nsecsElapsed()
    oldcnrg = oldff.energy(oldff.components().coulomb()).value()
    oldljnrg = oldff.energy(oldff.components().lj()).value()

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()
    newcnrg = newff.energy(newff.components().coulomb()).value()
    newljnrg = newff.energy(newff.components().lj()).value()

    t.start()
    parnrgs = parff.energies()
    parns = t.nsecsElapsed()
    parcnrg = parff.energy(parff.components().coulomb()).value()
    parljnrg = parff.energy(parff.components().lj()).value()

    if verbose:
        print("\nTotal energy")
        print("CLJFUNC:  %s  %s  %s  : %s ms" % (cnrg+ljnrg,cnrg,ljnrg,
                                                 0.000001*ns))
        print("OLD FF :  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("PAR FF :  %s  %s  %s  : %s ms" % (parcnrg+parljnrg,parcnrg,parljnrg,
                                                 0.000001*parns))

    assert_almost_equal( cnrg, newcnrg )
    assert_almost_equal( ljnrg, newljnrg )
    assert_almost_equal( oldcnrg, newcnrg )
    assert_almost_equal( oldljnrg, newljnrg )
    assert_almost_equal( parcnrg, newcnrg )
    assert_almost_equal( parljnrg, newljnrg )

    if verbose:
        print("\nMoving a water molecule...\n")

    water = mols[ MolIdx(0) ].molecule()
    water = water.move().translate( Vector(1,0,0) ).commit()

    t.start()
    oldff.update(water)
    old_change = t.nsecsElapsed()

    t.start()
    newff.update(water)
    new_change = t.nsecsElapsed()

    t.start()
    parff.update(water)
    par_change = t.nsecsElapsed()

    if verbose:
        print("ff.update(): OLD %s ms, NEW %s ms, PAR %s ms" % \
           (old_change*0.000001,new_change*0.000001,par_change*0.000001))

    t.start()
    oldnrgs = oldff.energies()
    oldns = t.nsecsElapsed()
    oldcnrg = oldff.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldff.energy( oldff.components().lj() ).value()

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()
    newcnrg = newff.energy( newff.components().coulomb() ).value()
    newljnrg = newff.energy( newff.components().lj() ).value()

    t.start()
    parnrgs = parff.energies()
    parns = t.nsecsElapsed()
    parcnrg = parff.energy(parff.components().coulomb()).value()
    parljnrg = parff.energy(parff.components().lj()).value()

    if verbose:
        print("\nChanged energy")
        print("OLD FF :  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("PAR FF :  %s  %s  %s  : %s ms" % (parcnrg+parljnrg,parcnrg,parljnrg,
                                                 0.000001*parns))

    assert_almost_equal( oldcnrg, newcnrg )
    assert_almost_equal( oldljnrg, newljnrg  )
    assert_almost_equal( oldcnrg, parcnrg )
    assert_almost_equal( oldljnrg, parljnrg )

    oldff.mustNowRecalculateFromScratch()
    newff.mustNowRecalculateFromScratch()
    parff.mustNowRecalculateFromScratch()

    t.start()
    oldnrgs = oldff.energies()
    oldns = t.nsecsElapsed()

    r_oldcnrg = oldff.energy( oldff.components().coulomb() ).value()
    r_oldljnrg = oldff.energy( oldff.components().lj() ).value()

    t.start()
    newnrgs = newff.energies()
    newns = t.nsecsElapsed()

    r_newcnrg = newff.energy( newff.components().coulomb() ).value()
    r_newljnrg = newff.energy( newff.components().lj() ).value()

    t.start()
    parnrgs = parff.energies()
    parns = t.nsecsElapsed()

    r_parcnrg = parff.energy(parff.components().coulomb()).value()
    r_parljnrg = parff.energy(parff.components().lj()).value()

    if verbose:
        print("\nRecalculated energy")
        print("OLD FF :  %s  %s  %s  : %s ms" % (r_oldcnrg+r_oldljnrg,r_oldcnrg,r_oldljnrg,
                                                 0.000001*oldns))
        print("NEW FF :  %s  %s  %s  : %s ms" % (r_newcnrg+r_newljnrg,r_newcnrg,r_newljnrg,
                                                 0.000001*newns))
        print("PAR FF :  %s  %s  %s  : %s ms" % (r_parcnrg+r_parljnrg,r_parcnrg,r_parljnrg,
                                                 0.000001*parns))

    assert_almost_equal( oldcnrg, r_oldcnrg )
    assert_almost_equal( oldljnrg, r_oldljnrg )
    assert_almost_equal( newcnrg, r_newcnrg )
    assert_almost_equal( newljnrg, r_newljnrg )
    assert_almost_equal( parcnrg, r_parcnrg )
    assert_almost_equal( parljnrg, r_parljnrg )

def test_sim(verbose = False):

    oldsys = System()
    newsys = System()
    parsys = System()

    oldsys.add(mols)
    newsys.add(mols)
    parsys.add(mols)

    oldsys.add(oldff)
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

    oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldsys.energy( oldff.components().lj() ).value()

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

    moves = RigidBodyMC(mols)
    moves.setGenerator( RanGenerator( 42 ) )
    moves.enableOptimisedMoves()

    t.start()
    moves.move(oldsys, nmoves, False)
    move_oldns = t.nsecsElapsed()

    old_naccepted = moves.nAccepted()
    old_nrejected = moves.nRejected()

    moves.setGenerator( RanGenerator( 42 ) )
    moves.clearStatistics()

    t.start()
    moves.move(newsys, nmoves, False)
    move_newns = t.nsecsElapsed()

    new_naccepted = moves.nAccepted()
    new_nrejected = moves.nRejected()

    moves.setGenerator( RanGenerator( 42 ) )
    moves.clearStatistics()

    t.start()
    moves.move(parsys, nmoves, False)
    move_parns = t.nsecsElapsed()

    par_naccepted = moves.nAccepted()
    par_nrejected = moves.nRejected()

    t.start()
    nrgs = oldsys.energies()
    oldns = t.nsecsElapsed()
    
    t.start()
    nrgs = newsys.energies()
    newns = t.nsecsElapsed()

    t.start()
    nrgs = parsys.energies()
    parns = t.nsecsElapsed()

    oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldsys.energy( oldff.components().lj() ).value()
    
    newcnrg = newsys.energy( newff.components().coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value()

    parcnrg = parsys.energy( parff.components().coulomb() ).value()
    parljnrg = parsys.energy( parff.components().lj() ).value()

    if verbose:
        print("\nMoves: OLD %s ms, NEW %s ms, PAR %s ms" % \
            (0.000001*move_oldns, 0.000001*move_newns, 0.000001*move_parns))
        print("OLD SYS:  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*oldns))
        print("nAccepted() = %s, nRejected() = %s" % (old_naccepted, old_nrejected))
        print("NEW SYS:  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*newns))
        print("nAccepted() = %s, nRejected() = %s" % (new_naccepted, new_nrejected))
        print("PAR SYS:  %s  %s  %s  : %s ms" % (parcnrg+parljnrg,parcnrg,parljnrg,
                                                 0.000001*parns))
        print("nAccepted() = %s, nRejected() = %s" % (par_naccepted, par_nrejected))
    
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

    r_oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    r_oldljnrg = oldsys.energy( oldff.components().lj() ).value()

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


if __name__ == "__main__":
    test_energy(True)
    test_sim(True)
