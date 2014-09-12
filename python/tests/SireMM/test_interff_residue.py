
from Sire.MM import *
from Sire.IO import *
from Sire.System import *
from Sire.Mol import *
from Sire.Units import *
from Sire.Qt import *
from Sire.Vol import *
from Sire.Move import *
from Sire.Maths import *
import Sire.Stream

import os

nmoves = 500

coul_cutoff = 10 * angstrom
lj_cutoff = 10 * angstrom

switchfunc = HarmonicSwitchingFunction(coul_cutoff,coul_cutoff,lj_cutoff,lj_cutoff)

if os.path.exists("../io/proteinbox.s3"):
    (mols, space) = Sire.Stream.load("../io/proteinbox.s3")
else:
    (mols, space) = Amber().readCrdTop("../io/proteinbox.crd", "../io/proteinbox.top")
    Sire.Stream.save( (mols,space), "../io/proteinbox.s3")

# space = Cartesian()

oldff = InterCLJFF("oldff")
oldff.setProperty("switchingFunction",switchfunc)
oldff.setProperty("space",space)
oldff.setShiftElectrostatics(True)

newff = InterFF("newff")
newff.setProperty("cljFunction", CLJShiftFunction())
newff.setProperty("switchingFunction",switchfunc)
newff.setProperty("space",space)
newff.disableParallelCalculation()

parff = InterFF("parff")
parff.setProperty("cljFunction", CLJShiftFunction())
parff.setProperty("space",space)
parff.setProperty("switchingFunction",switchfunc)
parff.enableParallelCalculation()

oldff.add(mols)
newff.add(mols)
parff.add(mols)

cljmols = CLJBoxes( CLJAtoms(mols.molecules()) )

def assert_equal( x, y ):
    if not (x == y):
        print("ERROR: %s is not equal to %s" % (x,y))
        assert(False)

def assert_almost_equal( x, y, delta = 0.5 ):
    if abs(x-y) > delta:
        print("ERROR: %s is not equal to %s within a delta of %s" % (x,y,delta))
        assert(False)

def test_energy(verbose = False):
    t = QElapsedTimer()

    cljfunc = CLJShiftFunction(coul_cutoff, lj_cutoff)
    cljfunc.setSpace(space)
    cljcalc = CLJCalculator()
    t.start()
    (rawcnrg,rawljnrg) = cljcalc.calculate(cljfunc,cljmols)
    ns = t.nsecsElapsed()

    if verbose:
        print("RAW CALC:  %s  %s  %s  : %s ms" % (rawcnrg+rawljnrg,rawcnrg,rawljnrg,
                                                  0.000001*ns))

    if verbose:
        print("Adding atoms to the forcefield...")

    oldff.energy()
    oldff.mustNowRecalculateFromScratch()
    newff.energy()
    newff.mustNowRecalculateFromScratch()
    parff.energy()
    parff.mustNowRecalculateFromScratch()

    t.start()
    old_nrg = oldff.energy()
    old_ns = t.nsecsElapsed()
    oldcnrg = oldff.energy(oldff.components().coulomb()).value()
    oldljnrg = oldff.energy(oldff.components().lj()).value()

    t.start()
    new_nrg = newff.energy()
    new_ns = t.nsecsElapsed()
    newcnrg = newff.energy(newff.components().coulomb()).value()
    newljnrg = newff.energy(newff.components().lj()).value()

    t.start()
    par_nrg = parff.energy()
    par_ns = t.nsecsElapsed()
    parcnrg = parff.energy(parff.components().coulomb()).value()
    parljnrg = parff.energy(parff.components().lj()).value()

    if verbose:
        print("\nTotal energy")
        print("OLD FF :  %s  %s  %s  : %s ms" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg,
                                                 0.000001*old_ns))
        print("NEW FF :  %s  %s  %s  : %s ms" % (newcnrg+newljnrg,newcnrg,newljnrg,
                                                 0.000001*new_ns))
        print("PAR FF :  %s  %s  %s  : %s ms" % (parcnrg+parljnrg,parcnrg,parljnrg,
                                                 0.000001*par_ns))

    assert_almost_equal( parcnrg, newcnrg )
    assert_almost_equal( parljnrg, newljnrg )
    assert_almost_equal( oldcnrg, newcnrg )
    assert_almost_equal( oldljnrg, newljnrg )

def test_moves(verbose = False):
    newsys = System()
    oldsys = System()
    optsys = System()
    parsys = System()

    res = MoleculeGroup("residues")
    protein = mols[ MolWithResID("ALA") ].molecule()

    for residue in protein.residues():
        res.add(residue)

    newsys.add(newff)
    oldsys.add(oldff)
    optsys.add(newff)
    parsys.add(parff)

    newsys.add(res)
    oldsys.add(res)
    optsys.add(res)
    parsys.add(res)

    moves = RigidBodyMC(res)
    moves.disableOptimisedMoves()

    moves.setMaximumTranslation( 0.2 * angstrom )
    moves.setMaximumRotation( 0.5 * degrees )

    opt_moves = RigidBodyMC(res)
    opt_moves.enableOptimisedMoves()

    opt_moves.setMaximumTranslation( 0.2 * angstrom )
    opt_moves.setMaximumRotation( 0.5 * degrees )

    t = QElapsedTimer()

    if verbose:
        print("Calculating initial energy...")

    t.start()
    oldnrg = oldsys.energy().value()
    oldns = t.nsecsElapsed()

    t.start()
    newnrg = newsys.energy().value()
    newns = t.nsecsElapsed()

    t.start()
    optnrg = optsys.energy().value()
    optns = t.nsecsElapsed()

    t.start()
    parnrg = parsys.energy().value()
    parns = t.nsecsElapsed()

    if verbose:
        print("\nTIMES: old = %s ms, new = %s ms, opt = %s ms, par = %s ms" % \
                        (oldns*0.000001,newns*0.000001,optns*0.000001,parns*0.000001))
        print("\nENERGIES: old = %s, new = %s, opt = %s, par = %s" % \
                        (oldnrg,newnrg,optnrg,parnrg))

    assert_almost_equal( oldnrg, newnrg, 0.5 )
    assert_almost_equal( optnrg, newnrg, 0.1 )
    assert_almost_equal( parnrg, newnrg, 0.1 )

    if verbose:
        print("\nPerforming simulation...")

    moves.clearStatistics()
    moves.setGenerator( RanGenerator(42) )
    t.start()
    moves.move(oldsys, nmoves, False)
    oldns = t.nsecsElapsed()
    old_naccept = moves.nAccepted()
    old_nreject = moves.nRejected()

    oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    oldljnrg = oldsys.energy( oldff.components().lj() ).value()

    oldsys.mustNowRecalculateFromScratch()    

    check_oldcnrg = oldsys.energy( oldff.components().coulomb() ).value()
    check_oldljnrg = oldsys.energy( oldff.components().lj() ).value()

    moves.clearStatistics()
    moves.setGenerator( RanGenerator(42) )
    t.start()
    moves.move(newsys, nmoves, False)
    newns = t.nsecsElapsed()
    new_naccept = moves.nAccepted()
    new_nreject = moves.nRejected()

    newcnrg = newsys.energy( newff.components().coulomb() ).value()
    newljnrg = newsys.energy( newff.components().lj() ).value()

    newsys.mustNowRecalculateFromScratch()    

    check_newcnrg = newsys.energy( newff.components().coulomb() ).value()
    check_newljnrg = newsys.energy( newff.components().lj() ).value()

    opt_moves.clearStatistics()
    opt_moves.setGenerator( RanGenerator(42) )
    t.start()
    opt_moves.move(optsys, nmoves, False)
    optns = t.nsecsElapsed()
    opt_naccept = moves.nAccepted()
    opt_nreject = moves.nRejected()

    optcnrg = optsys.energy( newff.components().coulomb() ).value()
    optljnrg = optsys.energy( newff.components().lj() ).value()

    optsys.mustNowRecalculateFromScratch()

    check_optcnrg = optsys.energy( newff.components().coulomb() ).value()
    check_optljnrg = optsys.energy( newff.components().lj() ).value()

    opt_moves.clearStatistics()
    opt_moves.setGenerator( RanGenerator(42) )
    t.start()
    opt_moves.move(parsys, nmoves, False)
    parns = t.nsecsElapsed()
    par_naccept = moves.nAccepted()
    par_nreject = moves.nRejected()

    parcnrg = parsys.energy( parff.components().coulomb() ).value()
    parljnrg = parsys.energy( parff.components().lj() ).value()

    parsys.mustNowRecalculateFromScratch()    

    check_parcnrg = parsys.energy( parff.components().coulomb() ).value()
    check_parljnrg = parsys.energy( parff.components().lj() ).value()

    if verbose:
        print("\nTIMES: old = %s ms, new = %s ms, opt = %s ms, par = %s ms" % \
                        (oldns*0.000001,newns*0.000001,optns*0.000001,parns*0.000001))

        print("\nNACCEPT: old = %s, new = %s, opt = %s, par = %s" % \
                        (old_naccept,new_naccept,opt_naccept,par_naccept))

        print("NREJECT: old = %s, new = %s, opt = %s, par = %s" % \
                        (old_nreject,new_nreject,opt_nreject,par_nreject))

        print("\nTotal energy")
        print("OLD FF :  %s  %s  %s" % (oldcnrg+oldljnrg,oldcnrg,oldljnrg))
        print("NEW FF :  %s  %s  %s" % (newcnrg+newljnrg,newcnrg,newljnrg))
        print("OPT FF :  %s  %s  %s" % (optcnrg+optljnrg,optcnrg,optljnrg))
        print("PAR FF :  %s  %s  %s" % (parcnrg+parljnrg,parcnrg,parljnrg))

        print("\nCheck energy")
        print("OLD FF :  %s  %s  %s" % (check_oldcnrg+check_oldljnrg,check_oldcnrg,check_oldljnrg))
        print("NEW FF :  %s  %s  %s" % (check_newcnrg+check_newljnrg,check_newcnrg,check_newljnrg))
        print("OPT FF :  %s  %s  %s" % (check_optcnrg+check_optljnrg,check_optcnrg,check_optljnrg))
        print("PAR FF :  %s  %s  %s" % (check_parcnrg+check_parljnrg,check_parcnrg,check_parljnrg))

    assert_equal( new_naccept, old_naccept )
    assert_equal( new_nreject, old_nreject )
    assert_equal( opt_naccept, old_naccept )
    assert_equal( opt_nreject, old_nreject )
    assert_equal( par_naccept, old_naccept )
    assert_equal( par_nreject, old_nreject )

    assert_almost_equal( oldcnrg, check_oldcnrg, 0.1 )
    assert_almost_equal( oldljnrg, check_oldljnrg, 0.1 )
    assert_almost_equal( newcnrg, check_newcnrg, 0.1 )
    assert_almost_equal( newljnrg, check_newljnrg, 0.1 )
    assert_almost_equal( optcnrg, check_optcnrg, 0.1 )
    assert_almost_equal( optljnrg, check_optljnrg, 0.1 )
    assert_almost_equal( parcnrg, check_parcnrg, 0.1 )
    assert_almost_equal( parljnrg, check_parljnrg, 0.1 )
    assert_almost_equal( newcnrg, oldcnrg, 0.5 )
    assert_almost_equal( newljnrg, oldljnrg, 0.5 )
    assert_almost_equal( optcnrg, newcnrg, 0.1 )
    assert_almost_equal( optljnrg, newljnrg, 0.1 )
    assert_almost_equal( parcnrg, newcnrg, 0.1 )
    assert_almost_equal( parljnrg, newljnrg, 0.1 )


if __name__ == "__main__":
    test_energy(True)
    test_moves(True)
