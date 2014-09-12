
from Sire.Mol import *
from Sire.MM import *
from Sire.IO import *
from Sire.Vol import *
from Sire.Qt import *
from Sire.Maths import *
from Sire.Units import *
from Sire.Base import *

coul_cutoff = 10 * angstrom
lj_cutoff = 10 * angstrom

dielectric = 78.0

combining_rules = "geometric"

(waterbox, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

switchfunc = HarmonicSwitchingFunction(coul_cutoff, coul_cutoff,
                                       lj_cutoff, lj_cutoff)

oldff = InterCLJFF("cljff")
oldff.add(waterbox)
oldff.setSpace(space)
oldff.setSwitchingFunction(switchfunc)
oldff.setCombiningRules(combining_rules)
oldff.setUseReactionField(True)
oldff.setReactionFieldDielectric(dielectric)

newff = InterFF("cljff")
cljfunc = CLJRFFunction()
cljfunc.setDielectric(dielectric)
newff.setCLJFunction(cljfunc)
newff.add(waterbox)
newff.setProperty("space",space)
newff.setProperty("switchingFunction", switchfunc)
newff.setProperty("combiningRules", wrap(combining_rules))
newff.enableParallelCalculation()

newff.energies()
oldff.energies()

t = QElapsedTimer()

def assert_almost_equal( x, y, delta = 0.5 ):
    if abs(x-y) > delta:
        print("ERROR: %s is not equal to %s within a delta of %s" % (x,y,delta))
        assert(False)

def assert_equal( x, y ):
    if not (x == y):
        print("ERROR: %s is not equal to %s" % (x, y))
        assert(False)

def _compareEnergies(oldff, newff, verbose):
     t.start()
     nrgs = oldff.energies()
     oldns = t.nsecsElapsed()

     oldcnrg = oldff.energy(oldff.components().coulomb()).value()
     oldljnrg = oldff.energy(oldff.components().lj()).value()

     t.start()
     newff.energies()
     newns = t.nsecsElapsed()

     newcnrg = newff.energy(newff.components().coulomb()).value()
     newljnrg = newff.energy(newff.components().lj()).value()

     if verbose:
         print("OLD: %s  %s  %s : %s ms" % (oldcnrg+oldljnrg, oldcnrg,
                                            oldljnrg, 0.000001*oldns))
         print("NEW: %s  %s  %s : %s ms" % (newcnrg+newljnrg, newcnrg,
                                            newljnrg, 0.000001*newns))

     assert_almost_equal(oldcnrg, newcnrg, 0.5)
     assert_almost_equal(oldljnrg, newljnrg, 0.5)

def _calc_energies(verbose):
     _compareEnergies(oldff, newff, verbose)

def test_geo_vac(verbose = False):
    oldff.setSpace(Cartesian())
    newff.setProperty("space", Cartesian())
    oldff.setCombiningRules("geometric")
    newff.setProperty("combiningRules", wrap("geometric"))

    print("\nGEOMETRIC - VACUUM")
    _calc_energies(verbose)

def test_ari_vac(verbose = False):
    oldff.setSpace(Cartesian())
    newff.setProperty("space", Cartesian())
    oldff.setCombiningRules("arithmetic")
    newff.setProperty("combiningRules", wrap("arithmetic"))

    print("\nARITHMETIC - VACUUM")
    _calc_energies(verbose)

def test_geo_box(verbose = False):
    oldff.setSpace(space)
    newff.setProperty("space", space)
    oldff.setCombiningRules("geometric")
    newff.setProperty("combiningRules", wrap("geometric"))

    print("\nGEOMETRIC - BOX")
    _calc_energies(verbose)

def test_ari_box(verbose = False):
    oldff.setSpace(space)
    newff.setProperty("space", space) 
    oldff.setCombiningRules("arithmetic")
    newff.setProperty("combiningRules", wrap("arithmetic"))
    
    print("\nARITHMETIC - BOX")
    _calc_energies(verbose)

if __name__ == "__main__":
    test_geo_vac(True)
    test_ari_vac(True)
    test_geo_box(True)
    test_ari_box(True)

