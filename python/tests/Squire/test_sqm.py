
from Sire.IO import *
from Sire.Squire import *
from Sire.MM import *
from Sire.Mol import *
from Sire.System import *
from Sire.Move import *
from Sire.Maths import *
from Sire.Base import *

from nose.tools import assert_almost_equal

import os

# do we have sqm? (look in AMBERHOME)
try:
    print("AMBERHOME = %s" % os.getenv("AMBERHOME"))
    sqm_exe = findExe("%s/bin/sqm" % os.getenv("AMBERHOME")).absoluteFilePath()
    print("Using SQM executable %s" % sqm_exe)
except:
    sqm_exe = None
    print("Skipping SQM tests as a sqm executable cannot be found.")

if sqm_exe:
    (waters, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

def test_qmmm(verbose = False):
    if not sqm_exe:
        return

    qm_water = waters[MolIdx(0)].molecule()
    mm_water = waters[MolIdx(1)].molecule()

    dist = Vector.distance( qm_water.evaluate().center(), mm_water.evaluate().center() )

    mm_water = mm_water.move().translate( (dist-3.0) * (qm_water.evaluate().center() -
                                                        mm_water.evaluate().center()).normalise()
                                         ).commit()

    dist = Vector.distance( qm_water.evaluate().center(), mm_water.evaluate().center() )

    qmmmff = QMMMFF("qmmmff")
    sqm = SQM()
    sqm.setExecutable(sqm_exe)
    sqm.setMethod("AM1/d")
    qmmmff.setQuantumProgram(sqm)

    qmmmff.add(qm_water, MGIdx(0))

    zero_energy = qmmmff.energy()

    if verbose:
        print("Zero energy = %s" % zero_energy)

    qmmmff.setZeroEnergy(zero_energy)

    qmmmff.add(mm_water, MGIdx(1))

    if verbose:
        print(qmmmff.energies())

    # I know that the AM1/d energy should be 2.7859 kcal mol-1
    assert_almost_equal( qmmmff.energy().value(), 2.7859, 3 )

if __name__ == "__main__":
    test_qmmm(True)

