
from Sire.IO import *
from Sire.MM import *
from Sire.System import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.FF import *
from Sire.Move import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Qt import *

import os

amber = Amber()

(molecules, space) = amber.readCrdTop("test/io/waterbox.crd", "test/io/waterbox.top")

system = System()

swapwaters = MoleculeGroup("swapwaters")
waters = MoleculeGroup("waters")

molnums = molecules.molNums();

for molnum in molnums:
    water = molecules[molnum].molecule()

    if water.residue().number() == ResNum(2025):
        center_water = water

swapwaters.add(center_water)

center_point = center_water.evaluate().center()

for molnum in molnums:
    if molnum != center_water.number():
        water = molecules[molnum].molecule()

        if Vector.distance(center_point, water.evaluate().center()) < 7.5:
            water = water.residue().edit().setProperty("PDB-residue-name", "SWP").commit()
            swapwaters.add(water)
        else:
            waters.add(water)

system.add(swapwaters)
system.add(waters)

gridff = GridFF("gridff")
gridff.setBuffer(2 * angstrom)
gridff.setGridSpacing( 1 * angstrom )
gridff.setLJCutoff( 7.5 * angstrom )
gridff.setCoulombCutoff( 30 * angstrom )
gridff.setShiftElectrostatics(True)
#gridff.setUseAtomisticCutoff(True)
#gridff.setUseReactionField(True)

gridff.add(swapwaters, MGIdx(0))
gridff.add(waters, MGIdx(1))
gridff.setSpace( Cartesian() )

gridff2 = GridFF("gridff2")
gridff2.setBuffer(2*angstrom)
gridff2.setGridSpacing( 1 * angstrom )
gridff2.setLJCutoff( 7.5 * angstrom )
gridff2.setCoulombCutoff( 30 * angstrom )
gridff2.setShiftElectrostatics(True)
#gridff2.setUseAtomisticCutoff(True)
#gridff2.setUseReactionField(True)

gridff2.add( swapwaters, MGIdx(0) )
gridff2.addFixedAtoms(waters.molecules())
gridff2.setSpace( Cartesian() )

cljff = InterGroupCLJFF("cljff")
cljff.setSwitchingFunction( HarmonicSwitchingFunction(30*angstrom,30*angstrom,7.5*angstrom,7.5*angstrom) )
cljff.add(swapwaters, MGIdx(0))
cljff.add(waters, MGIdx(1))
cljff.setShiftElectrostatics(True)
#cljff.setUseAtomisticCutoff(True)
#cljff.setUseReactionField(True)
cljff.setSpace( Cartesian() )

print(gridff.energies())
print(gridff2.energies())

print("\nEnergies")
print(gridff.energies())
print(gridff2.energies())
print(cljff.energies())
