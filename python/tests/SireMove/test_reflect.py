
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *
from Sire.Move import *
from Sire.System import *
from Sire.Units import *
from Sire.Maths import *

(mols, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")

radius = 10*angstrom

waters = MoleculeGroup("waters")

for molnum in mols.molNums():
    mol = mols[molnum].molecule()

    if mol.evaluate().center().length() < radius.value():
        waters.add(mol)

ff = InterCLJFF("cljff")
ff.add(waters)

system = System()
system.add(ff)
system.add(waters)

def test_central(verbose=False):
    testsys = System(system)

    moves = RigidBodyMC(waters)
    moves.setReflectionSphere(Vector(0), radius)

    if verbose:
        print("Performing 5000 moves 10 times...")

        PDB().write(testsys.molecules(), "test_central0000.pdb")

    for i in range(1,11):
        moves.move(testsys, 5000, False)

        if verbose:
            PDB().write(testsys.molecules(), "test_central%0004d.pdb" % i)

    # the reflection sphere should ensure that none of the 
    # water molecules diffuse outside the sphere
    mols = testsys.molecules()

    for molnum in mols.molNums():
        mol = mols[molnum]

        dist = mol.evaluate().center().length()

        if verbose:
            print("%s : %s A" % (molnum.value(), dist))

        assert( dist <= radius.value() )


def test_local(verbose = False):
    testsys = System(system)

    moves = RigidBodyMC(waters)

    local_radius = 1 * angstrom

    if verbose:
        print("Setting local reflection spheres...")

    for molnum in waters.molNums():
        center = waters[molnum].evaluate().center()
        moves.setReflectionSphere(molnum, center, local_radius)

        if verbose:
            print("%s : %s == %s, %s == %s" % (molnum.value(), \
                center, moves.reflectionSphereCenter(molnum), \
                local_radius.value(), moves.reflectionSphereRadius(molnum).value()))

        assert( moves.reflectionSphereCenter(molnum) == center )
        assert( moves.reflectionSphereRadius(molnum) == local_radius )

    if verbose:                               
        print("Performing 5000 moves 10 times...")
        PDB().write(testsys.molecules(), "test_local0000.pdb")

    for i in range(1,11):
        moves.move(testsys, 5000, False)

        if verbose:
            PDB().write(testsys.molecules(), "test_local%0004d.pdb" % i)

    # the reflection sphere should ensure that none of the
    # water molecules diffuse outside the sphere
    mols = testsys.molecules()

    for molnum in mols.molNums():
        mol = mols[molnum]
        oldmol = waters[molnum]

        dist = Vector.distance( mol.evaluate().center(), oldmol.evaluate().center() )

        if verbose:
            print("%s : %s A" % (molnum.value(), dist))

        assert( dist <= local_radius.value() )

if __name__ == "__main__":
    test_central(True)
    test_local(True)
