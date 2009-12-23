
from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *

import os

protomsdir = "%s/Work/ProtoMS" % os.getenv("HOME")

protoms = ProtoMS( "%s/protoms2" % protomsdir )

protoms.addParameterFile( "%s/parameter/amber99.ff" % protomsdir )
protoms.addParameterFile( "%s/parameter/gaff.ff" % protomsdir )
protoms.addParameterFile( "test/ff/sb2.ff" )

mol = PDB().readMolecule("test/io/sb2.pdb")

mol = mol.edit().rename("SB2").commit()

mol = protoms.parameterise(mol, ProtoMS.SOLUTE)

print mol.property("connectivity")

perturbations = mol.property("perturbations")

print perturbations

lam = perturbations.symbols().Lambda()

chgs = {}
ljs = {}

for i in range(0,101,10):
   mol = perturbations.perturb(mol, {lam:(0.01*i)})

   PDB().write(mol, "test_%003d.pdb" % i)

   chgs[i] = mol.property("charge")
   ljs[i] = mol.property("LJ")

keys = chgs.keys()
keys.sort()

print "initial", mol.property("initial_charge")

for key in keys:
    print key, chgs[key]

print "final", mol.property("final_charge")

for key in keys:
    print key, ljs[key]

