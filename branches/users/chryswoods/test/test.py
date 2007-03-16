
from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.MM import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.Units import *

import copy

timer = QTime()

#read in all of the molecules
print "Loading the molecules..."
timer.start()
mols = PDB().read("test/io/water.pdb")

ms = timer.elapsed()
print "... took %d ms" % ms

#specify the space in which the molecules are placed
space = Cartesian()

space = PeriodicBox(Vector(-18.3854,-18.66855,-18.4445), \
                    Vector( 18.3854, 18.66855, 18.4445))

#specify the type of switching function to use
switchfunc = HarmonicSwitchingFunction(80.0)
switchfunc = HarmonicSwitchingFunction(15.0, 14.5)

#create a forcefield for the molecules
cljff = InterCLJFF(space, switchfunc)

cljff2 = Tip4PFF(space, switchfunc)

ljff = InterLJFF(space, switchfunc)

ljff_a = InterLJFF(space, switchfunc)
ljff_b = InterLJFF(space, switchfunc)

ljff_a_b = InterGroupLJFF(space, switchfunc)

#parametise each molecule and add it to the forcefield
print "Parametising the molecules..."

chgs = AtomicCharges( [0.0, 0.52 * mod_electrons,
                            0.52 * mod_electrons,
                           -1.04 * mod_electrons] )

ljs = AtomicLJs( [ LJParameter( 3.15365 * angstrom, \
                                0.1550 * kcal_per_mol ), \
                   LJParameter.dummy(), \
                   LJParameter.dummy(), \
                   LJParameter.dummy() ] )

tip4p = False

timer.start()

rand = RanGenerator()

n_in_a = 0
n_in_b = 0

for mol in mols:
      
      mol.setProperty( "charges", chgs )
      mol.setProperty( "ljs", ljs )

      if (not tip4p):
          tip4p = mol

      #randomly divide the molecules into the two groups
      if (rand.rand() < 0.5):          
          ljff_a.add(mol, {ljff_a.parameters().lj() : "ljs"})
          ljff_a_b.addTo(ljff_a_b.groups().A(),
                         mol, {ljff_a_b.parameters().lj() : "ljs"})
                         
          n_in_a = n_in_a + 1
      else:
          ljff_b.add(mol, {ljff_b.parameters().lj() : "ljs"})
          ljff_a_b.addTo(ljff_a_b.groups().B(),
                         mol, {ljff_a_b.parameters().lj() : "ljs"})
      
          n_in_b = n_in_b + 1

      cljff.add(mol, {cljff.parameters().coulomb() : "charges",
                      cljff.parameters().lj() : "ljs"})
                      
      cljff2.add(mol, {cljff2.parameters().coulomb() : "charges",
                       cljff2.parameters().lj() : "ljs"})

      ljff.add(mol, {ljff.parameters().lj() : "ljs"})

ms = timer.elapsed()
print "... took %d ms" % ms
      
print "(%d molecules in group A, %d in group B)" % (n_in_a, n_in_b)
      
#now calculate the energy of the forcefield
print "Calculating the energy..."

timer.start()
nrg = cljff.energy()
ms = timer.elapsed()

print "InterCLJFF ",cljff.energy(), "kcal mol-1"
print "    Coulomb = ", cljff.energy(cljff.components().coulomb())
print "         LJ = ", cljff.energy(cljff.components().lj())

print "... took %d ms" % ms

timer.start()
nrg = ljff.energy()
ms = timer.elapsed()

print "InterLJFF ",ljff.energy(), "kcal mol-1"
print "         LJ = ", ljff.energy(ljff.components().lj())

print "... took %d ms" % ms

print "Calculating the energy..."

timer.start()
nrg = cljff2.energy()
ms = timer.elapsed()

print "Tip4PFF ",cljff2.energy(), "kcal mol-1"

print "... took %d ms" % ms

# sum up the partial forcefields
timer.start()
nrg = ljff_a.energy() + ljff_b.energy() + ljff_a_b.energy()
ms = timer.elapsed()

print "LJ_partials ",nrg,"kcal mol-1"
print "... took %d ms" % ms

timer.start()

nmoves = 1000
for i in range(0,nmoves):
    cljff.change( tip4p )
    nrg = cljff.energy()

ms = timer.elapsed()

print "InterCLJFF ",cljff.energy(), "kcal mol-1"
print "... took %d ms (%f moves per second)" % (ms, nmoves*1000.0/ms)

tip4p.translate( (1.0,0.0,0.0) )

timer.start()

cljff.change(tip4p)
ljff.change(tip4p)
cljff2.change(tip4p)

ms = timer.elapsed()

print "Changing took %d ms" % ms

timer.start()
nrg = cljff.energy()
ms = timer.elapsed()

print "InterCLJFF ",cljff.energy(), "kcal mol-1"
print "    Coulomb = ", cljff.energy(cljff.components().coulomb())
print "         LJ = ", cljff.energy(cljff.components().lj())

print "... took %d ms" % ms

timer.start()
nrg = ljff.energy()
ms = timer.elapsed()

print "InterLJFF ",ljff.energy(), "kcal mol-1"
print "         LJ = ", ljff.energy(ljff.components().lj())

print "... took %d ms" % ms

timer.start()
nrg = cljff2.energy()
ms = timer.elapsed()

print "Tip4PFF ",cljff2.energy(), "kcal mol-1"

print "... took %d ms" % ms

timer.start()

nmoves = 1000
for i in range(0,nmoves):
    tip4p.translate( (0.00001,0,0) )
    cljff.change( tip4p )
    nrg = cljff.energy()

ms = timer.elapsed()

print "%d moves of InterCLJFF took %d ms" % (nmoves, ms)

timer.start()

nmoves = 1000

delta = 0

old_nrg = ljff.energy()
old_version = ljff.version()

for i in range(0,nmoves):
    tip4p.translate( (0.00001,0,0) )
    
    old_ljff = copy.copy(ljff)
    
    ljff.change( tip4p )
    nrg = ljff.energy()

    if (nrg == old_nrg):
       print "Energies are wrongly the same!!!"
    
    if (ljff.version() == old_version):
       print "Versions are wrongly the same!!!"
    
    ljff = old_ljff
    
    nrg = ljff.energy()
    
    if (nrg != old_nrg):
       print "Energies are wrongly different!!!"

    if (ljff.version() != old_version):
       print "Versions are wrongly different!!!"

ms = timer.elapsed()

print "%d moves of InterLJFF took %d ms" % (nmoves, ms)

timer.start()

nmoves = 1000

delta = 0

old_nrg = ljff_a.energy() + ljff_b.energy() + ljff_a_b.energy()

for i in range(0,nmoves):
    tip4p.translate( (0.00001,0,0) )
    
    old_ljff_a = copy.copy(ljff_a)
    old_ljff_b = copy.copy(ljff_b)
    old_ljff_a_b = copy.copy(ljff_a_b)
    
    ljff_a.change( tip4p )
    ljff_b.change( tip4p )
    ljff_a_b.change( tip4p )
    
    nrg = ljff_a.energy() + ljff_b.energy() + ljff_a_b.energy()

    if (nrg == old_nrg):
       print "Energies are wrongly the same!!!"
    
    ljff_a = old_ljff_a
    ljff_b = old_ljff_b
    ljff_a_b = old_ljff_a_b
    
    nrg = ljff_a.energy() + ljff_b.energy() + ljff_a_b.energy()
    
    if (nrg != old_nrg):
       print "Energies are wrongly different!!!"

ms = timer.elapsed()

print "%d moves of LJ_partials took %d ms" % (nmoves, ms)
