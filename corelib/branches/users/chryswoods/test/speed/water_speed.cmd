
parfile1 solvents.ff
solvent1 ../io/water.pdb

streamwarning STDOUT
streamfatal STDOUT
streamheader STDOUT
streaminfo STDOUT

cutoff  15.0
feather 0.5

#SIRE_COMPATIBILITY_MODE on

boundary solvent

chunk1 singlepoint

