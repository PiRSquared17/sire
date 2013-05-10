#############################
##
## The Sire Tools module
##
## This contains the python scripts that provide
## useful tools and shortcuts in Sire
##
## (C) Christopher Woods
##

#ensure that the SireQt and SireError libraries are loaded as
#these are vital for the rest of the module
import Sire.Qt
import Sire.Error

from Sire.Units import *

# Create a "readParams" function that reads a text file and returns a dictionary
# of key - value pairs
def readParams( filename ):
    """Read a text file containing key-value pairs and return as a dictionary"""
    lines = open(filename, "r").readlines()

    params = {}

    print "Reading in parameters..."

    for line in lines:
        # format is "key = value", "#" is used for commenting
        s = line.find("=")
        if s != -1:
            key = line[0:s].lstrip().rstrip()
            value = line[s+1:].lstrip().rstrip()

            if key.find("#") != -1:
                #ignore this line, as it must be a weird comment!
                next

            if value.find("#") != -1:
                value = value[0:value.find("#")].lstrip().rstrip()

            # now see if we have to turn the value from a string into a python object
            try:
                exec("testval = %s" % value, globals())
                value = testval
            except:
                pass

            print "%s = %s" % (key, value)
            params[key] = value

    print "\n",

    return params

#now import every tool, so that someone writing "import Sire.Tools" will
#get access to everything

# Module to simplify loading and setting up systems from Amber
# coordinate and topology files
import AmberLoader

# Module to perform waterswap calculations for protein-ligand
# binding free energies
import WSRC