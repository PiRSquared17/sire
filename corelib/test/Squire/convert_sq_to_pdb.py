
import sys

lines = open(sys.argv[1], "r").readlines()

#6  C  0.40322105  -0.87854830  -0.42860640
#8 O -8.54617705 -4.36022630 9.32893690 -0.83400000

#ATOM      4   N1 LIG     1       1.606   0.314  -2.127

natom = 0

for line in lines:
    words = line.split()

    if len(words) == 5:
        natom += 1
        print("ATOM %5d %4s QMA     1    %8.3f%8.3f%8.3f" % \
                  (natom, words[1], float(words[2]), float(words[3]),
                   float(words[4])) )
    elif len(words) == 6:
        natom += 1
        print("ATOM %5d %4s MMA     2    %8.3f%8.3f%8.3f" % \
                  (natom, words[1], float(words[2]), float(words[3]),
                   float(words[4])) )
