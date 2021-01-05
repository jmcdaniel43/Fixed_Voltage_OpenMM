#!/usr/bin/env python
import sys

# this shifts nanotubes by a set amount in the z-direction...
match1=[ 'nancB' , 'nancC' , 'nanhB' , 'nanhC' ]
match2=[ 'nancD' , 'nanhD' ]
shift_z1=-1.5
shift_z2=-3.0

pdbfile="test_no_electrolyte.pdb"

# read pdb
with open(pdbfile) as f:
    for line in f:
        if any(s in line for s in match1):
            # shift z coord of atom and print ...
            sys.stdout.write("%4s  %5s %4s %4s%s%4s    %8.3f%8.3f%8.3f  1.00  0.00      %4s%2s\n" %( line[0:4], line[6:11], line[12:16], line[17:21], line[21:22], line[22:26] , float(line[30:38]) , float(line[38:46]) , float(line[46:54]) + shift_z1 , line[72:76], line[76:78] ) )
        elif any(s in line for s in match2):
            # shift z coord of atom and print ...
            sys.stdout.write("%4s  %5s %4s %4s%s%4s    %8.3f%8.3f%8.3f  1.00  0.00      %4s%2s\n" %( line[0:4], line[6:11], line[12:16], line[17:21], line[21:22], line[22:26] , float(line[30:38]) , float(line[38:46]) , float(line[46:54]) + shift_z2 , line[72:76], line[76:78] ) )
        else:
            sys.stdout.write( line )

#atomdata.append( [ line[0:4] , line[6:11] , line[12:16] , line[17:21] , line[22:26] , line[30:38] , line[38:46] , line[46:54] , line[72:76] , line[76:78] ] )
#fout.write("%4s  %5d %4s %4s%s%4d    %8.3f%8.3f%8.3f  1.00  0.00      %4s%2s\n" %(data[0], iatom, data[2], data[3],'A', ires, x , y , z , data[8], data[9] ) )




             

