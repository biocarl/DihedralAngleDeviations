#!/usr/bin/python
from Bio import SeqIO
import sys
if len(sys.argv) == 1:
    sys.exit()


file=sys.argv[1]


handle = open(file, "rU")

for record in SeqIO.parse(handle, "pdb-seqres") :
    print ">" + record.id + "\n" + record.seq            
handle.close()

