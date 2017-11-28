#!/usr/bin/python
from Bio import SeqIO
from random import randint
import sys

"""
Carl Hauck,
created 11.2017
https://github.com/biocarl/DihedralAngleDeviations
"""
'''
Extracts seq from pdb file and introduces n mutations randomly

'''
if len(sys.argv) <3:
    print 'Wrong arguments'
    sys.exit()

file=sys.argv[1]
n=int(sys.argv[2])

handle = open(file, "rU")

for record in SeqIO.parse(handle, "pdb-seqres") :
    sequence = str(record.seq)
    length = len(sequence)
    if(length <= n):
        print 'Error, you want to introduce more mutations than the sequence has residues.'
        sys.exit()

    for _ in range(n):
        rand1 = randint(0,length-1)
        rand2 = randint(0,length-1)
        sequence = sequence[:rand1]+ sequence[rand2] +sequence[rand1+1:]
    print ">" + record.id + "\n" + sequence

handle.close()

