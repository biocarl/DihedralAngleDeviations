#!/usr/bin/python
from Bio import SeqIO
import sys
if len(sys.argv) <3:
        sys.exit()


#Vars
fasta_origin=sys.argv[1]
fasta_mutate=sys.argv[2]
sequence_origin = []
sequence_mutate = []


#Split sequences of fasta into list, first entry
record_first = SeqIO.read(fasta_origin, "fasta")
sequence_origin = list(record_first.seq)

record_first = SeqIO.read(fasta_mutate, "fasta")
sequence_mutate = list(record_first.seq)

#Chancel if length is different, for now only point mutation is allowed. No deletion/insertion
if(len(sequence_mutate) != len(sequence_origin)) :
    print("Unequal fasta seq length!")
    sys.exit()


#Search for mutations and save their positions
position = -1
for x, y in zip(sequence_origin, sequence_mutate):
    position+=1
    if(x != y) :
        print x, y, ' at position:',position
