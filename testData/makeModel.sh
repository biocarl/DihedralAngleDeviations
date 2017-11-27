#!/bin/bash


reference_file=3BWM.pdb
sample_file=mutant.fasta


mutations=120
iterations=100

COUNTER=0
while [  $COUNTER -lt $iterations ]; do
echo The counter is $COUNTER
let COUNTER=COUNTER+1

#Produce input files for homoloy
python introduce_mutation.py $reference_file $mutations > mutant.fasta

#homology modeling
header=">P1;E2R2K7 \nsequence:E2RK7:::::::0.00: 0.00"

footer="*"

seq=$(
cat $sample_file | head -2  | tail -1
)


mv $reference_file 1m4r.pdb

printf "$header\n$seq$footer\n" > E2R2K7.txt

#model alignment
bash /home/chauck/miniconda2/bin/mod9.19 alinhar.py
#modelling
bash /home/chauck/miniconda2/bin/mod9.19 modelar.py

#renaming files
mv E2R2K7.B99990001.pdb model.pdb
mv 1m4r.pdb $reference_file

#delete files
rm *E2R2K7*
rm modelar.log alinhar.log

#analyse model
python analize_mutations_ori.py $reference_file  model.pdb >> deviations_ori.csv
python analize_mutations_mut.py $reference_file  model.pdb >> deviations_mut.csv


#clean
rm model.pdb

done

echo "DONE_JOB" > "DONE_JOB"
