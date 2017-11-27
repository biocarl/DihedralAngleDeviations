# DihedralAngleDeviations ~ *Impact of mutations on protein structure in homology modelling* 

Scripts for calculating/visualizing the change in dihedral angles of a sample relative to the template.

### *Part1*: Quantitative Analysis: How differ dihedral angles of models modelled through homology modelling when applying random mutations?

**Usage:**
*Analyse mutation effect quantitatively (plot,tables)*

**$** pymol analize_mutations.py 3BWM.pdb 3BWM_mutant.pdb

**Output:**
Mean deviation, stdtdev,
Ramachandran plot where (angles of reference: blue| angles of template red)
![plot](https://github.com/biocarl/DihedralAngleDeviations/raw/master/testData/results_example/rest_rama.png)


### *Part2*: Qualitative Analysis: See the different orientations of each mutated residue in the proposed model.

**Usage:**
*Analyse mutation effect visually: Open in Pymol and show mutation effect*

**$** pymol 3BWM.pdb 3BWM_mutant.pdb 
> run pymol_script.py 
> show_mutations 


###Possible workflow

**Extract sequence** 
$ python extract_sequence.py 3BWM.pdb > 3BWM.fasta

**Introduce mutation**
$ python introduce_mutation.py 3BWM.pdb 10 > 3BWM_mutant.fasta

**Show mutation**
$ python find_mutations.py 3BWM.fasta 3BWM_mutant.fasta

**Model mutated sequence with homology modelling e.g. Modeller or SWISS-MODEL**
3BWM_mutant.fasta -> 3BWM_mutant.pdb

**Analyse mutation effect quantitatively (plot,tables)**
$ pymol analize_mutations.py 3BWM.pdb 3BWM_mutant.pdb

**(Superimpose structures (->comparison.pdb) and merge as models)**
$ python superimpose.py 3BWM.pdb 3BWM_mutant.pdb

**Analyse mutation effect visually: Open in Pymol and show mutation effect**
$ pymol 3BWM.pdb 3BWM_mutant.pdb
> run pymol_script.py
> show_mutations

