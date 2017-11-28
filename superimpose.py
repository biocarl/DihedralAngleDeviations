import sys,os
from Bio.PDB import *

"""
Carl Hauck,
created 11.2017
https://github.com/biocarl/DihedralAngleDeviations
"""

if len(sys.argv) <3:
            sys.exit()

#Vars
pdb_origin=sys.argv[1]
pdb_mutate=sys.argv[2]

#parse atoms of the structures
parser = PDBParser()
structure_origin = parser.get_structure(pdb_origin.split('.')[0],pdb_origin)
atoms_origin = structure_origin.get_atoms()

structure_mutate = parser.get_structure(pdb_mutate.split('.')[0],pdb_mutate)
atoms_mutate = structure_mutate.get_atoms()

#Superimpose the two structures

#Approach: Find region which is identical and superimpose the structures based on that region, should work in the most cases

atoms_mutate_subset = list()
atoms_origin_subset = list()

#Define Offset, sometime the modelled proteins lack a residue and positions shifted. TODO detect that offset
offset = 1

range_min = 10 #n same residues
range = 0
start = 0
stop = 0
position = -1


#apply offset
list1 = structure_origin.get_residues()
list2 = structure_mutate.get_residues()
for i in xrange(offset):
    list1.next()
#Get range where residues are the same
for x, y in zip(list1,list2):
        position+=1
        print x.get_resname() , y.get_resname(), ' at position:',position
        if x.get_resname() == y.get_resname() :
            atoms_origin_subset.extend( x.get_atom() )
            atoms_mutate_subset.extend( y.get_atom() )
            range +=1
            if(range == range_min):
                stop = position
                start = position+1 - range_min
                break
        else: 
            range=0
            atoms_mutate_subset = list()
            atoms_origin_subset = list()

print 'Picked range for base of superposition: start',start,'stop',stop

#print select atoms
#for x, y in zip(atoms_origin_subset,atoms_mutate_subset):
#    print x.get_name() , y.get_name()

#Superpostion on subset and re-orientation on whole structure
sup = Superimposer()
# Specify the atom lists
sup.set_atoms(atoms_origin_subset,atoms_mutate_subset)

# Apply rotation/translation to the moving atoms
sup.apply(atoms_mutate)


#Save both models to a joined structure

#model1 = structure_origin.get_models().next() #first model
#model2 = structure_mutate.get_models().next() #first model

structures = [structure_origin,structure_mutate]

io = PDBIO()
io.set_structure(structures[0])
io.save('model1.pdb',write_end=False)

io = PDBIO()
io.set_structure(structures[1])
io.save('model2.pdb',write_end=False)

#Merge to one file with Pymol(?) file specification
with file('model1.pdb', 'r') as original1: data1 = original1.read()
with file('model2.pdb', 'r') as original2: data2 = original2.read()

with file('comparison.pdb', 'w') as modified: modified.write("MODEL        1\n"+data1+"ENDMDL\n"+ "MODEL        1\n"+data2+"ENDMDL\n")

#delete intermediate files
os.remove('model1.pdb')
os.remove('model2.pdb')

#OTHER
#resolution = structure.header['resolution']
#print(resolution)

