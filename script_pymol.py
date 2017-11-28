import sys
from pymol import cmd,stored

"""
Carl Hauck,
created 11.2017
https://github.com/biocarl/DihedralAngleDeviations
"""

'''
Offset = if there is a shift
'''
def show_mutations(offset=0):
	obj_list = cmd.get_names('objects')

	length = len(obj_list)
	if(length < 2):
		print 'Error! Not two models!'
		sys.exit()

	#Declare models
	origin = obj_list[0]
	mutate = obj_list[1]

	#Get residues of structures
	residues_dic = { 'original_list' : [], 'mutant_list' : []}
	#populate original list
	cmd.iterate("("+"model "+origin+" and name ca)","original_list.append((resi,resn))",space=residues_dic)
	print origin,residues_dic['original_list']
	#populate mutant list
	cmd.iterate("("+"model "+mutate+" and name ca)","mutant_list.append((resi,resn))",space=residues_dic)
	print mutate,residues_dic['mutant_list']

	#Superimpose structure
	cmd.align(origin,mutate)

	#Compare residues:check for amino acid length
	if(abs(len(residues_dic['original_list']) -  len(residues_dic['mutant_list']))>1):
		print 'Unequal length of amino acids'#on is tolerated because of offset
		sys.exit()

	#Compare residues: find mutations
	mutations = list()
	select_origin = "("+"model "+origin+" and resi 0" #select all atoms from model1 with the following residues ids
	select_mutant = "("+"model "+mutate+" and resi 0"

	for x, y in zip(residues_dic['original_list'][int(offset):],residues_dic['mutant_list']):
		if(x[1] != y[1]):
			#print x,y,' is different.'
			mutations.append((x,y))
			select_origin += '+' + str(int(x[0]))
			select_mutant += '+' + str(int(y[0]))

	select_origin += ")"
	select_mutant += ")"


	#prepare view
	cmd.show("cartoon")

	#Create surface
	cmd.show("surface","(Model "+origin+")")
	cmd.show("surface","(Model "+mutate+")")

	cmd.set("transparency", 0.8)
	cmd.set("surface_color", "white", "*")

	#create objects of affected residues
	cmd.create("original",select_origin,0,1)
	cmd.create("mutate",select_mutant,0,1)

	#Change vis of altered residues
	cmd.color("red","(Model original)")
	cmd.color("blue","(Model mutate)")
	cmd.hide("everything","(Model mutate)")
	cmd.hide("everything","(Model original)")
	cmd.show("sticks","(Model original)")
	cmd.show("sticks","(Model mutate)")

cmd.extend('show_mutations',show_mutations)
