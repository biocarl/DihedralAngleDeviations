import sys,os,math
import matplotlib as mpl
mpl.use('Agg') #to work on no-x-server

import matplotlib.pyplot as plt
import numpy as np

from Bio.PDB import *

def main():
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

    #Find residues which differ (origin/mutate)
    residues_mutate_subset = list()
    residues_origin_subset = list()

    #Define Offset, sometime the modelled proteins lack a residue and positions shifted. TODO detect that offset
    offset = 0
    position = -1

    #apply offset
    list1 = structure_origin.get_residues()
    list2 = structure_mutate.get_residues()
    for i in xrange(offset):
        list1.next()

    #Get all residue pairs which are distinct from each other
    for x, y in zip(list1,list2):
        position+=1
        #consider only true residues
        if x.get_resname() == y.get_resname() and not x.get_id()[0].strip() and not y.get_id()[0].strip() : #empty strings are considered as falsy in python
            #print 'saved',x.get_resname() , y.get_resname(), ' at position:',position
            #print 'raw',x, y, ' at position:',position
            #print 'raw',x.get_id()[1], y.get_id()[1], ' at position:',position
            residues_origin_subset.append((x,position+offset)) #saving tuple: residue object --- position in array != id-position
            residues_mutate_subset.append((y,position))
        if position == 210 :
            break

    if(len(residues_origin_subset)==0):
        #print "No mutations found! Cancel."
        sys.exit()
    #1' Dihedral angles
    #Get dihedral angles-origin
    ppb=PPBuilder()
    polypeptides = ppb.build_peptides(structure_origin)
    phi_psi_origin = polypeptides[0].get_phi_psi_list()
    #print '#angles',len(phi_psi_origin)

    #Get dihedral angles-mutate
    ppb=PPBuilder()
    polypeptides = ppb.build_peptides(structure_mutate)
    phi_psi_mutate = polypeptides[0].get_phi_psi_list()
    #print '#angles',len(phi_psi_mutate)


    #save dihedral angles for distinct residues to file and plot
    phi_list_origin = list()
    phi_list_mutate = list()
    psi_list_origin = list()
    psi_list_mutate = list()
    dev_phi = list() #for later statistics
    dev_psi = list()


    alert = 'checking position: \n'
    alert+=residues_origin_subset[0][0].get_resname()+"("+str(residues_origin_subset[0][0].get_id()[1])+")\t" #as_origin
    alert+=residues_mutate_subset[0][0].get_resname()+"("+str(residues_mutate_subset[0][0].get_id()[1])+")\n" #as_mutate
    alert+="<"+residues_origin_subset[0][0].get_resname()+"("+str(residues_origin_subset[0][1])+")\t" #as_origin
    alert+="<"+residues_mutate_subset[0][0].get_resname()+"("+str(residues_mutate_subset[0][1])+")\n" #as_mutate
    #print alert

    table = "as_origin\tphi\tpsi\tas_mutate\tphi\tpsi\n" #tab-delimited table
    for x, y in zip(residues_origin_subset,residues_mutate_subset):
        #tmp values of list entries


        if(phi_psi_origin[x[1]][0] is not None and phi_psi_mutate[y[1]][0] is not None and phi_psi_origin[x[1]][1] is not None and phi_psi_mutate[y[1]][1] is not None ) :
            phi_o = math.degrees(float(phi_psi_origin[x[1]][0])) #x[1] position in sequence array
            phi_m = math.degrees(float(phi_psi_mutate[y[1]][0]))
            psi_o = math.degrees(float(phi_psi_origin[x[1]][1]))
            psi_m = math.degrees(float(phi_psi_mutate[y[1]][1]))
            #phi_o = math.degrees(float(phi_psi_origin[x[0].get_id()[1]-1][0]))
            #phi_m = math.degrees(float(phi_psi_mutate[y[0].get_id()[1]-1][0]))
            #psi_o = math.degrees(float(phi_psi_origin[x[0].get_id()[1]-1][1]))
            #psi_m = math.degrees(float(phi_psi_mutate[y[0].get_id()[1]-1][1]))
        else:
            #print "Dehidral angle imcomplete! Skipping"
            continue


        #print phi_o,phi_m, psi_o, psi_m

        #String for file write
        table+=x[0].get_resname()+"("+str(x[0].get_id()[1])+")\t" #as_origin
        table+=str(phi_o)+"\t"+str(psi_o)+"\t"#phi\tpsi
        table+=y[0].get_resname()+"("+str(y[0].get_id()[1])+")\t" #as_mutate
        table+=str(phi_m)+"\t"+str(psi_m)+"\n"#phi\tpsi
        #List for plotting
        phi_list_origin.append(phi_o)
        phi_list_mutate.append(phi_m)
        psi_list_origin.append(psi_o)
        psi_list_mutate.append(psi_m)
        #for statistics
        dev_phi.append(abs(phi_o-phi_m))
        dev_psi.append(abs(psi_o-psi_m))

    #Plot ramachandran
    #plot_ramachandran(phi_list_origin,psi_list_origin,phi_list_mutate,psi_list_mutate,"phi","psi")


    #Calculating some statics
        #mean angle deviation !float
    dev_phi_mean =  sum(dev_phi)/float(len(dev_phi))
    dev_psi_mean =  sum(dev_psi)/float(len(dev_psi))
    table += "mean deviation phi\t"+str(dev_phi_mean)+"\tmean deviation psi\t"+str(dev_psi_mean)+"\n"

        #standart deviation over all deviations
    dev_phi_stddev = np.std(dev_phi)
    dev_psi_stddev = np.std(dev_psi)
    table += "stddev of deviation phi\t"+str(dev_phi_stddev)+"\tstddev of deviation psi\t"+str(dev_psi_stddev)+"\n"

    #out = "mean-deviation_phi"+"\tstddev_of_deviation_phi"+"\t"+"mean_deviation_psi"+"\tstddev_of_deviation_psi"+"\n"
    out = str(dev_phi_mean) +"\t"+str(dev_phi_stddev)+"\t"+str(dev_psi_mean)+"\t"+str(dev_psi_stddev)
    print out

    #write to file
    #with open("phi_psi_angles.csv", "w") as text_file:
    #        text_file.write(table)




def plot_ramachandran(x1,y1,x2,y2,x_label_in,y_label_in):

    # Generate plot (with the help of http://azevedolab.net/resources/rama_plot_py.pdf)
    plt.plot(x1, y1, 'b.', x2, y2, 'r.')
    plt.xlim(-180,180) # Sets axis limits
    plt.ylim(-180,180) # Sets axis limits
    plt.xticks(np.arange(-180.1,180.1,30)) # Sets ticks markers for x axis
    plt.yticks(np.arange(-180.1,180.1,30)) # Sets ticks makers for y axis
    plt.xlabel(x_label_in) # Adds axis label
    plt.ylabel(y_label_in) # Adds axis label
    plt.arrow(-180,0,360,0) # Creates an arrow
    plt.arrow(0,-180,0,360) # Creates an arrow
    # Show plot
    #plt.show()
    # Create a file for plot
    fig = plt.gcf() # Creates a figure
    fig.set_size_inches(7.0, 7.0) # Changes figure size
    fig.savefig('rama.png', dpi=300) # Saves figures

if __name__ == "__main__":
    main()
