from modeller import *
env = environ()
aln = alignment(env)
mdl = model(env, file='1m4r', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='1m4rA', atom_files='1m4r.pdb')
aln.append(file='E2R2K7.txt', align_codes='E2R2K7')
aln.align2d()
aln.write(file='E2R2K7-1m4r.ali', alignment_format='PIR')
