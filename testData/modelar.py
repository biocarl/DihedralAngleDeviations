from modeller import *
from modeller.automodel import *
env = environ()
a = automodel(env, alnfile='E2R2K7-1m4r.ali',
knowns='1m4rA', sequence='E2R2K7',
assess_methods=(assess.DOPE,
assess.GA341))
a.starting_model = 1
a.ending_model = 1
a.make()

