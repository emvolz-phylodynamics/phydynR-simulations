#!/usr/bin/python
import subprocess
from pylab import *
from multiprocessing import Pool 

#  python sir.py 5 5 5000 100

DIR = 'sim0'
NP = 40
nsims = 200
#~ NP = 4
#~ nsims = 3
wub = 9
wlb = 1

S0 = 5.e3
tfin = 100. 
initial_sample_time = 95. 
sampleSize = 500

try: 
	os.makedirs( DIR )
except:
	pass
#


def call_sim( w ): 
	wacute = w
	wrisk2 = wub - wacute 
	print subprocess.check_call( ['python', 'sir.py', repr(wacute), repr(wrisk2), repr(S0), repr(tfin+1.), DIR] )
#

def call_treesampler(x):
	try:
		subprocess.check_call( ['Rscript', 'sample.tree.R', DIR+'/'+x, repr(sampleSize) , repr(initial_sample_time) ] ) # path n st0
	except:
		pass
#

p = Pool( NP )
p.map( call_sim, linspace(wlb, wub-1., nsims).tolist() )


import os,re
x = [fn for fn in os.listdir(DIR) if 'traj' in fn]
m = [ re.findall('[0-9]+', xx )[0]  for xx in x]
p.map( call_treesampler, m )
