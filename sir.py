# acute/chronic + 2 risk levels + assortativity + aging 
# hosts age from younger level to older level 
#~ NOTE aging goes from rl 2 to rl1 

# usage: python sir.py RR_acute RR_highRiskGroup SusceptibleSize TimeFinal
# example: python sir.py 5 5 5000 100
from pylab import * 
import csv 
import numpy
import pdb 
import os, sys

print 'arguments '
print sys.argv 

wacute = float( sys.argv[1] )
wrisk2 = float( sys.argv[2] )
S0 = float(sys.argv[3] ) # 
tfin = float(sys.argv[4] )
#~ odir = sys.argv[5]
wchronic = 1.; 

# other parameters: 
prisk2 = .30 #approximate proportion new infections in risk level 2 (young age group)
assortprob = .80 # prob that transmission goes to same age group as donor 
gamma0 = 1.
gamma1 = 1/9. 
mu = 1/40.;
agerate = 1/5.
beta = .25 #R0=2



# output
pid = repr(os.getpid())

try:
	odir = sys.argv[5]
except IndexError:
	odir = '_'.join( [ 'sir' ] + sys.argv[1:]  )
#
try: 
	os.makedirs( odir)
except:
	pass
#
ofnstem = odir + '/' + pid
ofn1 = ofnstem + '-transm.csv'
ofn2 = ofnstem + '-rem.csv'
ofn3 = ofnstem + '-traj.csv'

###
# output parameters 
ofn4 = ofnstem + '-parms.csv'
csv.writer( open( ofn4, 'w' )).writerow( [ 'S0', 'tfin', 'wacute', 'wrisk2' ] )
csv.writer( open( ofn4, 'a' )).writerow( [ S0, tfin, wacute, wrisk2 ] )

###


out = open( ofn1, 'w')
writer = csv.writer( out ) # pid1, pid2, t
out2 = open( ofn2, 'w')
writer2 = csv.writer( out2 ) # pid, t, tiplab
out_traj = open( ofn3, 'w')
writer_traj = csv.writer( out_traj ) 

ic = 1;  

# s, i1, i2, i12, i22
x = [ S0, 1., 0. , 0., 0.]

pids1 = [ic]
pids2 = []
pids12 = [ic] #TODO ? 
pids22 = [] 

#events
PROG1 = 1
PROG2 = 2
DEATH1= 3
DEATH2= 4
DEATH_S = 5
TRANSM = 6
BIRTH = 7

PROG12 = 8
PROG22 = 9
DEATH12= 10
DEATH22 = 11

AGING1 = 12 # aging for early inf 
AGING2 = 13 # aging for later inf 


def state2rates(x):
	f = beta * (x[1] + x[2] + x[3] + x[4]) * x[0] / sum(x)
	r12 = x[1]*gamma0
	d1 = mu * x[1]
	d2 = mu * x[2] 
	r2d = x[2] * gamma1
	d_s = mu * x[0]
	birth= S0 * mu
	
	r12_2 = x[3] * gamma0
	d1_2 = x[3] * mu
	d2_2 = x[4] * mu 
	r2d_2 = x[4] * gamma1
	
	aging1 = agerate * x[3] # early inf
	aging2 = agerate * x[4] # later inf
	
	return [ r12, r2d, d1, d2, d_s, f, birth, r12_2, r2d_2, d1_2, d2_2 , aging1 , aging2 ]
#

def rates2event(rates):
	return numpy.random.choice( range(1,len(rates)+1) , p = array(rates) / sum(rates) )
#

def rates2deltat(rates):
	r = sum(rates)
	return exponential( 1./r )
#

def state2infectorType(pids1, pids2, pids12, pids22):
	inf1w = len(pids1) * wacute
	inf2w = len(pids2) * wchronic
	inf12w = len(pids12) * wacute * wrisk2
	inf22w = len(pids22) * wchronic * wrisk2
	infw = array([inf1w, inf2w, inf12w, inf22w ])
	return numpy.random.choice( range(1,5) , p = infw/ sum( infw ) ) 
#

t = 0.; 
it = 0
while t < tfin:
	it += 1
	if mod( it, 1e2) == 0:
		print t
		print x
		print '---'
		writer_traj.writerow( [t] + x )
	#
	if sum( x[1:] )<=0:
		# restart 
		print x
		print 'restarting simulation'
		out.close()
		out2.close()
		out_traj.close() 
		
		out = open( ofn1, 'w')
		writer = csv.writer( out ) # pid1, pid2, t
		out2 = open( ofn2, 'w')
		writer2 = csv.writer( out2 ) # pid, t, tiplab
		out_traj = open( ofn3, 'w')
		writer_traj = csv.writer( out_traj ) 
		
		ic = 1;  
		
		# s, i1, i2, i12, i22
		x = [ S0, 1., 0. , 0., 0.]
		
		pids1 = [ic]
		pids2 = []
		pids12 = [ic]
		pids22 = [] 
		
		t = 0.
		it = 0
		
		continue 
	#
	rates = state2rates(x)
	dt = rates2deltat(rates)
	event = rates2event(rates)
	t += dt
	if event == PROG1:
		ipid = numpy.random.choice( len(pids1) )
		pid = pids1[ipid]
		pids1.pop(ipid)
		pids2.append(pid)
		x[1]-=1
		x[2]+=1 
	elif event == PROG2: 
		ipid = numpy.random.choice( len(pids2) )
		pid = pids2[ipid]
		pids2.pop(ipid)
		writer2.writerow( [ pid , t, repr(pid)+'_chron' ] )
		x[2]-=1 
	elif event == DEATH1:
		ipid = numpy.random.choice( len(pids1) )
		pid = pids1[ipid]
		pids1.pop(ipid)
		writer2.writerow( [ pid, t, repr(pid)+'_acute' ] )
		x[1]-=1 
	elif event == DEATH2:
		ipid = numpy.random.choice( len(pids2) )
		pid = pids2[ipid]
		pids2.pop(ipid)
		writer2.writerow( [pid, t, repr(pid)+'_chron' ] )
		x[2]-=1 
		if pid==None:
			pdb.set_trace()
	if event == PROG12:## risk2
		ipid = numpy.random.choice( len(pids12) )
		pid = pids12[ipid]
		pids12.pop(ipid)
		pids22.append(pid)
		x[3]-=1
		x[4]+=1 
	elif event == PROG22:##risk2
		ipid = numpy.random.choice( len(pids22) )
		pid = pids22[ipid]
		pids22.pop(ipid)
		writer2.writerow( [ pid , t, repr(pid)+'_chron2' ] )
		x[4]-=1 
	elif event == DEATH12:##risk2
		ipid = numpy.random.choice( len(pids12) )
		pid = pids12[ipid]
		pids12.pop(ipid)
		writer2.writerow( [ pid, t, repr(pid)+'_acute2' ] )
		x[3]-=1 
	elif event == DEATH22:##risk2
		ipid = numpy.random.choice( len(pids22) )
		pid = pids22[ipid]
		pids22.pop(ipid)
		writer2.writerow( [pid, t, repr(pid)+'_chron2' ] )
		x[4]-=1 
	elif event == DEATH_S:
		x[0]-=1
	elif event == TRANSM:
		inftype = state2infectorType(pids1, pids2, pids12, pids22)
		if inftype==1:
			donor = numpy.random.choice(pids1 )
		elif inftype==2:
			donor = numpy.random.choice(pids2 )
		elif inftype == 3:
			donor = numpy.random.choice(pids12 )
		else:
			donor = numpy.random.choice(pids22 )
		#
		ic += 1
		x[0] -= 1.
		
		# assortativity 
		if inftype == 1 or inftype == 2:
			transtype = 1
		else: 
			transtype = 2
		#
		if rand() < (1.-assortprob):
			if rand() < prisk2:
				transtype = 2
			else:
				transtype = 1
			#
		#
		
		if transtype==1:
			x[1] += 1.
			pids1.append(ic )
		else:
			x[3] += 1.
			pids12.append(ic)
		#
		
		writer.writerow( [donor, ic , t ] )
	elif event == BIRTH:
		x[0] += 1.
	elif event == AGING1:
		#~ stage 1 rl2 to stage 1 rl1
		ipid = numpy.random.choice( len(pids12) )
		pid = pids12[ipid]
		pids12.pop(ipid)
		pids1.append(pid)
		x[3]-=1
		x[1]+=1 
	elif event == AGING2:
		#~ stage 1 rl2 to stage 1 rl1
		ipid = numpy.random.choice( len(pids22) )
		pid = pids22[ipid]
		pids22.pop(ipid)
		pids2.append(pid)
		x[4]-=1
		x[2]+=1 
	#
	
	
	if sum(x[1:]) == 0.:
		out.close()
		out2.close()
		exit
#

# write events for extant at tfin
for pid in pids1:
	writer2.writerow( [pid, tfin, repr(pid) + '_acute'] )
for pid in pids2:
	writer2.writerow( [pid, tfin, repr(pid) + '_chron'] )
#
for pid in pids12:
	writer2.writerow( [pid, tfin, repr(pid) + '_acute2'] )
for pid in pids22:
	writer2.writerow( [pid, tfin, repr(pid) + '_chron2'] )
#

writer_traj.writerow( [t] + x )

out.close()
out2.close()
out_traj.close() 
