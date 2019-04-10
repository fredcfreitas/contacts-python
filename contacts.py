#!/usr/bin/python2
#
######################################################
# Contacts.py is a simple scritp written to analyze  #
# simulations run in GROMACS. You must supply a      #
# contact file, with 10-12 parameters, the .tpr file #
# used for the simulation, and an trajectory (.xtc)  #
# file. You must also have Gromacs (4.x) installed   #
# on your machine. This is a straightforward script  #
# you can modify in any way you see fit. You must    #
# observe GNU license to use it.                     #
# Written by Paul Whitford, 11/02/2009.              #
# Debugged by Ronaldo Oliveira, 05/15/10             #
# Translated to python by Frederico Campos Freitas   #
######################################################



# CONTFILE is the file that defines the contacts.  Specific formatting must be
# followed: Copy the "pairs" terms from your C-Alpha Structure-based topology file.
# remove the "1" and reformat each line so it is space delimited.  If you have
# a blank line in the contact file, the program will probably crash.
# example formatting can be found at http://sbm.ucsd.edu/contact_ex

import sys
import numpy as np
import scipy as sc
#import progressbar
#from time import sleep
from itertools import islice
from scipy import stats
import subprocess

CUTOFF = 1.2
SKIPFRAMES = 1 #1 means no skkiped frames. 2 will skip each 1 frame and so on.
t0 = 99984  #initial time to extract trajectory
ttf = 100100 #1E20 #final time to extract trajectory
DDT = 10000 #10000  #time increment to generate temporary pdb files
GROMACSpath = '' #gromacs executable path files



##################################################################################################
# Function to convert binary trajectory file into readable temporary pieces
#
##################################################################################################
def ConvertReadable(gmxpath,filetpr,filextc,frameskip,Ti,Tf):

	runtrjconv = "echo 0 | " + gmxpath + "trjconv -b " + str(Ti) + " -e " + str(Tf) + " -nice 0 -skip " + str(frameskip) + " -s " + filetpr + " -o teste-" + str(Ti) + ".pdb -f " + filextc + " " #bash command to be runned
	print runtrjconv
	tstatus = subprocess.check_output(['bash','-c', runtrjconv]) # run trjconv to every timestep
	return tstatus
##################################################################################################

##################################################################################################
# Function to delete converted trajectory temporary files
#
##################################################################################################
def DeleteTemporary(Ti):

	deletetemp = "rm teste-" + str(Ti) + ".pdb" #bash command to be runned
   	subprocess.check_output(['bash','-c', deletetemp]) # run trjconv to every timestep
   	return
##################################################################################################


##################################################################################################
# Function to calculate the contacts
#
##################################################################################################
def DoContacts(fcref,ttraj,weigthfile):
	Q = 0
	Qopt = 0
	for a in range(len(fcref)):
		Ci = int(fcref.item((a,0))-1) #correction due indexation of python
		Cj = int(fcref.item((a,1))-1) #correction due indexation of python
		Rij = np.sqrt((((ttraj.item(Ci,0))-((ttraj.item(Cj,0))))**2)+(((ttraj.item(Ci,1))-((ttraj.item(Cj,1))))**2)+(((ttraj.item(Ci,2))-((ttraj.item(Cj,2))))**2))
		if Rij <= fcref.item((a,2)) * CUTOFF:
			Q = Q+1 #usual contact calculation
			Qopt = Qopt+1*weigthfile[a] #calculating Optimized contacts
	return Q,Qopt
##################################################################################################

def main():

	if len(sys.argv) > 2: ## To open just if exist  file in argument
		TOPOLTPR = sys.argv[1]
		TRAJXTC = sys.argv[2]
		CONTFILE = sys.argv[3]
		WEIGHT = sys.argv[4]
	else:
		print ('One (or more) input file(s) is(are) missing. Please insert files using contacts_XX.py file.TPR file.XTC file.cont')
		sys.exit()
	try:
		#txtc = np.genfromtxt(TOPOLTPR, dtype=float) #get values of TOPOLTPR from numerical array.
		#ttraj = np.genfromtxt(TRAJXTC, skip_header=3, skip_footer=2, dtype=float, usecols=(5,6,7))
		tcontfile = np.genfromtxt(CONTFILE, dtype=float)
		Aweigthfile = np.genfromtxt(WEIGHT, dtype=float)
	except (IOError) as errno:
		print ('I/O error. %s' % errno)
		sys.exit()
	print 'Reading a contact file'


	Rfcref = ((np.sqrt((6.0/5.0)*((tcontfile[:,4])/(tcontfile[:,3]))))*10)[np.newaxis, :].T # Distance between two contacts from file.cont
	Afcref = np.concatenate((tcontfile[:,[0,1]], Rfcref), axis=1) #create array with Iaa and Jaa indices and Raa from file.cont

	contacts = []
	optcontacts = []

	to = t0
	te = t0 + DDT

	utimes = [] #timesteps already read in simulation
	setimes = ''
	X = []
	Y = []
	Z = []

	while (to < ttf):
		try:
			ttstatus = ConvertReadable(GROMACSpath,TOPOLTPR,TRAJXTC,SKIPFRAMES,to,te)
			print 'Olhe'
			print ttstatus
		except (IOError) as errnoa:
			print ('I/O error. %s' % errnoa)
			sys.exit()
		except ValueError:
			print 'Value'
			#ConvertReadable(GROMACSpath,TOPOLTPR,TRAJXTC,SKIPFRAMES,to=to-1.5*DDT,te='')
		except:
			print 'Undefined'
			print ttstatus
			#teste = ttstatus[110][0]
			to=to-2*DDT
			te = -1 #last frame to read from trajectory
			ttstatus = ConvertReadable(GROMACSpath,TOPOLTPR,TRAJXTC,SKIPFRAMES,to,te)
		try:
			postemp = open('teste-' + str(to) + '.pdb', 'r') #open temporary Translated trajectory file
			tempfile = postemp.readlines() #split it in a list of lines
		except (IOError) as errnoa:
			print ('Ehhh I/O error. %s' % errnoa)
			sys.exit()

		for line in tempfile:
			if 't=' in line:
				setimes = float(line[26:1000])
				#print setimes
			if setimes not in utimes:
				repos = True
				if ('ATOM' in line) and (repos):
					X.append(float(line[30:37]))
					Y.append(float(line[38:45]))
					Z.append(float(line[46:53]))
				elif 'TER' in line:
					utimes.append(setimes)
					Attraj = np.transpose(np.array([X,Y,Z])) #reconstruced array with positions of all atoms in each time
					Q,Qopt = DoContacts(Afcref,Attraj,Aweigthfile)
					contacts.append([Q]) #inside time "for loop"
					optcontacts.append([Qopt]) #inside time "for loop"
					#print Attraj
					#np.savetxt('setimes' + str(setimes) +'.dat',Attraj)
					X = [] #position vector of each time
					Y = []
					Z = []
					repos = False
		DeleteTemporary(to)
		to = to + DDT
		te = te + DDT
		np.savetxt('contacts.dat',contacts)
		np.savetxt('opt-contacts.dat',optcontacts)


if __name__ == "__main__": main()
