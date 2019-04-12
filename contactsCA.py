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
#import scipy as sc
#import progressbar
#from time import sleep
#from itertools import islice
#from scipy import stats
import subprocess
import multiprocessing
from functools import partial


CUTOFF = 1.2
SKIPFRAMES = 1 #1 means no skkiped frames. 2 will skip each 1 frame and so on.
t0 = 0 #initial time to extract trajectory
ttf = 1E20 #final time to extract trajectory
DDT = 10000  #time increment to generate temporary pdb files
GROMACSpath = '' #gromacs executable path files



##################################################################################################
# Function to convert binary trajectory file into readable temporary pieces
#
##################################################################################################
def ConvertReadable(gmxpath,filetpr,filextc,frameskip,Ti,Tf):

	runtrjconv = "echo 0 | " + gmxpath + "trjconv -b " + str(Ti) + " -e " + str(Tf) + " -nice 0 -skip " + str(frameskip) + " -s " + filetpr + " -o teste-" + str(Ti) + ".pdb -f " + filextc + " " #bash command to be runned
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
# Function to call contacts calculation and parallelize it
#
##################################################################################################
def CallDoContacts(cfcref,cttraj,cweigthfile,ca):
	numcores = (multiprocessing.cpu_count())
	pool = multiprocessing.Pool(processes=numcores)
	C1Contacts = partial(DoContacts, cfcref)
	C2Contacts = partial(C1Contacts, cttraj)
	C3Contacts = partial(C2Contacts, cweigthfile)
	Qvec = pool.map(C3Contacts, ca)
	pool.close()
	pool.join()
	return Qvec
#
##################################################################################################

##################################################################################################
# Function to calculate the contacts
#
##################################################################################################

def DoContacts(fcref,ttraj,weigthfile,a):
	QQ = 0
	QQopt = 0
	Ci = int(fcref.item((a,0))-1) #correction due indexation of python
	Cj = int(fcref.item((a,1))-1) #correction due indexation of python
	Rij = np.sqrt((((ttraj.item(Ci,0))-((ttraj.item(Cj,0))))**2)+(((ttraj.item(Ci,1))-((ttraj.item(Cj,1))))**2)+(((ttraj.item(Ci,2))-((ttraj.item(Cj,2))))**2))
	if Rij <= fcref.item((a,2)) * CUTOFF:
		QQ = 1 #usual contact calculation
		QQopt = 1*weigthfile[a] #calculating Optimized contacts
	return QQ,QQopt
#
##################################################################################################


def main():

	if len(sys.argv) > 4: ## To open just if exist  file in argument
		TOPOLTPR = sys.argv[1]
		TRAJXTC = sys.argv[2]
		CONTFILE = sys.argv[3]
		WEIGHT = sys.argv[4]
	elif (len(sys.argv) == 4):
		TOPOLTPR = sys.argv[1]
		TRAJXTC = sys.argv[2]
		CONTFILE = sys.argv[3]
	else:
		print ('One (or more) input file(s) is(are) missing. Please insert files using (at least): ./contacts_XX.py file.TPR file.XTC file.cont')
		sys.exit()
	try:
		tcontfile = np.genfromtxt(CONTFILE, dtype=float)
		if len(sys.argv) > 4:
			Aweigthfile = np.genfromtxt(WEIGHT, dtype=float)
		if len(sys.argv) == 4:
			Aweigthfile = np.ones(len(tcontfile))
			print ('Without weigth file.')
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
	numa = 0 #variable to count number of atoms
	flagnuma = True

	while (to < ttf):

		try:
			ConvertReadable(GROMACSpath,TOPOLTPR,TRAJXTC,SKIPFRAMES,to,te)
		except (IOError) as errnoa:
			print ('I/O error. %s' % errnoa)
			sys.exit()
		except ValueError:
			print 'There is something wrong.'
		except:
			print 'You made a bad choice for initial (or final) time. But there is no problem.'
			to=to-2*DDT
			te = -1 #last frame to read from trajectory
			ConvertReadable(GROMACSpath,TOPOLTPR,TRAJXTC,SKIPFRAMES,to,te)
		try:
			postemp = open('teste-' + str(to) + '.pdb', 'r') #open temporary Translated trajectory file
			tempfile = postemp.readlines() #split it in a list of lines
		except (IOError) as errnoa:
			print ('I/O error. %s' % errnoa)
			sys.exit()

		for line in tempfile:

			if to==t0: #node to calculate number of atoms
				if ('ATOM' in line) and (flagnuma):
					numa +=1
				elif 'TER' in line:
					flagnuma = False

			if 't=' in line:
				setimes = float(line[26:1000])
			if setimes not in utimes:
				repos = True #open to read positions
				if ('ATOM' in line) and (repos):
					X.append(float(line[30:37]))
					Y.append(float(line[38:45]))
					Z.append(float(line[46:53]))
				elif 'TER' in line:
					repos = False #close to read positions
					utimes.append(setimes)
					Attraj = np.transpose(np.array([X,Y,Z])) #reconstruced array with positions of all atoms in each time
					aa = range(len(Afcref))
					BQQopt = CallDoContacts(Afcref,Attraj,Aweigthfile,aa)
					TQQopt = np.sum(BQQopt, axis=0)
					Q = TQQopt[0]
					Qopt = TQQopt[1]
					contacts.append([Q]) #inside time "for loop"
					optcontacts.append([Qopt]) #inside time "for loop"
					X = [] #position vector of each time
					Y = []
					Z = []
		DeleteTemporary(to)
		to = to + DDT
		te = te + DDT
		np.savetxt('contacts.dat',contacts)
		np.savetxt('opt-contacts.dat',optcontacts)
#	print numa #to print the number of elements analyzed.


if __name__ == "__main__": main()
