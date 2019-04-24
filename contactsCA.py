#!/usr/bin/env python3
#
################################################################################
# Contacts.py is a simple scritp written to analyze simulations run in GROMACS.#
# You must supply a contact file, with 10-12 parameters, the .tpr file used for#
# the simulation, and an trajectory (.xtc) file. You must also have Gromacs	   #
# (4.x) installed on your machine. This is a straightforward script you can	   #
# modify in any way you see fit. You must observe GNU license to use it.       #
# Written by Paul Whitford, 11/02/2009.				                           #
# Debugged by Ronaldo Oliveira, 05/15/10                                       #
# Translated to python by Frederico Campos Freitas 03/06/2018   			   #
################################################################################



# CONTFILE is the file that defines the contacts. Specific formatting must be
# followed: Copy the "pairs" terms from your C-Alpha Structure-based topology
# file. Remove the "1" and reformat each line so it is space delimited.
# If you have a blank line in the contact file, the program will probably crash.
# Example formatting can be found at http://sbm.ucsd.edu/contact_ex

import sys
import subprocess
import time
import multiprocessing
import multiprocessing.pool
import concurrent.futures
#from functools import partial
from bisect import bisect_left
import numpy as np
#import scipy as sc
#import progressbar
#from time import sleep
#from itertools import islice
#from scipy import stats



CUTOFF = 1.2
SKIPFRAMES = 1 #1 means no skkiped frames. 2 will skip each 1 frame and so on.
t0 = 0 #initial time to extract trajectory
ttf = 1E20 #final time to extract trajectory
DDT = 4000 #20000  #time increment to generate temporary pdb files
GROMACSpath = '' #gromacs executable path files

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
	Process = NoDaemonProcess


################################################################################
# Function to convert binary trajectory file into readable temporary pieces	   #
#																			   #
################################################################################
def ConvertReadable(gmxpath,filetpr,filextc,frameskip,Ti,Tf):
	runtrjconv = "echo 0 | " + gmxpath + "trjconv -b " + str(Ti) + " -e " + str(Tf) + " -nice 0 -skip " + str(frameskip) + " -s " + filetpr + " -o teste-" + str(Ti) + ".pdb -f " + filextc + " " #bash command to be runned
	tstatus = subprocess.Popen(['bash','-c', runtrjconv], stdout=subprocess.PIPE).communicate() # run trjconv to every timestep
	return tstatus
################################################################################

################################################################################
# Function to delete converted trajectory temporary files										                     #
#												                                                                         #
################################################################################
def DeleteTemporary(Ti):
	deletetemp = "rm teste-" + str(Ti) + ".pdb" #bash command to be runned
	tstatus = subprocess.Popen(['bash','-c', deletetemp], stdout=subprocess.PIPE).communicate() # run trjconv to every timestep
	return tstatus
################################################################################


################################################################################
# Function to call contacts calculation and parallelize it					   #
#												 							   #
################################################################################
def CallDoContacts(cfcref,cttraj,cweigthfile,ca):
	pool = multiprocessing.Pool(multiprocessing.cpu_count())
	#C1Contacts = partial(DoContacts, cfcref)
	#C2Contacts = partial(C1Contacts, cttraj)
	#C3Contacts = partial(C2Contacts, cweigthfile)
	#Qvec = pool.map(C3Contacts, ca)
	Qvec = [pool.apply_async(DoContacts, args=(cfcref, cttraj, cweigthfile, a)) for a in ca]
	pool.close()
	#pool.join()
	return Qvec

################################################################################

################################################################################
# Function to calculate the contacts								 		   #
#							 												   #
################################################################################

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

################################################################################

################################################################################
# Function to get atom coordinates from pdb file and calculate contacts		   #
#							 												   #
################################################################################
def read_and_calculate(setimes, utimes, contacts, optcontacts, Afcref, Aweigthfile, X, Y, Z, tempfile):

	end5 = time.time()
	with concurrent.futures.ProcessPoolExecutor() as executor:
		for line in tempfile:
			if 't=' in line:
				setimes = float(line[27:100])
			#if setimes not in utimes:
			if (bisect_left(utimes, setimes) >= len(utimes)):
				repos = True #open to read positions
				if ('ATOM' in line) and (repos):
					X = np.append(X, (line[30:37]))
					Y = np.append(Y, (line[38:45]))
					Z = np.append(Z, (line[46:53]))
				elif 'TER' in line:
					end6 = time.time()
					repos = False #close to read positions
					utimes = np.append(utimes, setimes)
					Attraj = np.transpose(np.array([X,Y,Z], dtype=float)) #reconstruced array with positions of all atoms in each time
					aa = list(range(len(Afcref)))
					BQQopt = CallDoContacts(Afcref, Attraj, Aweigthfile, aa)
					TQQopt = np.sum(BQQopt, axis=0)
					Q = TQQopt[0]
					Qopt = TQQopt[1]
					contacts = np.append(contacts, Q) #inside time "for loop"
					optcontacts = np.append(optcontacts, Qopt) #inside time "for loop"
					X = np.array([]) #position vector of each time
					Y = np.array([])
					Z = np.array([])

	end7 = time.time()

	return contacts, optcontacts

################################################################################




def main():

	start1 = time.time()

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
		print('One (or more) input file(s) is(are) missing. Please insert files using (at least): ./contacts_XX.py file.TPR file.XTC file.cont [weigthfile]')
		sys.exit()
	try:
		tcontfile = np.genfromtxt(CONTFILE, dtype=float)
		if len(sys.argv) > 4:
			Aweigthfile = np.genfromtxt(WEIGHT, dtype=float)
		if len(sys.argv) == 4:
			Aweigthfile = np.ones(len(tcontfile))
			print('Without weigth file.')
	except (IOError) as errno:
		print(('I/O error. %s' %errno))
		sys.exit()


	print('Reading a contact file')

#	if '.tpr' in TOPOLTPR:
#		print('TPR file')
#	elif: '.ndx' in TOPOLTPR:
#		print('ndx file')
#	else:
#		'I did not understand you'

	end1 = time.time()


	Rfcref = ((np.sqrt((6.0/5.0)*((tcontfile[:,4])/(tcontfile[:,3]))))*10)[np.newaxis, :].T # Distance between two contacts from file.cont
	Afcref = np.concatenate((tcontfile[:,[0,1]], Rfcref), axis=1) #create array with Iaa and Jaa indices and Raa from file.cont



	contacts = np.array([])
	optcontacts = np.array([])

	to = t0
	te = t0 + DDT

	utimes = [] #timesteps already read in simulation
	setimes = ''
	X = np.array([])
	Y = np.array([])
	Z = np.array([])



	end2 = time.time()

	while (to < ttf):

		end3 = time.time()

		try:
			ConvertReadable(GROMACSpath,TOPOLTPR,TRAJXTC,SKIPFRAMES,to,te)
		except (IOError) as errnoa:
			print(('I/O error. %s' % errnoa))
			sys.exit()
		except ValueError:
			print('There is something wrong.')
		except:
			print('You made a bad choice for initial (or final) time. But there is no problem.')
			to=to-2*DDT
			te = -1 #last frame to read from trajectory
			ConvertReadable(GROMACSpath,TOPOLTPR,TRAJXTC,SKIPFRAMES,to,te)
		try:
			postemp = open('teste-' + str(to) + '.pdb', 'r') #open temporary Translated trajectory file
			tempfile = postemp.readlines() #split it in a list of lines
		except (IOError) as errnoa:
			print(('I/O error. %s' % errnoa))
			sys.exit()

		end4 = time.time()


		contacts, optcontacts = read_and_calculate(setimes, utimes, contacts, optcontacts, Afcref, Aweigthfile, X, Y, Z, tempfile)
		#pool = MyPool(multiprocessing.cpu_count())

		#[contacts, optcontacts] = [pool.apply(read_and_calculate, args=(setimes, utimes, contacts, optcontacts, Afcref, Aweigthfile, X, Y, Z, line)) for line in tempfile]

		#pool.close()
		#pool.join()

		DeleteTemporary(to)
		to = to + DDT
		te = te + DDT

		#print(start1)
		#print(end1)
		#print(end2)
		#print(end3)
		#print(end4)
		#print(end5)
		#print(end6)
		#print(end7)

		np.savetxt('contacts.dat', contacts, fmt='%d')
		np.savetxt('opt-contacts.dat', optcontacts, fmt='%10.3f')


if __name__ == "__main__": main()
