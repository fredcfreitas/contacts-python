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
# Written by Paul Whitford.                          #
# Debugged by Ronaldo Oliveira, 15/05/10             #
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
DDT = 10  #time increment to generate temporary pdb files
GROMACSpath = '' #gromacs executable path files



##################################################################################################
# Function to convert binary trajectory file into readable temporary pieces
#
##################################################################################################
def ConvertReadable(gmxpath,filetpr,filextc,frameskip,Ti,Tf):

	runtrjconv = "echo 0 | " + gmxpath + "trjconv -b " + str(Ti) + " -e " + str(Tf) + " -nice 0 -skip " + str(frameskip) + " -s " + filetpr + " -o teste-" + str(Ti) + ".pdb -f " + filextc + " " #bash command to be runned
   	subprocess.check_output(['bash','-c', runtrjconv]) # run trjconv to every timestep
   	return
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
		#ttpr = np.genfromtxt(TRAJXTC, skip_header=3, skip_footer=2, dtype=float, usecols=(5,6,7))
		tcontfile = np.genfromtxt(CONTFILE, dtype=float)
		weigthfile = np.genfromtxt(WEIGHT, dtype=float)
	except (IOError) as errno:
		print ('I/O error. %s' % errno)
		sys.exit()
	print 'Reading a contact file'


	Rfcref = ((np.sqrt((6.0/5.0)*((tcontfile[:,4])/(tcontfile[:,3]))))*10)[np.newaxis, :].T # Distance between two contacts from file.cont
	fcref = np.concatenate((tcontfile[:,[0,1]], Rfcref), axis=1) #create array with Iaa and Jaa indices and Raa from file.cont
	#np.savetxt('fcref',fcref)
	#np.savetxt('Rfcref',Rfcref)
	Q = 0.0
	Qopt = 0.0

	contacts = []
	optcontacts = []

	to = 0 #initial time to extract trajectory
	ttf = 100 #final time to extract trajectory
	te = DDT
	while (to < ttf):
		ConvertReadable(GROMACSpath,TOPOLTPR,TRAJXTC,SKIPFRAMES,to,te)
		#analysis function
		DeleteTemporary(to)
		to = to + DDT
		te = te + DDT



	#for t in totaltime:
	# for a in range(len(fcref)):
	# 	Ci = int(fcref.item((a,0))-1) #correction due indexation of python
	# 	Cj = int(fcref.item((a,1))-1) #correction due indexation of python
	# 	Rij = np.sqrt((((ttpr.item(Ci,0))-((ttpr.item(Cj,0))))**2)+(((ttpr.item(Ci,1))-((ttpr.item(Cj,1))))**2)+(((ttpr.item(Ci,2))-((ttpr.item(Cj,2))))**2))
	# 	if Rij <= fcref.item((a,2)) * CUTOFF:
	# 		Q = Q+1 #usual contact calculation
	# 		Qopt = Qopt+1*weigthfile[a] #calculating Optimized contacts
	#contacts.append([Q]) #inside time "for loop"
	#optcontacts.append([Qopt]) #inside time "for loop"

	#print Q
	#print Qopt

	#np.savetxt('contacts.dat',contacts)
	#np.savetxt('opt-contacts.dat',optcontacts)

	#teste-"Tf".pdb #file to be deleted



	#print ttpr
	#print txtc
	#g = ttpr.item((7, 4	)) #get the item aij with 0 in the begining
	#print 'A variavel 8,5 do arquivo file.cont e ' + str(g) + '.'


#	runtrjconv = "echo 0 | trjconv -b 0.3 -e 0.5 -nice 0 -skip 1 -s run.0.tpr -o teste-.pdb -f run.0.xtc" #bash command to be runned
#	output = subprocess.check_output(['bash','-c', runtrjconv]) #bash command being runned
#	lines = ttpr.readlines() #to split ttpr in lines
#Qmax = np.int(np.max(Q))
	#Q = np.asarray([float(line.rstrip()) for line in islice(ttpr, CONNUMaa, None)])
#	print lines[100] #show 100th line
#	linhas = []
#	for item in lines:
#		linhas.append(float(lines))
#	linhas = np.array(lines) + 0.
#	print linhas[100]

if __name__ == "__main__": main()
