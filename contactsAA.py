#!/usr/bin/env python3
#
###################################################################################################
# Contacts.py is a simple scritp written to analyze simulations run in GROMACS. You must supply   #
# a contact file, with 10-12 parameters, the .tpr file used for the simulation, and an trajectory #
# (.xtc) file. You must also have Gromacs (4.x) installed on your machine. This is a straightfor- #
# ward script you can modify in any way you see fit. You must observe GNU license to use it.      #
# Written by Paul Whitford, 11/02/2009.				                       		                  #
# Debugged by Ronaldo Oliveira, 05/15/10                                                          #
# Translated to python by Frederico Campos Freitas 03/06/2019			                          #
###################################################################################################



# CONTFILE is the file that defines the contacts.  Specific formatting must be
# followed: Copy the "pairs" terms from your C-Alpha Structure-based topology file.
# remove the "1" and reformat each line so it is space delimited.  If you have
# a blank line in the contact file, the program will probably crash.
# example formatting can be found at http://sbm.ucsd.edu/contact_ex
