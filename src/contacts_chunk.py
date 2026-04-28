#!/usr/bin/env python3
#
################################################################################
# Script to calculate contacts from *.xtc gromacs files using mdtraj library   #
# To better results, use the topology file given by SMOG or SMOG2 (*.gro)      #
################################################################################
#
#USAGE: python3 contacts_chunk.py #MODEL (CA or AA) #PDB_BASE #XTC_FILE \
# #file.cont_FILE #OUTPUT_filename
# file.cont is the GROMACS forcefield pairs section (SMOG-like)

import sys
import numpy as np
import mdtraj as md


# Threshold to be used to consider a native contact when comparing the distances
THRESHOLD = 1.5

# Number of frames to be read at a time
CHUNK = 10000

# Skip every STRIDE frames when analyzing
STRIDE = 1

def evaluate_r_initial(contacts, model="AA", precision=np.double):
    """
    Function to evaluate initial pairwise distances accordingly the model \
    simulated.
    Input:
     contacts_list - Array with the definitions from pairs section of the \
     forcefield (TPR file.)
     model - AA = All-Atom; CA = Carbon_alpha Coarse-Grained
    Output:
     r_initial - vector with the initial distance of each pair.
    """
    contacts = np.asarray(contacts, dtype=precision)
    if model == "CA":
        r_initial = np.power(np.divide(np.multiply(contacts[:, -1], 1.2), \
                                                   contacts[:, -2]), \
                             np.divide(1, 2))
    elif model == "AA":
        r_initial = np.power(np.divide(np.multiply(contacts[:, -1], 2), \
                                                   contacts[:, -2]),\
                             np.divide(1, 6))
    else:
        print("You have not provided a model.")
        sys.exit()
    return r_initial


def evaluating_contacts_chunk(pdb_file, xtc_file, pairs_indexes, r_initial, \
                              threshold=1.5, chunk=10000, stride=1):
    """
    Function to evaluate the number of contacts for each given timestep.
    Input:
     pdb_file - File with your structure (PDB or GRO files for instance).
     xtc_file - Trajectory.
     pairs_indexes - Numpy array Nx2 with the pairs to be used to evaluate \
     the contacts. (The first two columns of the pairs section in the TPR file \
     without the header).
     r_initial - Initial distance for each given pair to be used as a reference.
     threshold - Value to be used as a threshold to evaluate the contacts.
     chunk - Size of each chunk in which the trajectory will be analyzed.
    Output: Nx1 numpy array with the total number of contacts for each \
     timestep.
    """
    contacts = []
    for chunk_trajectory in md.iterload(xtc_file, top=pdb_file, chunk=chunk):
        trajectory = md.compute_distances(chunk_trajectory, pairs_indexes)
        print((chunk_trajectory))
        contacts.append(np.sum(np.less_equal(trajectory, np.multiply(r_initial,\
                                                           threshold)), axis=1))

    contacts = np.concatenate((contacts))
    return contacts[::stride]

def main():

    #cores = multiprocessing.cpu_count()

    pdb_file = sys.argv[2]

    xtc_file = sys.argv[3]

    pairs_contacts = sys.argv[4]

    output_file = sys.argv[5]

    contacts = np.genfromtxt(pairs_contacts)

    # This MUST be done due to python indexes, that starts from zero
    pairs_indexes = np.subtract(contacts[:, 0:2], 1)

    r_initial = evaluate_r_initial(contacts, model=str(sys.argv[1]))

    final_contacts = evaluating_contacts_chunk(pdb_file, xtc_file, \
                                               pairs_indexes, r_initial,\
                                               threshold=THRESHOLD, \
                                               chunk=CHUNK, stride=STRIDE)

    np.savetxt(output_file, final_contacts, fmt="%d")

    return 0

if __name__ == "__main__":
    main()
