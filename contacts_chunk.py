#!/usr/bin/env python3
#usage: python3 contacts_chunk.py #PDB_BASE #XTC_FILE #FILE.count_FILE #OUTPUT


#import os
import sys
#import time
#import glob
import numpy as np
#import multiprocessing
import mdtraj as md


def evaluating_contacts_chunk(pdb_file, xtc_file, pairs_indexes, r_initial, \
                              threshold=1.5, chunk=10000):
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
    return contacts

def main():

    #cores = multiprocessing.cpu_count()

    pdb_file = sys.argv[1]

    xtc_file = sys.argv[2]

    pairs_contacts = sys.argv[3]

    output_file = sys.argv[4]

    contacts = np.genfromtxt(pairs_contacts)

    # This MUST be done due to python indexes, that starts from zero
    pairs_indexes = np.subtract(contacts[:, 0:2], 1)

    r_initial = np.power(np.divide(np.multiply(contacts[:, 4], 2), contacts[:,\
                                                           3]), np.divide(1, 6))

    final_contacts = evaluating_contacts_chunk(pdb_file, xtc_file, \
                                               pairs_indexes, r_initial)

    np.savetxt(output_file, final_contacts, fmt="%d")

    return

if __name__ == "__main__":
    main()
