#!/usr/bin/env python3
#usage: python3 contacts_chunk.py #PDB_BASE #XTC_FILE #FILE.count_FILE #OUTPUT


#import os
import sys
#import time
#import glob
import numpy as np
#import multiprocessing
import mdtraj as md

THRESHOLD = 1.5
CHUNK = 1000

def main():

    #cores = multiprocessing.cpu_count()

    pdb_file = sys.argv[1]

    xtc_file = sys.argv[2]

    pairs_contacts = sys.argv[3]

    output_file = sys.argv[4]

    contacts = np.genfromtxt(pairs_contacts)

    # This MUST be done due to python indexes, that starts from zero
    contacts_list = np.subtract(contacts[:, 0:2], 1)

    r_initial = np.power(np.divide(np.multiply(contacts[:, 4], 2), contacts[:, 3]), np.divide(1, 6))

    # trajectory = md.compute_distances(md.load(xtc_file, top=pdb_file), contacts_list)
    contacts = []
    for chunk_trajectory in md.iterload(xtc_file, top=pdb_file, chunk=CHUNK):
        trajectory = md.compute_distances(chunk_trajectory, contacts_list)
        print((chunk_trajectory))
        contacts.append(np.sum(np.less_equal(trajectory, np.multiply(r_initial, THRESHOLD)), axis=1))

    contacts = np.concatenate((contacts))

    np.savetxt(output_file, contacts, fmt="%d")

    return

if __name__ == "__main__": 
    main()
