#!/usr/bin/env python3
#usage: python3 contact_probability.py #

import sys
import numpy as np
import mdtraj as md

def fromlistto1d(inputlist):
    """
    Function to change the input (a list, an array or a sequence) in a column \
    vector.
    """
    inputlist = np.asarray(inputlist)
    return np.reshape(inputlist, (np.shape(inputlist)[0], 1))

def gen_contact_probability(pdb_file, xtc_file, pairs_contacts, threshold=1.5,\
                            chunk=10000):
    """
    Function to evaluate the contact probability per residue.
    Input:
     pdb_file - File with your structure (PDB or GRO files for instance).
     xtc_file - Trajectory.
     pairs_contacts - File with the pairs section extracted from the \
     forcefield. (The pairs section of the TPR file without the header).
     threshold - Value to be used as a threshold to evaluate the contacts.
     chunk - Size of each chunk in which the trajectory will be analyzed.
    Output: Nx2 numpy array with the contact probability for each pair.
    """
    contacts = np.genfromtxt(pairs_contacts)

    # This MUST be done due to python indexes, that starts from zero
    pairs_indexes = np.subtract(contacts[:, 0:2], 1)

    r_initial = np.power(np.divide(np.multiply(contacts[:, 4], 2), contacts[:, \
                                           3]), np.divide(1, 6))

    cutoff = np.multiply(threshold, r_initial)

    probability = np.zeros(np.shape(r_initial))

    n_frames = 0
    for chunk_trajectory in md.iterload(xtc_file, top=pdb_file, chunk=chunk):
        trajectory = md.compute_distances(chunk_trajectory, pairs_indexes)
        print((chunk_trajectory))
        # Getting the number of frames of each chunk and adding to the total
        n_frames += np.shape(trajectory)[0]
        # for idx in range(np.shape(pairs_indexes)[0]):
        for idj in range(np.shape(trajectory)[0]):
            below_threshold = np.less_equal(trajectory[idj], cutoff)
            probability += np.multiply(below_threshold, 1)

    final_probability = np.divide(probability, n_frames)

    # To certify it has no zero in the final_probability
    assert not np.any(np.less(final_probability, 0))

    residues = np.add(np.arange(np.shape(r_initial)[0]), 1)

    return np.concatenate((fromlistto1d(residues), \
                           fromlistto1d(final_probability)), axis=1)

def main():

    # pdb_file = sys.argv[1]
    #
    # xtc_file = sys.argv[2]
    #
    # pairs_contacts = sys.argv[3]
    #
    # output_file = sys.argv[4]

    final_contacs_prob = gen_contact_probability(sys.argv[1], sys.argv[2], \
                                                 sys.argv[3])

    np.savetxt(sys.argv[4], final_contacs_prob)
    return

if __name__ == "__main__":
    main()
