#!/usr/bin/env python3
# USAGE: python3 phi_values.py #MODEL(AA or CA) #PDB_BASE #XTC_FILE \
# #file.cont_FILE #OUTPUT *boundaries* [should be in the order: unfolded \
# (lower and upper limit), transition state (lower and upper limit); folded \
# (lower and upper limit) in all cases.

import sys
import itertools
import numpy as np
import mdtraj as md

def evaluate_r_initial(contacts, model="AA"):
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
    if model == "CA":
        r_initial = np.power(np.divide(np.multiply(contacts[:, 4], 1.2), \
                                                   contacts[:, 3]), \
                             np.divide(1, 2))
    elif model == "AA":
        r_initial = np.power(np.divide(np.multiply(contacts[:, 4], 2), \
                                                   contacts[:, 3]),\
                             np.divide(1, 6))
    else:
        print("You have not provided an unimplemented model.")
        exit()
    return r_initial


def return_ordered(num_1, num_2):
    """Function to get two numbers and return the pair in ascending order"""
    if num_1 == num_2:
        print("Both are the same. \n Please restart with the rigth values.")
        pair = num_1, num_2
    elif num_1 < num_2:
        pair = num_1, num_2
    else:
        pair = num_2, num_1
    return pair


def test_overlap_pairs(pair_1, pair_2):
    """
    Function to test if there is an overlap between two given pairs
    Input: Two pairs.
    Output: Boolean - True if there is NO overlap.
    """
    return np.logical_or(np.less(pair_1[1], pair_2[0]), np.greater(pair_1[0], \
                                                                   pair_2[1]))


def checking_intervals_overlaping(input_array):
    """
    Function to check if it has an overlap in any pair of a given boundary \
    array.
    Input: Nx2 array with boundaries
    Output: Boolean - True if there is NO overlap in any interval.
    """
    results = []
    for i in itertools.combinations(range(np.shape(input_array)[0]), 2):
        results.append(test_overlap_pairs(input_array[i[0]], input_array[i[1]]))
    return np.all(results)


def fromlistto1d(inputlist):
    """
    Function to change the input (a list, an array or a sequence) in a column \
    vector.
    """
    inputlist = np.asarray(inputlist)
    return inputlist.reshape(-1, 1)


def phi_i(pdb_file, xtc_file, pairs_indexes, r_initial, boundaries, \
          threshold=1.5, chunk=10000):
    """
    Function to evaluate the phi-value for each atom/residue.
    Input:
     pdb_file - File with your structure (PDB or GRO files, for instance).
     xtc_file - Trajectory.
     pairs_indexes - Numpy array Nx2 with the pairs to be used to evaluate \
     the contacts. (The first two columns of the pairs section in the TPR file \
     without the header).
     r_initial - Initial distance for each given pair to be used as a reference.
     threshold - Value to be used as a threshold to evaluate the contacts.
     chunk - Size of each chunk in which the trajectory will be analyzed.
    Output
     phi - Numpy array with phi_values as a function of i - (atom/residue)
    """

    cutoff = np.multiply(threshold, r_initial)

    results = np.zeros((np.shape(boundaries)[0], np.shape(r_initial)[0]))

    # Correcting the numbering of atoms/residues involved.
    atoms_indexes = np.unique(pairs_indexes)
    atoms_involved = np.add(atoms_indexes, 1)

    # The last number of frames will store the total number for sanity check.
    n_frames = np.zeros(np.shape(boundaries)[0] + 1)
    for chunk_trajectory in md.iterload(xtc_file, top=pdb_file, chunk=chunk):
        trajectory = md.compute_distances(chunk_trajectory, pairs_indexes)
        print(chunk_trajectory)
        # Getting the number of frames of each chunk and adding to the total
        n_frames[-1] += np.shape(trajectory)[0]
        below_threshold = np.less_equal(trajectory, cutoff)
        # Generate a matrix with 1 where contacts are formed.
        num_below_threshold = np.multiply(below_threshold, 1)
        # (number of contact per timestep)
        contacts_time = np.sum(num_below_threshold, axis=1)
        for i, pair in enumerate(boundaries):
            idx = np.logical_and(np.greater_equal(contacts_time, pair[0]), \
                                 np.less_equal(contacts_time, pair[1]))
            n_frames[i] += idx.sum()
            results[i] += np.sum(num_below_threshold[idx], axis=0)

    # To normalize all the probabilities in each dimension after all pieces are\
    # read
    for i in range(np.shape(boundaries)[0]):
        results[i] = np.divide(results[i], n_frames[i])

    # Checking if all individual probabilities are normalized.
    assert np.less_equal(np.max(results), 1)

    # Initiaizing the phi-values
    pij_transition_unfolded = np.zeros(np.shape(atoms_involved))
    pij_folded_unfolded = np.zeros(np.shape(atoms_involved))
    # pij_folded = np.zeros(np.shape(atoms_involved))

    # Evaluating both parts of the phi fraction for each atom/residue.
    for i, atom in enumerate(atoms_indexes):
        idx_atom = np.isin(pairs_indexes, atom).any(axis=1)
        pij_transition_unfolded[i] += \
                    np.sum(np.subtract(results[1, idx_atom], \
                                       results[0, idx_atom]))
        pij_folded_unfolded[i] += \
                    np.sum(np.subtract(results[2, idx_atom], \
                                       results[0, idx_atom]))

    # Sanity check of number of frames read
    assert np.less_equal(np.sum(n_frames[:-1]), n_frames[-1])

    phi = np.nan_to_num(np.divide(pij_transition_unfolded, pij_folded_unfolded))

    # # This return will give phi(i) where i is the atom/residue number given \
    # # by the forcefield.
    return np.concatenate((fromlistto1d(atoms_involved), \
                           fromlistto1d(phi)), axis=1)

    # atoms_indexes = np.add(np.arange(np.shape(atoms_involved)[0]), 1)
    # # This return will give phi(i) where i is the i-th atom/residue.
    # return np.concatenate((fromlistto1d(atoms_indexes), \
    #                        fromlistto1d(phi)), axis=1)



def main():

    # pdb_file = sys.argv[1]
    #
    # xtc_file = sys.argv[2]
    #
    pairs_contacts = sys.argv[4]
    #
    # output_file = sys.argv[4]

    contacts = np.genfromtxt(pairs_contacts)

    # This MUST be done due to python indexes, that starts from zero
    pairs_indexes = np.subtract(contacts[:, 0:2], 1)

    r_initial = evaluate_r_initial(contacts, model=str(sys.argv[1]))

    # Defining the transition state
    first_lower, first_upper = return_ordered(float(sys.argv[6]), \
                                              float(sys.argv[7]))

    # Defining the unfolded state
    second_lower, second_upper = return_ordered(float(sys.argv[8]), \
                                                 float(sys.argv[9]))

    # Defining the folded state
    third_lower, third_upper = return_ordered(float(sys.argv[10]), \
                                              float(sys.argv[11]))

    # Building the boundaries array
    boundaries = np.asarray([[first_lower, first_upper], [second_lower, \
                             second_upper], [third_lower, third_upper]])

    # Checking if the intervals do not overlap each other
    assert checking_intervals_overlaping(boundaries)


    final_phi_i = phi_i(sys.argv[2], sys.argv[3], pairs_indexes, r_initial, \
                        boundaries)

    np.savetxt(sys.argv[5], final_phi_i)
    return 0

if __name__ == "__main__":
    main()
