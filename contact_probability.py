#!/usr/bin/env python3
# usage: python3 contact_probability.py #MODEL(AA or CA) #PDB_BASE #XTC_FILE \
# #file.cont_FILE #OUTPUT

import sys
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


def fromlistto1d(inputlist):
    """
    Function to change the input (a list, an array or a sequence) in a column \
    vector.
    """
    inputlist = np.asarray(inputlist)
    return np.reshape(inputlist, (np.shape(inputlist)[0], 1))


def gen_contact_probability(pdb_file, xtc_file, pairs_indexes, r_initial, \
                            threshold=1.5, chunk=10000):
    """
    Function to evaluate the contact probability per residue.
    Input:
     pdb_file - File with your structure (PDB or GRO files for instance).
     xtc_file - Trajectory.
     pairs_indexes - Numpy array Nx2 with the pairs to be used to evaluate \
     the contacts. (The first two columns of the pairs section in the TPR file \
     without the header).
     r_initial - Initial distance for each given pair to be used as a reference.
     threshold - Value to be used as a threshold to evaluate the contacts.
     chunk - Size of each chunk in which the trajectory will be analyzed.
    Output: Nx3 numpy array with the contact probability for each atom/residue \
    at each given total contacts value.
    """

    cutoff = np.multiply(threshold, r_initial)

    results = np.zeros((np.shape(pairs_indexes)[0], np.shape(r_initial)[0]))

    # Correcting the numbering of atoms/residues involved.
    atoms_indexes = np.unique(pairs_indexes)
    atoms_involved = np.add(atoms_indexes, 1)

    # Correcting the contacts value and its correspondent index.
    contacts_indexes = np.arange(np.shape(pairs_indexes)[0])
    contacts_involved = np.add(contacts_indexes, 1)

    # The last number of frames will store the total number for sanity check.
    n_frames = np.zeros(np.shape(contacts_indexes)[0] + 1)
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
        # Iterating over the number of contacts found. i-index, j-actual
        for i, j in zip(contacts_indexes, contacts_involved):
            idx = np.equal(contacts_time, j)
            n_frames[i] += idx.sum()
            results[i] += np.sum(num_below_threshold[idx], axis=0)

    # To normalize all the probabilities in each dimension after all pieces are\
    # read.
    for i in contacts_indexes:
        results[i] = np.nan_to_num(np.divide(results[i], n_frames[i]))

    # Checking if all individual probabilities are normalized.
    assert np.less_equal(np.max(results), 1)

    # Initiaizing the P(Q,i)
    p_q_i = np.zeros((np.shape(pairs_indexes)[0], \
                        np.shape(atoms_involved)[0]))
    # The probability for each atom is given multiplying the probability of all\
    # pairs with this atom.
    for i, atom in enumerate(atoms_indexes):
        idx_atom = np.isin(pairs_indexes, atom).any(axis=1)
        p_q_i[:, i] += np.sum(results[:, idx_atom], axis=1)
        # normalization over the number of contacts with each given atom/residue
        p_q_i[:, i] = np.divide(p_q_i[:, i], idx_atom.sum())
    # Sanity check of number of frames read
    assert np.less_equal(np.sum(n_frames[:-1]), n_frames[-1])

    # Formatted array. First column and row are the contacts and atoms involved\
    # respectively.
    p_q_i_formatted = \
    np.hstack((fromlistto1d(np.hstack((np.asarray([0]), contacts_involved))), \
               np.vstack((atoms_involved, p_q_i))))

    return p_q_i, p_q_i_formatted


def main():

    pairs_contacts = sys.argv[4]

    contacts = np.genfromtxt(pairs_contacts)

    # This MUST be done due to python indexes, that starts from zero
    pairs_indexes = np.subtract(contacts[:, 0:2], 1)

    r_initial = evaluate_r_initial(contacts, model=str(sys.argv[1]))

    raw_prob, formatted_prob = \
    gen_contact_probability(sys.argv[2], sys.argv[3], pairs_indexes, r_initial)

    np.savetxt("raw-" + str(sys.argv[5]), raw_prob)
    np.savetxt("formatted-" + str(sys.argv[5]), formatted_prob)
    return 0

if __name__ == "__main__":
    main()
