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
    return inputlist.reshape(-1, 1)



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
    Output:
     p_q_i - QxN numpy array with the contact probability at each given total \
     contacts value for each atom/residue.
     contacs_list - numpy 1-D (Q) array with the number of contacts.
     atoms_list - numpy 1-D (N) array with the atoms/residues involved.
    """

    cutoff = np.multiply(r_initial, threshold)

    # Correcting the numbering of atoms/residues involved.
    atoms_indexes = np.unique(pairs_indexes)
    atoms_involved = np.add(atoms_indexes, 1)

    # Correcting the contacts value and its correspondent index.
    contacts_indexes = np.arange(np.shape(pairs_indexes)[0] + 1)

    # Initializing the contacts involved. Is expected the total number of \
    # contacts is formed at least in the first frame.
    contacts_involved = np.asarray([np.shape(pairs_indexes)[0]])

    # Initializing the results array
    results = np.zeros((np.shape(contacts_indexes)[0], np.shape(r_initial)[0]))

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
        # Evaluating the contacts formed.
        contacts_involved = np.unique(np.concatenate((contacts_time, \
                                                      contacts_involved)))
        # Iterating over the number of contacts found.
        # n_frames receive number of frames found with Q contacts
        for i in contacts_indexes:
            idx = np.equal(contacts_time, i)
            n_frames[i] += idx.sum()
            results[i] += np.sum(num_below_threshold[idx], axis=0)

    # To normalize all the probabilities in each dimension after all pieces are\
    # read.
    for i in contacts_involved:
        results[i] = np.nan_to_num(np.divide(results[i], n_frames[i]))

    # Extracting nonzero results
    results_nz = results[contacts_involved]

    # Checking if all individual probabilities are normalized.
    assert np.less_equal(np.max(results_nz), 1)

    # Initiaizing the P(Q,i)
    p_q_i = np.zeros((np.shape(contacts_involved)[0], \
                        np.shape(atoms_involved)[0]))
    # The probability for each atom is given multiplying the probability of all\
    # pairs with this atom.
    for i, atom in enumerate(atoms_indexes):
        idx_atom = np.isin(pairs_indexes, atom).any(axis=1)
        p_q_i[:, i] += np.sum(results_nz[:, idx_atom], axis=1)
        # normalization over the number of contacts with each given atom/residue
        p_q_i[:, i] = np.nan_to_num(np.divide(p_q_i[:, i], idx_atom.sum()))
    # Sanity check of number of frames read
    assert np.less_equal(np.sum(n_frames[:-1]), n_frames[-1])

    return p_q_i, contacts_involved, atoms_involved


def main():

    pairs_contacts = sys.argv[4]

    contacts = np.genfromtxt(pairs_contacts)

    # This MUST be done due to python indexes, that starts from zero
    pairs_indexes = np.subtract(contacts[:, 0:2], 1)

    r_initial = evaluate_r_initial(contacts, model=str(sys.argv[1]))

    raw_prob, contacts, atoms = \
    gen_contact_probability(sys.argv[2], sys.argv[3], pairs_indexes, r_initial)

    np.savetxt("raw-" + str(sys.argv[5]), raw_prob)

    np.savetxt("Q-involved-"+ str(sys.argv[5]), contacts,
               newline="\n", header="# contacts involved", fmt="%d")

    np.savetxt("atoms-involved-"+ str(sys.argv[5]), atoms, \
               newline="\n", header="# atoms involved", fmt="%d")

    return 0

if __name__ == "__main__":
    main()
