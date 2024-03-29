#| ### This code calculate the contact probability using mdtraj and numpy libraries.
#|
#!/usr/bin/env python3
# usage: python3 contact_probability.py #MODEL(AA or CA) #PDB_BASE #XTC_FILE \
# #file.cont_FILE #OUTPUT
#
# #file.cont_FILE is the pairs section extracted from SMOG *.top file

import sys
import numpy as np
import mdtraj as md

#------------------------------------------------------------------------------


def evaluate_r_initial(contacts, model="AA"):
    """
    Function to evaluate initial pairwise distances found in the initial model.
    Input:
     contacts_list - Array with the definitions from pairs section of the \
     forcefield (SMOG TOP file.)
     model - AA = All-Atom; CA = Carbon_alpha Coarse-Grained
    Output:
     r_initial - vector with the initial distances of each pair.
    """
    if model == "CA":
        r_initial = np.power(
            np.divide(np.multiply(contacts[:, 4], 1.2), contacts[:, 3]), \
                np.divide(1, 2)
        )
    elif model == "AA":
        r_initial = np.power(
            np.divide(np.multiply(contacts[:, 4], 2), contacts[:, 3]), \
                np.divide(1, 6)
        )
    else:
        print("You have not provided a recognizable model.")
        sys.exit()
    return r_initial


def fromlistto1d(inputlist):
    """
    Function to change the input (a list, an array or a sequence) into a \
    column vector.
    """
    inputlist = np.asarray(inputlist)
    return inputlist.reshape(-1, 1)


def gen_contact_probability(
    pdb_file, xtc_file, pairs_indexes, r_initial, threshold=1.5, chunk=10000\
    ):
    """
    Function to evaluate the contact probability per atom/residue.
    Input:
     pdb_file - File with your structure (PDB or GRO file).
     xtc_file - Trajectory.
     pairs_indexes - Numpy array Nx2 with the pairs to be used to evaluate \
      the contacts. (The first two columns of the pairs section in the TOP \
      file without the header and corrected to python indexation).
     r_initial - Initial distance for each given pair.
     threshold - Value to be used as a threshold to evaluate the contacts.
     chunk - Size of each chunk in which the trajectory will be analyzed.
    Output:
     p_q_i - QxN numpy array with the contact probability at each given total \
      contacts value for each atom/residue.
     contacts_list - numpy 1-D (Q) array with the number of contacts.
     atoms_list - numpy 1-D (N) array with the atoms/residues involved.
    """
    # Calculating the cutoff distances for each pair
    cutoff = np.multiply(r_initial, threshold)
    # Correcting the  atoms/residue numbering.
    atoms_indexes = np.unique(pairs_indexes)
    atoms_involved = np.add(atoms_indexes, 1)
    # Correcting the contact index to start from 1.
    contacts_indexes = np.arange(np.shape(pairs_indexes)[0] + 1)
    # Initializing the contacts involved. Is expected the total number of \
    # contacts is formed at least in the first frame.
    contacts_involved = np.asarray([np.shape(pairs_indexes)[0]])
    # Initializing the results array
    results = np.zeros((np.shape(contacts_indexes)[0], np.shape(r_initial)[0]))
    # n_frames will store the total number of frames read as a sanity check.
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
        # Adding the number of contacts evaluated to the list of previous \
        # number of formed contacts list.
        contacts_involved = np.unique(
            np.concatenate((contacts_time, contacts_involved))
        )
        # Iterating over the number of contacts found.
        # n_frames receive number of frames found with Q contacts
        for i in contacts_indexes:
            idx = np.equal(contacts_time, i)
            n_frames[i] += idx.sum()
            results[i] += np.sum(num_below_threshold[idx], axis=0)
    # To normalize all the probabilities in each dimension after all pieces \
    # are read.
    for i in contacts_involved:
        results[i] = np.nan_to_num(np.divide(results[i], n_frames[i]))
    # Extracting nonzero results
    results_nz = results[contacts_involved]
    # Checking if all individual probabilities are normalized.
    assert np.less_equal(np.max(results_nz), 1)
    # Initializing the P(Q,i)
    p_q_i = np.zeros((np.shape(contacts_involved)[0], \
                      np.shape(atoms_involved)[0]))
    # The probability for each atom is given multiplying the probability of \
    # all pairs with this atom.
    for i, atom in enumerate(atoms_indexes):
        idx_atom = np.isin(pairs_indexes, atom).any(axis=1)
        p_q_i[:, i] += np.sum(results_nz[:, idx_atom], axis=1)
        # normalization over the number of contacts with each given \
        # atom/residue
        p_q_i[:, i] = np.nan_to_num(np.divide(p_q_i[:, i], idx_atom.sum()))
    # Sanity check of number of frames read
    assert np.less_equal(np.sum(n_frames[:-1]), n_frames[-1])
    return p_q_i, contacts_involved, atoms_involved


#------------------------------------------------------------------------------

#| #### Here starts the main part. Be sure to use the right files.


def main():
    """
    Starting the main function. Take a look at the comments on top regarding \
    input files.
    """
    pairs_contacts = sys.argv[4]
    # comments are inserted to avoid problems with common file header
    contacts = np.loadtxt(pairs_contacts, comments=['#', '[', ';'])
    # This MUST be done due to python indexes, that starts from zero
    pairs_indexes = np.subtract(contacts[:, 0:2], 1)
    # Getting distances from TOP file (from original structure)
    r_initial = evaluate_r_initial(contacts, model=str(sys.argv[1]))
    # Calling the contact probability function
    raw_prob, contacts, atoms = gen_contact_probability(
        sys.argv[2], sys.argv[3], pairs_indexes, r_initial
    )
    # Saving the contact probability array
    np.savetxt("raw-" + str(sys.argv[5]), raw_prob)
    # Saving the list of contacts.
    np.savetxt(
        "Q-involved-" + str(sys.argv[5]),
        contacts,
        newline="\n",
        header="# contacts involved",
        fmt="%d"
    )
    # Saving the list of atoms/residues involved
    np.savetxt(
        "atoms-involved-" + str(sys.argv[5]),
        atoms,
        newline="\n",
        header="# atoms involved",
        fmt="%d"
    )
    return 0

#------------------------------------------------------------------------------
#| Uncomment the first lines below with the 'sys.argv' to run the notebook
#| instead of the script.
#| Do not forget to properly replace the options. Example to run the notebook
#| with shared files:

#| `sys.argv = "contact_probability.py AA ../share/ci2-AA-120-run.gro ../share/ci2-AA-120-run.xtc ../share/ci2-AA-contacts.dat output".split()`
#sys.argv = "contact_probability.py #MODEL(AA or CA) #PDB_BASE #XTC_FILE #file.cont_FILE #OUTPUT".split()
if __name__ == "__main__":
    main()
