#!/usr/bin/env python3
# usage: python3 contact_probability.py #MODEL(AA or CA) #PDB_BASE #XTC_FILE \
# #file.cont_FILE #OUTPUT

import sys
import contact_probability as cp
import numpy as np
import mdtraj as md




def test_evaluate_r_initial():#pairs_section, initial_distances):
    """
    Function to test the evaluation of  initial pairwise distances \
    accordingly the simulated model.
    Input:
     contacts_list - Array with the definitions from pairs section of the \
     forcefield (TPR file.)
     model - AA = All-Atom; CA = Carbon_alpha Coarse-Grained
    Output:
     r_initial - vector with the initial distance of each pair.
    """
    pairs_section_filepath = '../share/ci2-AA-contacts.dat'
    pairs_section = np.genfromtxt(pairs_section_filepath)
    pairs_indexes = np.subtract(pairs_section[:, 0:2], 1)
    pdb_file = '../share/ci2-adjusted.pdb'
    pdb_loaded = md.load(pdb_file)
    initial_distances = md.compute_distances(pdb_loaded, pairs_indexes)
    r_initial = cp.evaluate_r_initial(pairs_section)
    test = (np.subtract(r_initial, initial_distances) < 0.001).all()
    assert test


#raw_prob, contacts, atoms = gen_contact_probability('ci2-AA-EM-end.gro', 'ci2-AA-120-run.xtc', pairs_indexes, r_initial)

