#!/usr/bin/env python3
# usage: pylint

#import sys
import numpy as np
import mdtraj as md
import contact_probability as cp
import contacts_chunk as cc


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
    threshold_distances = 0.001
    pairs_section_filepath = 'share/ci2-AA-contacts.dat'
    pairs_section = np.genfromtxt(pairs_section_filepath)
    pairs_indexes = np.subtract(pairs_section[:, 0:2], 1)
    pdb_file = 'share/ci2-adjusted.pdb'
    pdb_loaded = md.load(pdb_file)
    initial_distances = md.compute_distances(pdb_loaded, pairs_indexes)
    r_initial = cp.evaluate_r_initial(pairs_section)
    test = (np.subtract(r_initial, \
        initial_distances) <= threshold_distances).all()
    assert test

def test_gen_contact_probability():
    """
    Function to test the gen_contact_probability if below given thresholds.
    """
    threshold_prob = 0.01
    threshold_contacts = 0
    threshold_atoms = 0
    pairs_section_filepath = 'share/ci2-AA-contacts.dat'
    pairs_section = np.genfromtxt(pairs_section_filepath)
    r_initial = cp.evaluate_r_initial(pairs_section)
    pairs_indexes = np.subtract(pairs_section[:, 0:2], 1)
    defined_raw_output = np.genfromtxt('share/raw-output-ci2.dat')
    defined_atoms_involved = \
        np.genfromtxt('share/atoms-involved-output-ci2.dat')
    defined_contacts_involved = \
        np.genfromtxt('share/Q-involved-output-ci2.dat')
    raw_prob, contacts, atoms = cp.gen_contact_probability(
        'share/ci2-AA-120-run.gro', 'share/ci2-AA-120-run.xtc',
        pairs_indexes, r_initial
        )
    test_raw = (np.subtract(raw_prob, \
        defined_raw_output) <= threshold_prob).all()
    test_contacts = (np.subtract(contacts, \
        defined_contacts_involved) <= threshold_contacts).all()
    test_atoms = (np.subtract(atoms, \
        defined_atoms_involved) <= threshold_atoms).all()
    test = test_raw and test_contacts and test_atoms
    assert test


def test_contacts_chunk():
    """
    Function to test if contacts chunk is correct within the error margin.
    """
    error_margin_contacts = 1
    pairs_section_filepath = 'share/ci2-AA-contacts.dat'
    pairs_section = np.genfromtxt(pairs_section_filepath)
    r_initial = cp.evaluate_r_initial(pairs_section)
    pairs_indexes = np.subtract(pairs_section[:, 0:2], 1)
    defined_contacts_output = np.genfromtxt('share/contacts-output-ci2.dat')
    contacts = cc.evaluating_contacts_chunk('share/ci2-adjusted.pdb', \
        'share/ci2-AA-120-run.xtc', pairs_indexes, r_initial)
    test = (np.absolute(np.subtract(contacts, defined_contacts_output))\
            <= error_margin_contacts).all()
    assert test