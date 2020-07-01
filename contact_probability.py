#!/usr/bin/env python3
#usage: python3 contact_probability.py #

import os
import sys
import glob
import numpy as np
import mdtraj as md
import multiprocessing

THRESHOLD = 1.5
CHUNK = 10000

#input: pdb, xtc, file.cont output_name
#
# read contact list
# evaluate r0 for each pair
# compute distances for each pair for each step
# evaluate contact_probability
# end

cores = multiprocessing.cpu_count()

pdb_file = sys.argv[1]

xtc_file = sys.argv[2]

pairs_contacts = sys.argv[3]

output_file = sys.argv[4]

contacts = np.genfromtxt(pairs_contacts)

# This MUST be done due to python indexes, that starts from zero
pairs_indexes = np.subtract(contacts[:, 0:2], 1)

r_initial = np.power(np.divide(np.multiply(contacts[:, 4], 2), contacts[:, \
                                           3]), np.divide(1, 6))

cutoff = np.multiply(THRESHOLD, r_initial)

probability = np.zeros(np.shape(r_initial))

n_frames = 0


for chunk_trajectory in md.iterload(xtc_file, top=pdb_file, chunk=CHUNK):
    trajectory = md.compute_distances(chunk_trajectory, pairs_indexes)
    print((chunk_trajectory))
    # Getting the number of frames of each chunk and adding to the total
    n_frames += np.shape(trajectory)[0]
    # for idx in range(np.shape(pairs_indexes)[0]):
    for idj in range(np.shape(trajectory)[0]):
        below_threshold = np.less_equal(trajectory[idj], cutoff)
        probability += np.multiply(below_threshold, 1)


final_probability = np.divide(probability, n_frames)
assert not np.any(np.less(final_probability, 0))




#output: contact_probability
