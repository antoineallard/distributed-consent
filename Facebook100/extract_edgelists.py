# -​*- coding: utf-8 -*​-
# @author: Antoine Allard <antoineallard.info>

import glob
import numpy as np
import os
import scipy.io
import scipy.sparse
import tarfile


for filename in glob.glob('*.mat'):

    if filename == 'schools.mat':
        continue

    rootname = filename.replace('.mat', '').replace(' ', '')
    edgelist_filename = rootname + '.txt'
    archive_filename = rootname + '.txt.tar.xz'

    if not os.path.isfile(archive_filename):

        print(rootname)

        # Extracts the edgelist from the MATLAB file.
        mat = scipy.io.loadmat(filename)
        idx = scipy.sparse.find(mat['A'])
        edgelist = np.transpose(np.array(idx)[:-1,:])
        np.savetxt(edgelist_filename, edgelist, fmt="%15d", delimiter=" ", header=" SourceVertex    TargetVertex")

        # Compresses the edgelist.
        with tarfile.open(archive_filename, "w:xz") as tar:
            tar.add(edgelist_filename)

        # Removes the edgelist.
        os.remove(edgelist_filename)
