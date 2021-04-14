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

    # Extracts the edgelist from the MATLAB file.
    mat = scipy.io.loadmat(filename)
    idx = scipy.sparse.find(mat['A'])
    edgelist = np.transpose(np.array(idx)[:-1,:])
    np.savetxt(filename.replace('.mat', '.txt'), edgelist, fmt="%15d", delimiter=" ", header=" SourceVertex    TargetVertex")

    with tarfile.open(filename.replace('.mat', '.txt.tar.xz'), "w:xz") as tar:
        tar.add(filename.replace('.mat', '.txt'))

    os.remove(filename.replace('.mat', '.txt'))
