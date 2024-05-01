from os.path import join, isfile
import scipy.sparse
import numpy as np
import math

node_num_cutoff = [1000, 2500, 5000, 10000]
cutoff_perc = [0.1, 0.15, 0.2, 0.25, 0.3]


def calc_3D_features(data_path, chr_name, site_coord, reso):
    # Loading the dists file
    dists_file = join(data_path, f"{chr_name}.npz")
    assert isfile(dists_file), f"{dists_file} not found"

    dists = scipy.sparse.load_npz(dists_file)
    dists = dists.todense()
    dists = np.array(dists)

    # Finding the cutoffs
    dist_vec = dists.flatten()
    dist_vec = dist_vec[np.nonzero(dist_vec)]
    dist_vec.sort()
    cutoff = []
    for c in cutoff_perc:
        cutoff.append(dist_vec[math.floor(len(dist_vec) * c)])
    del dist_vec

    # Symmetrizing the triangular matrix
    dists = dists + dists.T - np.diag(np.diag(dists))

    feats = {}
    node_name = round(site_coord / reso)
    if node_name < len(dists):
        dist_row = np.sort(dists[node_name, :])
        dist_row = dist_row[np.nonzero(dist_row)]
        for c in range(len(cutoff)):
            feats[f"density_{100 * cutoff_perc[c]:.0f}"] = sum(dist_row < cutoff[c])

        for c in node_num_cutoff:
            if len(dist_row) > 0:
                mean_dist = dist_row[1:c].mean()
                feats[f"dist_close_{c}"] = mean_dist

    return feats
