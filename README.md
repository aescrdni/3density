Supplementary code to calculate Hi-C based density features for CRISPR efficiency prediction. Written in Python 3.9. The code is available under MIT license.

# Data structure and assumptions
The main function, ```calc_3D_features```, requires upper triangular distance matrices (whose generation is described in the manuscript) in scipy sparse files; the code assumes that the distance matrix for each chromosome is stored in ```<data_path>/<chr_name>.npz```, with ```data_path``` and ```chr_name``` given as inputs.
The function receives the genomic coordinates of a CRISPR target site and returns the values of 9 3D features.

# Input
1. ```data_path```
The path to the distance matrix files.

3. ```chr_name```
The name of the target site's chromosome (for which a sparse distance matrix exists in ```data_path```).

4. ```site_coord```
Genomic coordinate of the target site's cut site on the (+) strand.

5. ```reso```
Resolution of the distance matrix (e.g. 10000).

# Output
```feats``` a dict with 9 3D features estimating the input target site's spatial density.

# Credits
Authors: Shaked Bergman and Tamir Tuller.
