This program takes as input either a 2D or a 3D cubic lattice of *classical* dipoles with up to two basis vectors and computes:
1. Frenkel eigenstates
2. Band structure
3. Macroscopic polarization

The user can set the type of interaction in the calc_interaction_X where X is either 2D or 3D to be dipole-dipole, nearest-neighbors, or to have a cutoff for convergence (especially important in 3D).
The program uses the tigh-binding model to compute the above quantities.
