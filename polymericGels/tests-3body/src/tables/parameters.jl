"""
    Parameters for the tables
"""

# General parameters
M=7;
N = 2^M;
sig = 0.4;
rc=1.5*sig;
rmin = sig/2;
rmax = rc;

# Parameters for two body interaction
eps = 1;

# Parameters for threebody interaction
eps_ij = 1.0;
eps_ik = 1.0;
eps_jk = 1.0;
thi = 180/(4*N)
thf = 180 - thi;
w=1.1;

