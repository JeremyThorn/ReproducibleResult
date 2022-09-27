import numpy as np
import matplotlib.pyplot as plt
import os
np.random.seed() 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#             RUN PARAMETERS             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
run_name = "synthetic_data"
n_atoms = 9 
n_dims = 3 
n_grid_points = 4000 
n_samples = 2000
cutoff = 1.0 
cutoff_width = 5e-3
grid_extent = 1.0
sig = 5e-3
##########################################

cell = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])

def P(x, sig, atoms, neigh_vecs, cutoff):
    p = np.zeros_like(x) + 1e-12
    for atomi in atoms:
        for vec in neigh_vecs:
            for atomj in atoms:
                if not (atomi == vec+atomj).all():
                    diff = atomi - (vec + atomj)
                    dist = np.linalg.norm(diff)
                    if dist <= cutoff-cutoff_width:
                        p += gaussian(x, dist, sig)
                    if dist <= cutoff and dist > cutoff-cutoff_width:
                        p += gaussian(x, dist, sig)*0.5*(1.0 + np.cos(np.pi/cutoff_width*(dist-cutoff)))
    p /= np.trapz(p, x)
    return p
 
def gaussian(x, m, s):
    g = np.exp(-0.5*((x-m)/(s))**2)
    return g

def get_neighbour_vecs(cell, cutoff, n_atoms, n_dims):
    neighbour_degree = np.ceil(cutoff)
    neighbour_direcs = np.arange(-neighbour_degree, neighbour_degree+1.0, 1)
    n_neighbours = ((2.0*neighbour_degree+1)**3).astype(int)
    extended_cell = np.einsum("ij,k->ikj", cell, neighbour_direcs)
    neighbour_vecs = np.zeros((n_neighbours, n_dims))
    c = 0
    for v0 in extended_cell[0]:
        for v1 in extended_cell[1]:
            for v2 in extended_cell[2]:
                neighbour_vecs[c,:] = v0 + v1 + v2
                c += 1
    return neighbour_vecs

def fill_cell_crystal(cell, n_atoms, n_dims):
    dirs = [1.0/4.0, -1.0/4.0]
    atoms = np.zeros((n_atoms, n_dims))
    c = 0
    for i in dirs:
        for j in dirs:
            for k in dirs:
                atoms[c,0] = i + np.random.normal(scale=0.03)
                atoms[c,1] = j + np.random.normal(scale=0.03)
                atoms[c,2] = k + np.random.normal(scale=0.03)
                c += 1

    atoms = np.einsum("ij,li->lj", cell, atoms)
    disp = 0.5*np.sum(cell, axis=0)
    atoms = atoms+disp
    return atoms

neigh_vecs = np.array(get_neighbour_vecs(cell, cutoff, n_atoms, n_dims))
x = np.linspace(0.0, grid_extent, n_grid_points)

ps = np.zeros((n_samples, n_grid_points))

print("Generating Synthetic Data...")
for i in range(n_samples):
    print("Generating sample ", i)
    atoms = fill_cell_crystal(cell, n_atoms, n_dims)
    ps[i] = P(x, sig, atoms, neigh_vecs, cutoff)

try:
    os.mkdir(run_name)
except:
    pass

np.savetxt(os.path.join(run_name,"rmd2_qfile.dat"), ps)

with open(os.path.join(run_name, "rmd2_infile.dat"), "w") as infofile:
    infofile.write(str(n_atoms) + "\n")
    infofile.write(str(n_dims) + "\n")
    infofile.write(str(n_grid_points) + "\n")
    infofile.write(str("...\n")) # for n_iters
    infofile.write(str(n_samples) + "\n")
    infofile.write(str(grid_extent) + "\n")
    infofile.write(str(cutoff) + "\n")
    infofile.write(str("...\n")) # for dt
    infofile.write(str(sig) + "\n")
    infofile.write(str(cutoff_width) + "\n")
    for row in cell:
        rowstr = ""
        for i in row:
            rowstr += f'{i:0.12f}'+ " "
        infofile.write(rowstr + "\n")
    for i in range(n_atoms):
        randpos = np.random.uniform(size=3)
        rowstr = ""
        for coord in randpos:
            rowstr += f'{coord:0.12f}'+ " "
        infofile.write(rowstr + "\n")
