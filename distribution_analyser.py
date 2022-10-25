import os
import numpy as np
import matplotlib.pyplot as plt
import sys

run_data_dir, t1, t2 = sys.argv[1:]

with open(os.path.join(run_data_dir, "rmd2_infile.dat")) as infile:
    lines = infile.readlines()
    n_grid_points = int(lines[2])
    grid_extent = float(lines[5])

t1 = int(t1)
t2 = int(t2)

p_files = sorted(os.listdir(os.path.join(run_data_dir,"p_hist")))
p_files = [os.path.join(run_data_dir+"/p_hist", p_file) for p_file in p_files]
p_hist = np.zeros(shape=(t2-t1,n_grid_points ))
for i, p_file in enumerate(p_files[t1:t2]):
    a = np.loadtxt(p_file)
    p_hist[i,:] = a
q_hist = np.loadtxt(os.path.join(run_data_dir, "rmd2_qfile.dat"))

D_KL_RMDs = np.zeros(shape=(n_grid_points))
D_KL_Gausss = np.zeros(shape=(n_grid_points))
print("Calculating D_KLs...")
for slice_idx in range(n_grid_points): 

    p_slice = p_hist[:, slice_idx]
    q_slice = q_hist[:, slice_idx]
    q_mean = np.mean(q_slice)
    q_std = np.std(q_slice)
    

    biggest_x = 0

    for x in p_slice:
        if x > biggest_x:
            biggest_x = x
    for x in q_slice:
        if x > biggest_x:
            biggest_x = x

    x = np.linspace(0, 12, 1000)

    p_kde = np.full_like(x, 1e-12)
    q_kde = np.full_like(x, 1e-12)

    sig = biggest_x/20

    for m in p_slice:
        p_kde += 1/(sig*np.sqrt(np.pi*2))*np.exp(-0.5*((x-m)/(sig))**2)

    for m in q_slice:
        q_kde += 1/(sig*np.sqrt(np.pi*2))*np.exp(-0.5*((x-m)/(sig))**2)

    p_kde /= np.trapz(p_kde, x)
    q_kde /= np.trapz(q_kde, x)

    
    theory_kde = 1e-12 + 1/(q_std*np.sqrt(np.pi*2))*np.exp(-0.5*((x-q_mean)/(q_std))**2)
    theory_kde /= np.trapz(theory_kde, x)

    D_KL_RMD = np.sum(p_kde * np.log(p_kde/q_kde), axis=0)
    D_KL_Gauss = np.sum(theory_kde * np.log(theory_kde/q_kde), axis=0)
    D_KL_RMDs[slice_idx] = D_KL_RMD
    D_KL_Gausss[slice_idx] = D_KL_Gauss


mean_kl_rmd = np.mean(D_KL_RMDs)
std_kl_rmd = np.std(D_KL_RMDs)
mean_kl_gauss = np.mean(D_KL_Gausss)
std_kl_gauss = np.std(D_KL_Gausss)

print("\t\t\t Mean \t\t\t STD")
print("D_KL(RMD):\t\t",mean_kl_rmd, "\t", std_kl_rmd)
print("D_KL(Theory):\t\t", mean_kl_gauss, "\t",std_kl_gauss)

"""
plt.plot(x, p_kde, label="RMD")
plt.plot(x, q_kde, label="Synthetic Target (ST)")
plt.plot(x, theory_kde, label="Gaussian")
plt.plot([], [], ' ', label="$D_{KL}(RMD, ST)$ = " + f"{D_KL_RMD:.2f}")
plt.plot([], [], ' ', label="$D_{KL}(Gaussian, ST)$ = " + f"{D_KL_Gauss:.2f}")
plt.xlabel("$g(r)|_{r=r'}$")
plt.ylabel("Frequency")
plt.legend()
plt.show()
"""
