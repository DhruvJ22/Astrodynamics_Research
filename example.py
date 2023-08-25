"""
Created on Tue Feb  8 22:13:15 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

This is an exmaple file to test the CR3BP functions
"""

import matplotlib.pyplot as plt
from cr3bp_char_quant import sys_chars
from cr3bp_lib_calc import lib_pt_loc
from cr3bp_PO_main import cr3bp_model
from cr3bp_PO_main import periodic_orbit
# from PO_master import periodic_orbit
sys_p1p2 = sys_chars('Earth','Moon')
mu = sys_p1p2.mu
lib_loc = lib_pt_loc(sys_p1p2)

ic = [1.05903, -0.067492, -0.103524, -0.170109, 0.0960234, -0.135279]
cr3bp_obj = cr3bp_model(sys_p1p2, ic)
JC0 = cr3bp_obj.JC(ic)

li = lib_loc[:, :]  # 0 for L1 and  1 for L2....
print("Earth-Moon Li:", li)

print("Jacobi constant:", JC0)

# Propagate the arbitrary state for time = 0 to tf
cr3bp_obj.tf = 10
results = cr3bp_obj.propagate()

po = periodic_orbit(sys_p1p2, ic)


# Plot P1, P2, Libration points and configuration space state history of the arbirtray state
pltnum = 1
plt.figure(pltnum)
ax = plt.axes(projection="3d")
ax.set_title("CR3BP EM, trajectory, T = 10[nd], tol = 1e-12")
ax.plot3D(results["states"][:, 0], results["states"][:, 1], results["states"][:, 2])
ax.scatter(
    results["states"][0, 0],
    results["states"][0, 1],
    results["states"][0, 2],
    color="black",
    label="t=0",
)
ax.scatter(li[:, 0], li[:, 1], li[:, 2], color="red", label="Li")
ax.scatter(-mu, 0, 0, color="blue", label="Earth")
ax.scatter(1 - mu, 0, 0, color="grey", label="Moon")
ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")])
ax.set_ylabel("y [nd]")
ax.set_xlabel("x [nd]")
ax.set_zlabel("z [nd]")
plt.legend()
pltnum = pltnum + 1
