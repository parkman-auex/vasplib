#!/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import z2pack

# Creating the System. Note that the SCF charge file does not need to be
# copied, but instead can be referenced in the .files file.
# The k-points input is appended to the .in file
# The command (mpirun ...) will have to be replaced to match your system.
system = z2pack.fp.System(
    input_files=[
        "input/CHGCAR", "input/INCAR", "input/POSCAR", "input/POTCAR",
        "input/wannier90.win"
    ],
    kpt_fct=z2pack.fp.kpoint.vasp,
    kpt_path="KPOINTS",
    command="mpirun -np 16 /opt/vasp.5.4.4/bin/vasp_ncl_w2.1 >& log"
)

if not os.path.exists('./results'):
    os.mkdir('./results')
if not os.path.exists('./plots'):
    os.mkdir('./plots')
if not os.path.exists('./WCC_data'):
    os.mkdir('./WCC_data')
# Running the WCC calculation - standard settings
# kz0
result_kz0 = z2pack.surface.run(
    system=system,
    num_lines=11,
    min_neighbour_dist=0.001,
    surface=lambda s, t: [s , t, 0],
    save_file='./results/result_kz0.p',
    load=True
)
# # kzP
# result_kzP = z2pack.surface.run(
#     system=system,
#     num_lines=11,
#     min_neighbour_dist=0.001,
#     surface=lambda s, t: [s , t, 0.5],
#     save_file='./results/result_kzP.p',
#     load=True
# )
# # kx0
# result_kx0 = z2pack.surface.run(
#     system=system,
#     num_lines=11,
#     min_neighbour_dist=0.001,
#     surface=lambda s, t: [0 , s, t],
#     save_file='./results/result_kx0.p',
#     load=True
# )
# # kxP
# result_kxP = z2pack.surface.run(
#     system=system,
#     num_lines=11,
#     min_neighbour_dist=0.001,
#     surface=lambda s, t: [0.5 , s, t],
#     save_file='./results/result_kxP.p',
#     load=True
# )
# # ky0
# result_ky0 = z2pack.surface.run(
#     system=system,
#     num_lines=11,
#     min_neighbour_dist=0.001,
#     surface=lambda s, t: [t , 0, s],
#     save_file='./results/result_ky0.p',
#     load=True
# )
# # kyP
# result_kyP = z2pack.surface.run(
#     system=system,
#     num_lines=11,
#     min_neighbour_dist=0.001,
#     surface=lambda s, t: [t , 0.5, s],
#     save_file='./results/result_kyP.p',
#     load=True
# )

# print(result_0.t) # prints the positions of the lines
# print(result_0.pol) # prints the sum of WCC for each line

# Plotting WCC evolution
# gen ax
ax1 = plt.subplot(111) #kz0
# ax2 = plt.subplot(322) #kzP
# ax3 = plt.subplot(323) #kx0
# ax4 = plt.subplot(324) #kxP
# ax5 = plt.subplot(325) #ky0
# ax6 = plt.subplot(326) #kyP
#
ax1.set_ylim([0.0, 1])
# ax2.set_ylim([0.0, 1])
# ax3.set_ylim([0.0, 1])
# ax4.set_ylim([0.0, 1])
# ax5.set_ylim([0.0, 1])
# ax6.set_ylim([0.0, 1])
#
ax1.set_xticks([0,0.5,1])
# ax2.set_xticks([0,0.5,1])
# ax3.set_xticks([0,0.5,1])
# ax4.set_xticks([0,0.5,1])
# ax5.set_xticks([0,0.5,1])
# ax6.set_xticks([0,0.5,1])
#fig, ax = plt.subplots(1, 2, sharey=True, figsize=(6, 6))
# data function

def write_WCC(results,filename):
    WCC_f = open(filename, "w")
    WCC_f.write('#      k        sum(wcc(:,ik))  wcc(i, ik)(i=1, NumOccupied)\n')
    for line in results.lines:
        WCC_f.write('{0:.6} '.format(line.t))
        for x in line.wcc:
            WCC_f.write('{0:.6} '.format(x));
        WCC_f.write('\n')
    WCC_f.close()    

# data_gen
write_WCC(result_kz0,'./WCC_data/WCC_kz0_FILENAME.dat')
# write_WCC(result_kzP,'./WCC_data/WCC_kzP_FILENAME.dat')
# write_WCC(result_kx0,'./WCC_data/WCC_kx0_FILENAME.dat')
# write_WCC(result_kyP,'./WCC_data/WCC_kxP_FILENAME.dat')
# write_WCC(result_ky0,'./WCC_data/WCC_ky0_FILENAME.dat')
# write_WCC(result_kyP,'./WCC_data/WCC_kyP_FILENAME.dat')
#plot_data
z2pack.plot.wcc(result_kz0, axis=ax1, wcc_settings={'s': 5.0, 'lw': 1.0, 'facecolor': 'r', 'edgecolors': 'r'},gaps=False)
# z2pack.plot.wcc(result_kzP, axis=ax2, wcc_settings={'s': 5.0, 'lw': 1.0, 'facecolor': 'r', 'edgecolors': 'r'},gaps=False)
# z2pack.plot.wcc(result_kx0, axis=ax3, wcc_settings={'s': 5.0, 'lw': 1.0, 'facecolor': 'r', 'edgecolors': 'r'},gaps=False)
# z2pack.plot.wcc(result_kxP, axis=ax4, wcc_settings={'s': 5.0, 'lw': 1.0, 'facecolor': 'r', 'edgecolors': 'r'},gaps=False)
# z2pack.plot.wcc(result_ky0, axis=ax5, wcc_settings={'s': 5.0, 'lw': 1.0, 'facecolor': 'r', 'edgecolors': 'r'},gaps=False)
# z2pack.plot.wcc(result_kyP, axis=ax6, wcc_settings={'s': 5.0, 'lw': 1.0, 'facecolor': 'r', 'edgecolors': 'r'},gaps=False)






#ax.set_yticks([0,0.5,1])
#z2pack.plot.wcc(result_1, axis=ax[1],wcc_settings={'s': 5.0, 'lw': 1.0, 'facecolor': 'r', 'edgecolors': 'r'},gaps=False)


# save
plt.savefig('plots/wilson_FILENAME.pdf', bbox_inches='tight')


# print(
#     'Z2 topological invariant at kx = 0: {0}'.format(
#         z2pack.invariant.z2(result_0)
#     )
# )

