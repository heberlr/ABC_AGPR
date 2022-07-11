'''
    File name: kernel_comparison.py
    Author: Joao Vitor de Oliveira Silva
    Date created: 2/17/2022
    Python Version: 3.6
    Description: Code to present GP kernel comparison as surrogate models of a Verhulst tumor growth model.
'''

import numpy as np
from matplotlib import pyplot as plt
import matplotlib

# Latex font rendering
# matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
# matplotlib.rc('text', usetex = True)


# Loading data
data = np.loadtxt("modeldata.dat", delimiter = '\t')
t = data[:, 0]
qoi = data[:, 1]
qoi_RBF, qoi_RBF_std = np.loadtxt("sqrexp_qoi.txt"), np.loadtxt("sqrexp_qoi_std.txt")
qoi_Matern_15, qoi_Matern_15_std = np.loadtxt("Matern_qoi.txt"), np.loadtxt("Matern_qoi_std.txt")
qoi_Dot, qoi_Dot_std = np.loadtxt("Dot_qoi.txt"), np.loadtxt("Dot_qoi_std.txt")
qoi_Per, qoi_Per_std = np.loadtxt("Expsinsquared_qoi.txt"), np.loadtxt("Expsinsquared_qoi_std.txt")

plt.figure(figsize=(8,6.3))

# Old plot
# plt.errorbar(t, qoi, 200, marker='o', color='k', label='Data', linestyle = "None", alpha=0.8)
#plt.plot(t, qoi_RBF, marker='^', color='#5c607d', label='Sqr. exp', markersize = 10, alpha=0.92)
#plt.plot(t, qoi_Matern_15, marker='v', color='#0b57b0', label='Matérn ($\\nu = 1.5)$', linestyle = "dashed", markersize = 10, alpha=0.62)
#plt.plot(t, qoi_Dot, marker='h', color='#0b8ab0', label='Dot', markersize = 9, linestyle = "dotted", alpha=0.62)
#plt.plot(t, qoi_Per, marker='d', color='#1529bf', label='Periodic', markersize = 9, alpha=0.92)
# qoi_RBF_std = np.array(len(t)*[qoi_RBF_std])
# qoi_Matern_15_std = np.array(len(t)*[qoi_Matern_15_std])
# qoi_Dot_std = np.array(len(t)*[qoi_Dot_std])
# qoi_Per_std = np.array(len(t)*[qoi_Per_std])

# New plot
plt.scatter(0,1000.0,200,color='w')
plt.errorbar(t, qoi, 200, marker='o', color='k', label='Data', linestyle = "None", alpha=0.8)
plt.plot(t, qoi_RBF, marker='^', color='#8B4513', label='Sqr. exp', markersize = 10, alpha=0.95)
plt.fill_between(t,  qoi_RBF - 1.96 * qoi_RBF_std, qoi_RBF + 1.96 * qoi_RBF_std, color='#8B4513', alpha=0.12)
plt.plot(t, qoi_Matern_15, marker='v', color='#2E8B57', label='Matérn ($\\nu = 1.5)$', linestyle = "dashed", markersize = 10, alpha=0.62)
plt.fill_between(t,  qoi_Matern_15 - 1.96 * qoi_Matern_15_std, qoi_Matern_15 + 1.96 * qoi_Matern_15_std, color='#2E8B57', alpha=0.12)
plt.plot(t, qoi_Dot, marker='h', color='red', label='Dot', markersize = 9, linestyle = "dotted", alpha=0.62)
# plt.fill_between(t,  qoi_Dot - 1.96 * qoi_Dot_std, qoi_Dot + 1.96 * qoi_Dot_std, color='red', alpha=0.12)
plt.plot(t, qoi_Per, marker='d', color='#1E90FF', label='Periodic', markersize = 9, alpha=0.92)
# plt.fill_between(t,  qoi_Per - 1.96 * qoi_Per_std, qoi_Per + 1.96 * qoi_Per_std, color='#1E90FF', alpha=0.12)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('Time (hours)', size = 16)
plt.ylabel('Number of Cells', size = 16)
# plt.xlim([0,290])
#plt.ylim([-800, 11000])
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
          fancybox=True, shadow=True, ncol=3, fontsize = 16)

plt.savefig("kernel_comp.jpg", dpi=120,bbox_inches='tight')
# plt.show()
