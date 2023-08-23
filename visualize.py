import numpy as np

import matplotlib.pyplot as plt
import os 

import random


def read_log_file(log_file):
    k_list = []
    n_list = []
    v_list = []
    t_new_list = []
    t_oldkn2_list = []
    t_oldk2n_list = []
    f = open(log_file, 'r')
    for line in f.readlines():
        if len(line) < 2 or line[:3] not in ['new', 'old']:
            continue
        ele = line.split()
        for x in ele:
            if x[:3] in ['new', 'old']:
                method = x
            elif x[:2] == 'n:':
                n = eval(x[2:])
            elif x[:2] == 'k:':
                k = eval(x[2:])
            elif x[:4] == 'tot:':
                tot = eval(x[4:])
        if method == 'new':
            k_list.append(k)
            n_list.append(n)
            t_new_list.append(tot)
        elif method == 'oldkn2':
            t_oldkn2_list.append(tot)
        elif method == 'oldk2n':
            t_oldk2n_list.append(tot)
    f.close()
    if n_list[0] == n_list[1]:
        v_list = k_list
    else:
        v_list = n_list
    return n_list, k_list, v_list, t_new_list, t_oldkn2_list, t_oldk2n_list



log_file = "log.txt"
n_list, k_list, v_list, t_new_list, t_oldkn2_list, t_oldk2n_list = read_log_file(log_file)

plt.rc('font', family='serif', size=24)
fig, ax = plt.subplots(figsize=(10,7))

if k_list == v_list:
    ax.set_xlabel("k")
else:
    ax.set_xlabel("n")
ax.set_ylabel("s")

p1 = ax.plot(v_list, t_new_list, label="New",linewidth=3, markersize=10, alpha=0.8,marker='o')#,color='steelblue')
p2 = ax.plot(v_list, t_oldkn2_list, label="Old_kn2",linewidth=3, markersize=10, alpha=0.8,marker='^')#,color='indianred')
p3 = ax.plot(v_list, t_oldk2n_list, label="Old_k2n",linewidth=3, markersize=10, alpha=0.8,marker='s')#,color='indianred')

ax.legend(fancybox=True, framealpha=0.5, fontsize = 20)

fig.tight_layout()

# plt.show()

plt.savefig("visual.pdf",format='pdf')