import numpy as np

import matplotlib.pyplot as plt
import os 

import random


# color = ['r', 'b', 'g', 'darkorange', 'black', 'yellow', 'aqua', 'purple', 'm', 'Magenta']








fontsize = 20

log_file = "log/res.txt"

def get_template(log_file):
    template = []
    f = open(log_file, 'r')
    for line in f.readlines():
        if len(line) < 2 or line[:3] not in ['new', 'old']:
            continue
        ele = line.split()
        for x in ele:
            if x in ['new', 'old']:
                method = x
            if x[:2] == 'n:':
                n = eval(x[2:])
            if x[:2] == 'k:':
                k = eval(x[2:])
            if x[:4] == 'st1:':
                st1 = eval(x[4:])
            if x[:4] == 'st2:':
                st2 = eval(x[4:])
            if x[:4] == 'tot:':
                tot = eval(x[4:])
        # print(ele)
        template.append({"method":method, "n":n, "k":k, "st1":st1, "st2":st2, "tot":tot})
    f.close()
    return template



def drawline(template, _label = 'tmp'):

    x = []
    y = []

    for n in n_list:
        for k in k_list:
            for tup in template:
                if tup['n'] == n and tup['k'] == k and tup['method'] == method:
                    x.append(tup[xlabel])
                    y.append(tup[ylabel])

    plt.plot(x, y, label=_label, lw=3)









plt.figure(figsize=(8, 6))
plt.ylabel("avg time (s)", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)







tmp_coincide = False

plt.axes().set_ylim(0, 120)
# plt.axes().set_xlim(0, 250)
k_list = [500, 1000, 1500, 2500, 5000, 10000]
n_list = [100]
xlabel = 'k'


plt.xlabel("k", fontsize=fontsize)

template = get_template(log_file)

method = 'old'
ylabel = 'st1'
drawline(template, "cal weight")
ylabel = 'st2'
drawline(template, "filter cluster")
ylabel = 'tot'
drawline(template, "total")
























# set loc=9 for upper center legend
plt.legend(loc=0, fontsize=fontsize)
plt.tight_layout()
tmp = plt.gcf()
plt.show()
tmp.savefig("pic/tmp.png")