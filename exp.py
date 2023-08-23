import os, sys
from Bio import Phylo
from Bio.Phylo.Consensus import _BitString

def _bitstrs(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString("".join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs


def compare(tree1, tree2):
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        return False
    # true if _BitStrings are the same
    #print(_bitstrs(tree1))
    #print(_bitstrs(tree2))
    if _bitstrs(tree1) == _bitstrs(tree2):
        return True
    else:
        return False

k_list = [1000]
n_list = [10000]
different_trees = False
loop = 3

path = "/opt/homebrew/Cellar/boost/1.82.0_1/include/"

f = open("log.txt", 'a')
f.write("loop: {}\n".format(loop))

os.system("g++ -std=c++17 databuilder_paper.cpp -o EXEC_exp")
os.system("g++ -std=c++17 main_new.cpp -I {} -o EXEC_new".format(path))
os.system("g++ -std=c++17 main_old.cpp -I {} -o EXEC_old".format(path))

for n in n_list:
    for k in k_list:
        oldkn2 = [0,0,0]
        oldk2n = [0,0,0]
        new = [0,0,0]
        for l in range(loop):
            if different_trees:
                os.system("./EXEC_exp -d -n {} -k {}".format(n, k))
            else:
                os.system("./EXEC_exp -n {} -k {}".format(n, k))

            # for FDCT_new
            st = os.popen("./EXEC_new").read().split()
            st = [eval(_) for _ in st]
            new[0] += st[0] / loop
            new[1] += (st[1]-st[0]) / loop
            new[2] += st[1] / loop

            # for FDCT_old
            st = os.popen("./EXEC_old").read().split()
            st = [eval(_) for _ in st]
            oldkn2[0] += st[0] / loop
            oldkn2[1] += (st[1]-st[0]) / loop
            oldkn2[2] += st[1] / loop
            
            st = os.popen("./EXEC_old -k").read().split()
            st = [eval(_) for _ in st]
            oldk2n[0] += st[0] / loop
            oldk2n[1] += (st[1]-st[0]) / loop
            oldk2n[2] += st[1] / loop
            
            tree1 = Phylo.read("oup.txt", "newick")
            tree2 = Phylo.read("oup2.txt", "newick")
            tree3 = Phylo.read("oup3.txt", "newick")
            assert os.popen("diff oup.txt oup2.txt").read() == '' or compare(tree1, tree2)
            assert os.popen("diff oup.txt oup3.txt").read() == '' or compare(tree1, tree3)
            assert os.popen("diff oup2.txt oup3.txt").read() == '' or compare(tree2, tree3)

            # assert output of both method is the same
            #assert os.popen("diff oup.txt oup2.txt").read() == ''

        f.write("new n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}\n".format(n, k, new[0], new[1], new[2]))
        f.write("oldkn2 n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}\n".format(n, k, oldkn2[0], oldkn2[1], oldkn2[2]))
        f.write("oldk2n n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}\n".format(n, k, oldk2n[0], oldk2n[1], oldk2n[2]))
        print("new n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}".format(n, k, new[0], new[1], new[2]))
        print("oldkn2 n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}".format(n, k, oldkn2[0], oldkn2[1], oldkn2[2]))
        print("oldk2n n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}".format(n, k, oldk2n[0], oldk2n[1], oldk2n[2]))
f.close()