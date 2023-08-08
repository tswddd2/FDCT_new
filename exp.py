import os, sys

n_list = [500, 1000, 2000]
k_list = [10]
loop = 1

if len(sys.argv) == 2:
    path = sys.argv[1]
else: path = "~/boost_1_55_0/"

f = open("log.txt", 'a')
f.write("loop: {}\n".format(loop))

os.system("g++ databuilder_paper.cpp -o EXEC_exp")
os.system("g++ main_new.cpp -I {} -o EXEC_new".format(path))
os.system("g++ main_old.cpp -I {} -o EXEC_old".format(path))

for n in n_list:
    for k in k_list:
        os.system("./EXEC_exp -n {} -k {}".format(n, k))
        old = [0,0,0]
        new = [0,0,0]
        for l in range(loop):

            # for FDCT_new
            st = os.popen("./EXEC_new").read().split()
            st = [eval(_) for _ in st]
            new[0] += st[0] / loop
            new[1] += (st[1]-st[0]) / loop
            new[2] += st[1] / loop

            # for FDCT_old
            st = os.popen("./EXEC_old").read().split()
            st = [eval(_) for _ in st]
            old[0] += st[0] / loop
            old[1] += (st[1]-st[0]) / loop
            old[2] += st[1] / loop

            # assert output of both method is the same
            assert os.popen("diff oup.txt oup2.txt").read() == ''

        f.write("new n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}\n".format(n, k, new[0], new[1], new[2]))
        f.write("old n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}\n".format(n, k, old[0], old[1], old[2]))
        print("new n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}".format(n, k, new[0], new[1], new[2]))
        print("old n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}".format(n, k, old[0], old[1], old[2]))
f.close()