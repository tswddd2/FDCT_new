import os

n_list = [500, 1000, 1500, 2500, 5000, 10000]
k_list = [100]
loop = 1

# path = r"C:/dev/boost_1_58_0/"
path = "boost_1_55_0/"

f = open("log.txt", 'a')
f.write("loop: {}\n".format(loop))

os.system("g++ databuilder.cpp -o EXEC_exp")
os.system("g++ main_new.cpp -I {} -o EXEC_new".format(path))
os.system("g++ main_old.cpp -I {} -o EXEC_old".format(path))

for n in n_list:
    for k in k_list:
        os.system("EXEC_exp -n {} -k {}".format(n, k))
        old = [0,0,0]
        new = [0,0,0]
        for l in range(loop):

            # for FDCT_new
            st = os.popen("EXEC_new").read().split()
            st = [eval(_) for _ in st]
            new[0] += st[0] / loop
            new[1] += (st[1]-st[0]) / loop
            new[2] += st[1] / loop

            # for FDCT_old
            st = os.popen("EXEC_old").read().split()
            st = [eval(_) for _ in st]
            old[0] += st[0] / loop
            old[1] += (st[1]-st[0]) / loop
            old[2] += st[1] / loop

        f.write("new n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}\n".format(n, k, new[0], new[1], new[2]))
        f.write("old n:{} k:{} st1:{:.2f} st2:{:.2f} tot:{:.2f}\n".format(n, k, old[0], old[1], old[2]))

f.close()