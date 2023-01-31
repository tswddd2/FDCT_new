#include <iostream>
#include <fstream>
#include <unistd.h>
#include <random>
#include <ctime>
using namespace std;

#define Lrand() (rand() * 30000 + rand())
#define Frand() ((rand()+0.)/RAND_MAX)
#define F(now) now->father
#define S(now) now->sons
#define L(now) now->label

int n = 100, m = 0, k = 10;
float r = 0;
string out_file = "inp.txt";
bool bin = true, sim = true;



void readArg(int argc, char **argv) {
    if (argc == 1)
        return ;
    bin = false;
    int t_opt;
    m = -1;
    while ((t_opt = getopt(argc, argv, "bn:r:m:k:o:d") ) != -1)
        switch (t_opt)
        {
            case 'b':
                bin = true;         break;
            case 'n':
                n = atoi(optarg);   break;
            case 'm':
                m = atoi(optarg);   break;
            case 'k':
                k = atoi(optarg);   break;
            case 'r':
                r = stof(optarg);   break;
            case 'o':
                out_file = optarg;  break;
            case 'd':
                sim = false;         break;
        }
    if (m == -1 && sim)
        m = min(20, n / 30);
    if (bin)
        r = 0;
}

struct TreeNode
{
	int label = -1, nid;
	vector<TreeNode*> sons;
    TreeNode* father = NULL;
    bool removed = false;
};

class Tree
{
    public:
    vector<TreeNode*> nodes;
    TreeNode* root;
    
    TreeNode* at(int nid) {
        return nodes[nid];
    }

    TreeNode* add_node(TreeNode* father = NULL) {
        TreeNode* now = new TreeNode();
        now->nid = nodes.size();
        F(now) = father;
        if (father != NULL)
            S(father).push_back(now);
        nodes.push_back(now);
        if (nodes.size() == 1)
            root = nodes[0];
        return now;
    }

    TreeNode* remove_node(TreeNode* now) {
        vector<TreeNode*>::iterator it;
        for (it = S(F(now)).begin(); it != S(F(now)).end(); it++)
            if (*it == now) {
                S(F(now)).erase(it);
                break;
            }

        for (auto son : S(now)) {
            F(son) = F(now);
            S(F(now)).push_back(son);
        }
        now->removed = true;
        return F(now);
    }

    TreeNode* remove_subtree(TreeNode* now) {
        vector<TreeNode*>::iterator it;
        for (it = S(F(now)).begin(); it != S(F(now)).end(); it++)
            if (*it == now) {
                S(F(now)).erase(it);
                break;
            }

        if (S(F(now)).size() > 1)
            return F(now);
        else
            return remove_node(F(now));
    }

    void add_subtree_bin(TreeNode* now, TreeNode* father) {
        TreeNode* son = S(father)[0], *inter;
        S(father).erase(S(father).begin());
        inter = add_node(father);
        S(inter).push_back(son);
        S(inter).push_back(now);
        F(son) = F(now) = inter;
    }

    void add_subtree(TreeNode* now, TreeNode* father, bool bin = false) {
        if (bin)
            add_subtree_bin(now, father);
        else {
            F(now) = father;
            S(father).push_back(now);
        }
    }
};

Tree *trees, base;

void buildRandomTree(Tree &t, int n, bool bin) {
    t.add_node();
    for (int i = 1; i < n; i++) {
        TreeNode* now;
        do {
            now = t.nodes[Lrand() % t.nodes.size()];
        } while (S(now).size() > 0);

        t.add_node(now);
        t.add_node(now);
    }

    int label = 0;
    for (auto now : t.nodes) {
        if (S(now).size() == 0)
            L(now) = label++;
        else if (!bin && F(now) != NULL && Frand() < r)
            t.remove_node(now);
    }

}

void cloneTree(Tree &t, TreeNode* clone, TreeNode* now = NULL) {
    // cout << clone->nid << '\n' << "sons ";
    // for (auto c : S(clone))
    //     cout << c->nid << ' ';
    // cout << '\n';
    if (now == NULL) {
        now = t.add_node();
        L(now) = L(clone);
    }
    for (auto clone_son : S(clone)) {
        TreeNode* son = t.add_node(now);
        L(son) = L(clone_son);
        cloneTree(t, clone_son, son);
    }
}

void mutation(Tree &t, bool bin = false) {
    if (Lrand() % 2 == 0 && !bin) {
        for (int i = 0; i < 5; i++) {
            TreeNode *now = t.at(Lrand() % t.nodes.size());
            if (L(now) == -1 && F(now) != NULL && !now->removed) {
                t.remove_node(now);
                return;
            }
        }
    } else {
        vector<TreeNode*> ancestor;
        for (int i = 0; i < 5; i++) {
            TreeNode *now = t.at(Lrand() % t.nodes.size());
            if (!now->removed && F(now) != NULL && F(F(now)) != NULL) {
                TreeNode* f = t.remove_subtree(now);
                ancestor.clear();
                while (f != NULL && ancestor.size() < 5) {
                    ancestor.push_back(f);
                    f = F(f);
                }
                t.add_subtree(now, ancestor[Lrand()%ancestor.size()], bin);
                return;
            }
        }
    }
}

string printTree(TreeNode* now) {
    if (L(now) > -1)
        return to_string(L(now)+1);

    string s = "(";
    for (auto son : S(now)) {
        if (s.length() > 1)
            s += ',';
        s += printTree(son);
    }
    s += ')';
    return s;
}



int main(int argc, char **argv) {
    int seed = time(0);
    srand(seed);
    readArg(argc, argv);

    trees = new Tree[k];

    ofstream fout(out_file);

    if (sim) {
        buildRandomTree(base, n, bin);
        for (int ti = 0; ti < k; ti++) {
            cloneTree(trees[ti], base.root);
            for (int i = 0; i < m; i++)
                mutation(trees[ti], bin);
        }
    } else {
        for (int ti = 0; ti < k; ti++) {
            buildRandomTree(trees[ti], n, bin);
        }
    }
    // fout << n << '\n';
    for (int ti = 0; ti < k; ti++)
        fout << printTree(trees[ti].root) << '\n';
    // cout << seed << '\n';
    return 0;
}
