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
bool sim = true;
vector<int> inter;
int dbg = 0;


void readArg(int argc, char **argv) {
    if (argc == 1)
        return ;
    int t_opt;
    m = -1; r = 0.2;
    while ((t_opt = getopt(argc, argv, "bn:r:m:k:o:d") ) != -1)
        switch (t_opt)
        {
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
        m = max(20, n/20);
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

    void remove_subtree(TreeNode* now) {
        vector<TreeNode*>::iterator it;
        for (it = S(F(now)).begin(); it != S(F(now)).end(); it++)
            if (*it == now) {
                S(F(now)).erase(it);
                break;
            }
    }

    void add_subtree(TreeNode* now, TreeNode* father) {
        F(now) = father;
        S(father).push_back(now);
    }
};

Tree *trees, base;

void buildRandomTree(Tree &t, int n) {
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
        else if (F(now) != NULL && Frand() < r)
            t.remove_node(now);
    }

}

void cloneTree(Tree &t, TreeNode* clone, TreeNode* now = NULL) {
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
    vector<TreeNode*> ancestor;
    for (int i = 0; i < 5; i++) {
        TreeNode *u = t.at(inter[Lrand() % inter.size()]);
        TreeNode *v = t.at(inter[Lrand() % inter.size()]);
        if (u->removed || v->removed)
            cout << "ERROR\n";
        if (u!=v && F(u) != NULL) {
            TreeNode* f = v;
            while (f != u && f != NULL)
                f = F(f);
            if (f == u) {
                f = u; u = v; v = f;
            }
            t.remove_subtree(u);
            t.add_subtree(u, v);
            return;
        }
    }
}

string printTree(TreeNode* now) {

    if (L(now) > -1)
        return to_string(L(now)+1);

    string s = "", res;
    bool brc = false;
    for (auto son : S(now)) {
        res = printTree(son);
        if (res.length() > 0) {
            if (s.length() > 0) {
                s += ',';
                brc = true;
            }
            s += res;
        }
    }
    if (brc) s = "("+s+")";

    if (dbg) cout << now->nid << '\t' << s << endl;
    return s;
}



int main(int argc, char **argv) {
    int seed = time(0);
    // int seed = 42;
    srand(seed);
    readArg(argc, argv);

    trees = new Tree[k];

    ofstream fout(out_file);

    if (sim) {
        buildRandomTree(base, n);
        for (int ti = 0; ti < k; ti++) {
            cloneTree(trees[ti], base.root);
            inter.clear();
            for (auto node: trees[ti].nodes)
                if (!node->removed && F(node)!=NULL && node->label==-1)
                    inter.push_back(node->nid);
            for (int i = 0; i < m; i++)
                mutation(trees[ti]);
            fout << printTree(trees[ti].root) << '\n';
        }
    } else
        for (int ti = 0; ti < k; ti++) {
            buildRandomTree(trees[ti], n);
            fout << printTree(trees[ti].root) << '\n';
        }
    return 0;
}
