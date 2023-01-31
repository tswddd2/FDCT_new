#include <iostream>
#include <algorithm>
// using namespace std;

struct radix_node_t {
    int key;
    void *val;
    radix_node_t *prev;

    void assign(int _key, void *_val) {
        key = _key, val = _val;
    }
};

struct radix_t {
    int node_num, key_num;
    int n = 0, k = 0;
    radix_node_t **tails, **nodes;
    radix_node_t **out;

    radix_t(int node_num, int key_num) : node_num(node_num), key_num(key_num), 
        out(new radix_node_t*[node_num]), tails(new radix_node_t*[key_num+1]), 
        nodes(new radix_node_t*[node_num]) {
            for (int i = 0; i < node_num; i++) {
                nodes[i] = new radix_node_t();
                out[i] = new radix_node_t();
            }
        }

    void clear() { n = k = 0; }

    void add(int key, void *val) {

        if (n == node_num) { cout << "Error radix sort n range\n"; exit(1); }

        if (key > k) k = key;
        nodes[n++]->assign(key, val);

        if (k > key_num)
            { cout << "Error radix sort k range\n"; exit(1); }
    }

    // stable O(k) sort
    void sort(bool desc = false) {
        static int i, j, key;

        for (i = 0; i <= k; i++) tails[i] = NULL;
        for (i = 0; i < n; i++) nodes[i]->prev = NULL;

        for (i = 0; i < n; i++) {
            key = nodes[i]->key;

            if (tails[key] == NULL)
                tails[key] = nodes[i];
            else {
                nodes[i]->prev = tails[key];
                tails[key] = nodes[i];
            }
        }

        j = n-1;
        static radix_node_t *node;
        for (i = 0; i <= k; i++) {
            node = desc ? tails[i] : tails[k-i];
            while (node != NULL) {
                *out[j] = *node;
                node = node->prev;
                j--;
            }
        }
    }

    // unstable O(nlogn) quick sort
    void quicksort(bool desc = false) {
        static auto comp_desc = [](radix_node_t *x, radix_node_t *y) {return x->key > y->key;};
        static auto comp = [](radix_node_t *x, radix_node_t *y) {return x->key < y->key;};
        static int i;
        for (i = 0; i < n; i++) *out[i] = *nodes[i];
        if (desc) std::sort(out, out+n, comp_desc);
        else std::sort(out, out+n, comp);
    }
    
};