/*
 * freqdiff2.h
 *
 *  Created on: 2021
 *      Author: Victor
 */

#ifndef FREQDIFF_H_
#define FREQDIFF_H_

#include <iostream>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <cassert>
#include <boost/dynamic_bitset.hpp>
using std::cout;
using std::endl;

#include "taxas_ranges.h"
#include "lca_preprocessing.h"
#include "radix_sort.h"
#include "utils.h"

#include "Tree.h"

#include <chrono>
using namespace std::chrono;
#define Now high_resolution_clock::now
#define duration_sec(start,end) (float)(duration_cast<milliseconds>(end-start).count() / 1000.)

struct node_bitvec_t {
	Tree::Node* node;
	boost::dynamic_bitset<>* bitvector;

	node_bitvec_t() : node(NULL), bitvector(NULL) {}
	node_bitvec_t(Tree::Node* node, boost::dynamic_bitset<>* bitvector) : node(node), bitvector(bitvector) {}
	//~node_bitvec_t() { delete bitvector; } FIXME
};
bool operator < (const node_bitvec_t& bv1, const node_bitvec_t& bv2) {
	return *(bv1.bitvector) < *(bv2.bitvector);
}

// cp_root (node id -> node*): show root of nodes
// cp_rmqs (node id -> rmq): only root id have rmq structure, store weight in (relative) depth order
// depths (node id -> depth): use to process cp_rmq, but not start from 0, need to converge
struct subpath_query_info_t {
	size_t nodes_num;
	Tree::Node** cp_roots;
	gen_rmq_t** cp_rmqs;
	int* depths;

	subpath_query_info_t(size_t nodes_num) : nodes_num(nodes_num), cp_roots(new Tree::Node*[nodes_num]),
			cp_rmqs(new gen_rmq_t*[nodes_num]), depths(new int[nodes_num]) {}
	~subpath_query_info_t() {
		delete[] cp_roots;
		for (size_t i = 0; i < nodes_num; i++) {
			delete cp_rmqs[i];
		}
		delete[] cp_rmqs;
		delete[] depths;
	}
};

subpath_query_info_t* preprocess_subpaths_queries(Tree* tree);
struct prob_set {
	Tree *tree1, *tree2;
	taxas_ranges_t *t1_tr, *t2_tr;
	lca_t *t2_lcas;
	subpath_query_info_t *subpq_t2;

	prob_set(Tree *t1, Tree *t2) : tree1(t1), tree2(t2), t1_tr(build_taxas_ranges(t1)), 
		t2_tr(build_taxas_ranges(t2)), t2_lcas(lca_preprocess(t2)), 
		subpq_t2(preprocess_subpaths_queries(t2)) {}
};

int* start,* stop;
int* e,* m;
std::vector<Tree::Node*>* rsort_lists;

bool* marked;

Tree::Node** _left,** _right;
size_t* orig_pos_in_parent;

size_t* counter;
bool* BT;
int* leaf_p_index;

int* vleft,* vright;
int* pointer;
int* levels,* ids;
int* parent;
bool* exists;
Tree::Node** tree_nodes;




bool g_tmp = true, check_ans = true;

































// working here

struct label_node {
	int id;
	label_node *prev_node, *father;
	int label = 0, lval = 0, rval = 0;
};

struct label_prob_set {
	std::vector<Tree*> trees;
	lca_t **lcas;
	int n, k;

	// record whether given node id is visited and store corresponding address of label node
	// see reference
	label_node **vis;
	std::vector<label_node> *inp;

	radix_t *radix;

	label_prob_set(std::vector<Tree*> &_trees) {
		k = _trees.size();
		n = _trees[0]->get_leaves_num();
		trees.resize(k);
		lcas = new lca_t*[k];
		inp = new std::vector<label_node>[k];

		int i, j;
		Tree::Node *node;
		for (i = 0; i < k; i++) {
			trees[i] = _trees[i];
			inp[i].reserve(2*n);
			for (j = 0; j < trees[i]->get_nodes_num(); j++)
				inp[i].push_back({j, NULL});
			for (j = 0; j < trees[i]->get_nodes_num(); j++) {
				node = trees[i]->get_node(j);
				if (node->parent == NULL) inp[i][j].father = NULL;
				else inp[i][j].father = &inp[i][node->parent->id];
			}
			lcas[i] = lca_preprocess(trees[i]);
		}

		vis = new label_node*[2*n];
		radix = new radix_t(2*n*k, 2*n*k);
	}
};

// build left / right node list given taxa range [st, ed)
std::vector<label_node>* build_sublist(label_prob_set *prob, std::vector<label_node> *lnodes, int st, int ed) {
	Tree::Node *prev, *node, *lca_node;
	label_node **vis = prob->vis;
	std::vector<label_node> *res;
	int n = ed - st, k = prob->k;
	int i, j;

	res = new std::vector<label_node>[k];

	for (i = 0; i < k; i++) {
		res[i].reserve(2*n);

		// map node id to label node address
		for (j = 0; j < lnodes[i].size(); j++)
			vis[lnodes[i][j].id] = &lnodes[i][j];

		prev = NULL;
		for (j = 0; j < lnodes[i].size(); j++) {

			// find leave node in taxa range [st, ed)
			node = prob->trees[i]->get_node(lnodes[i][j].id);
			if (!node->is_leaf()) continue;
			if (!(node->taxa >= st && node->taxa < ed)) continue;

			// add leave node
			res[i].push_back({node->id, vis[node->id]});

			// add lca node not added before
			if (prev != NULL) {
				lca_node = prob->trees[i]->get_node(lca(prob->lcas[i], node->id, prev->id));
				if (vis[lca_node->id] != NULL) {
					res[i].push_back({lca_node->id, vis[lca_node->id]});
					vis[lca_node->id] = NULL;
				}
			}
			prev = node;
		}

		label_node *plnode;
		for (j = 0; j < lnodes[i].size(); j++)
			vis[lnodes[i][j].id] = NULL;
		for (j = 0; j < res[i].size(); j++)
			vis[res[i][j].id] = &res[i][j];
		for (j = 0; j < res[i].size(); j++) {
			plnode = res[i][j].prev_node->father;
			while (plnode != NULL && vis[plnode->id] == NULL)
				plnode = plnode->father;
			if (plnode == NULL)
				res[i][j].father = NULL;
			else
				res[i][j].father = vis[plnode->id];
		}

	}
	return res;
}

void label_nodes(label_prob_set *prob, std::vector<label_node> *lnodes, int t_start, int t_end) {

	bool tmp = false;
	int n = t_end - t_start;
	int k = prob->k;

	int i, j;
	if (t_end - t_start == 1) {
		for (int i = 0; i < k; i++)
			lnodes[i][0].label = 1;
		return ;
	}


	int mid = (t_start+t_end)/2;
	std::vector<label_node> *llnodes = build_sublist(prob, lnodes, t_start, mid);
	std::vector<label_node> *rlnodes = build_sublist(prob, lnodes, mid, t_end);
	label_nodes(prob, llnodes, t_start, mid);
	label_nodes(prob, rlnodes, mid, t_end);


	// compute lval, rval
	label_node *node;
	for (i = 0; i < k; i++) {
		for (j = 0; j < llnodes[i].size(); j++)
			llnodes[i][j].prev_node->lval = llnodes[i][j].label;
		for (j = 0; j < rlnodes[i].size(); j++)
			rlnodes[i][j].prev_node->rval = rlnodes[i][j].label;
		for (j = 0; j < lnodes[i].size(); j++) {
			node = &lnodes[i][j];
			if (node->lval != 0)
				while (node->father != NULL && node->father->lval == 0) {
					node->father->lval = node->lval;
					node = node->father;
				}
			node = &lnodes[i][j];
			if (node->rval != 0)
				while (node->father != NULL && node->father->rval == 0) {
					node->father->rval = node->rval;
					node = node->father;
				}
		}
	} 

	// radix sort based on lval
	radix_t *radix = prob->radix;
	radix->clear();
	for (i = 0; i < k; i++)
		for (j = 0; j < lnodes[i].size(); j++)
			radix->add(lnodes[i][j].lval, &lnodes[i][j]);
	radix->sort();

	// radix sort based on rval
	int radix_n = radix->n;
	radix->clear();
	for (i = 0; i < radix_n; i++) {
		node = (label_node*)radix->out[i]->val;
		radix->add(node->rval, node);
	}
	radix->sort();

	label_node *prev = NULL;
	int cnt = 0;
	for (i = 0; i < radix->n; i++) {
		node = (label_node*)radix->out[i]->val;
		if (prev == NULL || prev->lval != node->lval 
			|| prev->rval != node->rval)
			cnt++;
		node->label = cnt;
		prev = node;
	}
}

void calc_w_knlogn(std::vector<Tree*>& trees) {
	size_t n = Tree::get_taxas_num();
	size_t k = trees.size();
	int* count = new int[2*k*n];

	std::fill(count, count+2*k*n, 0);
	label_prob_set *prob = new label_prob_set(trees);
	label_nodes(prob, prob->inp, 0, n);

	int i, j;
	for (i = 0; i < trees.size(); i++)
		for (j = 0; j < trees[i]->get_nodes_num(); j++)
			count[prob->inp[i][j].label]++;
	for (i = 0; i < trees.size(); i++)
		for (j = 0; j < trees[i]->get_nodes_num(); j++) {
			Tree::Node* node = trees[i]->get_node(j);
			node->weight = count[prob->inp[i][j].label];
		}
}














// for node_id in T1:
// start[node_id] = smallest rank (i.e. leftmost leaf) in T2 among \Lambda(T1[node_id])
// stop[node_id] = biggest rank (i.e. rightmost leaf) in T2 among \Lambda(T1[node_id])
void compute_start_stop(Tree* tree1, Tree* tree2, int* t2_leaves_ranks) {
	// for each cluster in t1, find its leftmost and rightmost leaves in t2
	for (int i = tree1->get_nodes_num()-1; i >= 0; i--) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf()) {
			start[i] = stop[i] = t2_leaves_ranks[node->taxa];
		} else {
			start[i] = INT32_MAX;
			stop[i] = 0;
			for (Tree::Node* child : node->children) {
				if (start[i] > start[child->id])
					start[i] = start[child->id];
				if (stop[i] < stop[child->id])
					stop[i] = stop[child->id];
			}
		}
	}
}


subpath_query_info_t* preprocess_subpaths_queries(Tree* tree) {
	subpath_query_info_t* subpath_query_info = new subpath_query_info_t(Tree::get_taxas_num()*2);

	// pointer from each node to the root of its cp
	Tree::Node** cp_roots = subpath_query_info->cp_roots;
	cp_roots[0] = tree->get_root();
	for (size_t i = 1; i < tree->get_nodes_num(); i++) {
		Tree::Node* node = tree->get_node(i);
		if (node->pos_in_parent == 0) {
			cp_roots[i] = cp_roots[node->parent->id];
		} else {
			cp_roots[i] = node;
		}
	}

	gen_rmq_t** cp_rmqs = subpath_query_info->cp_rmqs;
	gen_rmq_t* nullp = NULL;
	std::fill(cp_rmqs, cp_rmqs+Tree::get_taxas_num()*2, nullp);
	for (size_t i = 0; i < tree->get_nodes_num(); i++) {
		if (cp_rmqs[cp_roots[i]->id] == NULL) {
			cp_rmqs[cp_roots[i]->id] = new gen_rmq_t;
		}
		cp_rmqs[cp_roots[i]->id]->v.push_back(-tree->get_node(i)->weight);
	}
	for (size_t i = 0; i < tree->get_nodes_num(); i++) {
		if (cp_rmqs[i] != NULL && cp_rmqs[i]->v.size() > 1) {
			general_rmq_preprocess(cp_rmqs[i]);
		}
	}

	int* depths = subpath_query_info->depths;
	depths[0] = 0; // root depth
	for (size_t i = 1; i < tree->get_nodes_num(); i++) {
		depths[i] = 1 + depths[tree->get_node(i)->parent->id];
	}

	return subpath_query_info;
}

// return maximum weight of path from ancestor to descendant
// the node of ancestor and descendant should belong to base tree of subpq_info 
int max_subpath_query(subpath_query_info_t* subpq_info, Tree::Node* ancestor, Tree::Node* descendant) {
	Tree::Node* curr = descendant;
	int res = 0;
	while (subpq_info->cp_roots[curr->id]->id != subpq_info->cp_roots[ancestor->id]->id) {
		gen_rmq_t* currpath_rmq = subpq_info->cp_rmqs[subpq_info->cp_roots[curr->id]->id];
		int query_endp = subpq_info->depths[curr->id] - subpq_info->depths[subpq_info->cp_roots[curr->id]->id];
		if (query_endp == 0) {
			res = std::min(res, currpath_rmq->v[0]);
		} else {
			res = std::min(res, currpath_rmq->v[general_rmq(currpath_rmq, 0, query_endp)]);
		}
		curr = subpq_info->cp_roots[curr->id]->parent;
	}
	gen_rmq_t* currpath_rmq = subpq_info->cp_rmqs[subpq_info->cp_roots[curr->id]->id];
	int query_startp = subpq_info->depths[ancestor->id] - subpq_info->depths[subpq_info->cp_roots[curr->id]->id];
	int query_endp = subpq_info->depths[curr->id] - subpq_info->depths[subpq_info->cp_roots[curr->id]->id];
	if (query_startp < query_endp) {
		res = std::min(res, currpath_rmq->v[general_rmq(currpath_rmq, query_startp+1, query_endp)]);
	}
	return -res;
}

























Tree* contract_tree_fast(Tree *tree, lca_t *lcas, std::vector<int>& marked, subpath_query_info_t *subpq) {
	if (marked.empty()) return NULL;

	// count: contracted node size (with rebundant)
	int count = 0;

	// record leaves and lca between leaves in dfs order
	// levels, ids: depth and node id of contracted nodes
	levels[count] = tree->get_leaf(marked[0])->depth;
	ids[count] = tree->get_leaf(marked[0])->id;
	count++;
	for (size_t i = 1; i < marked.size(); i++) {
		int lca_id = lca(lcas, tree->get_leaf(marked[i-1])->id, tree->get_leaf(marked[i])->id);
		levels[count] = tree->get_node(lca_id)->depth;
		ids[count] = tree->get_node(lca_id)->id;
		count++;
		levels[count] = tree->get_leaf(marked[i])->depth;
		ids[count] = tree->get_leaf(marked[i])->id;
		count++;
	}


	// compute vleft: previous node with smaller depth in dfs order
	std::stack<int> Sl;
	for (int i = 0; i < count; i++) {
		while (!Sl.empty() && levels[i] <= levels[Sl.top()]) {
			Sl.pop();
		}

		if (!Sl.empty()) {
			vleft[i] = Sl.top();
		} else {
			vleft[i] = -1;
		}
		Sl.push(i);
	}

	// compute vright: next node with smaller depth in dfs order
	std::stack<int> Sr;
	for (int i = count-1; i >= 0; i--) {
		while (!Sr.empty() && levels[i] <= levels[Sr.top()]) {
			Sr.pop();
		}

		if (!Sr.empty()) {
			vright[i] = Sr.top();
		} else {
			vright[i] = -1;
		}
		Sr.push(i);
	}

	// if (g_tmp) cout << "cok compute vleft & vright\n";


	// compute parent according to vleft, vright
	// clear existing rebundant node
	int root_pos = -1;
	std::fill(exists, exists+count, true);
	for (int i = 0; i < count; i++) pointer[i] = i;
	for (int i = 0; i < count; i++) {
		parent[i] = -1;
		if (vleft[i] == -1 && vright[i] == -1) { //root
			if (root_pos == -1) {
				root_pos = i;
			}
		} else if (vleft[i] == -1) {
			parent[i] = vright[i];
		} else if (vright[i] == -1) {
			parent[i] = pointer[vleft[i]];
		} else {
			if (levels[vleft[i]] >= levels[vright[i]]) {
				parent[i] = pointer[vleft[i]];
				if (levels[vleft[i]] == levels[vright[i]]) { // deals with non-binary nodes
					pointer[vright[i]] = pointer[vleft[i]];
					exists[vright[i]] = false;
				}
			} else {
				parent[i] = pointer[vright[i]];
			}
		}
	}
	// if (g_tmp) cout << "cok compute parent, pointer\n";



	// construct new_tree
	// tree_nodes: pointer of node from new_tree via dfs index
	Tree* new_tree = new Tree(count*2);
	Tree::Node* nullp = NULL;
	std::fill(tree_nodes, tree_nodes+count, nullp);
	for (int i = 0; i < count; i++) {
		// create node
		if (i%2 == 0 || exists[i]) {
			if (tree_nodes[i] == NULL) {
				if (i%2 == 0) {
					tree_nodes[i] = new_tree->add_node(marked[i/2]);
				} else {
					tree_nodes[i] = new_tree->add_node();
				}
			}
			tree_nodes[i]->weight = tree->get_node(ids[i])->weight;
			tree_nodes[i]->orig_w = tree->get_node(ids[i])->orig_w;
			tree_nodes[i]->secondary_id = ids[i];
		}

		// connect to parent
		if (tree_nodes[i] != NULL && parent[i] != -1) {
			if (tree_nodes[parent[i]] == NULL) {
				tree_nodes[parent[i]] = new_tree->add_node();
			}
			tree_nodes[parent[i]]->add_child(tree_nodes[i]);
		}
	}
	// if (g_tmp) cout << "cok build tree\n";



	// insert special nodes
	if (subpq != NULL) {
		size_t newtree_nodes = new_tree->get_nodes_num();
		for (size_t i = 0; i < newtree_nodes; i++) {
			Tree::Node* curr_node = new_tree->get_node(i);

			if (curr_node->is_root()) continue;

			Tree::Node* desc_par = tree->get_node(curr_node->secondary_id)->parent;
			Tree::Node* anc_par = tree->get_node(curr_node->parent->secondary_id);
			int msq = max_subpath_query(subpq, anc_par, desc_par);
			if (msq > 0) {
				Tree::Node* sp_node = new_tree->add_node();
				sp_node->weight = sp_node->orig_w = msq;
				sp_node->secondary_id = anc_par->id;
				curr_node->parent->null_child(curr_node->pos_in_parent);
				curr_node->parent->add_child(sp_node);
				sp_node->add_child(curr_node);
			}
		}
	}

	new_tree->fix_tree(tree_nodes[root_pos]);
	for (int i = 0; i < new_tree->get_nodes_num(); i++) {
		Tree::Node *node = new_tree->get_node(i);
		Tree::Node *fnode = tree->get_node(node->secondary_id);
		if (fnode->spoil || fnode->leaf_size > node->leaf_size)
			node->spoil = 1;
	}
	// if (g_tmp) cout << "out contract\n";
	return new_tree;
}





















int *minv, *maxv, *fillv, *str_n;
int *ancest;
int *origw_to_w, *subtree_cnt;
struct Interval {
	int l, r;
} static *intervals;
int *inext, *iid;
radix_t *radix;



bool fake_data = false;

inline int unionset_find(int *ancest, int id) {
	static int res, next;
	res = id;
	while (ancest[res] != res)
		res = ancest[res];

	while (id != res) {
		next = ancest[id];
		ancest[id] = res;
		id = next;
	}

	return res;
}

void filter_clusters_nlogn(prob_set *prob, Tree::Node* t1_root, Tree* tree2, bool* to_del) {

	Tree *tree1 = prob->tree1;
	taxas_ranges_t *t1_tr = prob->t1_tr;

	int taxa_st = t1_tr->intervals[t1_root->id].start,
		taxa_ed = t1_tr->intervals[t1_root->id].end;
	int t1_st = t1_root->id,
		t1_ed = t1_st + t1_root->node_size;

	bool tmp = false;

	if (tmp) cout << "in " << t1_root->id << endl;


	// construct heavy path of this subtree 
	std::vector<Tree::Node*> path;
	Tree::Node *node = t1_root;
	while (!node->is_leaf()) node = node->children[0];
	while (true) {
		path.push_back(node);
		if (node == t1_root) break;
		node = node->parent;
	}





	// n: leaf size; m: side root size; pn: path size; 
	// leaf p index (taxa id => side root id): map leaves to corresponding sidetree
	// use in compress subtree weight, but not use in uncompress

	int n = tree2->get_leaves_num();
	int m = 0, pn = path.size();
	std::vector<Tree::Node*> s1_roots, s2_roots;

	// compute all side root and leaf_p_index
	// pnode->tree id: used in weight assignment during weight compress
	leaf_p_index[path[0]->taxa] = -1;
	Tree::Node *pnode;
	int pi, i, j, nid;
	for (pi = 0; pi < pn; pi++) {
		pnode = path[pi];
		pnode->tree_id = -1;
		for (i = 1; i < pnode->children.size(); i++) {
			node = pnode->children[i];
			nid = node->id;
			for (j = t1_tr->intervals[nid].start; j <= t1_tr->intervals[nid].end; j++) {
				leaf_p_index[t1_tr->taxas[j]] = m;
			} 

			m++;
			s1_roots.push_back(node);
		}
	}

	if (tmp) cout << "ok side tree root\n";

	// store weight of tree1
	// build orig_w to w mapping for tree2 and sub_t2
	int *layer_t1_w = new int[t1_ed-t1_st];
	for (i = t1_st; i < t1_ed; i++)
		layer_t1_w[i-t1_st] = tree1->get_node(i)->weight;

	for (i = 0; i < tree2->get_nodes_num(); i++) {
		origw_to_w[tree2->get_node(i)->orig_w] = tree2->get_node(i)->weight;
	}

	// sub_t2: store contracted trees
	Tree **sub_t2 = new Tree*[m];
	std::vector<int> *marked = new std::vector<int>[m];

	// marked (side root id => leaf id set): store leaves in t1 subtrees in t2 left to right order
	// use in contract
	taxas_ranges_t *t2_tr = build_taxas_ranges(tree2);
	for (i = 0; i < n; i++)
		if (leaf_p_index[t2_tr->taxas[i]] > -1)
			marked[leaf_p_index[t2_tr->taxas[i]]].push_back(t2_tr->taxas[i]);

	// TODO: find out the reason that why different leave order work
	if (check_ans && false) {
		cout << "----check ans----\n";
		taxas_ranges_t *orig_tr = prob->t2_tr;
		int orig_n = orig_tr->taxas_num;
		int *subseq = new int[n];
		bool error = false;
		for (int i = 0; i < n; i++)
			cout << t2_tr->taxas[i] << ' '; cout << endl;
		
		int taxa, j = 0;
		for (int i = 0; i < orig_n; i++) {
			taxa = orig_tr->taxas[i];
			// cout << t1_tr->intervals[tree1->get_leaf(taxa)->id].start << ' ' << taxa_st << ' ' << taxa_ed << endl;
			if (t1_tr->intervals[tree1->get_leaf(taxa)->id].start >= taxa_st && t1_tr->intervals[tree1->get_leaf(taxa)->id].start <= taxa_ed)
				subseq[j++] = taxa;
		}
		for (int i = 0; i < n; i++) {
			cout << subseq[i] << ' ';
			if (subseq[i] != t2_tr->taxas[i]) error = true;
		} cout << endl;

		if (error) {
			cout << "ERROR sub_t2 taxas order" << endl;
			exit(1);
		}
	}
	if (tmp) cout << "ok marked\n";




	// contract tree + radix sort
	Tree::Node *snode;
	radix->clear();
	int cnt1 = 0, cnt2 = 0;
	for (i = 0; i < m; i++) {
		node = s1_roots[i];
		sub_t2[i] = contract_tree_fast(prob->tree2, prob->t2_lcas, marked[i], prob->subpq_t2);

		nid = node->id;

		for (j = nid; j < nid + node->node_size; j++) {
			// cout << "t1 " << j << endl;
			snode = tree1->get_node(j);
			snode->tree_id = i;
			radix->add(snode->weight, snode);
		}
		for (j = 0; j < sub_t2[i]->get_nodes_num(); j++) {
			snode = sub_t2[i]->get_node(j);
			snode->weight = origw_to_w[snode->orig_w];
			snode->tree_id = i;
			radix->add(snode->weight, snode);
		}
	}
	radix->sort();
	// exit(0);

	if (tmp) cout << "ok contract\n";


	// compress weight in sub_t1 and sub_t2
	std::fill(subtree_cnt, subtree_cnt+m, 0);
	Tree::Node *prev = NULL;
	for (int i = 0; i < radix->n; i++) {
		node = (Tree::Node*)radix->out[i]->val;
		if (prev == NULL || prev->tree_id != node->tree_id || prev->orig_w != node->orig_w)
			subtree_cnt[node->tree_id] ++;
		node->weight = subtree_cnt[node->tree_id];
		prev = node;
	}

	if (tmp) cout << "ok compress weight" << endl;




	// solve recursively
	for (i = 0; i < m; i++) {
		filter_clusters_nlogn(prob, s1_roots[i], sub_t2[i], to_del);
	}

	// recover weight
	for (int i = t1_st; i < t1_ed; i++)
		tree1->get_node(i)->weight = layer_t1_w[i-t1_st];

	if (tmp) cout << "mid " << t1_root->id << endl;
	// return ;




	// work on nodes on centroid path and tree2
	// calculate stratify length

	// str_n (cpath node id -> rel taxa id): store (relative) t1_tr taxa order of endpoint of path node
	// rel taxa id: regard first taxa of centroid path leave as 0
	// cout << pn << ' ' << path.size() << endl;
	for (i = 0; i < pn; i++) {
		str_n[i] = t1_tr->intervals[path[i]->id].end - taxa_st;
	} 
	// return ;



	// minv, maxv, fillv (node id -> rel taxa id): store corresponding value of node in tree2
	// intervals: simply [min(minv, fillv), maxv]
	// ancest: union set of rel taxa id, use to compute fillv
	// iraidx: sort interval in weight, use to solve skyline problem

	radix->clear();

	int sid, anc_id;
	for (i = 0; i < n; i++)
		ancest[i] = i;
	for (i = tree2->get_nodes_num()-1; i >= 0; i--) {
		node = tree2->get_node(i);

		// calculate minv, maxv
		if (node->is_leaf()) {
			minv[i] = maxv[i] = t1_tr->intervals[tree1->get_leaf(node->taxa)->id].start - taxa_st;
			fillv[i] = -1;
		} else {
			for (j = 0; j < node->children.size(); j++) {
				sid = node->children[j]->id;
				if (j == 0)
					minv[i] = minv[sid], maxv[i] = maxv[sid], fillv[i] = fillv[sid];
				else {
					minv[i] = std::min(minv[i], minv[sid]);
					maxv[i] = std::max(maxv[i], maxv[sid]);
					fillv[i] = std::max(fillv[i], fillv[sid]);
				}
			}

			// union the taxas for fillv
			anc_id = unionset_find(ancest, minv[node->children[0]->id]);
			for (j = 1; j < node->children.size(); j++) {
				sid = node->children[j]->id;
				ancest[unionset_find(ancest, minv[sid])] = anc_id;
			}
		}

		// calculate fillv
		anc_id = unionset_find(ancest, minv[i]);
		while (fillv[i] < n-1 && unionset_find(ancest, fillv[i]+1) == anc_id)
			fillv[i]++;

		intervals[i].l = std::max(minv[i], fillv[i]+1);
		intervals[i].r = (node->spoil ? n-1 : maxv[i]-1);
		if (intervals[i].l <= intervals[i].r) {
			radix->add(node->weight, intervals+i);
			// cout << node->weight << ' ' << intervals[i].l << ' ' << intervals[i].r << endl;
		}
	}
	radix->sort(true);

	if (tmp) cout << "ok interval\n";
	// return ;



	// solve max-manhattan skyline problem
	int ti = 0;
	for (i = 0; i < n; i++) {
		inext[i] = str_n[ti];
		if (i == str_n[ti]) {
			iid[i] = path[ti]->id;
			ti++;
		} else iid[i] = -1;
	}
	if (tmp) cout << "ok iid\n";

	// return ;

	static Interval *interval;

	inext[n] = n;
	for (i = 0; i < radix->n; i++) {
		interval = (Interval*)radix->out[i]->val;
		// cout << radix->out[i]->key << ' ' << interval->l << ' ' << interval->r << endl;
		j = unionset_find(inext, interval->l);
		while (j <= interval->r) {
			if (iid[j] >= 0 && radix->out[i]->key >= tree1->get_node(iid[j])->weight) {
				to_del[iid[j]] = true;
			}
			inext[j] = unionset_find(inext, j+1);
			j = inext[j];
		}
	}
	if (tmp) cout << "ok to_del\n";

	if (tmp) cout << "out " << t1_root->id << endl;
}














// TODO: change n_to_w and w_to_n to int[k]
// TODO: use nlogn sort other than radix sort if k > n
void filter(Tree* tree1, Tree* tree2, 
		bool* to_del) {
	auto st = Now();

	// backup node weights as orig_w
	Tree::Node *node;
	int i;
	for (i = 0; i < tree1->get_nodes_num(); i++)
		tree1->get_node(i)->orig_w = tree1->get_node(i)->weight;
	for (i = 0; i < tree2->get_nodes_num(); i++)
		tree2->get_node(i)->orig_w = tree2->get_node(i)->weight;

	prob_set *prob = new prob_set(tree1, tree2);

	// compress node weights in nlogn (in case k >> n)
	if (radix->k > radix->n) {
		radix->clear();
		for (i = 0; i < tree1->get_nodes_num(); i++)
			radix->add(tree1->get_node(i)->weight, tree1->get_node(i));
		for (i = 0; i < tree2->get_nodes_num(); i++)
			radix->add(tree2->get_node(i)->weight, tree2->get_node(i));
		radix->quicksort();

		int cnt = 0;
		for (i = 0; i < radix->n; i++) {
			node = (Tree::Node*)radix->out[i]->val;
			if (i == 0 || radix->out[i-1]->key != radix->out[i]->key)
				cnt++;
			node->weight = cnt;
		}
	}

	// start filter clusters
	filter_clusters_nlogn(prob, tree1->get_root(), tree2, to_del);

	// recover node weights
	for (i = 0; i < tree1->get_nodes_num(); i++)
		tree1->get_node(i)->weight = tree1->get_node(i)->orig_w;
	for (i = 0; i < tree2->get_nodes_num(); i++)
		tree2->get_node(i)->weight = tree2->get_node(i)->orig_w;
	auto ed = Now();
	if (g_tmp) cout << duration_sec(st, ed) << endl;
}






































void compute_m(Tree::Node* node, int* e, int* m, std::vector<Tree::Node*>* rsort_lists) {
	if (node->is_leaf()) {
		m[node->id] = e[node->taxa];
	}
	for (Tree::Node* child : node->children) {
		compute_m(child, e, m, rsort_lists);
		if (m[node->id] > m[child->id]) {
			m[node->id] = m[child->id];
		}
	}
	if (!node->is_root()) {
		rsort_lists[m[node->id]].push_back(node);
	}
}


// See Section 2.4 of paper [XXX]
void merge_trees(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, lca_t* t2_lcas) {
	//compute e
	for (size_t j = 0; j < Tree::get_taxas_num(); j++) {
		e[t1_tr->taxas[j]] = j;
	}
	std::fill(m, m+tree2->get_nodes_num(), INT32_MAX);
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) rsort_lists[i].clear();

	compute_m(tree2->get_root(), e, m, rsort_lists);

	// sort tree2
	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		tree2->get_node(i)->clear_children();
	}
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		for (auto it = rsort_lists[i].begin(); it != rsort_lists[i].end(); it++) {
			(*it)->parent->add_child(*it);
		}
	}

	taxas_ranges_t* t2_tr = build_taxas_ranges(tree2);
	int* t2_ranks = get_taxas_ranks(t2_tr);
	compute_start_stop(tree1, tree2, t2_ranks);
	delete[] t2_ranks;

	// calc x_left and x_right
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		Tree::Node* curr = tree2->get_leaf(i);
		Tree::Node* parent = curr->parent;
		while (parent != NULL && *(parent->children.begin()) == curr) {
			curr = parent;
			parent = curr->parent;
		}
		_left[i] = curr;

		curr = tree2->get_leaf(i);
		parent = curr->parent;
		while (parent != NULL && *(parent->children.rbegin()) == curr) {
			curr = parent;
			parent = curr->parent;
		}
		_right[i] = curr;
	}

	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		orig_pos_in_parent[i] = tree2->get_node(i)->pos_in_parent;
	}

	for (int i = tree1->get_nodes_num()-1; i >= 1; i--) {
		Tree::Node* a = tree2->get_leaf(t2_tr->taxas[start[i]]);
		Tree::Node* b = tree2->get_leaf(t2_tr->taxas[stop[i]]);
		if (a == b) continue; // tree1->node i is a leaf

		Tree::Node* ru = tree2->get_node(lca(t2_lcas, a->id, b->id));

		Tree::Node* a_left = _left[a->taxa];
		Tree::Node* b_right = _right[b->taxa];

		size_t du_pos = (a_left->depth > ru->depth) ? orig_pos_in_parent[a_left->id] : 0;
		size_t eu_pos = (b_right->depth > ru->depth) ? orig_pos_in_parent[b_right->id] : ru->get_children_num()-1;
		if (du_pos == 0 && eu_pos == ru->get_children_num()-1) continue;

		Tree::Node* newnode = tree2->add_node();
		newnode->weight = tree1->get_node(i)->weight;
		newnode->leaf_size = tree1->get_node(i)->leaf_size;
		for (size_t j = du_pos; j <= eu_pos; j++) {
			if (ru->children[j] != NULL) {
				newnode->add_child(ru->children[j]);
				// lazy child deletion
				ru->null_child(j);
			}
		}
		ru->set_child(newnode, du_pos);
	}
	tree2->fix_tree();

	delete t2_tr;
}






























Tree* freqdiff(std::vector<Tree*>& trees) {

	// freopen("oup.txt", "w", stdout);

	start = new int[Tree::get_taxas_num()*2];
	stop = new int[Tree::get_taxas_num()*2];
	e = new int[Tree::get_taxas_num()];
	m = new int[Tree::get_taxas_num()*2];
	rsort_lists = new std::vector<Tree::Node*>[Tree::get_taxas_num()];

	marked = new bool[Tree::get_taxas_num()*2];
	std::fill(marked, marked+Tree::get_taxas_num()*2, false);
	bool* to_del_t = new bool[Tree::get_taxas_num()*2];
	bool* to_del_ti = new bool[Tree::get_taxas_num()*2];

	_left = new Tree::Node*[Tree::get_taxas_num()];
	_right = new Tree::Node*[Tree::get_taxas_num()];
	orig_pos_in_parent = new size_t[Tree::get_taxas_num()*2];

	BT = new bool[Tree::get_taxas_num()*2];
	counter = new size_t[Tree::get_taxas_num()*2];
	leaf_p_index = new int[Tree::get_taxas_num()];

	// tree contract
	vleft = new int[Tree::get_taxas_num()*2];
	vright = new int[Tree::get_taxas_num()*2];
	pointer = new int[Tree::get_taxas_num()*2];
	levels = new int[Tree::get_taxas_num()*2];
	ids = new int[Tree::get_taxas_num()*2];
	parent = new int[Tree::get_taxas_num()*2];
	exists = new bool[Tree::get_taxas_num()*2];
	tree_nodes = new Tree::Node*[Tree::get_taxas_num()*2];


	int n = Tree::get_taxas_num();
	int k = trees.size();
	minv = new int[3*n]; maxv = new int[3*n]; fillv = new int[3*n];
	ancest = new int[n]; inext = new int[2*n]; iid = new int[2*n];
	str_n = new int[2*n];
	intervals = new Interval[3*n];
	radix = new radix_t(5*n, k);
	origw_to_w = new int[k+1];
	subtree_cnt = new int[n*2];









	// weights[i][node id] = weights of cluster of node (with node id) in tree i
	g_tmp = false;
	auto st = Now();

	if (g_tmp) cout << "start\n";
	calc_w_knlogn(trees);
	if (g_tmp) cout << "ok weight" << endl;

	auto ed = Now();
	cout << duration_sec(st, ed) << endl;

	// freopen("oup.txt", "w", stdout);
	// for (int i = 0; i < k; i++) {
	// 	for (int j = 0; j < trees[i]->get_nodes_num(); j++)
	// 		cout << trees[i]->get_node(j)->weight << ' ';
	// 	cout << endl;
	// }

	// exit(0);



	// for (int i = 0; i < trees.size(); i++) {
	// 	for (int j = 0; j < trees[i]->get_nodes_num(); j++) {
	// 		cout << trees[i]->get_node(j)->weight << ' ';
	// 	} cout << '\n'; }


	// if (fake_data) {
	// 	for (int i = 0; i < trees.size(); i++) {
	// 		for (int j = 0; j < trees[i]->get_nodes_num(); j++)
	// 			trees[i]->get_node(j)->weight = (7*i+13*j)%97+1;
	// 	}
	// }







	for (size_t i = 0; i < trees.size(); i++)
		trees[i]->reorder();

	// subpath_query_info_t** subpath_query_ti = NULL;
	// subpath_query_ti = new subpath_query_info_t*[trees.size()];
	// for (size_t i = 0; i < trees.size(); i++) {
	// 	subpath_query_ti[i] = preprocess_subpaths_queries(trees[i]);
	// }

	Tree* T = new Tree(trees[0]);
	for (size_t i = 1; i < trees.size(); i++) {
		if (g_tmp) cout << i << endl;

		Tree* Ti = new Tree(trees[i]);


		taxas_ranges_t* tr_Ti = build_taxas_ranges(Ti);

		taxas_ranges_t* tr_T = build_taxas_ranges(T);

		lca_t* lca_T = lca_preprocess(T);

		if (g_tmp) cout << "ok prepare\n";




		// if (check_ans) {


		// 	std::fill(to_del_ti, to_del_ti+Ti->get_nodes_num(), false);
		// 	filter(Ti, T, to_del_ti);

		// 	freopen("oup.txt", "w", stdout);
		// 	for (int j = 0; j < Ti->get_nodes_num(); j++)
		// 		cout << to_del_ti[j] << ' '; cout << endl;
		// 	exit(0);
		// }



		// filter clusters
		std::fill(to_del_ti, to_del_ti+Ti->get_nodes_num(), false);
		filter(Ti, T, to_del_ti);
		if (g_tmp) cout << "ok del_ti 1" << endl;




		std::fill(to_del_t, to_del_t+T->get_nodes_num(), false);
		filter(T, Ti, to_del_t);
		if (g_tmp) cout << "ok del_t 2" << endl;

		// cout << Ti->get_nodes_num() << ' ' << T->get_nodes_num() << endl;

		// freopen("oup.txt", "w", stdout);
		// for (int i = 0; i < Ti->get_nodes_num(); i++)
		// 	cout << to_del_ti[i] << ' '; cout << endl;
		// for (int i = 0; i < T->get_nodes_num(); i++)
		// 	cout << to_del_t[i] << ' '; cout << endl;


		// exit(0);

		Ti->delete_nodes(to_del_ti);
		T->delete_nodes(to_del_t);
		if (g_tmp) cout << "ok delete nodes" << endl;

		delete lca_T;
		delete tr_T;
		delete tr_Ti;

		lca_T = lca_preprocess(T);
		tr_Ti = build_taxas_ranges(Ti);

		merge_trees(Ti, T, tr_Ti, lca_T);
		if (g_tmp) cout << "ok merge" << endl;

		delete lca_T;
		delete tr_Ti;

		// reset secondary ids
		for (size_t j = 0; j < T->get_nodes_num(); j++) {
			T->get_node(j)->secondary_id = T->get_node(j)->id;
		}

		delete Ti;
	}

	std::fill(to_del_t, to_del_t+T->get_nodes_num(), false);
	taxas_ranges_t* tr_T = build_taxas_ranges(T);
	for (size_t i = 0; i < trees.size(); i++) {
		filter(T, trees[i], to_del_t);
		delete trees[i];
		if (g_tmp) cout << "ok final delete tree " << i << endl;
	}
	T->delete_nodes(to_del_t);
	delete tr_T;

	delete[] start;
	delete[] stop;
	delete[] e;
	delete[] m;
	delete[] rsort_lists;

	delete[] marked;
	delete[] to_del_ti;
	delete[] to_del_t;

	delete[] orig_pos_in_parent;
	delete[] _left;
	delete[] _right;


	return T;
}


#endif /* FREQDIFF_H_ */
