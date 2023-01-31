#ifndef NEX_PARSER_H_
#define NEX_PARSER_H_

#include <istream>
#include <iostream>
#include <cassert>

#include "Tree.h"

using std::cout;
using std::endl;

std::vector<Tree*> parse_nex_short(std::istream& in) {
	std::string s;
	// Parse trees
	std::vector<Tree*> trees;
	while (getline(in, s)) {
		// cout << s << endl;
		if (s.length() <= 2)
			break;
		trees.push_back(new Tree(s));
	}
	// for (auto tree: trees)
	// 	cout << tree->to_string() << endl;
	return trees;
}

std::vector<Tree*> parse_nex(std::istream& in) {
	std::string s;
	while (getline(in, s)) {
		if (s.find("BEGIN TREES;") != std::string::npos)
			break;
	}
	getline(in, s); // "translate"

	// Parse taxas
	while (getline(in, s)) {
		if (s.find(";") != std::string::npos)
			break;
	}
	// Parse trees
	std::vector<Tree*> trees;
	while (getline(in, s)) {
		if (s.find("END;") != std::string::npos)
			break;
		trees.push_back(new Tree(s));
	}
	return trees;
}


#endif
