#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

#include "nex_parser.h"

#include "freqdiff.h"
// #include "wrong.h"
// #include "freqdiff2.h"
// #include "freqdiff2_backup.h"
#include "local_consensus.h"
#include "Tree.cpp"

using namespace std::chrono;

#define LOCAL_TEST

#define Now high_resolution_clock::now
#define duration_sec(start,end) (float)(duration_cast<milliseconds>(end-start).count() / 1000.)

using namespace std;

bool readArg(int argc, char **argv) {
    if (argc == 1)
        return false;
    int t_opt;
    while ((t_opt = getopt(argc, argv, "k") ) != -1)
        if(t_opt == 'k')
            return true;
    return false;
}

int main(int argc, char** argv) {

	#ifdef LOCAL_TEST

	string algo("freq");
	// ifstream fin("input.nex");
	// vector<Tree*> trees = parse_nex(fin);

	ifstream fin("inp.txt");
	vector<Tree*> trees = parse_nex_short(fin);

    bool alg_choice = readArg(argc, argv);

	#else

	if (argc != 3) {
		cout << "Usage: exec [freq|minrs|minis] filename" << endl;
		cout << "More info to come" << endl;
		return 0;
	}

	string algo(argv[1]);
	ifstream fin(argv[2]);
	vector<Tree*> trees = parse_nex(fin);

    bool alg_choice = false;

	#endif



	auto st = Now();

	Tree* consensus;
	// try 
	// {
		if (algo == "freq") {
			// consensus = freqdiff(trees, Tree::get_taxas_num() > 1000);
			consensus = freqdiff(alg_choice, trees);
		} else if (algo == "minrlc_exact") {
			consensus = minRLC_exact(trees);
		} else if (algo == "minilc_exact") {
			consensus = minILC_exact(trees);
		} else if (algo == "aho-build") {
			consensus = ahoBuild(trees);
		} else {
			std::cerr << "Algorithm not supported" << std::endl;
			return 1;
		}
	// }
	// catch ( exception &e ) 
	// {
	// 	cout << "Caught " << e.what( ) << endl;
	// 	cout << "Type " << typeid( e ).name( ) << endl;
	// };

	auto ed = Now();
	cout << duration_sec(st, ed) << endl;

	if (consensus == NULL) {
		std::cout << "No valid consensus found." << std::endl;
	} else {
        alg_choice ? freopen("oup3.txt", "w", stdout) : freopen("oup2.txt", "w", stdout);
		cout << consensus->to_newick() << endl;
		fclose(stdout);
	}

	delete consensus;

	return 0;
}
