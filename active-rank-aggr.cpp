
#include <iomanip>
#include <iostream>

#include <string.h>
#include <time.h>
#include <cmath>

#include <sstream>

#include <algorithm>

#include <queue>
#include <vector>
#include <map>
	
#include "kutils.hpp"
#include "kalgorithms.hpp"
#include "kstructure.h"

using namespace std;
using namespace klib;

typedef int element_type;
typedef int rank_type;
typedef pair<double, Permutation<element_type, rank_type> *> result;

int main(int argc, char const *argv[])
{
	typedef pair<vector<Permutation<element_type, rank_type>*>, std::map<string, int> > query_set_type;
	
	string qresult_text = argv[1];//"../../../Machine Learning/Tp Machine Learning/results/test Marcos Andre/5fold-change/qresults_50_1";
	string qresult_image = argv[2];//"../../../Machine Learning/Tp Machine Learning/results/test Marcos Andre/resultado/result_50_1.txt";
	
	cout << qresult_text << endl;

	query_set_type aggrs;
	query_set_type text = load_format_trec<element_type, rank_type, PermutationTree<element_type, rank_type> >(qresult_text);
	query_set_type img = load_format_trec<element_type, rank_type, PermutationTree<element_type, rank_type> >(qresult_image);

	std::map<string, int>::iterator i, t;
	int idx;
	for (i = img.second.begin(), t = text.second.begin(), idx = 0; i != img.second.end() && t != text.second.end(); ++i, ++t, idx++)
	{
		cout << i->first << " : " << t->first << endl;
		vector<Permutation<element_type, rank_type>*> perms;
		perms.push_back(img.first[i->second]);
		perms.push_back(text.first[t->second]);
		Permutation<element_type, rank_type> * sigma = new PermutationTree<element_type, rank_type>();
		union_all(*sigma, perms);

		IGraph * tournament = create_majority_graph(perms, perms.size(), sigma->size(), sigma->size(), *sigma);
		
		//result mix = mixed<element_type, rank_type>(*tournament, perms, *sigma);
		//aggrs.first.push_back(mix.second);

		aggrs.first.push_back(kFAS_pivot<element_type, rank_type>(*tournament, *sigma, true));
		aggrs.second[t->first] = idx;

		delete tournament;
		delete sigma;
	}

	string out_file = argv[3];//"../../../Machine Learning/Tp Machine Learning/results/test Marcos Andre/resultado/qresults-aggr_50_1.txt";
	save_format_trec<element_type, rank_type>(out_file, aggrs);

	return EXIT_SUCCESS;
}