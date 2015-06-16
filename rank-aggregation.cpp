
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

typedef char element_type;
typedef int rank_type;
typedef pair<double, Permutation<element_type, rank_type> *> result;

int main(int argc, char const *argv[])
{
	int seed = time(NULL);//1433258641;//time(NULL);//1433249523;
	srand(seed);

	//const int k = 11; // number of list
	//const int n = 9; // number of candidates
	//const int m = n; // top-m rank
	cout << "seed : " << seed << endl;
	for (int n = 3; n < 9; ++n)
	{
		int m = n;	

		for (int k = 10; k < 100; k+=10)
		{
			vector<Permutation<element_type, rank_type>* > perms;
			
			double *sum_dist = new double[k];
			memset(sum_dist, 0, sizeof(double)*k);
			
			cout << "Generating syntetic rankings..." << endl;
			for (int i = 0; i < k; ++i)
			{	
			 	perms.push_back(gen_rank<char, int, PermutationTree<element_type, rank_type> >(1, n, m, 'A'));
			 	perms[i]->print(); cout << endl;
			}

			cout << "# rankings : " << k << endl;
			cout << "# candidates : " << n << endl;
			cout << "# top-m size : " << m << endl;

			Permutation<element_type, rank_type> * sigma = new PermutationTree<element_type, rank_type>();
			union_all(*sigma, perms);
			
			sigma->print(); cout << endl<< endl;

			IGraph * tournament = create_majority_graph(perms, k, m, sigma->size(), *sigma);
			((Graph_Adj_Matrix*)tournament)->printAsMatrix();
			cout << endl;

			stringstream filename;
			
			filename << "rank-results/complete list/results_" << "k" << k << "_n" << n << "_m" << m;
			
			testing_all(*tournament, perms, *sigma, 'A', filename.str());
		}

	}

	/*
	cout << "------------------------------------------------------------" << endl;
	cout << "Executing Optimal-Kemeny-Consensus..." << endl;
	cout << "------------------------------------------------------------" << endl;
	
	Permutation<element_type, rank_type> * opt = kemeny_consensus(sigma, perms);

	opt->print();
	cout << "\nSum Kendall rank(optimal) : " << _kemeny_rule<element_type, rank_type>(*opt, perms) << endl;

	Permutation<element_type, rank_type> * opt_kemenized = locally_kemenysation(*tournament, *opt, 'A');
	opt_kemenized->print();
	cout << "\nSum Kendall rank(optimal kemenized) : " << _kemeny_rule<element_type, rank_type>(*opt_kemenized, perms) << endl;

	cout << "------------------------------------------------------------" << endl;
	cout << "Executing Pink-A-List 2-aproximated algorithm..." << endl;
	cout << "------------------------------------------------------------" << endl;

	Permutation<element_type, rank_type> * pal = pick_a_list<element_type, rank_type>(perms);
	union_perm(*pal, *sigma);
	pal->print();
	cout << "\nSum Kendall rank(pick a list) : " << _kemeny_rule<element_type, rank_type>(*pal, perms) << endl;

	
	pal->print(); cout << "\n\n";
	Permutation<element_type, rank_type> * kemenized = locally_kemenysation(*tournament, *pal, 'A');
	kemenized->print();
	cout << "\nSum Kendall rank(pal kemenized) : " << _kemeny_rule<element_type, rank_type>(*kemenized, perms) << endl;

	cout << "------------------------------------------------------------" << endl;
	cout << "Executing FAS-pivot 2-aproximated algorithm..." << endl;
	cout << "------------------------------------------------------------" << endl;
	Permutation<element_type, rank_type> * fas_pivot = kFAS_pivot<element_type, rank_type>(*tournament, 'A');
	fas_pivot->print();
	cout << "\nSum Kendall rank(FAS Heuristic) : " << _kemeny_rule<element_type, rank_type>(*fas_pivot, perms) << endl;

	Permutation<element_type, rank_type> * fas_pivot_kemenized = kFAS_pivot<element_type, rank_type>( *tournament, 'A', true);
	fas_pivot_kemenized->print();
	cout << "\nSum Kendall rank(FAS Heuristic) : " << _kemeny_rule<element_type, rank_type>(*fas_pivot_kemenized, perms) << endl;

	cout << "------------------------------------------------------------" << endl;
	cout << "Executing FHP_greedy algorithm..." << endl;
	cout << "------------------------------------------------------------" << endl;
	Permutation<element_type, rank_type> * fhp_greedy = FHP_greedy<element_type, rank_type>(*tournament, 'A');
	fhp_greedy->print();
	cout << "\nSum Kendall rank(FHP Heuristic) : " << _kemeny_rule<element_type, rank_type>(*fhp_greedy, perms) << endl;

	Permutation<element_type, rank_type> * fhp_greedy_kemenized = FHP_greedy<element_type, rank_type>( *tournament, 'A', true);
	fhp_greedy_kemenized->print();
	cout << "\nSum Kendall rank(FHP Heuristic) : " << _kemeny_rule<element_type, rank_type>(*fhp_greedy_kemenized, perms) << endl;

	cout << "------------------------------------------------------------" << endl;
	cout << "Executing Mixed algorithm..." << endl;
	cout << "------------------------------------------------------------" << endl;
	result mixed_result = mixed<element_type, rank_type>(*tournament, perms, 'A');
	mixed_result.second->print();
	cout << "\nSum Kendall rank(FHP Heuristic) : " << mixed_result.first << endl;

	((Graph_Adj_Matrix*)tournament)->printAsMatrix();
	*/
	return EXIT_SUCCESS;
}