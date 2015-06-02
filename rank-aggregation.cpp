
#include <iomanip>
#include <iostream>

#include <string.h>
#include <time.h>
#include <cmath>

#include <algorithm>

#include <queue>
#include <vector>
#include <map>
	
#include "kutils.hpp"
#include "kalgorithms.hpp"
#include "kstructure.h"

using namespace std;
using namespace klib;

#include "Graph_Hamilton_Path.cpp"
#include "Graph_Search.cpp"
#include "Graph_Minimum_Spanning_Tree.cpp"
#include "Single_Source_Shortest_Path.cpp"
#include "Graph_Feedback_Arc_Set.cpp"

char* union_all(map<char, int> ranks[], int n, int &size){
	map<char, int> rank = ranks[0];
	for (int i = 1; i < n; ++i)
	{	
		
		for (map<char, int>::iterator it = ranks[i].begin(); it !=  ranks[i].end(); ++it)
		{
			rank[it->first] = 0;
		}
	}

	int i = 0;
	char *c = new char[rank.size()];
	for (map<char, int>::iterator it = rank.begin(); it !=  rank.end(); ++it)
	{
		c[i] = it->first;
		i++;
	}

	size = i;

	return c;
}

int main(int argc, char const *argv[])
{
	int seed = time(NULL);
	srand(seed);

	cout << "seed : " << seed << endl;
	
	const int k = 10; // number of list
	const int n = 8; // number of candidates
	const int m = n; // top-m rank

	typedef char element_type;
	typedef int rank_type;

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

	cout << "------------------------------------------------------------" << endl;
	cout << "Executing Optimal-Kemeny-Consensus..." << endl;
	cout << "------------------------------------------------------------" << endl;
	IGraph * tournament = create_majority_graph(perms, k, m, n, 'A');
	Permutation<element_type, rank_type> * opt = kemeny_consensus(perms[0], perms);
	
	opt->print();
	cout << "\nSum Kendall rank(optimal) : " << _kemeny_rule<element_type, rank_type>(*opt, perms) << endl;

	Permutation<element_type, rank_type> * opt_kemenized = locally_kemenysation(*tournament, *opt, 'A');
	opt_kemenized->print();
	cout << "\nSum Kendall rank(optimal kemenized) : " << _kemeny_rule<element_type, rank_type>(*opt_kemenized, perms) << endl;
	cout << "------------------------------------------------------------" << endl;
	cout << "Executing Pink-A-List 2-aproximated algorithm..." << endl;
	cout << "------------------------------------------------------------" << endl;
	for (int i = 0; i < k; ++i)
	{
		for (int j = i+1; j < k; ++j)
		{
			double dist = kendall_dist<element_type, rank_type>(*perms[i], *perms[j]);
			sum_dist[i] += dist;
			sum_dist[j] += dist;
		}
	}

	int min_dist = 0;
	for (int i = 1; i < k; ++i)
	{
		if(sum_dist[i] < sum_dist[min_dist]){
			min_dist = i;
		}
	}

	perms[min_dist]->print();
	cout << "\nSum Kendall rank(" << min_dist + 1 << ") : " << sum_dist[min_dist] << endl;

	Permutation<element_type, rank_type> * kemenized = locally_kemenysation(*tournament, *perms[min_dist], 'A');
	kemenized->print();
	cout << "\nSum Kendall rank(" << min_dist + 1 << ") : " << _kemeny_rule<element_type, rank_type>(*kemenized, perms) << endl;

	cout << "------------------------------------------------------------" << endl;
	cout << "Executing FAS-pivot 2-aproximated algorithm..." << endl;
	cout << "------------------------------------------------------------" << endl;
	Permutation<element_type, rank_type> * fas_pivot = kFAS_pivot<element_type, rank_type>(*tournament, 'A');
	fas_pivot->print();
	cout << "\nSum Kendall rank(FAS Heuristic) : " << _kemeny_rule<element_type, rank_type>(*fas_pivot, perms) << endl;

	Permutation<element_type, rank_type> * fas_pivot_kemenized = kFAS_pivot<element_type, rank_type>( *tournament, 'A', true);
	fas_pivot_kemenized->print();
	cout << "Sum Kendall rank(FAS Heuristic) : " << _kemeny_rule<element_type, rank_type>(*fas_pivot_kemenized, perms) << endl;

	cout << "------------------------------------------------------------" << endl;
	cout << "Executing FHP_greedy algorithm..." << endl;
	cout << "------------------------------------------------------------" << endl;
	Permutation<element_type, rank_type> * fhp_greedy = FHP_greedy<element_type, rank_type>(*tournament, 'A');
	fhp_greedy->print();
	cout << "\nSum Kendall rank(FAS Heuristic) : " << _kemeny_rule<element_type, rank_type>(*fhp_greedy, perms) << endl;

	Permutation<element_type, rank_type> * fhp_greedy_kemenized = FHP_greedy<element_type, rank_type>( *tournament, 'A', true);
	fhp_greedy_kemenized->print();
	cout << "Sum Kendall rank(FAS Heuristic) : " << _kemeny_rule<element_type, rank_type>(*fhp_greedy_kemenized, perms) << endl;	

	return EXIT_SUCCESS;
}