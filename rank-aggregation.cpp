
#include <iomanip>
#include <iostream>

#include <string.h>
#include <time.h>
#include <cmath>

#include <algorithm>

#include <queue>
#include <vector>
#include <map>

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

	srand(time(NULL));
	
	int k = 6; // number of list
	int n = 4;  // number of candidates
	int m = n; // top-m rank

	map<char, int> *ranks = new map<char, int>[k];
	vector<Permutation<char, int>* > perms;
	char **ranks_arr = new char*[k];
	map<char, int> opt;
	
	double *sum_dist = new double[k];
	memset(sum_dist, 0, sizeof(double)*k);
	
	for (int i = 0; i < k; ++i)
	{	
		ranks_arr[i] = new char[m];
	 	ranks[i] = gen_rank(1, n, m, ranks_arr[i]);
	 	perms.push_back(gen_rank<char, int, PermutationTree<char, int> >(1, n, m, 'A'));
	 	cout << ranks_arr[i] << endl;

	 	perms[i]->print(); cout << endl;
	 	cout << "VV :"  << (*perms[i])('A') << " : " << (*perms[i])('J') << endl;
	}

	IGraph * tournament = create_majority_graph(perms, k, m, n, 'A');
	kemeny_consensus(perms[0], perms)->print();

    cout << "Gen all..." <<endl;

	for (int i = 0; i < k; ++i)
	{
		for (int j = i+1; j < k; ++j)
		{
			double dist = kendall_dist(ranks[i], ranks[j]);
			sum_dist[i] += dist;
			sum_dist[j] += dist;

			//cout << i << " - " << j << " : " << dist << endl;
		}
	}

	int min_dist = 0;
	for (int i = 1; i < k; ++i)
	{
		//print_rank(ranks[i]);
		//cout << "Sum Kendall rank(" << i + 1 << ") : " << sum_dist[i] << endl;
		if(sum_dist[i] < sum_dist[min_dist]){
			min_dist = i;
		}
	}

	print_rank(ranks[min_dist]);
	cout << "Sum Kendall rank(" << min_dist + 1 << ") : " << sum_dist[min_dist]/k << endl;

	int cl;
	char *c = union_all(ranks, k, cl);
	cout << "size : " << cl << endl;
	print_array(c, cl);
	cout << "OPTIMAL :: " << endl;
	opt = kemeny_consensus(c, cl, ranks, k);
	int * out = new int[cl];
	for (map<char, int>::iterator i = opt.begin(); i != opt.end(); ++i)
	{
		cout << i->second << " " << i->first << endl;
		out[i->second-1] = (int) i->first - 'A';
	}
	print_rank(opt);
	cout << "Sum Kendall rank(optimal) : " << _kemeny_rule(opt, ranks, k) << endl;

	cout << "Heuristic :: " << endl;

	/*int **G = new int*[cl];
	for (int i = 0; i < cl; ++i)
	{
		G[i] = new int[cl];
		memset(G[i], 0, sizeof(int)*cl);
	}*/

	//pair<int, int> * outdegree = create_majority_graph(G, ranks_arr, count, k, cl);
	IGraph * graph = create_majority_graph(ranks_arr, k, m, cl);

	((Graph_Adj_Matrix*)graph)->printAsMatrix();
	cout << "Optimal with locally Kemenysation" << endl;
	print_array(out, cl);
	DVhamiltonPathForTournament(*graph, out);
	print_array(out, cl);
	opt = rank_from_array( out, cl);
	print_rank(opt);
	cout << "Sum Kendall rank(optimal with locally Kemenysation) : " << _kemeny_rule(opt, ranks, k) << endl;


	cout <<  "aqui : " ;
	int* igr = DVhamiltonPathForTournament(*graph);
	print_array(igr, cl);
	map<char, int> rank_heuI = rank_from_array( igr, cl);
	print_rank(rank_heuI);
	cout << "Sum Kendall rank(Heuristic) : " << _kemeny_rule(rank_heuI, ranks, k) << endl;

	cout << "FAS-pivot" << endl;
	int * FAS = graph::feedback_arc_set_pivot(*graph);

	map<char, int> rank_heuFAS = rank_from_array( FAS, cl);
	print_rank(rank_heuFAS);
	cout << "Sum Kendall rank(FAS Heuristic) : " << _kemeny_rule(rank_heuFAS, ranks, k) << endl;

	// Applying locally Kemenysation
	cout << "FAS-pivot`with locally Kemenysation" << endl;
	DVhamiltonPathForTournament(*graph, FAS);
	rank_heuFAS = rank_from_array( FAS, cl);
	print_rank(rank_heuFAS);
	cout << "Sum Kendall rank(FAS Heuristic) : " << _kemeny_rule(rank_heuFAS, ranks, k) << endl;

	cout << "topological" << endl;
	int * topsort = hamiltonPathForTournament(*graph);
	map<char, int> rank_heutopsort = rank_from_array( topsort, cl);
	print_rank(rank_heutopsort);
	cout << "Sum Kendall rank(topsort) : " << _kemeny_rule(rank_heutopsort, ranks, k) << endl;


	DVhamiltonPathForTournament(*graph, topsort);
	rank_heutopsort = rank_from_array( topsort, cl);
	print_rank(rank_heutopsort);
	cout << "Sum Kendall rank(topsort) : " << _kemeny_rule(rank_heutopsort, ranks, k) << endl;


	/*char * rh = hamiltonian_path_tournament(G, outdegree, cl);

	
	print_array(rh, cl);
	map<char, int> rank_heu = rank_from_array( rh, cl);

	print_rank(rank_heu);
	cout << "Sum Kendall rank(Heuristic) : " << _kemeny_rule(rank_heu, ranks, k) << endl;
	
	for (int i = 0; i < cl; ++i)
	{
		for (int j = 0; j < cl; ++j)
		{
			cout << setw(5) << G[i][j] << " ";	
		}

		cout << endl;
	}*/

	((Graph_Adj_Matrix*)graph)->printAsMatrix();

	return 0;
}