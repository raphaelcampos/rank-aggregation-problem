
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
	for (int n = 5; n < 11; ++n)
	{
		int m = n-2;	

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

			IGraph * tournament = create_majority_graph(perms, k, m, sigma->size(), 'A');
			((Graph_Adj_Matrix*)tournament)->printAsMatrix();
			cout << endl;

			stringstream filename;
			
			filename << "rank-results/partial list/results_" << "k" << k << "_n" << n << "_m" << m;
			
			testing_all(*tournament, perms, *sigma, 'A', filename.str());
		}

	}

	return EXIT_SUCCESS;
}