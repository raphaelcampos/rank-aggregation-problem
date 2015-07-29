
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

typedef unsigned short int element_type;
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
	//for (int i = 0; i < 30; ++i)
	for (int n = 4; n <= 7000; n = (n<10)?n+1:n*2)
	{
		int m = n;	

		for (int k = 10; k <= 10 ; k+=10)
		{
			vector<Permutation<element_type, rank_type>* > perms;
			
			double *sum_dist = new double[k];
			memset(sum_dist, 0, sizeof(double)*k);
			
			cout << "Generating syntetic rankings..." << endl;
			for (int i = 0; i < k; ++i)
			{	
			 	perms.push_back(gen_rank<element_type, rank_type, PermutationTree<element_type, rank_type> >(1, n, m, 0));
			}

			cout << "# rankings : " << k << endl;
			cout << "# candidates : " << n << endl;
			cout << "# top-m size : " << m << endl;

			Permutation<element_type, rank_type> * sigma = new PermutationTree<element_type, rank_type>();
			union_all(*sigma, perms);
			cout << sigma->size() << endl;
			cout << (sizeof(sigma) + sigma->size()*(sizeof(rank_type) + sizeof(element_type)))/pow(2.0,30) << "GB" << endl;
			cout << (sigma->size()*sigma->size()*(sizeof(double)+sizeof(IGraph::vertex*)) + sigma->size()*sizeof(IGraph::vertex))/pow(2.0,30) << "GB" << endl;
			
			//return EXIT_SUCCESS;
			//sigma->print(); cout << endl<< endl;

			timestamp_t t0 = get_timestamp();
			IGraph * tournament = create_majority_graph(perms, perms.size(), sigma->size(), sigma->size(), (unsigned short int)0);		
			cout << (get_timestamp() - t0) / 1000000.0L << endl;

			stringstream filename;
			filename << "rank-results/time/results_" << "k" << k << "_n" << n << "_m" << m;
			
			std::ofstream outfile;
			outfile.open(filename.str().c_str(), std::ios_base::app);
		
			outfile << n << " " << (get_timestamp() - t0) / 1000000.0L << " ";
			outfile.close();
			testing_all(*tournament, perms, *sigma, *sigma, filename.str());
			delete tournament;
			delete sigma;
		}

	}
	
	return EXIT_SUCCESS;
}