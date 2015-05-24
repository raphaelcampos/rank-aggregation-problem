
#include <iomanip>
#include <iostream>

#include <string.h>
#include <time.h>
#include <cmath>

#include <algorithm>

#include <queue>
#include <vector>
#include <map>

#include "IGraph.h"
#include "Graph_Adj_Matrix2.hpp"
#include "Graph_Hamilton_Path.cpp"
#include "Graph_Search.cpp"
#include "Graph_Minimum_Spanning_Tree.cpp"
#include "Single_Source_Shortest_Path.cpp"

using namespace std;


void print_array(int a[], int n){
	for (int l = 0; l < n; ++l)
	{
		cout << a[l] << ",";
	}
	cout << endl;
	
}

void print_array(pair<int, int> a[], int n){
	for (int l = 0; l < n; ++l)
	{
		cout << a[l].first << " - " << a[l].second << ",";
	}
	cout << endl;
	
}

void print_array(char a[], int n){
	for (int l = 0; l < n; ++l)
	{
		cout << a[l] << ",";
	}
	cout << endl;
	
}

int merge_and_count(int *A, int *buffer, int ini, int fim, int meio){

	memset(buffer, 0, sizeof(int)*(fim-ini+1));

	int count = 0;

	int i = ini;
	int j = meio + 1;
	
	int k = 0;
	while(i < meio + 1 && j <= fim){
		if(A[i] <= A[j]){
			buffer[k] = A[i];
			i++;
		}else if(A[i] > A[j]){
			buffer[k] = A[j];
			count += meio - i + 1;
			j++;
		}
		k++;
	}

	if(i < meio + 1){
		memcpy(buffer + k, A + i, sizeof(int)*(meio-i+1));
	}else{
		memcpy(buffer + k, A + j, sizeof(int)*(fim-j+1));
	}
	
	memcpy(A + ini, buffer, sizeof(int)*(fim-ini+1));

	return count;
}

int sort_and_count(int A[], int buffer[], int ini, int fim){

	if(fim + 1 - ini <= 1){
		return 0;
	}else{
		int meio = (ini + fim)/2;
		return sort_and_count(A, buffer, ini, meio) +
		 sort_and_count(A, buffer, meio + 1, fim) + 
		 merge_and_count(A, buffer, ini, fim, meio);
		 ;
	}

}


map<char, int> rank_from_array(int a[], int n){
	std::map<char, int> rank;
	for (int i = 0; i < n; ++i)
	{
		rank[(char)'A' + a[i]] = i + 1;
	}
	return rank;
}

map<char, int> rank_from_array(char a[], int n){
	std::map<char, int> rank;
	for (int i = 0; i < n; ++i)
	{
		rank[a[i]] = i + 1;
	}
	return rank;
}

map<char, int> gen_rank(int min, int max, int k, char r[]){
	std::map<char, int> rank;

	int rank_pos = 1;
	int gen = 0;
	std::vector<char> v;

	for (int i = 0; i < max; ++i)
	{
		v.push_back((char)((int)'A' + i));
	}

	while(v.size() > 0 && rank_pos <= k){
		std::vector<char>::iterator it = v.begin();
		gen = (int)(rand()%(v.size()));
		it += gen;
		rank[*it] = rank_pos;
		r[rank_pos-1] = *it;
		rank_pos++;
		v.erase(it);
	}

	return rank;
}

void print_rank(std::map<char, int> &rank){
	std::map<char,int>::iterator it;
	for (it = rank.begin(); it!=rank.end(); ++it)
	{
		cout << it->first << " : " << it->second << endl;
	}
}

double kendall_dist(map<char,int> A, map<char,int> B){
	std::map<char,int>::iterator it;

	int total_size = B.size()+A.size();
	int *l = new int[total_size];
	int *buffer = new int[total_size];

	int size = 0;
	for (it = B.begin(); it!=B.end(); ++it)
	{
		if(A.count(it->first)){
			l[it->second-1] = A[it->first]; 
			A.erase(it->first);
		}else{
			l[it->second-1] = B.size() + it->second;
		}
		size++;
	}

	if(A.size()>0){
		for (it = A.begin(); it!= A.end(); ++it)
		{
			l[B.size() + it->second - 1] = A[it->first];		
			size++;
		}
	}

	for (int i = B.size(); i < total_size-1; ++i)
	{
		if(l[i] == 0 && l[i+1] != 0){
			l[i] = l[i+1];
			l[i+1] = 0;
		}
	}

	int s = sort_and_count(l, buffer, 0, size-1);

	delete [] buffer;
	delete [] l;
	
	return (2*s/(double)(size*(size-1))); 
}

double _kemeny_rule(map<char, int> &optimal, map<char, int> ranks[], int n){
	double sum = 0;

	#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		sum += kendall_dist(optimal, ranks[i]);
	}

	return sum;
}

map<char, int> kemeny_consensus(char ini_rank[], int n, map<char, int> ranks[], int rs){

	map<char, int> opt;

	std::sort (ini_rank, ini_rank + n);

	double min = 10000000;

	
	do {
		map<char, int> r = rank_from_array(ini_rank, n);

	 	double v = _kemeny_rule(r, ranks, rs);
	    if (min > v){
	    	opt = r;
	    	min = v;
	    }
	} while ( std::next_permutation(ini_rank, ini_rank + n) );

	return opt;
}

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

/*void DFS(int i, int **G, int *visited, int n)
{
    int j;
    printf("\n%C", 'A' + i);
    visited[i]=1;
    for(j=0;j<n;j++)
        if(!visited[j]&&G[i][j]==1)
            DFS(j, G, visited, n);
}*/

class GreaterOutDegree
{
  	bool reverse;
	public:
	  	GreaterOutDegree(const bool& revparam=false)
		    {reverse=revparam;}
		  bool operator() (const pair<int,int>& lhs, const pair<int,int>&rhs) const
		  {
		    if (reverse) return (lhs.first > rhs.first);
		    else return (lhs.first<rhs.first);
		  }
};

/**
 * [create_majority_graph description]
 * @param  ranks [description]
 * @param  rs    [description]
 * @param  k     [description]
 * @param  cl    [description]
 * @return       [description]
 */
IGraph * create_majority_graph(char * ranks[], int rs, int k, int cl){

	IGraph *g = new Graph_Adj_Matrix(cl);

	for (int i = 0; i < rs; ++i)
	{
		for (int j = 0; j < k; ++j)
		{
			for (int l = j + 1; l < k; ++l)
			{	
				int u = (int)ranks[i][j] - ((int)'A');
				int v = (int)ranks[i][l] - ((int)'A');

				g->addEdge(u, v,  g->getEdge(u, v) + 1);
				
			}
		}
	}

	IGraph *gmst = new Graph_Adj_Matrix(8);
	IGraph::vertex * s = &(*gmst->begin());
	
	// A
	gmst->addEdge(0,1,1);
	gmst->addEdge(0,2,7);
	// B
	gmst->addEdge(1,0,1);
	gmst->addEdge(1,2,4);
	gmst->addEdge(1,3,9);
	gmst->addEdge(1,4,6);
	// C
	gmst->addEdge(2,0,7);
	gmst->addEdge(2,1,4);
	gmst->addEdge(2,7,8);
	// D
	gmst->addEdge(3,2,1);
	gmst->addEdge(3,4,11);
	gmst->addEdge(3,6,2);
	gmst->addEdge(3,7,5);
	// E
	gmst->addEdge(4,1,6);
	gmst->addEdge(4,3,11);
	gmst->addEdge(4,5,3);
	// F
	gmst->addEdge(5,4,3);
	gmst->addEdge(5,6,10);
	// G
	gmst->addEdge(6,3,2);
	gmst->addEdge(6,5,10);
	gmst->addEdge(6,7,13);
	// H
	gmst->addEdge(7,2,8);
	gmst->addEdge(7,3,5);
	gmst->addEdge(7,6,13);
	
	MST_prim(*gmst, *s);
	graph::dijkstra(*gmst, *s);
	graph::bellman_ford(*gmst, *s);

	IGraph * tour = new Graph_Adj_Matrix(g->numVertices());

	complete2Tournament(*g, *tour);

	s = &(*tour->begin());
	//BFS(*tour, *s);

	DFS(*tour, *s);
	s = &(*(tour->begin()++));
	DFS(*tour, *s);

	//cout << (int)((Graph_Adj_Matrix*)tour)->thereIsUniversalSink() << endl;
	((Graph_Adj_Matrix*)tour)->printAsMatrix();

	return tour;
}

pair<int, int> * create_majority_graph(int **G, char **ranks, int rs, int k, int n){

	pair<int, int> *outdegree = new pair<int, int>[n];
	int *distance = new int[n];
	int *visited = new int[n];
	queue<int> q;

	memset(outdegree, 0, sizeof(pair<int, int>)*n);
	memset(distance, 0, sizeof(int)*n);
	memset(visited, 0, sizeof(int)*n);

	for (int i = 0; i < rs; ++i)
	{
		for (int j = 0; j < k; ++j)
		{
			for (int l = j + 1; l < k; ++l)
			{	
				int u = (int)ranks[i][j] - ((int)'A');		
				int v = (int)ranks[i][l] - ((int)'A');
				G[u][v] += 1;
			}
		}
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			cout << setw(5) << G[i][j] << " ";	
		}

		cout << endl;
	}

	for (int u = 0; u < n; ++u)
	{
		for (int v = u + 1; v < n; ++v)
		{
			//if(G[u][v] == 0 && G[v][u] == 0) continue;

			if(G[u][v] > G[v][u]){
				G[u][v] = 1;
				G[v][u] = 0;

				outdegree[u].first += 1;
				outdegree[u].second = u;
			}else{
				G[u][v] = 0;
				G[v][u] = 1;

				outdegree[v].first += 1;
				outdegree[v].second = v;			
			}		
		}
	}

	return outdegree;
}

char * hamiltonian_path_tournament(int **G, pair<int, int> *outdegree, int n){

	char * path = new char[n];

	queue< pair<int, int> > Q;
	int max = 0;
	
	for (int i = 1; i < n; ++i)
	{
		if(outdegree[max].first < outdegree[i].first){
			max = i;
		}
	}
	Q.push(outdegree[max]);

	int it = 0;
	int v;
	while(!Q.empty() && it < n){
		v = Q.front().second;
		Q.pop();

		outdegree[v].first = 0;
		max = 0;
		int left_adj = -1;
		for (int u = 0; u < n; ++u)
		{
			if(G[u][v]){
				G[u][v] = 0;
				outdegree[u].first--;
			}
			if(G[v][u] && outdegree[u].first > outdegree[max].first){
				max = u;
			}
		}

		if(outdegree[max].first){
			Q.push(outdegree[max]);
		}

		path[it] = (char)('A' + v);
		it++;
	}

	for (int u = 0; u < n; ++u)
	{
		if(G[v][u]){
			path[it] = (char)('A' + u);
			break;
		}
	}

	return path;
}

int main(int argc, char const *argv[])
{

	srand(time(NULL));
	
	int count = 3;
	int k = 8; // top-k rank
	map<char, int> *ranks = new map<char, int>[count];
	char **ranks_arr = new char*[count];
	map<char, int> opt;
	
	double *sum_dist = new double[count];
	memset(sum_dist, 0, sizeof(double)*count);
	
	
	for (int i = 0; i < count; ++i)
	{
		
		ranks_arr[i] = new char[k];
	 	ranks[i] = gen_rank(1, k, k, ranks_arr[i]);
	 	//cout << ranks_arr[i] << endl;
	}

    cout << "Gen all..." <<endl;

	for (int i = 0; i < count; ++i)
	{
		for (int j = i+1; j < count; ++j)
		{
			double dist = kendall_dist(ranks[i], ranks[j]);
			sum_dist[i] += dist;
			sum_dist[j] += dist;

			//cout << i << " - " << j << " : " << dist << endl;
		}
	}

	int min_dist = 0;
	for (int i = 1; i < count; ++i)
	{
		//print_rank(ranks[i]);
		//cout << "Sum Kendall rank(" << i + 1 << ") : " << sum_dist[i] << endl;
		if(sum_dist[i] < sum_dist[min_dist]){
			min_dist = i;
		}
	}

	print_rank(ranks[min_dist]);
	cout << "Sum Kendall rank(" << min_dist + 1 << ") : " << sum_dist[min_dist] << endl;

	int cl;
	char *c = union_all(ranks, count, cl);
	cout << "size : " << cl << endl;
	print_array(c, cl);
	cout << "OPTIMAL :: " << endl;
	opt = kemeny_consensus(c, cl, ranks, count);
	print_rank(opt);
	cout << "Sum Kendall rank(optimal) : " << _kemeny_rule(opt, ranks, count) << endl;

	cout << "Heuristic :: " << endl;

	int **G = new int*[cl];
	for (int i = 0; i < cl; ++i)
	{
		G[i] = new int[cl];
		memset(G[i], 0, sizeof(int)*cl);
	}

	pair<int, int> * outdegree = create_majority_graph(G, ranks_arr, count, k, cl);
	IGraph * graph = create_majority_graph(ranks_arr, count, k, cl);
	int* igr = DVhamiltonPathForTournament(*graph);

	cout <<  "aqui : " ;
	print_array(igr, cl);
	map<char, int> rank_heuI = rank_from_array( igr, cl);
	print_rank(rank_heuI);
	cout << "Sum Kendall rank(Heuristic) : " << _kemeny_rule(rank_heuI, ranks, count) << endl;

	
	char * rh = hamiltonian_path_tournament(G, outdegree, cl);
	
	print_array(rh, cl);
	map<char, int> rank_heu = rank_from_array( rh, cl);

	print_rank(rank_heu);
	cout << "Sum Kendall rank(Heuristic) : " << _kemeny_rule(rank_heu, ranks, count) << endl;

	for (int i = 0; i < cl; ++i)
	{
		for (int j = 0; j < cl; ++j)
		{
			cout << setw(5) << G[i][j] << " ";	
		}

		cout << endl;
	}

	((Graph_Adj_Matrix*)graph)->printAsMatrix();

	return 0;
}