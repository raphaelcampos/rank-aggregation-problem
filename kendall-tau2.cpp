
#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <time.h>
#include <algorithm>

using namespace std;

void print_array(int a[], int n){
	for (int l = 0; l < n; ++l)
	{
		cout << a[l] << ",";
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


unordered_map<char, int> rank_from_array(char a[], int n){
	std::unordered_map<char, int> rank;
	for (int i = 0; i < n; ++i)
	{
		rank[a[i]] = i + 1;
	}
	return rank;
}

unordered_map<char, int> gen_rank(int min, int max, int k){
	std::unordered_map<char, int> rank;

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
		rank_pos++;
		v.erase(it);
	}

	return rank;
}

void print_rank(std::unordered_map<char, int> &rank){
	std::unordered_map<char,int>::iterator it;
	for (it = rank.begin(); it!=rank.end(); ++it)
	{
		cout << it->first << " : " << it->second << endl;
	}
}

double kendall_dist(unordered_map<char,int> A, unordered_map<char,int> B){
	std::unordered_map<char,int>::iterator it;

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

double _kemeny_rule(unordered_map<char, int> &optimal, unordered_map<char, int> ranks[], int n){
	double sum = 0;

	#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		sum += kendall_dist(optimal, ranks[i]);
	}

	return sum;
}

unordered_map<char, int> kemeny_consensus(char ini_rank[], int n, unordered_map<char, int> ranks[], int rs){

	unordered_map<char, int> opt;

	std::sort (ini_rank, ini_rank + n);

	double min = 10000000;

	
	do {
		unordered_map<char, int> r = rank_from_array(ini_rank, n);

	 	double v = _kemeny_rule(r, ranks, rs);
	    if (min > v){
	    	opt = r;
	    	min = v;
	    }
	} while ( std::next_permutation(ini_rank, ini_rank + n) );

	return opt;
}

char* union_all(unordered_map<char, int> ranks[], int n, int &size){
	unordered_map<char, int> rank = ranks[0];
	for (int i = 1; i < n; ++i)
	{	
		
		for (unordered_map<char, int>::iterator it = ranks[i].begin(); it !=  ranks[i].end(); ++it)
		{
			rank[it->first] = 0; 
		}
	}

	int i = 0;
	char *c = new char[rank.size()];
	for (unordered_map<char, int>::iterator it = rank.begin(); it !=  rank.end(); ++it)
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
	
	int count = 50;
	unordered_map<char, int> *ranks = new unordered_map<char, int>[count];
	unordered_map<char, int> opt;
	
	double *sum_dist = new double[count];
	memset(sum_dist, 0, sizeof(double)*count);
	
	
	for (int i = 0; i < count; ++i)
	{
	 	ranks[i] = gen_rank(1,10,10);
	}

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

	return 0;
}