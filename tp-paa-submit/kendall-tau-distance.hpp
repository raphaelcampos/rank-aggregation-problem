#ifndef __I_K_KENDALL_DIST__
#define __I_K_KENDALL_DIST__

#include "kstructure.h"

namespace klib{

	template<
		typename rank_type
	>
	int merge_and_count(rank_type *A, rank_type *buffer, int ini, int fim, int meio){

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

	template<
		typename rank_type
	>
	int sort_and_count(rank_type A[], rank_type buffer[], int ini, int fim){

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

	template<
		typename element_type,
		typename rank_type,
		typename perm_type
	>
	double kendall_dist(perm_type &A, perm_type &B){

		int n = B.size();
		int *l = new rank_type[n];
		int *buffer = new rank_type[n];

		for (typename perm_type::iterator it = B.begin(); it != B.end(); ++it)
		{
			l[it->second - 1] = A(it->first);
		}

		int s = klib::sort_and_count(l, buffer, 0, n-1);

		delete [] buffer;
		delete [] l;
		
		return (2*s/(double)(n*(n-1))); 
	}

}

#endif