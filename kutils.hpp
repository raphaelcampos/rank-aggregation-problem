#ifndef __I_K_UTILS__
#define __I_K_UTILS__

#include "kstructure.h"

namespace klib{
	template<
		typename element_type
	>
	int * permutation_to_vertex_perm(Permutation<element_type, int> &p, element_type first){
		int * perm = new int[p.size()];

		for (typename Permutation<element_type, int>::iterator it = p.begin(); it != p.end(); ++it)
		{
			perm[it->second - 1] = it->first - first;
		}

		return perm;
	}

	template<
		typename element_type
	>
	Permutation<element_type, int>  * vertex_perm_to_permutation(int * vertex_ids, int n, element_type first){
		Permutation<element_type, int> * perm = new PermutationTree<element_type, int>;

		for (int i = 0 ;  i < n; i++)
		{
			perm->addElement(vertex_ids[i] + first, i + 1);
		}

		return perm;
	}

	template<
		typename element_type,
		typename rank_type,	
		class ForwardIt1,
		class ForwardIt2
	>
	void iter_swap(ForwardIt1 a, ForwardIt2 b)
	{
	   rank_type tmp; 
	   
	   tmp = a->second;
	   a->second = b->second;
	   b->second = tmp;
	}

	template<
		typename element_type,
		typename rank_type,
		class BidirIt
	>
	void reverse(BidirIt first, BidirIt last)
	{
	    while ((first != last) && (first != --last)) {
	        klib::iter_swap<element_type, rank_type>(first++, last);
	    }
	}

	template<
		typename element_type,
		typename rank_type,
		typename It
	>
	bool next_permutation(It begin, It end)
	{
	        if (begin == end)
	                return false;

	        It i = begin;
	        ++i;
	        if (i == end)
	                return false;

	        i = end;
	        --i;

	        while (true)
	        {
	                It j = i;
	                --i;

	                if ((*i).second < (*j).second)
	                {
	                        It k = end;

	                        while (!((*i).second < (*--k).second))
	                                /* pass */;

	                        klib::iter_swap<element_type, rank_type>(i, k);

	                        klib::reverse<element_type, rank_type>(j, end);
	                        return true;
	                }

	                if (i == begin)
	                {
	                        klib::reverse<element_type, rank_type>(begin, end);
	                        return false;
	                }
	        }
	}
}

#endif