#ifndef __I_K_ALGORITHM__
#define __I_K_ALGORITHM__

#include "kstructure.h"
#include "kendall-tau-distance.hpp"

namespace klib{

	#include "Graph_Hamilton_Path.cpp"

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

	template<
		typename element_type,
		typename rank_type,
		typename permutation_impl_type
	>
	Permutation<element_type, rank_type> * gen_rank(int min, int max, int m, element_type first){
		Permutation<element_type, rank_type> * rank = new permutation_impl_type();

		int rank_pos = 1;
		int gen = 0;
		std::vector<element_type> v;

		for (int i = 0; i < max; ++i)
		{
			v.push_back(first + i);
		}

		while(v.size() > 0 && rank_pos <= m){
			std::vector<char>::iterator it = v.begin();
			gen = (int)(rand()%(v.size()));
			it += gen;
			rank->addElement(*it, rank_pos);

			rank_pos++;
			v.erase(it);
		}

		return rank;
	}

	template<
		typename element_type,
		typename rank_type
	>
	IGraph * create_majority_graph(const vector<Permutation<element_type, rank_type> *> &ranks, int k, int m, int n, element_type first){
		typedef Permutation<element_type, rank_type> perm_type;

		IGraph *g = new Graph_Adj_Matrix(n);
		bool *showed = new bool[n];

		//O(km²)
		for (int i = 0; i < ranks.size(); ++i)
		{
			memset(showed, 0, sizeof(bool)*n);
			perm_type *rank = ranks[i];
			for (typename perm_type::iterator j = rank->begin(); j != rank->end(); ++j)
			{
				typename perm_type::iterator l = j;
				for (l++; l != rank->end(); ++l)
				{	

					int u = (int)j->first - ((int)first);
					int v = (int)l->first - ((int)first);

					if(j->second < l->second){
						g->addEdge(u, v,  g->getEdge(u, v) + 1.0/k);
					}else{
						g->addEdge(v, u,  g->getEdge(v, u) + 1.0/k);
					}
					showed[u] = true;
					showed[v] = true;
				}
			}
			/*if(m < cl){
				for (int v = 0; v < cl; ++v)
				{
					if(!showed[v]){
						for (int u = 0; u < cl; ++u)
						{
							if(u != v)
								g->addEdge(u, v,  g->getEdge(u, v) + 1.0/rs);
						}
					}
				}
			}	*/	
		}
		cout << "graph created" << endl;
		((Graph_Adj_Matrix*)g)->printAsMatrix();

		IGraph * tour = new Graph_Adj_Matrix(n);
		cout << "instance" << endl; 
		complete2Tournament(*g, *tour);
		
		cout << "WEIGHTED : " << endl;
		((Graph_Adj_Matrix*)tour)->printAsMatrix();

		return tour;
	}

	template<
		typename element_type,
		typename rank_type,
		typename perm_type
	>
	double _kemeny_rule(perm_type &optimal, std::vector<perm_type *> &ranks){
		double sum = 0;

		#pragma omp parallel for
		for (int i = 0; i < ranks.size(); ++i)
		{	
			cout << "KKK" ; optimal.print(); cout << endl; ranks[i]->print(); cout << endl;
			sum += kendall_dist<element_type, rank_type>(optimal, *ranks[i]);
		}

		return sum;
	}

	template<
		typename element_type,
		typename rank_type
	>
	Permutation<element_type, rank_type> * kemeny_consensus(Permutation<element_type, rank_type> * sigma, vector<Permutation<element_type, rank_type> *> &ranks){
		
		Permutation<element_type, rank_type> *opt;
		vector<element_type> ini_rank;
		int n = sigma->size();
		
		ini_rank.resize(n);
		for (typename Permutation<element_type, rank_type>::iterator j = sigma->begin(); j != sigma->end(); ++j){
			ini_rank[j->second-1] = j->first;
		}

		std::sort (ini_rank.begin(), ini_rank.end());
		double min = 10000000;
		do {
			Permutation<element_type, rank_type> * r = new PermutationTree<element_type, rank_type>(ini_rank);

		 	double v = _kemeny_rule<element_type, rank_type>(*r, ranks);
		    if (min > v){
		    	opt = r;
		    	min = v;
		    }
		} while ( std::next_permutation(ini_rank.begin(), ini_rank.end()) );

		return opt;
	}

}

#endif