#ifndef __I_K_ALGORITHM__
#define __I_K_ALGORITHM__

#include "kstructure.h"
#include "kendall-tau-distance.hpp"
#include "kutils.hpp"

#include "Graph_Hamilton_Path.cpp"
#include "Graph_Feedback_Arc_Set.cpp"

using namespace graph;

namespace klib{

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
		}
		
		IGraph * tour = new Graph_Adj_Matrix(n);
		complete2Tournament(*g, *tour);
		
		
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
		    	delete opt;
		    	opt = r;
		    	min = v;
		    }else{
		    	delete r;
		    }
		} while ( std::next_permutation(ini_rank.begin(), ini_rank.end()) );

		return opt;
	}

	template<
		typename element_type,
		typename rank_type
	>
	Permutation<element_type, rank_type> * locally_kemenysation(IGraph &tour, Permutation<element_type, rank_type> &pi, element_type first){
		
		int * vertex_ids_perm = permutation_to_vertex_perm(pi, first);

		DVhamiltonPathForTournament(tour, vertex_ids_perm);

		return vertex_perm_to_permutation(vertex_ids_perm, tour.numVertices(), first);

	}

	template<
		typename element_type,
		typename rank_type
	>
	Permutation<element_type, rank_type> * kFAS_pivot(IGraph &tour, element_type first, bool kemenyzation = false){
		
		int * vertex_ids_perm = feedback_arc_set_pivot(tour);

		if(kemenyzation){
			DVhamiltonPathForTournament(tour, vertex_ids_perm);
		}
			
		return vertex_perm_to_permutation(vertex_ids_perm, tour.numVertices(), first);
	}

	template<
		typename element_type,
		typename rank_type
	>
	Permutation<element_type, rank_type> * FHP_greedy(IGraph &tour, element_type first, bool kemenyzation = false){
		
		int * vertex_ids_perm = hamiltonPathForTournament(tour);

		if(kemenyzation){
			DVhamiltonPathForTournament(tour, vertex_ids_perm);
		}
			
		return vertex_perm_to_permutation(vertex_ids_perm, tour.numVertices(), first);

	}

}

#endif