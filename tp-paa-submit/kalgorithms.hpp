#ifndef __I_K_ALGORITHM__
#define __I_K_ALGORITHM__

#include <fstream>

#include "kstructure.h"
#include "kendall-tau-distance.hpp"
#include "kutils.hpp"

#include "Graph_Hamilton_Path.cpp"
#include "Graph_Feedback_Arc_Set.cpp"

#include <sys/time.h>
#include <unistd.h>

using namespace graph;
using namespace std;

namespace klib{

	typedef unsigned long long timestamp_t;

	static timestamp_t get_timestamp ()
	{
	  struct timeval now;
	  gettimeofday (&now, NULL);
	  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
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
			typename std::vector<element_type>::iterator it = v.begin();
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
			if(m < n){
				for (int v = 0; v < n; ++v)
				{
					if(!showed[v]){
						for (int u = 0; u < n; ++u)
						{
							if(u != v)
								g->addEdge(u, v,  g->getEdge(u, v) + 1.0/k);
						}
					}
				}
			}		
		}
		
		IGraph * tour = new Graph_Adj_Matrix(n);
		complete2Tournament(*g, *tour);
		
		return tour;
	}

	template<
		typename element_type,
		typename rank_type
	>
	IGraph * create_majority_graph(const vector<Permutation<element_type, rank_type> *> &ranks, int k, int m, int n, Permutation<element_type, int> &table){
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
					int u = table(j->first) - 1;
					int v = table(l->first) - 1;
					if(j->second < l->second){
						g->addEdge(u, v,  g->getEdge(u, v) + 1.0/k);
					}else{
						g->addEdge(v, u,  g->getEdge(v, u) + 1.0/k);
					}
					showed[u] = true;
					showed[v] = true;
				}
			}
			if(m < n){
				for (int v = 0; v < n; ++v)
				{
					if(!showed[v]){
						for (int u = 0; u < n; ++u)
						{
							if(u != v)
								g->addEdge(u, v,  g->getEdge(u, v) + 1.0/k);
						}
					}
				}
			}		
		}
		
		IGraph * tour = new Graph_Adj_Matrix(n);
		complete2Tournament(*g, *tour);
		
		delete g;
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

		return sum/ranks.size();
	}

	template<
		typename element_type,
		typename rank_type
	>
	Permutation<element_type, rank_type> * kemeny_consensus(Permutation<element_type, rank_type> * sigma, vector<Permutation<element_type, rank_type> *> &ranks){
		
		Permutation<element_type, rank_type> *opt = NULL;
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
		    	if(opt != NULL) delete opt;
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
	Permutation<element_type, rank_type> * locally_kemenysation(IGraph &tour, Permutation<element_type, rank_type> &pi, Permutation<element_type, int>& table){
		
		int * vertex_ids_perm = permutation_to_vertex_perm(pi, table);

		DVhamiltonPathForTournament(tour, vertex_ids_perm);

		return vertex_perm_to_permutation(vertex_ids_perm, tour.numVertices(), table);

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
	Permutation<element_type, rank_type> * kFAS_pivot(IGraph &tour, Permutation<element_type, int>& table, bool kemenyzation = false){
		
		int * vertex_ids_perm = feedback_arc_set_pivot(tour);

		if(kemenyzation){
			DVhamiltonPathForTournament(tour, vertex_ids_perm);
		}
			
		return vertex_perm_to_permutation(vertex_ids_perm, tour.numVertices(), table);
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

	template<
		typename element_type,
		typename rank_type
	>
	Permutation<element_type, rank_type> * FHP_greedy(IGraph &tour, Permutation<element_type, int>& table, bool kemenyzation = false){
		
		int * vertex_ids_perm = hamiltonPathForTournament(tour);
		if(kemenyzation){
			DVhamiltonPathForTournament(tour, vertex_ids_perm);
		}
			
		return vertex_perm_to_permutation(vertex_ids_perm, tour.numVertices(), table);

	}

	template<
		typename element_type,
		typename rank_type
	>
	Permutation<element_type, rank_type> * pick_a_list(vector<Permutation<element_type, rank_type> *> &perms){
		
		const int k = perms.size();

		double *sum_dist = new double[k];
		memset(sum_dist, 0, sizeof(double)*k);

		Permutation<element_type, rank_type> * sigma = new PermutationTree<element_type, rank_type>();
		union_all(*sigma, perms);
		 
		for (int i = 0; i < k; ++i)
		{
			Permutation<element_type, rank_type> * perm = new PermutationTree<element_type, rank_type>();
			
			for (typename Permutation<element_type, rank_type>::iterator it = perms[i]->begin(); it != perms[i]->end(); ++it)
			{
				perm->addElement(it->first, it->second);
			}	
			union_perm(*perm, *sigma);
			

			sum_dist[i] = _kemeny_rule<element_type, rank_type>(*perm, perms);

			delete perm;

			
		}

		int min_dist = 0;
		for (int i = 1; i < k; ++i)
		{
			if(sum_dist[i] < sum_dist[min_dist]){
				min_dist = i;
			}
		}

		Permutation<element_type, rank_type> * result = new PermutationTree<element_type, rank_type>();;
		for (typename Permutation<element_type, rank_type>::iterator it = perms[min_dist]->begin(); it != perms[min_dist]->end(); ++it)
		{
			result->addElement(it->first, it->second);
		}	
		union_perm(*result, *sigma);
		return result;
	}

	template<
		typename element_type,
		typename rank_type
	>
	pair<double, Permutation<element_type, rank_type>* > mixed(IGraph &tournament, vector<Permutation<element_type, rank_type> *> &perms, element_type first){
		
		typedef pair<double, Permutation<element_type, rank_type> *> result;

		const int k = perms.size();
		std::vector<result> results;

		Permutation<element_type, rank_type> * pal = pick_a_list<element_type, rank_type>(perms);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*pal, perms), pal));

		Permutation<element_type, rank_type> * kemenized = locally_kemenysation(tournament, *pal, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*kemenized, perms), kemenized));
		
		Permutation<element_type, rank_type> * fas_pivot = kFAS_pivot<element_type, rank_type>(tournament, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fas_pivot, perms), fas_pivot));

		Permutation<element_type, rank_type> * fas_pivot_kemenized = kFAS_pivot<element_type, rank_type>( tournament, first, true);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fas_pivot_kemenized, perms), fas_pivot_kemenized));

		Permutation<element_type, rank_type> * fhp_greedy = FHP_greedy<element_type, rank_type>(tournament, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fhp_greedy, perms), fhp_greedy));

		Permutation<element_type, rank_type> * fhp_greedy_kemenized = FHP_greedy<element_type, rank_type>( tournament, first, true);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fhp_greedy_kemenized, perms), fhp_greedy_kemenized));

		return *min_element(results.begin(), results.end());
	}

	template<
		typename element_type,
		typename rank_type
	>
	pair<double, Permutation<element_type, rank_type>* > mixed(IGraph &tournament, vector<Permutation<element_type, rank_type> *> &perms, Permutation<element_type, int> &first){
		
		typedef pair<double, Permutation<element_type, rank_type> *> result;

		const int k = perms.size();
		std::vector<result> results;

		Permutation<element_type, rank_type> * pal = pick_a_list<element_type, rank_type>(perms);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*pal, perms), pal));

		Permutation<element_type, rank_type> * kemenized = locally_kemenysation(tournament, *pal, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*kemenized, perms), kemenized));
		
		Permutation<element_type, rank_type> * fas_pivot = kFAS_pivot<element_type, rank_type>(tournament, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fas_pivot, perms), fas_pivot));

		Permutation<element_type, rank_type> * fas_pivot_kemenized = kFAS_pivot<element_type, rank_type>( tournament, first, true);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fas_pivot_kemenized, perms), fas_pivot_kemenized));

		Permutation<element_type, rank_type> * fhp_greedy = FHP_greedy<element_type, rank_type>(tournament, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fhp_greedy, perms), fhp_greedy));

		Permutation<element_type, rank_type> * fhp_greedy_kemenized = FHP_greedy<element_type, rank_type>( tournament, first, true);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fhp_greedy_kemenized, perms), fhp_greedy_kemenized));

		return *min_element(results.begin(), results.end());
	}

	template<
		typename element_type,
		typename rank_type
	>
	void testing_all(IGraph &tournament, vector<Permutation<element_type, rank_type> *> &perms, Permutation<element_type, rank_type> &sigma, element_type first, string filename){
		
		typedef pair<double, Permutation<element_type, rank_type> *> result;

		const int k = perms.size();
		std::vector<result> results;

		if(sigma.size() < 10){
			Permutation<element_type, rank_type> * opt = kemeny_consensus(&sigma, perms);
			opt->print(); cout << endl;
			results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*opt, perms), opt));
		}
		Permutation<element_type, rank_type> * pal = pick_a_list<element_type, rank_type>(perms);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*pal, perms), pal));
		pal->print(); cout << endl;

		Permutation<element_type, rank_type> * kemenized = locally_kemenysation(tournament, *pal, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*kemenized, perms), kemenized));
		kemenized->print(); cout << endl;


		Permutation<element_type, rank_type> * fas_pivot = kFAS_pivot<element_type, rank_type>(tournament, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fas_pivot, perms), fas_pivot));
		fas_pivot->print(); cout << endl;


		Permutation<element_type, rank_type> * fas_pivot_kemenized = kFAS_pivot<element_type, rank_type>( tournament, first, true);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fas_pivot_kemenized, perms), fas_pivot_kemenized));

		Permutation<element_type, rank_type> * fhp_greedy = FHP_greedy<element_type, rank_type>(tournament, first);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fhp_greedy, perms), fhp_greedy));

		Permutation<element_type, rank_type> * fhp_greedy_kemenized = FHP_greedy<element_type, rank_type>( tournament, first, true);
		results.push_back(make_pair(_kemeny_rule<element_type, rank_type>(*fhp_greedy_kemenized, perms), fhp_greedy_kemenized));

		results.push_back(mixed(tournament, perms, first));

		std::ofstream outfile;
		outfile.open(filename.c_str(), std::ios_base::app);
		
		for (int i = 0; i < results.size(); ++i)
		{
			outfile << setprecision(10) << results[i].first << " ";
		}
		outfile << endl;
		outfile.close();
	}

	template<
		typename element_type,
		typename rank_type
	>
	void testing_all(IGraph &tournament, vector<Permutation<element_type, rank_type> *> &perms, Permutation<element_type, rank_type> &sigma,  Permutation<element_type, int> &first, string filename){
		
		typedef double result;

		const int k = perms.size();
		const int n = sigma.size();
		
		std::vector<result> results;

		timestamp_t t0;

		if(sigma.size() < 13){
			t0 = get_timestamp();
			Permutation<element_type, rank_type> * opt = kemeny_consensus(&sigma, perms);
			results.push_back((get_timestamp() - t0) / 1000000.0L);
		}else{
			results.push_back(100000);
		}
		t0 = get_timestamp();
		Permutation<element_type, rank_type> * pal = pick_a_list<element_type, rank_type>(perms);
		results.push_back((get_timestamp() - t0) / 1000000.0L);

		Permutation<element_type, rank_type> * kemenized = locally_kemenysation(tournament, *pal, first);
		results.push_back((get_timestamp() - t0) / 1000000.0L);
		
		t0 = get_timestamp();
		Permutation<element_type, rank_type> * fas_pivot = kFAS_pivot<element_type, rank_type>(tournament, first);
		results.push_back((get_timestamp() - t0) / 1000000.0L);

		Permutation<element_type, rank_type> * fas_pivot_kemenized = kFAS_pivot<element_type, rank_type>( tournament, first, true);
		results.push_back((get_timestamp() - t0) / 1000000.0L);

		t0 = get_timestamp();
		Permutation<element_type, rank_type> * fhp_greedy = FHP_greedy<element_type, rank_type>(tournament, first);
		results.push_back((get_timestamp() - t0) / 1000000.0L);

		Permutation<element_type, rank_type> * fhp_greedy_kemenized = FHP_greedy<element_type, rank_type>( tournament, first, true);
		results.push_back((get_timestamp() - t0) / 1000000.0L);

		t0 = get_timestamp();
		mixed(tournament, perms, first);
		results.push_back((get_timestamp() - t0) / 1000000.0L);

		std::ofstream outfile;
		outfile.open(filename.c_str(), std::ios_base::app);
		
		for (int i = 0; i < results.size(); ++i)
		{
			outfile << setprecision(10) << results[i] << " ";
		}
		outfile << endl;
		outfile.close();
	}

}

#endif