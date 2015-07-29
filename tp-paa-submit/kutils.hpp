#ifndef __I_K_UTILS__
#define __I_K_UTILS__

#include <fstream>

#include "kstructure.h"
#include "utils.h"

namespace klib{
	template<
		typename arr_type
	>
	void print_array(arr_type a[], int n){
		for (int l = 0; l < n; ++l)
		{
			cout << a[l] << ",";
		}
		cout << endl;
		
	}

	template<
		typename element_type,
		typename rank_type
	>
	void union_perm(Permutation<element_type, rank_type> &sigma, Permutation<element_type, rank_type> &tau){
		
		for (typename Permutation<element_type, rank_type>::iterator it = tau.begin(); it != tau.end(); ++it)
		{
			sigma.addElement(it->first, sigma(it->first));
		}
	}

	template<
		typename element_type,
		typename rank_type
	>
	void union_all(Permutation<element_type, rank_type> &sigma, vector<Permutation<element_type, rank_type> *>&perms){
		//Permutation<element_type, rank_type> * perm = perms[0];
		for (int i = 0; i < perms.size(); ++i)
		{	
			for (typename Permutation<element_type, rank_type>::iterator it = perms[i]->begin(); it != perms[i]->end(); ++it)
			{
				sigma.addElement(it->first, sigma(it->first));
			}
		}
	}


	template<
		typename element_type,
		typename rank_type,
		typename perm_imp_type
	>
	pair<std::vector<Permutation<element_type, rank_type>* >, std::map<string, int> >  load_format_trec(string qresult){
		std::string line;
		std::map<string, int> map_queries;
		std::vector<Permutation<element_type, rank_type>* > rankings;
		
		ifstream infile;
    	infile.open(qresult.c_str());

    	while(infile.good() && std::getline( infile, line )){
    		std::vector<std::string> strs = split(line, ' ');

    		element_type elem = atoi(&(strs[2]).c_str()[0]);
    		rank_type pos = atoi(&(strs[3]).c_str()[0]);

    		std::map<string, int>::iterator it = map_queries.find(strs[0]);
            if(it != map_queries.end()){
                rankings[it->second]->addElement(elem, pos);
            }else{
            	Permutation<element_type, rank_type>* p = new perm_imp_type();
    			p->addElement(elem, pos);
    			rankings.push_back(p);
    			map_queries[strs[0]] = map_queries.size() - 1;
            }
    	}

    	infile.close();
    	return make_pair(rankings, map_queries);
	}

	template<
		typename element_type,
		typename rank_type
	>
	void save_format_trec(string qresult, pair< std::vector<Permutation<element_type, rank_type>* >, std::map<string, int> > &results){
		
		ofstream outfile;
    	outfile.open(qresult.c_str());

    	for (std::map<string, int>::iterator m = results.second.begin(); m != results.second.end(); ++m)
    	{
    		int i = m->second;
    		for (typename Permutation<element_type, rank_type>::iterator it = results.first[i]->begin(); it != results.first[i]->end(); ++it)
    		{
    			outfile << m->first << " 0 " << it->first << " " << it->second << " " << results.first[i]->size()-it->second << " STANDARD" << endl;
    		}
     	}

    	outfile.close();
	}

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
	int * permutation_to_vertex_perm(Permutation<element_type, int> &p, Permutation<element_type, int> &table){
		int * perm = new int[p.size()];

		for (typename Permutation<element_type, int>::iterator it = p.begin(); it != p.end(); ++it)
		{
			perm[it->second - 1] = table(it->first) - 1;
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
		typename element_type
	>
	Permutation<element_type, int>  * vertex_perm_to_permutation(int * vertex_ids, int n, Permutation<element_type, int> &table){
		Permutation<element_type, int> * perm = new PermutationTree<element_type, int>;

		for (int i = 0 ;  i < n; i++)
		{
			element_type elem =  table[vertex_ids[i]+1];
			perm->addElement(elem, i + 1);
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