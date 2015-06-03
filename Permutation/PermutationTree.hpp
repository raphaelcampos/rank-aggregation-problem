#ifndef __I_PERMUT_TREE__
#define __I_PERMUT_TREE__

#include <iostream>
#include <map>

#include "Permutation.h"

using namespace std;

template<
    typename element_type,
    typename rank_type
>
class PermutationTree : public Permutation<element_type, rank_type>{
    typedef std::map<element_type, rank_type> storage_type;
    public:
        class iterator_imp : public Permutation<element_type, rank_type>::iterator_imp{
            public:
                void next(){
                    ptr_++;
                }

                void prev(){
                    ptr_--;
                }
                
                typename Permutation<element_type, rank_type>::iterator_imp * clone(){
                    return new iterator_imp(ptr_);
                }
                
                pair<const element_type, rank_type> * get() const{  
                    return &(*ptr_);
                }

                bool isEqual(const typename Permutation<element_type, rank_type>::iterator_imp& other) const{
                    return get() == other.get();
                }

                iterator_imp(typename map<element_type, rank_type >::iterator ptr){
                    this->ptr_ = ptr;
                }

            private:
                typename map< element_type, rank_type >::iterator ptr_;
        };

        PermutationTree(){
            highest = 0;
        }

        PermutationTree(const std::vector<element_type> &t){
            for (int i = 0; i < t.size(); ++i)
            {
                perm[t[i]] = i + 1; 
            }
        }

        /**
         * [addElement description]
         * @param e [description]
         * @param r [description]
         */
        void addElement(element_type e, rank_type r){
            perm[e] = r;
            if(highest < r){
                highest = r;
            }
        }
        /**
         * given the vertice id v, remove it from
         * the graph
         * @param v     vertice id
         */
        void removeElement(element_type e){

        }

        int size(){
            return perm.size();
        }

        /**
         * O(log n)
         * @param  e element in the list
         * @return element's rank in the list
         */
        rank_type operator()(element_type e){
            it = perm.find(e);
            if(it != perm.end()){
                return it->second;
            }else{
                return highest + 1;
            }
        }

        void print(){
            cout << "{" ;
            int i = 0;
            for (it = perm.begin(); it != perm.end() && i < perm.size() - 1; ++it, i++)
            {
                cout << it->first << ":" << it->second << ",";
            }

            cout << it->first << ":" << it->second << "}";
        }

        typename Permutation<element_type, rank_type>::iterator begin(){
            iterator_imp * first = new iterator_imp(this->perm.begin());
            
            typename Permutation<element_type, rank_type>::iterator f(first);

            return f;
        }

        typename Permutation<element_type, rank_type>::iterator end(){
            iterator_imp * last = new iterator_imp(this->perm.end());

            typename Permutation<element_type, rank_type>::iterator l(last);

            return l;
        }

    private:
        storage_type perm;
        typename storage_type::iterator it;
        rank_type highest;



};

#endif 