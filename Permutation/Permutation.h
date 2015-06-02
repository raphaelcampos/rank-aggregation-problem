#ifndef __I_PERMUT__
#define __I_PERMUT__

using namespace std;

template<
	typename element_type,
	typename rank_type
>
class Permutation{

	public:
        class iterator_imp{
            public:

                virtual ~iterator_imp(){}
                
                virtual void next(){}

                virtual void prev(){}

                virtual Permutation::iterator_imp* clone(){}
                
                virtual pair<const element_type, rank_type> * get() const{}

                virtual bool isEqual(const Permutation::iterator_imp& other) const{}
        };

    public:
        class iterator{
            public:
                typedef iterator self_type;
                typedef pair<const element_type, rank_type> value_type;
                typedef value_type& reference;
                typedef value_type* pointer;
                typedef int difference_type;
                self_type operator--(){
                    imp->prev();
                    return *this;
                }
                self_type operator++(){
                    imp->next();
                    return *this;
                }
                self_type operator++(int junk){
                    imp->next();
                    return *this;
                }
                reference operator*(){
                    return *(imp->get());
                }
                pointer operator->(){
                    return imp->get();
                }
                bool operator==(const self_type& rhs){
                    return imp->isEqual(*rhs.imp);
                }
                bool operator!=(const self_type& rhs){
                    return !imp->isEqual(*rhs.imp);
                }

                iterator( const iterator& other ){
                	this->	imp = other.imp->clone();
			  	}

                iterator(iterator_imp * imp){
                    this->imp = imp;
                }

                ~iterator(){
                    delete imp;
                }

            private:
                iterator_imp * imp;
        };

	public:


		virtual ~Permutation() {}
        
        virtual void addElement(element_type e, rank_type r) = 0;
        /**
         * given the vertice id v, remove it from
         * the graph
         * @param v     vertice id
         */
        virtual void removeElement(element_type e) = 0;

        virtual int size() = 0;
        
        /**
         * [getVertex description]
         * @param  v [description]
         * @return   [description]
         */
        virtual rank_type operator()(element_type e) = 0;

        virtual void print() = 0;

        /**
         * iterator at first position
         * it allows iterations on 
         * the vertice's storage
         * @return IGraph::vertex_iterator
         */
        virtual iterator begin() = 0;

        /**
         * iterator for the end
         * it determine the end of iteration
         * @return IGraph::vertex_iterator
         */
        virtual iterator end() = 0;

        
};

#endif 