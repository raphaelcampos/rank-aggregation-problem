#ifndef __I_GRAPH__
#define __I_GRAPH__

using namespace std;

class IGraph {
    public:
        class vertex;
        typedef pair<double, vertex*> ve;

    public:
        class vertex_iterator_imp{
            public:

                virtual ~vertex_iterator_imp(){}
                
                virtual void next(){}

                virtual IGraph::vertex_iterator_imp* clone(){}
                
                virtual IGraph::vertex * get() const{}

                virtual bool isEqual(const IGraph::vertex_iterator_imp& other) const{}
        };

    public:
        class vertex
        {
            protected:
                class iterator_imp{
                    public:
                        virtual void next() {}
                        virtual iterator_imp* clone() {}
                        virtual ve * get() const {}
                        virtual bool isEqual(const iterator_imp& other) const {}
                };

            public:
                class iterator{
                    public:
                        typedef iterator self_type;
                        typedef ve value_type;
                        typedef value_type& reference;
                        typedef value_type* pointer;
                        typedef int difference_type;
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

                        iterator(iterator_imp * imp){
                            this->imp = imp;
                        }

                    private:
                        iterator_imp * imp;
                };

                virtual iterator begin() = 0;
                virtual iterator end() = 0;

                int indegree;
                int outdegree;
                int id;

        };


        class vertex_iterator{
            public:
                typedef vertex_iterator self_type;
                typedef vertex value_type;
                typedef value_type& reference;
                typedef value_type* pointer;
                typedef int difference_type;
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

                vertex_iterator(vertex_iterator_imp * imp){
                    this->imp = imp;
                }

            private:
                vertex_iterator_imp * imp;
        };

    public:
        virtual ~IGraph() {}
        
        virtual void addVertex(int v) = 0;
        /**
         * given the vertice id v, remove it from
         * the graph
         * @param v     vertice id
         */
        virtual void removeVertex(int v) = 0;
        
        /**
         * Add an new edge to the graph
         * @param v source vertex id
         * @param u destination vertex id
         * @param w edge weight
         */
        virtual void addEdge(int v, int u, double w) = 0;
        
        /**
         * Remove a vertex from the graph
         * @param v [description]
         * @param u [description]
         */
        virtual void removeEdge(int v, int u) = 0;
    
        /**
         * [getVertex description]
         * @param  v [description]
         * @return   [description]
         */
        virtual vertex * getVertex(int v) = 0;
        
        /**
         * [getEdge description]
         * @param  v [description]
         * @param  u [description]
         * @return   [description]
         */
        virtual double getEdge(int v, int u) = 0;
        
        /**
         * [updateEdgeWeight description]
         * @param v [description]
         * @param u [description]
         * @param w [description]
         */
        virtual void updateEdgeWeight(int v, int u, double w) = 0;
        
        /**
         * [edgeExists description]
         * @param  v [description]
         * @param  u [description]
         * @return   [description]
         */
        virtual bool edgeExists(int v, int u) = 0;
        
        /**
         * [vertexExists description]
         * @param  v [description]
         * @return   [description]
         */
        inline virtual bool vertexExists(int v) = 0;

        /**
         * 
         * @param v [description]
         */
        virtual void nextAdjVertex(int v) = 0;
        
        /**
         * true if the graph is complete, false otherwise
         * @return bool
         */
        virtual bool isComplete() = 0;
        
        /** 
         * returns total number of vertices
         * @return int
         */
        virtual int numVertices() = 0;

        /** 
         * returns total number of edges
         * @return int
         */
        virtual int numEdges() = 0;

        /**
         * iterator at first position
         * it allows iterations on 
         * the vertice's storage
         * @return IGraph::vertex_iterator
         */
        virtual vertex_iterator begin() = 0;

        /**
         * iterator for the end
         * it determine the end of iteration
         * @return IGraph::vertex_iterator
         */
        virtual vertex_iterator end() = 0;

    protected:    
        bool directed;
};

#endif 