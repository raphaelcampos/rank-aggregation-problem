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
                virtual void next() = 0;
                virtual vertex_iterator_imp* clone() = 0;
                virtual IGraph::vertex * get() const = 0;
                virtual bool isEqual(const vertex_iterator_imp& other) const = 0;
        };

    public:
        class vertex
        {
            protected:
                class iterator_imp{
                    public:
                         void next();
                         iterator_imp* clone();
                         ve * get() const;
                         bool isEqual(const iterator_imp& other) const;
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
                            iterator_imp i = *rhs.imp;
                            return imp->isEqual(i);
                        }
                        bool operator!=(const self_type& rhs){
                             iterator_imp i = *rhs.imp;
                            return !imp->isEqual(i);
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
                    //vertex_iterator_imp * i = *rhs.imp;
                    return imp->isEqual(*rhs.imp);
                }
                bool operator!=(const self_type& rhs){
                    //vertex_iterator_imp i = *rhs.imp;
                    return !imp->isEqual(*rhs.imp);
                }

                vertex_iterator(vertex_iterator_imp * imp){
                    //cout << imp->get()->id << endl;
                    this->imp = imp;
                }

            private:
                vertex_iterator_imp * imp;
        };

    public:
        virtual ~IGraph() {}
        virtual void addVertex(int v) = 0;
        virtual void removeVertex(int v) = 0;
        virtual void addEdge(int v, int u, double w) = 0;
        virtual void removeVertex(int v, int u) = 0;
        virtual void inDegree(int v) = 0;
        virtual void outDegree(int v) = 0;
        virtual void getVertex(int v) = 0;
        virtual double getEdge(int v, int u) = 0;
        virtual void updateEdgeWeight(int v, int u, double w) = 0;
        virtual bool edgeExists(int v, int u) = 0;
        virtual bool vertexExists(int v) = 0;
        virtual void nextAdjVertex(int v) = 0;
        virtual bool isComplete() = 0;
        virtual int numVertices() = 0;

        virtual vertex_iterator begin() = 0;
        virtual vertex_iterator end() = 0;
};

#endif 