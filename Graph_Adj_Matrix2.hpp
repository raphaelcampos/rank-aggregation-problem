#ifndef __GRAPH_ADJ_MATRIX__
#define __GRAPH_ADJ_MATRIX__

#include <string.h>
#include <iostream>
#include <iomanip>

#include <vector>

#include "IGraph.h"

using namespace std;

class Graph_Adj_Matrix : public IGraph { 
    public:
        class vertex : public IGraph::vertex
        {
            friend Graph_Adj_Matrix;
            typedef pair<double, IGraph::vertex*> ve;

            protected:
                class iterator_imp : public IGraph::vertex::iterator_imp{
                    public:
                        void next(){
                            if(ptr_ == end_)
                                return;

                            ptr_++;
                            if(ptr_->second == NULL && ptr_ != end_)
                                next();
                        }
                        
                        IGraph::vertex::iterator_imp* clone(){
                            return this;
                        }
                        
                        IGraph::ve * get() const{
                            IGraph::ve * v = &(*ptr_);

                            return v;
                        }

                        bool isEqual(const IGraph::vertex::iterator_imp& other) const {
                            return get() == other.get();
                        }

                        iterator_imp(vector< Graph_Adj_Matrix::vertex::ve >::iterator ptr, vector< Graph_Adj_Matrix::vertex::ve >::iterator end){
                            ptr_ = ptr;
                            end_ = end;
                        }

                    private:
                        vector< Graph_Adj_Matrix::vertex::ve >::iterator ptr_;
                        vector< Graph_Adj_Matrix::vertex::ve >::iterator end_;
                };


            public:
                vertex(){
                    indegree = 0;
                    outdegree = 0;
                }
                vertex(int n){
                    indegree = 0;
                    outdegree = 0;
                    for (int j = 0; j < n; ++j)
                    {
                        pair<double, vertex*> adj;
                        adj.first = 0;
                        adj.second = NULL;
                        
                        this->adj.push_back(adj);
                    }
                }

                IGraph::vertex::iterator begin(){ 
                    iterator_imp * first = new iterator_imp(this->adj.begin(), this->adj.end());
                    
                    if(first->get()->second == NULL)
                        first->next();

                    IGraph::vertex::iterator f(first);
                    return f;
                }

                IGraph::vertex::iterator end(){
                    iterator_imp * last = new iterator_imp(this->adj.end(), this->adj.end());
                    IGraph::vertex::iterator l(last);
                    return l;
                }

            public:
                vector<ve> adj;
         };

        Graph_Adj_Matrix();
        Graph_Adj_Matrix(const Graph_Adj_Matrix &);
        Graph_Adj_Matrix(int n, int ini = 0);
        
        ~Graph_Adj_Matrix();
        void addVertex(int v);
        void removeVertex(int v);
        void addEdge(int v, int u, double w);
        void removeEdge(int v, int u);
        IGraph::vertex * getVertex(int v);
        double getEdge(int v, int u);
        void updateEdgeWeight(int v, int u, double w);
        bool edgeExists(int v, int u);
        bool vertexExists(int v);
        void nextAdjVertex(int v);
        bool isComplete();
        int numVertices();
        void printAsMatrix();

        IGraph::vertex_iterator begin();

        IGraph::vertex_iterator end();

        Graph_Adj_Matrix & operator=(const Graph_Adj_Matrix&);

    public:
        class vertex_iterator_imp : public IGraph::vertex_iterator_imp{
            public:
                void next(){
                    ptr_++;
                }
                
                IGraph::vertex_iterator_imp * clone(){
                    return this;
                }
                
                IGraph::vertex * get() const{  
                    return &(*ptr_);
                }

                bool isEqual(const IGraph::vertex_iterator_imp& other) const{
                    return get() == other.get();
                }

                vertex_iterator_imp(vector< Graph_Adj_Matrix::vertex >::iterator ptr){
                    this->ptr_ = ptr;
                }

            private:
                vector< Graph_Adj_Matrix::vertex >::iterator ptr_;
        };

    private:
        vector< Graph_Adj_Matrix::vertex > vertices;
        int n;
        int m;
};

Graph_Adj_Matrix::Graph_Adj_Matrix(){
    n = vertices.size();
    m = 0;
}

Graph_Adj_Matrix::Graph_Adj_Matrix(const Graph_Adj_Matrix & g){

    this->n = g.n;
    this->m = g.m;
}

Graph_Adj_Matrix::Graph_Adj_Matrix(int n, int ini){
    this->n = n;
    
    for (int i = 0; i < n; ++i)
    {
        Graph_Adj_Matrix::vertex v(n);
        v.id = i;
        vertices.push_back(v);
    }
}

Graph_Adj_Matrix::~Graph_Adj_Matrix(){
}

void Graph_Adj_Matrix::addVertex(int v){

}

void Graph_Adj_Matrix::removeEdge(int v, int u){

}


void Graph_Adj_Matrix::updateEdgeWeight(int v, int u, double w){
    if(edgeExists(v, u))
        vertices[v].adj[u].first = w;
}

void Graph_Adj_Matrix::addEdge(int v, int u, double w){
    if(edgeExists(v, u)) return updateEdgeWeight(v, u, w);

    if(vertexExists(v) && vertexExists(u)){
        vertices[v].adj[u].first = w;
        vertices[v].adj[u].second = &vertices[u];

        vertices[v].outdegree++;
        
        vertices[v].adj[u].second->indegree++;
        
        m++;
    }
}

void Graph_Adj_Matrix::removeVertex(int v){
    if(!vertexExists(v)) return;

    vertices[v].outdegree = 0;
    vertices[v].indegree = 0;
            
    for (int u = 0; u < n; ++u)
    {
        if(edgeExists(v, u)){
            vertices[v].adj[u].second->indegree--;
        }

        if(edgeExists(u, v)){
            vertices[u].outdegree--; 
            vertices[u].adj[v].first = 0;
            vertices[u].adj[v].second = NULL;
        }
    }
}

IGraph::vertex * Graph_Adj_Matrix::getVertex(int v){
    if(vertexExists(v))
        return &vertices[v];
    return NULL;
}

double Graph_Adj_Matrix::getEdge(int v, int u){
    return vertices[v].adj[u].first;
}

bool Graph_Adj_Matrix::edgeExists(int v, int u){

    return (n > v && n > u && vertices[v].adj[u].second != NULL);
}

bool Graph_Adj_Matrix::vertexExists(int v){
    return (n > v);
}

void Graph_Adj_Matrix::nextAdjVertex(int v){

}

bool Graph_Adj_Matrix::isComplete(){
    return m >= (n)*(n-1);
}

int Graph_Adj_Matrix::numVertices(){
    return vertices.size();
}

Graph_Adj_Matrix & Graph_Adj_Matrix::operator=(const Graph_Adj_Matrix& g){
}

void Graph_Adj_Matrix::printAsMatrix(){
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << setw(5) << vertices[i].adj[j].first << " ";  
        }

        cout << endl;
    }
}

IGraph::vertex_iterator  Graph_Adj_Matrix::begin(){
    vertex_iterator_imp * first = new vertex_iterator_imp(this->vertices.begin());
    IGraph::vertex_iterator_imp * ff = first;
    first->get();
    ff->get();
    IGraph::vertex_iterator f(first);

    return f;
}

IGraph::vertex_iterator  Graph_Adj_Matrix::end(){

    Graph_Adj_Matrix::vertex_iterator_imp * last = new Graph_Adj_Matrix::vertex_iterator_imp(this->vertices.end());

    IGraph::vertex_iterator l(last);

    return l;
}

#endif 