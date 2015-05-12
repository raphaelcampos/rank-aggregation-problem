#ifndef __GRAPH_ADJ_MATRIX__
#define __GRAPH_ADJ_MATRIX__

#include <string.h>
#include <iostream>
#include <iomanip>

#include <vector>

#include "IGraph.h"

using namespace std;

template<class T>
class Graph_Adj_Matrix;

template<class U>
class vertex
{
    friend Graph_Adj_Matrix<U>;
    typedef pair<U, vertex*> ve;

    public:
        class iterator
        {
            public:
                typedef iterator self_type;
                typedef ve value_type;
                typedef ve& reference;
                typedef ve* pointer;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                iterator(pointer ptr) : ptr_(ptr) { }
                self_type operator++() { self_type i = *this; ptr_++; return i; }
                self_type operator++(int junk) { ptr_++; return *this; }
                reference operator*() { return ptr_->vertices[pos]; }
                pointer operator->() { return &(ptr_->vertices[pos]); }
                bool operator==(const self_type& rhs) { return ptr_ == rhs.ptr_; }
                bool operator!=(const self_type& rhs) { return ptr_ != rhs.ptr_; }
            private:
                pointer ptr_;
                int pos;
        };
        vertex(){
            indegree = 0;
            outdegree = 0;
        }


    public:
        vector<ve> adj;
        int indegree;
        int outdegree;
 };

template<class T>
class Graph_Adj_Matrix : public IGraph<T> { 
    public:
    
        class vertex_iterator : public forward_iterator_tag
        {
            public:
                typedef vertex_iterator self_type;
                typedef vertex<T> value_type;
                typedef value_type& reference;
                typedef value_type* pointer;
                typedef std::forward_iterator_tag iterator_category;
                typedef int difference_type;
                vertex_iterator(vector< vertex<T> > *ptr, int pos = 0) : ptr_(ptr) { this->pos = pos; }
                self_type operator++() { self_type i = *this; pos++; return i; }
                self_type operator++(int junk) { pos++; return *this; }
                reference operator*() { return (*ptr_)[pos]; }
                pointer operator->() { return &(*ptr_)[pos]; }
                bool operator==(const self_type& rhs) { return pos == rhs.pos; }
                bool operator!=(const self_type& rhs) { return pos != rhs.pos; }
            private:
                vector< vertex<T> > * ptr_;
                int pos;
        };

        Graph_Adj_Matrix();
        Graph_Adj_Matrix(const Graph_Adj_Matrix &);
        Graph_Adj_Matrix(int n, int ini = 0);
        
        ~Graph_Adj_Matrix();
        void addVertex(int v);
        void removeVertex(int v);
        void addEdge(int v, int u, T w);
        void removeVertex(int v, int u);
        void inDegree(int v);
        void outDegree(int v);
        void getVertex(int v);
        T getEdge(int v, int u);
        void updateEdgeWeight(int v, int u, T w);
        bool edgeExists(int v, int u);
        bool vertexExists(int v);
        void nextAdjVertex(int v);
        bool isComplete();
        void printAsMatrix();

        Graph_Adj_Matrix<T>::vertex_iterator begin(){
            vertex_iterator first(&this->vertices);
            return first;
        }

        Graph_Adj_Matrix<T>::vertex_iterator end(){
            vertex_iterator last(&this->vertices, this->vertices.size()-1);
            return last;
        }

        Graph_Adj_Matrix & operator=(const Graph_Adj_Matrix&);

    private:
        vector< vertex<T> > vertices;
        int n;
        int m;
};

template <class T>
Graph_Adj_Matrix<T>::Graph_Adj_Matrix(){
    n = vertices.size();
    m = 0;
}

template <class T>
Graph_Adj_Matrix<T>::Graph_Adj_Matrix(const Graph_Adj_Matrix<T> & g){

    this->n = g.n;
    this->m = g.m;
}

template <class T>
Graph_Adj_Matrix<T>::Graph_Adj_Matrix(int n, int ini){
    this->n = n;
    
    for (int i = 0; i < n; ++i)
    {
        vertex<T> v;
        for (int j = 0; j < n; ++j)
        {
            
            pair<T, vertex<T>*> adj;
            adj.first = 0;
            adj.second = NULL;
            
            v.adj.push_back(adj);
        }

        vertices.push_back(v);
    }
}

template <class T>
Graph_Adj_Matrix<T>::~Graph_Adj_Matrix(){
}

template <class T>
void Graph_Adj_Matrix<T>::addVertex(int v){

}


template <class T>
void Graph_Adj_Matrix<T>::updateEdgeWeight(int v, int u, T w){
    if(edgeExists(v, u))
        vertices[v].adj[u].first = w;
}

template <class T>
void Graph_Adj_Matrix<T>::removeVertex(int v){

}

template <class T>
void Graph_Adj_Matrix<T>::addEdge(int v, int u, T w){
    if(edgeExists(v, u)) return updateEdgeWeight(v, u, w);

    if(vertexExists(v) && vertexExists(u)){
        vertices[v].adj[u].first = w;
        vertices[v].adj[u].second = &vertices[u];

        vertices[v].outdegree++;
        
        vertices[v].adj[u].second->indegree++;
        
        m++;
    }
}

template <class T>
void Graph_Adj_Matrix<T>::removeVertex(int v, int u){

}

template <class T>
void Graph_Adj_Matrix<T>::inDegree(int v){

}

template <class T>
void Graph_Adj_Matrix<T>::outDegree(int v){

}

template <class T>
void Graph_Adj_Matrix<T>::getVertex(int v){

}

template <class T>
T Graph_Adj_Matrix<T>::getEdge(int v, int u){
    return vertices[v].adj[u].first;
}

template <class T>
bool Graph_Adj_Matrix<T>::edgeExists(int v, int u){

    return (n > v && n > u && vertices[v].adj[u].second != NULL);
}

template <class T>
bool Graph_Adj_Matrix<T>::vertexExists(int v){
    return (n > v);
}

template <class T>
void Graph_Adj_Matrix<T>::nextAdjVertex(int v){

}

template <class T>
bool Graph_Adj_Matrix<T>::isComplete(){
    return m >= n*(n-1)/2;
}

template <class T>
Graph_Adj_Matrix<T> & Graph_Adj_Matrix<T>::operator=(const Graph_Adj_Matrix<T>& g){
}

template <class T>
void Graph_Adj_Matrix<T>::printAsMatrix(){
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << setw(5) << vertices[i].adj[j].first << " ";  
        }

        cout << endl;
    }
}

#endif 