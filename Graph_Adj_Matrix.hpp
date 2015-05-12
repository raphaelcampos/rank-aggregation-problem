#ifndef __GRAPH_ADJ_MATRIX__
#define __GRAPH_ADJ_MATRIX__

#include <string.h>
#include <iostream>
#include <iomanip>

#include <vector>

#include "IGraph.h"

using namespace std;

template<class T>
class Graph_Adj_Matrix : public IGraph<T> { 

    public:
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

        Graph_Adj_Matrix & operator=(const Graph_Adj_Matrix&);

    private:
        T ** adjM;
        T ** W;
        vector< pair<int, int> > indegree;
        vector< pair<int, int> > outdegree;
        int n;
        int m;

        template<class U>
        U** allocate(U **M, int size, int ini = 0){
            // If m != null deallocate it
            if(M != NULL) 
                this->deallocate(M, size);

            M = new U*[size];
            for (int i = 0; i < size; ++i)
            {
                M[i] = new U[size];
                memset(M[i], ini, sizeof(U)*size);
            }
            

            return M;
        }

        template<class U>
        U** copy(U **M, U **s, int size){

            U ** local = new U*[size];

            for (int i = 0; i < size; ++i)
            {
                local[i] = new U[size];
                memcpy(&(local[i][0]), &(s[i][0]), sizeof(U)*size);
            }

            // If m != null deallocate it
            if(M != NULL)
                deallocate(m, size);

            M = local;
            return M;
        }

        template<class U>
        void deallocate(U **M, int size){
            for (int i = 0; i < size; ++i)
            {
                delete [] M[i];
            }
        }

};

template <class T>
Graph_Adj_Matrix<T>::Graph_Adj_Matrix(){
    n = 0;
    m = 0;
    this->adjM = NULL;
    this->W = NULL;
}

template <class T>
Graph_Adj_Matrix<T>::Graph_Adj_Matrix(const Graph_Adj_Matrix<T> & g){

    this->n = g.n;
    this->m = g.m;

    adjM = copy(this->adjM, g.adjM);
    W = copy(this->W, g.W);
}

template <class T>
Graph_Adj_Matrix<T>::Graph_Adj_Matrix(int n, int ini){
    this->n = n;
    this->m = 0;
    this->adjM = NULL;
    this->W = NULL;
    
    adjM = allocate(this->adjM, n, ini);
    W = allocate(this->W, n, 0);

    outdegree.resize(n);
    indegree.resize(n);
}

template <class T>
Graph_Adj_Matrix<T>::~Graph_Adj_Matrix(){
    for (int i = 0; i < this->n; ++i)
    {
        delete [] adjM[i];
        delete [] W[i];
    }
}

template <class T>
void Graph_Adj_Matrix<T>::addVertex(int v){
    if(!vertexExists(v)){
        int tam = n + (v - n + 1);
        T **aux = new T*[tam];
        
        for (int i = 0; i < tam; ++i)
        {
            aux[i] = new T[tam];

            memcpy(&(aux[i][0]), &(adjM[i][0]), sizeof(T)*(n));
            aux[i][n] = 0;
        }
        
        // clear adjM
        for (int i = 0; i < n; ++i)
        {
            delete [] adjM[i];
        }

        adjM = aux;
        n = tam;
    }
}


template <class T>
void Graph_Adj_Matrix<T>::updateEdgeWeight(int v, int u, T w){
    if(edgeExists(v, u))
        W[v][u] = w;
}

template <class T>
void Graph_Adj_Matrix<T>::removeVertex(int v){

}

template <class T>
void Graph_Adj_Matrix<T>::addEdge(int v, int u, T w){
    if(edgeExists(u, v)) return updateEdgeWeight(u, v, w);


    if(vertexExists(v) && vertexExists(u)){
        W[v][u] = w;
        adjM[v][u] = 1;

        outdegree[v].first++;
        outdegree[v].second = v;

        indegree[u].first++;
        indegree[u].second = u;
        
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
    if(edgeExists(v, u))
        return W[v][u];
    else
        return 0;
}

template <class T>
bool Graph_Adj_Matrix<T>::edgeExists(int v, int u){
    return (n > v && n > u && adjM[u][v] != 0);
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
    if(this != &g){
        Graph_Adj_Matrix<T> tmp(g);

        if(this->adjM)
            for (int i = 0; i < this->n; ++i)
            {
                delete [] adjM[i];
            }

        this->adjM = tmp.adjM;
        this->n = tmp.n;
        this->m= tmp.m;
    }

    return *this;
}

template <class T>
void Graph_Adj_Matrix<T>::printAsMatrix(){
    cout << "asdasd" << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << setw(5) << adjM[i][j] << " ";  
        }

        cout << endl;
    }
}

#endif 