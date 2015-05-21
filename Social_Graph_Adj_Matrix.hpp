#ifndef __SOCIAL_GRAPH_ADJ_MATRIX__
#define __SOCIAL_GRAPH_ADJ_MATRIX__

#include <string.h>
#include <iostream>
#include <iomanip>

#include <vector>

#include "Graph_Adj_Matrix2.hpp"
#include "IGraph.h"

using namespace std;

class Social_Graph_Adj_Matrix : public Graph_Adj_Matrix {
    
    public:
        enum { GROUP_A, GROUP_B };

        class vertex : public Graph_Adj_Matrix::vertex
        {
            //friend Graph_Adj_Matrix;
            typedef pair<double, IGraph::vertex*> ve;

            public:
                vertex()
                :Graph_Adj_Matrix::vertex(){
                    init();
                }
                vertex(int n)
                :Graph_Adj_Matrix::vertex(){
                    init();
                    for (int j = 0; j < n; ++j)
                    {
                        pair<double, vertex*> adj;
                        adj.first = 0;
                        adj.second = NULL;
                        
                        this->adj.push_back(adj);
                    }
                }

            private:
                void init(){
                    inGroup = 0;
                    outGroup = 0;
                    group = GROUP_B;
                }

            public:
                vector<ve> adj;
                int inGroup;
                int outGroup;
                int group;
         };

        Social_Graph_Adj_Matrix();
        Social_Graph_Adj_Matrix(const Social_Graph_Adj_Matrix &);
        Social_Graph_Adj_Matrix(int n, int ini = 0);
        
        ~Social_Graph_Adj_Matrix();

        void addVertex(int v);
        void removeVertex(int v);
        
        void addEdge(int v, int u, double w);
        void removeEdge(int v, int u);
        void updateEdgeWeight(int v, int u, double w);
        
        IGraph::vertex * getVertex(int v);
        double getEdge(int v, int u);
        
        void partionate();
        double calcCNormalized(int eAA, int eAB);
        void putVertexInA(IGraph::vertex &v);
        inline double rank(int id);
    
    public:
        int *group;
        int *inGroup;
        int *outGroup;
        int eGroupA;
        int eGroupAB;
        vector<IGraph::vertex*> groupA;
        vector<IGraph::vertex*> groupB;
};

Social_Graph_Adj_Matrix::Social_Graph_Adj_Matrix()
:Graph_Adj_Matrix(false){
    eGroupAB = 0;
    eGroupA = 0;
}

Social_Graph_Adj_Matrix::Social_Graph_Adj_Matrix(const Social_Graph_Adj_Matrix & g)
:Graph_Adj_Matrix(g){
    eGroupAB = 0;
    eGroupA = 0;
}

Social_Graph_Adj_Matrix::Social_Graph_Adj_Matrix(int n, int ini)
:Graph_Adj_Matrix(n, ini, false){
    eGroupAB = 0;
    eGroupA = 0;
    inGroup = new int[numVertices()];
    outGroup = new int[numVertices()];
    group = new int[numVertices()];

    for (int i = 0; i < numVertices(); ++i)
    {
        group[i] = GROUP_B;
        groupB.push_back(getVertex(i));
    }

    memset(outGroup, 0, sizeof(int)*numVertices());
    memset(inGroup, 0, sizeof(int)*numVertices());
}

Social_Graph_Adj_Matrix::~Social_Graph_Adj_Matrix(){
    delete [] group;
    delete [] inGroup;
    delete [] outGroup;
}

void Social_Graph_Adj_Matrix::addVertex(int v){
    Graph_Adj_Matrix::addVertex(v);
}

void Social_Graph_Adj_Matrix::removeEdge(int v, int u){
    Graph_Adj_Matrix::removeEdge(v, u);
}


void Social_Graph_Adj_Matrix::updateEdgeWeight(int v, int u, double w){
    Graph_Adj_Matrix::updateEdgeWeight(v, u, w);
}

void Social_Graph_Adj_Matrix::addEdge(int v, int u, double w){
    IGraph::vertex *sv = (IGraph::vertex*)getVertex(v);
    IGraph::vertex *su = (IGraph::vertex*)getVertex(u);

    if(edgeExists(v, u)) return;

    if(v == 1 || u == 1)
        cout << v << "," << u << endl;

    if(group[sv->id] == group[su->id]){
        inGroup[sv->id]++;
        inGroup[su->id]++;
        if(group[sv->id] == GROUP_A){
            eGroupA++;
        }
    }else{
        outGroup[sv->id]++;
        outGroup[su->id]++;
        eGroupAB++;
    }

    Graph_Adj_Matrix::addEdge(v, u, w);
}

/**
 * Puts vertex u in the group A
 * O(n) where n is the number of vertices adj of u 
 * @param u vertex
 */
void Social_Graph_Adj_Matrix::putVertexInA(IGraph::vertex &u){
    // already in A
    if(group[u.id] == GROUP_A) return;

    // u becomes A
    group[u.id] = GROUP_A;

    //update eAA e eAB
    eGroupA += outGroup[u.id];
    eGroupAB += inGroup[u.id] - outGroup[u.id];
    
    // swap
    int aux = outGroup[u.id];
    outGroup[u.id] = inGroup[u.id];
    inGroup[u.id] = aux;
    
    // spread the news to the neighbors
    for (IGraph::vertex::iterator e = u.begin(); e != u.end() ; ++e)
    {   
        IGraph::vertex * v = e->second;
        if(group[v->id] == GROUP_A){
            inGroup[v->id]++; 
            outGroup[v->id]--;
        }else{
            inGroup[v->id]--;
            outGroup[v->id]++;
        }
    }
}

inline double Social_Graph_Adj_Matrix::rank(int id){

    /*int tmp_eAA = eGroupA;
    int tmp_eAB = eGroupAB;        

    tmp_eAA += outGroup[id];
    tmp_eAB += inGroup[id] - outGroup[id];
   
    return calcCNormalized(tmp_eAA, tmp_eAB);*/
    //return (double)outGroup[id]/(double)(inGroup[id]+1);
    /*if(outGroup[id] != 0 && inGroup[id] != 0){
        return (double)outGroup[id]/(double)inGroup[id];
    }else if(outGroup[id] == 0){
        return -100000000;
    }else if(inGroup[id] == 0){
        return 10000000;
    }*/
}

/**
 * O(1)
 * @param  eAA [description]
 * @param  eAB [description]
 * @return     [description]
 */
double Social_Graph_Adj_Matrix::calcCNormalized(int eAA, int eAB){
    double E = numEdges()/2;
    return (eAA/(double)(eAA+eAB)) - ((E - eAA)/(double)(E + eAB));
}

void Social_Graph_Adj_Matrix::partionate(){
    typedef pair<double, IGraph::vertex*> item;
    priority_queue<item> Q;
    groupB.clear();
    groupA.clear();

    for (IGraph::vertex_iterator it = begin(); it != end(); ++it)
    {
        int id = it->id;
        if(group[id] == GROUP_A){
            groupA.push_back(&(*it));
        }else{
            groupB.push_back(&(*it));
            Q.push(make_pair(Social_Graph_Adj_Matrix::rank(id), &(*it)));
        }
    }

    for (int i = 0; i < groupB.size(); ++i)
    {
        int id = groupB[i]->id;
        cout << id << " : " << Social_Graph_Adj_Matrix::rank(id) << endl;
    }


    double CNant, CNdep;
    do{
        int tmp_eAA = eGroupA;
        int tmp_eAB = eGroupAB;
        CNant = calcCNormalized(eGroupA, eGroupAB);
        IGraph::vertex * v = Q.top().second;        

        tmp_eAA += outGroup[v->id];
        tmp_eAB += inGroup[v->id] - outGroup[v->id];
       
        CNdep = calcCNormalized(tmp_eAA, tmp_eAB);

        Q.pop();

        if(group[v->id] == GROUP_B && CNant < CNdep){
            cout << v->id << " : " << CNant << " < " << CNdep << endl;    
            putVertexInA(*v);
            groupA.push_back(&(*v));
            
            // spread the news to the neighbors
            for (IGraph::vertex::iterator e = v->begin(); e != v->end() ; ++e)
            {   
                IGraph::vertex * u = e->second;
                if(group[u->id] == GROUP_B){
                    Q.push(make_pair(Social_Graph_Adj_Matrix::rank(u->id), &(*u)));
                }
            }
        }
    }while(!Q.empty() && CNant < CNdep);

    cout << "EM A" << endl;
    for (int i = 0; i < groupA.size(); ++i)
    {
        int id = groupA[i]->id;
        cout << id << " : " << Social_Graph_Adj_Matrix::rank(id) << endl;
    }

}

void Social_Graph_Adj_Matrix::removeVertex(int v){
    Graph_Adj_Matrix::removeVertex(v);
}

IGraph::vertex * Social_Graph_Adj_Matrix::getVertex(int v){
    return Graph_Adj_Matrix::getVertex(v);
}

double Social_Graph_Adj_Matrix::getEdge(int v, int u){
    return Graph_Adj_Matrix::getEdge(v, u);
}
#endif 