#ifndef __SOCIAL_GRAPH_ADJ_MATRIX__
#define __SOCIAL_GRAPH_ADJ_MATRIX__

#include <string.h>
#include <iostream>
#include <iomanip>

#include <vector>
#include <set>

#include "Graph_Adj_Matrix2.hpp"
#include "IGraph.h"
#include "PriorityQueue.cpp"
#include "dlib/matrix.h"

using namespace std;
using namespace dlib;

class Social_Graph_Adj_Matrix : public Graph_Adj_Matrix {
    
    public:
        enum { GROUP_A, GROUP_B };

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
        
        void partion();
        void partionate();
        double calcCNormalizedVertex(int v);
        double calcCNormalized(int eAA, int eAB);

        double clusteringCoefficient(const set<int> &V);
        double modularity();
        double mixingTime(double epsilon = 10e-3);

        void putVertexInA(IGraph::vertex &v);
        inline double rank(int id);

        void loadSybils(string filename, bool isOrdered = true);
    
        std::vector<double> calcMetrics();

        std::vector<double> printMetrics();

    public:
        int *group;
        int *inGroup;
        int *outGroup;
        int eGroupA;
        int eGroupAB;
        set<int> groupA;
        set<int> groupB;
        set<int> sybils;
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
        groupB.insert(i);
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

    if(Graph_Adj_Matrix::edgeExists(v, u)){
        return;
    }

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
    double E = numEdges()/2;
    return -inGroup[id] + outGroup[id];
    if(outGroup[id] != 0 && inGroup[id] != 0){
        return (double)outGroup[id]/(double)inGroup[id];
        //return ((double)outGroup[id]/(double)inGroup[id]) - ((E - outGroup[id])/(E + (double)inGroup[id] -  outGroup[id]));
        
    }else if(outGroup[id] == 0){
        return -100000000;
    }else if(inGroup[id] == 0){
        return 10000000;
    }
}

/**
 * O(1)
 * @param  eAA [description]
 * @param  eAB [description]
 * @return      Normalized Condutancy
 */
double Social_Graph_Adj_Matrix::calcCNormalized(int eAA, int eAB){
    double E = numEdges()/2;
    return (eAA/(double)(eAA+eAB)) - ((E - eAA)/(double)(E + eAB));
}

/**
 * O(1)
 * @param  v  vertex id
 * @return     Normalized Condutancy
 */
double Social_Graph_Adj_Matrix::calcCNormalizedVertex(int v){
    int tmp_eAA = eGroupA;
    int tmp_eAB = eGroupAB;

    tmp_eAA += outGroup[v];
    tmp_eAB += inGroup[v] - outGroup[v];
   
    return calcCNormalized(tmp_eAA, tmp_eAB);
}

void Social_Graph_Adj_Matrix::partion(){
    typedef pair<double, IGraph::vertex*> item;
    IndexedPriorityQueue<double, std::greater<double> > Q(numVertices());
    groupB.clear();
    groupA.clear();

    // O(Vlg(V))
    for (IGraph::vertex_iterator it = begin(); it != end(); ++it)
    {
        int id = it->id;
        if(group[id] == GROUP_A){
            groupA.insert(it->id);
            groupB.insert(it->id);
        }else{
            groupB.insert(it->id); // log(v)
            Q.insert(id, Social_Graph_Adj_Matrix::calcCNormalizedVertex(id)); // log(v)
        }
    }

    //
    // V²logV + O(Vlog(V)) + O(E)
    //
    double CNant, CNdep;
    do{
        IGraph::vertex * v = getVertex(Q.minIndex());
        
        CNant = calcCNormalized(eGroupA, eGroupAB);
        CNdep = calcCNormalizedVertex(v->id);

        Q.deleteMin(); // O(log(v))

        if(group[v->id] == GROUP_B && CNant < CNdep){
            //cout << v->id << " : " << CNant << " < " << CNdep << endl;    
            putVertexInA(*v);
            groupA.insert(v->id); // O(log(v))
            
            for (int i = 0; i < numVertices(); ++i)
            {
                if(group[i] != GROUP_B) continue;
                v = getVertex(i);

                Q.changeKey(v->id, Social_Graph_Adj_Matrix::calcCNormalizedVertex(v->id));
            }
        }
    }while(!Q.isEmpty() && CNant < CNdep);
}

void Social_Graph_Adj_Matrix::partionate(){
    typedef pair<double, IGraph::vertex*> item;
    IndexedPriorityQueue<double, std::greater<double> > Q(numVertices());
    groupB.clear();
    groupA.clear();

    for (IGraph::vertex_iterator it = begin(); it != end(); ++it)
    {
        int id = it->id;
        if(group[id] == GROUP_A){
            groupA.insert(it->id);
            groupB.insert(it->id);
        }else{
            groupB.insert(it->id);
            Q.insert(id, Social_Graph_Adj_Matrix::rank(id));
        }
    }

    double CNant, CNdep;
    do{
        IGraph::vertex * v = getVertex(Q.minIndex());
        
        CNant = calcCNormalized(eGroupA, eGroupAB);
        CNdep = calcCNormalizedVertex(v->id);

        Q.deleteMin();

        if(group[v->id] == GROUP_B && CNant < CNdep){
            //cout << v->id << " : " << CNant << " < " << CNdep << endl;    
            putVertexInA(*v);
            groupA.insert(v->id);

            // spread the news to the neighbors
            for (IGraph::vertex::iterator e = v->begin(); e != v->end() ; ++e)
            {   
                IGraph::vertex * u = e->second;
                if(group[u->id] == GROUP_B){
                    Q.changeKey(u->id, Social_Graph_Adj_Matrix::rank(u->id));
                }
            }
            
        }
    }while(!Q.isEmpty() && CNant < CNdep);
    
    set<int> groupC;
    //sort(groupA.begin(), groupA.end());
    //sort(groupB.begin(), groupB.end());

    set_difference(groupB.begin(),groupB.end(),groupA.begin(),groupA.end(), inserter(groupC, groupC.begin()));

    //groupC.resize(it - groupC.begin());
    cout << "EM A " << groupA.size() << endl;
    cout << "EM B " << groupB.size() << endl;
    for (set<int>::iterator i = groupB.begin(); i != groupB.end(); ++i)
    {
        int id = (*i);
        //if(!id) continue; 
        //cout << id << endl;
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

void Social_Graph_Adj_Matrix::loadSybils(string filename, bool isOrdered){
    std::ifstream is(filename.c_str());
    std::string str;

    set<int>::iterator it = sybils.begin();
    while(std::getline(is,str)) {
        std::istringstream ss(str);
        int i;

        ss >> i;

        if(isOrdered){
            sybils.insert(it, i);
            it++;
        }else{
            sybils.insert(i);
        }
    }
}

std::vector<double>  Social_Graph_Adj_Matrix::calcMetrics(){
    set<int> honests;
    set<int> result;
    
    set_difference(groupB.begin(),groupB.end(),sybils.begin(),sybils.end(), inserter(honests, honests.begin()));

    result.clear();
    set_difference(groupB.begin(),groupB.end(),groupA.begin(),groupA.end(), inserter(result, result.begin()));

    groupB = result;
    //groupC.resize(it - groupC.begin());
    cout << "EM A " << groupA.size() << endl;
    cout << "EM B " << groupB.size() << endl;

    result.clear();
    set_intersection(groupB.begin(),groupB.end(),sybils.begin(),sybils.end(), inserter(result, result.begin()));

    double fscc = result.size()/(double)sybils.size();

    result.clear();
    set_intersection(groupA.begin(),groupA.end(), honests.begin(), honests.end(), inserter(result, result.begin()));
    double fhcc = result.size()/(double)honests.size();
    
    std::vector<double> metrics;

    metrics.push_back(numEdges()/((double)numVertices()));
    metrics.push_back(modularity());
    metrics.push_back(eGroupAB/(numEdges()/2.0 - eGroupA - eGroupAB));
    metrics.push_back((double)eGroupAB/eGroupA);
    metrics.push_back(clusteringCoefficient(groupB));
    metrics.push_back(clusteringCoefficient(groupA));
    metrics.push_back(fscc);
    metrics.push_back(fhcc);
    metrics.push_back(1 - fhcc);
    metrics.push_back(1 - fscc);
    metrics.push_back(mixingTime());

    return metrics; 
}

std::vector<double> Social_Graph_Adj_Matrix::printMetrics(){

    std::vector<double> metrics = calcMetrics();
    
    cout << "(a) Grau médio : " << metrics[0] << endl;
    cout << "(b) Modularidade : " << metrics[1] << endl;
    cout << "(c) Condutancia Sybil : " << metrics[2] << endl;
    cout << "(d) Condutancia Honesta : " << metrics[3] << endl;
    cout << "(e) Coef agrup Sybil : " << metrics[4] << endl;
    cout << "(f) Coef agrup Honesta : " << metrics[5] << endl;
    cout << "(g) Fracao sybil corretamente classif : " << metrics[6] << endl;
    cout << "(h) Fracao honestos corretamente classif : " << metrics[7] << endl;
    cout << "(i) Fracao falso positivos : " << metrics[8] << endl;
    cout << "(j) Fracao falso negativos : " << metrics[9] << endl;
    
    cout << "(k) Mixing time : " << metrics[10] << endl;

    return metrics;
}

/**
 * O(V³)
 * @param  V [description]
 * @return   [description]
 */
double Social_Graph_Adj_Matrix::clusteringCoefficient(const set<int> &V){
    double sum = 0;
    int total = 0;
    for (set<int>::iterator i = V.begin(); i != V.end(); ++i)
    {
        int numEdges = 0;
        double coef = 0;

        for (set<int>::iterator j = V.begin(); j != V.end(); ++j)
        {
            if(i == j) continue;
            if(!edgeExists(*i,*j)) continue;
            for (set<int>::iterator k = V.begin(); k != V.end(); ++k)
            {
                if(edgeExists(*j, *k) && edgeExists(*k, *i)) numEdges++;
            }
        }
        total += inGroup[*i];
        coef = 2.0*numEdges/(double)(inGroup[*i]*(inGroup[*i]-1));
        sum += coef;
    }
    
    return sum/V.size();
}

double Social_Graph_Adj_Matrix::modularity(){
    double mA, mB, ls, ds, L = numEdges()/2.0;
    
    ls = eGroupA; // number of edges in A
    ds = (2.0*eGroupA) + eGroupAB;
    mA = (ls/L) - (ds/(2.0*L))*(ds/(2.0*L));
    
    ls = L - eGroupA - eGroupAB; // number of edges in B
    ds = 2*ls + eGroupAB;
    mB = ((ls)/L) - (ds/(2.0*L))*(ds/(2.0*L));
    
    return mA + mB;
}

double Social_Graph_Adj_Matrix::mixingTime(double epsilon){
    typedef std::map<unsigned long, double> sample;
    int n = this->numVertices();
    matrix<double> P = zeros_matrix<double>(n,n);
    matrix<double> ds(1, n);
    matrix<double> iD(1, n);
    matrix<double> delta = ones_matrix<double>(1, n);
    bool *seen = new bool[n];
    
    for (IGraph::vertex_iterator u = this->begin(); u != this->end(); ++u)
    {

        ds(0, u->id) = u->outdegree/((double)numEdges());
        iD(0, u->id) = 1.0/n;
        seen[u->id] = false;
        sample samp;
        for (IGraph::vertex::iterator e = u->begin(); e != u->end(); ++e)
        {
            IGraph::vertex * v = e->second;
            P(u->id, v->id) = 1.0/(u->outdegree);
        }
    }

    matrix<double> Pt = P;

    int t = 0;
    
    double maxi = 10;
    
    /*matrix<double> eigenvector = real_eigenvalues(P);
    std::sort(&eigenvector(0), &eigenvector(0)+eigenvector.size());
    double mi = max(abs(eigenvector(1)), abs(eigenvector(eigenvector.size()-2)));
    cout << "MI : " << mi << endl;
    cout << "upper bound : " << (log(n) + log(1/epsilon))/(1-mi) << endl;
    cout << "lower bound : " << (mi/(2-2*mi))*log(1/(2*epsilon)) << endl;*/
    
    while(maxi > epsilon){
        
        double max_dist = 0;
        for (int i = 0; i < n; ++i)
        {
            double dist = sum(abs(rowm(Pt, i) - ds))/2;
            delta(0, i) = dist;
            iD(0, i) = t;
            
            if(delta(0, i) > max_dist){
                max_dist = delta(0, i);
            }
        }
        maxi = max_dist;
        
        Pt = Pt * P;

        t++;
    }

    return max(iD);

}

#endif 