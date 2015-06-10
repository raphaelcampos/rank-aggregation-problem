
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <string.h>
#include <time.h>
#include <cmath>

#include <algorithm>

#include <queue>
#include <vector>
#include <map>

#include "IGraph.h"
#include "Graph_Adj_Matrix2.hpp"
#include "Social_Graph_Adj_Matrix.hpp"

using namespace std;

typedef std::vector<int> nl_t;
typedef std::vector<nl_t> nll_t;

void print_array(int a[], int n){
	for (int l = 0; l < n; ++l)
	{
		cout << a[l] << ",";
	}
	cout << endl;
	
}

Social_Graph_Adj_Matrix * load_from_file(string filename){
	std::ifstream is(filename.c_str());
    std::string str;

    nll_t nll;

	while(std::getline(is,str)) {
        std::istringstream ss(str);
        nl_t nl;
        int i;

        if(ss >> i){
        	nl.push_back(i);
	        while(ss >> i) {
	        	nl.push_back(i);
	        }

	        nll.push_back(nl);
	    }
    }

    Social_Graph_Adj_Matrix *G = new Social_Graph_Adj_Matrix(nll.size());
    
    int E = 0;
    for (int i = 0; i < nll.size(); ++i)
    {
    	for (int j = 1; j < nll[i].size(); ++j)
    	{
    		G->addEdge(nll[i][0], nll[i][j], 1);
    		E++;
    	}
    }

    cout << "E = "  << E << endl; 

    return G;

}

template<
	typename storage_type
>
void save_region(string filename, storage_type V){
	ofstream out;
	out.open(filename.c_str());
	for (typename storage_type::iterator i = V.begin(); i != V.end(); ++i)
    {
    	out << *i << endl;
    }
    out.close();
}

void metrics(int argc, char const *argv[]){
	Social_Graph_Adj_Matrix *graphA = load_from_file(argv[1]);
	Social_Graph_Adj_Matrix *graphB = load_from_file(argv[3]);
	

	graphA->loadSybils(argv[2]);
	
	set_difference(graphA->groupB.begin(), graphA->groupB.end(), graphA->sybils.begin(), graphA->sybils.end(), inserter(graphA->groupA, graphA->groupA.begin()));

    set_difference(graphA->groupA.begin(), graphA->groupA.end(), graphA->sybils.begin(), graphA->sybils.end(), inserter(graphA->groupB, graphA->groupB.begin()));

    for (set<int>::iterator i = graphA->groupA.begin(); i != graphA->groupA.end(); ++i)
    {
    	graphA->putVertexInA(*graphA->getVertex(*i));
    }

	std::vector<double> metrics = graphA->printMetrics();

	graphB->loadSybils(argv[4]);
	
	
    set_difference(graphB->groupB.begin(), graphB->groupB.end(), graphB->sybils.begin(), graphB->sybils.end(), inserter(graphB->groupA, graphB->groupA.begin()));

    set_difference(graphB->groupA.begin(), graphB->groupA.end(), graphB->sybils.begin(), graphB->sybils.end(), inserter(graphB->groupB, graphB->groupB.begin()));

    for (set<int>::iterator i = graphB->groupA.begin(); i != graphB->groupA.end(); ++i)
    {
    	graphB->putVertexInA(*graphB->getVertex(*i));
    }

	metrics = graphB->printMetrics();
}

int main(int argc, char const *argv[])
{

	Social_Graph_Adj_Matrix *graphA = load_from_file(argv[1]);
	Social_Graph_Adj_Matrix *graphB = load_from_file(argv[3]);
	
	int ids[] = {5,7,10,13,18,20,23,35,37,45,49,53,59,67,79,82,90,92,96,99};
	for (int i = 0; i < 20; ++i)
	{
		graphA->putVertexInA(*graphA->getVertex(ids[i]));
		graphB->putVertexInA(*graphB->getVertex(ids[i]));
	}

	graphA->loadSybils(argv[2]);
	graphA->partion();
	std::vector<double> metrics = graphA->printMetrics();
	save_region("regionHonestGA.txt", graphA->groupA);
	save_region("regionSybilGA.txt", graphA->groupB);
	save_region("metricsGA.txt", metrics);

	delete graphA;
	graphB->loadSybils(argv[4]);
	graphB->partion();
	metrics = graphB->printMetrics();
	save_region("regionHonestGB.txt", graphB->groupA);
	save_region("regionSybilGB.txt", graphB->groupB);
	save_region("metricsGB.txt", metrics);
}