
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
    
    for (int i = 0; i < nll.size(); ++i)
    {
    	for (int j = 1; j < nll[i].size(); ++j)
    	{
    		G->addEdge(nll[i][0], nll[i][j], 1);
    	}
    }

    return G;

}

int main(int argc, char const *argv[])
{
	char * f = "instancia/inGA.txt";
	Social_Graph_Adj_Matrix *graph = load_from_file(argv[1]);
	
	int ids[] = {5,7,10,13,18,20,23,35,37,45,49,53,59,67,79,82,90,92,96,99};
	for (int i = 0; i < 20; ++i)
	{
		//graph->putVertexInA(*graph->getVertex(rand()%100));
		graph->putVertexInA(*graph->getVertex(ids[i]));
	}
	graph->partion();
	graph->loadSybils(argv[2]);

	//delete [] graph;

	Social_Graph_Adj_Matrix *graph1 = load_from_file(argv[1]);
	srand(time(NULL));
	//int ids[] = {5,7,10,13,18,20,23,35,37,45,49,53,59,67,79,82,90,92,96,99};
	for (int i = 0; i < 20; ++i)
	{
		//graph->putVertexInA(*graph->getVertex(rand()%100));
		graph1->putVertexInA(*graph1->getVertex(ids[i]));
	}
	graph1->partionate();

	cout << "Grau médio : " << graph->numEdges()/(2.0*graph->numVertices()) << endl;

	//delete [] graph1;
}