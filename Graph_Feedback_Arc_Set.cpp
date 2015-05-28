#include <iostream>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <vector>

#include "IGraph.h"

using namespace std;

namespace graph{

	void FAS_pivot(IGraph &graph, vector<int> &vertices){
		
		if(vertices.size() <= 1) return;
		// pick random pivot
		int i = vertices[rand()%(vertices.size())];
		
		vector<int> vL;
		vector<int> vR;
		
		for (int j = 0; j < vertices.size(); ++j)
		{
			if(i == vertices[j]) continue;

			if(graph.edgeExists(vertices[j], i)){
				vL.push_back(vertices[j]);
			}else if(graph.edgeExists(i, vertices[j])){
				vR.push_back(vertices[j]);
			}
		}
		
		FAS_pivot(graph, vL);
		FAS_pivot(graph, vR);
		
		vertices.clear();
		for (int k = 0; k < vL.size(); ++k)
		{
			vertices.push_back(vL[k]);
		}

		vertices.push_back(i);

		for (int k = 0; k < vR.size(); ++k)
		{
			vertices.push_back(vR[k]);
		}

	}


	int * feedback_arc_set_pivot(IGraph &graph){
		srand(time(NULL));

		int n = graph.numVertices();
		
		int *v = new int[n];
		vector<int> vertices;
		for (IGraph::vertex_iterator v = graph.begin(); v != graph.end(); ++v)
		{
			vertices.push_back(v->id);
		}
		FAS_pivot(graph, vertices);
		// /quickSort(graph, vertices, 0, n-1);
		for (int i = 0; i < vertices.size(); ++i)
		{
			v[i] = vertices[i];
			cout << (char)((int)'A' + vertices[i]) << " ";
		}
		cout << endl;

		return v;
	}


}