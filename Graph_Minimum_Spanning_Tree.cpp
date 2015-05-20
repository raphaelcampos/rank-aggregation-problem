#include <string.h>
#include <queue>
#include <algorithm>
#include <functional>

#include "IGraph.h"

/**
 * O(V + Vlog(V) + Elog(V))
 * @param G [description]
 * @param s [description]
 */
void MST_prim(IGraph &G, IGraph::vertex &s){
	int n = G.numVertices();
	bool *isQ = new bool[n];
	IGraph::vertex ** pi = new IGraph::vertex *[n];

	memset(isQ, 1, sizeof(bool)*n);

	typedef pair<double, int> vi;
	
 	vector<vi> vertices;

	for (int i = 0; i < n; ++i)
	{
		vertices.push_back(make_pair(1000000, i));
		pi[i] = NULL;
	}

	vertices[s.id].first = 0;
	priority_queue< vi, vector<vi>, std::greater<vi> > Q(vertices.begin(), vertices.begin() + n);

	int it = 0;
	while(it < n){

		vi e = Q.top(); // O(1)
		Q.pop(); // O(lg(v))

		isQ[e.second] = false;
		IGraph::vertex *u = G.getVertex(e.second);
		for (IGraph::vertex::iterator v = u->begin(); v != u->end(); ++v)
		{
			double w = G.getEdge(u->id, v->second->id);
			if(isQ[v->second->id] &&  w < vertices[v->second->id].first){// O(1)
				pi[v->second->id] = u;
				Q.push(make_pair(w, v->second->id)); // O(log(v))
				vertices[v->second->id].first = w;
			}
		}

		it++;

	}

	cout << "PRIMs : ";
	for (int i = 0; i < n; ++i)
	{
		if(pi[i] != NULL)
			cout << pi[i]->id << " ";
		else
			cout << "nil ";
	}
	cout << endl;

	delete [] isQ;

}