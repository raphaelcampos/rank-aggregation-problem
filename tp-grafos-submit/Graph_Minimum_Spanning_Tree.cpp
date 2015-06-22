#include <string.h>
#include <queue>
#include <algorithm>
#include <functional>

#include "IGraph.h"
#include "PriorityQueue.cpp"

/**
 * O(V + Vlog(V) + Elog(V))
 * @param G [description]
 * @param s [description]
 */
void MST_prim(IGraph &G, IGraph::vertex &s, double inf = 1e100){
	int n = G.numVertices();
	bool *isQ = new bool[n];
	IGraph::vertex ** pi = new IGraph::vertex *[n];

	typedef pair<double, int> vi;
	
 	vector<vi> vertices;

 	IndexedPriorityQueue<double, less<double> > Q(n);

	for (int i = 0; i < n; ++i)
	{
		vertices.push_back(make_pair(inf, i));
		Q.insert(i, inf);
		pi[i] = NULL;
		isQ[i] = true;
	}

	vertices[s.id].first = 0;
	Q.changeKey(s.id, 0);

	int it = 0;
	while(!Q.isEmpty()){

		int id = Q.minIndex(); // O(1)
		Q.deleteMin(); // O(lg(v))

		cout << id << endl;

		isQ[id] = false;
		IGraph::vertex *u = G.getVertex(id);
		for (IGraph::vertex::iterator v = u->begin(); v != u->end(); ++v)
		{
			double w = G.getEdge(id, v->second->id);
			if(isQ[v->second->id] &&  w < vertices[v->second->id].first){// O(1)
				pi[v->second->id] = u;
				Q.changeKey(v->second->id, w); // O(log(v))
				vertices[v->second->id].first = w;
			}
		}

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