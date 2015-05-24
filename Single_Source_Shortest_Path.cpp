#include "IGraph.h"
#include "PriorityQueue.cpp"

void dijkstra(IGraph &G, IGraph::vertex &s){
	int n = G.numVertices();
	bool *isQ = new bool[n];
	IGraph::vertex ** pi = new IGraph::vertex *[n];

	typedef pair<double, int> vi;
	
 	vector<vi> d;

 	IndexedPriorityQueue<double, less<double> > Q(n);

	for (int i = 0; i < n; ++i)
	{
		d.push_back(make_pair(1000000, i));
		Q.insert(i, 1000000);
		pi[i] = NULL;
		isQ[i] = true;
	}

	d[s.id].first = 0;
	Q.changeKey(s.id, 0);

	//priority_queue< vi, vector<vi>, std::greater<vi> > Q(vertices.begin(), vertices.begin() + n);

	int it = 0;
	while(!Q.isEmpty()){

		int id = Q.minIndex(); // O(1)
		Q.deleteMin(); // O(lg(v))

		isQ[id] = false;
		IGraph::vertex *u = G.getVertex(id);
		for (IGraph::vertex::iterator v = u->begin(); v != u->end(); ++v)
		{
			double w = G.getEdge(id, v->second->id);
			if(w + d[id].first < d[v->second->id].first){// O(1)
				pi[v->second->id] = u;
				Q.changeKey(v->second->id, w + d[id].first); // O(log(v))
				d[v->second->id].first = w + d[id].first;
			}
		}

	}

	cout << "DIJKSTRA : ";
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