#include "IGraph.h"
#include "PriorityQueue.cpp"

#include <cmath>

namespace graph{

	void dijkstra(IGraph &G, IGraph::vertex &s, double inf = 1e100){
		int n = G.numVertices();
		bool *isQ = new bool[n];
		IGraph::vertex ** pi = new IGraph::vertex *[n];

		typedef pair<double, int> vi;
		
	 	vector<vi> d;

	 	IndexedPriorityQueue<double, less<double> > Q(n);

		for (int i = 0; i < n; ++i)
		{
			d.push_back(make_pair(inf, i));
			Q.insert(i, inf);
			pi[i] = NULL;
		}

		d[s.id].first = 0;
		Q.changeKey(s.id, 0);

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
	}

	bool bellman_ford(IGraph &G, IGraph::vertex &s, double inf = 1e100){
		int n = G.numVertices();
		bool *isQ = new bool[n];
		IGraph::vertex ** pi = new IGraph::vertex *[n];
		vector<double> d;

		for (int i = 0; i < n; ++i)
		{
			d.push_back(inf);
			pi[i] = NULL;
		}

		d[s.id] = 0;

		// performs relaxation |V| - 1 times
		// for all edges in G 
		for (int i = 0; i < n-1; ++i)
		{
			// complexity O(E)
			for(IGraph::vertex_iterator u = G.begin(); u != G.end(); u++){
				for (IGraph::vertex::iterator e = u->begin(); e != u->end(); e++)
				{
					double w = e->first;
					IGraph::vertex *v = e->second;
					//relaxation
					if(d[v->id] > d[u->id] + w){
						d[v->id] = d[u->id] + w;
						pi[v->id] = &(*u);
					}
				}
			}
		}

		// complexity O(E)
		for(IGraph::vertex_iterator u = G.begin(); u != G.end(); u++){
			for (IGraph::vertex::iterator e = u->begin(); e != u->end(); e++)
			{
				double w = e->first;
				IGraph::vertex *v = e->second;
				//relaxation
				if(d[v->id] > d[u->id] + w){
					return false;
				}
			}
		}

		cout << "BELLMAN : ";
		for (int i = 0; i < n; ++i)
		{
			if(pi[i] != NULL)
				cout << pi[i]->id << " ";
			else
				cout << "nil ";
		}
		cout << endl;

		return true;

	}

}