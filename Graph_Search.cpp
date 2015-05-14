#include <string.h>
#include <queue>

#include "IGraph.h"

void BFS(IGraph &G, IGraph::vertex &s){
	enum {BLACK, WHITE, GRAY};
	int n = G.numVertices();
	int *color = new int[n];
	int *d = new int[n];
	IGraph::vertex **pi = new IGraph::vertex*[n];
	
	for (IGraph::vertex_iterator v = G.begin(); v != G.end(); ++v)
	{
		color[v->id] = WHITE;
		pi[v->id] = NULL;
		d[v->id] = 1000000;
	}

	color[s.id] = GRAY;
	d[s.id] = 0;

	queue<IGraph::vertex * > Q;
	Q.push(&s); 
	while(!Q.empty()){

		IGraph::vertex * u = Q.front();
		for (IGraph::vertex::iterator e = u->begin(); e != u->end(); e++)
		{
			IGraph::vertex * v = e->second;	
			if(color[v->id] == WHITE){
				color[v->id] = GRAY;
				d[v->id] = d[u->id] + 1;
				pi[v->id] = u;
				Q.push(&(*v)); 
			}
		}
		color[u->id] = BLACK;
		Q.pop();
	}

	cout << "S : " << s.id << endl;
	for (int i = 0; i < n; ++i)
	{
		if(pi[i] != NULL)
			cout << pi[i]->id << " ";
		else
			cout << "nil ";
	}
	cout << endl;

	delete [] color;
	delete [] d;
	delete [] pi;
}