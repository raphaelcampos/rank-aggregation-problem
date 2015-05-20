#include <string.h>
#include <queue>

#include "IGraph.h"

enum {BLACK, WHITE, GRAY};

/**
 * O(V + E)
 * @param G  Graph
 * @param s source vertex
 */
void BFS(IGraph &G, IGraph::vertex &s){
	
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

void DFS_visit(IGraph &G, IGraph::vertex &v, int color[], int d[], int f[], IGraph::vertex **pi, int &dtime){

	dtime++;
	color[v.id] = GRAY;
	d[v.id] = dtime;
	
	cout << v.id << endl;	

	for (IGraph::vertex::iterator e = v.begin(); e != v.end(); e++){
		IGraph::vertex *u = e->second;
		if(color[u->id] == WHITE){	
			pi[u->id] = &v;
			DFS_visit(G, *u, color, d, f, pi, dtime);
		}
	}

	dtime++;
	color[v.id] = BLACK;
	f[v.id] = dtime;
}

/**
 * complexity : O(V + E)
 * @param G Graph
 * @param s source vertex
 */
void DFS(IGraph &G, IGraph::vertex &s){

	int n = G.numVertices();
	int *color = new int[n];
	int *d = new int[n];
	int *f = new int[n];
	IGraph::vertex **pi = new IGraph::vertex*[n];
	
	for (IGraph::vertex_iterator v = G.begin(); v != G.end(); ++v)
	{
		color[v->id] = WHITE;
		pi[v->id] = NULL;
		d[v->id] = 0;
		f[v->id] = 0;
	}

	int dtime = 0;

	DFS_visit(G, s, color, d, f, pi, dtime);
	for (IGraph::vertex_iterator v = G.begin(); v != G.end(); ++v)
	{
		if(color[v->id] == WHITE){
			DFS_visit(G, *v, color, d, f, pi, dtime);
		}
	}

	cout << "S : " << s.id << endl;
	for (int i = 0; i < n; ++i)
	{
		if(pi[i] != NULL)
			cout << pi[i]->id << "(" << d[i] << "/" << f[i] <<  ") ";
		else
			cout << "nil" << "(" << d[i] << "/" << f[i] <<  ") ";
	}
	cout << endl;

	delete [] color;
	delete [] d;
	delete [] f;
	delete [] pi;
}