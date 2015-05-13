#include "IGraph.h"

void complete2Tournament(IGraph &comp, IGraph &tour){

	if(comp.isComplete()){
		for (IGraph::vertex_iterator v = comp.begin(); v != comp.end(); v++){
			for (IGraph::vertex::iterator u = v->begin(); u != v->end(); u++){
		
				double w1 = comp.getEdge(v->id, u->second->id);
				double w2 = comp.getEdge(u->second->id, v->id);
				
				if(w1 > w2)
					tour.addEdge(v->id, u->second->id, 1);
				else
					tour.addEdge(u->second->id, v->id, 1);
			}
		}

	}
}

int * hamiltonPathForTournament(IGraph &tour){
	int * path = new int[tour.numVertices()];

	int n = tour.numVertices(); 

	queue< IGraph::vertex*> Q;
	int max = 0;

	IGraph::vertex * v = NULL;
	for (IGraph::vertex_iterator it = tour.begin(); it != tour.end(); it++)
	{
		if(max < it->outdegree){
			max = it->outdegree;
			
			v = &(*it);
		}
	}

	Q.push(v);
	
	int it = 0;
	while(!Q.empty() && it < n){

		v = Q.front();
		Q.pop();
		
		if(v->begin() != v->end())
			if(IGraph::vertex * aux = dynamic_cast<IGraph::vertex*>(v->begin()->second)){
				
				max = aux->outdegree - (int)tour.edgeExists(aux->id, v->id);
				for (IGraph::vertex::iterator u = v->begin(); u != v->end(); ++u)
				{	
					if(max < u->second->outdegree - (int)tour.edgeExists(u->second->id, v->id)){
						max = u->second->outdegree - (int)tour.edgeExists(u->second->id, v->id);
						aux = u->second;
					}
				}

				Q.push(aux);
			}
		
		path[it] = v->id;

		tour.removeVertex(v->id);
		it++;
	}

	return path;
}