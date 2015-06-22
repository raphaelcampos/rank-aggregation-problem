#ifndef __I_G_HAMILTON_PATH__
#define __I_G_HAMILTON_PATH__

#include <cmath>
#include "IGraph.h"

void complete2Tournament(IGraph &comp, IGraph &tour, bool unweighted = true){

	//if(comp.isComplete()){
	//	cout << "complete2Tournament" << endl;
		for (IGraph::vertex_iterator v = comp.begin(); v != comp.end(); v++){
			for (IGraph::vertex::iterator u = v->begin(); u != v->end(); u++){
				if(!(tour.edgeExists(v->id,u->second->id) || tour.edgeExists(u->second->id, v->id))){
					double w1 = comp.getEdge(v->id, u->second->id);
					double w2 = comp.getEdge(u->second->id, v->id);
					
					double wVal = (unweighted)? 1 : abs(w1 - w2); 

					if(w1 > w2)
						tour.addEdge(v->id, u->second->id, (unweighted)? 1 : w1);
					else if(w1 < w2)
						tour.addEdge(u->second->id, v->id, (unweighted)? 1 : w2);
					else{
						tour.addEdge(u->second->id, v->id, (unweighted)? 1 : w2);
						tour.addEdge(v->id, u->second->id, (unweighted)? 1 : w1);
					}
				}
			}
		}
	//}
}

int * hamiltonPathForTournament(IGraph &tour){
	enum {BLACK, WHITE, GRAY};
	
	int * path = new int[tour.numVertices()];
	int * color = new int[tour.numVertices()];
	int * indegree = new int[tour.numVertices()];

	queue< IGraph::vertex* > Q;
	int n = tour.numVertices();

	int min = 100000;
	IGraph::vertex * u = NULL;
	for (IGraph::vertex_iterator it = tour.begin(); it != tour.end(); it++)
	{
		color[it->id] = WHITE;
		indegree[it->id] = it->indegree;
		if(min > it->indegree){
			min = it->indegree;
			u = &(*it);
		}
	}

	color[u->id] = GRAY;
	Q.push(u);
	int i = 0;
	while(!Q.empty()){
		u = Q.front();
		Q.pop();

		int min = 100000;
		IGraph::vertex * s = NULL;
		for (IGraph::vertex::iterator e = u->begin(); e != u->end(); ++e)
		{
			IGraph::vertex *v = e->second;
			if(color[v->id] == WHITE){
				indegree[v->id]--;
				
				if(min > indegree[v->id]){
					min = indegree[v->id];	
					s = &(*v);
				}
			}
		}
		if(s != NULL){
			color[s->id] = GRAY;
			Q.push(s);
		}

		color[u->id] = BLACK;
		path[i] = u->id;
		i++;
	}

	return path;
}

/*int * hamiltonPathForTournament(IGraph &tour){
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
}*/

int tour_ham_merge_and_count(IGraph &tour, int *A, int *buffer, int ini, int fim, int meio){

	memset(buffer, 0, sizeof(int)*(fim-ini+1));

	int count = 0;

	int i = ini;
	int j = meio + 1;
	
	int k = 0;
	while(i < meio + 1 && j <= fim){
		if(tour.edgeExists(A[i], A[j])){
			buffer[k] = A[i];
			i++;
		}else{
			buffer[k] = A[j];
			count += meio - i + 1;
			j++;
		}
		k++;
	}

	if(i < meio + 1){
		memcpy(buffer + k, A + i, sizeof(int)*(meio-i+1));
	}else{
		memcpy(buffer + k, A + j, sizeof(int)*(fim-j+1));
	}
	
	memcpy(A + ini, buffer, sizeof(int)*(fim-ini+1));

	return count;
}

int tour_ham_sort_and_count(IGraph &tour, int A[], int buffer[], int ini, int fim){

	if(fim + 1 - ini <= 1){
		return 0;
	}else{
		int meio = (ini + fim)/2;
		return tour_ham_sort_and_count(tour, A, buffer, ini, meio) +
		 tour_ham_sort_and_count(tour, A, buffer, meio + 1, fim) + 
		 tour_ham_merge_and_count(tour, A, buffer, ini, fim, meio);
		 ;
	}
}

void DVhamiltonPathForTournament(IGraph &tour, int * path){
	int * buffer = new int[tour.numVertices()];

	tour_ham_sort_and_count(tour, path, buffer, 0, tour.numVertices() - 1);
}

int * DVhamiltonPathForTournament(IGraph &tour){
	int * path = new int[tour.numVertices()];
	int * buffer = new int[tour.numVertices()];

	for (IGraph::vertex_iterator it = tour.begin(); it != tour.end(); it++)
	{
		path[it->id] = it->id;
	}

	tour_ham_sort_and_count(tour, path, buffer, 0, tour.numVertices() - 1);

	return path;
}

#endif