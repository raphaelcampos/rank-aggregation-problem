#include "IGraph.h"

template<class T>
IGraph<T> complete2Tournament(IGraph<T> &comp){
	IGraph<T> tour(comp.numVertices());

	if(comp.isComplete()){
		for(IGraph<int>::vertex_iterator v = comp.begin(); v != comp.end(); v++){

			for (IGraph<int>::vertex::iterator u = v->begin(); u != v->end(); ++u)
			{
				
			}

		}

	}

	return tour;
}

template<class T1, class T2>
T1 * hamiltonPathForTournament(IGraph<T2> tour){

}