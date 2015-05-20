
#include <iomanip>
#include <iostream>

#include <string.h>
#include <time.h>
#include <cmath>

#include <algorithm>

#include <queue>
#include <vector>
#include <map>

#include "IGraph.h"
#include "Graph_Adj_Matrix2.hpp"
#include "Social_Graph_Adj_Matrix.hpp"

using namespace std;


void print_array(int a[], int n){
	for (int l = 0; l < n; ++l)
	{
		cout << a[l] << ",";
	}
	cout << endl;
	
}

int main(int argc, char const *argv[])
{
	Social_Graph_Adj_Matrix G(8);

	G.addEdge(0, 3, 1);
	G.addEdge(0, 5, 1);
	//G.addEdge(0, 6, 1);
	//G.addEdge(0, 7, 1);

	G.addEdge(1, 2, 1);
	G.addEdge(1, 3, 1);
	G.addEdge(1, 4, 1);
	//G.addEdge(1, 6, 1);
	//G.addEdge(1, 7, 1);

	G.addEdge(2, 4, 1);
	G.addEdge(2, 5, 1);
	G.addEdge(2, 6, 1);
	
	G.addEdge(3, 7, 1);
	
	G.addEdge(4, 5, 1);
	//G.addEdge(4, 6, 1);

	G.addEdge(6, 7, 1);

	G.printAsMatrix();

	cout << G.eGroupA << " " << G.eGroupAB << endl;

	G.putVertexInA(*G.getVertex(0));

	cout << G.eGroupA << " " << G.eGroupAB << endl;

	G.putVertexInA(*G.getVertex(1));

	cout << G.eGroupA << " " << G.eGroupAB << endl;

	print_array(G.inGroup, 8);
	print_array(G.outGroup, 8);

	G.partionate();

	Social_Graph_Adj_Matrix G1(7);

	G1.addEdge(0, 1, 1);
	G1.addEdge(0, 2, 1);
	//G.addEdge(0, 6, 1);
	//G.addEdge(0, 7, 1);

	//G1.addEdge(1, 2, 1);
	G1.addEdge(1, 3, 1);
	//G.addEdge(1, 6, 1);
	//G.addEdge(1, 7, 1);

	G1.addEdge(2, 3, 1);
	G1.addEdge(2, 4, 1);
	//G.addEdge(2, 6, 1);
	
	G1.addEdge(3, 4, 1);
	G1.addEdge(3, 5, 1);
	G1.addEdge(3, 6, 1);
	
	G1.addEdge(4, 5, 1);
	G1.addEdge(4, 6, 1);

	G1.addEdge(5, 6, 1);

	G1.putVertexInA(*G1.getVertex(1));
	G1.putVertexInA(*G1.getVertex(0));

	G1.printAsMatrix();

	G1.partionate();
}