#!/bin/bash
rm -r tp-grafos-submit
mkdir tp-grafos-submit

echo "copying dlib..."
cp -r dlib tp-grafos-submit/.
echo "copying graph lib..."
cp Graph_* tp-grafos-submit/.
cp IGraph.h tp-grafos-submit/.
cp PriorityQueue.cpp tp-grafos-submit/.
echo "copying social graph"
cp Social_Graph_Adj_Matrix.hpp tp-grafos-submit/.
echo "copying main file..."
cp tp-graph.cpp  tp-grafos-submit/.
echo "copying bash files..."
cp compile.sh tp-grafos-submit/.
cp execute.sh tp-grafos-submit/.