#!/bin/bash
rm -r tp-paa-submit
mkdir tp-paa-submit

echo "copying Permutation..."
cp -r Permutation tp-paa-submit/.
echo "copying graph lib..."
cp Graph_* tp-paa-submit/.
cp IGraph.h tp-paa-submit/.
cp PriorityQueue.cpp tp-paa-submit/.
echo "copying klib..."
cp kutils.hpp tp-paa-submit/.
cp kalgorithms.hpp tp-paa-submit/.
cp kstructure.h tp-paa-submit/.
cp kendall-tau-distance.hpp tp-paa-submit/.
echo "copying utils..."
cp utils.h tp-paa-submit/.
echo "copying mains file..."
cp rank-aggregation.cpp tp-paa-submit/.
cp time-measure-aggr.cpp tp-paa-submit/.
echo "copying bash file..."
cp execute_rank_aggre_30times.sh tp-paa-submit/.
cp execute-time-measuring.sh tp-paa-submit/.
