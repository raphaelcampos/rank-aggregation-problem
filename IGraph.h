#ifndef __I_GRAPH__
#define __I_GRAPH__

template<class T>
class IGraph { 

	public:
        virtual ~IGraph() {}
        virtual void addVertex(int v) = 0;
        virtual void removeVertex(int v) = 0;
        virtual void addEdge(int v, int u, T w) = 0;
        virtual void removeVertex(int v, int u) = 0;
        virtual void inDegree(int v) = 0;
        virtual void outDegree(int v) = 0;
        virtual void getVertex(int v) = 0;
        virtual T getEdge(int v, int u) = 0;
        virtual void updateEdgeWeight(int v, int u, T w) = 0;
        virtual bool edgeExists(int v, int u) = 0;
        virtual bool vertexExists(int v) = 0;
        virtual void nextAdjVertex(int v) = 0;
        virtual bool isComplete() = 0;
};

#endif 