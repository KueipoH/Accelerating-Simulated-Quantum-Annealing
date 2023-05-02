// file.cpp
#include "file.h"

ReadGraph::ReadGraph(const string& filename) {
    ifstream inFile(filename);
    if (inFile.is_open()) {
        string line;
        getline(inFile, line);
        while (getline(inFile, line)) {
            stringstream ss(line);
            int u, v;
            ss >> u >> v;
            edges.push_back({u, v});
            numNodes = max(numNodes, max(u, v));
        }
        inFile.close();
    }
}

vector<Edge> ReadGraph::getEdges() const {
    return edges;
}

int ReadGraph::getNumNodes() const {
    return numNodes;
}