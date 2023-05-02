// file.h
#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

struct Edge {
    int u, v;
};

class ReadGraph{
public:
    ReadGraph(const string& filename);
    vector<Edge> getEdges() const;
    int getNumNodes() const;

private:
    vector<Edge> edges;
    int numNodes = 0;
};
