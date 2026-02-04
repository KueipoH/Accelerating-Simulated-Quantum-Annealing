#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <random>
#include <string>
#include <chrono>
#include <unordered_map>
#include <cmath>
#include "file.h"
#include "utils.h"
#include "gset.h"

using namespace std;
using namespace std::chrono;
extern unordered_map<string, string> my_graph;

vector<vector<int>> edgeListToQUBO(vector<Edge>& edges, int numNodes);
vector<int> spin_to_binary(const vector<int>& spin_vector);
int choice_spin(const vector<vector<int>>&, const vector<vector<int>>&);
int ising_energy(const vector<vector<int>>&, const vector<int>&);

random_device rd;
mt19937 gen1(rd());
uniform_real_distribution<double> uniform_dist(0.0, 1.0);

int main(){

    string filename = "G1";

    // Coupling Matrix <- Couplings Data
    ReadGraph J_Data(filename);
    vector<Edge> edges = J_Data.getEdges();
    int numNodes = J_Data.getNumNodes();

    // To Qubo Form (kept for final energy calculation)
    vector<vector<int>> graph = edgeListToQUBO(edges, numNodes);

    int N = graph.size();
    int trotter_M = 32;
    int MC_Step = 1024*7;

    // [Opt 1] Build adjacency list from QUBO (sparse representation)
    // Instead of iterating all N nodes for local_field update, only visit neighbors
    vector<vector<pair<int,int>>> adj(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (graph[i][j] != 0) {
                adj[i].push_back({j, graph[i][j]});
            }
        }
    }

    // [Opt 2] Flat 1D arrays for better cache locality
    vector<int> r_spin(N * trotter_M);
    vector<int> local_field(N * trotter_M, 0);

    // Random spin init {-1, +1}
    {
        random_device rd1;
        mt19937 rng(rd1());
        uniform_int_distribution<int> dis(0, 1);
        for (int i = 0; i < N * trotter_M; i++) {
            r_spin[i] = dis(rng) * 2 - 1;
        }
    }

    // Construct initial local field using adjacency list
    for (int i = 0; i < N; i++) {
        for (int m = 0; m < trotter_M; m++) {
            int sum = 0;
            for (const auto& p : adj[i]) {
                sum += p.second * r_spin[p.first * trotter_M + m];
            }
            local_field[i * trotter_M + m] = sum;
        }
    }

    float delta_H;

    cout<<"Gset: "<< filename <<endl;
    cout<<"Graph Size: "<< N << endl;
    cout<<"Number of Trotters: "<< trotter_M<<endl;
    cout<<"MC Steps: "<< MC_Step <<endl;
    auto start = high_resolution_clock::now();

    // [Opt 3] Precompute buffer for str_c values (avoid log/cosh in hot loop)
    vector<float> str_c_vals(trotter_M);

    for (int t_mc = 1; t_mc <= MC_Step; t_mc++) {
        // [Opt 4] Precompute T once per MC step (was recomputed N*M times)
        float T = 1.0f / (t_mc * (1.0f - 1.0f/8.0f) / static_cast<float>(MC_Step));

        // Precompute str_c for each Trotter index at this MC step
        float mystery = 32.0f;
        float t_ratio = static_cast<float>(t_mc / MC_Step); // integer division preserved
        for (int m = 0; m < trotter_M; m++) {
            str_c_vals[m] = T / 2.0f * logf(coshf(mystery * (1.0f - t_ratio) / (m * T)));
        }

        float neg_inv_T = -1.0f / T;

        for (int i = 0; i < N; i++) {
            int base_i = i * trotter_M;
            for (int m = 0; m < trotter_M; m++) {
                int idx = base_i + m;

                // [Opt 5] Simplified Trotter boundary (no 3-way branch)
                int trotter_coupling = 0;
                if (m > 0) trotter_coupling += r_spin[idx - 1];
                if (m < trotter_M - 1) trotter_coupling += r_spin[idx + 1];

                delta_H = r_spin[idx] * (local_field[idx] - str_c_vals[m] * trotter_coupling);

                float rate = uniform_dist(gen1);
                if (expf(delta_H * neg_inv_T) > rate) {
                    r_spin[idx] = -r_spin[idx];
                    // [Opt 1] Sparse local field update via adjacency list
                    int flip_val = 2 * r_spin[idx];
                    for (const auto& p : adj[i]) {
                        local_field[p.first * trotter_M + m] += p.second * flip_val;
                    }
                }
            }
        }
    }

    // Convert flat r_spin back to 2D for choice_spin (one-time cost)
    vector<vector<int>> r_spin_2d(N, vector<int>(trotter_M));
    for (int i = 0; i < N; i++)
        for (int m = 0; m < trotter_M; m++)
            r_spin_2d[i][m] = r_spin[i * trotter_M + m];

    int maxcut = choice_spin(r_spin_2d, graph);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();
    cout << "Time taken by program: " << duration << " milliseconds" << endl;
    cout<<"Max Cut Value: "<< maxcut <<"\n";
    cout<<"The Best Value: "<< my_graph[filename] <<"\n";
    cout<<"Accuracy: "<< maxcut /  static_cast<float>(stoi(my_graph[filename])) <<endl;
};
