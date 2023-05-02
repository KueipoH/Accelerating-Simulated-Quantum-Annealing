#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_map>
//#include <omp.h>

random_device rd1;
mt19937 gen(rd1());
uniform_int_distribution<int> dis(0, 1);

vector<int> spin_to_binary(const vector<int>& spin_vector) {
    vector<int> binary_vector(spin_vector.size());
    for (size_t i = 0; i < spin_vector.size(); ++i) {
        binary_vector[i] = (spin_vector[i] + 1) / 2;
    }
    return binary_vector;
}

vector<vector<int>> edgeListToQUBO(vector<Edge>& edges, int numNodes) {
    vector<vector<int>> Qubo(numNodes, vector<int>(numNodes));
    
    for (auto& edge : edges) {
        Qubo[edge.u-1][edge.v-1] -= 1;
        Qubo[edge.v-1][edge.u-1] -= 1;
        Qubo[edge.u-1][edge.u-1] += 1;
        Qubo[edge.v-1][edge.v-1] += 1;
    }
    return Qubo;
}

float str_c(float M, int t_mc, int MC_Step){
    float mystery = 32.;
    float T = 1./(t_mc*(1. - 1./8.)/ static_cast<float>(MC_Step));
    return T/2.*log(cosh(mystery*(1.-t_mc/MC_Step)/(M* T)));
} 

vector<vector<int>> random_2d_spin(int N, int trotter_M){
    vector<vector<int>> spin_init(N, vector<int>(vector<int>(trotter_M)));
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < trotter_M; k++) {
            spin_init[i][k] = dis(gen) * 2 - 1;
        }
    }
    return spin_init;
}

vector<vector<int>> local_field_init(int N, int trotter_M) {
    return vector<vector<int>>(N, vector<int>(trotter_M, 0));
}

int ising_energy(vector<vector<int>>graph, vector<int>spin){
    int calc_energy = 0;
    vector<int> result(graph[0].size());

    int num_threads = 8;

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for(auto it = 0; it < spin.size(); it++){
        int temp = 0;
        for(auto jt = 0; jt < graph[it].size(); jt++){
            temp += spin[jt]*graph[it][jt];
        }
        result[it] = temp;
        calc_energy += result[it]*spin[it];
    }
    return -1*calc_energy;
};

int choice_spin(vector<vector<int>>r_spin, vector<vector<int>>graph){
    int N = r_spin.size();
    int trotter_M = r_spin[0].size();
    unordered_map<int, int> e_map;
    for(int M = 0; M < trotter_M; M++){
        vector<int> local_spin(N, 0);
        for(int i = 0; i < N; i++){
            local_spin[i] = r_spin[i][M];
        }
        e_map[M] = -1*ising_energy(graph, spin_to_binary(local_spin) );
    }
    int max_key = -1;
    int max_value = -1;

    for (const auto &pair : e_map) {
        if (pair.second > max_value) {
            max_key = pair.first;
            max_value = pair.second;
        }
    }
    return max_value;
} 

void Unit_Test_str(int MC_Step, int trotter_M, float mystery){
    float T;
    for(auto t_mc= 1; t_mc < MC_Step ; t_mc++){
        for(auto M = 1; M < trotter_M; M++){
            T = 1./(t_mc*(1. - 1./8.)/ static_cast<float>(MC_Step));
            cout<<"MC: "<<t_mc<<" M: "<<M<<" T: "<< T <<" str: "<<T/2.*log(cosh(mystery*(1.-t_mc/MC_Step)/(M* T)))<<endl;
        }
    } 
}