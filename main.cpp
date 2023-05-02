#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <random>
#include <string>
#include <chrono>
#include <unordered_map>
#include "file.h"
#include "utils.h"
//#include <omp.h>

using namespace std;
using namespace std::chrono;

vector<vector<int>> edgeListToQUBO(vector<Edge>& edges, int numNodes);
vector<vector<int>> random_2d_spin(int , int );
vector<vector<int>> local_field_init(int , int );
float str_c(float, int, int);
void Unit_Test_str(int MC_Step, int M, float mystery);

vector<int> spin_to_binary(const vector<int>& spin_vector);
int choice_spin(vector<vector<int>>, vector<vector<int>>);
int ising_energy(vector<vector<int>>, vector<int>);

random_device rd;
mt19937 gen1(rd());
uniform_real_distribution<double> uniform_dist(0.0, 1.0);


int main(){
    unordered_map<string, string> my_graph;
    my_graph["G1"] = "11624";
    my_graph["G2"] = "11620"; 
    my_graph["G3"] = "11622";  
    my_graph["G4"] = "11646"; 
    my_graph["G5"] = "11631";
    my_graph["G14"] = "3064";
    my_graph["G15"] = "3050";
    my_graph["G16"] = "3052";
    my_graph["G17"] = "3047";
    my_graph["G36"] = "7680";
    my_graph["G47"] = "6657";
    my_graph["G51"] = "3848";
    my_graph["G55"] = "10289";
    my_graph["G58"] = "19293";
    my_graph["G60"] = "14188";
    my_graph["G63"] = "27045";
    my_graph["G70"] = "9522";
    string filename = "G1";
    
    // Coupling Matrix ← Couplings Data
    ReadGraph J_Data(filename);
    vector<Edge> edges = J_Data.getEdges();
    int numNodes = J_Data.getNumNodes();

    // To Qubo Form
    vector<vector<int>> graph = edgeListToQUBO(edges, numNodes);
    
    int N = graph.size();
    int trotter_M = 32; 
    auto MC_Step = 1024*7;

    /*** This function will help you adjust the Transverse field strength parameter to optimize your Ising Energy.
    // Refer to the Paper 4 Performance Evaluation: G0 is initially set to 8
    Unit_Test_str(MC_Step, trotter_M, 8);
    ***/

    //Randomly initialize σ ∈ {−1, 1}
    vector<vector<int>> local_field = local_field_init(N, trotter_M); 

    //localfield <- 0
    vector<vector<int>> r_spin = random_2d_spin(N, trotter_M); 
   
    // Construct local-field energy
    for(int i = 0; i < N; i++){
        for(int M = 0; M < trotter_M; M++){
            for( int j = 0; j < N; j++){
                local_field[i][M] += graph[i][j]*r_spin[j][M];
            }
        }
    }
    
    float rate, delta_H, T;

    cout<<"Gset: "<< filename <<endl;
    cout<<"Graph Size: "<< N << endl;
    cout<<"Number of Trotters: "<< trotter_M<<endl;
    cout<<"MC Steps: "<< MC_Step <<endl;
    auto start = high_resolution_clock::now();
    for(auto t_mc= 1; t_mc <= MC_Step; t_mc++){
        for(int i= 0; i < N; i++){
            for(auto M = 0; M < trotter_M; M++){            

                if(M!=0 && M!=trotter_M-1 ){
                    delta_H = r_spin[i][M]*(local_field[i][M] - str_c(M, t_mc, MC_Step)*(r_spin[i][M+1] + r_spin[i][M-1]));
                }
                else if(M==0){
                    delta_H = r_spin[i][M]*(local_field[i][M] - str_c(M, t_mc, MC_Step)*(r_spin[i][M+1]));
                }
                else if(M==trotter_M-1){
                    delta_H = r_spin[i][M]*(local_field[i][M] - str_c(M, t_mc, MC_Step)*(r_spin[i][M-1]));
                }

                rate = uniform_dist(gen1);
                T = 1./(t_mc*(1. - 1./8.)/ static_cast<float>(MC_Step));
                if(exp(-1*delta_H / T ) > rate ){
                    r_spin[i][M] = -1 * r_spin[i][M];
                    for(int j = 0; j < N; j++ ){
                        local_field[j][M] += 2 * graph[i][j]*r_spin[i][M]; 
                    }
                }
            }
        }
    }

    int maxcut = choice_spin(r_spin, graph);
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();
    cout << "Time taken by program: " << duration << " milliseconds" << endl;    
    cout<<"Max Cut Value: "<< maxcut <<"\n";
    cout<<"The Best Value: "<< my_graph[filename] <<"\n";
    cout<<"Accuracy: "<< maxcut /  static_cast<float>(stoi(my_graph[filename])) <<endl;
};  
