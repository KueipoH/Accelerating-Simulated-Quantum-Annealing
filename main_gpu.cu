/*
 * CUDA implementation of Simulated Quantum Annealing using cuBLAS GEMM + Tensor Cores
 * Based on Algorithm 2 - uses FP16 tensor cores for local_field matrix operations
 *
 * Compile: nvcc -O3 -arch=sm_70 main_gpu.cu file.cpp -lcublas -o sqa_gpu
 * (Use sm_80 for Ampere, sm_90 for Hopper)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <string>
#include <chrono>
#include <unordered_map>
#include <cmath>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuda_fp16.h>

#include "file.h"
#include "gset.h"

using namespace std;
using namespace std::chrono;

extern unordered_map<string, string> my_graph;

// Error checking macros
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " << cudaGetErrorString(err) << endl; \
        exit(1); \
    } \
} while(0)

#define CUBLAS_CHECK(call) do { \
    cublasStatus_t status = call; \
    if (status != CUBLAS_STATUS_SUCCESS) { \
        cerr << "cuBLAS error at " << __FILE__ << ":" << __LINE__ << " - " << status << endl; \
        exit(1); \
    } \
} while(0)

// Build QUBO matrix from edges
vector<vector<int>> edgeListToQUBO(vector<Edge>& edges, int numNodes) {
    vector<vector<int>> Qubo(numNodes, vector<int>(numNodes, 0));
    for (auto& edge : edges) {
        Qubo[edge.u-1][edge.v-1] -= static_cast<int>(edge.w);
        Qubo[edge.v-1][edge.u-1] -= static_cast<int>(edge.w);
        Qubo[edge.u-1][edge.u-1] += static_cast<int>(edge.w);
        Qubo[edge.v-1][edge.v-1] += static_cast<int>(edge.w);
    }
    return Qubo;
}

// CUDA kernel: Convert int to half precision
__global__ void convert_int_to_half(const int* src, half* dst, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = __float2half(static_cast<float>(src[idx]));
    }
}

// CUDA kernel: Convert half to float
__global__ void convert_half_to_float(const half* src, float* dst, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = __half2float(src[idx]);
    }
}

// CUDA kernel: Metropolis update for a block of spins
// Each thread handles one (i, m) pair within the block
__global__ void metropolis_update_kernel(
    int* d_spin,              // [N x M] spin configuration (row-major: spin[i*M + m])
    float* d_local_field,     // [N x M] local field
    const float* d_str_c,     // [M] precomputed coupling strengths
    float neg_inv_T,          // -1/T for Metropolis
    int N, int M,
    int block_start,          // Starting node index for this block
    int block_size,           // Number of nodes in this block
    unsigned long long seed,  // Base seed for RNG
    int t_mc,                 // Current MC step (for seed variation)
    int* d_flipped            // [block_size x M] output: which spins flipped
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int num_elements = block_size * M;

    if (tid >= num_elements) return;

    int local_i = tid / M;
    int m = tid % M;
    int i = block_start + local_i;

    if (i >= N) return;

    int spin_idx = i * M + m;

    // Xorshift64 RNG with unique seed per thread
    unsigned long long rng = seed ^ (tid * 2654435761ULL) ^ (t_mc * 1099511628211ULL);
    rng ^= rng >> 12;
    rng ^= rng << 25;
    rng ^= rng >> 27;
    float rand_val = ((rng * 0x2545F4914F6CDD1DULL) >> 40) * (1.0f / 16777216.0f);

    // Trotter coupling (interaction with adjacent Trotter slices)
    int trotter_coupling = 0;
    if (m > 0) trotter_coupling += d_spin[spin_idx - 1];
    if (m < M - 1) trotter_coupling += d_spin[spin_idx + 1];

    // Delta energy
    float delta_H = d_spin[spin_idx] * (d_local_field[spin_idx] - d_str_c[m] * trotter_coupling);

    // Metropolis acceptance
    int flip = (expf(delta_H * neg_inv_T) > rand_val) ? 1 : 0;
    d_flipped[tid] = flip;

    if (flip) {
        d_spin[spin_idx] = -d_spin[spin_idx];
    }
}

// CUDA kernel: Prepare the "judged_spin" matrix for GEMM
// judged_spin[k, m] = 2 * spin[block_start + k, m] if flipped, else 0
__global__ void prepare_judged_spin_kernel(
    const int* d_spin,
    const int* d_flipped,
    half* d_judged_spin,  // [blk_sz x M] in half precision for GEMM
    int N, int M,
    int block_start, int block_size
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int num_elements = block_size * M;

    if (tid >= num_elements) return;

    int local_k = tid / M;
    int m = tid % M;
    int i = block_start + local_k;

    if (i >= N) {
        d_judged_spin[tid] = __float2half(0.0f);
        return;
    }

    if (d_flipped[tid]) {
        // Spin was flipped, contribute 2 * new_spin
        float val = 2.0f * d_spin[i * M + m];
        d_judged_spin[tid] = __float2half(val);
    } else {
        d_judged_spin[tid] = __float2half(0.0f);
    }
}

// CUDA kernel: Add GEMM result to local_field
__global__ void add_to_local_field_kernel(
    float* d_local_field,
    const float* d_gemm_output,
    int n
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        d_local_field[idx] += d_gemm_output[idx];
    }
}

// Compute Ising energy for one Trotter slice
int compute_energy_slice(const int* spin, const vector<vector<int>>& J, int N, int M, int m) {
    int energy = 0;
    for (int i = 0; i < N; i++) {
        int si = (spin[i * M + m] + 1) / 2;  // Convert {-1,1} to {0,1}
        for (int j = 0; j < N; j++) {
            int sj = (spin[j * M + m] + 1) / 2;
            energy += J[i][j] * si * sj;
        }
    }
    return -energy;
}

// Find best spin configuration across all Trotter slices
int choice_spin_gpu(const int* spin, const vector<vector<int>>& J, int N, int M) {
    int max_cut = -1;
    for (int m = 0; m < M; m++) {
        int cut = -compute_energy_slice(spin, J, N, M, m);
        if (cut > max_cut) max_cut = cut;
    }
    return max_cut;
}

int main() {
    string filename = "G1";

    // Read graph
    ReadGraph J_Data(filename);
    vector<Edge> edges = J_Data.getEdges();
    int numNodes = J_Data.getNumNodes();

    // Build QUBO matrix
    vector<vector<int>> graph = edgeListToQUBO(edges, numNodes);

    int N = numNodes;
    int M = 32;              // Trotter slices
    int MC_Step = 1024 * 7;  // Monte Carlo steps
    int blk_sz = 64;         // Block size (K in the algorithm)

    cout << "========================================" << endl;
    cout << "[GPU - cuBLAS GEMM + Tensor Core]" << endl;
    cout << "========================================" << endl;
    cout << "Gset: " << filename << endl;
    cout << "Graph Size (N): " << N << endl;
    cout << "Trotter Slices (M): " << M << endl;
    cout << "MC Steps: " << MC_Step << endl;
    cout << "Block Size (K): " << blk_sz << endl;

    // Initialize cuBLAS
    cublasHandle_t cublas_handle;
    CUBLAS_CHECK(cublasCreate(&cublas_handle));

    // Enable Tensor Core operations
    CUBLAS_CHECK(cublasSetMathMode(cublas_handle, CUBLAS_TENSOR_OP_MATH));

    // Flatten QUBO matrix to 1D (row-major)
    vector<int> J_flat(N * N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            J_flat[i * N + j] = graph[i][j];
        }
    }

    // Initialize spins randomly {-1, +1}
    vector<int> h_spin(N * M);
    random_device rd;
    mt19937 rng(rd());
    uniform_int_distribution<int> dis(0, 1);
    for (int i = 0; i < N * M; i++) {
        h_spin[i] = dis(rng) * 2 - 1;
    }

    // Compute initial local field: local_field[i,m] = sum_j J[i,j] * spin[j,m]
    // This is equivalent to: local_field = J @ spin (matrix multiply)
    vector<float> h_local_field(N * M, 0.0f);
    for (int i = 0; i < N; i++) {
        for (int m = 0; m < M; m++) {
            float sum = 0.0f;
            for (int j = 0; j < N; j++) {
                sum += J_flat[i * N + j] * h_spin[j * M + m];
            }
            h_local_field[i * M + m] = sum;
        }
    }

    // Allocate device memory
    int *d_spin, *d_J_int, *d_flipped;
    float *d_local_field, *d_str_c, *d_gemm_output;
    half *d_J_half, *d_judged_spin;

    size_t spin_bytes = N * M * sizeof(int);
    size_t field_bytes = N * M * sizeof(float);
    size_t J_bytes = N * N * sizeof(int);
    size_t J_half_bytes = N * blk_sz * sizeof(half);  // Only need J[:, block] for each GEMM
    size_t judged_bytes = blk_sz * M * sizeof(half);
    size_t flipped_bytes = blk_sz * M * sizeof(int);

    CUDA_CHECK(cudaMalloc(&d_spin, spin_bytes));
    CUDA_CHECK(cudaMalloc(&d_J_int, J_bytes));
    CUDA_CHECK(cudaMalloc(&d_local_field, field_bytes));
    CUDA_CHECK(cudaMalloc(&d_str_c, M * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_flipped, flipped_bytes));
    CUDA_CHECK(cudaMalloc(&d_judged_spin, judged_bytes));
    CUDA_CHECK(cudaMalloc(&d_gemm_output, field_bytes));

    // Allocate J in half precision for GEMM (full matrix, columns will be selected)
    half *d_J_full_half;
    CUDA_CHECK(cudaMalloc(&d_J_full_half, N * N * sizeof(half)));

    // Copy initial data to device
    CUDA_CHECK(cudaMemcpy(d_spin, h_spin.data(), spin_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_J_int, J_flat.data(), J_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_local_field, h_local_field.data(), field_bytes, cudaMemcpyHostToDevice));

    // Convert J to half precision
    int threads = 256;
    int blocks_J = (N * N + threads - 1) / threads;
    convert_int_to_half<<<blocks_J, threads>>>(d_J_int, d_J_full_half, N * N);
    CUDA_CHECK(cudaDeviceSynchronize());

    // Precompute str_c values buffer
    vector<float> h_str_c(M);

    // Kernel launch config
    int blocks_spin = (blk_sz * M + threads - 1) / threads;
    int blocks_field = (N * M + threads - 1) / threads;

    // RNG seed
    unsigned long long base_seed = rd() | ((unsigned long long)rd() << 32);

    cout << "Starting MC simulation with cuBLAS GEMM..." << endl;
    auto start = high_resolution_clock::now();

    // GEMM parameters for local_field update:
    // local_field += J[:, block] @ judged_spin[block, :]
    // Dimensions: [N x M] += [N x K] @ [K x M]
    // cuBLAS uses column-major, so we need to be careful with transpositions
    //
    // Using cublasGemmEx with Tensor Cores:
    // C = alpha * A @ B + beta * C
    // For row-major data treated as column-major:
    // C^T = alpha * B^T @ A^T + beta * C^T
    //
    // We want: local_field[N,M] += J[N,K] @ judged[K,M]
    // In column-major view with our row-major data:
    // local_field_cm = judged_cm^T @ J_cm^T

    float alpha = 1.0f;
    float beta = 1.0f;  // Accumulate into output

    for (int t_mc = 1; t_mc <= MC_Step; t_mc++) {
        // Compute temperature
        float T = 1.0f / (t_mc * (1.0f - 1.0f/8.0f) / static_cast<float>(MC_Step));
        float neg_inv_T = -1.0f / T;

        // Precompute str_c for each Trotter index
        float mystery = 32.0f;
        float t_ratio = static_cast<float>(t_mc / MC_Step);  // integer division preserved
        for (int m = 0; m < M; m++) {
            h_str_c[m] = T / 2.0f * logf(coshf(mystery * (1.0f - t_ratio) / (m * T)));
        }
        CUDA_CHECK(cudaMemcpy(d_str_c, h_str_c.data(), M * sizeof(float), cudaMemcpyHostToDevice));

        // Process spins in blocks
        for (int block_start = 0; block_start < N; block_start += blk_sz) {
            int K = min(blk_sz, N - block_start);

            // 1. Metropolis update for this block
            CUDA_CHECK(cudaMemset(d_flipped, 0, flipped_bytes));
            metropolis_update_kernel<<<blocks_spin, threads>>>(
                d_spin, d_local_field, d_str_c,
                neg_inv_T, N, M,
                block_start, K,
                base_seed, t_mc,
                d_flipped
            );

            // 2. Prepare judged_spin matrix (only flipped spins contribute)
            CUDA_CHECK(cudaMemset(d_judged_spin, 0, judged_bytes));
            prepare_judged_spin_kernel<<<blocks_spin, threads>>>(
                d_spin, d_flipped, d_judged_spin,
                N, M, block_start, K
            );

            // 3. Update local_field using GEMM: local_field += J[:, block] @ judged_spin
            // Using cublasGemmEx with FP16 inputs, FP32 accumulation, Tensor Cores
            //
            // Matrix dimensions:
            //   J[:, block_start:block_start+K] is [N x K]
            //   judged_spin is [K x M]
            //   Result is [N x M]
            //
            // cuBLAS column-major: we treat row-major as transposed
            // So we compute: result = judged_spin^T @ J_block^T (in column-major view)
            // Which gives us: [M x N] = [M x K] @ [K x N]

            // Get pointer to the relevant columns of J
            half* d_J_block = d_J_full_half + block_start;  // J[:, block_start], stride N

            // Zero out gemm output before accumulation
            CUDA_CHECK(cudaMemset(d_gemm_output, 0, field_bytes));

            // GEMM: d_gemm_output = d_J_block @ d_judged_spin
            // [N x M] = [N x K] @ [K x M]
            // In cuBLAS column-major: C = A @ B
            // m=M, n=N, k=K, A=judged (K x M in col-maj = M x K in row-maj)
            // B=J_block (K x N in col-maj = N x K in row-maj)
            // C=output (N x M in col-maj = M x N in row-maj)

            CUBLAS_CHECK(cublasGemmEx(
                cublas_handle,
                CUBLAS_OP_N,        // op(A) = A
                CUBLAS_OP_N,        // op(B) = B
                M, N, K,            // m, n, k (output is m x n)
                &alpha,
                d_judged_spin, CUDA_R_16F, M,      // A: K x M (col-maj), lda=M
                d_J_block, CUDA_R_16F, N,          // B: K x N (col-maj), lda=N  (strided access)
                &beta,
                d_gemm_output, CUDA_R_32F, M,      // C: M x N (col-maj), ldc=M
                CUBLAS_COMPUTE_32F,
                CUBLAS_GEMM_DEFAULT_TENSOR_OP
            ));

            // 4. Add GEMM result to local_field (need to transpose the result)
            // Actually, let's use a simpler approach - update directly with atomic adds
            // For now, we'll use a direct kernel to update local_field

            add_to_local_field_kernel<<<blocks_field, threads>>>(
                d_local_field, d_gemm_output, N * M
            );

            CUDA_CHECK(cudaDeviceSynchronize());
        }
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start).count();

    // Copy results back
    CUDA_CHECK(cudaMemcpy(h_spin.data(), d_spin, spin_bytes, cudaMemcpyDeviceToHost));

    // Find best cut
    int maxcut = choice_spin_gpu(h_spin.data(), graph, N, M);

    cout << "========================================" << endl;
    cout << "Time taken: " << duration << " milliseconds" << endl;
    cout << "Max Cut Value: " << maxcut << endl;
    cout << "Best Known Value: " << my_graph[filename] << endl;
    cout << "Accuracy: " << maxcut / static_cast<float>(stoi(my_graph[filename])) << endl;
    cout << "========================================" << endl;

    // Cleanup
    CUBLAS_CHECK(cublasDestroy(cublas_handle));
    CUDA_CHECK(cudaFree(d_spin));
    CUDA_CHECK(cudaFree(d_J_int));
    CUDA_CHECK(cudaFree(d_local_field));
    CUDA_CHECK(cudaFree(d_str_c));
    CUDA_CHECK(cudaFree(d_flipped));
    CUDA_CHECK(cudaFree(d_judged_spin));
    CUDA_CHECK(cudaFree(d_gemm_output));
    CUDA_CHECK(cudaFree(d_J_full_half));

    return 0;
}
