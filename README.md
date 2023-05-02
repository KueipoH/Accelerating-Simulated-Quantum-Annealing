# Welcome to the Simulated Quantum Annealing repository

Simulated Quantum Annealing (SQA) is inspired by Quantum Annealing (QA) and aims to emulate the quantum tunneling effect on classical computers using a path-integral Monte Carlo simulation. This simulation introduces the strength of couplings between replicas.

In this project, I have referred to __Accelerating Simulated Quantum Annealing with GPU and Tensor Cores__ as a foundation. To quickly grasp the concept of SQA, I have initially developed a CPU-based version for validation. Following this, I will proceed to implement a CUDA version based on the literature.

## Getting Started

This project is a __CPU-based__ implementation of the research paper. _Accelerating Simulated Quantum Annealing with GPU and Tensor Cores_ by Author Yi-Hua Chung, Cheng-Jhih Shih, and Shih-Hao Hung. [Link to paper](https://link.springer.com/chapter/10.1007/978-3-031-07312-0_9).

Firstly, it is crucial to understand the concept of the Local-field matrix's dimensions. Upon initial reading, one might easily be misled by the Trotter diagrams in terms of comprehending the Local-field matrix's dimensional space. It is also important to consider why the author has chosen to implement the outer loop for spins and the inner loop for Trotters. This choice is the rationale behind updating the local-field energy. Finally, when it comes to implementation details, the crux lies in selecting the final spin. While theory provides a foundation, the secret is hidden within my very own code.

## Requirements & Command

- c++ 14 or higher
```cpp
g++ -std=c++14 main.cpp file.cpp -o sqa && ./sqa
```

## License
The code in this repository is licensed under [Apache License 2.0](https://github.com/NVIDIA/cuda-quantum/blob/main/LICENSE).
Contributing a pull request to this repo requires accepting the Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us the rights to use your contribution. A CLA-bot will automatically determine whether you need to provide a CLA and decorate the PR appropriately. Simply follow the instructions provided by the bot. You will only need to do this once.


## Results

Achieving over 98% accuracy on the G-set requires tuning the number of Trotter slices, Monte Carlo steps, and transverse field _G0_ param.
- A larger MC-STEP usually produces a better result closer to the global optimum.


## Feedback

Please let me know your feedback and ideas, the Cuda version is coming soon.


## References

- __Accelerating Simulated Quantum Annealing with GPU and Tensor Cores__, High Performance Computing: 37th International Conference, ISC High Performance 2022, Hamburg, Germany, May 29–June 2, 2022, Proceedings. 
- Dataset. [G-Set](https://web.stanford.edu/~yyye/yyye/Gset/)

