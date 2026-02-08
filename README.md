# QKD_LDPC_V

This repository contains a C++ implementation of protocol using LDPC (Low-Density Parity-Check) codes for Quantum Key Distribution (QKD). It supports various decoding algorithms, code rate adaptation, privacy maintenance, and matrices in various formats.

## Building and Running Without Docker

To build and run the project locally without Docker, ensure you have CMake (version 3.30 or higher), g++, and build-essential installed on your Ubuntu system (or equivalent dependencies on other OS).

1. Clone the repository:
   ```
   git clone https://github.com/ColdCloudd/QKD_LDPC_V.git
   cd QKD_LDPC_V
   ```

2. Create a build directory and compile:
   ```
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j   # Use parallel build for speed
   ```

3. Run the executable (from the build directory):
   ```
   ./QKD_LDPC
   ```
   - The program expects configuration files in the `configs` directory, parity check matrices in the `sparse_matrices` directory, and outputs to `results`.

## Building and Running with Docker Compose

This project uses a multi-stage Dockerfile for building an optimized runtime image and docker-compose.yml for easy container management.

1. Clone the repository (if not already done):
   ```
   git clone https://github.com/ColdCloudd/QKD_LDPC_V.git
   cd QKD_LDPC_V
   ```

2. Build the Docker image:
   ```
   sudo docker compose build
   ```

3. Run the container:
   ```
   sudo docker compose run --rm qkd_ldpc
   ```
   - Use Ctrl-C to interrupt program execution.

## Usage

The program is a command-line tool for simulating LDPC-based error correction in QKD. It reads configurations from files in `configs` (if the folder contains several configuration files, the program will process them all sequentially) and processes parity check matrices from `sparse_matrices`, saving results to `results`.

To view the full help message with all configuration parameters:
```
./QKD_LDPC --help   #without Docker
```
```
sudo docker compose run --rm qkd_ldpc --help   #with Docker Compose
```
(or `-help`).
