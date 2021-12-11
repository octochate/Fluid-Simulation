#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
#include <time.h>

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}

__global__ void add_source(float* x, float* s, float dt, int size)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    if (idx < size) {
        x[idx] += dt * s[idx];
    }
}

__device__ void set_bnd(int N, int b, float* x, int index, int elementsPerThread)
{
    int i = index + 1;
    if (i > N + 1) {
        return;
    }

    int size = (N + 2) * (N + 2);

    while (i < (index + elementsPerThread) && i <= N + 1) {
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
        i++;
    }

    __syncthreads();

    if (index == 0) {
        x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }
}

__global__ void lin_solve(int N, int b, float* x, float* x0, float a, float c, int elementPerThread)
{
    int i, j, k, idxNew;
    int idx = (blockIdx.x * blockDim.x + threadIdx.x) * elementPerThread;

    if (idx >= N * N) {
        return;
    }

    for (k = 0; k < 20; k++)
    {
        idxNew = idx;
        while (idxNew < idx + elementPerThread && idxNew < N * N) {
            j = idxNew / N + 1;
            i = idxNew % N + 1;
            x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
            idxNew++;
        }
        __syncthreads();
        set_bnd(N, b, x, idx, elementPerThread);
    }
}

__global__ void advect(int N, int b, float* d, float* d0, float* u, float* v, float dt, int elementPerThread)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x) * elementPerThread;

    if (idx >= N * N) {
        return;
    }

    int idxNew = idx;
    while (idxNew < idx + elementPerThread && idxNew < N * N) {
        int j = idxNew / N + 1;
        int i = idxNew % N + 1;

        int i0, j0, i1, j1;
        float x, y, s0, t0, s1, t1, dt0;

        dt0 = dt * N;
        x = i - dt0 * u[IX(i, j)]; y = j - dt0 * v[IX(i, j)];
        if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; i0 = (int)x; i1 = i0 + 1;
        if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; j0 = (int)y; j1 = j0 + 1;
        s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;

        d[IX(i, j)] =  s0* (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
            s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);

        idxNew++;
    }
    
    __syncthreads();
    set_bnd(N, b, d, idx, elementPerThread);
}

__global__ void project1(int N, float* u, float* v, float* p, float* div, int elementPerThread)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x) * elementPerThread;

    if (idx >= N * N) {
        return;
    }

    int idxNew = idx;
    while (idxNew < idx + elementPerThread && idxNew < N * N) {
        int j = idxNew / N + 1;
        int i = idxNew % N + 1;

        div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
        p[IX(i, j)] = 0;
        idxNew++;
    }

    __syncthreads();
    set_bnd(N, 0, div, idx, elementPerThread);
    set_bnd(N, 0, p, idx, elementPerThread);
}

__global__ void project3(int N, float* u, float* v, float* p, int elementPerThread)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x) * elementPerThread;

    if (idx >= N * N) {
        return;
    }

    int idxNew = idx;
    while (idxNew < idx + elementPerThread && idxNew < N * N) {
        int j = idxNew / N + 1;
        int i = idxNew % N + 1;

        u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
        v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);

        idxNew++;
    }
    __syncthreads();
    set_bnd(N, 1, u, idx, elementPerThread);
    set_bnd(N, 2, v, idx, elementPerThread);
}

void get_from_UI(float* d, float* u, float* v, int force, int source, int N)
{
    int i, j, size = (N + 2) * (N + 2);

    for (i = 0; i < size; i++)
    {
        u[i] = v[i] = d[i] = 0.0f;
    }

    i = (N+2) / 2;
    j = (N+2) / 2;

    if (i < 1 || i > N || j < 1 || j > N)
        return;

    u[IX(i, j)] = force * (1.0);
    v[IX(i, j)] = force * (1.0);
    d[IX(i, j)] = source * 1.0;

    return;
}

static void printFrameMatrices(float* dens, float* u, float* v, int N)
{
    printf("Density Matrix:\t\t\t\t\t\t\tVelocity U Matrix:\t\t\t\t\t\tVelocity V Matrix:\n");
    for (int i = 0; i <= N + 1; i++) {
        for (int j = 0; j <= N + 1; j++)
        {
            printf("%f, ", dens[IX(i, j)]);
        }
        printf("\t");
        for (int j = 0; j <= N + 1; j++)
        {
            printf("%f, ", u[IX(i, j)]);
        }
        printf("\t");
        for (int j = 0; j <= N + 1; j++)
        {
            printf("%f, ", v[IX(i, j)]);
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]) {

    // declare variables
    float dt = 0.1f, diff = 0.1f, visc = 0.1f;
    float force = 100.0f, source = 100.0f;
    int N = 5, num_threads = 2, num_iterations = 10, display_output = 0, use_streams = 0;

    // load parameters from command line
    if (argc != 6) {
        fprintf(stderr, "usage : %s <grid size> <number of threads> <number of iterations> <display output> <use CUDA streams>\n", argv[0]);
        fprintf(stderr, "where:\n"); \
        fprintf(stderr, "\t grid_size                : grid resolution\n");
        fprintf(stderr, "\t number_of_threads        : number of threads running per kernel (1 block/kernel enforced)\n");
        fprintf(stderr, "\t number_of_iterations     : number of iterations to run the simulation\n");
        fprintf(stderr, "\t display_output           : Whether to display the opdated grids after every iteration\n");
        fprintf(stderr, "\t use_CUDA_streams         : Whether to use CUDA streams or have all kernels run sequentially\n");
        exit(1);
    }


    N = atoi(argv[1]);
    num_threads = atoi(argv[2]);
    num_iterations = atoi(argv[3]);
    display_output = atoi(argv[4]);
    use_streams = atoi(argv[5]);

    // validate parameters
    if (N <= 0) {
        fprintf(stderr, "grid_size must be an integer larger than 0\n");
        exit(1);
    }
    if (num_threads <= 0 || num_threads > 1024) {
        fprintf(stderr, "number_of_threads must be an integer in the range [1, 1024]\n");
        exit(1);
    }
    if (num_iterations <= 0) {
        fprintf(stderr, "number_of_iterations must be an integer larger than 0\n");
        exit(1);
    }
    if (display_output != 0 && display_output != 1) {
        fprintf(stderr, "display_output must be 1 if TRUE and 0 if FALSE\n");
        exit(1);
    }
    if (use_streams != 0 && use_streams != 1) {
        fprintf(stderr, "use_streams must be 1 if TRUE and 0 if FALSE\n");
        exit(1);
    }

    
    int size = (N + 2) * (N + 2);

    // unified memory: pointers to data
    float* u_cuda, * v_cuda, * u_prev_cuda, * v_prev_cuda;
    float* p_cuda, * div_cuda;
    float* dens_cuda, * dens_prev_cuda;

    // Allocate space for arrays
    cudaMallocManaged(&u_cuda, sizeof(float) * size);
    cudaMallocManaged(&v_cuda, sizeof(float) * size);
    cudaMallocManaged(&p_cuda, sizeof(float) * size);
    cudaMallocManaged(&div_cuda, sizeof(float) * size);
    cudaMallocManaged(&u_prev_cuda, sizeof(float) * size);
    cudaMallocManaged(&v_prev_cuda, sizeof(float) * size);
    cudaMallocManaged(&dens_cuda, sizeof(float) * size);
    cudaMallocManaged(&dens_prev_cuda, sizeof(float) * size);

    // Initialize all arrays to 0
    cudaMemset(u_cuda, 0, sizeof(float) * size);
    cudaMemset(v_cuda, 0, sizeof(float) * size);
    cudaMemset(p_cuda, 0, sizeof(float) * size);
    cudaMemset(div_cuda, 0, sizeof(float) * size);
    cudaMemset(u_prev_cuda, 0, sizeof(float) * size);
    cudaMemset(v_prev_cuda, 0, sizeof(float) * size);
    cudaMemset(dens_cuda, 0, sizeof(float) * size);
    cudaMemset(dens_prev_cuda, 0, sizeof(float) * size);


    // Parallelization scheme for add_source
    int num_threads_source = (N + 2);
    int num_blocks_source = (N + 2);
    
    // Parallelization scheme for only using one block
    int elementsPerThread = (N * N + 1) / num_threads;

    // Instantiate cuda streams
    cudaStream_t stream1, stream2, stream3;
    cudaStreamCreate(&stream1);
    if (use_streams) {
        cudaStreamCreate(&stream2);
        cudaStreamCreate(&stream3);
    }
    else {
        stream2 = stream1;
        stream3 = stream1;
    }
    


    // update loop
    for (int i = 0; i < num_iterations; i++) {

        // get input from 'GUI'
        get_from_UI(dens_prev_cuda, u_prev_cuda, v_prev_cuda, force, source, N);


        //// Add GUI inputs to grids
        add_source <<< num_blocks_source, num_threads_source >>> (u_cuda, u_prev_cuda, dt, size);
        add_source <<< num_blocks_source, num_threads_source >>> (v_cuda, v_prev_cuda, dt, size);
        add_source <<< num_blocks_source, num_threads_source >>> (dens_cuda, dens_prev_cuda, dt, size);
        cudaDeviceSynchronize();

        SWAP(u_prev_cuda, u_cuda);
        SWAP(v_prev_cuda, v_cuda);
        SWAP(dens_prev_cuda, dens_cuda);

        //// diffuse step
        float a = dt * visc * N * N;
        lin_solve << < 1, num_threads, 0, stream1 >> > (N, 1, u_cuda, u_prev_cuda, a, 1 + 4 * a, elementsPerThread); // diffuse u_cuda
        lin_solve << < 1, num_threads, 0, stream2 >> > (N, 2, v_cuda, v_prev_cuda, a, 1 + 4 * a, elementsPerThread); // diffuse v_cuda
        a = dt * diff * N * N;
        lin_solve << < 1, num_threads, 0, stream3 >> > (N, 0, dens_cuda, dens_prev_cuda, a, 1 + 4 * a, elementsPerThread); // diffuse dens_cuda
        cudaDeviceSynchronize();


        // projection step (no swapping beforehand)
        project1 << < 1, num_threads >> > (N, u_cuda, v_cuda, p_cuda, div_cuda, elementsPerThread);
        lin_solve << < 1, num_threads >> > (N, 0, p_cuda, div_cuda, 1, 4, elementsPerThread);
        project3 << < 1, num_threads >> > (N, u_cuda, v_cuda, p_cuda, elementsPerThread);
        cudaDeviceSynchronize();

        SWAP(u_prev_cuda, u_cuda);
        SWAP(v_prev_cuda, v_cuda);

        // advection step (velocity grid)
        advect << < 1, num_threads, 0, stream1 >> > (N, 1, u_cuda, u_prev_cuda, u_prev_cuda, v_prev_cuda, dt, elementsPerThread);
        advect << < 1, num_threads, 0, stream2 >> > (N, 2, v_cuda, v_prev_cuda, u_prev_cuda, v_prev_cuda, dt, elementsPerThread);
        cudaDeviceSynchronize();


        // projection step (no swapping beforehand)
        project1 << < 1, num_threads >> > (N, u_cuda, v_cuda, p_cuda, div_cuda, elementsPerThread);
        lin_solve << < 1, num_threads >> > (N, 0, p_cuda, div_cuda, 1, 4, elementsPerThread);
        project3 << < 1, num_threads >> > (N, u_cuda, v_cuda, p_cuda, elementsPerThread);
        cudaDeviceSynchronize();
    

        // advection step (density grid)
        SWAP(dens_prev_cuda, dens_cuda);
        advect << < 1, num_threads >> > (N, 0, dens_cuda, dens_prev_cuda, u_cuda, v_cuda, dt, elementsPerThread);
        cudaDeviceSynchronize();

        // Display results if needed
        printf("Finished iteration %d:\n", i);
        if (display_output) {
            printFrameMatrices(dens_cuda, u_cuda, v_cuda, N);
            puts("\n");
        }
    }

    // Free Device space
    cudaFree(u_cuda);
    cudaFree(v_cuda);
    cudaFree(u_prev_cuda);
    cudaFree(v_prev_cuda);
    cudaFree(dens_cuda);
    cudaFree(dens_prev_cuda);
    cudaFree(p_cuda);
    cudaFree(div_cuda);
    

    // Destroy cuda streams
    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaStreamDestroy(stream3);

    return 0;
}
