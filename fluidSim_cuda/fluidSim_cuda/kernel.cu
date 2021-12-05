#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>
#include <cuda_runtime_api.h>


#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}

__global__ void add_source(float* x, float* s, float dt, int size)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    if (idx < size) {
        x[idx] += dt * s[idx];
    }
}

__device__ void set_bnd(int N, int b, float* x, int index)
{
    if (index > N) {
        return;
    }
    int i = index;

    x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
    x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
    x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
    x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
    
    if (index == 0) {
        x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }
}
__global__ void lin_solve(int N, int b, float* x, float* x0, float a, float c)
{
    int i, j, k;
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);
    if (idx > (N + 2) * (N + 2)) {
        return;
    }

    j = idx / N;
    i = idx % N;

    for (k = 0; k < 20; k++)
    {
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
        __syncthreads();
         set_bnd(N, b, x, idx);
        __syncthreads();
    }
}


__global__ void advect(int N, int b, float* d, float* d0, float* u, float* v, float dt)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    if (idx > (N + 2) * (N + 2)) {
        return;
    }

    int j = idx / N;
    int i = idx % N;

    int i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * N;
        x = i - dt0 * u[IX(i, j)]; y = j - dt0 * v[IX(i, j)];
    if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; i0 = (int)x; i1 = i0 + 1;
    if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; j0 = (int)y; j1 = j0 + 1;
    s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
    d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
        s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);

    __syncthreads();
    set_bnd(N, b, d, idx);
}


__global__ void project1(int N, float* u, float* v, float* p, float* div)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    if (idx > (N + 2) * (N + 2)) {
        return;
    }

    int j = idx / N;
    int i = idx % N;

        div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
        p[IX(i, j)] = 0;

        __syncthreads();
        set_bnd(N, 0, div, idx); 
        __syncthreads();
        set_bnd(N, 0, p, idx);
}

__global__ void project3(int N, float* u, float* v, float* p, float* div)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    if (idx > (N + 2) * (N + 2)) {
        return;
    }

    int j = idx / N;
    int i = idx % N;

        u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
    v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);

    __syncthreads();
    set_bnd(N, 1, u, idx); 
    __syncthreads();
    set_bnd(N, 2, v, idx);
}

void get_from_UI(float* d, float* u, float* v, int force, int source, int N)
{
    int i, j, size = (N + 2) * (N + 2);

    for (i = 0; i < size; i++)
    {
        u[i] = v[i] = d[i] = 0.0f;
    }

    i = N / 2;
    j = N / 2;

    if (i < 1 || i > N || j < 1 || j > N)
        return;

    u[IX(i, j)] = force * (1.0);
    v[IX(i, j)] = force * (1.0);
    d[IX(i, j)] = source * 1.0;

    return;
}

static void printFrameMatrices(float* dens, float* u, float* v, int N)
{
        printf("Density Matrix:\t\t\t\t\tVelocity U Matrix:\t\t\t\t\tVelocity V Matrix:\n");
        for (int j = 0; j < 1; j++)
        {
            printf("%f, ", dens[j * N]);
            printf("%f, ", dens[1 + j * N]);
            printf("%f, ", dens[2 + j * N]);
            printf("%f", dens[3 + j * N]);
            printf("\t\t");
            printf("%f, ", u[j * N]);
            printf("%f, ", u[1 + j * N]);
            printf("%f, ", u[2 + j * N]);
            printf("%f", u[3 + j * N]);
            if (j == 1)
            {
                printf("\t");
            }
            else
            {
                printf("\t\t");
            }
            if (j != 1)
            {
                printf("%f, ", v[j * N]);
                printf("%f, ", v[1 + j * N]);
                printf("%f, ", v[2 + j * N]);
                printf("%f", v[3 + j * N]);
            }
            else
            {
                printf("\t%f, ", v[j * N]);
                printf("%f, ", v[1 + j * N]);
                printf("%f, ", v[2 + j * N]);
                printf("%f", v[3 + j * N]);
            }
            printf("\n");
        }
        printf("\n");
}

int main(int argc, char* argv[]) {

    int N = 4;
    float dt = 0.1f, diff = 0.0f, visc = 0.0f;
    float force = 0.5f, source = 100.0f;

    // Allocate space for host
    int size = (N + 2) * (N + 2);

    // Device copies of data
    float* u_cuda, * v_cuda, * u_prev_cuda, * v_prev_cuda;
    float* dens_cuda, * dens_prev_cuda;

    // Allocate space for device copies
    cudaMallocManaged((void**)&u_cuda, sizeof(float));
    cudaMallocManaged((void**)&v_cuda, sizeof(float));
    cudaMallocManaged((void**)&u_prev_cuda, sizeof(float));
    cudaMallocManaged((void**)&v_prev_cuda, sizeof(float));
    cudaMallocManaged((void**)&dens_cuda, sizeof(float));
    cudaMallocManaged((void**)&dens_prev_cuda, sizeof(float));

    // have copy of host data
    get_from_UI(dens_cuda, u_cuda, v_cuda, force, source, N);


    // Parallelize computation :: TODO update values to match
    int num_threads = (N + 2);
    int num_blocks = (N + 2);

    cudaStream_t stream1, stream2, stream3;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);
    cudaStreamCreate(&stream3);

    // Velocity timestep parallelization
    add_source <<< num_blocks, num_threads, 0, stream1>>> (u_cuda, u_prev_cuda, dt, size);
    add_source <<< num_blocks, num_threads, 0, stream2>>> (v_cuda, v_prev_cuda, dt, size);
    add_source <<< num_blocks, num_threads, 0, stream3>>> (dens_cuda, dens_prev_cuda, dt, size);
    cudaDeviceSynchronize();

    SWAP(u_prev_cuda, u_cuda);
    SWAP(v_prev_cuda, v_cuda);
    SWAP(dens_prev_cuda, dens_cuda);

    // diffuse
    float a = dt * visc * N * N;
    lin_solve << < num_blocks, num_threads >> > (N, 1, u_prev_cuda, u_cuda, a, 1 + 4 * a); // diffuse u_cuda
    cudaDeviceSynchronize();

    lin_solve << < num_blocks, num_threads >> > (N, 2, v_prev_cuda, v_cuda, a, 1 + 4 * a); // diffuse v_cuda
    a = dt * diff * N * N;
    lin_solve << < num_blocks, num_threads >> > (N, 0, dens_prev_cuda, dens_cuda, a, 1 + 4 * a); // diffuse dens_cuda
    cudaDeviceSynchronize();

    project1 << < num_blocks, num_threads >> > (N, u_cuda, v_cuda, u_prev_cuda, v_prev_cuda);
    cudaDeviceSynchronize();

    lin_solve << < num_blocks, num_threads >> > (N, 0, u_prev_cuda, v_prev_cuda, 1, 4);
    cudaDeviceSynchronize();

    project3 << < num_blocks, num_threads >> > (N, u_cuda, v_cuda, u_prev_cuda, v_prev_cuda);
    cudaDeviceSynchronize();

    SWAP(dens_prev_cuda, dens_cuda);
    SWAP(u_prev_cuda, u_cuda);
    SWAP(v_prev_cuda, v_cuda);


    advect << < num_blocks, num_threads >> > (N, 1, u_cuda, u_prev_cuda, u_prev_cuda, v_prev_cuda, dt);
    advect << < num_blocks, num_threads >> > (N, 2, v_cuda, v_prev_cuda, u_prev_cuda, v_prev_cuda, dt);
    cudaDeviceSynchronize();

    project1 << < num_blocks, num_threads >> > (N, u_cuda, v_cuda, u_prev_cuda, v_prev_cuda);
    cudaDeviceSynchronize();
    lin_solve << < num_blocks, num_threads >> > (N, 0, u_prev_cuda, v_prev_cuda, 1, 4);
    cudaDeviceSynchronize();
    project3 << < num_blocks, num_threads >> > (N, u_cuda, v_cuda, u_prev_cuda, v_prev_cuda);
    cudaDeviceSynchronize();

    // Density timestep parallelization
    advect << < num_blocks, num_threads >> > (N, 0, dens_cuda, dens_prev_cuda, u_cuda, v_cuda, dt);
    cudaDeviceSynchronize();

    // Copy result back to host
    printf("%f", dens_cuda[0]);

    // Free Device space
    cudaFree(u_cuda);
    cudaFree(v_cuda);
    cudaFree(u_prev_cuda);
    cudaFree(dens_cuda);
    cudaFree(dens_prev_cuda);
    cudaFree(v_prev_cuda);

    return 0;
}
