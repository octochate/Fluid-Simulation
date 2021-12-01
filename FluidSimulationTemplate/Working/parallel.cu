#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

__device__ void add_source(int N, float* x, float* s, float dt, int size)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    if (idx < size) {
        x[idx] += dt * s[idx];
    }
}

__device__ void set_bnd(int N, int b, float* x)
{
    int i;

    for (i = 1; i <= N; i++)
    {
        x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
    }
    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
    x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}
__device__ void lin_solve(int N, int b, float* x, float* x0, float a, float c)
{
    int i, j, k;
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);
    j = idx / N;
    i = idx % N;

    for (k = 0; k < 20; k++)
    {
         x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
         __syncthreads();
         if (idx == 0) {
             set_bnd(N, b, x);
         }
    }
}

__global__ void diffuse(int N, int b, float* x, float* x0, float diff, float dt)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);
    float a = dt * diff * N * N;
    lin_solve(N, b, x, x0, a, 1 + 4 * a);
}

__global__ void advect(int N, int b, float* d, float* d0, float* u, float* v, float dt)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * N;
    FOR_EACH_CELL
        x = i - dt0 * u[IX(i, j)]; y = j - dt0 * v[IX(i, j)];
    if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; i0 = (int)x; i1 = i0 + 1;
    if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; j0 = (int)y; j1 = j0 + 1;
    s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
    d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
        s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
    END_FOR
        set_bnd(N, b, d);
}


__global__ void project(int N, float* u, float* v, float* p, float* div)
{
    int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    int i, j;

    FOR_EACH_CELL
        div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
    p[IX(i, j)] = 0;
    END_FOR
        set_bnd(N, 0, div); set_bnd(N, 0, p);

    lin_solve(N, 0, p, div, 1, 4);

    FOR_EACH_CELL
        u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
    v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
    END_FOR
        set_bnd(N, 1, u); set_bnd(N, 2, v);
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

    u[IX(i, j)] = force * (1);
    v[IX(i, j)] = force * (1);
    d[IX(i, j)] = source;

    return;
}

int main(int argc, char* argv[]) {
       
    int N = 4;
    float dt = 0.1f, diff = 0.0f, visc = 0.0f;
    float force = 0.5f, source = 100.0f;

    // Host copies of data
    float* u, * v, * u_prev, * v_prev;
    float* dens, * dens_prev;

    // Allocate space for host
    int size = (N + 2) * (N + 2);
    u = (float*)malloc(size * sizeof(float));
    v = (float*)malloc(size * sizeof(float));
    u_prev = (float*)malloc(size * sizeof(float));
    v_prev = (float*)malloc(size * sizeof(float));
    dens = (float*)malloc(size * sizeof(float));
    dens_prev = (float*)malloc(size * sizeof(float));

    // Device copies of data
    float* u_cuda, * v_cuda, * u_prev_cuda, * v_prev_cuda;
    float* dens_cuda, * dens_prev_cuda;
    
    // Allocate space for device copies
    cudaMalloc((void**)&u_cuda, sizeof(float));
    cudaMalloc((void**)&v_cuda, sizeof(float));
    cudaMalloc((void**)&u_prev_cuda, sizeof(float));
    cudaMalloc((void**)&v_prev_cuda, sizeof(float));
    cudaMalloc((void**)&dens_cuda, sizeof(float));
    cudaMalloc((void**)&dens_prev_cuda, sizeof(float));

    // have copy of host data
    get_from_UI(dens, u_cuda, v_cuda, force, source, N);

    // Copy data from host to device
    cudaMemcpy(u_cuda, u, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(v_cuda, u, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(u_prev_cuda, u, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(v_prev_cuda, u, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dens, u, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dens_prev_cuda, u, size * sizeof(float), cudaMemcpyHostToDevice);

    // Parallelize computation :: TODO update values to match
    int num_threads = 1024;
    int num_blocks = 1;
    
    // Velocity timestep parallelization
    add_source(N, u_cuda, u_prev_cuda, dt, size);
    add_source(N, v_cuda, v_prev_cuda, dt, size);
    add_source(N, dens_cuda, dens_prev_cuda, dt, size);
    cudaDeviceSynchronize();

    SWAP(u_prev_cuda, u_cuda);
    SWAP(v_prev_cuda, v_cuda);
    SWAP(dens_prev_cuda, dens_cuda);
    cudaDeviceSynchronize();

    diffuse(N, 1, u_cuda, u_prev_cuda, visc, dt);
    diffuse(N, 2, v_cuda, v_prev_cuda, visc, dt);
    diffuse(N, 0, dens_cuda, dens_prev_cuda, diff, dt);
    cudaDeviceSynchronize();

    project(N, u_cuda, v_cuda, u_prev_cuda, v_prev_cuda);
    SWAP(dens_prev_cuda, dens_cuda);
    cudaDeviceSynchronize();

    SWAP(u_prev_cuda, u_cuda);
    SWAP(v_prev_cuda, v_cuda);
    cudaDeviceSynchronize();

    advect(N, 1, u_cuda, u_prev_cuda, u_prev_cuda, v_prev_cuda, dt);
    advect(N, 2, v_cuda, v_prev_cuda, u_prev_cuda, v_prev_cuda, dt);
    cudaDeviceSynchronize();

    project(N, u_cuda, v_cuda, u_prev_cuda, v_prev_cuda);
    cudaDeviceSynchronize();

    // Density timestep parallelization
    advect(N, 0, dens_cuda, dens_prev_cuda, u_cuda, v_cuda, dt);
    cudaDeviceSynchronize();

    // Copy result back to host

    // Free Device space
    cudaFree(u_cuda);
    cudaFree(v_cuda);
    cudaFree(u_prev_cuda);
    cudaFree(dens_cuda);
    cudaFree(dens_prev_cuda);
    cudaFree(v_prev_cuda);

    // Free Host space
    free(u); 
    free(v); 
    free(u_prev);
    free(v_prev); 
    free(dens); 
    free(dens_prev);

	return 0;
}
