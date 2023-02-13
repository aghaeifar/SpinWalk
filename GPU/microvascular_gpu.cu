/* --------------------------------------------------------------------------
 * Project: Microvascular
 * File: microvascular_gpu.cu
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 10.02.2023
 * Descrip  : simulating BOLD in microvascular network
 * -------------------------------------------------------------------------- */

// compile :  nvcc microvascular_gpu.cu -Xptxas -v -O3  -arch=compute_86 -code=sm_86  -Xcompiler -fopenmp

#include <sstream>
#include <random>

#include <thrust/random.h>

#include "../common/miscellaneous.h"
#include "../common/helper_cuda.h"
#include "../common/rotation.h"



#define GAMMA           267515315. // rad/s.T
#define MAX_N_FIELDMAP  3
#define CONFIG_FILE     "../inputs/config.ini"
#define THREADS_PER_BLOCK       64


__global__ void simulation_kernel_cuda(simulation_parameters *param_orig, 
                                       float *d_pFieldMap, 
                                       bool *d_pMask,
                                       float *d_random_walk_xyz_init_scaled,  // 3 * param.n_spins
                                       float *spin_mxyz);
__global__ void generate_initial_position_cuda(float *, simulation_parameters *, bool *);
__global__ void scale_initial_positions(float *, float *, float, int);

using namespace std;

int main()
{
    map<string, vector<string> > filenames = {{"fieldmap", vector<string>()},
                                              {"mask", vector<string>()},
                                              {"output", vector<string>()} }; 

    std::vector<float> sample_length_all;

    simulation_parameters param;
    float *pFieldMap = NULL;
    bool *pMask = NULL;

    // ========== read config file ==========
    if(read_config(CONFIG_FILE, param, sample_length_all, filenames) == false)
    {
        std::cout << "Reading config file failed. Aborting...!" << std::endl;
        return 1;
    }
    if(param.n_fieldmaps > MAX_N_FIELDMAP)
    {
        std::cout << "Due to some memory issues, we don't support more than 3 field-maps at the moment" << std::endl;
        return 1;
    }
    param.seed = std::random_device{}();

    // ========== load field-maps ==========
    std::ifstream in(filenames.at("fieldmap")[0], std::ios::in | std::ios::binary);
    in.read((char *)&param.fieldmap_size[0], sizeof(int) * 3);
    in.read((char *)&param.sample_length, sizeof(float));
    in.close();
    param.scale2grid = (param.fieldmap_size[0] - 1.) / param.sample_length;  
    param.matrix_length = param.fieldmap_size[0] * param.fieldmap_size[1] * param.fieldmap_size[2];
    if (param.fieldmap_size[0] != param.fieldmap_size[1] || param.fieldmap_size[0] != param.fieldmap_size[2])
    {
        std::cout << "Simulator does not support non-isotropic sample for now! we will fix it in far future :)";
        return 1;
    }

    // concatenate all fieldmaps
    pFieldMap = new float[param.matrix_length * param.n_fieldmaps];
    for (int i = 0; i < param.n_fieldmaps; i++)
    {
        std::ifstream in_field(filenames.at("fieldmap")[i], std::ios::in | std::ios::binary);
        in_field.seekg(sizeof(int) * 3 + sizeof(float) * 2); // skip header
        in_field.read((char *)&pFieldMap[i * param.matrix_length], sizeof(float) * param.matrix_length);
        in_field.close();
    }

    // read mask
    pMask = new bool[param.matrix_length];
    std::ifstream in_mask(filenames.at("mask")[0], std::ios::in | std::ios::binary);
    in_mask.seekg(sizeof(int) * 3 + sizeof(float) * 2); // skip header
    in_mask.read((char *)&pMask[0], sizeof(bool) * param.matrix_length);
    in_mask.close();

    // ========== Dump Settings ==========
    if(param.enDebug)
    {
        std::cout << "Dumping what were read:" << std::endl;
        for (int i = 0; i < param.n_fieldmaps; i++)
            std::cout << filenames.at("fieldmap")[i] << std::endl;
        
        for (int i = 0; i < param.n_sample_length; i++)
            std::cout << sample_length_all[i] << std::endl;

        param.dump();
    }

    // ========== Preparation ==========
    param.n_timepoints = param.TR / param.dt; // includes start point
    if (param.n_dummy_scan % 2 != 0)
    {
        std::cout << "Number of dummy scans must be even!";
        return 1;
    }

    // alpha/2 RF pulse (along x-axis) + TE=TR/2 relaxation
    float m0_init[3] = {0., -sinf(param.FA/2.), cosf(param.FA/2.)}; // alpha/2 RF pulse (along x-axis) + TE=TR/2 relaxation
    relax(param.e12, param.e22, m0_init);

    // ========== simulating steady-state signal ==========
    if(param.enSteadyStateSimulation)
    {
        float m0t[3], m1t[3], s_cycl;
        std::copy(m0_init, m0_init + 3, m0t);
        std::cout << "M0 after alpha/2 pulse & TR/2 relaxation = [" << m0t[0] << " " << m0t[1] << " " << m0t[2] << "]" << std::endl;
        for (int rep = 0; rep < 3; rep++)
            for (int dummy_scan = 0; dummy_scan < param.n_dummy_scan; dummy_scan++)
            {
                s_cycl = (dummy_scan % 2 == 0) ? -param.s : param.s;
                xrot(s_cycl, param.c, m0t, m1t);
                if (dummy_scan != param.n_dummy_scan - 1)
                    relax(param.e1, param.e2, m1t);
                else
                {
                    relax(param.e12, param.e22, m1t);
                    std::cout << "Magnetization after " <<param. n_dummy_scan * (rep + 1) << " RF shots = [" << m1t[0] << " " << m1t[1] << " " << m1t[2] << "]" << std::endl;
                    relax(param.e12, param.e22, m1t);
                }
                std::copy(m1t, m1t + 3, m0t);
            }
    }

    // ========== outputs ==========
    int device_count;
    checkCudaErrors(cudaGetDeviceCount(&device_count));
    std::cout << "Number of GPU(s): " << device_count << std::endl;
    
    cudaEvent_t start_0, end_0;
    checkCudaErrors(cudaEventCreate (&start_0));
    checkCudaErrors(cudaEventCreate (&end_0));

    param.n_spins /= device_count; // spins will be distributed to multiple GPUs (if there is). We hope it is divisible 
    float *p_spin_xyz = new float[3 * param.n_spins * param.n_fieldmaps * param.n_sample_length * device_count];

    checkCudaErrors(cudaEventRecord (start_0));
    #pragma omp parallel for
    for(int d=0; d<device_count; d++)
    {
        checkCudaErrors(cudaSetDevice(d));
        cudaStream_t   stream;
        cudaStreamCreate(&stream);

        simulation_parameters *d_param, param_local;
        float *d_pFieldMap, *d_random_walk_xyz_init, *d_random_walk_xyz_init_scaled;
        float*d_mxyz;
        bool *d_pMask;

        checkCudaErrors(cudaMalloc((void**)&d_param,                sizeof(simulation_parameters)));
        checkCudaErrors(cudaMalloc((void**)&d_pFieldMap,            sizeof(float) * param.matrix_length * param.n_fieldmaps));    
        checkCudaErrors(cudaMalloc((void**)&d_random_walk_xyz_init, sizeof(float) * 3 * param.n_spins));
        checkCudaErrors(cudaMalloc((void**)&d_random_walk_xyz_init_scaled, sizeof(float) * 3 * param.n_spins));
        checkCudaErrors(cudaMalloc((void**)&d_mxyz,                 sizeof(float) * 3 * param.n_spins * param.n_fieldmaps));
        checkCudaErrors(cudaMalloc((void**)&d_pMask,                sizeof(bool) * param.matrix_length));

        checkCudaErrors(cudaMemcpyAsync(d_pFieldMap, pFieldMap, sizeof(float) * param.matrix_length * param.n_fieldmaps, cudaMemcpyHostToDevice, stream)); 
        checkCudaErrors(cudaMemcpyAsync(d_pMask, pMask, sizeof(bool) * param.matrix_length, cudaMemcpyHostToDevice, stream));
        checkCudaErrors(cudaMemcpyAsync(d_param, &param, sizeof(simulation_parameters), cudaMemcpyHostToDevice, stream));
        memcpy(&param_local, &param, sizeof(simulation_parameters));

        // ========== run ==========
        int numBlocks = (param.n_spins + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        float elapsedTime;
        cudaEvent_t start;
        cudaEvent_t end;
        checkCudaErrors(cudaEventCreate (&start));
        checkCudaErrors(cudaEventCreate (&end));
        checkCudaErrors(cudaEventRecord (start));
        // generate initial spatial position for spins, based on sample_length_ref
        generate_initial_position_cuda<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_random_walk_xyz_init, d_param, d_pMask);
        gpuCheckKernelExecutionError( __FILE__, __LINE__);

        float sample_length_ref = param.sample_length;
        
        for (int sl = 0; sl < param.n_sample_length; sl++)
        {        
            float sample_length_scale = sample_length_all[sl] / sample_length_ref;
            if (param.n_sample_length > 1)
                printf("%d, %d ) Simulating sample size = %.2f um, scale to reference = %.2f\n", d, sl, sample_length_all[sl]*1e6, sample_length_scale);

            param_local.sample_length = sample_length_all[sl];    
            param_local.scale2grid = (param_local.fieldmap_size[0] - 1.) / param_local.sample_length;    
            cudaMemcpy(d_param, &param_local, sizeof(simulation_parameters), cudaMemcpyHostToDevice);

            scale_initial_positions<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_random_walk_xyz_init_scaled, d_random_walk_xyz_init, sample_length_scale, param.n_spins);
            gpuCheckKernelExecutionError( __FILE__, __LINE__);

            simulation_kernel_cuda<<<numBlocks, THREADS_PER_BLOCK, 0, stream>>>(d_param, 
                                                                                d_pFieldMap, 
                                                                                d_pMask, 
                                                                                d_random_walk_xyz_init_scaled, 
                                                                                d_mxyz);
            gpuCheckKernelExecutionError( __FILE__, __LINE__);                                                       
            checkCudaErrors(cudaMemcpyAsync(p_spin_xyz + 3*param.n_spins*param.n_fieldmaps*sl + 3*param.n_spins*param.n_fieldmaps*param.n_sample_length*d
                                            , d_mxyz, sizeof(float)*3*param.n_spins*param.n_fieldmaps, cudaMemcpyDeviceToHost, stream));
        }    

        checkCudaErrors(cudaEventRecord(end));
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start, end));
        std::cout << d << ") Simulation took " << elapsedTime/1000 << " seconds" << std::endl;

        // ========== clean up GPU ==========
        checkCudaErrors(cudaFree(d_param));
        checkCudaErrors(cudaFree(d_pFieldMap));
        checkCudaErrors(cudaFree(d_pMask));
        checkCudaErrors(cudaFree(d_random_walk_xyz_init));
        checkCudaErrors(cudaFree(d_random_walk_xyz_init_scaled));
        checkCudaErrors(cudaFree(d_mxyz));
        checkCudaErrors(cudaEventDestroy(start));
        checkCudaErrors(cudaEventDestroy(end));
        checkCudaErrors(cudaStreamDestroy(stream));
    }

    float elapsedTime;
    checkCudaErrors(cudaEventRecord(end_0));
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaEventElapsedTime(&elapsedTime, start_0, end_0));
    std::cout << "Entire simulation over " << device_count << " GPU(s) took " << elapsedTime/1000 << " seconds" << std::endl;
    checkCudaErrors(cudaEventDestroy(start_0));
    checkCudaErrors(cudaEventDestroy(end_0));
    
    // ========== save results ========== 
    output_header oh(param.n_spins, param.n_fieldmaps, param.n_sample_length, device_count);
    std::fstream out_spin_xyz(filenames.at("output")[0], std::ios::out | std::ios::binary);
    out_spin_xyz.write((char*)&oh, sizeof(output_header));
    out_spin_xyz.write((char*)&p_spin_xyz[0], sizeof(float) * 3 * param.n_spins * param.n_fieldmaps * param.n_sample_length * device_count);
    out_spin_xyz.close();

    // ========== clean up CPU ==========
    delete[] p_spin_xyz;
    delete[] pFieldMap;
    delete[] pMask;
}

__global__ void simulation_kernel_cuda(simulation_parameters *param_orig, 
                                       float *d_pFieldMap, 
                                       bool *d_pMask,
                                       float *d_random_walk_xyz_init_scaled,  // 3 * param.n_spins
                                       float *d_mxyz)
{
    int spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    __shared__ simulation_parameters param;
    if(threadIdx.x == 0)
        memcpy(&param, param_orig, sizeof(simulation_parameters));
    __syncthreads();

    // if (spin_no == 1000)
    //     param.ddump();
    if (spin_no >= param.n_spins)
        return;

    int16_t n_timepoints_local;
    int32_t seed  = param.seed + spin_no;

    float accumulated_phase[MAX_N_FIELDMAP] = {0.}; // please note array size is hard-coded! works for now, but limits flexibility.
    float field[MAX_N_FIELDMAP];
    float m0[3*MAX_N_FIELDMAP], m1[3*MAX_N_FIELDMAP]; // m0[0,1,2] * nfieldmap - please note array size is hard-coded! works for now, but limits flexibility.
    float xyz[3], xyz_new[3];
    for(int i=0; i<3; i++)
        xyz[i] = d_random_walk_xyz_init_scaled[3*spin_no + i];

    thrust::minstd_rand  gen(seed);
    thrust::normal_distribution<float> dist_random_walk_xyz(0.f, sqrt(6 * param.diffusion_const * param.dt)); // running duration: 280 sec
    gen.discard(seed);

    // alpha/2 RF pulse (along x-axis) + TR/2 relaxation
    for(int i=0; i<param.n_fieldmaps; i++)
    {
        m0[3*i+0] = 0;
        m0[3*i+1] = -param.s2 * param.e22;
        m0[3*i+2] =  1. + param.e12 * (param.c2 - 1.);
    }

    bool is_lastdummy = false;
    for (int dummy_scan = 0; dummy_scan < param.n_dummy_scan + 1; dummy_scan++)
    {
        is_lastdummy = (dummy_scan == param.n_dummy_scan);
        n_timepoints_local = is_lastdummy ? (param.n_timepoints) / 2 : (param.n_timepoints); // random walk till TR/2 in the final execution of the loop
        
        // alpha RF pulse (along x-axis) + TR or TR/2 relaxation
        float s_cycl = (dummy_scan % 2 == 0) ? -param.s : param.s; // PI phase cycling, starts with -FA (since we have +FA/2 above)
        for (int i = 0; i< param.n_fieldmaps; i++)
            xrot(s_cycl, param.c, m0 + 3*i, m1 + 3*i);

        // copy m1 to m0
        for(int i=0; i<3*param.n_fieldmaps; i++)
            m0[i] = m1[i];

        // random walk with boundries and accomulate phase
        int32_t ind=0, ind_old=-1;
        int16_t current_timepoint = 0;
        for(int i=0; i<param.n_fieldmaps; i++)
            accumulated_phase[i] = 0;
        while (current_timepoint < n_timepoints_local)
        {
            for (int i=0; i<3; i++)
            {
                xyz_new[i] = xyz[i] + dist_random_walk_xyz(gen); // new spin position after random-walk
                if (xyz_new[i] < 0)
                    xyz_new[i] += param.sample_length;
                else if (xyz_new[i] > param.sample_length)
                    xyz_new[i] -= param.sample_length;
            }
            // subscripts to linear indices
            ind = sub2ind(ROUND(xyz_new[0]*param.scale2grid), ROUND(xyz_new[1]*param.scale2grid), ROUND(xyz_new[2]*param.scale2grid), param.fieldmap_size[0], param.fieldmap_size[1]);
            //accumulate phase
            if(ind != ind_old) // used this trick to access fewer to the global memorsy which is slow. Helpful for big samples!
            {               
                if (d_pMask[ind] == true) // check doesn't cross a vessel 
                    continue;
                for (int i = 0; i<param.n_fieldmaps; i++)            
                    field[i] = d_pFieldMap[ind + i*param.matrix_length];
                ind_old = ind; 
            }

            for (int i = 0; i<param.n_fieldmaps; i++)            
                accumulated_phase[i] += field[i];

            for (int i = 0; i < 3; i++)
                xyz[i] = xyz_new[i];
 
            current_timepoint++;            
        }

        for (int i=0; i<param.n_fieldmaps; i++)
        {            
            accumulated_phase[i] *= param.B0 * GAMMA * param.dt; // Fieldmap per Tesla to radian     
            zrot(accumulated_phase[i], m0 + 3*i, m1 + 3*i); // dephase
            relax(is_lastdummy ? param.e12 : param.e1, is_lastdummy ? param.e22 : param.e2, m1 + 3*i);
        }
        // copy m1 to m0 for the next iteration
        for(int i=0; i<3*param.n_fieldmaps; i++)
            m0[i] = m1[i];
    }

    for (int i=0; i<param.n_fieldmaps; i++)
    {  
        int ind = 3*param.n_spins*i + 3*spin_no;
        d_mxyz[ind + 0] = m1[3*i + 0];
        d_mxyz[ind + 1] = m1[3*i + 1];
        d_mxyz[ind + 2] = m1[3*i + 2];
    }
}

__global__ void scale_initial_positions(float *d_scaled_xyz, float *d_initial_xyz, float scale, int size)
{
    int n = blockIdx.x * blockDim.x + threadIdx.x ;
    if(n < size)
    {
        int ind = 3*n;
        d_scaled_xyz[ind+0] = d_initial_xyz[ind+0] * scale;
        d_scaled_xyz[ind+1] = d_initial_xyz[ind+1] * scale;
        d_scaled_xyz[ind+2] = d_initial_xyz[ind+2] * scale;
    }
}

// generate random initial position
__global__ void generate_initial_position_cuda(float *d_spin_position_xyz, simulation_parameters *param, bool *pMask)
{
    // prepare random generator engine
    int spin_no = blockIdx.x * blockDim.x + threadIdx.x ;
    if(spin_no >= param->n_spins)
        return;

    thrust::minstd_rand  gen(param->seed + spin_no);
    thrust::uniform_real_distribution<float> dist_initial_point((float)param->sample_length * 0.1, (float)param->sample_length * 0.9);
    gen.discard(param->seed + spin_no);
    
    float scale2grid = (param->fieldmap_size[0]-1.) / param->sample_length;
    int32_t index = 0;
    float *spin_position_xyz = d_spin_position_xyz + 3*spin_no;
    do
    {
        for (int i = 0; i < 3; i++)
            spin_position_xyz[i] = dist_initial_point(gen);
        index = sub2ind(ROUND(spin_position_xyz[0]*scale2grid), ROUND(spin_position_xyz[1]*scale2grid), ROUND(spin_position_xyz[2]*scale2grid), param->fieldmap_size[0], param->fieldmap_size[1]);
    }while (pMask[index] == true);
}

