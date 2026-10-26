#include <stdio.h>

#include "../include/datatypes.h"
#include "../definitions.h"
#include "../include/utils.h"

int block_start(int coord, int n, int p) 
{
    return (int)(((long)coord * n) / p);
}

void dims_create(int dims[3], ParamBag *params)
{
    const int x_max = params->NX;
    const int y_max = params->NY;
    const int z_max = params->NZ;

    const int N_proc = params->number_of_processes;
    int A = 1000000000;

    for (int x = 1; x <= x_max; x++)
    {
        for (int y = 1; y <= y_max; y++)
        {
            for (int z = 1; z <= z_max; z++)
            {
                const int N_proc_local = x*y*z;

                if (N_proc_local != N_proc)
                    continue;
                
                const int x_length = (params->NX + x - 1) / x;
                const int y_length = (params->NY + y - 1) / y;
                const int z_length = (params->NZ + z - 1) / z;

                const int A_local = 2*(x_length*y_length + x_length*z_length + y_length*z_length);

                if (A_local < A)
                {
                    A = A_local;
                    dims[0] = x;
                    dims[1] = y;
                    dims[2] = z;
                }
            }
        }
    }
}

void initialize_MPI(ParamBag *params)
{
    MPI_Comm_size(MPI_COMM_WORLD, &params->number_of_processes);

    int dims[3] = {0};

    dims_create(dims, params);

    MPI_Comm_rank(MPI_COMM_WORLD, &params->process_rank);

    if (dims[0] == 0)
    {
        if (params->process_rank == 0)
        {
            printf("Couldn't find decomposition!\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (params->process_rank == 0)
    {
        printf("Domain decomposition: [%d][%d][%d]\n", dims[0], dims[1], dims[2]);
    }

    int xperiodic = 0;
    int yperiodic = 0;
    int zperiodic = 0;
#ifdef XPERIODIC
    xperiodic = 1;
#endif
#ifdef YPERIODIC
    yperiodic = 1;
#endif
#ifdef ZPERIODIC
    zperiodic = 1;
#endif
    int periods[3] = {xperiodic, yperiodic, zperiodic};
    int reorder = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &params->comm_cart);
    MPI_Comm_rank(params->comm_cart, &params->process_rank);
    MPI_Cart_get(params->comm_cart, 3, dims, periods, params->process_coords);
    MPI_Cart_shift(params->comm_cart, 0, 1, &params->process_left, &params->process_right);
    MPI_Cart_shift(params->comm_cart, 1, 1, &params->process_bottom, &params->process_top);
    MPI_Cart_shift(params->comm_cart, 2, 1, &params->process_back, &params->process_front);

    params->i_start = block_start(params->process_coords[0], params->NX, dims[0]);
    params->i_end = block_start((params->process_coords[0] + 1), params->NX, dims[0]);

    params->j_start = block_start(params->process_coords[1], params->NY, dims[1]);
    params->j_end = block_start((params->process_coords[1] + 1), params->NY, dims[1]);

    params->k_start = block_start(params->process_coords[2], params->NZ, dims[2]);
    params->k_end = block_start((params->process_coords[2] + 1), params->NZ, dims[2]);

    params->NX_proc = params->i_end - params->i_start;
    params->NY_proc = params->j_end - params->j_start;
    params->NZ_proc = params->k_end - params->k_start;

    params->LX = params->NX * params->dx;
    params->LY = params->NY * params->dx;
    params->LZ = params->NZ * params->dx;
}

void initialize_HDF5(ParamBag *params)
{
    params->fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(params->fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    params->scalar_space = H5Screate(H5S_SCALAR);

    int NX = params->NX;
    int NY = params->NY;
    int NZ = params->NZ;
    int NX_proc = params->NX_proc;
    int NY_proc = params->NY_proc;
    int NZ_proc = params->NZ_proc;

    // Space occupied in file
    hsize_t dims_file[3] = {NX, NY, NZ};
    params->filespace = H5Screate_simple(3, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc[3] = {NX_proc + 4, NY_proc + 4, NZ_proc + 4};
    params->memspace = H5Screate_simple(3, dims_proc, NULL);

    params->dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    params->dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(params->dxpl_id, H5FD_MPIO_COLLECTIVE);
}

void initialize_fields(double *u, double *v, double *w, double *pressure, ParamBag *params)
{
    UNPACK_GRID

    FOR_DOMAIN
    {
        UNPACK_POSITIONS

        const int idx0 = IDX(i, j, k);

#ifdef INI_STATIONARY
        u[idx0] = 0.0;
        v[idx0] = 0.0;
        w[idx0] = 0.0;
#endif

#ifdef INI_PARABOLIC
        u[idx0] = 4.0 * INI_PARABOLIC * y_center * (LY - y_center);
        v[idx0] = 0.0;
        w[idx0] = 0.0;
#endif

        pressure[idx0] = 0.0;
    }
}