#include "../include/datatypes.h"
#include "../include/utils.h"

void output(double *field, char *name, ParamBag *params)
{
    UNPACK_GRID

    char filename[32];

    double time = params->time;

    // Create file
    sprintf(filename, "%s_%d.h5", name, params->n_output);
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    H5Pclose(fapl_id);

    // Write time
    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t dset_scalar = H5Dcreate2(file_id, "time", H5T_NATIVE_DOUBLE, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_scalar, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &time);
    H5Dclose(dset_scalar);
    H5Sclose(scalar_space);

    // Space occupied in file
    hsize_t dims_file[3] = {NX, NY, NZ};
    hid_t filespace = H5Screate_simple(3, dims_file, NULL);

    // Space occupied in processor memory
    hsize_t dims_proc[3] = {NX_proc + 2, NY_proc + 2, NZ_proc + 2};
    hid_t memspace = H5Screate_simple(3, dims_proc, NULL);

    // Create dataset
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hid_t dset_id = H5Dcreate2(file_id, "data", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Pclose(dcpl_id);
    H5Sclose(filespace);

    // File hyperslab
    hsize_t start_file[3] = {i_start, j_start, k_start};
    hsize_t count[3] = {NX_proc, NY_proc, NZ_proc};
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start_file, NULL, count, NULL);

    // Process hyperslab
    hsize_t start_proc[3] = {1, 1, 1};
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start_proc, NULL, count, NULL);

    // Write data
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl_id, field);

    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Pclose(dxpl_id);

    // Close file
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
}