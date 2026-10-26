#ifndef DATATYPES_H
#define DATATYPES_H

#include <mpi.h>
#include <hdf5.h>

typedef struct ParamBag
{
    // Simulation variables
    double time;
    double time_output;
    double time_log;
    int n_output;

    // Parameters
    int NX;   
    int NY;   
    int NZ;
    double dx;
    double t_final;
    double dt;
    double dt_output;
    double dt_log;
    double Re;
    double Fx;
    double Fy;
    double Fz;
    double poisson_epsilon;

    // MPI
    int number_of_processes; 
    int process_rank;     
    int process_coords[3]; 
    int process_left;
    int process_right;
    int process_bottom;
    int process_top;
    int process_back;
    int process_front;
    int i_start;              
    int i_end;                
    int NX_proc;         
    int j_start;           
    int j_end;              
    int NY_proc;            
    int k_start;             
    int k_end;              
    int NZ_proc;
    int LX;
    int LY;
    int LZ;            
    MPI_Comm comm_cart;

    // HDF5
    hid_t fapl_id;
    hid_t scalar_space;
    hid_t filespace;
    hid_t memspace;
    hid_t dcpl_id;
    hid_t dxpl_id;
} ParamBag;

typedef struct BufferBag
{
    int buffer_number_x;
    int buffer_number_y;
    int buffer_number_z;

    double *send_buffer_x; 
    double *send_buffer_y;
    double *send_buffer_z;
    double *recv_buffer_x;
    double *recv_buffer_y;
    double *recv_buffer_z;
} BufferBag;

#endif