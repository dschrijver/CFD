#include <mpi.h>
#include <stdlib.h>

#include "../include/datatypes.h"
#include "../include/utils.h"

void communicate(double *field, BufferBag *buffers, ParamBag *params)
{
    UNPACK_GRID
    UNPACK_BUFFERS

    MPI_Comm comm_cart = params->comm_cart;

    MPI_Status status_first;
 
    // Communicate to left
    for (int j = j_start; j < j_end; j++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            send_buffer_x[NZ_proc*(j - j_start) + (k - k_start)] = field[IDX(i_start, j, k)];
        }
    }
    MPI_Sendrecv(send_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_left, 0, recv_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_right, MPI_ANY_TAG, comm_cart, &status_first);
    for (int j = j_start; j < j_end; j++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            field[IDX(i_end, j, k)] = recv_buffer_x[NZ_proc*(j - j_start) + (k - k_start)];
        }
    }

    // Communicate to right
    for (int j = j_start; j < j_end; j++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            send_buffer_x[NZ_proc*(j - j_start) + (k - k_start)] = field[IDX(i_end - 1, j, k)];
        }
    }
    MPI_Sendrecv(send_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_right, 0, recv_buffer_x, buffer_number_x, MPI_DOUBLE, params->process_left, MPI_ANY_TAG, comm_cart, &status_first);
    for (int j = j_start; j < j_end; j++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            field[IDX(i_start - 1, j, k)] = recv_buffer_x[NZ_proc*(j - j_start) + (k - k_start)];
        }
    } 

    // Communicate to bottom
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            send_buffer_y[NZ_proc*(i - i_start) + (k - k_start)] = field[IDX(i, j_start, k)];
        }
    }
    MPI_Sendrecv(send_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_bottom, 0, recv_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_top, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            field[IDX(i, j_end, k)] = recv_buffer_y[NZ_proc*(i - i_start) + (k - k_start)];
        }
    }

    // Communicate to top
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            send_buffer_y[NZ_proc*(i - i_start) + (k - k_start)] = field[IDX(i, j_end - 1, k)];
        }
    }
    MPI_Sendrecv(send_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_top, 0, recv_buffer_y, buffer_number_y, MPI_DOUBLE, params->process_bottom, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start; i < i_end; i++)
    {
        for (int k = k_start; k < k_end; k++)
        {
            field[IDX(i, j_start - 1, k)] = recv_buffer_y[NZ_proc*(i - i_start) + (k - k_start)];
        }
    }

    // Communicate to back
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            send_buffer_z[NY_proc*(i - i_start) + (j - j_start)] = field[IDX(i, j, k_start)];
        }
    }
    MPI_Sendrecv(send_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_back, 0, recv_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_front, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            field[IDX(i, j, k_end)] = recv_buffer_z[NY_proc*(i - i_start) + (j - j_start)];
        }
    }

    // Communicate to front
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            send_buffer_z[NY_proc*(i - i_start) + (j - j_start)] = field[IDX(i, j, k_end - 1)];
        }
    }
    MPI_Sendrecv(send_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_front, 0, recv_buffer_z, buffer_number_z, MPI_DOUBLE, params->process_back, MPI_ANY_TAG, comm_cart, &status_first);
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            field[IDX(i, j, k_start - 1)] = recv_buffer_z[NY_proc*(i - i_start) + (j - j_start)];
        }
    }
}