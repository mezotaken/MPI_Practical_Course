#include <iostream>
#include <mpi.h>
#include <assert.h>

using namespace std;

int main (int argc, char** argv)
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size = -1;
    assert(MPI_Comm_size(MPI_COMM_WORLD, &world_size) == MPI_SUCCESS);

    // Get the rank of the process
    int world_rank = -1;
    assert(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank) == MPI_SUCCESS);

    // Print off a hello world message
    cout << "Hello world from processor rank " << world_rank << " out of " << world_size << " processors\n",

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}
