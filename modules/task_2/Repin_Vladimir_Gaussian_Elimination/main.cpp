// Copyright mezotaken
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>


int main(int argc, char* argv[]) {
  int status,        // MPI functions return status
    procId,        // Current process ID
    nProc;        // Number of processes

  double* matr = NULL;  // matrix n * n+1 for parallel algorithm
  double* sqmatr = NULL;  // matrix n * n+1 for sequential algorithm
  double* matrcopy = NULL;
  int mSize = 0;      // Size of processed matrix
  double wtime = 0;    // Work Time
  double div = 0;      // Multiplier for rows
  double* buf = NULL;     // Buffer for receiving rows
  double* row = NULL;    // Buffer for pivot row
  double* x = NULL;    // For parallel solution
  double* partx = NULL;  // parts of parallel solution
  int tasksize = 0;    // Aux for calc counts
  int shift = 0;      // Aux for calc displs
  int* counts = NULL;    // Number of rows for each process
  int* displs = NULL;    // Displacement of rows for each process
  int control_proc = 0;  // ID of controlling process

              // MPI initializing
  status = MPI_Init(&argc, &argv);
  if (status != MPI_SUCCESS) {
    std::cout << "Error while MPI Initializing";
    return -1;
  }
  // Getting process ID
  status = MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  if (status != MPI_SUCCESS) {
    std::cout << "Error while getting process ID";
    return -1;
  }
  // Getting number of processes
  status = MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  if (status != MPI_SUCCESS) {
    std::cout << "Error while getting number of processes";
    return -1;
  }

  // Calculating number of rows given to each process
  // and allocating auxilary data for parallel calculating
  mSize = atoi(argv[1]);
  tasksize = mSize / nProc;

  counts = new int[nProc];
  displs = new int[nProc];
  row = new double[mSize + 1];

  for (int i = 0; i < nProc; i++) {
    counts[i] = tasksize;
    if (i < mSize%nProc)
      counts[i]++;
    counts[i] *= (mSize + 1);
    displs[i] = shift;
    shift += counts[i];
  }

  buf = new double[counts[procId]];
  partx = new double[counts[procId] / (mSize + 1)];

  if (procId == 0) {
    x = new double[mSize];
    srand((unsigned int)time(NULL));
    // Initializing matrixes
    matrcopy = new double[mSize*(mSize + 1)];
    matr = new double[mSize*(mSize + 1)];
    sqmatr = new double[mSize*(mSize + 1)];
    for (int i = 0; i < mSize; i++)
      for (int j = 0; j < mSize + 1; j++)
        matrcopy[i*(mSize + 1) + j] = sqmatr[i*(mSize + 1) + j] =
        matr[i*(mSize + 1) + j] = (std::rand() % 20000) / 100.0 - 100.0;

    // If matrix is small enough, then print it
    if (mSize < 11) {
      std::cout << "Matrix" << std::endl;
      for (int i = 0; i < mSize; i++) {
        for (int j = 0; j < mSize + 1; j++)
          std::cout << matr[i*(mSize + 1) + j] << " ";
        std::cout << std::endl;
      }
    }

    // Sequential part
    wtime = MPI_Wtime();

    for (int i = 0; i < mSize; i++) {
      // Divide pivot row by pivot value
      div = sqmatr[i*(mSize + 1) + i];
      sqmatr[i*(mSize + 1) + i] = 1;
      for (int j = i + 1; j < mSize + 1; j++)
        sqmatr[i*(mSize + 1) + j] /= div;

      // Substract pivot row from every other row
      for (int j = 0; j < mSize; j++)
        if (j != i) {
          div = sqmatr[j*(mSize + 1) + i];
          for (int k = i; k < mSize + 1; k++)
            sqmatr[j*(mSize + 1) + k] -= sqmatr[i*(mSize + 1) + k] * div;
        }
    }

    wtime = MPI_Wtime() - wtime;
    // End of sequential part

    // Output of sequential results
    std::cout << "-------------------------------" << std::endl;
    std::cout << "Sequential absolute error: " << std::endl;
    double error = 0;
    double part = 0;
    for (int i = 0; i < mSize; i++) {
      part = 0;
      for (int j = 0; j < mSize; j++)
        part += matrcopy[i*(mSize + 1) + j] * sqmatr[j*(mSize + 1) + mSize];
      error += (part - matrcopy[i*(mSize + 1) + mSize])*
               (part - matrcopy[i*(mSize + 1) + mSize]);
    }
    std::cout << sqrt(error);
    std::cout << std::endl;
    std::cout << "Sequential time = " << wtime << std::endl << std::endl;
    // Parallel part
    wtime = MPI_Wtime();
  }

  // Scatter required rows between processes
  MPI_Scatterv(matr, counts, displs, MPI_DOUBLE, buf,
      counts[procId], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Aux data correction after sending
  for (int i = 0; i < nProc; i++) {
    displs[i] /= (mSize + 1);
    counts[i] /= (mSize + 1);
  }

  for (int i = 0; i < mSize; i++) {
    // Choosing control process depending on current row
    if (i == displs[control_proc + 1])
      control_proc++;

    // Divide pivot row by pivot value
    if (procId == control_proc) {
      div = buf[(i - displs[control_proc])*(mSize + 1) + i];
      for (int j = i; j < mSize + 1; j++)
        row[j - i] = buf[(i - displs[control_proc])*(mSize + 1) + j] /= div;
    }

    // Broadcasting pivot row to all processes
    MPI_Bcast(row, mSize - i + 1, MPI_DOUBLE, control_proc, MPI_COMM_WORLD);

    // Substract pivot row from every other row
    for (int j = 0; j < counts[procId]; j++)
      if (j + displs[procId] != i) {
        div = buf[j*(mSize + 1) + i];
        for (int k = i; k < mSize + 1; k++)
          buf[j*(mSize + 1) + k] -= row[k - i] * div;
      }
  }

  // Preparing parts of solution for send
  for (int j = 0; j < counts[procId]; j++)
    partx[j] = buf[j*(mSize + 1) + mSize];

  // Gather parts of solution in 0 process
  MPI_Gatherv(partx, counts[procId], MPI_DOUBLE, x, counts,
                    displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  // Output parallel results in process 0
  if (procId == 0) {
    wtime = MPI_Wtime() - wtime;
    double error = 0;
    double part = 0;
    for (int i = 0; i < mSize; i++) {
      part = 0;
      for (int j = 0; j < mSize; j++)
        part += matrcopy[i*(mSize + 1) + j] * x[j];
      error += (part - matrcopy[i*(mSize + 1) + mSize])*
               (part - matrcopy[i*(mSize + 1) + mSize]);
    }
    std::cout << sqrt(error);
    std::cout << std::endl;
    std::cout << "Parallel time: " << wtime << std::endl;
    std::cout << "-------------------------------" << std::endl;
  }

  // Free allocated memory
  if (procId == 0) {
    delete[] sqmatr;
    delete[] matr;
    delete[] matrcopy;
    delete[] x;
  }
  delete[] buf;
  delete[] counts;
  delete[] displs;
  delete[] row;
  delete[] partx;

  status = MPI_Finalize();
  return 0;
}
