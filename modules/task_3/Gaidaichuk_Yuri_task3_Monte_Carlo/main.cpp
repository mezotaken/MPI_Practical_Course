// copyright            : (C) 2018 by Yury
#include <assert.h>
#include <mpi.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>

#define leftBoard 0
#define rightBoard 1

double givenFunction(double a1) {
  return a1*a1;
}

double myRand(double fMin, double fMax) {
    double f = static_cast<double>rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double sequenceMonteCarloSum(int setPower) {
  int seed = static_cast<int>(MPI_Wtime());
  std::srand(seed);
  double sum = 0;
  for (int i = 0; i < setPower; i++) {
    double randomArg = myRand(leftBoard, rightBoard);
    double part = givenFunction(randomArg);
    sum += part;
  }
  return sum;
}

double partialMonteCarloSum(int partialSetPower, int rank) {
  // initialize specific seed
  int seed = static_cast<int>(MPI_Wtime()) * rank;
  std::srand(seed);
  double partialSum = 0;
  for (int i = 0; i < partialSetPower; i++) {
    double randomArg = myRand(leftBoard, rightBoard);
    partialSum += givenFunction(randomArg);
  }
  return partialSum;
}

int main(int argc, char * argv[]) {
  int status = 0, rank = 0, procNum = 0;
  int setPower = atoi(argv[1]);
  double localMonteCarloIntegral = 0;
  double globalMonteCarloIntegral = 0;
  double parTime1 = 0, parTime2 = 0, seqTime1 = 0, seqTime2 = 0;

  status = MPI_Init(& argc, & argv);
  // assert(status == MPI_SUCCESS);
  if (status != MPI_SUCCESS) {
    std::cout << "ERROR\n";
    return 1;
  }

  status = MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  // assert(status == MPI_SUCCESS);
  if (status != MPI_SUCCESS) {
    std::cout << "ERROR\n";
    return 1;
  }

  status = MPI_Comm_size(MPI_COMM_WORLD, & procNum);
  // assert(status == MPI_SUCCESS);
  if (status != MPI_SUCCESS) {
    std::cout << "ERROR\n";
    return 1;
  }
  // Parallel version
  // чтобы были разные наборы
  MPI_Barrier(MPI_COMM_WORLD);
  parTime1 = MPI_Wtime();
  localMonteCarloIntegral = partialMonteCarloSum(setPower / procNum, rank);
  MPI_Reduce(&localMonteCarloIntegral, &globalMonteCarloIntegral, 1,
    MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  parTime2 = MPI_Wtime();
  if (rank == 0) {
    std::cout << "Parallel result = " <<
      globalMonteCarloIntegral / setPower << "\n";
    std::cout << "Time consumed = "<< parTime2 - parTime1 << " sec\n";
    seqTime1 = MPI_Wtime();
    std::cout << "Sequence result = " <<
      sequenceMonteCarloSum(setPower) / setPower << "\n";
    seqTime2 = MPI_Wtime();
    std::cout << "Time consumed = "<< seqTime2 - seqTime1 << " sec\n";
  }
  status = MPI_Finalize();
  assert(status == MPI_SUCCESS);
  return 0;
}
