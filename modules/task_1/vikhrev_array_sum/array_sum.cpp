#include <mpi.h>
#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <time.h>

#define MainProc 0
#define OutputSize 26

int main(int argc, char* argv[])
{
	// mpi variables
	int status, rank, size, nProc = 0;
	double t1 = 0, t2 = 0;
	//usual variables
	int  cols,
		 rows,
		 arrSize,
		 N; // number of elems that will be given to one process
	double *arr = NULL, //matrix
		   *buff, // buffer for messages
		    partialSum = 0,
		    arrSum = 0;

	//read argv
	if (argc > 2)
	{
		rows = atoi(argv[1]);
		cols =  atoi(argv[2]);
	}
	else
		rows = cols = 10;
	arrSize = rows * cols;

	//mpi part
	status = MPI_Init(&argc, &argv);
	if (status != MPI_SUCCESS) { return -1; }

	status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (status != MPI_SUCCESS) { return -1; }

	status = MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (status != MPI_SUCCESS) { return -1; }

	if (size > 64) return -1; // limit

	nProc = (arrSize >= size) ? size : arrSize;
	N = arrSize / nProc;

	if (rank == MainProc) // sending
	{
		//init arr
		arr = new double[arrSize];
		std::srand((unsigned)time(NULL));
		for (int i = 0; i < arrSize; i++)
			arr[i] = (std::rand()%100)/100.0;

		if (arrSize < OutputSize)
		{
			for (int i = 0; i < arrSize; i++)
			{
				std::cout << arr[i] << " ";
				if((i+1) % cols == 0)
					std::cout << std::endl;
			}
		}
		std::cout << std::endl;
		buff = arr;

		//sequential part
		t1 = MPI_Wtime();
		for (int i = 0; i < arrSize; i++)
			arrSum += arr[i];
		t2 = MPI_Wtime();
		std::cout << "Sequential Time: " << t2 - t1 << std::endl;
		std::cout << "Sequential Sum = " << arrSum << std::endl;
		//end sequential part

		t1 = MPI_Wtime(); // start
		for (int i = 1; i < nProc - 1; i++)
			status = MPI_Send(&arr[i*N], N, MPI_DOUBLE, i, i, MPI_COMM_WORLD);

		if (size != 1) // sending last part separately in case if arrSize % size != 0
		    MPI_Send(&arr[N*(nProc - 1)], arrSize - N * (nProc - 1), MPI_DOUBLE, nProc - 1, nProc - 1, MPI_COMM_WORLD);


	}
	else //receiving
	{
		if (rank == nProc - 1)
		{
			N = arrSize - N * (nProc - 1);
			buff = new double[N];
		    MPI_Recv(buff, N, MPI_DOUBLE, MainProc, rank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}
		else
		{
			buff = new double[N];
			MPI_Recv(buff, N, MPI_DOUBLE, MainProc, rank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}
	}

	//common part
	for (int i = 0; i < N; i++)
		partialSum += buff[i];

	//sum
	if (rank == MainProc)
	{
		arrSum = 0;
		arrSum += partialSum;

		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&partialSum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			arrSum += partialSum;
		}
		t2 = MPI_Wtime();
		std::cout << std::endl;
		std::cout << "Parallel Time: " << t2 - t1 << std::endl;
		std::cout <<"Parallel Sum = " << arrSum << std::endl;
	}
	else
	{
		MPI_Send(&partialSum, 1, MPI_DOUBLE, MainProc, 0, MPI_COMM_WORLD);
	}

	//del
	if (rank == 0)
		delete[] arr;
	else
		delete[] buff;

	status = MPI_Finalize();
	if (status != MPI_SUCCESS) { return -1; }
	return 0;
}
