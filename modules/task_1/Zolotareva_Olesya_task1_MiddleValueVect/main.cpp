#include <math.h> 
#include <iostream>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <assert.h>
using namespace std;

int main(int argc, char *argv[])
{

    srand((int)time(0));

    int* x;
    int N;

    if (argc>1) {
        N= atoi(argv[1]);
    }
    else
        N = 50;


    int ProcRank, ProcNum,flag;

    MPI_Status Status;

	x = new int[N];

	for (int i = 0; i < N; i++) x[i] = rand() % 2000 + 1;
	if (N < 50)
	{
		cout << "x:  ";
		for (int i = 0; i < N; i++) cout << x[i] << " ";
	}

    MPI_Init(&argc,&argv);
	MPI_Initialized(&flag);

	if (!flag)
	{
		cout << "Error";
		delete[] x;
		return 0;
	}

    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&ProcRank);
		int ParallelSum = 0;
	//double tempsum = 0;
		double ParallelTimeStart=0;
		double ParallelTimeEnd = 0;
		int tempSum = 0;
    
	if (ProcRank == 0)
	{
		//serial v

		int SerialSum = 0;
		double SerialTimeStart = MPI_Wtime();

		for (int i = 0; i < N; i++)
		{
			SerialSum += x[i];
		}
		double SerialMidscore = (double)SerialSum / N;
		double SerialTimeEnd = MPI_Wtime();

		//double SerialTime = MPI_Wtime() - SerialTimeStart;

		cout << "Midscore (line) = " << SerialMidscore<<endl; 
		cout << "Time (line) = " << SerialTimeEnd - SerialTimeStart<<endl;

		//parallel v
		if (ProcNum > 1)
		{
			ParallelTimeStart = MPI_Wtime();
			int p = N / ProcNum ;
			for (int i = 1; i < ProcNum - 1; i++)
			{

				// MPI_Send(&blockSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send((x + i * p), p, MPI_INT, i, 0, MPI_COMM_WORLD);
			}

			// int blockSize = N - p*(ProcNum - 1);
			 //MPI_Send(&blockSize, 1, MPI_INT, ProcNum-1, 0, MPI_COMM_WORLD);
			MPI_Send((x + (ProcNum - 1)*p),( N- (ProcNum - 1)*p), MPI_INT, (ProcNum - 1), 0, MPI_COMM_WORLD);
			for (int i = 0; i < p; i++)
			{
				ParallelSum += x[i];
			}

			for (int i = 1; i < ProcNum; i++)
			
			{
				MPI_Recv(&tempSum, 1, MPI_INT, i , 0, MPI_COMM_WORLD, &Status);
				ParallelSum += tempSum;
			}

			ParallelTimeEnd = MPI_Wtime();
			cout << "\nMidScore(parallel)=" << (double)ParallelSum / N<<endl;
			cout << "\n Time(parallel)=" << (ParallelTimeEnd - ParallelTimeStart)<<endl;
		}
	}
    else
	{
		if (ProcRank != ProcNum-1) {


			int blockSize = N / ProcNum ;
			int *temp;
			temp = new int[blockSize];
			//MPI_Recv(&blockSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
			MPI_Recv(temp, blockSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
			for (int i = 0; i < blockSize; i++)
			{
				tempSum += temp[i];
			}
			MPI_Send(&tempSum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			//delete []temp;
		}
		else {
		int blockSize = N - (ProcNum - 1)*(N / ProcNum);
		int *temp;
		temp = new int[blockSize];
		MPI_Recv(temp, blockSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status);
		for (int i = 0; i < blockSize; i++)
			tempSum += temp[i];
		MPI_Send(&tempSum,1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}}
    MPI_Finalize();

    delete []x;
}
