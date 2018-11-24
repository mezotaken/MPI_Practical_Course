#include <iostream>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <assert.h>


using namespace std;

int main(int argc, char* argv[])
{
	srand((int)time(NULL));

	int rank, CountP, flag;
	if (argc < 2)
	{
		cout << "Error";
		return -1;
	}

	int sizev = atoi(argv[1]);
	if (sizev <= 0)
	{
		cout << "Error";
		return -1;
	}


	double *v1 = NULL;
	double *v2 = NULL;

	//MPI_Status status; 
	double Line_MultVec = 0;
	double MultVec = 0;
	// double tempSum = 0;
	double tmpMult = 0;
	double Time_begin = 0;
	double Time_end = 0;
	int i = 0;

	MPI_Init(&argc, &argv);
	MPI_Initialized(&flag);

	if (!flag)
	{
		cout << "Error";
		return -1;
	}
	MPI_Comm_size(MPI_COMM_WORLD, &CountP);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double *tmpv1 = new double[sizev / CountP];
	double *tmpv2 = new double[sizev / CountP];

	if (rank == 0)
	{
		v1 = new double[sizev];
		v2 = new double[sizev];

		for (i = 0; i < sizev; i++)
		{
			v1[i] = rand() % 10 - 5;
			v2[i] = rand() % 10 - 5;
		}

		if (sizev < 21)
		{
			cout << "v1:";
			for (i = 0; i < sizev; i++)
			{
				cout << " " << v1[i];
			}
			cout << "\nv2:";
			for (i = 0; i < sizev; i++)
			{
				cout << " " << v2[i];
			}
			cout << endl;
		}

		//  --------------  Line  --------------------
		Time_begin = MPI_Wtime();
		for (int j = 0; j < sizev; j++)
			Line_MultVec += v1[j] * v2[j];
		Time_end = MPI_Wtime();
		cout << "Line_Result = " << Line_MultVec << endl;
		cout << "Line_Time = " << Time_end - Time_begin << endl;

		// --------------  Parrallel  -------------

		if (CountP > 1)
		{
			//int offset1 = sizev/CountP;
			//int offset2 = (CountP-1)*(sizev/CountP);

			Time_begin = MPI_Wtime();
		}

	}
	/*else
		if (rank != CountP-1)
		{
			v1 = new double[sizev/(CountP)];
			v2 = new double[sizev/(CountP)];
			//MPI_Recv(v1,sizev/(CountP),MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
			//MPI_Recv(v2,sizev/(CountP),MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
			for (int i = 0; i < sizev/(CountP); i++)
				tempSum += v1[i]*v2[i];
			//MPI_Send(&tempSum, 1, MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		} */


	MPI_Scatter(v1, sizev / CountP, MPI_DOUBLE, tmpv1, sizev / CountP, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(v2, sizev / CountP, MPI_DOUBLE, tmpv2, sizev / CountP, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (i = 0; i < sizev / CountP; i++)
		tmpMult += tmpv1[i] * tmpv2[i];

	MPI_Reduce(&tmpMult, &MultVec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		for (i = sizev - sizev % CountP; i < sizev; i++)
		{
			MultVec += v1[i] * v2[i];
		}
		Time_end = MPI_Wtime();

		cout << "\nParallel_Result = " << MultVec << endl;
		cout << "Parallel_Time = " << Time_end - Time_begin << endl;
	}





	/*MPI_Send((v1+offset2), sizev - offset2, MPI_DOUBLE, (CountP-1), 0, MPI_COMM_WORLD);
	MPI_Send((v2+offset2), sizev - offset2, MPI_DOUBLE, (CountP-1), 0, MPI_COMM_WORLD);
	for (i = 1; i < CountP-1; i++)
	{
		MPI_Send(v1+i*(offset1), offset1, MPI_DOUBLE,i,0,MPI_COMM_WORLD);
		MPI_Send(v2+i*(offset1), offset1, MPI_DOUBLE,i,0,MPI_COMM_WORLD);
	}
	//MPI_Send((v1+offset2), sizev - offset2, MPI_DOUBLE, (CountP-1), 0, MPI_COMM_WORLD);
	//MPI_Send((v2+offset2), sizev - offset2, MPI_DOUBLE, (CountP-1), 0, MPI_COMM_WORLD);

	for (i = 0; i < (offset1); i++)
		MultVec += v1[i]*v2[i];

	for (i = 1; i < CountP; i++)
	{
		MPI_Recv(&tempSum,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
		MultVec += tempSum;
	}
	//MPI_Reduce(&tempSum, &MultVec, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



	//for (i = 0; i < (sizev/CountP); i++)
   //     MultVec += v1[i]*v2[i];
	Time_end = MPI_Wtime();

	cout << "\nParallel_Result = " << MultVec << endl;
	cout << "Parallel_Time = " << Time_end - Time_begin << endl;
*/
	MPI_Finalize();

	delete[]v1;
	delete[]v2;
	delete[]tmpv1;
	delete[]tmpv2;

	return 0;
}
