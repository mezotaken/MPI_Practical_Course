#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main (int argc, char* argv[])
{
    int status,        //MPI functions return status
        procId,        //Current process ID
        nProc;        //Number of processes
    int partres = 0,//Part of result
        result = 0; //Number of wrong-ordered pairs of elements in array
    int tasksize;    //Number of pairs to check for each process
    double* vect = NULL;    //Processed vector
    int vSize;        //Size of processed vector
    double wtime = 0;    //Work Time
    double* buf = NULL;    //Buffer for message receiving

    //MPI initializing
    status = MPI_Init (&argc, &argv);
    if (status != MPI_SUCCESS) { cout << "Error while MPI Initializing";  return -1; }
    //Getting process ID
    status = MPI_Comm_rank (MPI_COMM_WORLD, &procId);
    if (status != MPI_SUCCESS) { cout << "Error while getting process ID";  return -1; }
    //Getting number of processes
    status = MPI_Comm_size (MPI_COMM_WORLD, &nProc);
    if (status != MPI_SUCCESS) { cout << "Error while getting number of processes";  return -1; }

    //Calculating size of pairs given to each process
    if (argc < 2) { cout << "Too few arguments";  return -1; }
    vSize = atoi(argv[1]);
    tasksize = (vSize - 1) / nProc;


    if (procId == 0)
    {

        srand((unsigned int)time(NULL));

        //Initializing vector
        vect = new double[vSize];
        for (int i = 0; i < vSize; i++)
            vect[i] = (double)(rand() % 20000) / 100.0 - 100.0;
        //Sequential part
        wtime = MPI_Wtime();
        for (int i = 0; i < vSize - 1; i++)
            if (vect[i] > vect[i + 1])
                result++;
        wtime = MPI_Wtime() - wtime;
        //End of sequential part
        cout << "-------------------------------" << endl;
        cout << "Sequential result = " << result << endl;
        cout << "Sequential time = " << wtime << endl << endl;
        result = 0;



        //Parallel part
        wtime = MPI_Wtime();

        //Sending messages with required part of array
        for (int i = 0; i < nProc - 1; i++)
            MPI_Send(&vect[i*tasksize], tasksize + 1, MPI_DOUBLE, i + 1, i + 1, MPI_COMM_WORLD);


        //Processing last part with 0 process
        for (int i = tasksize*(nProc - 1); i < vSize - 1; i++)
            if (vect[i] > vect[i + 1])
                partres++;
    }
    else
    {
        buf = new double[tasksize+1];
        MPI_Recv(buf,tasksize+1,MPI_DOUBLE,0,procId,MPI_COMM_WORLD,MPI_STATUSES_IGNORE);
        for (int i = 0; i < tasksize; i++)
            if (buf[i] > buf[i + 1])
                partres++;
    }

    //Collecting partial results into process 0
    status = MPI_Reduce(&partres, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) { cout << "Error while reducing";  return -1; }
    //End of parallel part


    //Output parallel results in process 0
    if (procId == 0)
    {
        wtime = MPI_Wtime() - wtime;
        cout << "Parallel result: " << result << endl;
        cout << "Parallel time: " << wtime << endl;
        cout << "-------------------------------" << endl;
    }

    if (procId == 0)
        delete[] vect;
    else
        delete[] buf;

    status = MPI_Finalize();
    if (status != MPI_SUCCESS)
    { cout << "Error while MPI Finilizing";  return -1; }

    return 0;
}
