#include <iostream>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <assert.h>

using namespace std;

// minimum search function 
    double minimum (double* vector, int size) {
    	double min = vector[0];
    	for (int i = 1; i < size; i++) {
    	if (vector[i] < min)
        	min = vector[i];
    	}
    	return min;
	}

int main(int argc, char *argv[])
{
	// initialize srand
    srand((int)time(0));
	
	// global variables
    double* Vector;
    double* MyVector;
    double* min_vector;
    int size;
    MPI_Status status;
	
	// initialize argument from cmd
    if (argc>1) {
        size = atoi(argv[1]);
    }
    else
        size = 100;

    Vector = new double[size]; // create the vector

    // fill the vector
    for (int i = 0; i < size; i++) {
        Vector[i] = 20+(((double)(rand()%750))/ ((double)(rand()%143) + 1));
    }

    // variables for parallel block
    double Min = 0;
    double startTime = 0;
    double endTime = 0;
    double Time = 0;
    double myMin = 0;
    int remainder = 0;
    int flag;
    int Id, numProcs;

    // variables for line block
    double Min_l = 0;
    double startTime_l = 0;
    double endTime_l = 0;
    double Time_l = 0;

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Initialized(&flag); // check
    if (!flag) {
        cout << "Error in MPI_Init!!!";
    }

    // communicators
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &Id);

    if (Id == 0) {	
    // LINE BLOCK
        startTime_l = MPI_Wtime();
        Min_l = minimum(Vector, size);
    	// Line Results
        endTime_l = MPI_Wtime();
        cout << endl << "===LINE===" << endl;
        cout << "Min: " << Min_l << endl;
        Time_l = endTime_l - startTime_l;
        cout << "Time: " << Time_l << " sec" << endl;
    // END LINE BLOCK

    // PARALLEL BLOCK
    	startTime = MPI_Wtime(); // started counting time in 0 proccess
    	min_vector = new double[numProcs];
    	for (int i = 1; i < numProcs; i++)
    		MPI_Send(Vector+size/numProcs*i, size/numProcs, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
   	min_vector[0] = minimum(Vector,size/numProcs);
   	for (int i = 1; i < numProcs; i++)
            MPI_Recv(&min_vector[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
        // check the remainder of vector
        remainder = size % numProcs;
        if (remainder) {
            for (int i = 0; i < remainder; i++)
                if (Min > Vector[size - remainder  + i])
                    Min = Vector[size - remainder  + i];
            if (Min < min_vector[0])
            	min_vector[0] = Min;
    	}
        // Parallel Results
        Min = minimum(min_vector, numProcs);
        delete[]min_vector;
        endTime = MPI_Wtime();
        cout << endl <<"===PARALLEL===" << endl;
        cout << "Min: " << Min << endl;
        Time = endTime - startTime;
        cout << "Time: " << Time << " sec" << endl;
    }
    else {
    	MyVector = new double[size/numProcs];
        MPI_Recv(MyVector, size/numProcs, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        myMin = minimum(MyVector, size/numProcs);
        delete[]MyVector;
        MPI_Send(&myMin, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    // END PARALLEL BLOCK
    delete[]Vector; // delete the vector
    return 0;
}
