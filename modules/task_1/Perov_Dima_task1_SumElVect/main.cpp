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

    double* myVector;
    int size;

    string sizeStr;
    sizeStr = argv[1];
    size = atoi(sizeStr.c_str());

    myVector = new double[size];

    // for parallel block
    double myResult = 0;
    double startTime = 0;
    double endTime;
    double myTime = 0;
    double mySum = 0;
    int remainderDiv = 0;
    int flag;
    int myId, numProcs;

    // for line block
    double myResult_ = 0;
    double startTime_ = 0;
    double endTime_ = 0;
    double myTime_ = 0;

    // Initialize the MPI environment

    MPI_Init(&argc, &argv);
    MPI_Initialized(&flag); // check
    if (!flag) {
        std::cout << "Error: MPI_Init";
    }

    // Description of the communicator
    // communicator manages groups of parallel processes
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs); // determining the number of processes in a group
    MPI_Comm_rank(MPI_COMM_WORLD, &myId); // determining the rank of a process in a group

	double* myBlock;
	myBlock = new double[size / numProcs];

    if (myId == 0) {

        for (int i = 0; i < size; i++) {
            myVector[i] = (double)(rand()%5) / ((double)(rand()%10) + 1);
        }
        //*LINE BLOCK*
        startTime_ = MPI_Wtime();
        for (int i = 0; i < size; i++) {
            myResult_ += myVector[i];
        }
        endTime_ = MPI_Wtime();
        cout << endl << "*LINE*" << endl;
        cout << "Result: " << myResult_ << endl;
        myTime_ = endTime_ - startTime_;
        cout << "Time: " << myTime_ << " sec" << endl;
        //*END

        //*PARALLEL BLOCK*
    startTime = MPI_Wtime(); // started counting time in 0 proccess

    }

    //MPI_Bcast(myVector, size , MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(myVector, size / numProcs, MPI_DOUBLE, myBlock, size / numProcs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = 0; i < size / numProcs; i++) {
		mySum += myBlock[i];
	}

	/*
    for (int i = size / numProcs * myId; i < size / numProcs + size / numProcs * myId; i++) {
        mySum += myVector[i];
    }
	*/

    MPI_Reduce(&mySum, &myResult, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myId == 0) {
         remainderDiv = size % numProcs;
         if (remainderDiv) {
            for (int i = 0; i < remainderDiv; i++)
                myResult += myVector[size - remainderDiv  + i];
         }

        //Parallel Results
        endTime = MPI_Wtime();
        cout << endl <<"*PARALLEL*" << endl;
        cout << "Result: " << myResult << endl;
        myTime = endTime - startTime;
        cout << "Time: " << myTime << " sec" << endl << endl;
    }


    MPI_Finalize();
    //*END PARALLEL BLOCK*

    delete[]myVector;

    return 0;
}
