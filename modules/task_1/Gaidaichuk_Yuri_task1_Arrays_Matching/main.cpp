#include <mpi.h>
#include <iostream>
#include <assert.h>
#include<cmath>
#include<cstdlib>

int computeForOneRrocessor(int x){
    int *arr1 = new int[x];
    int *arr2 = new int[x];

    std::srand ( x );

    for(int i = 0; i < x; i++){
        arr1[i] = std::rand() % 10;
    }
    for(int i = 0; i < x; i++){
        arr2[i] = std::rand() % 10;
    }
    //print
    if(x <= 10){
        std::cout << "Arr1 = ";
        for(int i = 0; i < x; i++){
            std::cout << arr1[i] <<" ";
        }
        std::cout<< '\n';
        std::cout << "Arr2 = ";
        for(int i = 0; i < x; i++){
            std::cout << arr2[i] <<" ";
        }
        std::cout<< '\n';
    }
    int matchCount = 0;
    for(int i = 0; i < x; i++){
        if(arr1[i] == arr2[i]){
            matchCount++;
        }
    }
    delete[]arr1;
    delete[]arr2;
    return matchCount;
}

int main (int argc, char* argv[])

{
    int status = 0, rank = 0, size = 0;
    int x = atoi(argv[1]);
    int matchCount = 0;

    double starttime1, endtime1, starttime2, endtime2;

    status = MPI_Init (&argc, &argv);//assert(status == MPI_SUCCESS);
    if (status != MPI_SUCCESS) {
        std::cout<<"ERROR\n";
        return 1;
    }

    status = MPI_Comm_rank (MPI_COMM_WORLD, &rank);//assert(status == MPI_SUCCESS);
    if (status != MPI_SUCCESS) {
        std::cout<<"ERROR\n";
        return 1;
    }

    status = MPI_Comm_size (MPI_COMM_WORLD, &size);//assert(status == MPI_SUCCESS);
    if (status != MPI_SUCCESS) {
        std::cout<<"ERROR\n";
        return 1;
    }

    if(size == 1){
        std::cout<<"One processor result = " << computeForOneRrocessor(x) << '\n';
        status = MPI_Finalize();
        assert(status == MPI_SUCCESS);

        return 0;
    }

    if(rank == 0)
    {
            int *arr1 = new int[x * size];
            int *arr2 = new int[x * size];
            std::srand ( x );
            for(int i = 0; i < x * size; i++){
                arr1[i] = std::rand() % 10;
            }
            for(int i = 0; i < x * size; i++){
                arr2[i] = std::rand() % 10;
            }

            ///////////////////////////////////////count1
            starttime1 = MPI_Wtime();
            for(int i = 0; i < x * size; i++){
                if(arr1[i] == arr2[i]){
                    matchCount++;
                }
            }

            endtime1   = MPI_Wtime();
            ///////////////////////////////////////

            //print
            if(x <= 10){
                std::cout << "Arr1 = ";
                for(int i = 0; i < x * size; i++){
                    std::cout << arr1[i] <<" ";
                }
                std::cout<< '\n';
                std::cout << "Arr2 = ";
                for(int i = 0; i < x * size; i++){
                    std::cout << arr2[i] <<" ";
                }
                std::cout<< '\n';
                //std::cout<< "Count = " << matchCount <<'\n';
            }

            ////////////////////////////////count2
            starttime2 = MPI_Wtime();

            for(int i = 1; i <= size - 1; i++){
                MPI_Send(arr1 + x * i, x, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(arr2 + x * i, x, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

            int parallelMatchCount = 0;


            for(int i = 0; i < x; i ++){
                if(arr1[i] == arr2[i])
                {
                    parallelMatchCount++;
                }
            }
            for(int i = 1; i <= size - 1; i++){
                int receivedPartialCount = 0;
                MPI_Recv(&receivedPartialCount, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                parallelMatchCount += receivedPartialCount;
                //std::cout << i <<" returned "<< receivedPartialCount << '\n';
            }
            endtime2 = MPI_Wtime();
            ////////////////////////////////

            if(parallelMatchCount == matchCount){
                std::cout << "Success. Count = "<< matchCount <<
                "\n Sequense returns " << matchCount <<"\n" <<
                "\n Parallel returns " << parallelMatchCount << " \n" <<
                "\n Sequense takes " << endtime1 - starttime1 << " seconds \n" <<
                "\n Parallel takes " << endtime2 - starttime2 << " seconds \n" <<
                '\n';
            }
            else{
                std::cout << "Smth has gone wrong"<< '\n';
            }
            delete[]arr1;
            delete[]arr2;
    }
    if(rank != 0)
    {
            int *partialArr1 = new int[x];
            int *partialArr2 = new int[x];

            int partialMatchCount = 0;

            MPI_Recv(partialArr1, x, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            MPI_Recv(partialArr2, x, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

            for(int i = 0; i < x; i ++){
                if(partialArr1[i] == partialArr2[i])
                {
                    partialMatchCount++;
                }
                //std::cout << "Arr["<<i<<"] = "<< partialArr[i] << '\n';
            }
            MPI_Send(&partialMatchCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    }
    /*
    std::cout << "Process #" << rank << '\n';
    std::cout << "Count process: " << size << '\n';
    std::cout << "Param: " << x << '\n';
    */

    status = MPI_Finalize();
    assert(status == MPI_SUCCESS);

    return 0;
}
