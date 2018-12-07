// Copyright ivanvikhrev
#include <assert.h>
#include <iostream>
#include <ctime>
#include "mpi.h"

#define MainProc 0

void swap(uint32_t *a, uint32_t *b) {
    uint32_t tmp = *a;
    *a = *b;
    *b = tmp;
}

uint32_t* merge(uint32_t *arr1, uint32_t num1, uint32_t *arr2, uint32_t num2) {
    uint32_t *res;
    uint32_t i = 0,
        j = 0,
        index = 0;

    res = new uint32_t[num1 + num2];

    while (i < num1 && j < num2) {
        if (arr1[i] < arr2[j])
            res[index++] = arr1[i++];
        else
            res[index++] = arr2[j++];
    }

    while (i < num1)
        res[index++] = arr1[i++];

    while (j < num2)
        res[index++] = arr2[j++];
    return res;
}

/* sort unsigned ints */
void rad_sort_u(uint32_t *from, uint32_t *to, unsigned bit) {
    if (!bit || to < from + 1) return;

    unsigned *ll = from, *rr = to - 1;

    for (;;) {
        /* find left most with bit, and right most without bit, swap */
        while (ll < rr && !(*ll & bit)) ll++;
        while (ll < rr && (*rr & bit)) rr--;
        if (ll >= rr) break;
        swap(ll, rr);
    }

    if (!(bit & *ll) && ll < to) ll++;
    bit >>= 1;

    rad_sort_u(from, ll, bit);
    rad_sort_u(ll, to, bit);
}

bool is_sorted(uint32_t* arr, uint32_t size) {
    bool flag = 1;
    for (uint32_t i = 0; i < size - 1; i++)
        if (arr[i] > arr[i + 1]) {
            flag = 0;
            break;
        }
    return flag;
}

int main(int argc, char** argv) {
    int status = 0, rank = 0, size = 0;  // MPI vars
    int *scounts = NULL,  // use in scatterv and gatherv
        *displs = NULL;   // use in scatterv and gatherv
    uint32_t arrSize = 0;
    uint32_t  *arr = NULL,  // data
              *buff = NULL,  // buffers for message exchanging
              *buff2 = NULL,
              *resultS = NULL,  // results
              *resultP = NULL;
    double t1 = 0, t2 = 0;
    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
     return -1;
    }

    status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (status != MPI_SUCCESS) {
        return -1;
    }

    status = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (status != MPI_SUCCESS) {
        return -1;
    }

    if (rank == MainProc) {
        if (argc > 1)
            arrSize = atoi(argv[1]);
        else
            arrSize = 10;
        std::cout << "  Array size : " << arrSize << std::endl;
        arr = new uint32_t[arrSize];
        resultS = new uint32_t[arrSize];
        resultP = new uint32_t[arrSize];
        std::srand((unsigned)time(NULL));
        for (uint32_t i = 0; i < arrSize; i++)
            arr[i] = resultS[i] = std::rand() % 201;
        if (arrSize < 50) {
            for (uint32_t i = 0; i < arrSize; i++)
                std::cout << arr[i] << " ";
            std::cout << std::endl;
        }
        t1 = MPI_Wtime();
        rad_sort_u(resultS, resultS + arrSize, INT_MIN);
        t2 = MPI_Wtime();
        std::cout << "  Sequential Time: " << t2 - t1 << std::endl;
        if (bool f = is_sorted(resultS, arrSize))
            std::cout << "    Array is sorted, flag = " << f << std::endl;
        else
            std::cout << "    Array isn`t sorted, flag = " << f << std::endl;
        if (arrSize < 50) {
            for (uint32_t i = 0; i < arrSize; i++)
                std::cout << resultS[i] << " ";
            std::cout << std::endl;
        }
        buff = new uint32_t[arrSize / size + arrSize % size];
        buff2 = new uint32_t[arrSize / size ];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    MPI_Bcast(&arrSize, 1, MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);
    scounts = new int[size];
    displs = new int[size];
    displs[0] = 0;
    scounts[0] = arrSize / size + arrSize % size;
    for (int i = 1; i < size; i++) {
        scounts[i] = arrSize / size;
        displs[i] = scounts[0] + (i-1)*scounts[i];
    }
    if (rank != 0)
        buff = new uint32_t[scounts[rank]];

    MPI_Scatterv(arr, scounts, displs, MPI_UNSIGNED, buff,
        scounts[rank], MPI_UNSIGNED, MainProc, MPI_COMM_WORLD);

    rad_sort_u(buff, buff + scounts[rank], INT_MIN);

    if (rank != 0)
        MPI_Send(buff, scounts[rank], MPI_UNSIGNED, MainProc,
            rank, MPI_COMM_WORLD);


    if (rank == 0) {
        for (uint32_t i = 0; i < scounts[0]; i++)
            resultP[i] = buff[i];
        for (uint32_t i = 1; i < size; i++) {
            MPI_Recv(buff2, arrSize / size, MPI_UNSIGNED, MPI_ANY_SOURCE,
                i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            resultP = merge(buff2, arrSize / size, resultP,
               scounts[0] + (i - 1)*scounts[1]);
        }
        t2 = MPI_Wtime();
        std::cout << "  Parallel Time: " << t2 - t1 << std::endl;
        if (bool f = is_sorted(resultS, arrSize))
            std::cout << "    Array is sorted, flag = " << f << std::endl;
        else
            std::cout << "    Array isn`t sorted, flag = " << f << std::endl;
        if (arrSize < 50) {
            for (uint32_t i = 0; i < arrSize; i++)
                std::cout << resultP[i] << " ";
            std::cout << std::endl;
        }
        delete[] buff2;
        delete[] resultS;
        delete[] resultP;
        delete[] arr;
    }
    delete[] buff;
    delete[] scounts;
    delete[] displs;
    status = MPI_Finalize();
    if (status != MPI_SUCCESS) {
        return -1;
    }
    return 0;
}
