// Copyright 2018 Nesterov Alexander
#include <mpi.h>
#include <vector>
#include <iostream>
#include <random>
#include <cstdlib>
#include <cmath>

int main(int argc, char* argv[]) {
    int status, rank, size;
    const size_t count = atoi(argv[1]);

    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
        std::cout << "MPI_Init failed" << std::endl;
        return -1;
    }

    status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (status != MPI_SUCCESS) {
        std::cout << "MPI_Comm_rank failed" << std::endl;
        return -1;
    }

    status = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (status != MPI_SUCCESS) {
        std::cout << "MPI_Comm_size failed" << std::endl;
        return -1;
    }

    double *first_vector = nullptr;
    double *second_vector = nullptr;
    double parallel_start = 0.0;
    if (rank == 0) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-100, 100);
        first_vector  = new double[count];
        second_vector = new double[count];
        for (size_t i = 0; i < count; i++) {
            first_vector[i]  = dis(gen);
            second_vector[i] = dis(gen);
        }

        double start = MPI_Wtime();
        double sequential_result = 0.0;
        for (size_t i = 0; i < count; i++) {
            sequential_result += first_vector[i] * second_vector[i];
        }
        std::cout << "Sequential time: " << MPI_Wtime() - start << std::endl;
        std::cout << "Result of scalar multiplication = ";
        std::cout << sequential_result << std::endl;
        parallel_start = MPI_Wtime();
    }

    const size_t sub_count = static_cast<size_t>(ceil(count / size));
    double *sub_first_vector = nullptr;
    sub_first_vector  = new double[sub_count];
    double *sub_second_vector = nullptr;
    sub_second_vector = new double[sub_count];
    for (size_t i = 0; i < sub_count; i++) {
        sub_first_vector[i]  = 0;
        sub_second_vector[i] = 0;
    }

    MPI_Scatter(first_vector, sub_count, MPI_DOUBLE,
        sub_first_vector, sub_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Scatter(second_vector, sub_count, MPI_DOUBLE,
        sub_second_vector, sub_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double sub_result = 0.0;
    for (size_t i = 0; i < sub_count; i++) {
        sub_result += sub_first_vector[i] * sub_second_vector[i];
    }

    double parallel_result = 0.0;
    MPI_Reduce(&sub_result, &parallel_result, 1, MPI_DOUBLE,
        MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << std::endl << "Parallel time: ";
        std::cout << MPI_Wtime() - parallel_start << std::endl;
        std::cout << "Result of scalar multiplication = ";
        std::cout << parallel_result << std::endl;
    }

    if (sub_first_vector  != nullptr) { delete[]sub_first_vector; }
    if (sub_second_vector != nullptr) { delete[]sub_second_vector; }
    if (first_vector      != nullptr) { delete[]first_vector; }
    if (second_vector     != nullptr) { delete[]second_vector; }

    status = MPI_Finalize();
    if (status != MPI_SUCCESS) {
        std::cout << "MPI_Finalize failed" << std::endl;
        return -1;
    }

    return 0;
}
