// copyright : (C) by diper1998
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <string>

const int LIMIT = 5;

int main(int argc, char *argv[]) {
    srand(static_cast<int>(time(0)));

    int Tacts = 10;

    int Producer = 5;
    int Consumer = 5;

    int CountProducer = 1;
    int CountConsumer = 1;

    int pdR = 0;
    int csR = 0;

    int sph1 = 0;
    int sph2 = 0;

    int flag;
    int myId, numProcs;

    int flagConsumer = 1;
    int flagProducer = 1;

    int borderCP;

    int criticalSection = 0;

    std::string sizeStr;

    sizeStr = argv[1];
    Tacts = atoi(sizeStr.c_str());

    sizeStr = argv[2];
    Producer = atoi(sizeStr.c_str());

    sizeStr = argv[3];
    Consumer = atoi(sizeStr.c_str());

    sizeStr = argv[4];
    CountProducer = atoi(sizeStr.c_str());

    sizeStr = argv[5];
    CountConsumer = atoi(sizeStr.c_str());

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    MPI_Initialized(&flag);  // check
    if (!flag) {
        std::cout << "Error: MPI_Init";
    }

    // Description of the communicator
    // communicator manages groups of parallel processes
    MPI_Comm_size(MPI_COMM_WORLD,
        &numProcs);  // determining the number of processes in a group
    MPI_Comm_rank(MPI_COMM_WORLD,
        &myId);  // determining the rank of a process in a group

    MPI_Status status;

    MPI_Request rqSendProducer;
    MPI_Request rqRecvProducer;
    MPI_Request rqSendConsumer;
    MPI_Request rqRecvConsumer;

    MPI_Request rqSendProducerFlag;
    MPI_Request rqRecvProducerFlag;
    MPI_Request rqSendConsumerFlag;
    MPI_Request rqRecvConsumerFlag;

    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_int_distribution<int> dis(0, 100);

    borderCP = CountProducer;

    if (CountConsumer + CountProducer + 1 != numProcs) {
        std::cout << "Error" << std::endl;
    } else {
        for (int i = 0; i < Tacts; i++) {
            if (myId == 0) {  // Commander
                std::cout << "Tact: " << i << std::endl;

                // Probability



                for (int j = 0; j < 12; j++) {
                    // pdR += dis(gen);
                    pdR += std::rand();
                }
                pdR = pdR % 10 + 1;

                for (int j = 0; j < 12; j++) {
                    csR += std::rand();
                }
                csR = csR % 10 + 1;

                // Worker selection

                if (flagProducer) {
                    if (pdR <= Producer) {
                        sph1 = std::rand() % CountProducer + 1;
                    } else {
                        sph1 = 0;
                    }
                }

                if (flagConsumer) {
                    if (csR <= Consumer) {
                        sph2 = std::rand() % CountConsumer +
                            CountProducer + 1;
                    } else {
                        sph2 = 0;
                    }
                }

                // Data transfer to Producers
                for (int j = 1; j <= CountProducer; j++) {
                    MPI_Isend(&sph1, 1, MPI_INT, j, 94, MPI_COMM_WORLD,
                        &rqSendProducer);
                    MPI_Wait(&rqSendProducer, &status);
                    MPI_Isend(&criticalSection, 1, MPI_INT, j, 95,
                        MPI_COMM_WORLD,
                        &rqSendProducer);
                    MPI_Wait(&rqSendProducer, &status);
                }

                // Getting data from working Producer
                if (sph1 != 0) {
                    MPI_Irecv(&criticalSection, 1, MPI_INT, sph1, 99,
                        MPI_COMM_WORLD,
                        &rqRecvProducer);
                    MPI_Wait(&rqRecvProducer, &status);

                    MPI_Irecv(&flagProducer, 1, MPI_INT, sph1, 97,
                        MPI_COMM_WORLD,
                        &rqRecvProducerFlag);
                    MPI_Wait(&rqRecvProducerFlag, &status);
                }

                // Data transfer to Consumers
                for (int j = borderCP + 1; j < numProcs; j++) {
                    MPI_Isend(&sph2, 1, MPI_INT, j, 93, MPI_COMM_WORLD,
                        &rqSendConsumer);
                    MPI_Wait(&rqSendConsumer, &status);
                    MPI_Isend(&criticalSection, 1, MPI_INT, j, 92,
                        MPI_COMM_WORLD,
                        &rqSendConsumer);
                    MPI_Wait(&rqSendConsumer, &status);
                }

                // Getting data from working Consumer
                if (sph2 != 0) {
                    MPI_Irecv(&criticalSection, 1, MPI_INT, sph2, 98,
                        MPI_COMM_WORLD,
                        &rqRecvConsumer);
                    MPI_Wait(&rqRecvConsumer, &status);

                    MPI_Irecv(&flagConsumer, 1, MPI_INT, sph2, 96,
                        MPI_COMM_WORLD,
                        &rqRecvConsumerFlag);
                    MPI_Wait(&rqRecvConsumerFlag, &status);
                }
            }  // end Commander

            // Getting data from Commander to Producer
            if (myId <= borderCP && myId != 0) {
                MPI_Irecv(&sph1, 1, MPI_INT, 0, 94, MPI_COMM_WORLD,
                    &rqRecvProducer);
                MPI_Wait(&rqRecvProducer, &status);

                MPI_Irecv(&criticalSection, 1,

                    MPI_INT, 0, 95, MPI_COMM_WORLD, &rqRecvProducer);
                MPI_Wait(&rqRecvProducer, &status);
            }

            // Working Producer
            if (sph1 != 0 && myId == sph1) {
                flagProducer = 1;

                while (criticalSection >= LIMIT && i < Tacts) {
                    std::cout << "PRODUCER[" << myId << "]: "
                              << criticalSection << " WAITING" << std::endl;
                    MPI_Isend(&criticalSection, 1, MPI_INT, 0, 99,
                        MPI_COMM_WORLD,
                        &rqSendProducer);
                    MPI_Wait(&rqSendProducer, &status);
                    flagProducer = 0;
                    MPI_Isend(&flagProducer, 1, MPI_INT, 0, 97,
                        MPI_COMM_WORLD,
                        &rqSendProducerFlag);
                    MPI_Wait(&rqSendProducerFlag, &status);

                    i++;

                    if (i < Tacts) {
                        MPI_Irecv(&sph1, 1, MPI_INT, 0, 94,
                            MPI_COMM_WORLD,
                            &rqRecvProducer);
                        MPI_Wait(&rqRecvProducer, &status);

                        MPI_Irecv(&criticalSection, 1, MPI_INT, 0, 95,
                            MPI_COMM_WORLD,
                            &rqRecvProducer);
                        MPI_Wait(&rqRecvProducer, &status);

                        if (criticalSection < LIMIT) {
                            flagProducer = 1;
                        }
                    }
                }

                if (flagProducer) {
                    criticalSection++;
                    std::cout << "PRODUCER[" << myId << "]: "
                        << criticalSection << std::endl;
                    MPI_Isend(&criticalSection, 1, MPI_INT, 0, 99,
                        MPI_COMM_WORLD,
                        &rqSendProducer);
                    MPI_Wait(&rqSendProducer, &status);
                    flagProducer = 1;
                    MPI_Isend(&flagProducer, 1, MPI_INT, 0, 97,
                        MPI_COMM_WORLD,
                        &rqSendProducerFlag);
                    MPI_Wait(&rqSendProducerFlag, &status);
                } else {
                    MPI_Isend(&criticalSection, 1, MPI_INT, 0, 99,
                        MPI_COMM_WORLD,
                        &rqSendProducer);
                    MPI_Wait(&rqSendProducer, &status);
                    flagProducer = 0;
                    MPI_Isend(&flagProducer, 1, MPI_INT, 0, 97,
                        MPI_COMM_WORLD,
                        &rqSendProducerFlag);
                    MPI_Wait(&rqSendProducerFlag, &status);
                }
            }

            // Getting data from Commander to Consumer
            if (myId > borderCP) {
                MPI_Irecv(&sph2, 1, MPI_INT, 0, 93, MPI_COMM_WORLD,
                    &rqRecvConsumer);
                MPI_Wait(&rqRecvConsumer, &status);

                MPI_Irecv(&criticalSection, 1, MPI_INT, 0, 92,
                    MPI_COMM_WORLD,
                    &rqRecvConsumer);
                MPI_Wait(&rqRecvConsumer, &status);
            }

            // Working Consumer
            if (sph2 != 0 && myId == sph2) {
                flagConsumer = 1;

                while (criticalSection <= 0 && i < Tacts) {
                    std::cout << "CONSUMER[" << myId << "]: "
                        << criticalSection << " WAITING" << std::endl;
                    MPI_Isend(&criticalSection, 1, MPI_INT, 0, 98,
                        MPI_COMM_WORLD,
                        &rqSendConsumer);
                    MPI_Wait(&rqSendConsumer, &status);
                    flagConsumer = 0;
                    MPI_Isend(&flagConsumer, 1, MPI_INT, 0, 96,
                        MPI_COMM_WORLD,
                        &rqSendConsumerFlag);
                    MPI_Wait(&rqSendConsumerFlag, &status);

                    i++;

                    if (i < Tacts) {
                        MPI_Irecv(&sph2, 1, MPI_INT, 0, 93, MPI_COMM_WORLD,
                            &rqRecvConsumer);
                        MPI_Wait(&rqRecvConsumer, &status);

                        MPI_Irecv(&criticalSection, 1, MPI_INT, 0, 92,
                            MPI_COMM_WORLD,
                            &rqRecvConsumer);
                        MPI_Wait(&rqRecvConsumer, &status);

                        if (criticalSection > 0) {
                            flagConsumer = 1;
                        }
                    }
                }

                if (flagConsumer) {
                    criticalSection--;
                    std::cout << "CONSUMER[" << myId << "]: "
                        << criticalSection << std::endl;
                    MPI_Isend(&criticalSection, 1, MPI_INT, 0, 98,
                        MPI_COMM_WORLD,
                        &rqSendConsumer);
                    MPI_Wait(&rqSendConsumer, &status);
                    flagConsumer = 1;
                    MPI_Isend(&flagConsumer, 1, MPI_INT, 0, 96,
                        MPI_COMM_WORLD,
                        &rqSendConsumerFlag);
                    MPI_Wait(&rqSendConsumerFlag, &status);
                } else {
                    MPI_Isend(&criticalSection, 1, MPI_INT, 0, 98,
                        MPI_COMM_WORLD,
                        &rqSendConsumer);
                    MPI_Wait(&rqSendConsumer, &status);
                    flagConsumer = 0;
                    MPI_Isend(&flagConsumer, 1, MPI_INT, 0, 96,
                        MPI_COMM_WORLD,
                        &rqSendConsumerFlag);
                    MPI_Wait(&rqSendConsumerFlag, &status);
                }
            }
        }
    }

    MPI_Finalize();
    return 0;
}
