#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h"
#include "timing.h"

main(int argc, char **argv)
{

    /*
      This is the Hello World program for CPSC424/524.

      Author: Andrew Sherman, Yale University

      Date: 9/15/2014

      Credits: This program is based on a program provided by Barry Wilkinson (UNCC), which
               had a similar communication pattern, but did not include any simulated work.
    */

    int i, rank, size, type = 99;
    int worktime, sparm, rwork(int, int);
    double wct0, wct1, total_time, cput;
    char message[100];

    MPI_Status status;

    MPI_Init(&argc, &argv); // Required MPI initialization call

    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get no. of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Which process am I?

    /* If I am the manager (rank 0) ... */
    if (rank == 0)
    {
        char **messages = (char **)malloc((size - 1) * sizeof(char *));
        for (int i = 0; i < size - 1; ++i)
        {
            messages[i] = (char *)malloc(100 * sizeof(char));
        }

        sparm = rwork(0, 0); // initialize the workers' work times

        /* Create the message using sprintf */
        sprintf(message, "Hello, from process %d.", rank);

        MPI_Barrier(MPI_COMM_WORLD);
        wct0 = MPI_Wtime();

        /* Send the message to all the workers, which is where the work happens */
        for (i = 1; i < size; i++)
        {
            MPI_Send(message, strlen(message) + 1, MPI_CHAR, i, type, MPI_COMM_WORLD);
            MPI_Send(&sparm, 1, MPI_INT, i, type, MPI_COMM_WORLD);
        }

        // Receive messages from workers as they complete
        for (int i = 1; i < size; i++)
        {
            MPI_Recv(message, 100, MPI_CHAR, MPI_ANY_SOURCE, type, MPI_COMM_WORLD, &status);
            sleep(3); // simulate post processing
            strcpy(messages[status.MPI_SOURCE - 1], message);
        }

        /* print messages */
        for (size_t i = 1; i < size; i++)
        {
            printf("Message from process %d: %s\n", i, messages[i - 1]);
        }

        wct1 = MPI_Wtime();
        // timing(&wct1, &cput); // get the end time using the original timing routine
        total_time = wct1 - wct0;
        printf("Message printed by manager: Total elapsed time is %f seconds.\n", total_time);

        for (int i = 0; i < size - 1; i++)
        {
            free(messages[i]);
        }
        free(messages);
    }

    /* Otherwise, if I am a worker ... */
    else
    {
        MPI_Barrier(MPI_COMM_WORLD); // wait for everyone to be ready before starting

        /* Receive messages from the manager */
        MPI_Recv(message, 100, MPI_CHAR, 0, type, MPI_COMM_WORLD, &status);
        MPI_Recv(&sparm, 1, MPI_INT, 0, type, MPI_COMM_WORLD, &status);

        worktime = rwork(rank, sparm); // Simulate some work

        sprintf(message, "Hello manager, from process %d after working %d seconds.",
                rank, worktime);

        /* Send message back to manager */
        MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, type, MPI_COMM_WORLD);
    }

    MPI_Finalize(); // Required MPI termination call
}
