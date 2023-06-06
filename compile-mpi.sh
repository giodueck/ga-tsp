#!/bin/bash

mpicc -D MPI -Wall -o ga-tsp-mpi main.c genetic.c tsp_parser.c tsp.c -lrt -lm $1
