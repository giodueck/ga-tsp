#!/bin/bash

# For debug mode add "-DDEBUG"
# This will cause the master node to print its PID, which can be attached to via "$ gdb --pid <PID>"
# TODO: Make it print the hostname as well
mpicc -Wall -o ga-tsp-mpi main.c genetic.c tsp_parser.c tsp.c -lrt -lm -DMPI $1
