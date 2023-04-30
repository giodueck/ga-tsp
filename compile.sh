#!/bin/bash

gcc -Wall -o ga-tsp main.c genetic.c tsp_parser.c tsp.c -lrt -lm $1
