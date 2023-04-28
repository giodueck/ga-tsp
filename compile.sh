#!/bin/bash

gcc -Wall -o ga-tsp main.c genetic.c tsp_parser.c -lrt -lm $1
