#!/bin/bash

gcc -o ga-tsp main.c genetic.c tsp_parser.c -lrt $1
