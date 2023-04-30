#pragma once

#include <stdlib.h>
#include <stdint.h>
#include "genetic.h"
#include "tsp_parser.h"

// Initializes a random solution
void generate_tsp_solution(ga_solution_t *sol, size_t i, size_t chrom_len, void *chrom_chunk, uint8_t *marks);

// Euclidean distance
double dist(const tsp_2d_node_t a, const tsp_2d_node_t b);

// Distance based fitness
int64_t fitness(ga_solution_t *sol);

// Cross two solutions and produce a child solution with traits from both parents 
void crossover(ga_solution_t *p1, ga_solution_t *p2, ga_solution_t *child, uint8_t *marks);

// Apply random swaps of genes dictated by some small chance
void mutate(ga_solution_t *sol, int per_Mi);

