#pragma once

#include <stddef.h>

/*
    Functions used to execute the genetic algorithm
*/

/* Solution format */
typedef struct ga_solution {
    size_t chrom_len;
    char dead, elite;
    unsigned int generation, fitness, last_fitness;
    void *chromosome;
} ga_solution_t;

/* ga_select criteria */
#define GA_MAXIMIZE 0
#define GA_MINIMIZE 1

// Creates a new randomly generated population
void ga_init(ga_solution_t *pop, size_t size, size_t chrom_len, void (*chrom_gen_func)(ga_solution_t *solution));

// Evaluates every solution in the population using the given function
void ga_eval(ga_solution_t *pop, size_t size, unsigned int (*fitness_func)(ga_solution_t *));

// Sorts population by fitness depending on the given criteria, then marks solutions not selected as dead
void ga_select(ga_solution_t *pop, size_t size, int criteria, int percent_dead, int percent_elite);

// Creates the next generation by replacing dead solutions
void ga_next_generation(ga_solution_t *pop, size_t size,
                        void (*crossing_func)(ga_solution_t *, ga_solution_t *, ga_solution_t *),
                        void (*mutation_func)(ga_solution_t *, ga_solution_t *));

