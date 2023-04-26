#pragma once

#include <stddef.h>

/*
    Functions used to execute the genetic algorithm
*/

/* Solution format */
typedef struct {
    size_t chrom_len;   // length in whatever datatype
    char dead, elite;
    unsigned int generation, fitness, last_fitness;
    void *chromosome;
} ga_solution_t;

/* ga_select criteria */
#define GA_MAXIMIZE 0
#define GA_MINIMIZE 1

// Creates a new randomly generated population
// chrom_gen_func is a function that initializes a block of memory used for storing the gene pool.
// The block chrom_chunk must be big enough to contain size * chrom_len bytes
void ga_init(ga_solution_t *pop, size_t size, size_t chrom_len, void *chrom_chunk,
             void (*chrom_gen_func)(ga_solution_t *solution, size_t i, size_t chrom_len, void *chrom_chunk));

// Evaluates every solution in the population using the given function
void ga_eval(ga_solution_t *pop, size_t size, unsigned int (*fitness_func)(ga_solution_t *));

// Sorts population by fitness depending on the given criteria, then marks solutions not selected as dead
// O(size log size) / O(n^2) (worst)
void ga_select(ga_solution_t *pop, size_t size, int criteria, int percent_dead, int percent_elite);

// Creates the next generation by replacing dead solutions
// mutation_chance works as a 1 in n. E.g. 1 in a billion
// O(size)
void ga_next_generation(ga_solution_t *pop, size_t size, int percent_dead,
                        int percent_cross, void (*crossing_func)(ga_solution_t *, ga_solution_t *, ga_solution_t *),
                        int mutation_chance, void (*mutation_func)(ga_solution_t *));

// Retrieves some fitness information about the population
// O(size)
void ga_gen_info(ga_solution_t *pop, size_t size, int percent_elite, int *best, int *worst_elite,
                 int *average, int *worst);

