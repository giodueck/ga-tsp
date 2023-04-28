#pragma once

#include <stddef.h>
#include <stdint.h>

/*
    Functions used to execute the genetic algorithm
*/

/* Solution format */
typedef struct {
    size_t chrom_len, gene_size;
    char dead, elite;
    unsigned int generation;
    int64_t fitness;
    void *chromosome;
} ga_solution_t;

/* ga_select criteria */
#define GA_MAXIMIZE 0
#define GA_MINIMIZE 1

// Creates a new randomly generated population
// chrom_gen_func is a function that initializes a block of memory used for storing the gene pool.
// The block chrom_chunk must be big enough to contain size * chrom_len bytes
void ga_init(ga_solution_t *pop,
             size_t size,
             size_t chrom_len,
             size_t gene_size,
             void *chrom_chunk,
             void (*chrom_gen_func)(ga_solution_t *solution, size_t i, size_t chrom_len, void *chrom_chunk));

// Evaluates every solution in the population using the given function
void ga_eval(ga_solution_t *pop, size_t size, int64_t (*fitness_func)(ga_solution_t *));

// Sorts population by fitness depending on the given criteria, then marks solutions not selected as dead
// O(size log size) / O(n^2) (worst)
void ga_select_trunc(ga_solution_t *pop, size_t size, int criteria, int percent_dead, int percent_elite);

// Creates the next generation by replacing dead solutions
// mutation_chance is a number in a million (actually 1024*1024)
// O(size)
int ga_next_generation(ga_solution_t *pop,
                       size_t size,
                       int percent_dead,
                       int percent_cross,
                       void (*crossing_func)(ga_solution_t *, ga_solution_t *, ga_solution_t *),
                       int mutation_per_Mi,
                       void (*mutation_func)(ga_solution_t *, int));

// Retrieves some fitness information about the population. Requires pop to be
// sorted by fitness
// O(size)
void ga_gen_info(ga_solution_t *pop,
                 size_t size,
                 int percent_elite,
                 int64_t *best,
                 int64_t *worst_elite,
                 int64_t *average,
                 int64_t *worst);

