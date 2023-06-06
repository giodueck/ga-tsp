#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

/*
    Functions used to execute the genetic algorithm
*/

/* Solution format */
typedef struct {
    size_t chrom_len, gene_size;
    char dead, elite;
    unsigned int generation;
    int64_t fitness;
    unsigned int fit_gen;   // auxiliary to help caching fitness
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
             void (*chrom_gen_func)(ga_solution_t *solution, size_t i, size_t chrom_len, void *chrom_chunk, uint8_t *marks));

// Evaluates every solution in the population using the given function
void ga_eval(ga_solution_t *pop, size_t size, int64_t (*fitness_func)(ga_solution_t *));

// Sorts population by fitness depending on the given criteria, then marks solutions not selected as dead
// O(size log size) / O(n^2) (worst)
void ga_select_trunc(ga_solution_t *pop, size_t size, int criteria, int percent_dead, int percent_elite, int64_t (*fitness_func)(ga_solution_t *));

// Creates the next generation by replacing dead solutions
// mutation_chance is a number in a million (actually 1024*1024)
// O(size)
/* int ga_next_generation_trunc(ga_solution_t *pop,
                             size_t size,
                             int percent_dead,
                             int percent_cross,
                             void (*crossing_func)(ga_solution_t *, ga_solution_t *, ga_solution_t *, uint8_t *, struct drand48_data *),
                             int mutation_per_Mi,
                             void (*mutation_func)(ga_solution_t *, int, struct drand48_data *),
                             struct drand48_data *rbuf); */

// Creates tournaments of size k where the fittest individuals get to procreate, while losers
// are replaced with offspring. If k >= 4, the parents are selected in one tournament and
// the least fit losers are replaced with the offspring, otherwise two tournaments are held
// which each yield one parent and one offspring. Offspring do not participate in the current tournament
int ga_next_generation_tournament(ga_solution_t *pop,
                                  size_t size,
                                  int k,
                                  int criteria,
                                  int64_t (*fitness_func)(ga_solution_t *i),
                                  void (*crossing_func)(ga_solution_t *, ga_solution_t *, ga_solution_t *, uint8_t *, struct drand48_data *),
                                  int mutation_per_Mi,
                                  void (*mutation_func)(ga_solution_t *, int, struct drand48_data *),
                                  struct drand48_data *rbuf);

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

// Retrieves some fitness information about the population.
// O(size)
void ga_gen_info_unsorted(ga_solution_t *pop,
                          size_t size,
                          int percent_elite,
                          int64_t *best,
                          int64_t *average,
                          int64_t *worst);
