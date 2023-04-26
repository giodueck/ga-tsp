#include "genetic.h"
#include <stdlib.h>
#include <string.h>

// Creates a new randomly generated population
// chrom_gen_func is a function that initializes a block of memory used for storing the gene pool.
// The block chrom_chunk must be big enough to contain size * chrom_len bytes
void ga_init(ga_solution_t *pop, size_t size, size_t chrom_len, void *chrom_chunk,
             void (*chrom_gen_func)(ga_solution_t *solution, size_t i, size_t chrom_len, void *chrom_chunk))
{
    for (size_t i = 0; i < size; i++)
    {
        pop[i] = (ga_solution_t) { .chrom_len = chrom_len, .dead = 0, .elite = 0, .generation = 0, .fitness = 0, .last_fitness = 0, NULL };
        chrom_gen_func(&(pop[i]), i, chrom_len, chrom_chunk);
    }
}

// Evaluates every solution in the population using the given function
void ga_eval(ga_solution_t *pop, size_t size, unsigned int (*fitness_func)(ga_solution_t *))
{
    for (size_t i = 0; i < size; i++)
        pop[i].fitness = fitness_func(&(pop[i]));
}

// Maximizing fitness means first in the array is biggest, return b compared to a
int ga_max(const void *a, const void *b)
{
    return ((ga_solution_t *) b)->fitness - ((ga_solution_t *) a)->fitness;
}

// Minimizing fitness means first in the array is smallest, return a compared to b
int ga_min(const void *a, const void *b)
{
    return ((ga_solution_t *) a)->fitness - ((ga_solution_t *) b)->fitness;
}

// Sorts population by fitness depending on the given criteria, then marks solutions not selected as dead
// O(size log size) / O(n^2) (worst)
void ga_select(ga_solution_t *pop, size_t size, int criteria, int percent_dead, int percent_elite)
{
    // Fittest in front
    qsort(pop, size, sizeof(ga_solution_t), criteria ? ga_min : ga_max);

    // Kill lowest fitness solutions
    int dead_count = size * percent_dead / 100;
    for (size_t i = 0; i < dead_count; i++)
    {
        pop[size - 1 - i].dead = 1;
    }

    // Highest fitness solutions are made elite
    if (percent_elite)
    {
        int elite_count = size * percent_elite / 100;
        for (size_t i = 0; i < elite_count; i++)
            pop[i].elite = 1;
        for (size_t i = elite_count; i < size - dead_count; i++)
            pop[i].elite = 0;
    }
}

// Creates the next generation by replacing dead solutions
// mutation_chance works as a 1 in n. E.g. 1 in a billion
// O(size)
void ga_next_generation(ga_solution_t *pop, size_t size, size_t chrom_len, int percent_dead,
                        int percent_cross, void (*crossing_func)(ga_solution_t *, ga_solution_t *, ga_solution_t *),
                        int mutation_chance, void (*mutation_func)(ga_solution_t *, size_t))
{
    size_t threshold = size * percent_dead / 100;
    size_t cross = size * percent_cross / 100;
    int rn;
    for (size_t i = threshold; i < size; i++)
    {
        // Reproduce
        if (cross)
        {
            int _i = rand() % threshold, _j = rand() % threshold;
            crossing_func(&(pop[_i]), &(pop[_j]), &(pop[i]));
            cross--;
        }
        else
        {
            int _i = rand() % threshold;
            memcpy(&(pop[i]), &(pop[_i]), sizeof(ga_solution_t) - sizeof(void *));
            memcpy(pop[i].chromosome, pop[_i].chromosome, chrom_len);
        }

        // Sometimes mutate
        if (rand() % mutation_chance == 0)
            mutation_func(&(pop[i]), chrom_len);
    }
}

// Retrieves some fitness information about the population
// O(size)
void ga_gen_info(ga_solution_t *pop, size_t size, int percent_elite, int *best, int *worst_elite,
                 int *average, int *worst)
{
    *best = pop[0].fitness;
    *worst = pop[size - 1].fitness;
    if (percent_elite)
        *worst_elite = pop[size * percent_elite / 100 - 1].fitness;
    int aux = 0;
    for (size_t i = 0; i < size; i++)
        aux += pop[i].fitness;
    *average = aux / size;
}

