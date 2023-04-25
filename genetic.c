#include "genetic.h"
#include <stdlib.h>

// Creates a new randomly generated population
// TODO: Think about allocation: tons of small blocks are inefficient, do one large one and divide it up
void ga_init(ga_solution_t *pop, size_t size, size_t chrom_len, void (*chrom_gen_func)(ga_solution_t *solution))
{
    for (size_t i = 0; i < size; i++)
    {
        pop[i] = (ga_solution_t) { .chrom_len = chrom_len, .dead = 0, .elite = 0, .generation = 0, .fitness = 0, .last_fitness = 0, NULL };
        chrom_gen_func(&(pop[i]));
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
void ga_next_generation(ga_solution_t *pop, size_t size,
                        void (*crossing_func)(ga_solution_t *, ga_solution_t *, ga_solution_t *),
                        void (*mutation_func)(ga_solution_t *, ga_solution_t *))
{

}

