#include "genetic.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// Creates a new randomly generated population
// chrom_gen_func is a function that initializes a block of memory used for storing the gene pool.
// The block chrom_chunk must be big enough to contain size * chrom_len bytes
void ga_init(ga_solution_t *pop,
             size_t size,
             size_t chrom_len,
             size_t gene_size,
             void *chrom_chunk,
             void (*chrom_gen_func)(ga_solution_t *solution, size_t i, size_t chrom_len, void *chrom_chunk, uint8_t *marks))
{
    uint8_t *marks = (uint8_t *) malloc(sizeof(uint8_t) * chrom_len);
    for (size_t i = 0; i < size; i++)
    {
        pop[i] = (ga_solution_t) { .chrom_len = chrom_len, .gene_size = gene_size, .dead = 0, .elite = 0, .generation = 0, .fitness = 0, .fit_gen = 0, NULL };
        chrom_gen_func(&(pop[i]), i, chrom_len, chrom_chunk, marks);
    }
    free(marks);
}

// Evaluates every solution in the population using the given function
void ga_eval(ga_solution_t *pop, size_t size, int64_t (*fitness_func)(ga_solution_t *))
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
void ga_select_trunc(ga_solution_t *pop, size_t size, int criteria, int percent_dead, int percent_elite, int64_t (*fitness_func)(ga_solution_t *))
{
    if (!size)
        return;

    // Evaluate every solution
    ga_eval(pop, size, fitness_func);

    // Fittest in front
    qsort(pop, size, sizeof(ga_solution_t), criteria ? ga_min : ga_max);

    // Kill lowest fitness solutions
    int dead_count = size * percent_dead / 100;
    for (size_t i = 0; i < size; i++)
    {
        pop[i].dead = i >= (size - dead_count);
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
                                  struct drand48_data *rbuf)
{
    /* Alg:
        Mark all solutions as not dead
        Let N = number of solutions replaced
        For 1..size/2k
            Hold 2 tournaments
                Select k live individuals
            Winners are parents, last-place-losers become offspring
            Mark contestants as dead to avoid repeated selection
        Increase generation
    */

    if (!size)
        return 0;

    for (size_t i = 0; i < size; i++)
        pop[i].dead = 0;

    if (k < 2)
        k = 2;

    int *contestants = (int *) malloc (sizeof(int) * k);
    int64_t *fits = (int64_t *) malloc(sizeof(int64_t) * k);
    uint8_t *marks = (uint8_t *) malloc(sizeof(uint8_t) * pop->chrom_len);
    long lrand;

    // Number of tournaments. Lower k means more individuals get replaced
    // per generation, but higher k means weak individuals win less often
    int N = size / (k * 2);

    for (int n = 0; n < N; n++)
    {
        int p1 = 0, p2 = 0;
        int c1 = 0, c2 = 0;
        int64_t low = 0;
        int64_t high = 0;

        // Select contestants
        for (int i = 0; i < k; i++)
        {
            lrand48_r(rbuf, &lrand);
            int pot = lrand % size;
            while (pop[pot].dead)
                pot = (pot + 1) % size;
            contestants[i] = pot;
            pop[pot].dead = 1;
        }

        // Evaluate
        for (int i = 0; i < k; i++)
            fits[i] = fitness_func(&pop[contestants[i]]);

        low = fits[0];
        high = low;
        p1 = contestants[0];
        c1 = p1;
        for (int i = 1; i < k; i++)
        {
            if (fits[i] < low)
            {
                low = fits[i];
                p1 = contestants[i];
            }
            if (fits[i] > high)
            {
                high = fits[i];
                c1 = contestants[i];
            }
        }

        // Select contestants (round 2)
        for (int i = 0; i < k; i++)
        {
            lrand48_r(rbuf, &lrand);
            int pot = lrand % size;
            while (pop[pot].dead)
                pot = (pot + 1) % size;
            contestants[i] = pot;
            pop[pot].dead = 1;
        }

        p2 = contestants[0];
        c2 = p2;
        low = fitness_func(&pop[contestants[0]]);
        high = low;

        // Evaluate (round 2)
        for (int i = 0; i < k; i++)
            fits[i] = fitness_func(&pop[contestants[i]]);

        // Pick parents and losers (offspring)
        low = fits[0];
        high = low;
        p2 = contestants[0];
        c2 = p2;
        for (int i = 1; i < k; i++)
        {
            if (criteria == GA_MINIMIZE)
            {
                if (fits[i] < low)
                {
                    low = fits[i];
                    p2 = contestants[i];
                }
                if (fits[i] > high)
                {
                    high = fits[i];
                    c2 = contestants[i];
                }
            }

            if (criteria == GA_MAXIMIZE)
            {
                if (fits[i] > low)
                {
                    low = fits[i];
                    p2 = contestants[i];
                }
                if (fits[i] < high)
                {
                    high = fits[i];
                    c2 = contestants[i];
                }
            }
        }

        // Create offspring
        crossing_func(&(pop[p1]), &(pop[p2]), &(pop[c1]), marks, rbuf);
        mutation_func(&(pop[c1]), mutation_per_Mi, rbuf);
        pop[c1].fit_gen = 0;
        fitness_func(&pop[c1]);
        
        crossing_func(&(pop[p2]), &(pop[p1]), &(pop[c2]), marks, rbuf);
        mutation_func(&(pop[c2]), mutation_per_Mi, rbuf);
        pop[c2].fit_gen = 0;
        fitness_func(&pop[c2]);
    } 
    free(contestants);
    free(fits);
    free(marks);

    for (size_t i = 0; i < size; i++)
        pop[i].generation++;

    return pop->generation;
}

// Retrieves some fitness information about the population. Requires pop to be
// sorted by fitness
// O(size)
void ga_gen_info(ga_solution_t *pop,
                 size_t size,
                 int percent_elite,
                 int64_t *best,
                 int64_t *worst_elite,
                 int64_t *average,
                 int64_t *worst)
{
    if (!size)
        return;
    
    if (best) *best = pop[0].fitness;
    
    if (worst) *worst = pop[size - 1].fitness;
    
    if (percent_elite && worst_elite && pop[size * percent_elite / 100 - 1].elite)
        *worst_elite = pop[size * percent_elite / 100 - 1].fitness;
    
    if (!average) return;
    int64_t aux = 0;
    for (size_t i = 0; i < size; i++)
        aux += pop[i].fitness;
    *average = aux / size;
}

void ga_gen_info_unsorted(ga_solution_t *pop,
                          size_t size,
                          int percent_elite,
                          int64_t *best,
                          int64_t *average,
                          int64_t *worst)
{
    if (!size)
        return;

    if (best) *best = pop[0].fitness;
    if (worst) *worst = pop[0].fitness;
    if (average) *average = 0;
    for (size_t i = 0; i < size; i++)
    {
        if (best && pop[i].fitness < *best)
            *best = pop[i].fitness;
        if (worst && pop[i].fitness > *worst)
            *worst = pop[i].fitness;
        if (average)
            *average += pop[i].fitness;
    }
    *average /= size;
}
