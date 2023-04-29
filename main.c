#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "genetic.h"
#include "tsp_parser.h"

#define SEL_TRUNCATE   0
#define SEL_TOURNAMENT 1

tsp_2d_t tsp;

/* Parameters */
int population_size = 2500;     // population size per thread
int max_gens = 3000;            // when to stop the algorithm
int gen_info_interval = 100;    // how often to print information about the population
int mutations = 1000;           // mutations / (1024*1024) = mutation chance
int sel_strat = 1;
    /* Truncation selection */
int percent_elite = 5;          // percentage of elite selection, makes gen_info more informative for tournament
int percent_dead = 50;          // how many solutions get replaced 
int percent_cross = 50;         // how many solutions are derived from crossover
    /* Tournament selection */
int tournament_size = 4;        // how many individuals get picked per tournament

// Initializes a random solution
void generate_tsp_solution(ga_solution_t *sol, size_t i, size_t chrom_len, void *chrom_chunk)
{
    uint32_t *chromosome = (uint32_t *) chrom_chunk + i * chrom_len;
    uint8_t *marks = (uint8_t *) malloc(sizeof(uint8_t) * chrom_len);
    memset(marks, 0, sizeof(uint8_t) * chrom_len);

    for (size_t j = chrom_len; j; j--)
    {
        size_t r = rand() % j;
        size_t l = 0;
        while (marks[l] || r)
        {
            if (!marks[l] && r)
                r--;
            l++;
        }
        chromosome[j - 1] = l;
        marks[l] = 1;
    }

    sol->chromosome = (void *) chromosome;
    free(marks);
}

// Euclidean distance
double dist(tsp_2d_node_t a, tsp_2d_node_t b)
{
    double n = a.x - b.x;
    double m = a.y - b.y;
    return sqrt(n*n + m*m);
}

// Distance based fitness
int64_t fitness(ga_solution_t *sol)
{
    double d = 0;
    for (int i = 0; i < sol->chrom_len; i++)
    {
        int j = (i + 1) % sol->chrom_len;
        d += dist(tsp.nodes[((uint32_t *) sol->chromosome)[i]],
                  tsp.nodes[((uint32_t *) sol->chromosome)[j]]);
    }
    return (int64_t) d;
}

void mutate(ga_solution_t *sol, int per_Mi);

// Cross two solutions and produce a child solution with traits from both parents 
void crossover(ga_solution_t *p1, ga_solution_t *p2, ga_solution_t *child)
{
    // Take half of the chromosome of one parent, then the remaining half of the other such that
    // nodes don't repeat
    
    int start = rand() % (p1->chrom_len / 2);
    int l = p1->chrom_len / 2;
    char *marks = (char *) malloc(sizeof(char) * p1->chrom_len);
    memset(marks, 0, p1->chrom_len);

    // Copy half from parent 1
    for (int i = 0; i < l; i++)
    {
        uint32_t n = ((uint32_t *) p1->chromosome)[start + i];
        ((uint32_t *) child->chromosome)[i] = n;
        marks[n] = 1;
    }

    // Copy remaining
    for (int i = 0; i < p1->chrom_len; i++)
    {
        uint32_t n = ((uint32_t *) p2->chromosome)[i];
        if (marks[n])
            continue;

        ((uint32_t *) child->chromosome)[l++] = n;
    }

    // If solutions are very similar apply some high mutation rate
    int diff = 0;
    for (int i = 0; i < p1->chrom_len; i++)
        if (((uint32_t *) p1->chromosome)[i] != ((uint32_t *) p2->chromosome)[i])
            diff++;

    // If parents are less than 10% different
    if (diff <= p1->chrom_len / 10)
        mutate(child, mutations * 15); // 15 times as likely to have mutations
    // ^ This is funnily what happens in real life, but it gives better results in this case
}

// Apply random swaps of genes dictated by some small chance
void mutate(ga_solution_t *sol, int per_Mi)
{
    // 1024*1024 - 1
    // This is close enough to 1 million and a good mask to efficiently get small rand() numbers
    per_Mi &= 0xFFFFF;
    for (size_t i = 0; i < sol->chrom_len; i++)
    {
        int r = rand() & 0xFFFFF;
        if (r < per_Mi)
        {
            uint32_t aux = ((uint32_t*)sol->chromosome)[i];
            uint32_t j = rand() % sol->chrom_len;
            ((uint32_t*)sol->chromosome)[i] = ((uint32_t*)sol->chromosome)[j];
            ((uint32_t*)sol->chromosome)[j] = aux;
        }
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [file.tsp] <optional random seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    tsp = tsp_2d_read(argv[1]);

    if (argc >= 3)
        srand(atoi(argv[2]));

    uint32_t *chromosome_chunk = (uint32_t *) malloc(sizeof(uint32_t) * tsp.dim * population_size);
    ga_solution_t *population = (ga_solution_t *) malloc(sizeof(ga_solution_t) * population_size);

    /* Initialize population */
    ga_init(population, population_size, tsp.dim, sizeof(uint32_t), chromosome_chunk, generate_tsp_solution);

    int gen = 0;
    
    /* Evolve for max_gens number of generations */
    for (int i = 0; i < max_gens; i++)
    {
        if (sel_strat == SEL_TRUNCATE)
            /* Sort and select elite and survivors according to the truncation selection method */
            ga_select_trunc(population, population_size, GA_MINIMIZE, percent_dead, percent_elite, fitness);
        
        if ((gen + 1) % gen_info_interval == 0 || !gen)
        {
            int64_t best, worst_elite = 0, avg, worst;
            // Just to sort population, doesn't make any changes
            if (sel_strat != SEL_TRUNCATE)
                ga_select_trunc(population, population_size, GA_MINIMIZE, percent_dead, percent_elite, fitness);

            ga_gen_info(population, population_size, percent_elite, &best, &worst_elite, &avg, &worst);
            printf("%4d:\tB: %5lu\t%3d%%: %5lu\tA: %5lu\tW: %5lu\n", gen + 1, best, percent_elite, worst_elite, avg, worst);
        }

        if (sel_strat == SEL_TRUNCATE)
            /* Replace dead population with new offspring from surviving individuals */
            gen = ga_next_generation_trunc(population, population_size, percent_dead, percent_cross, crossover, mutations, mutate);
        else if (sel_strat == SEL_TOURNAMENT)
            /* Do tournaments to define which solutions are selected to cross.
            If the percentage dead is half or more, all individuals reproduce.
            The strongest solution stays in the population if it is not topped.*/
            gen = ga_next_generation_tournament(population, population_size, tournament_size, GA_MINIMIZE, fitness, crossover, mutations, mutate);
    }
    
    free(population);
    free(chromosome_chunk);
    tsp_2d_free(tsp);

    return 0;
}
