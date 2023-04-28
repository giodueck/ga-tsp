#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "genetic.h"
#include "tsp_parser.h"

#define POPSIZE 1000

tsp_2d_t tsp;
int percent_elite = 5, percent_dead = 60, percent_cross = 0;
int mutations = 10000;
int gen = 0, max_gens = 1000;

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

double dist(tsp_2d_node_t a, tsp_2d_node_t b)
{
    double n = a.x - b.x;
    n *= n;
    double m = a.y - b.y;
    m *= m;
    return sqrt(n + m);
}

unsigned int fitness(ga_solution_t *sol)
{
    double d = 0;
    for (int i = 0; i < sol->chrom_len; i++)
    {
        int j = (i + 1) % sol->chrom_len;
        d += dist(tsp.nodes[((uint32_t *) sol->chromosome)[i]],
                  tsp.nodes[((uint32_t *) sol->chromosome)[j]]);
    }
    return (unsigned int) d;
}

void crossover(ga_solution_t *p1, ga_solution_t *p2, ga_solution_t *child);

void mutate(ga_solution_t *sol, int per_Mi)
{
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

    uint32_t *chromosome_chunk = (uint32_t *) malloc(sizeof(uint32_t) * tsp.dim * POPSIZE);
    ga_solution_t population[POPSIZE] = {0};

    /* Initialize population */
    ga_init(population, POPSIZE, tsp.dim, sizeof(uint32_t), chromosome_chunk, generate_tsp_solution);
    
    /* Evolve for max_gens number of generations */
    for (int i = 0; i < max_gens; i++)
    {
        ga_eval(population, POPSIZE, fitness);
        ga_select(population, POPSIZE, GA_MINIMIZE, percent_dead, percent_elite);

        // for (int i = 0; i < POPSIZE && !population[i].dead; i++)
        // {
        //     printf("%s%d: %lu\n", population[i].elite ? "*" : " ", i, population[i].fitness);
        // }

        if ((gen + 1) % 10 == 0)
        {
            int64_t best, worst_elite = 0, avg, worst;
            ga_gen_info(population, POPSIZE, percent_elite, &best, &worst_elite, &avg, &worst);
            printf("%3d:\tB: %5lu\t%3d%%: %5lu\tA: %5lu\tW: %5lu\n", gen + 1, best, percent_elite, worst_elite, avg, worst);
        }

        gen = ga_next_generation(population, POPSIZE, percent_dead, percent_cross, NULL /* crossover */, mutations, mutate);
    }
    
    free(chromosome_chunk);
    tsp_2d_free(tsp);

    return 0;
    
}
