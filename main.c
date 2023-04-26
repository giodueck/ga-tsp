#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "genetic.h"
#include "tsp_parser.h"

#define POPSIZE 5000

tsp_2d_t tsp;

void generate_tsp_solution(ga_solution_t *sol, size_t i, size_t chrom_len, void *chrom_chunk)
{
    uint32_t *chromosome = (uint32_t *) chrom_chunk + i * chrom_len;
    uint8_t *marks = (uint8_t *) malloc(sizeof(uint8_t) * chrom_len);
    memset(marks, 0, sizeof(uint8_t) * chrom_len);

    for (size_t j = chrom_len; j; j--)
    {
        size_t r = rand() % j;
        size_t l = 0;
        while (r && !marks[l])
        {
            if (!marks[l]) r--;
            l++;
        }
        chromosome[j - 1] = l;
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

void mutate(ga_solution_t *sol);

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [file.tsp] <optional random seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    tsp = tsp_2d_read(argv[1]);

    if (argc >= 3)
    {
        printf("seed = %d\n", atoi(argv[2]));
        srand(atoi(argv[2]));
    }

    // for (size_t i = 0; i < tsp.dim; i++)
    // {
    //     printf("%lu: %lf %lf\n", i+1, tsp.nodes[i].x, tsp.nodes[i].y);
    // }

    uint32_t *chromosome_chunk = (uint32_t *) malloc(sizeof(uint32_t) * tsp.dim * POPSIZE);
    ga_solution_t population[POPSIZE] = {0};

    ga_init(population, POPSIZE, tsp.dim, chromosome_chunk, generate_tsp_solution);

    ga_eval(population, POPSIZE, fitness);
    ga_select(population, POPSIZE, GA_MINIMIZE, 98, 1);

    for (int i = 0; i < POPSIZE && !population[i].dead; i++)
    {
        printf("%s%d: %u\n", population[i].elite ? "*" : " ", i, population[i].fitness);
    }
    
    free(chromosome_chunk);
    tsp_2d_free(tsp);

    return 0;
    
}
