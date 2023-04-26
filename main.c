#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "genetic.h"
#include "tsp_parser.h"

void generate_tsp_solution(ga_solution_t *sol, size_t i, size_t chrom_len, void *chrom_chunk)
{
}

int fitness(ga_solution_t *sol);

void crossover(ga_solution_t *p1, ga_solution_t *p2, ga_solution_t *child);

void mutate(ga_solution_t *sol);

int main(int argc, char **argv)
{
    tsp_2d_t tsp = tsp_2d_read(argv[1]);

    for (size_t i = 0; i < tsp.dim; i++)
    {
        printf("%lu: %lf %lf\n", i+1, tsp.nodes[i].x, tsp.nodes[i].y);
    }

    tsp_2d_free(tsp);

    return 0;
}
