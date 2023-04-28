#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "genetic.h"
#include "tsp_parser.h"

tsp_2d_t tsp;

/* Parameters */
int population_size = 1500;
int percent_elite = 5, percent_dead = 50, percent_cross = 50;
int mutations = 1000; // mutations / (1024*1024) = mutation chance
int max_gens = 2000, gen_info_interval = 50;

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
}

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

    uint32_t *chromosome_chunk = (uint32_t *) malloc(sizeof(uint32_t) * tsp.dim * population_size);
    ga_solution_t *population = (ga_solution_t *) malloc(sizeof(ga_solution_t) * population_size);

    /* Initialize population */
    ga_init(population, population_size, tsp.dim, sizeof(uint32_t), chromosome_chunk, generate_tsp_solution);

    int gen = 0;
    
    /* Evolve for max_gens number of generations */
    for (int i = 0; i < max_gens; i++)
    {
        ga_eval(population, population_size, fitness);
        ga_select_trunc(population, population_size, GA_MINIMIZE, percent_dead, percent_elite);
        
        if ((gen + 1) % gen_info_interval == 0 || !gen)
        {
            int64_t best, worst_elite = 0, avg, worst;
            ga_gen_info(population, population_size, percent_elite, &best, &worst_elite, &avg, &worst);
            printf("%3d:\tB: %5lu\t%3d%%: %5lu\tA: %5lu\tW: %5lu\n", gen + 1, best, percent_elite, worst_elite, avg, worst);
        }

        gen = ga_next_generation(population, population_size, percent_dead, percent_cross, crossover, mutations, mutate);
    }
    
    free(population);
    free(chromosome_chunk);
    tsp_2d_free(tsp);

    return 0;
    
}
