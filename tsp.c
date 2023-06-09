#include "tsp.h"
#include "tsp_parser.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern tsp_2d_t tsp;
extern int mutations;
extern struct drand48_data *rbufs;

// Initializes a random solution
void generate_tsp_solution(ga_solution_t *sol, size_t i, size_t chrom_len, void *chrom_chunk, uint8_t *marks)
{
    uint32_t *chromosome = (uint32_t *) chrom_chunk + i * chrom_len;
    // uint8_t *marks = (uint8_t *) malloc(sizeof(uint8_t) * chrom_len);
    memset(marks, 0, sizeof(uint8_t) * chrom_len);

    for (size_t j = chrom_len; j; j--)
    {
        long n;
        lrand48_r(&rbufs[0], &n);
        size_t r = n % j;
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
}

// Euclidean distance
double dist(const tsp_2d_node_t a, const tsp_2d_node_t b)
{
    double n = a.x - b.x;
    double m = a.y - b.y;
    return sqrt(n*n + m*m);
}

// Distance based fitness
int64_t fitness(ga_solution_t *sol)
{
    int64_t d = 0;
    if (sol->fit_gen)
        return sol->fitness;
    for (int i = 0; i < sol->chrom_len; i++)
    {
        int j = (i + 1) % sol->chrom_len;
        d += round(dist(tsp.nodes[((uint32_t *) sol->chromosome)[i]],
                        tsp.nodes[((uint32_t *) sol->chromosome)[j]]));
    }
    sol->fitness = d;
    sol->fit_gen = 1;
    return d;
}

// Cross two solutions and produce a child solution with traits from both parents 
void crossover(ga_solution_t *p1, ga_solution_t *p2, ga_solution_t *child, uint8_t *marks, struct drand48_data *rbuf)
{
    // Take half of the chromosome of one parent, then the remaining half of the other such that
    // nodes don't repeat
    
    long lrand;
    lrand48_r(rbuf, &lrand);
    int start = lrand % (p1->chrom_len / 2);
    int l = p1->chrom_len / 2;

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

    // If parents are less than 5% different
    if (diff <= p1->chrom_len / 20)
        mutate(child, mutations * 20, rbuf); // 15 times as likely to have mutations
    // ^ This is not what happens in real life, but it gives better results in this case
}

// Apply random swaps of genes dictated by some small chance
void mutate(ga_solution_t *sol, int per_Mi, struct drand48_data *rbuf)
{
    // 1024*1024 - 1
    // This is close enough to 1 million and a good mask to efficiently get small rand() numbers
    per_Mi &= 0xFFFFF;
    int per_Mi2 = 3 * per_Mi / 4 + 1;
    long n, n2;
    lrand48_r(rbuf, &n);
    while ((n & 0xFFFFF) < per_Mi)
    {
        n2 = n;
        lrand48_r(rbuf, &n);
        uint32_t i = n % sol->chrom_len;
        uint32_t aux = ((uint32_t*)sol->chromosome)[i];
        lrand48_r(rbuf, &n);

        uint32_t j;
        // Most times make j a neighbor
        if ((n2 & 0xFFFFF) < per_Mi2)
            j = (i + 1) % sol->chrom_len;
        else 
            j = n % sol->chrom_len;

        // Sometimes do 2-swap
        if ((n & 0xF) < 0xA)
        {
            ((uint32_t*)sol->chromosome)[i] = ((uint32_t*)sol->chromosome)[j];
            ((uint32_t*)sol->chromosome)[j] = aux;
        } else // other times to 3-swap
        {
            lrand48_r(rbuf, &n);
            uint32_t k = n % sol->chrom_len;
            ((uint32_t*)sol->chromosome)[i] = ((uint32_t*)sol->chromosome)[j];
            ((uint32_t*)sol->chromosome)[j] = ((uint32_t*)sol->chromosome)[k];
            ((uint32_t*)sol->chromosome)[k] = aux;
        }
    }
}

void verify_tsp_solutions(ga_solution_t *sol, size_t size, struct drand48_data *rbuf)
{
    uint8_t *marks = (uint8_t *) malloc(sizeof(uint8_t) * sol->chrom_len);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < sol->chrom_len; j++)
            marks[j] = 0;

        for (int j = 0; j < sol->chrom_len; j++)
        {
            if (marks[((uint32_t *) sol[i].chromosome)[j]])
            {
                generate_tsp_solution(&sol[i], 1, sol->chrom_len, sol[i].chromosome, marks);
                break;
            }
            marks[((uint32_t *) sol[i].chromosome)[j]] = 1;
        }
    }
    free(marks);
}

