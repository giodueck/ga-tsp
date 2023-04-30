#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "genetic.h"
#include "tsp_parser.h"
#include "tsp.h"

#define SEL_TRUNCATE   0
#define SEL_TOURNAMENT 1

tsp_2d_t tsp = {0};
FILE *csv = NULL;
pthread_t *threads = NULL;
int *thread_bounds = NULL;

/* Parameters */
int population_size = 2500;     // population size per thread
int max_gens = 3000;            // when to stop the algorithm
int gen_info_interval = 100;    // how often to print information about the population
int mutations = 1000;           // mutations / (1024*1024) = mutation chance
int num_threads = 1;            // number of islands to evolve
int island_cross_interval = 0;  // how often island populations are allowed to cross
int f_answer = 0;               // if 1 print shortest path found
int sel_strat = SEL_TOURNAMENT; // selection strategy 
    /* Truncation selection */
int percent_elite = 5;          // percentage of elite selection, makes gen_info more informative for tournament
int percent_dead = 50;          // how many solutions get replaced 
int percent_cross = 50;         // how many solutions are derived from crossover
    /* Tournament selection */
int tournament_size = 4;        // how many individuals get picked per tournament

/* CLI arguments 

    -a      print the shortest path found
    -c      cross percentage (trunc)
    -d      dead percentage (trunc)
    -e      elite percentage (trunc)
    -f      TSP file
    -g      generations
    -h      print help
    -i      gen. info interval
    -k      tournament size 
    -m      mutation rate
    -o      output gen info to file as CSV format
    -p      population size
    -r      PRNG seed
    -s      switch to truncation
    -t      island (thread) count
    -u      island crossover interval

    Of these only -h and -s don't take arguments
*/

void print_help(char **argv)
{
    const char *help_text = "\
Usage: %s [options] <file.tsp>\n\
\n\
  Options:\n\
    -a              Print the shortest path found after finishing evolution.\n\n\
    -c [0-100]      Percentage of new population generated from crossover.\n\
                    Only has an effect if -s is also given.\n\
                        Default: 50\n\n\
    -d [0-100]      Percentage of population not selected for generating new\n\
                    offspring. Only has an effect if -s is also given.\n\
                        Default: 50\n\n\
    -e [0-100]      Percentage of population covered by elitist selection. Only\n\
                    has an effect if -s is also given. Also affects display of\n\
                    generation statistics regardless of if -s is given.\n\
                        Default: 5\n\n\
    -f [filename]   Load TSP from the given file. Must be TSPLIB format. Can\n\
                    also be given without the -f option.\n\n\
    -g [integer]    Number of generations to evolve.\n\
                        Default: 3000\n\n\
    -h              Display this help.\n\n\
    -i [integer]    Number of generations between statistics prints. The\n\
                    population is sorted by fitness to find this information, which\n\
                    randomly affects tournament selection.\n\
                    If the number of islands is more than one, the information is\n\
                    printed after each crossing between islands instead.\n\
                    -1 to disable all output.\n\
                    0 to disable printing info before the algorithm finishes.\n\
                        Default: 100\n\n\
    -k [integer]    Number of individuals per tournament. Every tournament\n\
                    selects one parent and one individual to be replaced by\n\
                    offspring, so they are held in pairs. Only has an effect if\n\
                    -s is not given.\n\
                        Default: 4\n\n\
    -m [integer]    Mutation rate out of 0x0FFFFF, or 1024x1024-1.\n\
                    Default: 1000 (~0.1%)\n\n\
    -o [filename]   Output generation info to a CSV file.\n\n\
    -p [integer]    Total population size. If there are more than one island this\n\
                    population is divided evenly among them.\n\
                        Default: 2500\n\n\
    -r [integer]    Supply a seed to the random number generator.\n\
                    Default: 1\n\n\
    -s              Switch selection strategy to truncation with elitism selection.\n\
                        Default is tournament selection.\n\n\
    -t [integer]    Number of islands, each of which is handled by a thread.\n\
                        Default: 1\n\n\
    -u [integer]    Number of generations after which islands will have their\n\
                    populations crossed.\n\
                    If the interval is below 1, the populations will never cross.\n\
                        Default: 0\n\n";

    printf(help_text, argv[0]);
}

void parse_args(int argc, char **argv)
{
    const char *optstring = "ac:d:e:f:g:hi:k:m:o:p:r:st:u:";
    int opt = 0;

    while ((opt = getopt(argc, argv, optstring)) != -1)
    {
        switch (opt)
        {
            case 'a':
                f_answer = 1;
                break;
            case 'c':
                percent_cross = atoi(optarg);
                break;
            case 'd':
                percent_dead = atoi(optarg);
                break;
            case 'e':
                percent_elite = atoi(optarg);
                break;
            case 'f':
                tsp = tsp_2d_read(optarg);
                break;
            case 'g':
                max_gens = atoi(optarg);
                break;
            case 'h':
                print_help(argv);
                exit(EXIT_SUCCESS);
            case 'i':
                gen_info_interval = atoi(optarg);
                break;
            case 'k':
                tournament_size = atoi(optarg);
                break;
            case 'm':
                mutations = atoi(optarg);
                break;
            case 'o':
                csv = fopen(optarg, "wt");
                break;
            case 'p':
                population_size = atoi(optarg);
                break;
            case 'r':
                srand(atoi(optarg));
                break;
            case 's':
                sel_strat = SEL_TRUNCATE;
                break;
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'u':
                island_cross_interval = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Usage: '%s [options] <file.tsp>'\nSee '%s -h' for help\n", argv[0], argv[0]);
                exit(EXIT_FAILURE);
        }
    }
}

int serial_ga(ga_solution_t *population, int gens)
{
    int gen = population->generation;
    while (gens-- > 0)
    {
        if (sel_strat == SEL_TRUNCATE)
            /* Sort and select elite and survivors according to the truncation selection method */
            ga_select_trunc(population, population_size, GA_MINIMIZE, percent_dead, percent_elite, fitness);

        if (sel_strat == SEL_TRUNCATE)
            /* Replace dead population with new offspring from surviving individuals */
            gen = ga_next_generation_trunc(population, population_size, percent_dead, percent_cross, crossover, mutations, mutate);
        else if (sel_strat == SEL_TOURNAMENT)
            /* Do tournaments to define which solutions are selected to cross.
            If the percentage dead is half or more, all individuals reproduce.
            The strongest solution stays in the population if it is not topped.*/
            gen = ga_next_generation_tournament(population, population_size, tournament_size, GA_MINIMIZE, fitness, crossover, mutations, mutate);
    }

    return gen;
}

struct parallel_ga_arg {
    ga_solution_t *population;
    int gens, low, high;
};

// Executes GA in parallel for a chunk of population, defined by the indices in the range [low, high)
void *parallel_ga(void *_arg)
{
    struct parallel_ga_arg arg = *(struct parallel_ga_arg *) _arg;
    // struct drand48_data rd;
    // srand48_r(arg.population->generation + arg.low, &rd);
    while (arg.gens-- > 0)
    {
        if (sel_strat == SEL_TRUNCATE)
            /* Sort and select elite and survivors according to the truncation selection method */
            ga_select_trunc(arg.population + arg.low, arg.high - arg.low, GA_MINIMIZE, percent_dead, percent_elite, fitness);

        if (sel_strat == SEL_TRUNCATE)
            /* Replace dead population with new offspring from surviving individuals */
            ga_next_generation_trunc(arg.population + arg.low, arg.high - arg.low, percent_dead, percent_cross, crossover, mutations, mutate);
        else if (sel_strat == SEL_TOURNAMENT)
            /* Do tournaments to define which solutions are selected to cross.
            If the percentage dead is half or more, all individuals reproduce.
            The strongest solution stays in the population if it is not topped.*/
            ga_next_generation_tournament(arg.population + arg.low, arg.high - arg.low, tournament_size, GA_MINIMIZE, fitness, crossover, mutations, mutate);
    }

    return NULL;
}

void gen_info(ga_solution_t *pop, int island)
{
    int64_t best, worst_elite = 0, avg, worst;
    int gen = pop->generation;
    // Just to sort population, doesn't make any changes
    if (num_threads <= 1)
    {
        ga_select_trunc(pop, population_size, GA_MINIMIZE, percent_dead, percent_elite, fitness);
        ga_gen_info(pop, population_size, percent_elite, &best, &worst_elite, &avg, &worst);
    }
    else
    {
        gen = (pop + thread_bounds[island])->generation;
        ga_select_trunc(pop + thread_bounds[island], thread_bounds[island + 1] - thread_bounds[island], GA_MINIMIZE, percent_dead, percent_elite, fitness);
        ga_gen_info(pop + thread_bounds[island], thread_bounds[island + 1] - thread_bounds[island], percent_elite, &best, &worst_elite, &avg, &worst);
    }

    if (csv)
        fprintf(csv, "%d,%d,%lu,%d,%lu,%lu,%lu\n", island, gen, best, percent_elite, worst_elite, avg, worst);
    if (num_threads > 1)
        printf("I: %3d\tG: %6d:\tB: %5lu\t%3d%%: %5lu\tA: %5lu\tW: %5lu\n", island, gen, best, percent_elite, worst_elite, avg, worst);
    else 
        printf("G: %6d:\tB: %5lu\t%3d%%: %5lu\tA: %5lu\tW: %5lu\n", gen, best, percent_elite, worst_elite, avg, worst);
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    if (!tsp.dim)
    {
        if (optind >= argc)
        {
            fprintf(stderr, "Usage: '%s [options] <file.tsp>'\nSee '%s -h' for help\n", argv[0], argv[0]);
            exit(EXIT_FAILURE);
        }

        tsp = tsp_2d_read(argv[optind]);
    }

    uint32_t *chromosome_chunk = (uint32_t *) malloc(sizeof(uint32_t) * tsp.dim * population_size);
    ga_solution_t *population = (ga_solution_t *) malloc(sizeof(ga_solution_t) * population_size);

    /* Initialize population */
    ga_init(population, population_size, tsp.dim, sizeof(uint32_t), chromosome_chunk, generate_tsp_solution);

    int gen = 0;
    int island = 0;

    if (csv)
        fprintf(csv, "Island,Generation,Best,Elite%%,Elite,Average,Worst\n");

    if (num_threads > 1)
    {
        threads = (pthread_t *) malloc(sizeof(pthread_t) * num_threads);
        thread_bounds = (int *) malloc(sizeof(int) * (num_threads + 1));

        int low = 0;
        thread_bounds[0] = low;
        for (int i = 0; i < num_threads - 1; i++)
        {
            low += population_size / num_threads;
            thread_bounds[i + 1] = low;
        }
        thread_bounds[num_threads] = population_size;
    }

    /* Evolve for max_gens number of generations */
    while (gen < max_gens)
    {
        /* Single-threaded */
        if (num_threads <= 1)
        {
            if (gen_info_interval > 0)
            {
                gen_info(population, island);

                gen = serial_ga(population, (max_gens - gen - gen_info_interval >= 0) ? gen_info_interval : max_gens - gen);
            }
            else 
                gen = serial_ga(population, max_gens);
            continue;
        }

        /* Multi-threaded */
        struct parallel_ga_arg *args = (struct parallel_ga_arg *) malloc(sizeof(struct parallel_ga_arg) * num_threads);
        if (island_cross_interval <= 0)
        {
            for (int i = 0; i < num_threads; i++)
            {
                if (gen_info_interval > 0)
                    gen_info(population, i);
                args[i] = (struct parallel_ga_arg) { .population = population, .gens = max_gens, .low = thread_bounds[i], .high = thread_bounds[i + 1] };
                pthread_create(&threads[i], NULL, parallel_ga, &args[i]);
            }

            gen = max_gens;
        } else 
        {
            for (int i = 0; i < num_threads; i++)
            {
                if (gen_info_interval > 0)
                    gen_info(population, i);
                args[i] = (struct parallel_ga_arg) { .population = population, .gens = ((max_gens - gen - island_cross_interval >= 0) ? island_cross_interval : max_gens - gen) - 1, .low = thread_bounds[i], .high = thread_bounds[i + 1] };
                pthread_create(&threads[i], NULL, parallel_ga, &args[i]);
            }

            gen += island_cross_interval;
        }

        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], NULL);
        }
        free(args);

        if (island_cross_interval > 0)
            gen = serial_ga(population, 1);
    }

    /* Print last generation */
    if (gen_info_interval >= 0)
    {
        if (num_threads <= 1)
            gen_info(population, island);
        else 
        {
            if (gen_info_interval > 0)
            {
                for (int i = 0; i < num_threads; i++)
                    gen_info(population, i);
                printf("\n");
            }
            
            /* Print total stats */
            int aux = num_threads;
            num_threads = 0;
            gen_info(population, 0);
            num_threads = aux;
        }
    }

    /* Print best path */
    if (f_answer)
    {
        ga_select_trunc(population, population_size, GA_MINIMIZE, percent_dead, percent_elite, fitness);
        printf("\nBest path after %d generations: %lu\n", max_gens, population[0].fitness);
        for (int i = 0; i < tsp.dim; i++)
        {
            uint32_t n = ((uint32_t *)population[0].chromosome)[i];
            printf("%s%u ", (i) ? "-> " : "", n);
        }
        printf("\n");
    }
    
    free(population);
    free(chromosome_chunk);
    tsp_2d_free(tsp);
    if (csv)
        fclose(csv);
    if (threads)
    {
        free(threads);
        free(thread_bounds);
    }

    return 0;
}
