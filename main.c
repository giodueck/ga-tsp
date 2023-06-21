#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#else
 #ifdef MPI
#include <mpi.h>
 #else
#include <pthread.h>
 #endif
#endif

#include "genetic.h"
#include "tsp_parser.h"
#include "tsp.h"

#define SEL_TRUNCATE   0
#define SEL_TOURNAMENT 1

tsp_2d_t tsp = {0};
FILE *csv = NULL;

#ifndef _OPENMP
#ifndef MPI
pthread_t *threads = NULL;
#endif
#endif

int *thread_bounds = NULL;
struct drand48_data *rbufs = NULL;

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
    -f      TSP file, exclude duplications
    -g      generations
    -h      print help
    -i      gen. info interval
    -k      tournament size
    -l      TSP file, keep duplications
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
    -e [0-100]      Affects display of generation statistics, shows fitness of\n\
                    top percentage of solutions.\n\
                        Default: 5\n\n\
    -f [filename]   Load TSP from the given file. Must be TSPLIB format.\n\
                    Will exclude duplicates.\n\n\
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
                    offspring, so they are held in pairs.\n\
                        Default: 4\n\n\
    -l [filename]   Load TSP from the given file. Must be TSPLIB format. Unlike -f\n\
                    it will keep all duplicates. Can be used implicitly.\n\n\
    -m [integer]    Mutation rate out of 0x0FFFFF, or 1024x1024-1.\n\
                    Default: 1000 (~0.1%)\n\n\
    -o [filename]   Output generation info to a CSV file.\n\n\
    -p [integer]    Total population size. If there are more than one island this\n\
                    population is divided evenly among them.\n\
                        Default: 2500\n\n\
    -r [integer]    Supply a seed to the random number generator.\n\
                    Default: 1\n\n\
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
    const char *optstring = "ae:f:g:hi:k:l:m:o:p:r:t:u:";
    int opt = 0;

    while ((opt = getopt(argc, argv, optstring)) != -1)
    {
        switch (opt)
        {
            case 'a':
                f_answer = 1;
                break;
            case 'e':
                percent_elite = atoi(optarg);
                break;
            case 'f':
                tsp = tsp_2d_read_dedup(optarg);
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
            case 'l':
                tsp = tsp_2d_read(optarg);
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
        if (sel_strat == SEL_TOURNAMENT)
            /* Do tournaments to define which solutions are selected to cross.
            If the percentage dead is half or more, all individuals reproduce.
            The strongest solution stays in the population if it is not topped.*/
            gen = ga_next_generation_tournament(population, population_size, tournament_size, GA_MINIMIZE, fitness, crossover, mutations, mutate, &rbufs[0]);
    }

    return gen;
}

#ifdef MPI
int mpi_ga(ga_solution_t *population, int gens, int island_size)
{
    int gen = population->generation;
    while (gens-- > 0)
    {
        if (sel_strat == SEL_TOURNAMENT)
            /* Do tournaments to define which solutions are selected to cross.
            If the percentage dead is half or more, all individuals reproduce.
            The strongest solution stays in the population if it is not topped.*/
            gen = ga_next_generation_tournament(population, island_size, tournament_size, GA_MINIMIZE, fitness, crossover, mutations, mutate, &rbufs[0]);
    }

    return gen;
}
#endif

struct parallel_ga_arg {
    ga_solution_t *population;
    int gens, low, high, t;
};

// Executes GA in parallel for a chunk of population, defined by the indices in the range [low, high)
void *parallel_ga(void *_arg)
{
    struct parallel_ga_arg arg = *(struct parallel_ga_arg *) _arg;
    // struct drand48_data rd;
    // srand48_r(arg.population->generation + arg.low, &rd);
    while (arg.gens-- > 0)
    {
        if (sel_strat == SEL_TOURNAMENT)
            /* Do tournaments to define which solutions are selected to cross.
            If the percentage dead is half or more, all individuals reproduce.
            The strongest solution stays in the population if it is not topped.*/
            ga_next_generation_tournament(arg.population + arg.low, arg.high - arg.low, tournament_size, GA_MINIMIZE, fitness, crossover, mutations, mutate, &rbufs[arg.t]);
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
        if (pop->generation != max_gens)
            ga_gen_info_unsorted(pop, population_size, percent_elite, &best, &avg, &worst);
        else 
            ga_gen_info(pop, population_size, percent_elite, &best, &worst_elite, &avg, &worst);
    }
    else
    {
        gen = (pop + thread_bounds[island])->generation;
        ga_select_trunc(pop + thread_bounds[island], thread_bounds[island + 1] - thread_bounds[island], GA_MINIMIZE, percent_dead, percent_elite, fitness);
        if (pop->generation != max_gens)
            ga_gen_info_unsorted(pop + thread_bounds[island], thread_bounds[island + 1] - thread_bounds[island], percent_elite, &best, &avg, &worst);
        else
            ga_gen_info(pop + thread_bounds[island], thread_bounds[island + 1] - thread_bounds[island], percent_elite, &best, &worst_elite, &avg, &worst);
    }

    if (csv)
        fprintf(csv, "%d,%d,%lu,%d,%lu,%lu,%lu\n", island, gen, best, percent_elite, worst_elite, avg, worst);
    if (num_threads > 1)
        printf("I: %3d\tG: %6d:\tB: %5lu\t%3d%%: %5lu\tA: %5lu\tW: %5lu\n", island, gen, best, percent_elite, worst_elite, avg, worst);
    else 
        printf("G: %6d:\tB: %5lu\t%3d%%: %5lu\tA: %5lu\tW: %5lu\n", gen, best, percent_elite, worst_elite, avg, worst);
}

#ifdef MPI
#define FLAG_CONT 1
#define FLAG_TERM 0

#define FLAG_TAG  1
#define DATA_TAG  2

// Send island population's genetic information to the process with ID dest_proc as an array of chars.
// A flag char is sent first: 0 when pop is NULL, signifies the end of the program, >0 otherwise, pop is sent afterwards.
// A 0 flag should only be sent from the master, as slaves have no knowledge of how many generations have passed.
void send_island(int dest_proc, ga_solution_t *pop, int from, int up_to)
{
    char flag[1] = {FLAG_TERM};
    if (pop == NULL)
    {
        MPI_Send(flag, 1, MPI_CHAR, dest_proc, FLAG_TAG, MPI_COMM_WORLD);
        return;
    }
    flag[0] = FLAG_CONT;
    MPI_Send(flag, 1, MPI_CHAR, dest_proc, FLAG_TAG, MPI_COMM_WORLD);

    int proc_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    for (int i = from; i < up_to; i++)
    {
        MPI_Send(((uint32_t*)pop[i].chromosome), pop->chrom_len, MPI_UINT32_T, dest_proc, DATA_TAG, MPI_COMM_WORLD);
    }
}

// Receive an island population's genetic information as an array of chars from the process with ID src_proc and
// store it in the dest array.
// A flag char is received first: 0 signifies the end of the program, >0 signifies that the population is sent next.
// If a 0 flag is received, the slave will terminate execution
int receive_island(int src_proc, ga_solution_t *dest, int from, int up_to)
{
    char flag[1] = {FLAG_TERM};
    MPI_Recv(flag, 1, MPI_CHAR, src_proc, FLAG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (flag[0] == FLAG_TERM)
        return FLAG_TERM;

    int proc_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    for (int i = from; i < up_to; i++)
    {
        MPI_Recv(((uint32_t*)dest[i].chromosome), dest->chrom_len, MPI_UINT32_T, src_proc, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    return FLAG_CONT;
}

// Slave main function, loops until receive_island returns FLAG_TERM
void slave_main(int proc_id, int from, int up_to, int gens)
{
    int island_size = up_to - from;
    uint32_t *chromosome_chunk = (uint32_t *) malloc(sizeof(uint32_t) * tsp.dim * island_size);
    ga_solution_t *pop = (ga_solution_t *) malloc(sizeof(ga_solution_t) * island_size);

    ga_init(pop, island_size, tsp.dim, sizeof(uint32_t), chromosome_chunk, generate_tsp_solution);

    printf("Process %d in slave_main, island_size = %d, from %d up to %d\n", proc_id, island_size, from, up_to);
    while (receive_island(0, pop, 0, island_size) != FLAG_TERM)
    {
        // printf("slave_main:%d: FLAG_CONT\n", proc_id);
        // Evolve
        mpi_ga(pop, gens, island_size);
        // Send back to master
        send_island(0, pop, 0, island_size);
    }
    // printf("slave_main:%d: FLAG_TERM\n", proc_id);

    free(pop);
    free(chromosome_chunk);
}
#endif //MPI

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    #ifdef _OPENMP
    printf("Using OpenMP!\n\n");
    #endif
    #ifdef MPI

    int num_procs, proc_id;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    if (proc_id == 0)
        printf("Using OpenMPI!\n\n");
    if (proc_id > num_threads)
    {
        printf("Note: Too few islands, node %d idle.\n", proc_id);
        MPI_Finalize();
        return 0;
    } 
    if (num_threads >= num_procs)
    {
        printf("Error: Too few nodes, for N islands need N+1 nodes.\n");
        exit(EXIT_FAILURE);
    }
    #endif

    if (!tsp.dim)
    {
        if (optind >= argc)
        {
            #ifdef MPI
            if (proc_id == 0)
            #endif
            fprintf(stderr, "Usage: '%s [options] <file.tsp>'\nSee '%s -h' for help\n", argv[0], argv[0]);
            exit(EXIT_FAILURE);
        }

        // There is a bug in the deduplicated reading causing wrong fitness results.
        tsp = tsp_2d_read(argv[optind]);
    }

    uint32_t *chromosome_chunk;
    ga_solution_t *population;
    
    #ifdef MPI
    if (proc_id == 0) {
    #endif

    printf("Dim = %lu\n", tsp.dim);
    chromosome_chunk = (uint32_t *) malloc(sizeof(uint32_t) * tsp.dim * population_size);
    population = (ga_solution_t *) malloc(sizeof(ga_solution_t) * population_size);

    if (csv)
        fprintf(csv, "Island,Generation,Best,Elite%%,Elite,Average,Worst\n");

    #ifdef MPI
    }
    #endif

    if (num_threads > 1)
    {
        #ifndef _OPENMP
        #ifndef MPI
        threads = (pthread_t *) malloc(sizeof(pthread_t) * num_threads);
        #endif
        #endif
        thread_bounds = (int *) malloc(sizeof(int) * (num_threads + 1));

        int low = 0;
        thread_bounds[0] = low;
        for (int i = 0; i < num_threads - 1; i++)
        {
            low += population_size / num_threads;
            thread_bounds[i + 1] = low;
        }
        thread_bounds[num_threads] = population_size;
    } else 
    {
        num_threads = 1;
        thread_bounds = (int *) malloc(sizeof(int) * (num_threads + 1));
        thread_bounds[0] = 0;
        thread_bounds[1] = population_size;
    }

    #ifndef MPI
    rbufs = (struct drand48_data *) malloc(sizeof(struct drand48_data) * num_threads);
    for (int i = 0; i < num_threads; i++)
        srand48_r(rand(), &rbufs[i]);
    #else
    rbufs = (struct drand48_data *) malloc(sizeof(struct drand48_data));
    srand48_r(rand() + proc_id, rbufs);

    if (proc_id > 0)
    {
        int gens = max_gens;
        if (island_cross_interval > 0)
            gens = island_cross_interval;
        if (num_threads > 1)
            // While thread bounds for i = 0 would be for the first thread, the first thread here is i = 1
            slave_main(proc_id, thread_bounds[proc_id - 1], thread_bounds[proc_id], gens - 1);
        else
            slave_main(proc_id, 0, population_size, gens - 1);
        MPI_Finalize();
        free(rbufs);
        tsp_2d_free(tsp);

        return 0;
    }
    #endif

    /* Initialize population */
    ga_init(population, population_size, tsp.dim, sizeof(uint32_t), chromosome_chunk, generate_tsp_solution);

    int gen = 0;

    // DEBUG
    #ifdef DEBUG
    printf("PID %d ready for attach\n", getpid());
    int debug_temp = 0;
    while (!debug_temp)
        sleep(5);
    #endif

    /* Evolve for max_gens number of generations */
    while (gen < max_gens)
    {
        #ifndef MPI
        /* Single-threaded */
        if (num_threads <= 1)
        {
            if (gen_info_interval > 0)
            {
                gen_info(population, 0);

                gen = serial_ga(population, (max_gens - gen - gen_info_interval >= 0) ? gen_info_interval : max_gens - gen);
            }
            else 
                gen = serial_ga(population, max_gens);
            continue;
        }
        #endif

        /* Multi-threaded */
        struct parallel_ga_arg *args = (struct parallel_ga_arg *) malloc(sizeof(struct parallel_ga_arg) * num_threads);
        if (island_cross_interval <= 0)
        {
            #ifdef _OPENMP
            #pragma omp parallel for
            for (int i = 0; i < num_threads; i++)
            {
                if (gen_info_interval > 0)
                    gen_info(population, i);
                args[i] = (struct parallel_ga_arg) { .population = population, .gens = max_gens, .low = thread_bounds[i], .high = thread_bounds[i + 1] };
                parallel_ga(&args[i]);
            }
            #else
            #ifdef MPI
            
            // Send islands from population array
            for (int i = 0; i < num_threads; i++)
            {
                // printf("Sending to %d, island_size = %d, from %d up to %d\n", i + 1, thread_bounds[i+1] - thread_bounds[i], thread_bounds[i], thread_bounds[i+1]);
                gen_info(population, i);
                send_island(i + 1, population, thread_bounds[i], thread_bounds[i + 1]);
            }

            // Receive islands, put back into population array
            for (int i = 0; i < num_threads; i++)
            {
                receive_island(i + 1, population, thread_bounds[i], thread_bounds[i + 1]);
            }

            #else
            for (int i = 0; i < num_threads; i++)
            {
                if (gen_info_interval > 0)
                    gen_info(population, i);
                args[i] = (struct parallel_ga_arg) { .population = population, .gens = max_gens, .low = thread_bounds[i], .high = thread_bounds[i + 1] };
                pthread_create(&threads[i], NULL, parallel_ga, &args[i]);
            }
            #endif
            #endif

            gen = max_gens;
        }
        else 
        {
            #ifdef _OPENMP
            #pragma omp parallel for
            for (int i = 0; i < num_threads; i++)
            {
                if (gen_info_interval > 0)
                    gen_info(population, i);
                args[i] = (struct parallel_ga_arg) { .population = population, .gens = ((max_gens - gen - island_cross_interval >= 0) ? island_cross_interval : max_gens - gen) - 1, .low = thread_bounds[i], .high = thread_bounds[i + 1], .t = i };
                parallel_ga(&args[i]);
            }
            #else
            #ifdef MPI
            
            // Send islands from population array
            for (int i = 0; i < num_threads; i++)
            {
                // printf("Sending to %d, island_size = %d, from %d up to %d\n", i + 1, thread_bounds[i+1] - thread_bounds[i], thread_bounds[i], thread_bounds[i+1]);
                gen_info(population, i);
                send_island(i + 1, population, thread_bounds[i], thread_bounds[i + 1]);
            }

            // Receive islands, put back into population array
            for (int i = 0; i < num_threads; i++)
            {
                receive_island(i + 1, population, thread_bounds[i], thread_bounds[i + 1]);
            }

            #else
            for (int i = 0; i < num_threads; i++)
            {
                if (gen_info_interval > 0)
                    gen_info(population, i);
                args[i] = (struct parallel_ga_arg) { .population = population, .gens = ((max_gens - gen - island_cross_interval >= 0) ? island_cross_interval : max_gens - gen) - 1, .low = thread_bounds[i], .high = thread_bounds[i + 1], .t = i };
                pthread_create(&threads[i], NULL, parallel_ga, &args[i]);
            }
            #endif
            #endif

            gen += island_cross_interval;
        }

        // Cross islands
        #ifndef _OPENMP
        #ifndef MPI
        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], NULL);
        }
        #endif
        #endif
        free(args);

        #ifdef MPI
        verify_tsp_solutions(population, population_size, rbufs);
        // TODO MPI code
        // printf("Master in main\n");
        for (int i = 0; i < population_size; i++)
            population[i].generation += (island_cross_interval <= 0) ? max_gens : island_cross_interval - 1;
        
        if (island_cross_interval > 0)
            gen = serial_ga(population, 1);
        #else

        if (island_cross_interval > 0)
            gen = serial_ga(population, 1);
        #endif
    }

    /* Print last generation */
    if (gen_info_interval >= 0)
    {
        if (num_threads <= 1)
            gen_info(population, 0);
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
    #ifndef _OPENMP
    #ifndef MPI
    if (threads)
        free(threads);
    #endif
    #endif
    free(thread_bounds);
    free(rbufs);

    #ifdef MPI
    for (int i = 0; i < num_threads; i++)
    {
        // Send termination signal
        send_island(i + 1, NULL, 0, 0);
    }
    MPI_Finalize();
    #endif
    return 0;
}
