#include "tsp_parser.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*`
typedef struct tsp_2d_t {
    size_t dim;
    tsp_2d_node_t *nodes;
} tsp_2d_t;
*/

#define DELIM " :\n"

tsp_2d_t tsp_2d_read(const char *filename)
{
    FILE *fd = fopen(filename, "rt");
    tsp_2d_t tsp = {0};
    char buf[BUFSIZ], *tok;
    int c, i = 0;
    size_t dimension = 0, index;

    c = getc(fd);
    do
    {
        // Read line
        if (c != '\n')
        {
            buf[i++] = c;
            continue;
        }
        else 
        {
            buf[i] = '\0';
            i = 0;
        }

        // Interpret keyword
        tok = strtok(buf, DELIM);
        if (strcmp(tok, "DIMENSION") == 0)
        {
            tok = strtok(NULL, DELIM);
            dimension = atoi(tok);
            tsp.nodes = malloc(sizeof(tsp_2d_node_t) * dimension);
            tsp.dim = dimension;
            continue;
        }
        
        // Read coordinate
        if ((index = atoi(tok)) != 0)
        {
            if (!dimension)
            {
                fprintf(stderr, "TSPLIB format error: coordinates read before dimension\n");
                exit(1);
            }

            tok = strtok(NULL, " \n");
            tsp.nodes[index - 1].x = atof(tok);
            tok = strtok(NULL, " \n");
            tsp.nodes[index - 1].y = atof(tok);
        }
    } while ((c = getc(fd)) != EOF);

    fclose(fd);
    return tsp;
}

tsp_2d_t tsp_2d_read_dedup(const char *filename)
{
    tsp_2d_t otsp = tsp_2d_read(filename);
    char *marks = (char *) malloc(sizeof(char) * otsp.dim);
    int cnt = otsp.dim;

    for (int i = 1; i < otsp.dim; i++)
    {
        marks[i] = 0;
        for (int j = 0; j < i; j++)
        {
            if (otsp.nodes[i].x == otsp.nodes[j].x && otsp.nodes[i].y == otsp.nodes[j].y)
            {
                marks[i] = 1;
                cnt--;
                j = i;
            }
        }
    }

    tsp_2d_t ntsp = {0};
    ntsp.nodes = malloc(sizeof(tsp_2d_node_t) * cnt);
    ntsp.dim = cnt;

    for (int i = 0, j = 0; i < otsp.dim; i++)
    {
        if (!marks[i])
            ntsp.nodes[j++] = otsp.nodes[i];
    }

    tsp_2d_free(otsp);
    return ntsp;
}

void tsp_2d_free(tsp_2d_t tsp)
{
    free(tsp.nodes);
}
