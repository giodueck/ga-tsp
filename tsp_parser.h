#pragma once

#include <stddef.h>

/* Basic (incomplete) parser for TSPLIB (.tsp) files
    http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf
*/

typedef struct {
    double x, y;
} tsp_2d_node_t;

typedef struct {
    size_t dim;
    tsp_2d_node_t *nodes;
} tsp_2d_t;

tsp_2d_t tsp_2d_read(const char *filename);

tsp_2d_t tsp_2d_read_dedup(const char *filename);

void tsp_2d_free(tsp_2d_t tsp);

