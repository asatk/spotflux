#include <stdio.h>
//#include <optarg.h>
//#include <unistd.h>

#include "constants.h"

void parseargs(int argc, char **argv) {

    (void) argc;
    (void) argv;

}


void store(double **grid, FILE *f) {
    int i, j;
    for ( i = 0 ; i < ntheta ; i++ ) {
        for ( j = 0 ; j < nphi ; j++ ) {
            fprintf(f, "%.4le ", grid[i][j]);
        }
    }
//    fprintf(f,"\n\n");
}
/**
typedef struct Args {
    int ntheta;
    int nphi;
    char bprof;
    char rartype;
    int nrar;
    char simtype;
} *args;

args parse(int argc, char **argv) {



   args a;

    a = (args) malloc(sizeof(struct Args));

    a->ntheta = 128;
    a->nphi = 256;
    a->nt = 10000;
    a->bprof = 1;
    a->rartype = 0;
    a->simtype = 1;

    ntheta = a->ntheta;
    nphi = a->nphi;
    nt = a->nt;
    bprof = a->bprof;
    rartype = a->rartype;
    simtype = a->simtype;
    return NULL;
    return a;
}
 */
