#ifndef _PIPE_H
#define _PIPE_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#ifdef PARALLEL
#include "mpi.h"
#endif


/* Constant Defines used throughout code */
#define MAX_LINE_LEN 24277
#define MAX_PROTEIN_LEN 24177
#define MAX_NUM_PROTEINS 151824
#define MAX_NEIGHBOURS 3000
#define W 20       /* Length of segment */
#define SCORE 35  /* PAM Score threshold */

#define DB_DIR "pipe_db"

/*
 * Token_Pair type stores two char arrays.  They are the '\t' separated
 * strings on an input line
 */
typedef struct {
  char *col1;
  char *col2;
} Token_Pair;

/********* Data Structures **********/
/*
 * Adj_list type is defined to be a struct, storing:
 *   id   - index of a Protein_Node array element
 *          Made this an int to allow for > 64k proteins
 *   next - next element in adjacency list.  I guess NULL would mean
 *          end of list?
 */
typedef struct Adj_list_STRUCT {
  unsigned int id;
  struct Adj_list_STRUCT *next;
} Adj_list;

/*
 * Protein_node type:
 *   adj      - start of the adjacency list for this protein
 *   yname    - the "yname" string which identifies this protein
 *   sequence - the amino acid character sequence
 */
typedef struct {
  Adj_list *adj;
  char *yname;
  char *sequence;
} Protein_Node;


/* Function Prototypes */
void init_mapback();
int count_lines(char *);
int return_organism(int prot_id, const int ranges[], unsigned int num_org);
int proper_organism(int prot_id, const int ranges[], unsigned int num_org, const int valid[]);
int read_line(FILE *, Token_Pair *);
int compare_segment(char *, char *);
int min(int, int);
void slide(char *, char *, int, int *, int, int);
double timer();


#ifdef __MAIN__

#ifdef PARALLEL
#define timer MPI_Wtime
#else
double timer() {
	static struct timeval tv;
	static struct timezone tz;

	gettimeofday(&tv, &tz);
	return (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
}
#endif


char MAPBACK[256];

/* PAM120 Substitution matrix */
int PAM120[23][23] = {
      /*A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X */
/*A*/ { 3,-3,-1, 0,-3,-1, 0, 1,-3,-1,-3,-2,-2,-4, 1, 1, 1,-7,-4, 0, 0,-1,-1},
/*R*/ {-3, 6,-1,-3,-4, 1,-3,-4, 1,-2,-4, 2,-1,-5,-1,-1,-2, 1,-5,-3,-2,-1,-2},
/*N*/ {-1,-1, 4, 2,-5, 0, 1, 0, 2,-2,-4, 1,-3,-4,-2, 1, 0,-4,-2,-3, 3, 0,-1},
/*D*/ { 0,-3, 2, 5,-7, 1, 3, 0, 0,-3,-5,-1,-4,-7,-3, 0,-1,-8,-5,-3, 4, 3,-2},
/*C*/ {-3,-4,-5,-7, 9,-7,-7,-4,-4,-3,-7,-7,-6,-6,-4, 0,-3,-8,-1,-3,-6,-7,-4},
/*Q*/ {-1, 1, 0, 1,-7, 6, 2,-3, 3,-3,-2, 0,-1,-6, 0,-2,-2,-6,-5,-3, 0, 4,-1},
/*E*/ { 0,-3, 1, 3,-7, 2, 5,-1,-1,-3,-4,-1,-3,-7,-2,-1,-2,-8,-5,-3, 3, 4,-1},
/*G*/ { 1,-4, 0, 0,-4,-3,-1, 5,-4,-4,-5,-3,-4,-5,-2, 1,-1,-8,-6,-2, 0,-2,-2},
/*H*/ {-3, 1, 2, 0,-4, 3,-1,-4, 7,-4,-3,-2,-4,-3,-1,-2,-3,-3,-1,-3, 1, 1,-2},
/*I*/ {-1,-2,-2,-3,-3,-3,-3,-4,-4, 6, 1,-3, 1, 0,-3,-2, 0,-6,-2, 3,-3,-3,-1},
/*L*/ {-3,-4,-4,-5,-7,-2,-4,-5,-3, 1, 5,-4, 3, 0,-3,-4,-3,-3,-2, 1,-4,-3,-2},
/*K*/ {-2, 2, 1,-1,-7, 0,-1,-3,-2,-3,-4, 5, 0,-7,-2,-1,-1,-5,-5,-4, 0,-1,-2},
/*M*/ {-2,-1,-3,-4,-6,-1,-3,-4,-4, 1, 3, 0, 8,-1,-3,-2,-1,-6,-4, 1,-4,-2,-2},
/*F*/ {-4,-5,-4,-7,-6,-6,-7,-5,-3, 0, 0,-7,-1, 8,-5,-3,-4,-1, 4,-3,-5,-6,-3},
/*P*/ { 1,-1,-2,-3,-4, 0,-2,-2,-1,-3,-3,-2,-3,-5, 6, 1,-1,-7,-6,-2,-2,-1,-2},
/*S*/ { 1,-1, 1, 0, 0,-2,-1, 1,-2,-2,-4,-1,-2,-3, 1, 3, 2,-2,-3,-2, 0,-1,-1},
/*T*/ { 1,-2, 0,-1,-3,-2,-2,-1,-3, 0,-3,-1,-1,-4,-1, 2, 4,-6,-3, 0, 0,-2,-1},
/*W*/ {-7, 1,-4,-8,-8,-6,-8,-8,-3,-6,-3,-5,-6,-1,-7,-2,-6,12,-2,-8,-6,-7,-5},
/*Y*/ {-4,-5,-2,-5,-1,-5,-5,-6,-1,-2,-2,-5,-4, 4,-6,-3,-3,-2, 8,-3,-3,-5,-3},
/*V*/ { 0,-3,-3,-3,-3,-3,-3,-2,-3, 3, 1,-4, 1,-3,-2,-2, 0,-8,-3, 5,-3,-3,-1},
/*B*/ { 0,-2, 3, 4,-6, 0, 3, 0, 1,-3,-4, 0,-4,-5,-2, 0, 0,-6,-3,-3, 4, 2,-1},
/*Z*/ {-1,-1, 0, 3,-7, 4, 4,-2, 1,-3,-3,-1,-2,-6,-1,-1,-2,-7,-5,-3, 2, 4,-1},
/*X*/ {-1,-2,-1,-2,-4,-1,-1,-2,-2,-1,-2,-2,-2,-3,-2,-1,-1,-5,-3,-1,-1,-1,-2}};

#else
extern char MAPBACK[256];
extern int PAM120[23][23];
#endif


#endif /* _PIPE_H */
