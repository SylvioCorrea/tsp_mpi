#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "mpi.h"

//Defines the number of cities to be used in the problem
#define N_OF_CS 14

//These define tags for mpi communication
#define WORK	1
#define RESULT	2
#define DIE	3

//Defines the amount of cities a slave will be required to perform permutations on
#define GRAIN	10

typedef struct {
    char *name;
    double x;
    double y;
}City;

/* Used for process comunication.
   Master sends a partially filled city path that will
   be used by a slave as basis for possible permutations
   with the remaining cities. best_length holds the best
   distance received by the master until now.
   A slave will use the struct to send back the best
   path found along with its length.
   */
typedef struct {
    int path[N_OF_CS];
    double best_length;
}Message;


void print_cities(City cities[]);
void print_int_arr(int arr[]);
double distance(City *c1, City *c2);
void copy_path(int path1[], int path2[]);
void fill_distance_m(double distance_m[N_OF_CS][N_OF_CS], City cities[]);
double calc_length(int path[], double distance_m[N_OF_CS][N_OF_CS]);
void print_distance_m(double *distance_m, int n);
void print_distance_m2(double distance_m[N_OF_CS][N_OF_CS]);

void master_routine(Message *message_ptr, int available[],
                    int best_path[], int path_size, int *burst);
void tsp_aux(int path[], int path_size, int available[],
             double distance_m[N_OF_CS-1][N_OF_CS],
             int best_path[], double *best_length);

