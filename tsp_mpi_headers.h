#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include <time.h>

//Defines the number of cities to be used in the problem
#define N_OF_CS 14
#define WORK	0
#define RESULT	1
#define DIE		2
#define GRAIN	6

typedef struct {
    char *name;
    double x;
    double y;
}City;


void print_cities(City cities[]);
void print_int_arr(int arr[]);
double distance(City *c1, City *c2);
void copy_path(int path1[], int path2[]);
void fill_distance_m(double distance_m[N_OF_CS][N_OF_CS], City cities[]);
double calc_length(int path[], double distance_m[N_OF_CS][N_OF_CS]);
void print_distance_m(double *distance_m, int n);
void print_distance_m2(double distance_m[N_OF_CS][N_OF_CS]);

void tsp(double distance_m[N_OF_CS][N_OF_CS], int start);
void tsp_aux(int path[], int path_size, int available[],
             double distance_m[N_OF_CS-1][N_OF_CS],
             int best_path[], double *best_length);

