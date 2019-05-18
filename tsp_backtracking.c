#include "tsp_mpi_headers.h"


//Compile using options:
//  -lm to link math library
//  -fopenmp to link openMP library

//N_OF_CS and N_OF_THREADS are defined in the header


//This function solves the tsp problem. The caller can choose the starting city,
//though it should be noted the sollution is a hamiltonian cycle and as such
//the solution list should be considered a ring, that is, it can be used
//for any starting city.
void tsp(double distance_m[N_OF_CS][N_OF_CS], int start) {
    
    int path[N_OF_CS]; //Stores current path of the computation
    path[0] = start; //Puts the first city in the path
    
    //best_paths stores in it's rows the best paths found by each
    //iteration of the loop below. It has one less row than the number of
    //cities because the first city is already chosen at the start
    //of the algorithm
    int best_paths[N_OF_CS][N_OF_CS];
    double best_lengths[N_OF_CS]; //Current best path length of each iteration
    
    int available[N_OF_CS]; //keeps track of cities not yet visited
    
    int i;
    
    
    for(i=0; i<N_OF_CS; i++) {
        best_paths[i][0] = start; //All paths in the iteration start in the same city
        best_lengths[i] = DBL_MAX;
        available[i] = 1; //Mark all cities as available (notice below)
    }
    available[start] = 0; //The first city is already on the path
    
    #pragma omp parallel firstprivate(i, path, available)
    #pragma omp for schedule(dynamic)
    for(i=0; i<N_OF_CS; i++) {
        if(i != start) {
            path[1] = i;
            available[i] = 0;
            //2 is the current path size, since it contains the starting city and
            //the current city in the iteration
            tsp_aux(path, 2, available, distance_m, best_paths[i], &best_lengths[i]);
            available[i] = 1;
        }
    }
    
    //Check which of the returned solutions is the best
    double solution_length = best_lengths[0];
    int solution = 0;
    for(i=0; i<N_OF_CS; i++) { //We could start at i=1, but might as well print 0 too
        if(i!=start) { 
            printf("Solution %d: %f\n", i, best_lengths[i]);
            print_int_arr(best_paths[i]);
            if(best_lengths[i] < solution_length) {
                solution_length = best_lengths[i];
                solution = i;
            }
        }
    }
    
    //Print solution
    printf("Best path:");
    print_int_arr(best_paths[solution]);
    printf("Path length: %f\n", best_lengths[solution]);
}



void tsp_aux(int path[], int path_size, int available[],
             double distance_m[N_OF_CS-1][N_OF_CS],
             int best_path[], double *best_length) {
    //If the path contains all cities, the current branch of the computation
    //tree is finished. If this path is better than the previously recorded path,
    //it is saved along with its length.
    if(path_size == N_OF_CS) {
        //printf("Checkpoint 2\n");
        double path_length = calc_length(path, distance_m);
        if(path_length < *best_length) {
            *best_length = path_length;
            copy_path(path, best_path);
        }
    } else {
        int i;
        //Finds the next city from the list of
        //cities which haven't been visited yet
        for(i=0; i<N_OF_CS; i++) {
            if(available[i]) {
                //Mark city as visited
                available[i] = 0;
                //Includes the new city in the current path
                path[path_size] = i;
                //Go on with the recursion. Notice the increment on path_size.
                tsp_aux(path, path_size+1, available, distance_m, best_path, best_length);
                //At this point the above recursion is backtracking.
                //Mark this city as available again before continuing on this loop,
                //so that deeper levels of the recursion can utilize it.
                available[i] = 1;
            }
        }
    }
}










int main() {
    
    int my_rank;       // Process id
	int proc_n;        // Number of processes (command line -np)
	int message;       // Message buffer
	MPI_Status status; // Status struct

	MPI_Init(&argc , &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_n);
	
	City cities[N_OF_CS];
	cities[0].name = "a"; cities[0].x = 125; cities[0].y = 832;
	cities[1].name = "b"; cities[1].x = 18; cities[1].y = 460;
	cities[2].name = "c"; cities[2].x = 176; cities[2].y = 386;
	cities[3].name = "d"; cities[3].x = 472; cities[3].y = 1000;
	cities[4].name = "e"; cities[4].x = 110; cities[4].y = 57;
	cities[5].name = "f"; cities[5].x = 790; cities[5].y = 166;
	cities[6].name = "g"; cities[6].x = 600; cities[6].y = 532;
	cities[7].name = "h"; cities[7].x = 398; cities[7].y = 40;
	cities[8].name = "i"; cities[8].x = 83; cities[8].y = 720;
	cities[9].name = "j"; cities[9].x = 829; cities[9].y = 627;
	cities[10].name = "k"; cities[10].x = 155; cities[10].y = 567;
	cities[11].name = "l"; cities[11].x = 930; cities[11].y = 106;
	cities[12].name = "m"; cities[12].x = 710; cities[12].y = 266;
	cities[13].name = "n"; cities[13].x = 33; cities[13].y = 680;
	//cities[14].name = "o"; cities[14].x = 672; cities[14].y = 415;
	//print_cities(cities);
	
    
    
    if(my_rank == 0) {
    //================master======================
    
    
	
	
	
	} else {
	//================slave=======================
		while(1) {
			MPI_Recv(&partial_path, GRAIN, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
			
			if(MPI_ANY_TAG == WORK) {
				
				
				
				
				MPI_Send(&solution, )
			} else {
			
			}
		}
    }
    

    
    
    
    //Creates and fills a distance table with distances between all cities
    double distance_m[N_OF_CS][N_OF_CS];
    fill_distance_m(distance_m, cities);
    
    printf("Computing solution using %d threads\n", THREADS);
    double begin = omp_get_wtime();
    tsp(distance_m, 0);
    double end = omp_get_wtime();
    printf("Solution found in %.2f seconds\n", (end - begin));
    //test((double *)distance_m, 3, 5);
}
