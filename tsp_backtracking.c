#include "tsp_mpi_headers.h"


//Compile using options:
//  -lm to link math library

//N_OF_CS and N_OF_THREADS are defined in the header

/*
  master_routine partially fills the city path then sends it to a slave so that it can
  compute all possible permutations on top of it with the remaining cities. The slave
  will send back a message containing the best permutation found for the partial path
  given. The master then sends the next partial path for the slave to compute.
  
  Parameters:
  message_ptr: pointer to the message struct used for sending. It holds the partial path on
    top of which every slave is going to complete all possible permutations. Also holds
    best path length found until now (this number should start at DBL_MAX).
  available: array marking cities not yet in path.
  best_path: holds best path found during computations.
  burst: pointer. Should start at 1. Used to determine wether the first burst of jobs has
    been sent. Incremented after every send until all slaves receive at least one job.
*/
void master_routine(Message *message_ptr, int available[], int best_path[],
                    int path_size, int *burst) {
    int i;
    if(path_size < N_OF_CS - GRAIN) {
        //Still building the partial path
        for(i=0; i<N_OF_CS; i++) {
            if(available[i]) {
                //City not yet chosen. Add it to the path. Mark it as unavailable.
                message->path[path_size] = i;
                available[i] = 0;
                //Go on with the recursion
                master_routine(message_ptr, availablepath_size+1);
                //Back from recursion. City is available again.
                available[i] = 1;
            } //Else we try the next city in the loop.
        }
        
    } else {
        //Time to forward the rest of the job for one of the slaves
        if((*burst) < proc_n) {
            //First burst of jobs is sent with no need for slave request.
            MPI_Send(message_ptr, sizeof(Message), MPI_BYTE, (*burst), WORK, MPI_COMM_WORLD);
            (*burst)++;
        } else {
            //All slaves received jobs already. Time to colect results
            //before sending anything else.
            Message results;
            MPI_Status status;
            MPI_Recv(&results, sizeof(Message), MPI_BYTE, MPI_ANY_SOURCE, RESULT, &status);
            if(results.best_length < message_ptr->best_length) {
                //A better path has been found.
                //Save it's length.
                message_ptr->best_length = results.best_length;
                //Copy path.
                for(i=0; i<N_OF_CS; i++) {
                    best_path[i] = results.path[i];
                }
            } //Else ignore results received
            
            //Send new job to this slave
            MPI_Send(message_ptr, sizeof(Message), MPI_BYTE, status.MPI_SOURCE, WORK, MPI_COMM_WORLD);
        }
    }
}








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


//path
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
    
    int my_rank;                //Process id
	int proc_n;                 //Number of processes (command line -np)
	Message message;            //Message buffer (see header file)
	int available[N_OF_CS];     //Tells which cities have yet to appear in a permutation
	int best_path[N_OF_CS];
	int i;
	for(i=0; i<N_OF_CS; i++) {
	    message.path[i] = -1;
	    best_path[i] = -1;
	    available[i] = 1;
	}
    message.best_length = DBL_MAX;
	
	MPI_Status status; //Status struct
    
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
	//Creates and fills a distance table with distances between all cities
    double distance_m[N_OF_CS][N_OF_CS];
    fill_distance_m(distance_m, cities);
    
    
    if(my_rank == 0) {
    //================master======================
        //Computation starts defining city 0 as starting city
        message.path[0] = 0;
        //The city is marked as unavailable.
		available[0] = 0;
        int burst = 1;
        master_routine(&message, available, best_path, 0, &burst);
        
        void master_routine(Message *message_ptr, available, int best_path[], int path_size,
                    int *burst, MPI_Status *status_ptr) {

        
    
	
	
	
	} else {
	//================slave=======================
		double best_length;
		
		while(1) {
		    MPI_Recv(&message, sizeof(Message), MPI_BYTE, 1, MPI_ANY_TAG, &status);
			
			if(status.MPI_TAG == WORK) { //Received a job
				//Mark all unavailable cities
				for(i=0; i<N_OF_CS-GRAIN; i++) {
				    available[message.path[i]] = 0;
				}
				best_length = message.best_lenght;
				
				//Work on permutations
				tsp_aux(message.path, N_OF_CS-GRAIN, available, distance_m,
				        best_path, &best_length);
				
				for(i=0; i<N_OF_CS; i++) {
				    //Copy best path to message
				    message.path[i] = best_path[i];
				    //Reset available for next job
				    available[i] = 1;
				}
				//Copy best found length to message
				message.best_length = best_length;
				
				
				//Send results
				MPI_Send(&message, sizeof(Message), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
			} else { //No more work to do
			    break;
			}
		}
    
    
    
    }
    //================wrapping up=====================
    

    
    
    
    
    
    printf("Computing solution using %d threads\n", THREADS);
    double begin = omp_get_wtime();
    tsp(distance_m, 0);
    double end = omp_get_wtime();
    printf("Solution found in %.2f seconds\n", (end - begin));
    //test((double *)distance_m, 3, 5);
}
