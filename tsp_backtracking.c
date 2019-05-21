#include "tsp_mpi_headers.h"

int proc_n; //Set by command line option -np
int jobs = 0;


//Compile using options:
//  -lm to link math library


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
                message_ptr->path[path_size] = i;
                available[i] = 0;
                //Go on with the recursion
                master_routine(message_ptr, available, best_path, path_size+1, burst);
                //Back from recursion. City is available again.
                available[i] = 1;
            } //Else we try the next city in the loop.
        }
        
    } else {
        //Time to forward the rest of the job for one of the slaves
        if((*burst) < proc_n) {
            //First burst of jobs is sent with no need for slave request.
            MPI_Send(message_ptr, sizeof(Message), MPI_BYTE, (*burst), WORK, MPI_COMM_WORLD);
            printf("job %d sent\n", *burst);
            (*burst)++;
            jobs++;
        } else {
            //All slaves received jobs already. Time to colect results
            //before sending anything else.
            Message results;
            MPI_Status status;
            //printf("Waiting results.\n");
            MPI_Recv(&results, sizeof(Message), MPI_BYTE, MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &status);
            //printf("Results received.\n");
            if(results.best_length < message_ptr->best_length) {
                //A better path has been found.
                //Save it's length.
                message_ptr->best_length = results.best_length;
                printf("%.2f\n", message_ptr->best_length);
                //Copy path.
                for(i=0; i<N_OF_CS; i++) {
                    best_path[i] = results.path[i];
                }
            } //Else ignore results received
            
            //Send new job to this slave
            MPI_Send(message_ptr, sizeof(Message), MPI_BYTE, status.MPI_SOURCE, WORK, MPI_COMM_WORLD);
            jobs++;
        }
    }
    
    
}



//Function called by slaves to compute permutations using all cities not currently
//in the path. Solutions that improve on best_length have their length take it's
//place and their path saved on best_path.
void tsp_aux(int path[], int path_size, int available[],
             double distance_m[N_OF_CS-1][N_OF_CS],
             int best_path[], double *best_length) {
    //If the path contains all cities, the current branch of the computation
    //tree is finished. If this path is better than the previously recorded path,
    //it is saved along with its length.
    if(path_size == N_OF_CS) {
        //printf("Calculating path length.\n");
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
                tsp_aux(path, path_size+1, available, distance_m,
                        best_path, best_length);
                //At this point the above recursion is backtracking.
                //Mark this city as available again before continuing on this loop,
                //so that deeper levels of the recursion can utilize it.
                available[i] = 1;
            }
        }
    }
}










int main(int argc, char **argv) {
    
    int my_rank;             //Process id
    Message message;         //Message buffer (see header file)
    
    //Tells which cities have yet to appear in a permutation.
    int available[N_OF_CS];  
    int best_path[N_OF_CS];
    int i;
    
    //Fill paths with placeholder -1 value.
    //Mark all cities as available.
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
        printf("Computing solution using %d processes\n", proc_n);
        
        double t1,t2;
	t1 = MPI_Wtime();  // inicia a contagem do tempo
	
	//Computation starts defining city 0 as starting city
        message.path[0] = 0;
        //The city is marked as unavailable.
		available[0] = 0;
        int burst = 1;
        master_routine(&message, available, best_path, 1, &burst);
        
        //All work sent. Slaves are blocked on send with their last results.
        //Receive and send final message.
        printf("All work sent. Waiting final results.\n");
        int done = 0;
        Message results;
        while(done<proc_n-1) {
            
            MPI_Recv(&results, sizeof(Message), MPI_BYTE,
                     MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &status);
            if(results.best_length < message.best_length) {
                //A better path has been found.
                //Save it's length.
                message.best_length = results.best_length;
                printf("final check: %.2f\n", message.best_length);
                //Copy path.
                for(i=0; i<N_OF_CS; i++) {
                    best_path[i] = results.path[i];
                }
            } //Else ignore results received.
            
            //Send final message to this slave.
            MPI_Send(&message, sizeof(Message), MPI_BYTE,
                     status.MPI_SOURCE, DIE, MPI_COMM_WORLD);
            done++;
        }
        
        t2 = MPI_Wtime(); // termina a contagem do tempo
        //Print solution
        printf("Best path: ");
        for(i=0; i<N_OF_CS; i++) {
            printf("%d ", best_path[i]);
        }
        printf("\nLength: %.2f\n", message.best_length);
        printf("Time: %.2f seconds\n", t2-t1);
        printf("Jobs sent: %d\n", jobs);
        
	} else {
	//================slave=======================
		double best_length;
		
		while(1) {
		    MPI_Recv(&message, sizeof(Message), MPI_BYTE, 0,
		             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(status.MPI_TAG == WORK) { //Received a job
				//printf("[%d]: job received\n", my_rank);
                //Mark all unavailable cities
				for(i=0; i<N_OF_CS-GRAIN; i++) {
				    available[message.path[i]] = 0;
				}
				best_length = message.best_length;
				
				//Work on permutations
				tsp_aux(message.path, N_OF_CS-GRAIN, available, distance_m,
				        best_path, &best_length);
				//printf("[%d]permutations done\n", my_rank);
                if(best_length<message.best_length) { //Found a better path
                    for(i=0; i<N_OF_CS; i++) {
                        //Copy best path to message
                        message.path[i] = best_path[i];
                        //Reset available for next job
                        available[i] = 1;
                    }
                    //Copy best found length to message
                    message.best_length = best_length;
                    printf("[%d]: %.2f\n", my_rank, message.best_length);
                } //Else no relevant results to report, will send back the
                  //same message it received.
				
				
				//Send back
				MPI_Send(&message, sizeof(Message), MPI_BYTE,
				         0, RESULT, MPI_COMM_WORLD);
                //printf("[%d]: results sent.\n", my_rank);
			} else { //No more work to do
			    printf("[%d]finished\n", my_rank);
			    break;
			}
		}
    
    
    
    }
    //================wrapping up=====================
    MPI_Finalize();
}
