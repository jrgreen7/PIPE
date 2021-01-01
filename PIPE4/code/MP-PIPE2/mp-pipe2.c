/***********************************MP-PIPE2 README****************************
Running Modes: Both can be enabled/disabled at the top of the code and then
recompiled.
    1) LOOCV: Classic LOOCV as performed throughout the PIPE project
    2) Enable PIPE-Sites: This will also run PIPE-Sites on the PIPE Matrix
        and append the results to the output data. If a site cannot be found,
        -1 will be put as output. If this is disabled, -2 will be placed in
        the output (see output format below).
    3) Export Matrices: This will output the unfiltered matrices to the
        local directory.
    4) Localize DB: This will copy the PIPE DB directory to /tmp and read
        the files from here
    5) PRELOAD_DB: This will load the entire PIPE database into RAM on each
        slave before processing

Output Data Format:
    protein_a    protein_b    PIPE_score    Matrix_max    running_time
    sim_weighted_score    site1_height    site1a_start    site1a_end
    site1b_start    site1b_end    site2_height    site2a_start    site2a_end
    site2b_start    site2b_end    site3_height    site3a_start    site3a_end
    site3b_start    site3b_end

To Compile:
    lab         mpicc -O3 -fopenmp -Wall mp-pipe2.c -m64 -lm -o mp-pipe2
    vf cluster  mpcc -O3 -xopenmp -xtarget=ultraT2 -xcache=8/16/4:4096/64/16 -m64 -lm -lmpi mp-pipe2.c -o mp-pipe2
                mpicc -O3 -xopenmp -xtarget=ultraT2 -xcache=8/16/4:4096/64/16 mp-pipe2.c -m64 -o mp-pipe2
    M9000       mpicc -fast -xarch=sparcima -xcache=64/64/2:6144/256/12 -xchip=sparc64vii mp-pipe2.c -o mp-pipe2

To Run:
    mpirun -np <# processes> mp-pipe2 <protein pairs input file> <output file>
    <PIPE Interaction Graph File> <PIPE DB dir (no trailing '/')
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include <omp.h>
#include <assert.h>
#include <stdint.h>
#include <sys/mman.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>

#define ENABLE_PIPE_SITES
//#define PRELOAD_DB
//#define ENABLE_LOOCV
#define EXPORT_MATRICES
#define EXPORT_SW_MATRIX //Kevin: Added to export SW matrix
//#define LOCALIZE_DB

#define WORKTAG 1
#define DIETAG 2
#define SMALL_PACKET_TAG 3
#define DATA_SIZE sizeof(data)  //K: Double memory for SW mat??
static const int WINDOW_SIZE = 20;
static const int FILTER_SIZE = 3;
static const double THRESHOLD = 0;
static const int PACKET_SIZE = 1000;
int THREADS = 2;

//CAT'S constants
#define MAX_DB_FILE 461532
#define NUM_SET_BITS 6080
#define MAX_PROTEIN_LEN 4910
#define MAX_NEIGHBOURS 362

//******
//PIPE-Sites constants
static const float alpha = 0.5;
static const int minBase = 4;
static const int min_peak_height = 0;
//******
char* INTERACTION_GRAPH_FILE;
char* PAIRS_FILE;
char* PIPE_DB_DIR;
char* ORG_SETTINGS;

typedef struct {
    int a;
    int b;
    int max_score;
    double average_score;
    double compute_time;
    double weighted_score;
    double weighted_score_old;
    int site1_height;
    int site1a_start, site1a_end, site1b_start, site1b_end;
    int site2_height;
    int site2a_start, site2a_end, site2b_start, site2b_end;
    int site3_height;
    int site3a_start, site3a_end, site3b_start, site3b_end;
} data;

typedef struct {
    unsigned int first;
    unsigned int second;
} Interaction;

typedef struct {
    int row;
    int col;
} Mat_Point;

typedef struct {
    int first;
    int second;
} training_pair; 

void master(char*, char*);
void slave();
void do_work(char***, MPI_Datatype);
int run_pipe(Interaction*, data*, char***, unsigned int**,
                unsigned int**, double**, unsigned int*, unsigned int*, int, const training_pair [], unsigned int pairs_list_size,
                unsigned int num_org);
void load_packets(data***, int*, char*, int*, data**);
void merge_results(data**, data**, int*, int);
void process_results(data**, int, char***, char*);
void load_names(char*** names, int* num_proteins);
void build_derived_type(MPI_Datatype*);
FILE *safe_fopen(char*, char*);
void readSimilarity(unsigned int*, char*);
int eFilter(unsigned int**, unsigned int**, int, int, int);
int getInteractionRange(unsigned int**, int, int, int*, int*,
                            int*, int*, int*);
Mat_Point getMaxPoint(unsigned int**, int, int);
void check_matrix(unsigned int**, int, int);
uint64_t find_normalization_factor(const unsigned int * const, const unsigned int * const, const training_pair [], unsigned int);

int main(int argc, char** argv) {

    if(argc != 7) {
        printf("\nUSAGE: mpirun -np <# processes> mp-pipe2 ");
        printf("<protein pairs input file> <output file> ");
        printf("<PIPE Interaction Graph File> <Training Pair Indexes> <PIPE DB dir");
        printf("(no trailing '/')> <ORG_SETTINGS_FILE>\n\n");
        return 0;
    }

    INTERACTION_GRAPH_FILE = argv[3];
    PAIRS_FILE = argv[4];
    PIPE_DB_DIR = argv[5];
    ORG_SETTINGS = argv[6];
    double start, end;
    int rank;

    static struct timeval tv;
    static struct timezone tz;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    if (rank == 0) {

        gettimeofday(&tv, &tz);
        start = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;

        master(argv[1], argv[2]);

        gettimeofday(&tv, &tz);
        end = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;

        printf("\nStart time: %f\nEnd time: %f\nTotal time: %f\n\n",
                    start, end, end - start);

    } else {

        #ifdef LOCALIZE_DB
            char command[50];
            sprintf(command, "rm -rf /tmp/PIPE_DB");
            system(command);
            sprintf(command, "cp -r %s/ /tmp/PIPE_DB/", PIPE_DB_DIR);
            system(command);
            PIPE_DB_DIR = "/tmp/PIPE_DB";
        #endif

        slave();

        #ifdef LOCALIZE_DB
            system("rm -r /tmp/PIPE_DB");
        #endif
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}

void master(char* input_file, char* output_file) {

    int i, processes, remaining_packets;
    int small_packet_size = 0, results_size = 0,num_proteins;
    char** names;
    data* small_packet;
    data** packets;
    MPI_Status status;
    MPI_Datatype data_type;

    build_derived_type(&data_type);

    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    load_packets(&packets, &remaining_packets, input_file,
                    &small_packet_size, &small_packet);

    if(small_packet) {
        MPI_Recv(&results_size, 1, MPI_INT, MPI_ANY_SOURCE,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        MPI_Send(small_packet, small_packet_size, data_type,
                    status.MPI_SOURCE, SMALL_PACKET_TAG, MPI_COMM_WORLD);

        MPI_Send(&small_packet_size, 1, MPI_INT, status.MPI_SOURCE,
                    SMALL_PACKET_TAG, MPI_COMM_WORLD);

        printf("Small packet sent to slave %d. %d packets remaining.\n",
                    status.MPI_SOURCE, remaining_packets);
    }

    load_names(&names, &num_proteins);

    while(remaining_packets > 0) {

        MPI_Recv(&results_size, 1, MPI_INT, MPI_ANY_SOURCE,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (results_size > 0) {

            data* results = malloc(results_size * DATA_SIZE);

            if (results == NULL) {
                fprintf(stderr, "ERROR: Master cannot allocate ");
                fprintf(stderr, "memory for incoming results.\n");
                exit(-1);
            }

            MPI_Recv(results, results_size, data_type, status.MPI_SOURCE,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            process_results(&results, results_size, &names, output_file);
            free(results);
            results = NULL;
        }

        remaining_packets--;

        MPI_Send(packets[remaining_packets], PACKET_SIZE, data_type,
                    status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);

        free(packets[remaining_packets]);
        packets[remaining_packets] = NULL;

        printf("Packet sent to slave %d. %d packets remaining.\n",
                    status.MPI_SOURCE, remaining_packets);
    }

    printf("\nCollecting final results...\n");

    for (i = 0; i < (processes -1)*2; i++) {
        MPI_Recv(&results_size, 1, MPI_INT, MPI_ANY_SOURCE,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (results_size > 0) {
            data* results = malloc(results_size * DATA_SIZE);

            if (results == NULL) {
                fprintf(stderr, "ERROR: Master cannot allocate ");
                fprintf(stderr, "memory for incoming results.\n");
                exit(-1);
            }

            MPI_Recv(results, results_size, data_type, status.MPI_SOURCE,
                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            process_results(&results, results_size, &names, output_file);
            free(results);
            results = NULL;
        }
        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, DIETAG, MPI_COMM_WORLD);
    }

    for(i = 0; i < num_proteins; i++) {
        free(names[i]);
        names[i] = NULL;
    }
    free(names);
    names = NULL;

    free(packets);
    packets = NULL;
    printf("Finished.\n");
}

void slave() {

    int i, num_proteins;
    char** names;
    MPI_Datatype data_type;

    build_derived_type(&data_type);
    load_names(&names, &num_proteins);
    do_work(&names, data_type);

    for (i = 0; i < num_proteins; i++) {
        free(names[i]);
        names[i] = NULL;
    }
    free(names);
    names = NULL;
}

void free_matrix(unsigned int **matrix) {
    int i;
    for(i=0; i<MAX_PROTEIN_LEN; i++){
        free(matrix[i]);
        matrix[i] = NULL;
    }
    free(matrix);
    matrix = NULL;
}

// Kevin:  Added to Free the SW Matrix
void free_matrix_sw(double **matrix) {
	    int i;
	    for(i=0; i<MAX_PROTEIN_LEN; i++){
		free(matrix[i]);
		matrix[i] = NULL;
	    }
	    free(matrix);
	    matrix = NULL;
}

void do_work(char*** names, MPI_Datatype data_type) {

    int numProteins, numPairs, check = 1, work_remaining = 0, results_size = 0;
    data* results = (data*) malloc(2 * PACKET_SIZE * DATA_SIZE);
    data* total_work = (data*) malloc(PACKET_SIZE * DATA_SIZE);

    if (total_work == NULL) {
        fprintf(stderr, "ERROR: Slave cannot allocate memory to store work.\n");
        exit(-1);
    }

    data work;
    MPI_Status status;
    Interaction *InteractionsList;
    
    // Brad's organism stuff
    FILE * org_file;
    unsigned int num_org, pairs_list_size;
    int k;

	org_file = safe_fopen(ORG_SETTINGS, "r");
	if (org_file == NULL) {
		fprintf(stderr, "ERROR: could not open organism settings file %s\n", ORG_SETTINGS);
		exit(-1);
	}

    fscanf(org_file, "%u", &num_org);
    fscanf(org_file, "%u", &pairs_list_size);

    training_pair pairs_list[pairs_list_size];

    for (k = 0; k < pairs_list_size; k++) {
        fscanf(org_file, "%u,%u", &pairs_list[k].first, &pairs_list[k].second);
    }

    printf("Number of training pairs: %u\n", pairs_list_size);
    for (k = 0; k < pairs_list_size; k++) {
        printf("Valid training pairs: %u, %u\n", pairs_list[k].first, pairs_list[k].second);
    }
    fclose(org_file);

    // Read in interactions index file
    int i, j;
    FILE* filePtr = safe_fopen(PAIRS_FILE, "r");
    if (fscanf(filePtr, "%d", &numPairs) < 0) {
        fprintf(stderr, "ERROR: Cannot read interactions list file.\n");
        exit(-1);
    }
    InteractionsList = (Interaction *) malloc(sizeof(Interaction) * numPairs);
    memset(InteractionsList, 0, sizeof(Interaction) * numPairs);
    if (InteractionsList == NULL) {
        fprintf(stderr, "ERROR: Slave cannot allocate memory ");
        fprintf(stderr, "for interactions list.\n");
        exit(-1);
    }
    for(i=0; i<numPairs;i++) {
        if (fscanf(filePtr, "%u", &InteractionsList[i].first) < 0) {
            fprintf(stderr, "ERROR: Cannot read interactions list file.\n");
            exit(-1);
        }
        if (fscanf(filePtr, "%u", &InteractionsList[i].second) < 0) {
            fprintf(stderr, "ERROR: Cannot read interactions list file.\n");
            exit(-1);
        }
    }
    fclose(filePtr);

    //**************

    #ifdef PRELOAD_DB
        //unsigned int *db = (unsigned int*)memalign(sysconf(_SC_PAGESIZE),
                                //MAX_DB_FILE*(numProteins*sizeof(unsigned int)));
                                //victoria falls
        unsigned int *db = malloc(MAX_DB_FILE *
                                (numProteins * sizeof(unsigned int)));
        if (db == NULL) {
            fprintf(stderr, "ERROR: Cannot allocate memory "
            fprintf(stderr, "for PIPE database.\n");
            exit(-1);
        }

        for (i=0; i<numProteins; i++)
            readSimilarity(&(db[i*MAX_DB_FILE]), (*names)[i]);

    #endif

    #pragma omp parallel default(shared) private (work,i,j) num_threads(THREADS)
    {
        //CAT'S & Kevin SW Mat Output
        unsigned int **H, **one;
	double **SW; // K: SW matrix
        unsigned int *similarityA, *similarityB;

        H = (unsigned int**) malloc(MAX_PROTEIN_LEN *
                                        sizeof(unsigned int*));
        one = (unsigned int**) malloc(MAX_PROTEIN_LEN *
                                        sizeof(unsigned int*));
	SW = (double**) malloc(MAX_PROTEIN_LEN *
					sizeof(double*));

        if (H == NULL || one == NULL || SW == NULL) {
            fprintf(stderr, "ERROR: Cannot allocate memory ");
            fprintf(stderr, "for PIPE matrices.\n");
            exit(-1);
        }

        for(i=0; i<MAX_PROTEIN_LEN; i++){
            H[i] = (unsigned int*) malloc(MAX_PROTEIN_LEN *
                                                sizeof(unsigned int));
            one[i] = (unsigned int*) malloc(MAX_PROTEIN_LEN *
                                                sizeof(unsigned int));
	    SW[i] = (double*) malloc(MAX_PROTEIN_LEN *
			    			sizeof(double));

            if (H[i] == NULL || one[i] == NULL || SW[i] == NULL) {
                fprintf(stderr, "ERROR: Cannot allocate memory ");
                fprintf(stderr, "for PIPE matrices.\n");
                exit(-1);
            }

            for(j=0; j<(MAX_PROTEIN_LEN); j++) {
                H[i][j] = 0;
                one[i][j] = 0;
		SW[i][j] = 0.0;
            }
        }

        #ifndef PRELOAD_DB
            similarityA = (unsigned int*)malloc(MAX_DB_FILE *
                                                    sizeof(unsigned int));
            similarityB = (unsigned int*)malloc(MAX_DB_FILE *
                                                    sizeof(unsigned int));

            if (similarityA == NULL || similarityB == NULL) {
                fprintf(stderr, "ERROR: Cannot allocate memory ");
                fprintf(stderr, "for PIPE DB files.\n");
                exit(-1);
            }
            memset(similarityA, 0, MAX_DB_FILE*sizeof(unsigned int));
            memset(similarityB, 0, MAX_DB_FILE*sizeof(unsigned int));
        #endif

        //**********
        while (1) {
            #pragma omp critical (edit_results)
            {
                if (check) {
                    if (work_remaining == 0) {

                        MPI_Send(&results_size, 1, MPI_INT, 0, 0,
                                    MPI_COMM_WORLD);

                        if(results_size > 0)
                            MPI_Send(results, results_size, data_type, 0, 0,
                                        MPI_COMM_WORLD);

                        MPI_Recv(total_work, PACKET_SIZE, data_type, 0,
                                 MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                        if (status.MPI_TAG == DIETAG)
                            check = 0;
                        else if (status.MPI_TAG == SMALL_PACKET_TAG)
                            MPI_Recv(&work_remaining, 1, MPI_INT, 0,
                                        MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                        else
                            work_remaining = PACKET_SIZE;

                        results_size = 0;
                    }

                    if (check) {
                        work_remaining--;
                        work = total_work[work_remaining];
                    }
                }
            }

            if (check) {
                #ifdef PRELOAD_DB
                    similarityA = &(db[work.a * MAX_DB_FILE]);
                    similarityB = &(db[work.b * MAX_DB_FILE]);
                #endif

                if (run_pipe(InteractionsList, &work, names, H, one, SW,
                            similarityA, similarityB, numPairs, pairs_list, pairs_list_size, num_org)) {
                    #pragma omp critical (edit_results)
                    {
                        results_size++;
                        results[results_size - 1] = work;
                    }
                }
            } else
                break;
        }
        free_matrix(one);
        free_matrix(H);
	free_matrix_sw(SW);
        #ifndef PRELOAD_DB
            free(similarityA);
            similarityA = NULL;
            free(similarityB);
            similarityB = NULL;
        #endif
    }

    free(InteractionsList);

    MPI_Send(&results_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    if(results_size > 0)
        MPI_Send(results, results_size, data_type, 0, 0, MPI_COMM_WORLD);

    MPI_Recv(total_work, PACKET_SIZE, data_type, 0, MPI_ANY_TAG,
                MPI_COMM_WORLD, &status);

    free(total_work);
    total_work = NULL;
    #ifdef PRELOAD_DB
        free(db);
        db = NULL;
    #endif
}

int run_pipe(Interaction* InteractionsList, data* work, char*** names,
                unsigned int** H, unsigned int** one, double** SW,
                unsigned int* similarityA, unsigned int* similarityB,
                int numPairs, const training_pair pairs_list [], unsigned int pairs_list_size,
                unsigned int num_org) {

    /*** MAIN PIPE ALGORITHM ***/

    int i, j, row, col, seq1_length = 0, seq2_length = 0;
    int seq1_numProt = 0, seq2_numProt = 0;
    double start, end;
    const unsigned int *simPtr1, *simPtr2; //CAT'S
    unsigned int sim_cnt1; //CAT'S
    static struct timeval tv;
    static struct timezone tz;


    // COMMENT OUT IN PARALLEL
    printf("%s\t%s\t\t%d\t%d\n", (*names)[work->a],
              (*names)[work->b], work->a, work->b);

    gettimeofday(&tv, &tz);
    start = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;

    #ifndef PRELOAD_DB
        readSimilarity(similarityA, (*names)[work->a]);
        readSimilarity(similarityB, (*names)[work->b]);
    #endif

    simPtr1 = similarityA;
    seq1_length = *simPtr1++;
    simPtr2 = similarityB;
    seq2_length = *simPtr2++;

    simPtr1 += seq1_length * (1 + num_org);
    simPtr2 += seq2_length * (1 + num_org);

    seq1_numProt = *simPtr1++;
    seq2_numProt = *simPtr2++;

    if (seq1_numProt != seq2_numProt) {
        fprintf(stderr, "ERROR: Length of two db files are different: %d %d\n", seq1_numProt, seq2_numProt);
        exit(-1);
    }

    int a = 0, b = 0; // counters
    int first = 0, second = 0; 
    int numSimA = 0, numSimB = 0;
    int r = 0, c = 0;
    int offsetA, offsetB;
    if (seq1_length >= WINDOW_SIZE && seq2_length >= WINDOW_SIZE) {
        offsetA = 1 + (seq1_length * (num_org + 1)) + 1;
        offsetB = 1 + (seq2_length * (num_org + 1)) + 1;
        for (i = 0; i < numPairs; i++) {
            first = offsetA + seq1_numProt + similarityA[offsetA + InteractionsList[i].first];
            second = offsetB + seq2_numProt + similarityB[offsetB + InteractionsList[i].second];
            numSimA = similarityA[first];
            numSimB = similarityB[second];
            for (a = 0; a < numSimA; a++) {
                for (b = 0; b < (numSimB); b++) {
                    r = similarityA[first + a + 1];
                    c = similarityB[second + b + 1];
                    // TODO this is an issue with weird readings of some kind - shouldn't be necessary
                    if (r < (seq1_length - WINDOW_SIZE + 1) && c < (seq2_length - WINDOW_SIZE + 1) && r >= 0 && c >= 0) {
                        H[r][c] += 1;
                    }
                }
            }
        }
    }   

    #ifdef EXPORT_MATRICES
        char filename[128];
        sprintf(filename, "%s-%s.mat", (*names)[work->a], (*names)[work->b]);
        FILE* matrix_file = fopen(filename, "w");

        for (row = 0; row < (seq1_length - WINDOW_SIZE + 1); row++) {
            for (col = 0; col < (seq2_length - WINDOW_SIZE + 1); col++) {
                fprintf(matrix_file, "%d ", H[row][col]);
            }
            fprintf(matrix_file, "\n");
        }
        fclose(matrix_file);
    #endif

    int cnt=0, max=0;
    double sum=0.0, avg=0.0, avgPreFilter=0.0, avgPreFilterOld=0.0;
    unsigned int total_hits2 = 0;
    //uint16_t oldNorm;

    eFilter(H, one, seq1_length - WINDOW_SIZE + 1,
                seq2_length - WINDOW_SIZE + 1, FILTER_SIZE);

    //SIM WEIGHTING
    for (row = 0; row < (seq1_length - WINDOW_SIZE + 1); row++) {
        sim_cnt1 = similarityA[1 + (row * (num_org + 1))];
        for (col = 0; col < (seq2_length - WINDOW_SIZE + 1); col++) {

            //ORIGINAL PIPE2 SCORES
            sum += (double) one[row][col];
            if(max < H[row][col])
                max = H[row][col];
            cnt++;
            /*********************/

            int hits = H[row][col];
            total_hits2 += hits;
            unsigned int sim_cnt2 = similarityB[1+ (col * (num_org+1))];

            if (hits) {
                avgPreFilterOld += hits * ((uint64_t)1000) /
                                    ((uint64_t)sim_cnt1 * (uint64_t)sim_cnt2);
                
                uint64_t norm_fact = find_normalization_factor(&similarityA[1+(row*(num_org+1))], &similarityB[1+(col*(num_org+1))], pairs_list, pairs_list_size);

                // Should never happen, but can sometimes if files aren't copied properly. This prevents divide by zero error
                if (norm_fact == 0) {
                    norm_fact = 100000;
                    //fprintf(stderr, "Major Issues with %s, %s at row=%d,col=%d\n", (*names)[work->a],(*names)[work->b], row,col);
                    //exit(-1);
                }
                avgPreFilter += hits * ((uint64_t)1000) /
                                    ( norm_fact);

		// Kevin: Adding to the SW Landscape
		SW[row][col] = hits * ((double)1000) / ((double)norm_fact);

            } else {
                ;
            }
        }
    }
    
    // Kevin: Adding SW Matrix Output
    #ifdef EXPORT_SW_MATRIX
    char swfilename[128];
    sprintf(swfilename, "%s-%s_SW.mat", (*names)[work->a], (*names)[work->b]);
    FILE* swmatrix_file = fopen(swfilename, "w");
    for (row = 0; row < (seq1_length - WINDOW_SIZE + 1); row++) {
    	for (col = 0; col < (seq2_length - WINDOW_SIZE + 1); col++) {
		fprintf(swmatrix_file, "%f ", SW[row][col]);
	}
	fprintf(swmatrix_file, "\n");
    }
    fclose(swmatrix_file);
    #endif


    if (cnt != 0) {
        avg = sum / (double) cnt;
        avgPreFilter /= (double) cnt;
        avgPreFilterOld /= (double) cnt;
    } else
        printf("Error! Divide by 0 (proteins %s and %s).\n",
                    (*names)[work->a], (*names)[work->b]);

    #ifdef ENABLE_PIPE_SITES
        work->site1_height = -1;
        work->site1a_start = -1;
        work->site1a_end = -1;
        work->site1b_start = -1;
        work->site1b_end = -1;
        work->site2_height = -1;
        work->site2a_start = -1;
        work->site2a_end = -1;
        work->site2b_start = -1;
        work->site2b_end = -1;
        work->site3_height = -1;
        work->site3a_start = -1;
        work->site3a_end = -1;
        work->site3b_start = -1;
        work->site3b_end = -1;

        if(max >= min_peak_height) {
            if(getInteractionRange(H, seq1_length, seq2_length,
                                    &(work->site1_height),
                                    &(work->site1a_start),
                                    &(work->site1a_end),
                                    &(work->site1b_start),
                                    &(work->site1b_end))) {
                if(getInteractionRange(H, seq1_length, seq2_length,
                                        &(work->site2_height),
                                        &(work->site2a_start),
                                        &(work->site2a_end),
                                        &(work->site2b_start),
                                        &(work->site2b_end))) {
                    getInteractionRange(H, seq1_length, seq2_length,
                                        &(work->site3_height),
                                        &(work->site3a_start),
                                        &(work->site3a_end),
                                        &(work->site3b_start),
                                        &(work->site3b_end));
                }
            }
        }
    #else
        work->site1_height = -2;
        work->site1a_start = -2;
        work->site1a_end = -2;
        work->site1b_start = -2;
        work->site1b_end = -2;
        work->site2_height = -2;
        work->site2a_start = -2;
        work->site2a_end = -2;
        work->site2b_start = -2;
        work->site2b_end = -2;
        work->site3_height = -2;
        work->site3a_start = -2;
        work->site3a_end = -2;
        work->site3b_start = -2;
        work->site3b_end = -2;
    #endif

    work->max_score = max;
    work->average_score = avg;
    work->weighted_score = avgPreFilter;
    work->weighted_score_old = avgPreFilterOld;

    for(i=0; i<(seq1_length - WINDOW_SIZE + 1); i++){
        for(j=0; j<(seq2_length - WINDOW_SIZE + 1); j++) {
            H[i][j]=0;
            one[i][j]=0;
        }
    }

    gettimeofday(&tv, &tz);
    end = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
    work->compute_time = end - start;

    /*printf("%d\t%d\t%.8f\t%d\t%.8f\t%.8f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                work->a, work->b, work->average_score,
                work->max_score, work->compute_time, work->weighted_score,
                work->site1_height, work->site1a_start, work->site1a_end,
                work->site1b_start, work->site1b_end, work->site2_height,
                work->site2a_start, work->site2a_end, work->site2b_start,
                work->site2b_end, work->site3_height, work->site3a_start,
                work->site3a_end, work->site3b_start, work->site3b_end);
    exit(0);*/

    if(avg >= THRESHOLD)
        return 1;
    return 0;
}

void load_packets(data*** packets, int* remaining_packets, char* file_name,
                    int* small_packet_size, data** small_packet) {

    FILE* datafile = fopen(file_name, "r");
    int i, j;
    double total_work;

    if (fscanf(datafile, "%lg\n", &total_work) < 0) {
        fprintf(stderr, "ERROR: Cannot read input file.\n");
        exit(-1);
    }

    *remaining_packets = total_work/PACKET_SIZE;
    *small_packet_size = fmod(total_work, (double)PACKET_SIZE);
    *packets = malloc((total_work/PACKET_SIZE) * DATA_SIZE);

    if (*packets == NULL) {
        fprintf(stderr, "ERROR: Master cannot allocate memory ");
        fprintf(stderr, "to load work from file.\n");
        exit(-1);
    }

    *small_packet = NULL;

    for (i = 0; i < *remaining_packets; i++) {
        (*packets)[i] = malloc(PACKET_SIZE * DATA_SIZE);
        if ((*packets)[i] == NULL) {
            fprintf(stderr, "ERROR: Master cannot allocate memory ");
            fprintf(stderr, "to load work from file.\n");
            exit(-1);
        }

        for (j = 0; j < PACKET_SIZE; j++) {
            int a,b;
            if (fscanf(datafile, "%d\t%d\n", &a, &b) < 0) {
                fprintf(stderr, "ERROR: Cannot read input file.\n");
                exit(-1);
            }
            data new_data = { a, b, 0, 0.0 };
            (*packets)[i][j] = new_data;
        }
    }

    if(*small_packet_size) {
        (*small_packet) = malloc((*small_packet_size) * DATA_SIZE);

        if ((*small_packet) == NULL) {
            fprintf(stderr, "ERROR: Master cannot allocate memory ");
            fprintf(stderr, "to load work from file.\n");
            exit(-1);
        }
        for (j = 0; j < (*small_packet_size); j++) {
            int a,b;
            if (fscanf(datafile, "%d\t%d\n", &a, &b) < 0) {
                fprintf(stderr, "ERROR: Cannot read input file.\n");
                exit(-1);
            }
            data new_data = { a, b, 0, 0.0 };
            (*small_packet)[j] = new_data;
        }
    }

    fclose(datafile);
}

void merge_results(data** results, data** new_results,
                    int* results_size, int new_results_size) {

    int i, old_results_size = *results_size;

    *results_size += new_results_size;
    *results = realloc(*results, ((*results_size) * DATA_SIZE));

    for(i = 0; i < new_results_size; i++)
        (*results)[old_results_size + i] = (*new_results)[i];
}

void process_results(data** results, int results_size,
                        char*** names, char* output_file) {

    FILE* out_file = safe_fopen(output_file, "a");
    int i;

    for (i = 0; i < results_size; i++) {
        fprintf(out_file, "%s\t%s\t%.8f\t%d\t%.8f\t%.8f\t%.8f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                (*names)[(*results)[i].a], (*names)[(*results)[i].b],
                (*results)[i].average_score, (*results)[i].max_score,
                (*results)[i].compute_time, (*results)[i].weighted_score, (*results)[i].weighted_score_old,
                (*results)[i].site1_height, (*results)[i].site1a_start,
                (*results)[i].site1a_end, (*results)[i].site1b_start,
                (*results)[i].site1b_end, (*results)[i].site2_height,
                (*results)[i].site2a_start, (*results)[i].site2a_end,
                (*results)[i].site2b_start, (*results)[i].site2b_end,
                (*results)[i].site3_height, (*results)[i].site3a_start,
                (*results)[i].site3a_end, (*results)[i].site3b_start,
                (*results)[i].site3b_end);
    }
    fclose(out_file);
}

void load_names(char*** names, int* num_proteins) {

    int i, j, count, index;
    char tmp_name[128];
    FILE* datafile = fopen(INTERACTION_GRAPH_FILE, "r");

    if (fscanf(datafile, "%d", num_proteins) < 0) {
        fprintf(stderr, "ERROR: Cannot read Interaction Graph file.\n");
        exit(-1);
    }
    (*names) = (char**)malloc((*num_proteins)*sizeof(char*));

    if ((*names) == NULL) {
        fprintf(stderr, "ERROR: Cannot allocate memory to ");
        fprintf(stderr, "store protein names.\n");
        exit(-1);
    }

    for(i = 0; i < (*num_proteins); i++) {
        (*names)[i] = (char*)malloc(25*sizeof(char));
        if ((*names)[i] == NULL) {
            fprintf(stderr, "ERROR: Cannot allocate memory to ");
            fprintf(stderr, "store protein names.\n");
            exit(-1);
        }
        if (fscanf(datafile, "%s", tmp_name) < 0) {
            fprintf(stderr, "ERROR: Cannot read Interaction Graph file.\n");
            exit(-1);
        }
        strcpy((*names)[i], tmp_name);
        if (fscanf(datafile, "%d", &count) < 0) {
            fprintf(stderr, "ERROR: Cannot read Interaction Graph file.\n");
            exit(-1);
        }

        for (j = 0; j < count; j++) {
            if (fscanf(datafile, "%d", &index) < 0) {
                fprintf(stderr, "ERROR: Cannot read Interaction Graph file.\n");
                exit(-1);
            }
        }
    }
    fclose(datafile);
}

void build_derived_type(MPI_Datatype* message_type_ptr) {

    data* indata = malloc(DATA_SIZE);

    if (indata == NULL) {
        fprintf(stderr, "ERROR: Cannot allocate memory to ");
        fprintf(stderr, "build derived data type.\n");
        exit(-1);
    }

    int blocks[22];
    MPI_Aint displacements[22];
    MPI_Aint addresses[23];
    MPI_Datatype typelist[22];

    typelist[0] = MPI_INT;
    typelist[1] = MPI_INT;
    typelist[2] = MPI_INT;
    typelist[3] = MPI_DOUBLE;
    typelist[4] = MPI_DOUBLE;
    typelist[5] = MPI_DOUBLE;
    typelist[6] = MPI_DOUBLE;
    typelist[7] = MPI_INT;
    typelist[8] = MPI_INT;
    typelist[9] = MPI_INT;
    typelist[10] = MPI_INT;
    typelist[11] = MPI_INT;
    typelist[12] = MPI_INT;
    typelist[13] = MPI_INT;
    typelist[14] = MPI_INT;
    typelist[15] = MPI_INT;
    typelist[16] = MPI_INT;
    typelist[17] = MPI_INT;
    typelist[18] = MPI_INT;
    typelist[19] = MPI_INT;
    typelist[20] = MPI_INT;
    typelist[21] = MPI_INT;

    blocks[0] = blocks[1] = blocks[2] = blocks[3] = blocks[4] = 1;
    blocks[5] = blocks[6] = blocks[7] = blocks[8] = blocks[9] = 1;
    blocks[10] = blocks[11] = blocks[12] = blocks[13] = 1;
    blocks[14] = blocks[15] = blocks[16] = blocks[17] = 1;
    blocks[18] = blocks[19] = blocks[20] = blocks[21] = 1;

    MPI_Address(indata, &addresses[0]);
    MPI_Address(&(indata->a), &addresses[1]);
    MPI_Address(&(indata->b), &addresses[2]);
    MPI_Address(&(indata->max_score), &addresses[3]);
    MPI_Address(&(indata->average_score), &addresses[4]);
    MPI_Address(&(indata->compute_time), &addresses[5]);
    MPI_Address(&(indata->weighted_score), &addresses[6]);
    MPI_Address(&(indata->weighted_score_old), &addresses[7]);
    MPI_Address(&(indata->site1_height), &addresses[8]);
    MPI_Address(&(indata->site1a_start), &addresses[9]);
    MPI_Address(&(indata->site1a_end), &addresses[10]);
    MPI_Address(&(indata->site1b_start), &addresses[11]);
    MPI_Address(&(indata->site1b_end), &addresses[12]);
    MPI_Address(&(indata->site2_height), &addresses[13]);
    MPI_Address(&(indata->site2a_start), &addresses[14]);
    MPI_Address(&(indata->site2a_end), &addresses[15]);
    MPI_Address(&(indata->site2b_start), &addresses[16]);
    MPI_Address(&(indata->site2b_end), &addresses[17]);
    MPI_Address(&(indata->site3_height), &addresses[18]);
    MPI_Address(&(indata->site3a_start), &addresses[19]);
    MPI_Address(&(indata->site3a_end), &addresses[20]);
    MPI_Address(&(indata->site3b_start), &addresses[21]);
    MPI_Address(&(indata->site3b_end), &addresses[22]);

    displacements[0] = addresses[1] - addresses[0];
    displacements[1] = addresses[2] - addresses[0];
    displacements[2] = addresses[3] - addresses[0];
    displacements[3] = addresses[4] - addresses[0];
    displacements[4] = addresses[5] - addresses[0];
    displacements[5] = addresses[6] - addresses[0];
    displacements[6] = addresses[7] - addresses[0];
    displacements[7] = addresses[8] - addresses[0];
    displacements[8] = addresses[9] - addresses[0];
    displacements[9] = addresses[10] - addresses[0];
    displacements[10] = addresses[11] - addresses[0];
    displacements[11] = addresses[12] - addresses[0];
    displacements[12] = addresses[13] - addresses[0];
    displacements[13] = addresses[14] - addresses[0];
    displacements[14] = addresses[15] - addresses[0];
    displacements[15] = addresses[16] - addresses[0];
    displacements[16] = addresses[17] - addresses[0];
    displacements[17] = addresses[18] - addresses[0];
    displacements[18] = addresses[19] - addresses[0];
    displacements[19] = addresses[20] - addresses[0];
    displacements[20] = addresses[21] - addresses[0];
    displacements[21] = addresses[22] - addresses[0];

    MPI_Type_struct(22, blocks, displacements, typelist, message_type_ptr);

    MPI_Type_commit(message_type_ptr);
    free(indata);
    indata = NULL;
}

FILE *safe_fopen(char *filename, char *mode) { //CAT'S
    FILE *filePtr;

    filePtr = fopen(filename, mode);

    if (filePtr == NULL) {
        fprintf(stderr, "ERROR: could not open file %s for '%s'\n",
                            filename, mode);
        exit(-1);
        }

    return filePtr;
}

void readSimilarity(unsigned int *similarity, char *proteinName) { //CAT'S
    char fileName[128];
    FILE *filePtr;
    size_t numRead;

    sprintf(fileName, "%s/%s.db", PIPE_DB_DIR, proteinName);
    filePtr = safe_fopen(fileName, "r");

    numRead = fread(similarity, sizeof(unsigned int), MAX_DB_FILE, filePtr);

    if (numRead <= 0 || !feof(filePtr)) {
        fprintf(stderr, "ERROR: problem reading DB file %s (too big?)\n",
                    fileName);
        exit(1);
    }

    fclose(filePtr);
}

int eFilter(unsigned int** in, unsigned int** out, int rows,
                int cols, int FILTER) {

    int i, j, k = 0, x, y, dir = 1, zeroes = 0, RADIUS = FILTER/2, answer = 0;

    /* Start at position (-1,0) */
    zeroes = (FILTER*FILTER)-(RADIUS*(RADIUS+1));

    for(i=0;i<RADIUS;i++) {
        for(j=0;j<(RADIUS+1);j++)
            zeroes += !in[i][j];
    }

    for(i=0;i<rows;i++) {
        /* Direction = Down */
        /* Remove row at the top */
        y = i-(RADIUS+1);

        if(y >= 0) {
            for(x = (k-RADIUS); x <= (k+RADIUS); x++) {
                if(x >= 0 && x < cols)
                    zeroes -= !in[y][x];
                else
                    zeroes--;
            }
        } else
            zeroes -= FILTER;

        /* Add row at the bottom */
        y = i+RADIUS;
        if(y < rows){
            for(x = (k-RADIUS); x <= (k+RADIUS); x++) {
                if(x >= 0 && x < cols)
                    zeroes += !in[y][x];
                else
                    zeroes++;
            }
        } else
            zeroes+=FILTER;

        if(zeroes <= ((FILTER*FILTER)/2)) {
            out[i][k] = 1;
            answer = 1;
        } else
            out[i][k]=0;

        k += dir;

        while(k >= 0 && k < cols) {
            /* Remove column */
            x = k-(dir*(RADIUS+1));

            if(x >= 0 && x < cols) {
                for(y = (i-RADIUS); y <= (i+RADIUS); y++) {
                    if(y >= 0 && y < rows)
                        zeroes -= !in[y][x];
                    else
                        zeroes--;
                }
            } else
                zeroes-=FILTER;

            /* Add column */
            x = k + (dir*RADIUS);

            if(x >= 0 && x < cols) {
                for(y = (i-RADIUS); y <= (i+RADIUS); y++) {
                    if(y >= 0 && y < rows)
                        zeroes += !in[y][x];
                    else
                        zeroes++;
                }
            } else
                zeroes+=FILTER;

            if(zeroes <= ((FILTER*FILTER)/2)) {
                out[i][k]=1;
                answer=1;
            } else
                out[i][k]=0;
            k+=dir;
        }

        k -= dir;
        dir =~ dir+1;
    }

    return answer;
}

//PIPE-Sites
int getInteractionRange(unsigned int** H, int seqLenA, int seqLenB,
                            int* height, int* a_start, int* a_end,
                            int* b_start, int* b_end){

    Mat_Point curPeak;

    int rows = seqLenA - WINDOW_SIZE + 1;
    int cols = seqLenB - WINDOW_SIZE + 1;
    int curPeakMin, i, j;
    int curPeakLeft = 0;
    int curPeakRight = 0;
    int curPeakUp = 0;
    int curPeakDown = 0; 
    //curPeakLeft = curPeak.col;
    //curPeakRight = curPeak.col;
    //curPeakUp = curPeak.row;
    //curPeakDown = curPeak.row;
    int upDone = 0, downDone = 0, leftDone = 0, rightDone = 0;
    int allDone = 0, stepCount = 0, foundRange = 0;

    while (foundRange == 0){

        curPeak = getMaxPoint(H, rows, cols);

        if ( H[curPeak.row][curPeak.col] == 0){
            *a_start = -1;
            *a_end = -1;
             *b_start = -1;
            *b_end = -1;
            return 0;
        } else {
            curPeakMin = H[curPeak.row][curPeak.col] * alpha;
            curPeakLeft = curPeak.col;
            curPeakRight = curPeak.col;
            curPeakUp = curPeak.row;
            curPeakDown = curPeak.row;
            stepCount=0;
            upDone = 0;
            downDone = 0;
            leftDone = 0;
            rightDone = 0;
            allDone=0;

            while( allDone == 0 ){
                stepCount++;
                //Check Left Boundary
                if ( cols <= 0 || curPeakLeft == 0 ||
                        H[curPeak.row][curPeakLeft - 1 ] < curPeakMin ) {
                    leftDone = 1;
                } else {
                    curPeakLeft--;
                }
                //Check right Boundary
                if ( cols <= 0 || curPeakRight == cols-1 ||
                        H[curPeak.row][curPeakRight +1] < curPeakMin){
                    rightDone = 1;
                } else {
                    curPeakRight++;
                }
                //Check Up Boundary
                if ( rows <= 0 || curPeakUp == 0 ||
                        H[curPeakUp - 1][curPeak.col] < curPeakMin ){
                    upDone = 1;
                } else {
                    curPeakUp--;
                }
                //Check Down Boundary
                if ( rows <= 0 || curPeakDown == rows - 1 ||
                        H[curPeakDown + 1][curPeak.col] < curPeakMin ){
                    downDone = 1;
                } else {
                    curPeakDown++;
                }
                //check to see if we are done in all the directions.
                if( upDone==1 && downDone==1 && leftDone==1 && rightDone==1 ) {
                    allDone = 1;
                }
            }
            *a_start = curPeakUp;
            *a_end = curPeakDown + WINDOW_SIZE;
             *b_start = curPeakLeft;
            *b_end = curPeakRight + WINDOW_SIZE;
            *height = H[curPeak.row][curPeak.col];

            for( i = *a_start; i <= *a_end - WINDOW_SIZE; i++ ){
                for( j = *b_start; j <= *b_end - WINDOW_SIZE; j++ ){
                    H[i][j] = 0;
                }
            }

            if( minBase < stepCount ){
                foundRange=1;
            }
        }
    }

    return 1;
}

//PIPE-Sites
Mat_Point getMaxPoint(unsigned int** H, int rows, int cols) {

    Mat_Point returnPoint;
    returnPoint.row = 0;
    returnPoint.col = 0;
    int i, j, curMax = 0;

    for( i = 0; i < rows; i++ ){
        for( j = 0; j < cols; j++ ){
            if( curMax < H[i][j] ) {
                curMax = H[i][j];
                returnPoint.row = i;
                returnPoint.col = j;
            }
        }
    }
    return returnPoint;
}

void check_matrix(unsigned int** matrix, int rows, int cols) {
    int i, j;

    printf("Checking matrix..\n");

    for( i = 0; i < rows; i++ )
        for( j = 0; j < cols; j++ )
            if(matrix[i][j] != 0)
                printf("matrix[%d][%d] = %d\n", i, j, matrix[i][j]);
    return;
}


uint64_t find_normalization_factor(const unsigned int * const p1pointer, const unsigned int * const p2pointer, const training_pair pairs_list[], unsigned int size) {
    uint64_t norm_factor = 0; // change to one later for better norm theory
    int i;
    for (i=0; i < size; i++) {
        norm_factor += (*(p1pointer + 1 + pairs_list[i].first)) * (*(p2pointer + 1 + pairs_list[i].second));
    }
    if (norm_factor == 0) {
        printf("Norm factor = 0.......BRADSEARCHTERM\n");
    }
    return norm_factor;
}
