#define __MAIN__

#include "PIPE.h"

//#define PROTEIN_SEQ_FILE "./data/protein_sequences.txt"

char* PROTEIN_SEQ_FILE;
char* DBDIR;
char* ORG_SETTINGS;
char* PAIRS_FILE;

typedef struct {
    unsigned int count;
    unsigned int ids[MAX_NEIGHBOURS];
} NeighbourList;

int main(int argc, char *argv[]) {
	int i, j, mype, npes, numProteins,numProteins1, index, seqSize, k;
	int p1, p2, len1, len2;
	int tmp_int, cnt_int;
	int *scratchPad, *scratchPadInts;
	double startTime=0, inputTime=0, seqCmpTime, outputTime;
	char *ptr1, *ptr2;
	char dbName[128];
	char dbName2[128];
	int *seqLen;
    unsigned int neigh_id;
	char **sequences, **seqName;
	FILE *seqFile, *dbFile, *dbFile2;
    NeighbourList *Neighbours;
	//FILE *logFile;
    //char logFileName[20];
	Token_Pair proSeqPair;

	if (argc != 5) {
		printf("Usage: %s PROTEIN_SEQ_FILE DB_DIR ORG_SETTINGS PROTEIN_PAIRS_INDEX\n", argv[0]);
		return -1;
	}
	
	PROTEIN_SEQ_FILE = argv[1];
	DBDIR = argv[2];
	ORG_SETTINGS = argv[3];
	PAIRS_FILE = argv[4];
	
#ifdef PARALLEL
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
#else
	mype = 0;
	npes = 1;
#endif

	init_mapback();

	if (mype == 0) {
		printf("Starting %s...\n", argv[0]);
		printf("Reading in input sequences...");
		fflush(stdout);
		startTime = timer();
	}
    
    ////////////// BRAD's ORGANISM NORMALIZATION CHANGES ////////////////////////
    FILE *org_file;
    unsigned int num_org;
	org_file = fopen(ORG_SETTINGS, "r");
	if (org_file == NULL) {
		fprintf(stderr, "ERROR: could not open organism settings file %s\n", ORG_SETTINGS);
		exit(-1);
	}
    fscanf(org_file, "%u", &num_org);

    printf("There are %u organisms\n", num_org);

    //printf("there are %d organisms!\n", num_org);
    int ranges[num_org+1];
    ranges[0] = -1;
    for (i = 1; i <= num_org; i++) {
        fscanf(org_file, "%d", &ranges[i]);
        ranges[i] += ranges[i-1];
        printf("Ranges[%d] == %d \n", i, ranges[i]);
    }
    ranges[0] = 0;

    // Read in the organsims we want to do genTab for (and later predict)
    int temp_num;
    fscanf(org_file, "%u", &temp_num);
    int valid[num_org];
    memset(valid, 0, sizeof(valid));

    int temp_pos;
    for (i = 0; i < temp_num; i++) {
        fscanf(org_file, "%d", &temp_pos);
        valid[temp_pos] = 1;
    }

    int org_counts[num_org]; // for counting how many similar proteins are seen in each organism


    //// Adding logFile
    //sprintf(logFileName, "logFile_%d.log", mype);
    //logFile = fopen(logFileName,"w");
    //if (logFile == NULL) {
	//	fprintf(stderr, "ERROR: could not open log file for PID %d\n", mype);
	//	exit(-1);
    //}

	numProteins = count_lines(PROTEIN_SEQ_FILE);
    printf("numProteins == %d\n", numProteins);

	seqLen = (int *) malloc(sizeof(int) * numProteins);
	seqName = (char **) malloc(sizeof(char *) * numProteins);
	sequences = (char **) malloc(sizeof(char *) * numProteins);

	/* Open and read in sequence file.
	   int seqLen[] is the length of each sequence
	   char *sequences[] is the MAPBACK sequence string
	   */
	seqFile = fopen(PROTEIN_SEQ_FILE, "r");
	if (seqFile == NULL) {
		fprintf(stderr, "ERROR: could not open protein sequence file %s\n", PROTEIN_SEQ_FILE);
		exit(-1);
	}
	index = 0;
	while (read_line(seqFile, &proSeqPair)) {
		seqName[index] = (char *) malloc(sizeof(char) * (strlen(proSeqPair.col1)+1));
		strcpy(seqName[index], proSeqPair.col1);

		seqSize = strlen(proSeqPair.col2);
		/* printf("%d: %s length %d\n", index, proSeqPair.col1, seqSize); */
		sequences[index] = (char *) malloc(sizeof(char) * seqSize);
		/* Convert ascii values to table index for easier table lookup */
		for (i = 0; i < seqSize; i++) {
			sequences[index][i] = MAPBACK[(int)proSeqPair.col2[i]];
		}
		seqLen[index] = seqSize;

		if (seqSize > MAX_PROTEIN_LEN) {
			fprintf(stderr, "ERROR: %d length (%d) > MAX_PROTEIN_LEN\n", seqSize, MAX_PROTEIN_LEN);
			exit(-1);
		}

		index++;
	}

    // ADDING neighbours list code!
    int i1, j1;
    FILE* filePtr = fopen(PAIRS_FILE, "r");
	if (filePtr == NULL) {
		fprintf(stderr, "ERROR: could not open protein index file %s\n", PAIRS_FILE);
		exit(-1);
	}
    if (fscanf(filePtr, "%d", &numProteins1) < 0) {
        fprintf(stderr, "ERROR: Cannot read Interaction Graph file.\n");
        exit(-1);
    }

	if (numProteins != numProteins1) {
		fprintf(stderr, "New: Diff number of proteins in pairs file than lines in seqs file: %d in seqs, %d in pairs\n",numProteins,numProteins1);
		exit(-1);
	}

    Neighbours = (NeighbourList *) malloc(sizeof(NeighbourList) * numProteins1);
    memset(Neighbours, 0, sizeof(NeighbourList) * numProteins1);

    if (Neighbours == NULL) {
        fprintf(stderr, "ERROR: Slave cannot allocate memory ");
        fprintf(stderr, "for interaction graph.\n");
        exit(-1);
    }


    for(i1=0; i1<numProteins1;i1++) {
        char name[16];

        if (fscanf(filePtr, "%s", name) < 0) {
            fprintf(stderr, "ERROR: Cannot read Interaction Graph file.\n");
            exit(-1);
        }

        //fprintf(stderr, "Comparison:%s\t%s\n", name, seqName[i1]);
        //if (!strcmp(name, seqName[i1])) {
        //    fprintf(stderr, "Diff names %s %s\n", name, seqName[i1]);
        //    exit(-1);
        //}

        if (fscanf(filePtr, "%d", &Neighbours[i1].count) < 0) {
            fprintf(stderr, "ERROR: Cannot read Interaction Graph file.\n");
            exit(-1);
        }
        for (j1=0; j1<Neighbours[i1].count; j1++) {
            if (fscanf(filePtr, "%d", &Neighbours[i1].ids[j1]) < 0) {
                fprintf(stderr, "ERROR: Cannot read Interaction Graph file.\n");
                exit(-1);
            }
        }
    }
    fclose(filePtr);

	if (mype == 0) {
		inputTime = timer();
		printf("done.\n\n");
		printf("Sequence Comparison...\n");
		fflush(stdout);
	}
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	scratchPad = (int *) malloc((size_t) MAX_PROTEIN_LEN * MAX_NUM_PROTEINS * sizeof(int));
	scratchPadInts = (int *) malloc((size_t) MAX_PROTEIN_LEN * MAX_NUM_PROTEINS * sizeof(int));
	//bzero(scratchPad, 2 * MAX_PROTEIN_LEN * MAX_NUM_PROTEINS);
	size_t temp;
	for(temp = 0; temp < (size_t) MAX_PROTEIN_LEN * MAX_NUM_PROTEINS; temp++) {
		scratchPad[temp] = 0;
		scratchPadInts[temp] = 0;
		//if(temp%10000000 == 0)
		//	printf("%ld\n", temp);
	}

	
	for (p1 = mype ; p1 < numProteins; p1 += npes) {
        // If not interested in doing PIPE runs on this organism's sequences, skip
        if (proper_organism(p1, ranges, num_org, valid) != 1) {
            continue;
        }

		printf("PE %d: Processing %d (%s)\n", mype, p1, seqName[p1]);
		//fprintf(logFile,"PE %d: Processing %d (%s)\n", mype, p1, seqName[p1]);
	   	//fflush(logFile);
	   	fflush(NULL);

		ptr1 = sequences[p1];
		len1 = seqLen[p1];

		for (p2 = 0; p2 < numProteins; p2++) {
			ptr2 = sequences[p2];
			len2 = seqLen[p2];

			for (j = 0; j < len2 - W + 1; j++)
				slide(ptr1, ptr2 + j, min(len1 - W + 1, len2 - W - j + 1), scratchPad, 0, p2);

			for (j = 1; j < len1 - W + 1; j++)
				slide(ptr1 + j, ptr2, min(len2 - W + 1, len1 - W - j + 1), scratchPad, j, p2);

		}
		
        // Write modified original file
		sprintf(dbName, "%s/%s.db", DBDIR, seqName[p1]);
		dbFile = fopen(dbName, "w");
		if (dbFile == NULL) {
			fprintf(stderr, "ERROR: could not open file %s for writing\n", dbName);
			exit(-1);
		}
        // Write modified new file
		sprintf(dbName2, "%s/%s.db_star", DBDIR, seqName[p1]);
		dbFile2 = fopen(dbName2, "w");
		if (dbFile2 == NULL) {
			fprintf(stderr, "ERROR: could not open file %s for writing\n", dbName2);
			exit(-1);
		}

		tmp_int = (int) seqLen[p1];
		fwrite(&tmp_int, sizeof(int), 1, dbFile);
		fwrite(&tmp_int, sizeof(int), 1, dbFile2);

		for (i=0; i<seqLen[p1]; i++) {
			cnt_int = 0;
            // Zero out previous organism count trackers for each new window
            memset(org_counts, 0, sizeof(org_counts));
			for (j = 0; j < numProteins; j++) {
				if (scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] != 0) {
				  cnt_int++;
                  // Add count for which organism it is
                  org_counts[return_organism(j,ranges,num_org)] += 1;
                } 
			}
            // Write total number of prots with similar window
			fwrite(&cnt_int, sizeof(int), 1, dbFile);
			fwrite(&cnt_int, sizeof(int), 1, dbFile2);
            // Write similar prot counts for individual organisms
            for (j=0; j<num_org;j++) {
                fwrite(&org_counts[j], sizeof(int), 1, dbFile);
            }
            for (j=0; j<num_org;j++) {
                fwrite(&org_counts[j], sizeof(int), 1, dbFile2);
            }
        }
		
        // Instead of writing length, write number of proteins tested against
		//tmp_int = (int) seqLen[p1];
        tmp_int = numProteins;
		fwrite(&tmp_int, sizeof(int), 1, dbFile);

        // Swap order of iteration from for s in seq, for j in prot to other way
        for (j = 0; j < numProteins; j++) {

			cnt_int = 0;
            for (i=0; i<seqLen[p1]; i++) {
				if (scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] != 0) {
				  cnt_int++;
                } 
			}
			
            // Write total number of prots with similar window
            // but actually write number of windows in protein with similar window to other protein
			fwrite(&cnt_int, sizeof(int), 1, dbFile);
            for (i=0; i<seqLen[p1]; i++) {
				if (scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] != 0)
				{
					tmp_int = (int) i;
					fwrite(&tmp_int, sizeof(int), 1, dbFile);
				}
			}			
            for (i=0; i<seqLen[p1]; i++) {
				if (scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] != 0) {
                    // do plus 1 to all neighbours for same offset
                    for (k = 0; k < Neighbours[j].count; k++) {
                        neigh_id = Neighbours[j].ids[k];
                        scratchPadInts[(size_t) i*MAX_NUM_PROTEINS + neigh_id] += 1;
                    }
                    scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] = 0;
                }			
            }
		}
		fclose(dbFile);

        // Write modified DB file
		
        // Instead of writing length, write number of proteins tested against
        tmp_int = numProteins;
		fwrite(&tmp_int, sizeof(int), 1, dbFile2);

        // Swap order of iteration from for s in seq, for j in prot to other way
        for (j = 0; j < numProteins; j++) {
			cnt_int = 0;
            // write number of windows in protein with similar window to other protein
            for (i=0; i<seqLen[p1]; i++) {
				if (scratchPadInts[(size_t) i*MAX_NUM_PROTEINS + j] != 0) {
				  cnt_int++;
                } 
			}
			fwrite(&cnt_int, sizeof(int), 1, dbFile2);

            for (i=0; i<seqLen[p1]; i++) {
				if (scratchPadInts[(size_t) i*MAX_NUM_PROTEINS + j] != 0)
				{
                    // Write window location
					tmp_int = (int) i;
					fwrite(&tmp_int, sizeof(int), 1, dbFile2);
                    // Write number of times interaction occured
                    tmp_int = scratchPadInts[(size_t) i*MAX_NUM_PROTEINS + j];
					fwrite(&tmp_int, sizeof(int), 1, dbFile2);
                    // zero out scratchpad for next time
                    scratchPadInts[(size_t) i*MAX_NUM_PROTEINS + j] = 0;
				}
			}			
		}
		fclose(dbFile2);
	}

fclose(org_file);
//fclose(logFile);

#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(mype == 0) {
		seqCmpTime = timer();
		printf("done.\n\n");

		printf("Number of proteins:\t%d\n", numProteins);
		seqSize = 0;
		for (i = 0; i < numProteins; i++)
			seqSize += seqLen[i];
		printf("Average protein length:\t%d\n\n", seqSize / numProteins);
		fflush(stdout);
	}


	if (mype == 0) {
		outputTime = timer();
		printf("********** Timing results **********\n");
		printf("  Input Phase:             %f \n", inputTime - startTime);
		printf("  Seq. Comparison Phase:   %f \n", seqCmpTime - inputTime);
		printf("  Total Runtime:           %f \n", outputTime - startTime);
		printf("************************************\n");
		fflush(stdout);
	}

#ifdef PARALLEL
	MPI_Finalize();
#endif
	return 0;
}
