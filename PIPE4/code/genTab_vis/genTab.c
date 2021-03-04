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
	int i, j, mype, npes, numProteins, index, seqSize, k;
	int p1, p2, len1, len2;
	int tmp_int, cnt_int, cum_int;
	int *scratchPad;
	double startTime=0, inputTime=0, seqCmpTime, outputTime;
	char *ptr1, *ptr2;
	char dbName[128];
	int *seqLen;
    unsigned int neigh_id;
	char **sequences, **seqName;
	FILE *seqFile, *dbFile;
	Token_Pair proSeqPair;

	if (argc != 4) {
		printf("Usage: %s PROTEIN_SEQ_FILE DB_DIR ORG_SETTINGS\n", argv[0]);
		return -1;
	}
	
	PROTEIN_SEQ_FILE = argv[1];
	DBDIR = argv[2];
	ORG_SETTINGS = argv[3];
	
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
    fclose(org_file);
    int org_counts[num_org]; // for counting how many similar proteins are seen in each organism

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
	size_t temp;
	for(temp = 0; temp < (size_t) MAX_PROTEIN_LEN * MAX_NUM_PROTEINS; temp++) {
		scratchPad[temp] = 0;
	}

	
	for (p1 = mype ; p1 < numProteins; p1 += npes) {
        // If not interested in doing PIPE runs on this organism's sequences, skip
        if (proper_organism(p1, ranges, num_org, valid) != 1) {
            continue;
        }

		printf("PE %d: Processing %d (%s)\n", mype, p1, seqName[p1]);
	   	fflush(NULL);

		ptr1 = sequences[p1]; //ptr1 is current sequence 
		len1 = seqLen[p1]; //current length of sequence 

		for (p2 = 0; p2 < numProteins; p2++) {
			ptr2 = sequences[p2];
			len2 = seqLen[p2];

			for (j = 0; j < len2 - W + 1; j++)
				slide(ptr1, ptr2 + j, min(len1 - W + 1, len2 - W - j + 1), scratchPad, 0, p2); //assuming comparing 2 proteins window by window

			for (j = 1; j < len1 - W + 1; j++)
				slide(ptr1 + j, ptr2, min(len2 - W + 1, len1 - W - j + 1), scratchPad, j, p2); //comparing same proteins again with an offset of a window I think returns a score

		}
		
        // Write database file
		sprintf(dbName, "%s/%s.db", DBDIR, seqName[p1]);
		dbFile = fopen(dbName, "w");
		if (dbFile == NULL) {
			fprintf(stderr, "ERROR: could not open file %s for writing\n", dbName);
			exit(-1);
		}
		// TODO write new file that I will use
		FILE *dbFile_Tom;
		char src[] = "Tom_";
		strcat(src,seqName[p1]);
		dbFile_Tom = fopen(src , "w");
		if (dbFile_Tom == NULL) {
			fprintf(stderr, "ERROR: could not open file %s for writing\n", src);
			exit(-1);
		}
		// TODO

        // STORAGE 1: Write num of amino acids, and info about each window of amino acids
		tmp_int = (int) seqLen[p1];
		//fwrite args: tmp_int size of array to be written, size of each element to be written in bytes, number of elements of size to be written, pointer to file to wirtten 
		fwrite(&tmp_int, sizeof(int), 1, dbFile); // Write length (amino acids) of protein
		// only happens once per file obviously because it is creating new db file everytime 

        // Loop through every window in protein to store information (e.g frequency) of that particular window
		for (i=0; i<seqLen[p1]; i++) { // Lock the window 
			cnt_int = 0;
            memset(org_counts, 0, sizeof(org_counts)); // Zero out previous organism count trackers for each new window
			for (j = 0; j < numProteins; j++) { // Loop the protein 
				if (scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] != 0) {
				  cnt_int++;
                  org_counts[return_organism(j,ranges,num_org)] += 1; // Add count for which organism it is
                } 
			}
            // Write total number of prots with similar window
			fwrite(&cnt_int, sizeof(int), 1, dbFile);
            // Write similar prot counts for individual organisms
            for (j=0; j<num_org;j++) { 	
                fwrite(&org_counts[j], sizeof(int), 1, dbFile);
            }
        }
		
        // STORAGE 2: Write num of other proteins, and index where further information about each protein is stored
        tmp_int = numProteins;
		fwrite(&tmp_int, sizeof(int), 1, dbFile); // Write number of proteins tested against

        // Store lookup table so can find location of information about each protein using its index
        cum_int = 0; // cumulative integer
        for (j = 0; j < numProteins; j++) { // Loop through protein then window, instead of opposite done above
			cnt_int = 0;
            for (i=0; i<seqLen[p1]; i++) {
				if (scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] != 0) {
				  cnt_int++;
                } 
			}
            fwrite(&cum_int, sizeof(int), 1, dbFile); // Write cumulative sum, which is location of protein information
            cum_int += (1+cnt_int); // +1 because also storing number of windows for each protein as well as each window
        }

        // STORAGE 3: Write information about related to DB protein for each other protein (e.g. windows DB protein has that are similar to said protein) ***
	// TODO output file with genTab info I need
        for (j = 0; j < numProteins; j++) { // Loop through protein then window, instead of opposite done above
			cnt_int = 0;
            for (i=0; i<seqLen[p1]; i++) {
				if (scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] != 0) { //only doing comparision that you did slide with
				  cnt_int++;
                } 
			}
			fwrite(&cnt_int, sizeof(int), 1, dbFile); // Write number of windows in protein with similar window to other protein
			//TODO
			//int parse[] = {15};
			//fwrite(&parse, sizeof(int), 1, dbFile);
			// Write the index that protein is being compared to
			fwrite(&j, sizeof(int), 1, dbFile_Tom);
			// TODO
            for (i=0; i<seqLen[p1]; i++) {
				if (scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] != 0)
				{
					tmp_int = (int) i;
					fwrite(&tmp_int, sizeof(int), 1, dbFile); // Write window location in DB protein is similar to other protein 
					//TODO
					//int parse1[] = {14};
					//fwrite(&parse1, sizeof(int), 1, dbFile);
					// Write the window location of the protein that is simillar
					fwrite(&tmp_int, sizeof(int), 1, dbFile_Tom);
                    scratchPad[(size_t) i*MAX_NUM_PROTEINS + j] = 0;
				}
			}			
			// TODO signify new protein being compared
			int parse1[] = {15};
                        fwrite(&parse1, sizeof(int), 1, dbFile);
			// TODO
		}
		fclose(dbFile);
		// TODO
		fclose(dbFile_Tom);
	}


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
