#include "PIPE.h"


void init_mapback() {
	int i;
	for(i=0;i<256;i++)
		MAPBACK[i] = -1;
	MAPBACK['A'] = 0; 	MAPBACK['a'] = 0;
	MAPBACK['R'] = 1; 	MAPBACK['r'] = 1;
	MAPBACK['N'] = 2; 	MAPBACK['n'] = 2;
	MAPBACK['D'] = 3; 	MAPBACK['d'] = 3;
	MAPBACK['C'] = 4; 	MAPBACK['c'] = 4;
	MAPBACK['Q'] = 5; 	MAPBACK['q'] = 5;
	MAPBACK['E'] = 6; 	MAPBACK['e'] = 6;
	MAPBACK['G'] = 7; 	MAPBACK['g'] = 7;
	MAPBACK['H'] = 8; 	MAPBACK['h'] = 8;
	MAPBACK['I'] = 9; 	MAPBACK['i'] = 9;
	MAPBACK['L'] = 10;	MAPBACK['l'] = 10;
	MAPBACK['K'] = 11;	MAPBACK['k'] = 11;
	MAPBACK['M'] = 12;	MAPBACK['m'] = 12;
	MAPBACK['F'] = 13;	MAPBACK['f'] = 13;
	MAPBACK['P'] = 14;	MAPBACK['p'] = 14;
	MAPBACK['S'] = 15;	MAPBACK['s'] = 15;
	MAPBACK['T'] = 16;	MAPBACK['t'] = 16;
	MAPBACK['W'] = 17;	MAPBACK['w'] = 17;
	MAPBACK['Y'] = 18;	MAPBACK['y'] = 18;
	MAPBACK['V'] = 19;	MAPBACK['v'] = 19;
	MAPBACK['B'] = 20;	MAPBACK['B'] = 20;
	MAPBACK['Z'] = 21;	MAPBACK['Z'] = 21;
	MAPBACK['X'] = 22;	MAPBACK['x'] = 22;
	MAPBACK['U'] = 22;	MAPBACK['u'] = 22;	
	return;
}

// Following 2 functions are Brad's for dealing with multiple organisms.
int return_organism(int prot_id, const int ranges[], unsigned int num_org) {
    int i;
    for (i = 1; i <=num_org; i++) {
        if (prot_id <= ranges[i]) {
            return i-1;
        }
    }
    fprintf(stderr, "ERROR: could not find corresponding organism for %d prot\n", prot_id);
    exit(-1);
}


int proper_organism(int prot_id, const int ranges[], unsigned int num_org, const int valid[]) {
    int id = return_organism(prot_id, ranges, num_org);
    //printf("organism number is: %d for %d prot_id, and proper?: %d\n", id, prot_id, valid[id]);
    return valid[id];
}


/* count_lines():
 *   fname - name of file
 *   Return Value: number of lines in the file fname
 */
int count_lines(char *fname) {
	FILE *fp;
	static char line[MAX_LINE_LEN];
	int length;
	fp=fopen(fname, "r");
	length=0;
	while (!feof(fp))
	{ 
		/* Quit the loop if fgets() fails (e.g. blank line at end of file?) */
		if (!fgets(line, MAX_LINE_LEN, fp))
			break;

		/* Make sure that we read an entire line, and that there was a line*/
		if (memchr(line, '\n', MAX_LINE_LEN))
			length++;
		else {
			/* This is an error.  This else happens if a line is longer
			 * than line length.  If we don't know how long lines will be,
			 * we could measure line length as well as the number
			 * of lines.  Otherwise, we can just use MAX_LINE_LEN, set to a
			 * large enough value.
			 */
			fprintf(stderr, "ERROR: Lines in %s longer than buffer (%d byte)\n",
					fname, MAX_LINE_LEN);
			exit(1);
		}
	}
	fclose(fp);
	return length;
}



/* read_line():
 *   fid  - file struct, created by the fopen() function
 *   pair - where the two returned char arrays are stored
 *   Return Value: false (0) if at the End-Of-File
 */
int read_line(FILE *fid, Token_Pair *pair) {
	static char line[MAX_LINE_LEN];
	int col1_len;
	int col2_len;
	int tab_pos;

	/* Read a line */
	if (!fgets(line, MAX_LINE_LEN, fid))
		return 0; /* Return false */

	/* Make sure that we read an entire line */
	if (!memchr(line, '\n', MAX_LINE_LEN)) {
		fprintf(stderr, "ERROR: Line longer than buffer\n");
		exit(1);
	}

	/* Find the index of the TAB character */
	tab_pos = strcspn(line, "\t");
	if (tab_pos == strlen(line)) {
		fprintf(stderr, "ERROR: Input file not in valid (tab-separated) format\n");
		exit(1);
	}

	/* Break the string into two columns */
	col1_len = tab_pos + 1;
	pair->col1 = (char *)malloc(sizeof(char) * col1_len);
	strncpy(pair->col1, line, tab_pos);
	pair->col1[tab_pos] = '\0';  /* Char-arrays are null-terminated */

	col2_len = strlen(line)-tab_pos-1;
	pair->col2 = (char *)malloc(sizeof(char) * col2_len);
	strncpy(pair->col2, &line[tab_pos+1], col2_len);
	pair->col2[col2_len-1] = '\0';  /* Char-arrays are null-terminated */
	return 1;   /* return true */
}



/* Compares two segments of length W and returns if they match with
   a score greater or equal to SCORE */
inline int compare_segment(char *segment1, char *segment2)
{  
	int i=0,score=0;

	for(i=0;i<W;i++)
		score = score + PAM120[(int)segment1[i]][(int)segment2[i]];
	
	if(score>=SCORE) 
		return 1;
	else 
		return 0; 
}


inline int min(int a, int b) {
	if (a < b)
		return a;
	else
		return b;
}


inline void slide(char *segment1, char *segment2, int len, int *table, int offset, int protein) {

	int i, score=0;
	int nextIndex, prevIndex;

	/* Get initial score for first window */
	for (i = 0; i < W; i++)
		score += PAM120[(int)segment1[i]][(int)segment2[i]];

	if (score >= SCORE){
    //if((offset*MAX_NUM_PROTEINS + protein) >=MAX_PROTEIN_LEN * MAX_NUM_PROTEINS) printf("ERROR: %d, %d\n", offset*MAX_NUM_PROTEINS + protein, MAX_PROTEIN_LEN * MAX_NUM_PROTEINS);
		//table[offset*MAX_NUM_PROTEINS + protein] = 1;
		table[(size_t) offset*MAX_NUM_PROTEINS + protein]++;
	}

	/* -  nextIndex points to amino acid entering the window
	 * -  prevIndex points to amino acid leaving the window
	 * => Iterate over remaining windows, updating the score and counting hits.
	 */ 
	nextIndex = W;
	prevIndex = 0;
	for (i = 1; i < len; i++) {
		score += PAM120[(int)segment1[nextIndex]][(int)segment2[nextIndex]];
		score -= PAM120[(int)segment1[prevIndex]][(int)segment2[prevIndex]];

	    if (score >= SCORE){
      //if((offset*MAX_NUM_PROTEINS + protein) >=MAX_PROTEIN_LEN * MAX_NUM_PROTEINS) printf("ERROR: %d, %d\n", offset*MAX_NUM_PROTEINS + protein, MAX_PROTEIN_LEN * MAX_NUM_PROTEINS);
			//table[(offset+i)*MAX_NUM_PROTEINS + protein] = 1;
			table[(size_t) (offset+i)*MAX_NUM_PROTEINS + protein]++;
	    }

		nextIndex++;
		prevIndex++;
	}

	return;
}
