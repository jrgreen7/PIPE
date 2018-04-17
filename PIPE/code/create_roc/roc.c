#include <stdio.h>
#include <stdlib.h>

#define LINE_LENGTH 64


/* Comparison Function */
int compare(const void* _a, const void* _b)
{
  const float* a = (const float*) _a;
  const float* b = (const float*) _b;

  if(*a > *b) return 1;      // first item is bigger than the second one -> return 1
  else
     if(*a == *b) return 0;  // equality -> return 0
     else return -1;         // second item is bigger than the first one -> return -1
}


void loadFile(char *filename, float *array, int lines)
{
  char line[LINE_LENGTH];
  char *nameA, *nameB, *avg;
  FILE *filePtr;
  int i;

  nameA = (char*) malloc(16*sizeof(char));
  nameB = (char*) malloc(16*sizeof(char));
  avg = (char*) malloc(10*sizeof(char));

  filePtr = fopen(filename, "r");
  if (filePtr == NULL) {
      fprintf(stderr, "ERROR: could not open file %s for reading\n", filename);
      exit(-1);
  }

  for(i=0;i<lines;i++){
    fgets(line, LINE_LENGTH, filePtr);
    sscanf(line, "%s\t%s\t%s\n", nameA, nameB, avg);
    array[i] = atof(avg);    
  } 
  fclose(filePtr);

}


/* Count how many values in an array are over the cutoff */
int countOC(float *array, int num, float cutoff)
{
  int i, count=0;

  for(i=0; i<num; i++){
    if(array[i]>=cutoff) count++;
  }
  return count;
}


/* Calculate Sensitivity and Specificity*/
void calcSensSpec(float *p, float *n, int pnum, int nnum, float *se, float *sp)
{
  int i;

  for(i=0;i<pnum;i++){
    se[i] = (float) (countOC(p, pnum, p[i]))/pnum;
    sp[i] = (float) (nnum - countOC(n, nnum, p[i]))/nnum;
    printf("sens: %.8f, spec: %.8f, cutoff: %.8f\n", se[i], sp[i], p[i]);
  }  
}


void printArray(float *array, int num)
{
  int i;

  for(i=0;i<num;i++){
    printf("%.8f\n", array[i]);
  }
}


void outputROC(float *se, float *sp, int num) 
{
  FILE *out;
  int i;
  float temp=0.0;

  out = fopen("results.roc", "w");

  for(i=0;i<num;i++){    
    if(temp!=(float) 1-sp[i])
      fprintf(out, "%.8f\t%.8f\n", (float) 1-sp[i], se[i]);
    temp = (float) 1-sp[i];
  }
  fclose(out);
}


int main(int argc, char *argv[]) 
{
  float *pos, *neg;  
  float *sens, *spec;
  int numpos, numneg;

  if (argc != 5) {
    printf("Usage: %s <pos file> <#pos> <neg file> <#neg>\n", argv[0]);
    return -1;
  }
  
  numpos = atoi(argv[2]);
  numneg = atoi(argv[4]);
  
  pos = (float*) malloc(numpos*sizeof(float));
  neg = (float*) malloc(numneg*sizeof(float));
  sens = (float*) malloc(numpos*sizeof(float));
  spec = (float*) malloc(numpos*sizeof(float));

  loadFile(argv[1], pos, numpos);   //Open positive CSV file
  loadFile(argv[3], neg, numneg);   //Open negative CSV file

  qsort (pos, numpos, sizeof(float), compare); //Sort positives
  qsort (neg, numneg, sizeof(float), compare); //Sort negatives

  calcSensSpec(pos, neg, numpos, numneg, sens, spec);
  outputROC(sens, spec, numpos);

  return 0;
}
