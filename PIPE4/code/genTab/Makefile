CC = gcc 
MPICC = mpicc
CFLAGS = -O3 -funroll-loops -finline-functions -Wall  
MPIFLAGS = -DPARALLEL

GEN_OBJS = genTab.o routines.o

all: genTab 

routines.o: routines.c PIPE.h
	$(CC) -c $(CFLAGS) $<

genTab:	$(GEN_OBJS)
	$(MPICC) $(CFLAGS) $(MPIFLAGS) $(GEN_OBJS) -o $@ 
genTab.o: genTab.c PIPE.h
	$(MPICC) -c $(CFLAGS) $(MPIFLAGS) $<

clean: 
	rm -f genTab *.o

