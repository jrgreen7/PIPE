CC = gcc
LIBS = -lm 
USER_CFLAGS =  -O $(INCLUDE_DIR)

default: all
all: roc

clean: 
	rm -f roc

roc:	roc.c
	$(CC) $(USER_CFLAGS) roc.c -o roc $(LIBS)
