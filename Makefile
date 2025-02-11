CC = g++
CFLAGS =
LDFLAGS = -lpthread

.SUFFIXES : .c .o
.c.o :
	$(CC) -c $(CFLAGS) $<

ALL = fastperm functions

all: $(ALL)

fastperm : fastperm.o functions.o
	$(CC) -o $@ $^ $(LDFLAGS)

functions : functions.o
	$(CC) -o $@ $< $(LDFLAGS)

clean :
	rm -rf *.o $(ALL)