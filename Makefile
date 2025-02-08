CC = g++
CFLAGS =
LDFLAGS = -lpthread

.SUFFIXES : .c .o
.c.o :
	$(CC) -c $(CFLAGS) $<

ALL = main myfuncs fastread

all: $(ALL)

main : main.o myfuncs.o
	$(CC) -o $@ $^ $(LDFLAGS)

myfuncs : myfuncs.o
	$(CC) -o $@ $< $(LDFLAGS)

fastread : fastread.o myfuncs.o
	$(CC) -o $@ $^ $(LDFLAGS)

clean :
	rm -rf *.o $(ALL)