
OBJS=   tlsqr.o \
        tsnnls.o \
	lsqr.o

all: taucs_basic_objs $(OBJS)

clean: 
	rm -f $(OBJS)
	(cd taucs_basic ; make clean)

taucs_basic_objs:
	(cd taucs_basic ; make all)

tlsqr.o: tsnnls.h
	$(CC) $(CFLAGS) -c tlsqr.c

tsnnls.o: tsnnls.h
	$(CC) $(CFLAGS) -c tsnnls.c

lsqr.o: lsqr.h
	$(CC) $(CFLAGS) -c lsqr.c
