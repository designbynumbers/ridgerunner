
# ridgerunner makefile -- needs made liboctrope-1.0 in its dir
#

CC=gcc
CFLAGS= -O3 -Iliboctrope-1.0/linklibsrc -Iliboctrope-1.0/octropesrc -Itsnnls_dist/tsnnls -I.

OBJS=	display.o errors.o linklib_additions.o settings.o \
	dlen.o ridgerunner_main.o stepper.o

all: tsnnls ridgerunner

tsnnls:
	(cd tsnnls_dist/tsnnls/ ; make libtsnnls)
	cp tsnnls_dist/tsnnls/libtsnnls.a .

clean: 
	rm -f $(OBJS)

ridgerunner: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o ridgerunner -L. -Lliboctrope-1.0 -Ltsnnls_dist/external_libraries -loctrope -ltsnnls -lgslcblas -llapacklinux -lblaslinux -lf77blas -lcblas -latlas -lg2c -lm

display.o: display.c gradient/eqedge.h
	$(CC) $(CFLAGS) -c display.c 

errors.o: errors.h
	$(CC) $(CFLAGS) -c errors.c

linklib_additions.o: linklib_additions.c gradient/eqedge.h
	$(CC) $(CFLAGS) -c linklib_additions.c

settings.o: settings.h settings.c
	$(CC) $(CFLAGS) -c settings.c

dlen.o: gradient/dlen.h gradient/dlen.c
	$(CC) $(CFLAGS) -c gradient/dlen.c

ridgerunner_main.o: gradient/ridgerunner_main.c
	$(CC) $(CFLAGS) -c gradient/ridgerunner_main.c

stepper.o: gradient/stepper.h gradient/stepper.c
	$(CC) $(CFLAGS) -c gradient/stepper.c

