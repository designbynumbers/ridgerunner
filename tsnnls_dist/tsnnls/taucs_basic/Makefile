
OBJS=	taucs_logging.o \
	taucs_malloc.o \
	taucs_ccs_order.o \
	taucs_ccs_ops.o \
	taucs_vec_base.o \
	taucs_complex.o \
	colamd.o \
	amdbar.o \
	amdexa.o \
	amdtru.o \
	genmmd.o \
	taucs_timer.o \
	taucs_sn_llt.o \
	taucs_ccs_base.o

all: $(OBJS)

clean: 
	rm -f $(OBJS)

taucs_logging.o: taucs.h taucs_logging.c
	$(CC) $(CFLAGS) -c taucs_logging.c

taucs_malloc.o: taucs.h 
	$(CC) $(CFLAGS) -c taucs_malloc.c

taucs_ccs_order.o: taucs.h 
	$(CC) $(CFLAGS) -c taucs_ccs_order.c

taucs_ccs_ops.o: taucs.h 
	$(CC) $(CFLAGS) -c taucs_ccs_ops.c

taucs_vec_base.o: taucs.h 
	$(CC) $(CFLAGS) -c taucs_vec_base.c

taucs_complex.o: taucs.h 
	$(CC) $(CFLAGS) -c taucs_complex.c

colamd.o: taucs.h colamd.h
	$(CC) $(CFLAGS) -c colamd.c

amdbar.o: taucs.h 
	$(CC) $(CFLAGS) -c amdbar.c

amdexa.o: taucs.h 
	$(CC) $(CFLAGS) -c amdexa.c

amdtru.o: taucs.h 
	$(CC) $(CFLAGS) -c amdtru.c

genmmd.o: taucs.h 
	$(CC) $(CFLAGS) -c genmmd.c

taucs_timer.o: taucs.h 
	$(CC) $(CFLAGS) -c taucs_timer.c

taucs_sn_llt.o: taucs.h 
	$(CC) $(CFLAGS) -c taucs_sn_llt.c

taucs_ccs_base.o: taucs.h 
	$(CC) $(CFLAGS) -c taucs_ccs_base.c

