
OBJS=   tlsqr.o \
        tsnnls.o \
	lsqr.o

TAUCS_OBJS=taucs_basic/taucs_logging.o \
	taucs_basic/taucs_malloc.o \
	taucs_basic/taucs_ccs_order.o \
	taucs_basic/taucs_ccs_ops.o \
	taucs_basic/taucs_vec_base.o \
	taucs_basic/taucs_complex.o \
	taucs_basic/colamd.o \
	taucs_basic/amdbar.o \
	taucs_basic/amdexa.o \
	taucs_basic/amdtru.o \
	taucs_basic/genmmd.o \
	taucs_basic/taucs_timer.o \
	taucs_basic/taucs_sn_llt.o \
	taucs_basic/taucs_ccs_base.o

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

libtsnnls: taucs_basic_objs $(OBJS)
	ar cru libtsnnls.a $(OBJS) $(TAUCS_OBJS)
	ranlib libtsnnls.a
