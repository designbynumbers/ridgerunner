
OBJS=   tlsqr.o \
        tsnnls.o \
	lsqr.o
	
TAUCS_OBJS=tsnnls/taucs_basic/taucs_logging.o \
	tsnnls/taucs_basic/taucs_malloc.o \
	tsnnls/taucs_basic/taucs_ccs_order.o \
	tsnnls/taucs_basic/taucs_ccs_ops.o \
	tsnnls/taucs_basic/taucs_vec_base.o \
	tsnnls/taucs_basic/taucs_complex.o \
	tsnnls/taucs_basic/colamd.o \
	tsnnls/taucs_basic/amdbar.o \
	tsnnls/taucs_basic/amdexa.o \
	tsnnls/taucs_basic/amdtru.o \
	tsnnls/taucs_basic/genmmd.o \
	tsnnls/taucs_basic/taucs_timer.o \
	tsnnls/taucs_basic/taucs_sn_llt.o \
	tsnnls/taucs_basic/taucs_ccs_base.o

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
