#
# Created 3/22/07. Updated to sibling build June 2012.
#

AM_CPPFLAGS = $(TSNNLS_INCLUDE) $(OCT_INCLUDE) $(PLC_INCLUDE)
LDADD = $(TSNNLS_LTLIB) $(OCT_LTLIB) $(PLC_LTLIB) 
AM_LDFLAGS =  

AM_CFLAGS = -g -finline-functions -Wall

bin_PROGRAMS = ridgerunner strutplot residual 

ridgerunner_SOURCES = errors.c \
	dlen.c \
	linklib_additions.c \
	settings.c \
	ridgerunner_main.c \
        stepper.c \
	dlen.h \
	errors.h \
	display.c \
	ridgerunner.h \
	portability.h \
	free_edge.c \
	globals.c 

strutplot_SOURCES = strutplot.c \
	free_edge.c \
	globals.c \
	ordie.c 

#install various supporting PERL scripts in $(bindir)

dist_bin_SCRIPTS = shortestknot.pl povanimate gridrunner.pl residual.pl batchrunner rcond

# Now include some testing stuff
if HAVEGSL
TESTS = csteptest mrgradtest strutgradtest rownumcvctest
else
TESTS = csteptest rownumcvctest
endif

if HAVEGSL
check_PROGRAMS = csteptest mrgradtest lookahead strutgradtest rownumcvctest
else 
check_PROGRAMS = csteptest lookahead rownumcvctest
endif

lookahead_SOURCES =  errors.c \
	dlen.c \
	linklib_additions.c \
	settings.c \
        stepper.c \
	dlen.h \
	errors.h \
	display.c \
	ridgerunner.h \
	portability.h \
	free_edge.c \
	globals.c \
	mangle.h \
	mangle.c \
	lookahead.c 

csteptest_SOURCES =  errors.c \
	dlen.c \
	linklib_additions.c \
	settings.c \
        stepper.c \
	dlen.h \
	errors.h \
	display.c \
	ridgerunner.h \
	portability.h \
	free_edge.c \
	globals.c \
	mangle.h \
	mangle.c \
	csteptest.c 

mrgradtest_SOURCES =  errors.c \
	dlen.c \
	linklib_additions.c \
	settings.c \
        stepper.c \
	dlen.h \
	errors.h \
	display.c \
	ridgerunner.h \
	portability.h \
	free_edge.c \
	globals.c \
	mangle.h \
	mangle.c \
	mrgradtest.c 

strutgradtest_SOURCES =  errors.c \
	dlen.c \
	linklib_additions.c \
	settings.c \
	stepper.c \
	dlen.h \
	errors.h \
	display.c \
	ridgerunner.h \
	portability.h \
	free_edge.c \
	globals.c \
	mangle.h \
	mangle.c \
	strutgradtest.c 

rownumcvctest_SOURCES =  errors.c \
	dlen.c \
	linklib_additions.c \
	settings.c \
	stepper.c \
	dlen.h \
	errors.h \
	display.c \
	ridgerunner.h \
	portability.h \
	free_edge.c \
	globals.c \
	mangle.h \
	mangle.c \
	rownumcvctest.c 

residual_SOURCES =  errors.c \
	dlen.c \
	linklib_additions.c \
	settings.c \
        stepper.c \
	dlen.h \
	errors.h \
	display.c \
	ridgerunner.h \
	portability.h \
	free_edge.c \
	globals.c \
	mangle.h \
	mangle.c \
	residual.c 
