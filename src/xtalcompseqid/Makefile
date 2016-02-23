CPPFLAGS+=-isystem../../3rdparty/mpich/include/

LDFLAGS+=-L../libxtalcommon
LDLIBS+=-lxtalcommon
LDFLAGS+=-L../libxtalutil
LDLIBS+=-lxtalutil
LDFLAGS+=-L../libxtaldata
LDLIBS+=-lxtaldata
LDFLAGS+=-L../../3rdparty/gsl/build/.libs/
LDFLAGS+=-L../../3rdparty/gsl/build/cblas/.libs/
LDLIBS+=-lgsl -lgslcblas
LDLIBS+=-l:../../3rdparty/mpich/lib/libmpi.a
LDLIBS+=-l:../../3rdparty/mpich/lib/libmpicxx.a


include ../Makefile.inc

