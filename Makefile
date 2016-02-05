GSL      = 3rdparty/gsl/build/.libs/libgsl.a
GSLCBLAS = 3rdparty/gsl/build/cblas/.libs/libgslcblas.a
GSLSRC   = 3rdparty/gsl/configure

MPI      = 3rdparty/mpich/
MPISRC   = 3rdparty/mpich_src/configure
MPIBUILD = 3rdparty/mpich_src/build

all: $(GSL) $(GSLCBLAS) $(MPI)
	$(MAKE) -C ./src/


GSLVERSION=1.16
MPIVERSION=3.2

$(GSLSRC):
	mkdir -p 3rdparty/ && \
	   cd 3rdparty/ && \
	   wget http://ftpmirror.gnu.org/gsl/gsl-$(GSLVERSION).tar.gz && \
	   tar zxf gsl-$(GSLVERSION).tar.gz && \
	   mv gsl-$(GSLVERSION) gsl	   	

$(GSL): $(GSLCBLAS) $(GSLSRC)


$(GSLCBLAS): $(GSLSRC)
	mkdir -p 3rdparty/gsl/build/ && cd 3rdparty/gsl/build/ && ../configure --host=armv6 && $(MAKE) CFLAGS="$(CXXFLAGS)"




$(MPISRC):
	mkdir -p 3rdparty/ && \
	  cd 3rdparty/ && \
	  wget -O mpich.tar.gz http://www.mpich.org/static/downloads/${MPIVERSION}/mpich-${MPIVERSION}.tar.gz && \
	  tar xzf mpich.tar.gz && $(RM) mpich.tar.gz && mv mpich-${MPIVERSION} mpich_src

$(MPI): $(MPISRC)
	$(eval DST := $(CURDIR)/$(MPI) )
	mkdir -p $(MPIBUILD) &&\
	  cd $(MPIBUILD) && \
	  ../configure --prefix=$(DST) --disable-fortran && $(MAKE) CFLAGS="$(CXXFLAGS)" && $(MAKE) install


raspi: CXXFLAGS+=-g3 -Wall -pedantic -static-libgcc -static-libstdc++ -mfloat-abi=hard -mfpu=vfp -ffast-math -fopenmp -Ofast -flto
raspi: $(GSL) $(GSLCBLAS)
	$(MAKE) -C ./src/ CXXFLAGS="$(CXXFLAGS)" CFLAGS="$(CXXFLAGS)"


.SILENT:

clean:
	$(RM) -R 3rdparty/gsl/build/
	$(MAKE) -C ./src/ clean
	$(RM) -R $(MPI)
	$(RM) -R $(MPIBUILD)

distclean: clean
	$(RM) -r 3rdparty/gsl*
	$(RM) -r 3rdparty/mpi*

.PHONY: gsl

