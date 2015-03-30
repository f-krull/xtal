GSL      = 3rdparty/gsl/build/.libs/libgsl.a
GSLCBLAS = 3rdparty/gsl/build/cblas/.libs/libgslcblas.a
GSLSRC   = 3rdparty/gsl/configure

all: $(GSL) $(GSLCBLAS)
	$(MAKE) -C ./src/


GSLVERSION=1.16

$(GSLSRC):
	mkdir -p 3rdparty/ && \
	   cd 3rdparty/ && \
	   wget http://ftpmirror.gnu.org/gsl/gsl-$(GSLVERSION).tar.gz && \
	   tar zxf gsl-$(GSLVERSION).tar.gz && \
	   mv gsl-$(GSLVERSION) gsl	   	

$(GSL): $(GSLCBLAS) $(GSLSRC)


$(GSLCBLAS): $(GSLSRC)
	mkdir -p 3rdparty/gsl/build/ && cd 3rdparty/gsl/build/ && ../configure --host=armv6 && $(MAKE) CFLAGS="$(CXXFLAGS)"


raspi: CXXFLAGS+=-g3 -Wall -pedantic -static-libgcc -static-libstdc++ -mfloat-abi=hard -mfpu=vfp -ffast-math -fopenmp -Ofast -flto
raspi: $(GSL) $(GSLCBLAS)
	$(MAKE) -C ./src/ CXXFLAGS="$(CXXFLAGS)" CFLAGS="$(CXXFLAGS)"


.SILENT:

clean:
	$(RM) -R 3rdparty/gsl/build/
	$(MAKE) -C ./src/ clean

distclean: clean
	$(RM) -r 3rdparty/gsl*

.PHONY: gsl

