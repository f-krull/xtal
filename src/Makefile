
CXXFLAGS+=-g3 -Wall -pedantic -static-libgcc -static-libstdc++
CXXFLAGS+=-march=native # change, if running elsewhere (core2, sandybridge, ...)
CXXFLAGS+=-fopenmp
CXXFLAGS+=-Ofast
CXXFLAGS+=-flto
# needed for -flto
GCC_AR  ?= gcc-ar
#GCC_AR := $(AR)

LIBDIRS := $(wildcard libxtal*)
BUILDLIBDIRS := $(LIBDIRS:%=build-%)
CLEANLIBDIRS := $(LIBDIRS:%=clean-%)

BINDIRS := $(wildcard xtal*)
BUILDBINDIRS := $(BINDIRS:%=build-%)
CLEANBINDIRS := $(BINDIRS:%=clean-%)

all: $(BUILDLIBDIRS) $(BUILDBINDIRS)

$(BUILDBINDIRS): $(BUILDLIBDIRS)
	@echo [Compile] $(@:build-%=%)
	$(MAKE) -C $(@:build-%=%) CXXFLAGS="$(CXXFLAGS)" AR="$(GCC_AR)"

$(BUILDLIBDIRS):
	@echo [Compile] $(@:build-%=%)
	$(MAKE) -C $(@:build-%=%) CXXFLAGS="$(CXXFLAGS)" AR="$(GCC_AR)"

clean: $(CLEANLIBDIRS) $(CLEANBINDIRS)

$(CLEANLIBDIRS):
	$(MAKE) -C $(@:clean-%=%) clean

$(CLEANBINDIRS):
	$(MAKE) -C $(@:clean-%=%) clean

.PHONY: $(DIRS)
.PHONY: $(BUILDLIBDIRS)
.PHONY: $(CLEANLIBDIRS)
.PHONY: $(BUILDBINDIRS)
.PHONY: $(CLEANBINDIRS)

