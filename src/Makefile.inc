SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))


PARTNAME :=  $(lastword $(subst /, ,$(CURDIR)))



ifneq (,$(findstring lib,$(PARTNAME)))
all: $(PARTNAME).a
else
all: $(PARTNAME)
endif

$(PARTNAME): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(PARTNAME) $(OBJS) $(LDLIBS)

$(PARTNAME).a: $(OBJS)
	$(AR) rcs $(PARTNAME).a $(OBJS)



clean:
	$(RM) $(PARTNAME) *.a *.o .depend

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

ifneq ($(MAKECMDGOALS),clean)
   -include .depend
endif

