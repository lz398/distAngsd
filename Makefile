#modied from htslib makefile
FLAGS=-ggdb -std=gnu99

CFLAGS += $(FLAGS)


CSRC = $(wildcard *.c)
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

LDFLAGS=-lz -llzma -lbz2 -lpthread -lcurl -lgsl

all: distAngsd

.PHONY: all clean test

# BAMDIR=""
# BDIR=$(realpath $(BAMDIR))

# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif

ifdef EIGEN
$(info EIGEN defined)
EIGEN_INCDIR=$(realpath $(EIGEN))
endif


-include $(OBJ:.o=.d)
ifdef EIGEN
ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR)	-I$(EIGEN_INCDIR)	$*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR)	-I$(EIGEN_INCDIR)	$*.c	>$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR)	-I$(EIGEN_INCDIR)	$*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR)	-I$(EIGEN_INCDIR)	$*.cpp >$*.d

distAngsd: $(OBJ)
	$(CXX) $(FLAGS)  -o distAngsd *.o $(HTS_LIBDIR) $(LDFLAGS)
else
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(EIGEN_INCDIR) $*.c
	$(CC) -MM $(CFLAGS) -I$(EIGEN_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(EIGEN_INCDIR)	$*.cpp
	$(CXX) -MM $(CXXFLAGS)	-I$(EIGEN_INCDIR)	$*.cpp >$*.d

distAngsd: $(OBJ)
	$(CXX) $(FLAGS)  -o distAngsd *.o $(LDFLAGS)
endif
else
ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR)	$*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR)	$*.c	>$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR)	$*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR)	$*.cpp >$*.d

distAngsd: $(OBJ)
	$(CXX) $(FLAGS)  -o distAngsd *.o $(HTS_LIBDIR) $(LDFLAGS)
else
%.o: %.c
	$(CC) -c  $(CFLAGS) $*.c
	$(CC) -MM $(CFLAGS) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)	$*.cpp >$*.d

distAngsd: $(OBJ)
	$(CXX) $(FLAGS)  -o distAngsd *.o -lz $(LDFLAGS)
endif
endif

clean:
	rm  -f distAngsd *.o *.d

test:
	echo "Only subset of analyses is being tested";
	cd test; ./testAll.sh;
