#modied from htslib makefile
FLAGS=-ggdb -std=c++11

CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)

CSRC = $(wildcard *.c)
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

all: distangsd

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


-include $(OBJ:.o=.d)

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp >$*.d

distangsd: $(OBJ)
	$(CXX) $(FLAGS)  -o distangsd *.o $(HTS_LIBDIR) -lz -llzma -lbz2 -lpthread -lcurl -lgsl
else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

distangsd: $(OBJ)
	$(CXX) $(FLAGS)  -o distangsd *.o -lz -llzma -lbz2 -lpthread -lcurl -lhts -lgsl
endif

clean:
	rm  -f distangsd *.o *.d
	
test:
	echo "Only subset of analyses is being tested";
	cd test; ./testAll.sh;
