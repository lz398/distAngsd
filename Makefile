CC	      = clang++
CC_FLAGS	= -g3 -O3 -Wall -I/usr/local/include/ -L/usr/local/include/ -mmacosx-version-min=10.15
EXEC	= test_e
LD_FLAGS	= -L/usr/local/lib -lhts -lgsl -lgslcblas -lcblas -lm
SE_OBJECTS      = bfgs.o	GLtest.o	GL2Dtest.o	ExistingMethods.o	vcftest.o	GL-Reads.o

${EXEC}:  $(SE_OBJECTS)
	$(CC) $(CC_FLAGS) $(SE_OBJECTS) -o ${EXEC} $(LD_FLAGS)
bfgs.o: bfgs.cpp bfgs.h
	$(CC) $(CC_FLAGS) -c bfgs.cpp
GLtest.o: GLtest.cpp GLtest.h shared.h
	$(CC) $(CC_FLAGS) -c GLtest.cpp
GL2Dtest.o: GL2Dtest.cpp GL2Dtest.h GLtest.h shared.h bfgs.h
	$(CC) $(CC_FLAGS) -c GL2Dtest.cpp
ExistingMethods.o: ExistingMethods.cpp ExistingMethods.h GLtest.h shared.h
	$(CC) $(CC_FLAGS) -c ExistingMethods.cpp
vcftest.o: vcftest.cpp vcftest.h GLtest.h shared.h
	$(CC) $(CC_FLAGS) -c vcftest.cpp
GL-Reads.o: GL-Reads.cpp vcftest.h GLtest.h GL2Dtest.h ExistingMethods.h shared.h
	$(CC) $(CC_FLAGS) -c GL-Reads.cpp
clean:
	rm -f ${EXEC} ${SE_OBJECTS}