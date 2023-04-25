CC:=gcc
GSL_FLAGS:=-lgsl
CBLAS_FLAGS:=-lopenblas
CRITERION_FLAGS:=-lcriterion
LIBNAME:=libpov.so
OPTLEVEL:=-O2
FLAGS:=-Wall -Wextra $(OPTLEVEL)
TESTFLAGS := $(FLAGS) -Wno-unused-function -Wno-unused-but-set-variable
.PHONY=clean move_shared

all: random-sweep.o model.o $(LIBNAME) move_shared random-sweep-debug.o random-sweep-prof.o test/test-random-sweep test/test-random-sweep-memalloc test/random-sweep-profile

random-sweep.o: random-sweep.c random-sweep.h
	$(CC) $(FLAGS) -fpic -c random-sweep.c 

model.o: model.c
	$(CC) $(FLAGS) -fpic -c model.c

random-sweep-prof.o: random-sweep.c random-sweep.h
	$(CC) $(FLAGS) -pg -c random-sweep.c -o random-sweep-prof.o

$(LIBNAME): random-sweep.o model.o
	$(CC) $(FLAGS) -shared -o $(LIBNAME) random-sweep.o model.o $(GSL_FLAGS) $(CBLAS_FLAGS)

random-sweep-debug.o: random-sweep.c random-sweep.h
	$(CC) $(FLAGS) -c random-sweep.c -o random-sweep-debug.o -DDEBUG

test/test-random-sweep: test/test-random-sweep.c random-sweep-debug.o 
	$(CC) $(TESTFLAGS) -I. test/test-random-sweep.c random-sweep-debug.o -o test/test-random-sweep $(CRITERION_FLAGS) $(GSL_FLAGS) $(CBLAS_FLAGS) 

test/test-random-sweep-memalloc: test/test-random-sweep-memalloc.c random-sweep-debug.o 
	$(CC) $(TESTFLAGS) -I. test/test-random-sweep-memalloc.c random-sweep-debug.o -o test/test-random-sweep-memalloc $(GSL_FLAGS) $(CBLAS_FLAGS) 

test/random-sweep-profile: test/random-sweep-profile.c
	$(CC) $(TESTFLAGS) -I. -pg test/random-sweep-profile.c random-sweep-prof.o -o test/random-sweep-profile $(GSL_FLAGS) $(CBLAS_FLAGS)

move_shared:
	cp $(LIBNAME) $${PERSONAL_C_LIB} 

clean:
	rm *.o
	rm *.so
	rm test/test-random-sweep test/test-random-sweep-memalloc test/random-sweep-profile

