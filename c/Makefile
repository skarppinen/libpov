CC:=gcc
GSL_FLAGS:=-lgsl
CBLAS_FLAGS:=-lopenblas
LIBNAME:=libpov.so
OPTLEVEL:=-O2
FLAGS:=-Wall -Wextra -Wpedantic $(OPTLEVEL)
.PHONY=clean 

all: random-sweep.o model.o $(LIBNAME) 

random-sweep.o: random-sweep.c random-sweep.h
	$(CC) $(FLAGS) -fpic -c random-sweep.c 

model.o: model.c
	$(CC) $(FLAGS) -fpic -c model.c

$(LIBNAME): random-sweep.o model.o
	$(CC) $(FLAGS) -shared -o $(LIBNAME) random-sweep.o model.o $(GSL_FLAGS) $(CBLAS_FLAGS)

clean:
	rm *.o
	rm *.so

