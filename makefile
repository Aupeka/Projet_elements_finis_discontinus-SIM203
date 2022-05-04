#---[ Compiler & Flags]—————————————————

compiler = icc
flags = -std=c99 -O3 -mkl -qopenmp 
links = 

#---[ Makefile ]----------------------------------

headers = $(wildcard include/*.h)
sources = $(wildcard src/*.c)
objects = $(subst src/,obj/,$(sources:.c=.o))
paths += -I./include

executables: main

init: $(objects) $(headers) main_initial.c
	$(compiler) -o dg $(flags) $(objects) $(paths) main_initial.c $(links)

noblas: $(objects) $(headers) mainNoBlas.c
	$(compiler) -o dg $(flags) $(objects) $(paths) mainNoBlas.c $(links)

blas: $(objects) $(headers) mainBlas.c
	$(compiler) -o dg $(flags) $(objects) $(paths) mainBlas.c $(links)

obj/%.o: src/%.c $(wildcard $(subst src/, include/,$(<:.c=.h))) $(wildcard $(subst src/, include/,$(<)))
	$(compiler) -o $@ $(flags) -c $(paths) $<

clean:
	rm -f -r dg obj/*
