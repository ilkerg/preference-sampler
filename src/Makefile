#CC = clang
CC = gcc

#FSANITIZATION = -fsanitize=memory
FSANITIZATION =

#FPROFILER = -lprofiler
FPROFILER =

#FOPENMP = -fopenmp
FOPENMP =

#FVECTORIZATION = -Rpass=loop-vectorize
#FVECTORIZATION = -fopt-info-vec
FVECTORIZATION =

LDFLAGS = $(FOPENMP) $(FSANITIZATION) $(FPROFILER) $(shell gsl-config --libs)

CFLAGS = -DNDEBUG -ggdb -O3 -march=native -ffast-math $(FOPENMP) $(FVECTORIZATION)
#CFLAGS = -g -Wshadow -Werror -Wall -O3 -march=native -ffast-math -std=c99 $(FOPENMP) $(FSANITIZATION)
#CFLAGS = -g -O0 $(FOPENMP) $(FSANITIZATION)

OBJ_DIR = ../obj/
OBJ = model.o helpers.o sampler.o set_counter.o
OBJECTS = $(addprefix $(OBJ_DIR), $(OBJ))

all: sampler tags
	cp ./sampler ../bin/

../obj/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

sampler: $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o sampler

.PHONY: tags
tags:
	ctags -R .

.PHONY: clean
clean:
	rm -f $(OBJ_DIR)*.o sampler

