#CC = clang
CC = gcc

#FSANITIZATION = -fsanitize=memory
#FPROFILER = -lprofiler
FOPENMP = -fopenmp

#FVECTORIZATION = -Rpass=loop-vectorize
#FVECTORIZATION = -fopt-info-vec

#LDFLAGS = $(FOPENMP) $(FSANITIZATION) $(FPROFILER) $(shell /usr/local/bin/gsl-config --libs)
LDFLAGS = $(FOPENMP) $(FSANITIZATION) $(FPROFILER) $(shell /usr/bin/gsl-config --libs)

CFLAGS = -DNDEBUG -ggdb -O3 -march=native -ffast-math $(FOPENMP) $(FVECTORIZATION)
#CFLAGS = -g -Wshadow -Werror -Wall -O3 -march=native -ffast-math -std=c99 $(FOPENMP) $(FSANITIZATION)
#CFLAGS = -g -O0 -march=native $(FOPENMP) $(FSANITIZATION)

OBJ_DIR = ../obj/
OBJ = model.o helpers.o sampler.o set_counter.o
OBJECTS = $(addprefix $(OBJ_DIR), $(OBJ))

all: sampler

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

.PHONY: install
install:
	cp ./sampler ../bin/
