CC = clang
#FSANITIZATION = -fsanitize=memory
FSANITIZATION =
#FOPENMP = -fopenmp
FOPENMP =
#FVECTORIZATION = -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize
FVECTORIZATION =
LDFLAGS = $(shell pkg-config --libs gsl) $(FOPENMP) $(FSANITIZATION) \
		  -march=native -lprofiler
CFLAGS = -DNDEBUG -O3 -march=native -ffast-math $(FOPENMP) $(FVECTORIZATION)
#CFLAGS = -g -Wshadow -Werror -Wall -O3 -march=native -ffast-math -std=c99 $(FOPENMP) $(FSANITIZATION)
#CFLAGS = -g -O0 $(FOPENMP) $(FSANITIZATION)

OBJ_DIR = ../obj/
OBJ = model.o helpers.o simplex.o sampler.o
OBJECTS = $(addprefix $(OBJ_DIR), $(OBJ))

all: sampler tags
	cp ./sampler ../bin/

../obj/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

sampler: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o sampler

.PHONY: tags
tags:
	ctags -R .

.PHONY: clean
clean:
	rm -f $(OBJ_DIR)*.o sampler test

