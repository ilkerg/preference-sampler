CC = clang
#FSANITIZATION = -fsanitize=memory
FSANITIZATION =
#FOPENMP = -fopenmp
FOPENMP =
LDFLAGS = $(shell pkg-config --libs gsl) $(FOPENMP) $(FSANITIZATION) \
		  -march=native -lprofiler
LDTESTFLAGS = -lcheck
CFLAGS = -DNDEBUG -O3 $(FOPENMP) -march=native -ffast-math \
		 -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize
#CFLAGS = -g -Wshadow -Werror -Wall -O3 -march=native -ffast-math -std=c99 $(FOPENMP) $(FSANITIZATION)
#CFLAGS = -g -O0 -std=c99 $(FOPENMP) $(FSANITIZATION)

OBJ_DIR = ../obj/
OBJ = model.o helpers.o vector.o matrix.o simplex.o
OBJECTS = $(addprefix $(OBJ_DIR), $(OBJ))

all: sampler test tags
	cp ./sampler ../bin/
	cp ./test ../bin/

../obj/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

sampler: $(OBJECTS) $(OBJ_DIR)sampler.o
	$(CC) $(LDFLAGS) $(OBJECTS) $(OBJ_DIR)/sampler.o -o sampler

test: $(OBJECTS) $(OBJ_DIR)test.o
	$(CC) $(LDFLAGS) $(LDTESTFLAGS) $(OBJECTS) $(OBJ_DIR)/test.o -o test

.PHONY: tags
tags:
	ctags -R .

.PHONY: clean
clean:
	rm -f $(OBJ_DIR)*.o sampler test

