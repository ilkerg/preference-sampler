CC = clang
FSANITIZATION = -fsanitize=memory
FSANITIZATION =
FOPENMP = -fopenmp
#FOPENMP =
LDFLAGS = $(shell pkg-config --libs gsl) $(FOPENMP) $(FSANITIZATION) -lprofiler
LDTESTFLAGS = -lcheck
#CFLAGS = -DNDEBUG -O2 $(FOPENMP)
CFLAGS = -g -O0 $(FOPENMP) $(FSANITIZATION)

all: sampler test tags
	cp ./sampler ../bin/
	cp ./test ../bin/

model.o: model.c model.h
	$(CC) $(CFLAGS) -c $< -o $@

helpers.o: helpers.c helpers.h
	$(CC) $(CFLAGS) -c $< -o $@

vector.o: vector.c vector.h
	$(CC) $(CFLAGS) -c $< -o $@

matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c $< -o $@

simplex.o: simplex.c simplex.h
	$(CC) $(CFLAGS) -c $< -o $@

sampler: sampler.o model.o helpers.o simplex.o
	$(CC) $(LDFLAGS) $? -o sampler

test: vector.o matrix.o simplex.o helpers.o model.o test.o
	$(CC) $(LDFLAGS) $(LDTESTFLAGS) $? -o test

.PHONY: tags
tags:
	ctags -R .

.PHONY: clean
clean:
	rm -f *.o sampler test

