CC = clang
LDFLAGS = -lgsl -lgslcblas -lm #-lprofiler
LDTESTFLAGS = -lcheck
#CFLAGS = -DNDEBUG -O2
CFLAGS = -g -O0

all: sampler test tags
	cp ./sampler ../bin/
	cp ./test ../bin/

model.o: model.c model.h
	$(CC) $(CFLAGS) -c $< -o $@

helpers.o: helpers.c helpers.h
	$(CC) $(CFLAGS) -c $< -o $@

vector.o: vector.c vector.h
	$(CC) $(CFLAGS) -c $< -o $@

sampler: sampler.o model.o helpers.o
	$(CC) $(LDFLAGS) $? -o sampler

test: test.o vector.o helpers.o model.o
	$(CC) $(LDFLAGS) $(LDTESTFLAGS) $? -o test

.PHONY: tags
tags:
	ctags -R .

.PHONY: clean
clean:
	rm -f *.o sampler test

