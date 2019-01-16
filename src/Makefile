CC = clang
LDFLAGS = -lgsl -lgslcblas -lm -lcheck #-lprofiler
#CFLAGS = -DNDEBUG -O2
CFLAGS = -g -O0
common-objects = helpers.o model.o vector.o

all: sampler test tags
	cp ./sampler ../bin/
	cp ./test ../bin/

sampler: sampler.o $(common-objects)

test: test.o $(common-objects)

.PHONY: tags
tags:
	ctags -R .

.PHONY: clean
clean:
	rm -f *.o sampler test

