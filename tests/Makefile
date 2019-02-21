CC = clang

all: full_conditional_1

%.o: %.c
	$(CC) -g -O0 -c $< -o $@

full_conditional_1: full_conditional_1.o ../obj/model.o ../obj/helpers.o ../obj/set_counter.o
	$(CC) full_conditional_1.o ../obj/model.o ../obj/helpers.o ../obj/set_counter.o -lm -o $@

clean:
	rm -f *.o full_conditional_1

