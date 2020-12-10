CC = icc

CFLAGS_OMP = -g -O3 -xHost -qopenmp

EXECUTABLES = weighted_sampling uniform_sampling

all: $(EXECUTABLES)

weighted_sampling: weighted_sampling.o tools.o
	$(CC) $(CFLAGS_OMP) -o $@ $^ -lm

uniform_sampling: uniform_sampling.o tools.o
	$(CC) $(CFLAGS_OMP) -o $@ $^ -lm

.cpp.o:
	$(CC) $(CFLAGS_OMP) -c $<

clean:
	rm -f $(EXECUTABLES) *.o
