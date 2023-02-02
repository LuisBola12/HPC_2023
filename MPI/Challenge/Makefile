PROGRAM=challenge
CC=mpicc
COMMON_FLAGS=-g -o $(PROGRAM) -lm
SOURCE=challenge.c


multinode: $(SOURCE)
	$(CC) $(COMMON_FLAGS) $(SOURCE)

clean:
	rm $(PROGRAM)
