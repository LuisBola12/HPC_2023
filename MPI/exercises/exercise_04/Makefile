PROGRAM=solution_04
CC=mpicc
COMMON_FLAGS=-g -o $(PROGRAM)
SOURCE=solution_04.c


multinode: $(SOURCE)
	$(CC) $(COMMON_FLAGS) $(SOURCE)

clean:
	rm $(PROGRAM)
