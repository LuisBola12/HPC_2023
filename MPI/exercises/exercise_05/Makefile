PROGRAM=solution_05
CC=mpicc
COMMON_FLAGS=-g -o $(PROGRAM)
SOURCE=solution_05.c


multinode: $(SOURCE)
	$(CC) $(COMMON_FLAGS) $(SOURCE)

clean:
	rm $(PROGRAM)
