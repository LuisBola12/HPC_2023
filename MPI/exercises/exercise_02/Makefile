PROGRAM=solution_02
CC=mpicc
COMMON_FLAGS=-g -o $(PROGRAM)
SOURCE=solution_02.c


multinode: $(SOURCE)
	$(CC) $(COMMON_FLAGS) $(SOURCE)

clean:
	rm $(PROGRAM)
