PROGRAM=solution_03
CC=mpicc
COMMON_FLAGS=-g -o $(PROGRAM)
SOURCE=solution_03.c


multinode: $(SOURCE)
	$(CC) $(COMMON_FLAGS) $(SOURCE)

clean:
	rm $(PROGRAM)
