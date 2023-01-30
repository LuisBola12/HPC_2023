PROGRAM=solution_01
CC=mpicc
COMMON_FLAGS=-g -o $(PROGRAM)
SOURCE=solution_01.c


multinode: $(SOURCE)
	$(CC) $(COMMON_FLAGS) $(SOURCE)

clean:
	rm $(PROGRAM)
