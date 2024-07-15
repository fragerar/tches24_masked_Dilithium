COMMON_SOURCES = utils.c random.c
FLAGS = -Wall -O3 -march=native -Wextra

all: $(COMMON_SOURCES)  fastdilithium_gadgets.c dilithium_gadgets.c convtable.c convba_2014.c impconvba.c
	gcc $(FLAGS) -D TEST_FASTDIL $(COMMON_SOURCES) fastdilithium_gadgets.c dilithium_gadgets.c convtable.c convba_2014.c impconvba.c -o main


clean:
	rm -f main 