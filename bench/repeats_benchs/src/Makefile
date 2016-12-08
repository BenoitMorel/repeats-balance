CC=g++
DBG=  -O3 -g  
CFLAGS=  $(DBG) -W -Wall -ansi -pedantic -I../../lib/current
LDFLAGS= $(DBG) -L../../lib/current -I../../lib/current

EXEC=pernode_bench rootpath_bench randompath_bench
OBJ= common/common.o 


all: $(EXEC) 


pernode_bench: $(OBJ) pernode_bench/pernode_bench.o
	$(CC) $(LDFLAGS) -o $@/$@ $(OTHER_OBJ)  $^ -lpll

rootpath_bench: $(OBJ) rootpath_bench/rootpath_bench.o
	$(CC) $(LDFLAGS) -o $@/$@ $(OTHER_OBJ)  $^ -lpll

randompath_bench: $(OBJ) randompath_bench/randompath_bench.o
	$(CC) $(LDFLAGS) -o $@/$@ $(OTHER_OBJ)  $^ -lpll



%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf pernode_bench/pernode_bench *.o pernode_bench/*.o rootpath_bench/rootpath_bench rootpath_bench/*.o randompath_bench/randompath_bench randompath_bench/randompath_bench.o
