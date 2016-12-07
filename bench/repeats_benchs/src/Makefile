CC=g++
DBG=  -O3 -g  
CFLAGS=  $(DBG) -W -Wall -ansi -pedantic -I../lib/current
LDFLAGS= $(DBG) -L../lib/current -I../lib/current

EXEC=pernode_bench
OBJDIR=.
OBJ= common/common.o pernode_bench/pernode_bench.o


all: $(EXEC) 


pernode_bench: $(OBJ) pernode_bench/pernode_bench.o
	$(CC) $(LDFLAGS) -o $@/$@ $(OTHER_OBJ)  $^ -lpll


$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf pernode_bench/pernode_bench *.o pernode_bench/*.o
