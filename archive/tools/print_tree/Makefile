CC=g++
DBG=  -O3 -g 
LIBPLL=../../bench/lib/current 
CFLAGS=  $(DBG) -W -Wall -ansi -pedantic -I$(LIBPLL)
LDFLAGS= $(DBG) -L$(LIBPLL) -I$(LIBPLL)

EXEC=print_tree


all: $(EXEC) 


print_tree: main.o
	$(CC) $(LDFLAGS) -o $@ $(OTHER_OBJ)  $^ -lpll


%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf print_tree print_tree.o 
