CC=g++
DBG=  -O3 -g  
CFLAGS=  $(DBG) -W -Wall -ansi -pedantic -I../lib/current
LDFLAGS= $(DBG) -L../lib/current -I../lib/current

EXEC=print_tree


all: $(EXEC) 


print_tree: main.o
	$(CC) $(LDFLAGS) -o $@ $(OTHER_OBJ)  $^ -lpll


%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf print_tree main.o 
