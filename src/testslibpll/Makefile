CC=g++
DBG=  -O3 -g 
CFLAGS=  $(DBG) -W -Wall -ansi -pedantic 
LDFLAGS= $(DBG) -L/home/kaplipa/benoit/github/repeats-balance/src/testslibpll/ 
EXEC=test
SRC= $(wildcard *.cpp)
OTHER_OBJDIR=../generated
OTHER_OBJ=  $(addprefix $(OTHER_OBJDIR)/, parse_sequences.o  scan_sequences.o  parse_partitions.o scan_partitions.o	parse_tree.o scan_tree.o LoadBalancing.o)
OBJDIR=.
OBJ= libpll_helper.o


all: $(EXEC) 


test: $(OBJ) $(OBJDIR)/main.o
	$(CC) $(LDFLAGS) -o $@ $(OTHER_OBJ)  $^ -lpll


$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf $(EXEC) *.o
