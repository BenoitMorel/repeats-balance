CC=g++

CFLAGS=-g -W -Wall -ansi -pedantic 
LDFLAGS=-g
EXEC=srlb
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o) parse_sequences.o scan_sequences.o
BISON= bison
FLEX= flex

all: $(EXEC)

srlb: $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

parse%.c: parse%.y
		$(BISON) -d -o $@ $<

%.c: %.l
		$(FLEX) -o $@ $<

test:
	./${EXEC}

clean:
	rm -rf *.o $(EXEC) parse*.h


