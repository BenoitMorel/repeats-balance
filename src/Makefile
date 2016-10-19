CC=g++

CFLAGS=-g -W -Wall -ansi -pedantic 
LDFLAGS=-g
EXEC=srlb
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(EXEC)

srlb: $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

test:
	./${EXEC}

clean:
	rm -rf *.o $(EXEC)


