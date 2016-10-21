CC=g++

CFLAGS=-g -W -Wall -ansi -pedantic 
LDFLAGS=-g
EXEC=srlb
SRC= $(wildcard *.cpp)
OBJ= $(SRC:%.cpp=$(OBJDIR)/%.o) $(OBJDIR)/parse_sequences.o $(OBJDIR)/scan_sequences.o
OBJDIR=generated
BISON= bison
FLEX= flex

all: clean directories $(EXEC) 

directories:
	@mkdir -p $(OBJDIR)

srlb: $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

$(OBJDIR)/parse%.c: parse%.y
		$(BISON) -d -o $@ $<

$(OBJDIR)/%.c: %.l
		$(FLEX) -o $@ $<

test:
	./${EXEC}

clean:
	rm -rf $(EXEC) generated


