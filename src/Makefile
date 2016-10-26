CC=g++


DBG=  -O3 #-g -pg
CFLAGS= $(DBG) -W -Wall -ansi -pedantic 
LDFLAGS= $(DBG)
EXEC=srlb tests
SRC= $(wildcard *.cpp)
OBJDIR=generated
OBJ=  $(addprefix $(OBJDIR)/, parse_sequences.o  scan_sequences.o  parse_partitions.o scan_partitions.o	parse_tree.o scan_tree.o LoadBalancing.o)

BISON= bison
FLEX= flex

all: clean directories $(EXEC) 

directories:
	@mkdir -p $(OBJDIR)

srlb: $(OBJ) $(OBJDIR)/main.o
	$(CC) $(LDFLAGS) -o $@ $^

tests: $(OBJ) $(OBJDIR)/Tests.o
	$(CC) $(LDFLAGS) -o $@ $^


$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

$(OBJDIR)/parse%.c: parsers/parse%.y
		$(BISON) -d -o $@ $<

$(OBJDIR)/%.c: parsers/%.l
		$(FLEX) -o $@ $<

test:
	./tests

run:
	./srlb

time:
	time ./srlb

clean:
	rm -rf $(EXEC) generated

prof:
	gprof -b srlb gmon.out > analysis.txt
