CC=g++


DBG= -O3 -g -pg
CFLAGS= $(DBG) -W -Wall -ansi -pedantic 
LDFLAGS= $(DBG)
EXEC=srlb tests
SRC= $(wildcard *.cpp)
OBJ= $(OBJDIR)/parse_sequences.o  $(OBJDIR)/scan_sequences.o

OBJDIR=generated
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

$(OBJDIR)/parse%.c: parse%.y
		$(BISON) -d -o $@ $<

$(OBJDIR)/%.c: %.l
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
