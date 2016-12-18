CC=g++
DBG=  -O3 -g 
CFLAGS=  $(DBG) -W -Wall -ansi -pedantic -I../../src/ 
LDFLAGS= $(DBG) -I../../src/generated

EXEC=export_lb
OBJDEPS := $(wildcard ../../src/generated/*.o)
OBJDEPS := $(filter-out ../../src/generated/main.o, $(OBJDEPS))
OBJDEPS := $(filter-out ../../src/generated/Tests.o, $(OBJDEPS))

all: $(EXEC) 


export_lb: export_lb.o 
	$(CC) $(LDFLAGS) -o $@ $(OBJDEPS)  $^ 


%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf export_lb export_lb.o 
