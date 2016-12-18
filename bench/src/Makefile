CC=g++
DBG=  -O3 -g  
CFLAGS=  $(DBG) -W -Wall -ansi -pedantic -I../lib/current
LDFLAGS= $(DBG) -L../lib/current -I../lib/current
OBJDIR=obj

EXEC=main
OBJ= $(addprefix $(OBJDIR)/,common.o full_traversal.o partitioned_full_traversal.o  pernode.o main.o)


all: dir $(EXEC) 

dir:
	mkdir -p obj


main: $(OBJ) 
	$(CC) $(LDFLAGS) -o $@  $^ -lpll

$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -r $(OBJDIR)/*.o $(EXEC)
