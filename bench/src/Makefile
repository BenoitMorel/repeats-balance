CC=g++
DBG= -O0 -g 
CFLAGS=  $(DBG) -std=gnu++11 -W -Wall -ansi -pedantic -I../lib/current
LDFLAGS= $(DBG) -std=gnu++11 -L../lib/current -I../lib/current
OBJDIR=obj

EXEC=main
OBJ= $(addprefix $(OBJDIR)/,common.o derivatives.o sumtables.o full_traversal.o pernode.o main.o calibrate_test.o repeats_rates.o export_repeats.o common/Tree.o common/Partition.o common/Model.o common/LikelihoodEngine.o depreciated/PLLHelper.o)


all: dir $(EXEC) 

dir:
	mkdir -p obj && mkdir -p obj/common obj/depreciated


main: $(OBJ) 
	$(CC) $(LDFLAGS) -o $@  $^ -lpll -lcudart 

$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf $(OBJDIR)/*.o $(EXEC)
