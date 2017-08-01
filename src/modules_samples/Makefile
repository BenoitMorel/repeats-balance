libpath=../../lib/current

CC=g++
DBG= -O3 -g 
CFLAGS=  $(DBG) -std=gnu++11 -W -Wall -ansi -pedantic -I$(libpath)
LDFLAGS= $(DBG) -std=gnu++11 -L$(libpath) -I$(libpath) -L../common
EXTERNAL_LINKS=-lrepeatsbalance -lpll_algorithm -lpll_msa -lpll_optimize -lpll_binary -lpll_tree -lpll_util -lpll
OBJDIR=obj

EXEC=main
OBJ= $(addprefix $(OBJDIR)/, main.o terraces_experiments.o)


all: dir $(EXEC) 

dir:
	mkdir -p obj 

main: $(OBJ) 
	$(CC) $(LDFLAGS) -o $@  $^ $(EXTERNAL_LINKS)

$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $< $(EXTERNAL_LINKS)


clean:
	rm -rf $(OBJDIR)/*.o $(EXEC)
