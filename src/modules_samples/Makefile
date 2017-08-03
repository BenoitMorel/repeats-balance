libpath=../../lib/current

CC=g++
DBG= -O3 -g 
CFLAGS=  $(DBG) -W -Wall -ansi -pedantic -I$(libpath)
LDFLAGS= $(DBG) -L$(libpath) -I$(libpath) -L../common
EXTERNAL_LINKS= -lrepeatsbalance -lpll_algorithm -lpll_optimize -lpll_tree -lpll_util -lpll -lm
OBJDIR=obj

EXEC=main
OBJ= $(addprefix $(OBJDIR)/, main.o terraces_experiments.o)


all: dir $(EXEC) 

dir:
	mkdir -p obj 

main: $(OBJ) 
	$(CC) $(CFLAGS) -o $@  $^ $(EXTERNAL_LINKS) $(LDFLAGS)

$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $< $(EXTERNAL_LINKS) $(LDFLAGS)


clean:
	rm -rf $(OBJDIR)/*.o $(EXEC)
