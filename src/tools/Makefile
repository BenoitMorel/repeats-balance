libpllpath=../../lib/current

CC=g++
DBG= -O3 -g 
CFLAGS=  $(DBG) -std=gnu++11 -W -Wall -ansi -pedantic -I$(libpllpath)
LDFLAGS= $(DBG) -std=gnu++11 -L$(libpllpath) -I$(libpllpath) -L../common -static-libstdc++
OBJDIR=obj

EXEC=main
OBJ= $(addprefix $(OBJDIR)/, main.o analyse_partitions.o generate_partitions.o generate_repeats_dataset.o sites_histogram.o eval_split_loss.o)


all: dir $(EXEC) 

dir:
	mkdir -p obj &


main: $(OBJ) 
	$(CC) $(LDFLAGS) -o $@  $^ -lpll  -lrepeatsbalance

$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf $(OBJDIR)/*.o $(EXEC)
