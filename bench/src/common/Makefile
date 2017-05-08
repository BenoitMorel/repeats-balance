
DBG= -O3 -g 
CFLAGS=$(DBG) -fPIC -W -Wall -ansi -pedantic -I../../lib/current
LDFLAGS=$(DBG) -L../../lib/current -I../../lib/current
OBJDIR=obj
EXEC=librepeatsbalance.so
OBJ= $(addprefix $(OBJDIR)/, Model.o Partition.o Tree.o LikelihoodEngine.o PartitionIntervals.o LoadBalancer.o Timer.o)


all: dir $(EXEC) 

dir:
	mkdir -p obj


$(EXEC): $(OBJ) 
	$(CC) -shared $(LDFLAGS) -o $@  $^ -lpll -lcudart 

$(OBJDIR)/%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	rm -rf $(OBJDIR)/*.o $(EXEC)
