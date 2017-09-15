CC		= gcc
CFLAGS	= -g -Wall
OBJS	= opt_mul.o
PROGRAM	= mul

all:	$(PROGRAM)

$(PROGRAM):	$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROGRAM)

clean:
	rm -f *.o $(PROGRAM)

