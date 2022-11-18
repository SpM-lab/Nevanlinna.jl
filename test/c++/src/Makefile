TARGET = nevanlinna

SRCS = nevanlinna.cpp

OBJS    = $(SRCS:.cpp=.o)

#CC = icpc
CC = g++

CFLAGS = -I /usr/include/eigen3 -lgmp -lgmpxx

$(TARGET): $(OBJS)
	    $(CC) -o $@ $^ $(CFLAGS) 

.cpp.o:
	$(CC) -c $< $(CFLAGS)

clean:
	-rm -f $(OBJS) $(TARGET) *.d
