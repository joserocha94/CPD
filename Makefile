
CC = g++
CFLAGS = -g -O3 -fopenmp

TARGET = foxes-rabbits

all: $(TARGET)
 
$(TARGET): $(TARGET).cpp
		$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp
 
clean:
		$(RM) $(TARGET)
