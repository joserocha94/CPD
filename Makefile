CC = g++
CFLAGS = -g -O3 -fopenmp

TARGET = foxes-rabbits

SOURCES = foxes-rabbits-omp.cpp

all: $(TARGET)
 
$(TARGET): $(SOURCES)
		$(CC) $(CFLAGS) -o $(TARGET) $(SOURCES)
 
clean:
		rm $(TARGET)
