CC = mpic++
CFLAGS = -Wall -O3 --std=c++11

parallel.exe : parallel.o random.o
	$(CC) random.o parallel.o -o parallel.exe
parallel.o : parallel.cpp
	$(CC) -c parallel.cpp -o parallel.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o parallel.exe
cleanoutput :
	rm *.genetic
