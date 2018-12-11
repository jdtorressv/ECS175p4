all: p4

p4: main.o
	g++ main.o -o p4 -lglut -lGL 
main.o: main.cpp  
	g++ -c main.cpp -lglut -lGL 
clean:
	rm -rf *o p4
