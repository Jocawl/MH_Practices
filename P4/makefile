# Makefile
TARGETS= bin/practica4

CXXFLAGS=  -I include/ -I obj/

CC = g++ -g -O2

default : $(TARGETS)

$(TARGETS) : obj/random.o obj/QAP.o obj/entrega4.o
	$(CC) -o $@ $(CXXFLAGS) obj/random.o obj/QAP.o obj/entrega4.o

obj/entrega4.o: src/entrega4.cpp include/QAP.h
	$(CC) $(CXXFLAGS) -c -o obj/entrega4.o src/entrega4.cpp

obj/QAP.o: src/QAP.cpp include/QAP.h
	$(CC) $(CXXFLAGS) -c -o obj/QAP.o src/QAP.cpp

obj/random.o: src/random_ppio.c include/random_ppio.h
	$(CC) $(CXXFLAGS) -c -o obj/random.o src/random_ppio.c

obj/timer.o: src/timer.cpp include/timer.h
	$(CC) $(CXXFLAGS) -c -o obj/timer.o src/timer.cpp

clean:
	rm -f obj/*.o
	clear

mrpropper:
	rm -f $(TARGETS)
	clear
