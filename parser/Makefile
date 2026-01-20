# Note: Make mymol and crystdist before making copcalc

PARSERINPUT_PATH := `pwd`
COMMON_ROOT_PATH := $(PARSERINPUT_PATH)/..
CRYSTDIST_LIB_PATH = ../crystdist
MYMOL_LIB_PATH = ../mymol


CC := g++ -m64
CFLAGS := -c -g -Wall -Wextra -I$(COMMON_ROOT_PATH)



all: lib/parser.a

lib/parser.a: src/parseinput.o
	mkdir -p lib
	ar cq lib/parser.a src/parseinput.o	

src/parseinput.o: src/parseinput.cpp
	$(CC) $(CFLAGS) src/parseinput.cpp -o src/parseinput.o 


clean:
	rm -rf src/*.o bin lib
