
SMCV_PATH = `pwd`
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(SMCV_PATH)/../"
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(SMCV_PATH)/../"
# This should point to the directory containing the "crystdist" folder:
CRYSTDIST_ROOT_PATH = "$(SMCV_PATH)/../"
SMCV_ROOT_PATH = "$(SMCV_PATH)/../"

CC = g++ -m64
#CC = g++ -m32
CFLAGS = -c -g -Wall -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) -I$(CRYSTDIST_ROOT_PATH) -I$(SMCV_ROOT_PATH)

all: lib/smcv.a

lib/smcv.a: src/opdfile.o
	mkdir -p lib
	ar cq lib/smcv.a src/opdfile.o

debug: lib/smcv.a src/smcv.o 
	mkdir -p bin
	$(CC) src/smcv.o lib/smcv.a -o bin/smcv

src/smcv.o: src/smcv.cpp
	$(CC) $(CFLAGS) src/smcv.cpp -o src/smcv.o

src/opdfile.o: src/opdfile.cpp
	$(CC) $(CFLAGS) src/opdfile.cpp -o src/opdfile.o

clean:
	rm -rf src/*o lib bin

