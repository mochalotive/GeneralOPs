
# Note: Make crystdist before making restforce

RESTFORCE_PATH = `pwd`
# This should point to the directory containing the crystdist library (crystdist.a):
CRYSTDIST_LIB_PATH = ../crystdist/lib/
# This should point to the directory containing the smcv library (smcv.a):
SMCV_LIB_PATH = ../smcv/lib/
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(RESTFORCE_PATH)/../"
# This should point to the directory containing the "crystdist" folder:
CRYSTDIST_ROOT_PATH = "$(RESTFORCE_PATH)/../"
# This should point to the directory containing the "smcv" folder:
SMCV_ROOT_PATH = "$(RESTFORCE_PATH)/../"

#CC = g++ -m32 
#CC = ../../clang/bin/clang++ -std=c++11 -stdlib=libc++

CC = g++ -std=c++11

CFLAGS = -c -g -Wall -I$(COMMON_ROOT_PATH) -I$(CRYSTDIST_ROOT_PATH) -I$(SMCV_ROOT_PATH)

all: restforce

restforce: src/restforce.o $(CRYSTDIST_LIB_PATH)/crystdist.a $(SMCV_LIB_PATH)/smcv.a
	$(CC) src/restforce.o $(CRYSTDIST_LIB_PATH)/crystdist.a $(SMCV_LIB_PATH)/smcv.a -o restforce
	mkdir -p bin
	mv restforce bin

src/restforce.o: src/restforce.cpp
	$(CC) $(CFLAGS) src/restforce.cpp -o src/restforce.o

clean:
	rm -rf src/*.o bin
