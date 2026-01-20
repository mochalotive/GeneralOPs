
# Note: Make crystdist before making smcvstep

SMCVSTEP_PATH = `pwd`
# This should point to the directory containing the mymol library (mymol.a):
MYMOL_LIB_PATH = ../mymol/lib
# This should point to the directory containing the crystdist library (crystdist.a):
CRYSTDIST_LIB_PATH = ../crystdist/lib/
# This should point to the directory containing the smcv library (smcv.a):
SMCV_LIB_PATH = ../smcv/lib
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(SMCVSTEP_PATH)/../"
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(SMCVSTEP_PATH)/../"
# This should point to the directory containing the "crystdist" folder:
CRYSTDIST_ROOT_PATH = "$(SMCVSTEP_PATH)/../"
# This should point to the directory containing the "smcv" folder:
SMCV_ROOT_PATH = "$(SMCVSTEP_PATH)/../"

CC = g++ -m64
CFLAGS = -c -g -Wall -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) \
                  -I$(CRYSTDIST_ROOT_PATH) -I$(SMCV_ROOT_PATH)

all: bin/smcvstep 

bin/smcvstep: src/smcvstep.o $(MYMOL_LIB_PATH)/mymol.a $(CRYSTDIST_LIB_PATH)/crystdist.a \
                             $(SMCV_LIB_PATH)/smcv.a
	mkdir -p bin
	$(CC) src/smcvstep.o $(MYMOL_LIB_PATH)/mymol.a $(CRYSTDIST_LIB_PATH)/crystdist.a \
                             $(SMCV_LIB_PATH)/smcv.a -o bin/smcvstep

src/smcvstep.o: src/smcvstep.cpp
	$(CC) $(CFLAGS) src/smcvstep.cpp -o src/smcvstep.o

clean:
	rm -rf src/*.o bin
