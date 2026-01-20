
# Note: Make mymol and crystdist before making copcalc

PLUMED2OPD_PATH = `pwd`
# This should point to the directory containing the mymol library (mymol.a):
MYMOL_LIB_PATH = ../mymol/lib
# This should point to the directory containing the crystdist library (crystdist.a):
CRYSTDIST_LIB_PATH = ../crystdist/lib
# This should point to the directory containing the "parser" folder:
PARSER_LIB_PATH = ../parser/lib
# This should point to the directory containing the "common" folder:
COMMON_ROOT_PATH = "$(PLUMED2OPD_PATH)/../"
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(PLUMED2OPD_PATH)/../"
# This should point to the directory containing the "mymol" folder:
MYMOL_ROOT_PATH = "$(PLUMED2OPD_PATH)/../"
# This should point to the directory containing the "crystdist" folder:
CRYSTDIST_ROOT_PATH = "$(PLUMED2OPD_PATH)/../"
PARSER_ROOT_PATH = "$(PLUMED2OPD_PATH)/../"
# This should point to the dir


CC = g++ -m64
CFLAGS = -c -g -I$(COMMON_ROOT_PATH) -I$(MYMATH_ROOT_PATH) \
                  -I$(MYMOL_ROOT_PATH) -I$(CRYSTDIST_ROOT_PATH) \
                   -I$(PARSER_ROOT_PATH)
all: bin/plumed2opd

bin/plumed2opd: src/plumed2opd.o  $(PARSER_LIB_PATH)/parser.a
	mkdir -p bin
	$(CC)   src/plumed2opd.o  $(PARSER_LIB_PATH)/parser.a -o bin/plumed2opd



src/plumed2opd.o: src/plumed2opd.cpp
	$(CC) $(CFLAGS) src/plumed2opd.cpp -o src/plumed2opd.o




clean:
	rm -rf src/*.o bin
