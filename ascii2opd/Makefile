
# Note: Make mymol and crystdist before making copcalc

ASCII_PATH = `pwd`
# This should point to the directory containing the mymol library (mymol.a):
PARSER_LIB_PATH = ../parser/lib
# This should point to smcv root directory
SMCV_LIB_PATH = ../smcv/lib
# This should point to the directory containing the "mymath" folder:
MYMATH_ROOT_PATH = "$(ASCII_PATH)/../"
# This should point to directory containing the Parser folder
PARSER_ROOT_PATH = "$(ASCII_PATH)/../"
SMCV_ROOT_PATH = "$(ASCII_PATH)../"

CC = g++ -m64
CFLAGS = -c -g   -I$(MYMATH_ROOT_PATH) \
                  -I$(SMCV_ROOT_PATH)  \
                   -I$(PARSER_ROOT_PATH)
all: bin/ascii2opd

bin/ascii2opd: src/ascii2opd.o $(SMCV_LIB_PATH)/smcv.a $(PARSER_LIB_PATH)/parser.a
	mkdir -p bin
	$(CC)   src/ascii2opd.o $(PARSER_LIB_PATH)/parser.a $(SMCV_LIB_PATH)/smcv.a -o bin/ascii2opd

src/ascii2opd.o: src/ascii2opd.cpp
	$(CC) $(CFLAGS) src/ascii2opd.cpp -o src/ascii2opd.o

clean:
	rm -rf src/*.o bin
