EXE_INCLUDEPATH = -I../include -I../include
CFLAGS=-O3 -std=c++0x -m64 -DMKL_ILP64 -fopenmp -DDEBUG


platform=$(shell uname)
ifeq ($(platform),Darwin)
#mac
LIBS =   -L../lib -lcustom_mkl
CFLAGS += -DOSX
CC = g++-5
else

#Linux
LIBS =  -Wl,--no-as-needed  -L../lib -lfastpca_custommkl -liomp5 -ldl -lpthread -lm
CFLAGS+=-DMKL_CUSTOM_LIBRARY='"libfastpca_custommkl.so"'
CC = g++

endif


EXE_CPP_FILES := $(wildcard *.cpp)
EXE_OBJ_FILES := $(addprefix ../build/,$(notdir $(EXE_CPP_FILES:.cpp=.o)))
EXE_TARGET = ../build/fastpca.xx

R_TARGET = ../fastRPCA/inst/libs/libfastpca.so

all: EXE_target R_target

debug: CFLAGS += -DDEBUG
debug: EXE_target R_target

#Cleans all the objects
clean:
	rm -f $(EXE_OBJ_FILES) $(R_OBJ_FILES) $(EXE_TARGET) $(R_TARGET)

#Link the objects into the final EXE executable
EXE_target: $(EXE_OBJ_FILES) 
	$(CC) $(CFLAGS) -o $(EXE_TARGET) $(EXE_OBJ_FILES) $(LIBS)

#Link the objects into the final R library
R_target: $(EXE_OBJ_FILES) 
	$(CC) $(CFLAGS) -o $(R_TARGET) $(EXE_OBJ_FILES) $(LIBS) -shared

#Build the CPP files into objects
../build/%.o: %.cpp
	$(CC) $(CFLAGS) -fPIC -c -o $@ $< $(EXE_INCLUDEPATH) 
