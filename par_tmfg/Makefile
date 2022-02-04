CC = g++
CCFLAGS = -O3  -Wall -std=c++17 -mcx16#-DPROFILE -DPARLAY_SEQUENTIAL# -fcilkplus -DCILKP -DPARLAY_CILK 
PLFLAGS = -ldl -pthread $(SDSLLF) #-fcilkplus 
PARLAY_DIR=../parlaylib/include
BOOST_ROOT=./boost_1_77_0

INCLUDE = -I$(PARLAY_DIR) -I$(BOOST_ROOT)

SRC = tmfg.cpp partmfg.cpp partmfg_double.cpp 
NNOBJECTS = $(SRC:.cpp=.o) 

tmfg: $(NNOBJECTS) Makefile
	$(CC) $(PLFLAGS) $(NNOBJECTS) -o $@

tmfg.o: $(SRC) dbht.h linkage.h chain.h partmfg.h
	$(CC) $(INCLUDE) $(CCFLAGS) -c $(SRC)

clean:
	rm -f tmfg *.o
