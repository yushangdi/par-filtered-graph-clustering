CC = g++
CCFLAGS = -O3  -Wall -std=c++17 -mcx16#-DPROFILE -DPARLAY_SEQUENTIAL# -fcilkplus -DCILKP -DPARLAY_CILK 
PLFLAGS = -ldl -pthread $(SDSLLF) #-fcilkplus 
PARLAY_DIR=../parlaylib/include
BOOST_ROOT=/home/ubuntu/boost_1_78_0

INCLUDE = -I$(PARLAY_DIR) -I$(BOOST_ROOT)

SRC = tmfg.cpp partmfg.cpp partmfg_double.cpp 
NNOBJECTS = $(SRC:.cpp=.o) 

tmfg: $(NNOBJECTS) Makefile
	$(CC) $(PLFLAGS) $(NNOBJECTS) -o $@

tmfg.o: $(SRC) dbht.h linkage.h chain.h partmfg.h
	$(CC) $(INCLUDE) $(CCFLAGS) -c $(SRC)

debug: CCFLAGS= -g -Wall -std=c++17  -I$(PARLAY_DIR) -DDEBUG -mcx16
# -fcilkplus  -DCILKP -DPARLAY_CILK -DPARLAY_USE_STD_ALLOC  -DPARLAY_SEQUENTIAL 

debug: tmfg

clean:
	rm -f tmfg *.o
