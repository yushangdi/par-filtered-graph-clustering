CC = g++
CCFLAGS = -O3  -Wall -std=c++17 -mcx16#-DPROFILE -DPARLAY_SEQUENTIAL# -fcilkplus -DCILKP -DPARLAY_CILK 
PLFLAGS = -ldl -pthread $(SDSLLF) #-fcilkplus 
PARLAY_DIR=../../parlaylib/include

INCLUDE = -I$(PARLAY_DIR)

SRC = linkage.cpp
NNOBJECTS = $(SRC:.cpp=.o) 

linkage: $(NNOBJECTS) Makefile
	$(CC) $(PLFLAGS) $(NNOBJECTS) -o $@

linkage.o: $(SRC) linkage.h
	$(CC) $(INCLUDE) $(CCFLAGS) -c $(SRC)

clean:
	rm -f linkage *.o
