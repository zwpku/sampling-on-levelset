
CFLAGS=-g -O2 -DDEBUG 
CC=g++

LOCAL_DIR=/home/mi/weizhang/local/
INCLUDES=-I$(LOCAL_DIR)/include 

LDFLAGS=-Wl,-rpath,$(LOCAL_DIR)/lib -L$(LOCAL_DIR)/lib -lm -lpthread

SOURCES=ex2_stiff_potential_on_s2.cpp 

all : ex2 gen_path 

ex2 : ex2_stiff_potential_on_s2.cpp
	$(CC) $(CFLAGS) ex2_stiff_potential_on_s2.cpp com.c linpack.c ranlib.c -o ex2

gen_path : 
	$(CC) $^ gen_path.cpp com.c linpack.c ranlib.c -o gen_path

clean:
	rm -f ex2

