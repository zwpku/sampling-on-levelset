
CFLAGS=-g -O2 -DDEBUG 
CC=g++

LOCAL_DIR=/home/wzhang/local/
#LOCAL_DIR=/home/mi/weizhang/local/
INCLUDES=-I$(LOCAL_DIR)/include 

LDFLAGS=-Wl,-rpath,$(LOCAL_DIR)/lib -L$(LOCAL_DIR)/lib -lm -lpthread -llapack -lblas -lranlib -lconfig++ -lgfortran

SOURCES=ex2_stiff_potential_on_s2.cpp 

all : ex2 gen_path 

ex2 : ex2_stiff_potential_on_s2.cpp
	$(CC) $(CFLAGS) $(INCLUDES)  ex2_stiff_potential_on_s2.cpp -o ex2 $(LDFLAGS)

gen_path : gen_path.cpp
	$(CC) $(CFLAGS) $(INCLUDES) gen_path.cpp -o gen_path $(LDFLAGS)

mcmc_sampler_son : mcmc_sampler_levelset.cpp mcmc.h config_mcmc.cpp
	$(CC) $(CFLAGS) $(INCLUDES) mcmc_sampler_levelset.cpp config_mcmc.cpp -o mcmc_sampler_son $(LDFLAGS)

clean:
	rm -f ex2 gen_path mcmc_sampler_son 

