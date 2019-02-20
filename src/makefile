
CFLAGS=-g -O2 -DDEBUG 
CC=g++

LOCAL_DIR=/home/wzhang/local/
#LOCAL_DIR=/home/mi/weizhang/local/
INCLUDES=-I$(LOCAL_DIR)/include 

LDFLAGS=-Wl,-rpath,$(LOCAL_DIR)/lib -L$(LOCAL_DIR)/lib -lm -lpthread -llapack -lblas -lranlib -lconfig++ -lgfortran

SOURCES=ex2_stiff_potential_on_s2.cpp 

all : ex2 gen_path mcmc_sampler_son ex4_no_mcmc 

ex2 : ex2_stiff_potential_on_s2.cpp ex2.h config_ex2.cpp
	$(CC) $(CFLAGS) $(INCLUDES)  ex2_stiff_potential_on_s2.cpp config_ex2.cpp -o ex2 $(LDFLAGS)

gen_path : gen_path.cpp
	$(CC) $(CFLAGS) $(INCLUDES) gen_path.cpp -o gen_path $(LDFLAGS)

mcmc_sampler_son : mcmc_sampler_levelset.cpp ex4.h config_ex4.cpp ex4_utils.cpp ex4_params.cpp
	$(CC) $(CFLAGS) $(INCLUDES) mcmc_sampler_levelset.cpp config_ex4.cpp ex4_utils.cpp ex4_params.cpp -o mcmc_sampler_son $(LDFLAGS)

ex4_no_mcmc: ex4_projection_based_sampler.cpp ex4.h config_ex4.cpp ex4_utils.cpp ex4_params.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) ex4_projection_based_sampler.cpp config_ex4.cpp ex4_utils.cpp ex4_params.cpp -o ex4_no_mcmc $(LDFLAGS)

clean:
	rm -f ex2 gen_path mcmc_sampler_son ex4_no_mcmc

