
CFLAGS=-g -O2 -DDEBUG 
CC=g++

#the directory needs to be provided depending on the environment
LOCAL_DIR=/home/wzhang/local/
#LOCAL_DIR=/home/mi/weizhang/local/

INCLUDES=-I$(LOCAL_DIR)/include 

LDFLAGS=-Wl,-rpath,$(LOCAL_DIR)/lib -L$(LOCAL_DIR)/lib -lm -lpthread -llapack -lblas -lranlib -lconfig++ -lgfortran

SOURCES=ex2_stiff_potential_on_s2.cpp 

all : ex1_gen_path ex2_mcmc ex2_no_mcmc ex3


ex1_gen_path : ex1_gen_path.cpp
	$(CC) $(CFLAGS) $(INCLUDES) ex1_gen_path.cpp -o ex1_gen_path $(LDFLAGS)

ex2_mcmc : ex2_mcmc.cpp ex2.h config_ex2.cpp ex2_utils.cpp ex2_params.cpp
	$(CC) $(CFLAGS) $(INCLUDES) ex2_mcmc.cpp config_ex2.cpp ex2_utils.cpp ex2_params.cpp -o ex2_mcmc $(LDFLAGS)

ex2_no_mcmc: ex2_projection_based_sampler.cpp ex2.h config_ex2.cpp ex2_utils.cpp ex2_params.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) ex2_projection_based_sampler.cpp config_ex2.cpp ex2_utils.cpp ex2_params.cpp -o ex2_no_mcmc $(LDFLAGS)

ex3 : ex3_stiff_potential_on_s2.cpp ex3.h config_ex3.cpp
	$(CC) $(CFLAGS) $(INCLUDES)  ex3_stiff_potential_on_s2.cpp config_ex3.cpp -o ex3 $(LDFLAGS)

clean:
	rm -f ex1_gen_path ex2_mcmc ex2_no_mcmc ex3 

