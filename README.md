Taisa Kushner   
contact:    taisa.kushner@colorado.edu



-------------------  NEURAL NETWORKS  ------------------
code for learning neural networks to predict quantiles - CMSB paper (within LearnNetwork folder)
example command line arg to run:

python3 learnGlucoseInsulin.py TensorFlowInput_DallaMan_ptid1_first250.csv 37 37 DallaMan_pt1

Additional Information:
input training data, number of glucose and insulin inputs, output name
The format of the CSV file should be as follows
No header, use commas to separate values
Insulin(t-180),  Insulin(t-175), ... ,Insulin(t), Glucose(t-180), .. Glucose(t), Glucose(t+Tg)

dumptoSherlock.m -- matlab script to convert learned network into format required for use in MPC by Sherlock tool.



-------------------  MPC ------------------
MPC code (within MPC_Reachability folder)

utlizes four networks - mean, high, low as well as plant
all must be trained independently of one another.  -- see examples in network_files
Format required as specified in sherlock-network-format.pdf
dumptoSherlock.m will convert tensorflow model to this style


Instructions to Compile:
Please modify the file Makefile.locale to help find Gurobi.
For a Mac with Gurobi 8.1: 

> HOST_ARCH=mac64
> GUROBI_PATH=/Library/gurobi810

You should feel free to modify these two variables. The Makefile will look for Gurobi headers under

> $(GUROBI_PATH)/$(HOST_ARCH)/include

and libraries under

> $(GUROBI_PATH)/$(HOST_ARCH)/include


Once these are set, you should type

> make 

to compile.

## Instructions to run -- will run with defaults in network_files folder

> ./run_file 


Additional specifications to be made, if desired:

-predict_network                   mean network file name
-upper_bound_network           	   upper networkm file name
-lower_bound_network               lower network file name
-sim_network                       plant network file name
-init_trace_file                   initial trace to seed -- see example input_1.csv
-mean_gluc_val                     control mean - eg 140
-upper_gluc_val                    bound upper network - eg 210
"-lower_gluc_val                   bound lower network - eg 70
-total_time_steps                  total time steps - default 100
-mpc_time                          MPC horizon, in minutes, default 60min
-apply_bounds                      binary, apply the network bounds, default true
-min_insulin_level                 eg. 0.0
-insulin_granularity               eg 0.001
-total_insulin                     total insulin allowable over time frame
-verbosity                         binary, display output, default true
-O                                 Output file base name



-------------------  PLEASE REFERENCE ------------------
S. Dutta, T. Kushner, S. Sankaranarayanan, Robust Data-Driven Control of Artificial
Pancreas Systems using Neural Networks. Computational Methods in Systems Biology.
(CMSB). Lecture Notes in Computer Science, vol 11095. (2018) Springer

Souradeep Dutta, Susmit Jha, Sriram Sankaranarayanan, and Ashish Tiwari, 
Output Range Analysis for Deep Feedforward Neural Networks 
In Proceedings of NASA Formal Methods Symposium (NFM), Volume 10811
 of Lecture Notes In Computer Science pp. 121-138 (2018). 
