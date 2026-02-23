#!/bin/bash

module purge
module load nvhpc/24.5

set -x


nvfortran -I/scratch/work/penigaudn/divers_pr_fev26/spnhsi/compile.gpu_nvhpc_d/fxtran-acdc/include  -acc=gpu -cuda -target=gpu -gopt -O1 -gpu=cc70,cc80 -gopt -lrt -Minfo=accel,all,ccff -o testbatched_top.x testbatched_top.F90 /scratch/work/penigaudn/divers_pr_fev26/spnhsi/compile.gpu_nvhpc_d/fxtran-acdc/lib/libFxtranACDC.a  -L/opt/softs/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/targets/x86_64-linux/lib -cudalib=cublas -Wl,-rpath,/opt/softs/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/targets/x86_64-linux/lib -lblas 

nvfortran -I/scratch/work/penigaudn/divers_pr_fev26/spnhsi/compile.gpu_nvhpc_d/fxtran-acdc/include  -acc=gpu -cuda -target=gpu -gopt -O1 -gpu=cc70,cc80 -gopt -lrt -Minfo=accel,all,ccff -o testbatched.x testbatched.F90 /scratch/work/penigaudn/divers_pr_fev26/spnhsi/compile.gpu_nvhpc_d/fxtran-acdc/lib/libFxtranACDC.a  -L/opt/softs/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/targets/x86_64-linux/lib -cudalib=cublas -Wl,-rpath,/opt/softs/nvidia/hpc_sdk/Linux_x86_64/24.5/math_libs/12.4/targets/x86_64-linux/lib -lblas 
