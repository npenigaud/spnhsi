#!/bin/bash

source /home/afar/modules/use.sh
module load rocm
module load afar/22.2.0

set -x

flang -march=native -O2 -fbackslash -ffp-contract=off -I/home/penigaudn/install-newfxtranacdc/spnhsi_300_brouillon/compile.gpu_afar_d/fxtran-acdc/include -g -fopenmp -fconvert=big-endian -fPIC --offload-arch=gfx942 -lflang_rt.hostdevice   -I/home/afar/software/compilers/afar/rocm-afar-8873-drop-22.2.0/include/hipfort/amdgcn/ -I/home/afar/software/compilers/afar/rocm-afar-8873-drop-22.2.0/include/rocblas/  -L/home/afar/software/compilers/afar/rocm-afar-8873-drop-22.2.0/lib/ -Qunused-arguments -o testbatched.x testbatched.F90 /home/penigaudn/install-newfxtranacdc/spnhsi_300_brouillon/compile.gpu_afar_d/fxtran-acdc/lib/libFxtranACDC.a -lrocblas -lamdhip64 /home/afar/software/compilers/afar/rocm-afar-8873-drop-22.2.0/lib/libhipfort-amdgcn.a -lblas

flang -march=native -O2 -fbackslash -ffp-contract=off -I/home/penigaudn/install-newfxtranacdc/spnhsi_300_brouillon/compile.gpu_afar_d/fxtran-acdc/include -g -fopenmp -fconvert=big-endian -fPIC --offload-arch=gfx942 -lflang_rt.hostdevice   -I/home/afar/software/compilers/afar/rocm-afar-8873-drop-22.2.0/include/hipfort/amdgcn/ -I/home/afar/software/compilers/afar/rocm-afar-8873-drop-22.2.0/include/rocblas/  -L/home/afar/software/compilers/afar/rocm-afar-8873-drop-22.2.0/lib/ -Qunused-arguments -o testbatched_top.x testbatched_top.F90 /home/penigaudn/install-newfxtranacdc/spnhsi_300_brouillon/compile.gpu_afar_d/fxtran-acdc/lib/libFxtranACDC.a -lrocblas -lamdhip64 /home/afar/software/compilers/afar/rocm-afar-8873-drop-22.2.0/lib/libhipfort-amdgcn.a -lblas


