NVCC=nvcc

NVCC_FLAGS=-O3 -arch=sm_70

NVCC_LIBS=-lcublas -lm

 
### compile
GPU_CG: GPU_CG.cu

        ${NVCC} ${NVCC_FLAGS} -o $@ $^ ${NVCC_LIBS}

 
### clean executables
clean:

        $(RM) GPU_CG
