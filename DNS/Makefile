# Objects
OBJS_unifTop         = tools_spectral_cyl.o main_unifTop.o
OBJS_kedgeTop        = tools_spectral_cyl.o main_kedgeTop.o
OBJS_freeSurfTop     = tools_spectral_cyl.o main_freeSurfTop.o
OBJS_kedgeTop2DFD    = tools_FD_cyl.o       main_kedgeTop2DFD.o
OBJS_freeSurfTop2DFD = tools_FD_cyl.o       main_freeSurfTop2DFD.o

# Compiler and flags
FC  = ifort
VEC = -xCORE-AVX2
VEC = -xCORE-AVX2 -axCORE-AVX512
OPT = -fpe0 -O3 -mkl=parallel -mcmodel=medium -warn all #-diag-file-append=intel_build.log
LNK = -shared-intel
LIB =
FLG = $(OPT) $(VEC) $(LNK)

#MKLLIB = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

LIBPLP = -lplplotd -lplplotf77d -lplplotf77cd

.PHONY: clean

# default rule
main_freeSurfTop: $(OBJS_freeSurfTop)
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) $(OBJS_freeSurfTop) $(LIB) -o ../../bin/$@

main_freeSurfTop2DFD: $(OBJS_freeSurfTop2DFD)
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) $(OBJS_freeSurfTop2DFD) $(LIB) -o ../../bin/$@

main_kedgeTop: $(OBJS_kedgeTop)
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) $(OBJS_kedgeTop) $(LIB) -o ../../bin/$@

main_kedgeTop2DFD: $(OBJS_kedgeTop2DFD)
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) $(OBJS_kedgeTop2DFD) $(LIB) -o ../../bin/$@

main_unifTop: $(OBJS_unifTop)
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) $(OBJS_unifTop) $(LIB) -o ../../bin/$@

# object rules
%.o: %.f
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) -c $<

%.o: %.F90
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) -c $<

clean:
	rm -f intel_build.log
	rm -f *.out
	rm -f *.o
	rm -f *.mod
	rm -f *_genmod*
