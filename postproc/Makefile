# Objects
#OBJSdef = in_out.o myParaview.o
OBJSdef = myParaview.o
OBJS90 = myParaview90.o
OBJS77 = myParaview77.o
OBJScyl = cyl_paraview.o

# Compiler and flags
FC  = ifort
VEC = -xCORE-AVX2
VEC = -xCORE-AVX2 -axCORE-AVX512
#OPT = -fpe0 -O3 -mkl=parallel -mcmodel=medium -warn all # -diag-file-append=intel_build.log
OPT = -O3 -mkl=parallel -mcmodel=medium -warn all #-diag-file-append=intel_build.log# -warn all
LNK = -shared-intel
LIB =
FLG = $(OPT) $(VEC) $(LNK)

#MKLLIB = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

LIBPLP = -lplplotd -lplplotf77d -lplplotf77cd

.PHONY: clean

# default rule
myParaview: $(OBJSdef)
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) $(OBJSdef) $(LIB) -o ../../bin/$@

pv90: $(OBJS90)
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) $(OBJS90) $(LIB) -o ../../bin/$@

pv77: $(OBJS77)
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) $(OBJS77) $(LIB) -o ../../bin/$@

cyl_paraview: cyl_paraview.o
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) cyl_paraview.o $(LIB) -o ../../bin/$@

# object rules
in_out.o: in_out.F90
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) -c $<

myParaview.o: myParaview.F90
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) -c $<

myParaview90.o: myParaview90.F90
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) -c $<

myParaview77.o: myParaview77.f
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) -c $<

cyl_paraview.o: cyl_paraview.f
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) -c $<

# module rules
%.mod: %.o %.F90
	module purge; \
	module load intel/2019.4; \
	$(FC) $(FLG) -c $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.out
	rm -f *genmod*
	rm -f intel_build.log
