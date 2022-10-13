#Makefile for parallel compiling
#=====================================================
# Compiler options
#============================================================================
# IFORT
#FC = h5pfc -r8 -ip -ipo -O3 -fpp -xHost  -qno-openmp -g -traceback
FC = h5pfc -r8 -O3 -qno-openmp
#FC = h5pfc -r8 -O0  -check all -g -traceback

#FC += -openmp
#FC += -debug all -warn all -check all -g -traceback
#FC += -fpe0 -ftrapuv

# CRAY

#=======================================================================
# Library
#======================================================================
#LINKS = /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3 \
	-L$(FFTWLIB) -lfftw3 -lfftw3f
#LINKS = -L$(FFTWLIB) -lfftw3 -L$(MKLLIB) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  

#LINKS = -lfftw3 -lblas -llapack

#LINKS = -lfftw3 -llapack -lessl 

#LINKS = -L/usr/local/lib -llapack -lblas -lfftw3 -lfftw3_threads

LINKSAFTER = -lfftw3 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lhdf5_fortran -lhdf5  -lz -ldl -lm 

#LINKS = -L${FFTW_LIB} -lfftw3 -L${LAPACK_LIB} -llapack -L${ESSL_LIB} -lesslbg -L${BLAS_LIB} -lblas
#============================================================================ 
# make PROGRAM   
#============================================================================

PROGRAM = boutnp

OBJECTS = CalcMaxCFL.o CreateGrid.o SetTempBCs.o \
          CalcLocalDivergence.o CheckDivergence.o gcurv.o ExplicitTermsVX.o ExplicitTermsVY.o ExplicitTermsVZ.o \
          ExplicitTermsTemp.o ReadInitCond.o CreateInitialConditions.o \
          ImplicitAndUpdateVX.o ImplicitAndUpdateVY.o ImplicitAndUpdateVZ.o ImplicitAndUpdateTemp.o \
          OpenLogs.o MpiAuxRoutines.o \
          matrix_transpose.o mpi_routines.o AuxiliaryRoutines.o InitTimeMarcherScheme.o \
          papero.o param.o CorrectPressure.o InitPointP.o \
          SolveImpVXY_Z.o SolveImpVXYZ_X.o SolveImpVXYZ_Y.o CalcPointPVel.o UpdatePointParticlePosition.o \
          StatRoutines.o PointPartAuxRoutines.o PointPartIO.o DeallocateArrays.o \
          SolveImpVZ_Z.o SolveImpTemp_Z.o  DefineIBMObject.o ReadPPartInput.o \
          tripvmy_line.o TimeMarcher.o CorrectVelocity.o GlobalQuantities.o InitPressureSolver.o Dump3DFiles.o \
          SolvePressureCorrection.o CalcDissipationNu.o InterpolateGrid.o LocateLargeDivergence.o CalcPlateNu.o \
	  CalcVorticity.o CalcMaterialDerivative.o InitArrays.o ReadInputFile.o

MODULES = param.o
#============================================================================ 
# Linking    
#============================================================================

$(PROGRAM) : $(MODULES) $(OBJECTS)
	$(FC) $(OBJECTS) -o $@ $(LINKSAFTER) 

#============================================================================
#  Dependencies
#============================================================================

param.o: param.F90
	$(FC) -c $(OP_COMP) param.F90

%.o:   %.F $(MODULES)
	$(FC) -c $(OP_COMP) $< 

%.o:   %.F90 $(MODULES)
	$(FC) -c $(OP_COMP) $< 
        
        

#============================================================================
#  Clean up
#============================================================================

clean :
	rm *.o 
	rm *.mod

