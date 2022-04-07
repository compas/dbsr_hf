#  FORTRAN compiler:
FC = gfortran

#  FORTRAN options for compilation and link:
comp = -c -O2
link = -O2

# Binary name
N  = dbsr_hf

# locate your bin directory with executables if desired:
BINDIR = ../

# locate LAPACK and BLAS libraries in your computer:
# many new compilers contain these libraries and you just need to indicate flag "-mkl",
# or to put -llapack -lblas instead the LAPACK and BLAS variables used below

#LIBDIR = $(HOME)/T/LIBS_new
#LAPACK = $(LIBDIR)/lapack.a
#BLAS =   $(LIBDIR)/blas.a
LAPACK = -llapack
BLAS   = -lblas

# source files:
S  =    dbsr_lib_zcom.f90            dbsr_lib_zconfjj.f90      dbsr_lib_dbs.f90       \
  	     hf_MOD_dbsr_hf.f90           hf_MOD_df_orbitals.f90    hf_MOD_energy.f90      \
        dbsr_hf.f90                                                                   \
	     hf_boundary.f90              hf_check_tails.f90        hf_def_conf_jj.f90     \
        hf_def_conf_LS.f90           hf_energy_BR.f90          hf_energy_SE.f90       \
        hf_energy_VP.f90             hf_get_case.f90           hf_get_estimates.f90   \
        hf_get_spline_param.f90      hf_matrix.f90             hf_orthogonality.f90   \
        hf_plot_bsw.f90              hf_read_pqbs.f90          hf_read_GRASP.f90      \
        hf_rk.f90                    hf_rotate.f90             hf_sk.f90              \
        hf_solve_HF.f90              hf_state_energies.f90     hf_summry.f90          \

# replace_iargc.f90   ->  needed for some old compilers only

O = $(S:.f90=.o)

$(N): $(O)
	$(FC) -o $(N) $(O) $(link)  $(LAPACK)  $(BLAS)

clean:
	rm -f *.o *.mod

x:
	rm -f $(N)

.SUFFIXES:
.SUFFIXES: .f90 .mod .o

.f90.o:
	$(FC) $(comp) $<
