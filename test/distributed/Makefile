
FORT=mpiifort

FLAGS=-debug all
LDFLAGS=-lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core     -xHost -axAVX -qopenmp

OBJDIR=build
SRCDIR_FCHL=../../qml/fchl
SRCDIR_MPI=../../qml/distributed

FSOURCES=ffchl_module.f90 ffchl_kernels.f90 ffchl_scalar_kernels.f90 readers.f90 printers.f90 ffchl_reader.f90 ffchl_wrapper.f90
FOBJECTS=$(addprefix $(OBJDIR)/, $(FSOURCES:.f90=.o))

# Don't run this makefile in parrallel
.NOTPARALLEL:

# Don't remove object files
.PRECIOUS: ${OBJDIR}/%.o

# These are not real files
.PHONY: test objects



all: objects exes test_fchl_gatherv_kernel test_fchl_qm7_nn_1p

exes: test_fchl_gatherv_kernel.exe test_fchl_qm7_nn_1p.exe

test_fchl_read_job: test_fchl_read_job.exe
	./test_fchl_read_job.exe


test_fchl_gatherv_kernel: test_fchl_gatherv_kernel.exe log
	@./test_fchl_gatherv_kernel.exe
	@ printf "%-40s" "$@"
	@./check.py jobname_dev/_fchl_kernel log/test_gatherv_kernel


test_fchl_qm7_nn_1p: test_fchl_qm7_nn_1p.exe
	@./test_fchl_qm7_nn_1p.exe
	@ printf "%-40s" "$@"
	@./check.py jobname_dev/_fchl_alpha log/test_qm7_nn_1p_alpha



objects: ${OBJDIR} ${FOBJECTS}

${OBJDIR}/%.o: ${SRCDIR_FCHL}/%.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

${OBJDIR}/%.o: ${SRCDIR_MPI}/%.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

${OBJDIR}/%.o: %.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

%.exe: ${OBJDIR}/%.o ${FOBJECTS}
	${FORT} ${LDFLAGS} ${FOBJECTS} ${OBJDIR}/${@:.exe=.o} -module ${OBJDIR} -o $@


${OBJDIR}:
	mkdir -p $@

log:
	mkdir -p $@

clean:
	rm ${OBJDIR}/*.o ${OBJDIR}/*.mod
	rm *.exe

super_clean:
	rm -r log
	rm -r ${OBJDIR}

