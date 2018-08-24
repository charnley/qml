
FLAGS=-debug all
FORT=mpiifort

OBJDIR=build
SRCDIR_FCHL=../../qml/fchl
SRCDIR_MPI=../../qml/distributed

FSOURCES=ffchl_module.f90 ffchl_scalar_kernels.f90 readers.f90 ffchl_reader.f90 ffchl_wrapper.f90
FOBJECTS=$(addprefix $(OBJDIR)/, $(FSOURCES:.f90=.o))

# Don't run this makefile in parrallel
.NOTPARALLEL:

# Don't remove object files
.PRECIOUS: ${OBJDIR}/%.o

# These are not real files
.PHONY: test objects



all: objects test_fchl_gatherv_kernel test_fchl_guido test_fchl_gatherv_kernel_x


test_fchl_gatherv_kernel: test_fchl_gatherv_kernel.exe log
	@ printf "%-40s" "$@"
	@./test_fchl_gatherv_kernel.exe
	@./check.py log/test_gatherv_kernel data_check/qm7_fchl_kernel


test_fchl_gatherv_kernel_x: test_fchl_gatherv_kernel.exe
	@ printf "%-40s" "$@"
	@./test_fchl_gatherv_kernel.exe
	@ echo "no" >> log/test_gatherv_kernel
	@./check.py log/test_gatherv_kernel data_check/qm7_fchl_kernel


test_fchl_guido:
	@ printf "%-40s" "$@"
	@ echo "\033[32mPASSED\033[0m"




objects: ${OBJDIR} ${FOBJECTS}

${OBJDIR}/%.o: ${SRCDIR_FCHL}/%.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

${OBJDIR}/%.o: ${SRCDIR_MPI}/%.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

${OBJDIR}/%.o: %.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

%.exe: ${OBJDIR}/%.o ${FOBJECTS}
	${FORT} ${FOBJECTS} ${OBJDIR}/${@:.exe=.o} -module ${OBJDIR} -o $@


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

