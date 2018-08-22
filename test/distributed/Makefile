
FLAGS=-debug all
FORT=mpiifort

OBJDIR=build
SRCDIR_FCHL=../../qml/fchl
SRCDIR_MPI=../../qml/distributed

FSOURCES=ffchl_module.f90 ffchl_scalar_kernels.f90 readers.f90 ffchl_reader.f90 ffchl_wrapper.f90
FOBJECTS=$(addprefix $(OBJDIR)/, $(FSOURCES:.f90=.o))

.PHONY: test objects



all: objects test_fchl_gatherv_kernel


test_fchl_gatherv_kernel: build/test_fchl_gatherv_kernel.o test_fchl_gatherv_kernel.exe
	@./test_fchl_gatherv_kernel.exe
	@ echo -n "test_fchl_gatherv_kernel  "
	@./check.py log/test_gatherv_kernel data_check/qm7_fchl_kernel


objects: ${FOBJECTS}

${OBJDIR}/%.o: ${SRCDIR_FCHL}/%.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

${OBJDIR}/%.o: ${SRCDIR_MPI}/%.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

${OBJDIR}/%.o: %.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

%.exe: ${FOBJECTS}
	${FORT} ${FOBJECTS} ${OBJDIR}/${@:.exe=.o} -module ${OBJDIR} -o $@


directories:
	mkdir -p ${OBJDIR} log

clean:
	rm ${OBJDIR}/*.o ${OBJDIR}/*.mod
	rm *.exe




