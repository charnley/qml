
TARGET=qmlexe

FLAGS=-debug all
FORT=mpiifort

OBJDIR=build
SRCDIR=../fchl

FSOURCES=ffchl_module.f90 ffchl_scalar_kernels.f90 readers.f90 ffchl_reader.f90 ffchl_wrapper.f90 ffchl_program.f90
FOBJECTS=$(addprefix $(OBJDIR)/, $(FSOURCES:.f90=.o))

.PHONY: test

all: ${TARGET}

${TARGET}: ${FOBJECTS}
	${FORT} ${FOBJECTS} -module ${OBJDIR} -o $@

${OBJDIR}/%.o: ${SRCDIR}/%.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

${OBJDIR}/%.o: %.f90
	${FORT} -c -module ${OBJDIR} ${FLAGS} $< -o $@

test:
	./qmlexe
	@echo
	./check.py test_ffchl_kernel data/qm7_fchl_kernel

test_mpi:
	mpirun -n 3 --host localhost ./qmlexe
	./check.py test_ffchl_kernel data/qm7_fchl_kernel

directories:
	mkdir -p ${BUILD}

clean:
	rm build/*.o build/*.mod

