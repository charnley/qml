
TARGET=qmlexe

FLAGS=-debug all
FORT=mpiifort

OBJDIR=build
SRCDIR=../fchl

FSOURCES=ffchl_module.f90 ffchl_scalar_kernels.f90
FOBJECTS=$(addprefix $(OBJDIR)/, $(FSOURCES:.f90=.o))

all: ${TARGET}

${TARGET}: ${FOBJECTS}
	${FORT} ${FOBJECTS} -o $@

${OBJDIR}/%.o: ${SRCDIR}/%.f90
	${FORT} -c ${FLAGS} $< -o $@

${OBJDIR}/%.o: %.f90
	${FORT} -c ${FLAGS} $< -o $@

directories:
	mkdir -p ${BUILD}

clean:
	rm build/*.o

