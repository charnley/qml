
# Quantum Machine Learning (QML) using Message Passing Interface (MPI)

This is the drivers for distributing the QML methods over a cluster.
Playing with MPI for now, more will come.

## Dependencies

    fortran compiler
    MPI
    and distributed solver stuff


## Install and setup

### Using intels compilerset

    mpiifort -c module.f90 -o module.o

### Using gnu compilerset

    gfortran -c -o $@ $<
    mpifort -c -o $@ $<

TODO add some flags



