rm -rf  *.o
sw5cc -host -OPT:IEEE_arith=1 -I/usr/sw-mpp/mpi2/include -std=c99 -c master.c
sw5cc -slave -OPT:IEEE_arith=1 -std=c99 -c slave.c
mpicc master.o slave.o -lstdc++ -OPT:IEEE_arith=1 -o a.out
