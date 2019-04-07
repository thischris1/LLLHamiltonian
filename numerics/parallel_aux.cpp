#ifdef COMPILE_PARALLELIZED

#include<parallel_aux.hpp>

int LLL_MPIStart(int argc,char *argv[])
{
  int numprocs;
  int myid;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  glLogger.redirectCout(numprocs,myid);
  glLogger.info(" - MPI initialized. Running process %d of %d.\n",myid,numprocs);
  return 0;
}

int LLL_MPIEnd()
{
  MPI_Finalize();
  return 0;
}

#endif
