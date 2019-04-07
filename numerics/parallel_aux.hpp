#ifdef COMPILE_PARALLELIZED

#ifndef PARALLEL_AUX
#define PARALLEL_AUX

/*!
  \file parallel_aux.hpp
  \brief Run this at the very beginning and the very end of a parallelized 
  program which uses MPI.
*/

#include <mpi.h>
#include <logger.hpp>

/*! \fn int LLL_MPIStart(int argc,char *argv[]);
  \param int argc, char *argv[]   Let the main() function of the program
  have the header int main(int argc,char *argv[]) and gimme these two 
  values.
 */

int LLL_MPIStart(int argc,char *argv[]);

int LLL_MPIEnd();


#endif

#endif
