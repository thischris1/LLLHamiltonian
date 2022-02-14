#ifndef TYPES_HPP
#define TYPES_HPP

#include <complex>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>

#include <mkl_lapack.h>


typedef MKL_Complex16 dComplex;
#include <vector>
#include <geometry/CxPeriodicPosition.hpp>
#include <geometry/CxVortex.hpp>
#include <geometry/CxImpurity.hpp>


typedef std::vector<CxPeriodicPosition> posVector ; //! Vector of positions
typedef std::vector<CxPeriodicPosition>::iterator  posVectorIt;
typedef std::vector<CxVortex> vortexVector ; //! Vector of vortices
typedef std::vector<CxVortex>::iterator  vortexIt; //! iterator for vortexVector
typedef std::vector<CxImpurity> impVector;
typedef std::vector<CxImpurity>::iterator impVectorIt;
typedef double MDOUBLE;
typedef  std::complex<MDOUBLE> CMDCOMPLEX;
typedef std::vector<MDOUBLE> DOUBLE_VECTOR;
typedef std::vector<CMDCOMPLEX> COMPLEX_VECTOR;

#endif
