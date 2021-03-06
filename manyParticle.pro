######################################################################
# Automatically generated by qmake (2.01a) Mo. Nov 6 21:23:16 2017
######################################################################

TEMPLATE = lib
TARGET += parallel 
CONFIG -= qt
CONFIG += dll
CONFIG += debug_and_release
DEFINES +=  MKL MKL_ILP64 GSL_DLL
QMAKE_LIBDIR= /usr/local/lib /opt/intel/mkl/lib/intel64
DEPENDPATH = . geometry LLLlib numerics utils
INCLUDEPATH = . geometry utils LLLlib numerics  /opt/intel/mkl/include /home/chris/include
LIBS += -L/usr/local/lib -L/opt/intel/mkl/lib/intel64 -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread  -ldl 

# Input
HEADERS += geometry/C2DVector.hpp \
           geometry/CDoubleGrid.h \
           geometry/CGaussianImpurityArray.h \
           geometry/CGrid.h \
           geometry/CImpurityArray.h \
           geometry/CImpurityArrayValueProvider.h \
           geometry/CIntegerGrid.h \
           geometry/CScalarGrid.h \
           geometry/CVectorGrid.h \
           geometry/CxImpurity.hpp \
           geometry/CxPath.hpp \
           geometry/CxPeriodicPosition.hpp \
           geometry/CxPosition.hpp \
           geometry/CxPositionArray.h \
           geometry/CxVortex.hpp \
           geometry/CxVortexSystem.hpp \
           LLLlib/BarrYEl.hpp \
           LLLlib/Basis.hpp \
           LLLlib/densOperator.hpp \
           LLLlib/eigSt.hpp \
           LLLlib/ERRORS.h \
           LLLlib/jxOp.hpp \
           LLLlib/jyOp.hpp \
           LLLlib/LLLhamBarrier.hpp \
           LLLlib/LLLhamGaussian.hpp \
           LLLlib/LLLhamGaussianArray.h \
           LLLlib/LLLhamiltonian.hpp \
           LLLlib/LLLhamInplaneMagImp.hpp \
           LLLlib/LLLhamInplaneMagImpRectg.hpp \
           LLLlib/LLLhamMagImp.hpp \
           LLLlib/LLLhamMagImpY.hpp \
           LLLlib/LLLlib.h \
           LLLlib/LLLRandomHamiltonian.h \
           LLLlib/LLLRandomLandauMatrix.h \
           LLLlib/mpi_tag_code.hpp \
           LLLlib/myaccessories.hpp \
           LLLlib/OneParticleOperator.hpp \
           LLLlib/Operator.hpp \
           LLLlib/persist_eigSt.hpp \
           LLLlib/sort_eigSt.hpp \
           LLLlib/spinDensOperator.hpp \
           LLLlib/spinDownDensOperator.hpp \
           LLLlib/spinSquareOperator.hpp \
           LLLlib/spinXDensOperator.hpp \
           LLLlib/spinXOperator.hpp \
           LLLlib/spinZDensOperator.hpp \
           LLLlib/spinZOperator.hpp \
           LLLlib/State.hpp \
           LLLlib/SymmetricMatrix.h \
           LLLlib/th1_natureConst.hpp \
           LLLlib/TPositionArray.h \
           LLLlib/Types.hpp \
           LLLlib/valForColumn.hpp \
           LLLlib/yosBasis.hpp \
           LLLlib/yosEigenState.hpp \
           LLLlib/yosState.hpp \
           numerics/CGslInterpolator.h \
           numerics/ComplexSquareMatrix.h \
           numerics/HermitMatrix.h \
           numerics/Matrix.h \
           numerics/matrix_full_cplx.hpp \
           numerics/matrix_full_real.hpp \
           numerics/Matrix_sparse_complex_parallel.h \
           numerics/Matrix_sparse_cplx.hpp \
           numerics/Matrix_sparse_real.hpp \
           numerics/Matrix_sparse_real_parallel.hpp \
           numerics/parallel_aux.hpp \
           numerics/RealSparseSquareMatrix.h \
           numerics/RealSquareMatrix.h \
           numerics/SquareMatrix.h \
           utils/C2dDistribution.hpp \
           utils/CDistanceDistribution.h \
           utils/CDistanceDistributionHR.h \
           utils/CDistribution.hpp \
           utils/CFileParser.h \
           utils/CLine.h \
           utils/ConstantValueProvider.h \
           utils/CRandomizer.h \
           utils/CRandomizerGsl.h \
           utils/CRandomizerMkl.h \
           utils/CSingleImpurityProvider.h \
           utils/CxBadValueError.hpp \
           utils/CxErrors.hpp \
           utils/CxFileNotFoundError.hpp \
           utils/CxIndextoLargeError.hpp \
           utils/CxIndextoSmallError.hpp \
           utils/CxNullPointerError.hpp \
           utils/CxOutOfBoundsError.hpp \
           utils/CxVortexSystem.hpp \
           utils/IInterpolator.h \
           utils/IValueProvider.h \
           utils/logger.hpp \
           utils/types.h
SOURCES += geometry/CDoubleGrid.cpp \
           geometry/CGaussianImpurityArray.cpp \
           geometry/CGrid.cpp \
           geometry/CImpurityArray.cpp \
           geometry/CImpurityArrayValueProvider.cpp \
           geometry/CIntegerGrid.cpp \
           geometry/CScalarGrid.cpp \
           geometry/CVectorGrid.cpp \
           geometry/CxImpurity.cpp \
           geometry/CxPath.cpp \
           geometry/CxPeriodicPosition.cpp \
           geometry/CxPosition.cpp \
           geometry/CxPositionArray.cpp \
           geometry/CxVortex.cpp \
           LLLlib/BarrYEl.cpp \
           LLLlib/Basis.cpp \
           LLLlib/densOperator.cpp \
           LLLlib/eigSt.cpp \
           LLLlib/jxOp.cpp \
           LLLlib/jyOp.cpp \
           LLLlib/LLLhamBarrier.cpp \
           LLLlib/LLLhamGaussian.cpp \
           LLLlib/LLLhamGaussianArray.cpp \
           LLLlib/LLLhamiltonian.cpp \
           LLLlib/LLLhamInplaneMagImp.cpp \
           LLLlib/LLLhamInplaneMagImpRectg.cpp \
           LLLlib/LLLhamMagImp.cpp \
           LLLlib/LLLhamMagImpY.cpp \
           LLLlib/LLLRandomHamiltonian.cpp \
           LLLlib/LLLRandomLandauMatrix.cpp \
           LLLlib/myaccessories.cpp \
           LLLlib/OneParticleOperator.cpp \
           LLLlib/Operator.cpp \
           LLLlib/persist_eigSt.cpp \
           LLLlib/sort_eigSt.cpp \
           LLLlib/spinDensOperator.cpp \
           LLLlib/spinDownDensOperator.cpp \
           LLLlib/spinSquareOperator.cpp \
           LLLlib/spinXDensOperator.cpp \
           LLLlib/spinXOperator.cpp \
           LLLlib/spinZDensOperator.cpp \
           LLLlib/spinZOperator.cpp \
           LLLlib/valForColumn.cpp \
           LLLlib/yosBasis.cpp \
           LLLlib/yosEigenState.cpp \
           LLLlib/yosState.cpp \
           numerics/CGslInterpolator.cpp \
           numerics/ComplexSquareMatrix.cpp \
           numerics/HermitMatrix.cpp \
           numerics/Matrix.cpp \
           numerics/matrix_full_cplx.cpp \
           numerics/matrix_full_real.cpp \
           numerics/Matrix_sparse_complex_parallel.cpp \
           numerics/Matrix_sparse_cplx.cpp \
           numerics/Matrix_sparse_real.cpp \
           numerics/Matrix_sparse_real_parallel.cpp \
           numerics/parallel_aux.cpp \
           numerics/RealSparseSquareMatrix.cpp \
           numerics/RealSquareMatrix.cpp \
           numerics/SquareMatrix.cpp \
           utils/C2dDistribution.cpp \
           utils/CDistanceDistribution.cpp \
           utils/CDistanceDistributionHR.cpp \
           utils/CDistribution.cpp \
           utils/CFileParser.cpp \
           utils/CLine.cpp \
           utils/ConstantValueProvider.cpp \
           utils/CRandomizer.cpp \
           utils/CRandomizerGsl.cpp \
           utils/CRandomizerMkl.cpp \
           utils/CSingleImpurityProvider.cpp \
           utils/CxBadValueError.cpp \
           utils/CxErrors.cpp \
           utils/CxFileNotFoundError.cpp \
           utils/CxIndextoLargeError.cpp \
           utils/CxIndextoSmallError.cpp \
           utils/CxNullPointerError.cpp \
           utils/CxOutOfBoundsError.cpp \
           utils/IInterpolator.cpp \
           utils/IValueProvider.cpp \
           utils/logger.cpp
