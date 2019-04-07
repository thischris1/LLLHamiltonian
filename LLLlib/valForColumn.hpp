/*!
  \file valForColumn.hpp
*/

#ifndef VALFORCOLUMN_HPP
#define VALFORCOLUMN_HPP

#include <densOperator.hpp>
#include <spinDensOperator.hpp>
#include <spinDownDensOperator.hpp>
#include <spinZDensOperator.hpp>
#include <spinXDensOperator.hpp>
#include <Basis.hpp>
#include <yosBasis.hpp>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <myaccessories.hpp>
//#include <DensOp.hpp>
#include <jxOp.hpp>
#include <jyOp.hpp>


class valForColumn
{
private:
  static const int MAX_DESCR_LEN = 50;

  //DensOperator *densOp;
  DensOperator *pDensOp;
  jxOp *pJxOp;
  jyOp *pJyOp;
  spinDensOperator *spinDensOp;
  spinDownDensOperator *spinDownDensOp;
  spinZDensOperator *spinZDensOp;
  spinXDensOperator *spinXDensOp;

  int whatToCompute;
  yosBasis *yosbas;
  int matType;
  //double *wf;

  char descr[MAX_DESCR_LEN];



public:
  /*!
  \fn valForColumn(int whatToCompute, yosBasis *new_yosbas, int new_matType, eigSt &stateToEval,
  double aTob, double alpha1, double alpha2)
  \brief Constructor.
  \param whatToCompute : Which quantity to cmpute for this column.
  \param new_yosbas : Basis in which to evaluate the Operators
  \param new_matType : Type of matrices to use
  \param stateToEval : For which state should the operators be calculated
  \param aTob : aspect ratio
  \param alpha1 : Flux of solenoid1 in units of h/e
  \param alpha2 : Flux of solenoid2 in units of h/e
  */
  valForColumn(int whatToCompute,
	       yosBasis *new_yosbas,
	       int new_matType,
	       eigSt &stateToEval, 
	       double aTob,
	       double alpha1,
	       double alpha2);

  ~valForColumn();
  double getValForColumn(double x,double y);

  /*!
    \fn void writeLandauDiag (char *strFile)
    \brief writes diagonal elements of Landau-Matrix to file 'strFile'
    \param strFile : Name of file to write to
  */
  void writeLandauDiag (char *strFile);

  /*!
    \fn const char *getDescrOfColumn();
    \brief Returns a tring that describes this column.
  */
  const char *getDescrOfColumn();

};

#endif















