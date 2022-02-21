/*!
  \file Matrix_sparse_cplx.cpp
 \brief Implementation of Class Matrix_sparse_cplx for sparse matrix diagonalization
 \author Christian Mueller
*/
#include <numerics/Matrix_sparse_cplx.hpp>
#include <stdio.h>
#include <LLLlib/myaccessories.hpp>
#include <utils/logger.hpp>
#include <iostream>
#ifdef MKL
#include <mkl.h>
#endif
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
/*
  static members
 */
const double Matrix_sparse_cplx::DIAG_TOLERANCE=1.e-4; //! precision at diagonalization
const int Matrix_sparse_cplx::ADD_TO_DIAG_WORKSPACE_SPARSE = 140; //! NAG erfahrungswert
const int Matrix_sparse_cplx::DIAG_MAX_ITERATIONS = 10000; //! according to mail from karel 

/*************************************
 Functions related to the class (and to f02fjf_ of the NAG library)
**************************************/

// !!!!!!!!! See comments in the .hpp file !!!!!!!!!!!!
/*!
  \fn void MatrixDotVectorWrapper(int *iflag,int *n,double *z, double w[],double rwork[],int *lrwork,int iwork[],int *liwork)
*/
void MatrixDotVectorWrapper(int *iflag,int *n,double *z, double w[],double rwork[],int *lrwork,int iwork[],int *liwork)
{
  //  Matrix_sparse_cplx * tempPtr = dynamic_cast< Matrix_sparse_cplx *> (iwork);
  //  tempPtr->matrixDotVector(iflag,z,w);
    ((Matrix_sparse_cplx *) iwork)->matrixDotVector(iflag,z,w);
}


/*!
  \brief Compute scalar product of two vectors 
*/
double VecDotVec(int *iflag,int *n,double z[],double w[],double rwork[],int *lrwork,int iwork[],int *liwork)
{
  double sum=0;
  for(int i=0;i<*n;i++) sum+=z[i]*w[i];
  return sum;
  // A slower alternative...
  //return ((Matrix_sparse_cplx *) iwork)->dot(iflag,z,w);
}

/*!
  \brief Monitoring function called bz=y NAG routine
*/
void GoToMonitC(int *istate,int *nextit,int *nevals,int *nevecs,int *k,double f[],double d[])
{
 
  std::cerr << "Diagonalization progress for complex sparse matrix receives "<<	*istate << " , " << *nextit;
  std::cerr  <<" , "<< *nevals << ", " << *nevecs << " , " << *k << std::endl;
  return;

}


/**********************************************
 Methods
***********************************************/
/*!
  \fn Matrix_sparse_cplx::Matrix_sparse_cplx(int *new_row_NrNz, int *new_col_ind, std::complex <double> *new_mat_sp, int new_nz_elements, int new_n,int new_how_many_vects, int add_to_workspace):
  \brief standard constructor


*/




Matrix_sparse_cplx::Matrix_sparse_cplx(vector<MKL_INT> new_row_NrNz,
				       vector <MKL_INT> new_col_ind,
				       vector < std::complex <double> > new_mat_sp,
				       int new_nz_elements,
				       int new_n,
				       int new_how_many_vects,
				       int add_to_workspace):
  row_NrNz(new_row_NrNz),
  mat_sp(new_mat_sp),
  nz_elements(new_nz_elements), // copies the elements deep - make it a shallow copy
  dimension(new_n),
  how_many_eigvects(new_how_many_vects),
  m_m0(1.0),
  add_to_diag_workspace_sparse(add_to_workspace),
  offSetDiagonal(-5.0),
  emin(-3.0),
  emax(-1.0),
  useRCI(false)
{
  glLogger.info("Entering ctor of sparse cplx matrix");
  col_ind = new_col_ind;
  // look for emin emax file, just 2 lines with doubles in it
  std::string theLine;
  	std::ifstream energyLimitFile("energyLimits.dat");
  	if (!energyLimitFile.good())
  	{
  		ERROR(" Energy limit file not found");
  		emin = -3;
  		useRCI=false;

  	}
  	else
  	{
  		int lineCount = 0;
  		useRCI =true;
  		ERROR("Reading Energy Limit file")
  		while (std::getline(energyLimitFile,theLine ))
  			{
  				if (istarts_with(theLine,"#"))
  				{

  					continue;
  				}
  				if (lineCount == 0)
  				{
  					emin = boost::lexical_cast<double>(theLine);
  				}
				if (lineCount == 1)
  				
  				{
  					emax = boost::lexical_cast<double>(theLine);
  				}
				if (lineCount ==2)
				  {
				    m_m0 =  boost::lexical_cast<double>(theLine);
				  }
  				lineCount++;
  	}
  	}

}



/*!  \fn Matrix_sparse_cplx::~Matrix_sparse_cplx()
  \brief destrcutor
*/
Matrix_sparse_cplx::~Matrix_sparse_cplx()
{


}

/*!
  \brief Accessor 
  \fn void Matrix_sparse_cplx::set_k(int k)
\param k 

*/

void Matrix_sparse_cplx::set_k(int k)
{
  add_to_diag_workspace_sparse=k;
}


/*!
  \fn int Matrix_sparse_cplx::foomax(int a,int b)

*/
int Matrix_sparse_cplx::foomax(int a,int b)
{
  if(a>b) return a;
  else return b;
}


/* Compute y=A*x; a slight trouble is that only one half (upper, say) of the matrix
   is stored (say A=U+D+L=U+D+U^T; stored is just D,U). Thus, to get Ax you have
   to compute U*x and U^T*x.
   Written according to
http://www.netlib.org/linalg/html_templates/node98.html#SECTION00932100000000000000
*/
/*!
 \fn void Matrix_sparse_cplx::matrixDotVector(int *iflag,std::complex<double> *z,double w[])

\brief Calculates the Product Matrix * vector where the input vector is in inVector and the result vector is in outVector (this*inVector = outVector)
USes member variable calculateReal to determine which part is to be calculated
*/

void Matrix_sparse_cplx::matrixDotVector(int *iflag,
					 double *inVector,
					 double *outVector)
{
  glLogger.debug("Entering matrixDot Vector");

  int i_Nz=0;
  std::complex<double> matEl;
  int 
    cInd = 0,
    elsInRow = 0;
  /* 
     Set result to 0
  */
  for(int index=0;
      index < dimension*2;
      index++)
    {
      /*
	Initialize vector
      */
      outVector[index]= 0.0;  
      glLogger.debug("InVector [%d] = (%f)", index, inVector[index]);
    }


  
  for(int rowIndex=0;
      rowIndex < dimension;
      rowIndex++)
    {
      /*
	Loop over rows in complex matrix 
      */

      elsInRow=row_NrNz[rowIndex]; // Number of non-zero elements in row [rowIndex] of the complex matrix 
      glLogger.debug("Row no (%d) has (%d) non zero elements", rowIndex, elsInRow);
      /* D*x would be counted twice, once here and once at U^T*x, otherwise */    

      // First element is the diagonalterm
      // Contribution of diagonal matrix element

	  outVector[rowIndex]= outVector[rowIndex]+ 
	    (mat_sp[i_Nz].real()+offSetDiagonal)*inVector[rowIndex] - 
	    mat_sp[i_Nz].imag()*inVector[rowIndex+dimension];
	  
	  outVector[rowIndex+dimension] = outVector[rowIndex+dimension]+  
	    mat_sp[i_Nz].imag()*inVector[rowIndex] + 
	    (mat_sp[i_Nz].real()+offSetDiagonal)*inVector[rowIndex+dimension];
	  
	  glLogger.debug("Diagonalelement contribution to (%d) = (%f)", 
			 rowIndex, outVector[rowIndex]);
	  glLogger.debug("2. Diagonalelement contrbution to (%d). element  = (%f)",
			 rowIndex+dimension, outVector[rowIndex +dimension]);
	  i_Nz++;
	  //	}
      for(int colIndex=i_Nz; // Start at begin of a row (in internal representation)
	  colIndex < i_Nz+elsInRow-1;
	  colIndex++) // Loop over elements in row 
	{
	  glLogger.debug("Loop over columns in matrix col(%d)", colIndex);
	  /*
	    Loop over columns in matrix
	  */
	  cInd=col_ind[colIndex]; // The index of non-zero element
	  std::complex<double> matrixElement = mat_sp[colIndex];
	  glLogger.debug("The index of the non zero element is (%d) value is (%f), (%f)", 
			cInd, matrixElement.real(), matrixElement.imag());
	  
	  /*
	    Unmittelbare Elemente berechnen
	  */
	  // erste Reihe
	  outVector[rowIndex] =  outVector[rowIndex] +
	    //	    matrixElement.real()*inVector[rowIndex] -  matrixElement.imag()*inVector[rowIndex + dimension];
		    matrixElement.real()*inVector[cInd] -  matrixElement.imag()*inVector[cInd + dimension];
	  // Zweite Reihe
	  outVector[rowIndex+dimension] =  outVector[rowIndex+dimension] +
	    //	    matrixElement.imag() *inVector[rowIndex] +  matrixElement.real()*inVector[rowIndex + dimension];
	  matrixElement.imag() *inVector[cInd] +  matrixElement.real()*inVector[cInd + dimension];
	  glLogger.debug("invector elements used (%d)=(%f), (%d)=(%f)", 
			 cInd, inVector[cInd], 
			 cInd+dimension, inVector[cInd+dimension]);
	  glLogger.debug("Offdiagonal part, after calculating contribution to element (%d) =(%f) and (%d) = (%f)",
			  rowIndex, outVector[rowIndex], rowIndex+dimension, outVector[rowIndex+dimension]);

	  /*
	    Symmetrische matrixbeitraege fuer andere vectorkomponenten berechnen
	  */
  glLogger.debug("Now calculating contribution to other components of vector due to symmetry of matrix");
	  outVector[cInd]= outVector[cInd] +
	    matrixElement.real()*inVector[rowIndex] + matrixElement.imag()*inVector[rowIndex+dimension];
	  outVector[cInd+dimension] = outVector[cInd+dimension] -
	     matrixElement.imag()*inVector[rowIndex] + matrixElement.real()*inVector[rowIndex+dimension];

	  glLogger.debug("Use invector components (%d)=(%f), (%d)=(%f)", 
		rowIndex, inVector[rowIndex], rowIndex+dimension, inVector[rowIndex+dimension]);
		 
	  glLogger.debug("Calculated contributions to elements (%d)=(%f) and (%d)=(%f)",
			 cInd, outVector[cInd], cInd+dimension, outVector[cInd+dimension] );
 

	  
	} // colIndex loop 
     i_Nz=  i_Nz+elsInRow-1;
    } // rowIndex loop 
 for (int index = 0; index < dimension*2; index++)
    {
      glLogger.debug("OutVector [%d] = (%f)", index, outVector[index]);
    }  



}


/* Compute scalar product of two vectors 
   Not used at the moment: GoToDotC() is enough as it does glLogge.info("Eigvals+Eigvects:");
  for(int i=0;i<m;i++) 
    {
      glLogger.info("%f: ",eigval[i]);
      for(int j=0;j<n-1;j++) printf("%f, ",x[i*n+j]);
      glLogger.info("%f)",x[i*n+n-1]);
    } 
  
not need to 
   know anything about the matrix.
*/
double Matrix_sparse_cplx::dot(int *iflag,double z[],double w[])
{
  double sum=0;
  for(int i=0;i<dimension;i++) sum+=z[i]*w[i]; // Should work always?
  glLogger.info("Sparse matrix returns %f", sum);
  return sum;
}

/* Can be used to monitor progress of the diagonalization
   Also not used at the moment...
*/
void Matrix_sparse_cplx::monit(int *istate,int *nextit,int *nevals,int *nevecs,int *k,double f[],double d[])
{
  glLogger.error("Diagonalization progress for complex sparse matrix receives %d, %d, %d, %d, %d",
		*istate, *nextit, *nevals, *nevecs, *k);
  return;
}


/*!
  \fn int Matrix_sparse_cplx::diagonalize(std::complex<double> *eigVec, int novecs, std::complex<double> *eigval)
  \param eigVec  x       eigenvectors, stored one after another in the 1D array
           (size has to be at least (k+4)*dimension???)
  \param.novecs number of eiegnvectors in eigVec
   \param eigval  eigenvalues
\return some error code
\brief does the actual diagonalization (see NRC p. 456 ..)
*/
int Matrix_sparse_cplx::diagonalize(std::complex<double> *eigVec,
				    int novecs,
				    double *eigval)
{
  glLogger.error("Entering Matrix_sparse_cplx::diagonalize, try to find %d eigenvalues, dimension %d", novecs,dimension);
 /* std::vector<MKL_Complex16> entries;
  	std::vector<MKL_INT> columnVector;
  	std::vector<MKL_INT> rowVector;
*/
  //std::cout << "# Anzahl der non zero einträge " << mat_sp.size() << " Column vector size "<< col_ind.size() << " rows " << row_NrNz.size() << " dimension "<< dimension <<"\n";
  //dumpCSRMatrix(std::cerr);
  /*
 
  */
  bool success = false;
  int trialCount = 0;
  MKL_INT fpm[128]{ };
  ERROR ("Vor fesat init");
   fpm[0]=1;
 ::feastinit(fpm);
  ERROR("NACH FEAST INIT");
  fpm[0]=1;

  double eps = 0;
  MKL_INT loop = 0;
  // double emin = -3;
  //double emax = -1;
  MKL_INT m0 = m_m0*dimension; // was dimension before
  ERROR("Before init of eigenvectors");
  
  MKL_Complex16* eigenvectors = 0;
  ERROR("Before eigenvalues");
  double * eigenvalues=  new double[m0];
  ERROR("Before res");
    MKL_INT mode = 0;
      MKL_INT info = 0;
      char uplo ='U';
  // whil hier
  while (!success)
    {
      eigenvectors = new MKL_Complex16[m0 * dimension];
      eigenvalues=  new double[m0];
      vector<double> res(m0);
    
      MKL_INT l_dimension = dimension;
      ERROR("About to start zfeast");
      try{
      zfeast_hcsrev(&uplo, &l_dimension, reinterpret_cast<MKL_Complex16*>(&mat_sp[0]), &row_NrNz[0], &col_ind[0],
		fpm, &eps, &loop, &emin, &emax, &m0, eigenvalues,
		eigenvectors, &mode, &res[0], &info);
      }
      catch(const std::exception& ex)
	{
	  ERROR ("Exception from Zfeast");
	  glLogger.error(ex.what());
	}
      trialCount ++;
      std::cerr << "Info is " << info << "\n";
      ERROR("back from zfeast");
      if (info != 0) {
	ERROR("Some error happenend during zfeast.");
	glLogger.error("Info value is %d", info);
	if (info == 3)
	  {
	    
	    delete [] eigenvectors;
	    eigenvectors = 0;
	    delete [] eigenvalues;
	    eigenvalues = 0;
	    res.clear();
	    m0= (dimension-m0)*0.5 + m0;
	    glLogger.error("Enlarge m0 to %f", m0);
	  }
      }
      else {
	for (int evalIndex = 0; evalIndex < mode; evalIndex++)
	  {
	    *eigval = eigenvalues[evalIndex];
	    eigval++;
	    for (int eigenvecIndex = 0; eigenvecIndex < dimension; eigenvecIndex++)
	      {
		//*eigVec =  eigenvectors[eigenvecIndex];
	       
		eigVec++;
	      }
	  }
	success = true;
	delete [] eigenvectors;
	eigenvectors = 0;
	delete [] eigenvalues;
	eigenvalues = 0;
	res.clear();
      }
    
    } // while 
  // Copy eigenvalues back
  
  ERROR("Done with diagonalize");
  return info;
}




void Matrix_sparse_cplx::dumpMatrix(void) const
{
  return (dumpMatrix(std::cout));

}


void Matrix_sparse_cplx::dumpMatrix(std::ostream & outStream) const
{
  if (!outStream.good())
    {
      return;
    }
  outStream << "# Dumping upper right half sparse square matrix of dimension \n";
  outStream << "# first non comment line :dimension (integer)";
  outStream << "# rowIndex(int) colIndex(int) value(complex) \n";
  outStream << dimension<<"\n";
  int matPtr = 0;
  for (int rowIndex = 0; rowIndex < dimension; rowIndex ++)
    {
      //outStream << "row Nr ("<< rowIndex<<") has ("<< row_NrNz[rowIndex]<<") non zero elements \n";
      for (int colIndex = 0;
	   colIndex < row_NrNz[rowIndex];
	   colIndex++)
	{
	  outStream << rowIndex <<" "<<  col_ind[matPtr]<< "  "<< mat_sp[matPtr]<<"\n";
	  matPtr++;
	} // Loop over columns 
    }
  //! Loop over rows
  
}

int Matrix_sparse_cplx::diagLowestOnly(std::complex<double> *eigVec,
					int novecs,
					double *eigval)
{
  
  MKL_INT pm[128]{0,0,1};
  matrix_descr descrA;
  descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
  descrA.mode =SPARSE_FILL_MODE_UPPER;
  descrA.diag = SPARSE_DIAG_NON_UNIT;
  /*  
  mkl_sparse_d_emkl_sparse_d_ev ('S', pm, reinterpret_cast<MKL_Complex16*>(&mat_sp[0]),  descrA, MKL_INT k0, MKL_INT *k, double *E, double *X, double *res);
  */
  return (-1);

}
void Matrix_sparse_cplx::dumpCSRMatrix(std::ostream &outStream) const
{
	if (!outStream.good())
	    {
	      return;
	    }
	outStream << "# Dumping matrixelements (only nonzero parts)\n";
	outStream << "# Anzahl der non zero einträge " << mat_sp.size() << " Column vector size "<< col_ind.size() << " rows " << row_NrNz.size() << " dimension "<< dimension <<"\n";
	for (unsigned int matIndex = 0; matIndex <mat_sp.size();matIndex ++)
	{
		outStream << "Index " << matIndex << " Value " << mat_sp[matIndex] << "col Index = " << col_ind[matIndex]<< "\n";
	}
	outStream << "# Dumping row Vector (starting points of the rows ) \n";
	for (unsigned int rowIndex = 0; rowIndex < row_NrNz.size(); rowIndex ++ )
	{
		outStream << "Index " << rowIndex << " Value " << row_NrNz[rowIndex] << "\n";
	}
}
