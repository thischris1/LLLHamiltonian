/*! 
  \file Operator.cpp
  \brief Implementation of class Operator
  \author Karel Vyborny
 */

/*! \mainpage LLLham project

  Main idea: construct (1) a basis, (2) an arbitrary
  operator and (3) compute the operator's matrix elements in this
  basis and diagonalize the matrix. 

  (1) make your own basis by deriving a class from Basis (see
  e.g. yosBasis )

  (2) make your operator by deriving a class from Operator and
  implementing Operator::matEl_nonzero(int,int,std::complex<double>
 *), i.e. specifying its matrix element between two arbitrary
  elements of your basis (see e.g. LLLhamiltonian )

  (3) use Operator::diagonalize(...).

 */

#include <Operator.hpp>
#include <utils/logger.hpp>
#include <ERRORS.h>

#include <cassert>

#include <numerics/Matrix_sparse_cplx.hpp>
#include <numerics/RealSquareMatrix.h>

#include <boost/filesystem.hpp>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
//#include <armadillo>

#include <Eigen/PardisoSupport>
#include <Eigen/SparseLU>

#include <MatOp/SparseGenMatProd.h>

/*! \class Operator
  \brief This class can compute all matrix elements (matEls) of a given operator and
  find the eigvecs and eigvals of this matrix. The function which
  returns the matEl between i1-th and i2-th element of the basis
  (matEl_nonzero) has to be implemented by deriving this class (see
  e.g. class  LLLhamiltonian).  
 */


/*
  static members
 */
const double Operator::MIN_MATEL_SPARSE=1.e-7; 
const int Operator::FREQ_REPORT_MATEL=100;  
const int Operator::ImplicVal_ADD_TO_DIAG_WORKSPACE_SPARSE=140; 
const int Operator::ImplicVal_ADD_TO_DIAG_WORKSPACE_SPARSE_STEP=200; 
const long int Operator::MAX_SIZE_SPARSE= 25000000;
const long int Operator::MAX_SIZE_SPARSE_CPLX = 12500000;
const long int Operator::MAX_DIM_FOR_FULL=40000; 
const char Operator::FilenameFor_tmp_eigvecsApproxs[50]="/scratch/LLL_sparse_full_eigvecs.tmp";



/*! \fn Operator(Basis* new_basis,int new_type)
  Constructor.

  \param *new_basis     Pointer to the basis in which the matEls are
                        to be computed. Basis can be either directly 
			implemented (like yosBasis) or it can consist 
			of linear combinations of vectors of some other
			basis (see Basis).
  \param new_type       How to store the matrix elements (0=FULL_REAL,
                        1=SPARSE_REAL, 2=FULL_COMPLEX 3 SPARSE COMPLEX)
 */

// Expected type of the operator matrix





Operator::Operator(Basis* new_basis,int new_type):
		my_basis(new_basis),
				basisType(new_basis->getBasisType()),
				endBasis(new_basis),
				logfile(),
				inpf_firstGuess(NULL),
				outf_firstGuess(NULL),
				matrix_allocated(0),
				spectrum_allocated(0),
				m_writeMatrixOnly(false),
				row_NrNz(0),
				col_ind(0),
				mat_sp(0),
				nz_elements(0),
				x(0),
				mat_full_real(0),
				mat_full_cplx(0),
				x_cplx(0),
				ADD_TO_DIAG_WORKSPACE_SPARSE(ImplicVal_ADD_TO_DIAG_WORKSPACE_SPARSE),
				ADD_TO_DIAG_WORKSPACE_SPARSE_STEP(ImplicVal_ADD_TO_DIAG_WORKSPACE_SPARSE_STEP)
{

		ERROR ("CHECK FOR EISTNCE of WRITEMATRIX.dat");
	if (boost::filesystem::exists( "./WRITEMATRIX.dat" ) )
	  {
	    m_writeMatrixOnly = true;
	    ERROR("(ONLY) SAVING THE MATRIX");
	  }
	if (new_basis == 0)
	{
		std::cerr << "Operator::Operator "
				<< "receives a NULL Basis \n";
		return;
	}
	my_basis=new_basis;
	glLogger.info("Operator ctor with matType =(%d)", new_type);
	switch(new_type){
	case 5: matrix_type = SPARSE_COMPLEX_EIGEN;break;
	case 4: matrix_type = SPARSE_COMPLEX_EIGEN;break;
	case 3: matrix_type = SPARSE_COMPLEX; break;
	case 2: matrix_type=FULL_COMPLEX;break;
	case 1: matrix_type=SPARSE_REAL;break;
	case 0:break;
	default: matrix_type=FULL_REAL;
	}
	endBasis=my_basis->getEndBasis();  /* If the basis comprises of
      linear combinations, find what is the 'implemented' basis at the end. */
	matrix_allocated=0;
	// if File "WRITEMATRIX.dat"  exists, set the trigger to true

}

Operator::~Operator()
{
	if(matrix_allocated)
	{
		switch(matrix_type){
		case FULL_REAL:
			delete[] mat_full_real;
			break;
		case FULL_COMPLEX:
			delete [] mat_full_cplx;
			break;
		case SPARSE_REAL:
			delete [] mat_sp;
			break;
		case SPARSE_COMPLEX:
			delete []mat_sp;
			break;
		default:
			throw CxErrors(__FILE__, __LINE__);
		}
		matrix_allocated=0;
	}
}

//void Operator::setLogFile(FILE *new_logfile) {logfile=new_logfile;}

/*! \fn void setFileInp_FirstGuessOfEigvecs(FILE *inpf_firstGuess)
    \brief Diagonalization for sparse routines may be accelerated
    by giving a suitable initial guess for the eigenvectors. 
    By default this option is switched off at constructing Operator. 
    It may be switched on by calling this routine. The file is 
    expected to be of same type as files produced by 
    setFileOut_FirstGuessOfEigvecs.
 */
void Operator::setFileInp_FirstGuessOfEigvecs(FILE *new_inpf_firstGuess) 
{inpf_firstGuess=new_inpf_firstGuess;}

/*! \fn void setFileOut_FirstGuessOfEigvecs(FILE *outf_firstGuess)
    \brief If this method is called, the vectors which may be used as 
    initial guess next time will be written to file.
    See setFileInp_FirstGuessOfEigvecs.
 */
void Operator::setFileOut_FirstGuessOfEigvecs(FILE *new_outf_firstGuess)
{outf_firstGuess=new_outf_firstGuess;}

/*! \fn change_ADD_TO_DIAG_WORKSPACE_SPARSE(int new_ADD,int new_ADD_STEP)
  \brief Change parameters ADD_TO_DIAG_WORKSPACE_SPARSE and
  ADD_TO_DIAG_WORKSPACE_SPARSE_STEP. ADD_TO... is not changed 
  if new_ADD<1, the same for new_ADD_STEP.
 */
void Operator::change_ADD_TO_DIAG_WORKSPACE_SPARSE(int new_ADD,int new_ADD_STEP)
{
	if(new_ADD>0) ADD_TO_DIAG_WORKSPACE_SPARSE=new_ADD;
	if(new_ADD_STEP>0) ADD_TO_DIAG_WORKSPACE_SPARSE_STEP=new_ADD_STEP;
}




/*! \fn std::complex<double> Operator::matEl(int i_st1,int i_st2)
   Returns the value of matrix element (a); matEl_NonZero does the calculation
   and it tells moreover whether a is precisely zero or not.
 */
std::complex<double> Operator::matEl(int i_st1,int i_st2)
{
	std::complex<double> a(0.0,0.0);
	matEl_nonzero(i_st1,i_st2,&a);
	return a;
}

/*! \fn int Operator::matEl_nonzero(int i_st1,int i_st2,std::complex<double> *a)
  \param   i_st1,i_st2    The matrix element is to be computed between
  the i_st1-th and i_st2-th element of the basis.
  \param *a               Value of the computed matrix element.
  \returns 0 if the matrix element is exactly zero (e.g. due to
  selection rules), 1 if it is non-zero.

  In the basic concept, you should derive a class from Operator where
  you overwrite this method and define thus your own operator.

  The code which is implemented right here concerns the case when
  the basis elements are linear combinations of elements of some other
  basis (see Basis).
 */
int Operator::matEl_nonzero(int i_st1,int i_st2,std::complex<double> *a)
{
	if(matrix_allocated==0)
	{
		preprocessMatEls((int) matrix_type);
		if(matrix_allocated==0) {
			throw CxErrors(__FILE__, __LINE__);
		}
	}
	if(matrix_type!=FULL_REAL) {
		throw CxErrors(__FILE__, __LINE__);
	}
	*a=std::complex<double>(mat_full_real[i_st1*my_basis->dimension()+i_st2],0.0);
	return 1;
}


/*! \fn Operator::preprocessMatEls(int newBasis_matrix_type)
   Computes the matrix elements when the basis vectors consist 
   of linear combinations of vectors of some other basis. */
int Operator::preprocessMatEls(int newBasis_matrix_type)
{
	int dim,dim1;
	std::complex<double> xtmp;

	switch(my_basis->getBasisType()){
	case LINEAR_COMBINATION:
	{
		if(newBasis_matrix_type!=0) throw CxErrors(__FILE__, __LINE__);
		Operator op(my_basis->getMyBasis(),newBasis_matrix_type);
		op.preprocessMatEls(newBasis_matrix_type);
		dim=my_basis->dimension();
		dim1=my_basis->getMyBasis()->dimension();
		try
		{
			mat_full_real = new double [dim*dim];
		}
		catch(...)
		{
			throw (CxErrors("Matrix allocation failed at",
					__FILE__,
					__LINE__));
		}
		matrix_allocated=1;
		// Compute the matEl between i-th and j-th vector of the basis
		for(int i=0;i<dim;i++)
			for(int j=i;j<dim;j++)
			{
				xtmp=0;
				for(int i1=0;i1<dim1;i1++)
					for(int j1=0;j1<dim1;j1++)
						xtmp+=matEl(i1,j1)*my_basis->getCoef(i,i1)*my_basis->getCoef(j,j1);
				mat_full_real[i*dim+j]=xtmp.real();
				mat_full_real[j*dim+i]=xtmp.real();
			}
		return 1;
	}
	case YOS_BASIS:  // No preprocessing needed: just ask for matEls
		return 1;
	default:
		throw (CxErrors("preprocessing failed", __FILE__,__LINE__));
	}

	return 0;
}



/*! \fn int Operator::diagonalize(int eigsToFind,double *eigvec,double *eigval)
  \brief This is the main routine of Operator. It computes all matrix
  elements and finds the eigvals and eigvecs of the matrix.

  For the 

   \param eigsToFind   Nr. of eigvals & eigvecs to be found (may 
                       be not fully respected)
   \param eigvec       array (size at least dimension of matrix for real, 2*dimension of matrix for complex) to put computed eigenstates to. 
                       Complex eigenvectors will be saved by components, i.e. real, imag
   \param eigval       array (size at least dimension of matrix) to put computed eigenvalues to 

   \return always zero
 */
int Operator::diagonalize(int eigsToFind,std::vector<double>&eigvec,std::vector<double> &eigval)
{
	glLogger.error(" - Starting diagonalization Handling matrix by ");

	switch(matrix_type){
	case SPARSE_REAL: {
		glLogger.info("SPARSE_REAL.");
		diagonalizeSPARSEREAL(eigsToFind,&eigvec[0],&eigval[0]);
		break;}
	case FULL_REAL: {
		glLogger.info("FULL_REAL.");
		diagonalizeFULLREAL(eigsToFind,&eigvec[0],&eigval[0]);
		break;}
	case FULL_COMPLEX:
	{
		glLogger.error("FULL_COMPLEX.");
		diagonalizeFULLCOMPLEX(eigsToFind,&eigvec[0],&eigval[0]);
		break;
	}
	case SPARSE_COMPLEX:
	{
		glLogger.error("SPARSE COMPLEX Matrix type");
		diagonalizeSPARSECOMPLEX(eigsToFind,&eigvec[0], &eigval[0]);
		break;
	}
	case SPARSE_COMPLEX_EIGEN:
	{
		ERROR("SPARSE EIGEN TYPE USED");
		std::vector<double> eigValVec;
		std::vector< std::vector < std::complex <double> > > eigvecVec;

		diagonalizeSparseEigen(eigsToFind, eigvec, eigval);

		break;
	}
	default:
	{
		throw (CxErrors("Bad Matrix type"));
	}
	}
	glLogger.error("Operator Done with diagonalize");
	return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

int Operator::diagonalizeFULLREAL(int eigsToFind,double *eigvec,double *eigval)
{
	int itmp;
	long int iltmp;
	char stmp[10],stmp1[10];
	int dim;
	time_t t1,t2;double tdiff;
	// Allocate the arrays to store matrix elements in
	dim=my_basis->dimension();
	if(dim>MAX_DIM_FOR_FULL) {
		throw CxErrors("Matrix too large", __FILE__, __LINE__);
	}
	iltmp=dim*dim;
	//mat_full_real=new double[iltmp];
	//matrix_allocated=1;
	iltmp=dim*eigsToFind;
	x=new double [iltmp];

	RealSquareMatrix m_Matrix(dim);

	// Compute the matrix elements
	CMDCOMPLEX b; // temporary variable for keeping matrix element
	glLogger.info ("Computing matrix elements...");

	t1=time(0);
	for(int i1=0;i1<dim;i1++)
	{
		for(int i2=i1;i2<dim;i2++)
		{
#
			itmp=matEl_nonzero(i1,i2,&b);
			/*
			 *(mat_full_real+i1*dim+i2)=b.real();
			 * */
			m_Matrix.setElement(b, i1,i2);

			glLogger.info ("Final MatEl %d, %d: %f.\n",i1,i2,b.real());

		}
		if(!(i1%FREQ_REPORT_MATEL))
		{
			if(!i1) Npercent(-1.,stdout);
			Npercent(1.*i1/dim,stdout);
		}
	}
	t2=time(0);tdiff=difftime(t2,t1);
	Ntime(tdiff,stmp1);
	/*
  iltmp=dim*dim*sizeof(*mat_full_real);
  Nbytes(iltmp,stmp);
	 */
	glLogger.info ("...done (%s used, took %s).\n",stmp,stmp1);

	// Diagonalize it
	glLogger.info ("Entering diagonalization...");
	t1=time(0);
	//  diag_full_real_matrix(dim,eigsToFind,mat_full_real,eigval,x);
	DOUBLE_VECTOR eVals(dim);
	COMPLEX_VECTOR eVecs(dim*dim);
	m_Matrix.getSpectrum(eVals,eVecs);
	for (int eigValIndex = 0; eigValIndex < eigsToFind; ++eigValIndex) {
		eigval[eigValIndex] = eVals.at(eigValIndex);
		for (int vectorIndex = 0; vectorIndex < dim; ++vectorIndex) {
			*(eigvec+eigValIndex*dim+vectorIndex)=eVecs[eigValIndex*dim+vectorIndex].real();
		}
	}


	t2=time(0);tdiff=difftime(t2,t1);
	Ntime(tdiff,stmp1);
	glLogger.info ("...done, took %s.\n",stmp1);



	// deallocate arrays which keep matEls
	delete[] mat_full_real;
	matrix_allocated=0;
	// translate x into eigvec
	/*
  for(int i=0;i<eigsToFind;i++)
    for(int j=0;j<dim;j++)
	 *(eigvec+i*dim+j)=x[i*dim+j];
  // deallocate array which keeps the comptd eigvecs
 // delete[] x;
	 */
	return 0;
}
#pragma GCC diagnostic pop
int Operator::diagonalizeSPARSEREAL(int eigsToFind,double *eigvec,double *eigval)
{
	int itmp,novecs;
	long int iltmp;
	char stmp[10],stmp1[10];
	int dim,add_to_diag_workspace;
	time_t t1,t2;double tdiff;

	FILE *f_tmp_eigvecsApproxs;
	FILE *f_tmp=NULL; // aux. variable

	// Allocate the arrays to store matrix elements in
	dim=my_basis->dimension();
	try
	{
		row_NrNz = new int [dim+MAX_SIZE_SPARSE];
		col_ind=row_NrNz+dim;
		mat_sp= new double [MAX_SIZE_SPARSE];
		nz_elements=0;
		matrix_allocated=1;
		x= new double [dim*(eigsToFind+ADD_TO_DIAG_WORKSPACE_SPARSE)];
	}

	catch(...)
	{
		throw CxErrors("Could not allocate memory at",
				__FILE__,
				__LINE__);
	}
	// Compute the matrix elements
	std::complex<double> b; // temporary variable for keeping matrix element
	glLogger.info("Computing matrix elements...");

	t1=time(0);
	for(int i1=0;i1<dim;i1++)
	{
		row_NrNz[i1]=0;
		for(int i2=i1;i2<dim;i2++)
		{
			if(!matEl_nonzero(i1,i2,&b)) continue;
			if(fabs(b.real())<MIN_MATEL_SPARSE) continue;
			// Only non-zero matrix elements are strored
			if(nz_elements>=MAX_SIZE_SPARSE)
			{
				throw CxErrors("Too many elements for sparse matrix in diagonalizeReal",
						__FILE__,
						__LINE__);
			}
			row_NrNz[i1]++;
			col_ind[nz_elements]=i2;
			mat_sp[nz_elements]=b.real();
			nz_elements++;
			glLogger.info("The Hamilton Matrix is :");
			glLogger.info("%d, %d: %20.15f.\n",i1,i2,b.real());
		}
		if(!(i1%FREQ_REPORT_MATEL))
		{
			if(!i1) Npercent(-1.,stdout);
			Npercent(1.*i1/dim,stdout);
		}	
	}
	Npercent(1.1,stdout);
	t2=time(0);
	tdiff=difftime(t2,t1);
	Ntime(tdiff,stmp1);
	iltmp=dim*sizeof(*row_NrNz)+nz_elements*(sizeof(*col_ind)+sizeof(*mat_sp));
	Nbytes(iltmp,stmp);
	glLogger.info ("...done: %ld non-zero elements (sparsity %5g, %s used), took %s.\n",
			nz_elements,1.*nz_elements/dim/dim,stmp,stmp1);


	// create a matrix object and try to diagonalize it (itmp is the exit err code)
	itmp=-1;
	add_to_diag_workspace=ADD_TO_DIAG_WORKSPACE_SPARSE;
	glLogger.info("Entering diagonalization...");

	novecs=0;
	if(inpf_firstGuess)
	{
		novecs=getFirstGuessFromFile(dim,x)-1;
	}
	/* try to get initial approximations for the eigvecs;
     novecs = Nr. of vectors passed as initial guess to mat */
	Matrix_sparse_real mat(row_NrNz,col_ind,mat_sp,nz_elements,
			dim,eigsToFind,add_to_diag_workspace);
	while(itmp)  // this loop attempts to diagonalize repeatedly (if the first try fails)
	{
		t1=time(0);
		itmp=mat.diagonalize(x,novecs,eigval);
		t2=time(0);
		tdiff=difftime(t2,t1);
		Ntime(tdiff,stmp1);
		if(itmp) // if diagonalization failed, try to increase add_to_diag_workspace
		{
			glLogger.error("...failed with exit code %d after %s. ",itmp,stmp1);

			if(itmp==2 || itmp==3 || itmp==4)
			{
				glLogger.info("Trying to increase k from %d to %d.\n",
						add_to_diag_workspace,
						add_to_diag_workspace+ADD_TO_DIAG_WORKSPACE_SPARSE_STEP);

				if((f_tmp_eigvecsApproxs=
						fopen(FilenameFor_tmp_eigvecsApproxs,"w"))!=NULL)
				{ // write down the approxs. for eigvecs ...
					f_tmp=outf_firstGuess;
					outf_firstGuess=f_tmp_eigvecsApproxs;
					writeFirstGuessToFile(eigsToFind+add_to_diag_workspace,dim,x);
				}
				add_to_diag_workspace+=ADD_TO_DIAG_WORKSPACE_SPARSE_STEP;
				mat.set_k(add_to_diag_workspace);
				delete[] x;
				x = new double [dim*(eigsToFind+add_to_diag_workspace)];
				if(f_tmp_eigvecsApproxs!=NULL)
				{ // ... and read them again
					fclose(f_tmp_eigvecsApproxs);
					f_tmp_eigvecsApproxs=fopen(FilenameFor_tmp_eigvecsApproxs,"r");
					inpf_firstGuess=f_tmp_eigvecsApproxs;
					novecs=getFirstGuessFromFile(dim,x)-1;
					inpf_firstGuess=NULL;
					outf_firstGuess=f_tmp;
					fclose(f_tmp_eigvecsApproxs);
				}
				if(add_to_diag_workspace<dim)
				{
					glLogger.info("Entering diagonalization...");
					continue;
				}
			}
			throw CxErrors(__FILE__, __LINE__);
			throw CxErrors(__FILE__, __LINE__);
		}
	} // diagonalization successful
	glLogger.info("...done, took %s.\n",stmp1);

	// deallocate arrays which keep matEls
	delete[] row_NrNz;
	delete[] mat_sp;
	mat.~Matrix_sparse_real();
	matrix_allocated=0;
	// translate x into eigvec
	for(int i=0;i<eigsToFind;i++)
		for(int j=0;j<dim;j++)
			*(eigvec+i*dim+j)=x[i*dim+j];
	// write eigvecs into a file so as to use them later as initial guess
	if(outf_firstGuess!=NULL) writeFirstGuessToFile(eigsToFind+add_to_diag_workspace,dim,x);
	// deallocate array which keeps the comptd eigvecs
	delete[] x;

	return 0;
}


/*!
  \fn int Operator::diagonalizeSPARSECOMPLEX(int eigsToFind,double *eigvec, double *eigval)
  \brief diagonalize a sparse complex matrix 
  \param eigsToFind number of requested eigenvalues
  \param  eigvec array with coefficents of eigenvectors
  \param  eigval array with eigenvalues
 */
int Operator::diagonalizeSPARSECOMPLEX(int eigsToFind,
		double *eigvec,
		double *eigval)
{
	int
	itmp = 0;
	ofstream matFile;
	long  int iltmp = 0;
	char
	stmp[10],
	stmp1[10];

	int
	dim,
	add_to_diag_workspace;

	time_t
	t1,
	t2;

	double tdiff;
	dim=my_basis->dimension();
	glLogger.error("Starting with dimension (%d)", dim);
	if (m_writeMatrixOnly) {
	  // open file
	  ERROR ("Just Writing matrix");
	  matFile.open("matrix_sparse.dat");
	  matFile << "# dimension " << dim <<std::endl;
	}


	FILE *f_tmp_eigvecsApproxs = 0;
	FILE *f_tmp = 0; // aux. variable

	// Allocate the arrays to store matrix elements in

	int * colIndex = 0;
	typedef Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor >  SpMat;
	SpMat theMatrix(dim,dim);
	// Compute the matrix elements
	std::complex<double> locMatEl; // temporary variable for keeping matrix element
	glLogger.info("Computing matrix elements...");

	t1=time(0);
	int elementCount =1;
	std::vector< std::complex<double> > cplx_matEl;
	std::vector<MKL_INT> columnVector;
	std::vector<MKL_INT> rowVector;

	for(int rowIndex=0;rowIndex<dim;rowIndex++)
	{
		rowVector.push_back(elementCount);
		for(int columnIndex=rowIndex;columnIndex<dim;columnIndex++)
		{

			if(!matEl_nonzero(rowIndex,columnIndex,&locMatEl))
			{
				continue;
			}
			if(fabs(locMatEl.real())<MIN_MATEL_SPARSE)
			{
				continue;
			}
			// Only non-zero matrix elements are strored
			if (m_writeMatrixOnly) {
			  matFile << rowIndex << "  " << columnIndex << "  " << locMatEl << std::endl;
			  continue;
			} 

			cplx_matEl.push_back(locMatEl);
			columnVector.push_back(columnIndex+1);
			elementCount++;

			glLogger.info("The Hamilton Matrix is :");
			glLogger.info("%d, %d: %20.15f.\n",rowIndex,columnIndex,locMatEl.real());

		}

	}


	std::cerr<<elementCount<<std::endl;
	rowVector.push_back(cplx_matEl.size()+1);
	Npercent(1.1,stdout);
	t2=time(0);tdiff=difftime(t2,t1);
	Ntime(tdiff,stmp1);
	iltmp=dim*sizeof(*row_NrNz)+nz_elements*(sizeof(*colIndex)+sizeof(*mat_sp));
	Nbytes(iltmp,stmp);
	glLogger.error ("...done: %ld non-zero elements (sparsity %5g, %s used), took %s.\n",
			elementCount, elementCount/dim/dim,stmp,stmp1);

	if (m_writeMatrixOnly ) {
	  matFile.close();
	  ERROR("DONE WRITING, LEAVE NOW");
	  return (-1);
	}
	// create a matrix object and try to diagonalize it (itmp is the exit err code)
	itmp=-1;
	add_to_diag_workspace=ADD_TO_DIAG_WORKSPACE_SPARSE;
	glLogger.info("Entering diagonalization...");

	int novecs = 0;
	if(inpf_firstGuess!=NULL) {
		novecs=getFirstGuessFromFile(dim,x)-1;
	}
	/* try to get initial approximations for the eigvecs;
     novecs = Nr. of vectors passed as initial guess to mat */
	Matrix_sparse_cplx mat(rowVector,columnVector,
			cplx_matEl, nz_elements,
			dim,eigsToFind,
			add_to_diag_workspace);
	std::vector<MKL_Complex16>cEigVec;
	cEigVec.reserve(dim*eigsToFind);
	while(itmp)  // this loop attempts to diagonalize repeatedly (if the first try fails)
	{
		t1=time(0);
		//      double *eigVal = new double[eigsToFind];

		glLogger.info("Size of cEigVec novecs=(%d),dim= (%d), eigsToFInd=(%d), (%d)",
				novecs, dim, eigsToFind,
				dim*eigsToFind);
		itmp=mat.diagonalize(cEigVec,eigsToFind,eigval);
		t2=time(0);
		tdiff=difftime(t2,t1);
		Ntime(tdiff,stmp1);
		if(itmp) // if diagonalization failed, try to increase add_to_diag_workspace
		{
			glLogger.error("...failed with exit code %d after %s. ",itmp,stmp1);

			if(itmp==2 || itmp==3 || itmp==4)
			{
				glLogger.info("Trying to increase k from %d to %d.\n",
						add_to_diag_workspace,
						add_to_diag_workspace+ADD_TO_DIAG_WORKSPACE_SPARSE_STEP);

				if((f_tmp_eigvecsApproxs=
						fopen(FilenameFor_tmp_eigvecsApproxs,"w"))!=NULL)
				{ // write down the approxs. for eigvecs ...
					f_tmp=outf_firstGuess;
					outf_firstGuess=f_tmp_eigvecsApproxs;
					writeFirstGuessToFile(eigsToFind+add_to_diag_workspace,dim,x);
				}
				add_to_diag_workspace+=ADD_TO_DIAG_WORKSPACE_SPARSE_STEP;
				mat.set_k(add_to_diag_workspace);
				delete[] x;
				x = new double [dim*(eigsToFind+add_to_diag_workspace)];
				if(f_tmp_eigvecsApproxs!=NULL)
				{ // ... and read them again
					fclose(f_tmp_eigvecsApproxs);
					f_tmp_eigvecsApproxs=fopen(FilenameFor_tmp_eigvecsApproxs,"r");
					inpf_firstGuess=f_tmp_eigvecsApproxs;
					novecs=getFirstGuessFromFile(dim,x)-1;
					inpf_firstGuess=NULL;
					outf_firstGuess=f_tmp;
					fclose(f_tmp_eigvecsApproxs);
				}
				if(add_to_diag_workspace<dim)
				{
					glLogger.info("Entering diagonalization...");
					continue;
				}
			}
			throw CxErrors(__FILE__, __LINE__);
		}
	} // diagonalization successful
	glLogger.info("...done, took %s.\n",stmp1);

	// deallocate arrays which keep matEls
	cplx_matEl.clear();
	rowVector.clear();
	columnVector.clear();

	matrix_allocated=0;
	// translate x into eigvec
	/*
	 */
	//! \todo fix this (should not depend on memeory layout
	//memcpy (eigvec, cEigVec, sizeof(double)*2*dim*eigsToFind);
	// write eigvecs into a file so as to use them later as initial guess
	if(outf_firstGuess!=NULL)
	{
		writeFirstGuessToFile(eigsToFind+add_to_diag_workspace,dim,x);
	}
	// deallocate array which keeps the comptd eigvecs
	delete[] x;
	cEigVec.clear();

	return 0;
}


void Operator::diagonalizeSparseEigen(int eigsToFind,std::vector < double >   & eigvec,std::vector <double> &eigval)
{
	long l_dim =my_basis->dimension();
	Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor > SpMat ( l_dim,l_dim);
	Eigen::SparseMatrix<std::complex<double> > SpMatsa;

	Eigen::SparseMatrix<std::complex<double> > B (l_dim,l_dim);
//	arma::sp_cx_mat SpMat(l_dim,l_dim);
	std::complex<double> locMatEl; // temporary variable for keeping matrix element
	glLogger.info("Computing matrix elements...");



		for(long rowIndex=0;rowIndex<l_dim;rowIndex++)
		{

			for(long columnIndex=rowIndex;columnIndex<l_dim;columnIndex++)
			{
				if(!matEl_nonzero(rowIndex,columnIndex,&locMatEl))
				{
					continue;
				}
				if(fabs(locMatEl.real())<MIN_MATEL_SPARSE)
				{
					continue;
				}
				// Only non-zero matrix elements are strored
				SpMat.insert(rowIndex, columnIndex) = locMatEl;
				//SpMat(columnIndex,rowIndex)= conj(locMatEl);


				glLogger.info("The Hamilton Matrix is :");
				glLogger.info("%d, %d: %20.15f.\n",rowIndex,columnIndex,locMatEl.real());

			}

		}

		// Matrix is filled
		if (!SpMat.isCompressed()) {
					SpMat.makeCompressed();
					//B.makeCompressed();
				}
		// make it selfadjoined
		SpMatsa = SpMat.selfadjointView<Eigen::Upper>();

		Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<std::complex<double> > > eigensolvers(SpMatsa);
		Eigen::VectorXcd tempEigVal = eigensolvers.eigenvalues();
	    Eigen::MatrixXcd Eigen_eigvecs = eigensolvers.eigenvectors();
		//eigval.resize(tempEigVal.size());
		unsigned int max = eigsToFind;
		std::vector<double>::iterator it = eigvec.begin();
		if (tempEigVal.size() < eigsToFind) {
			max = tempEigVal.size();
		}

		//Eigen::VectorXcd::Map(eigval[0], tempEigval.size()) = tempEigval;
		for (unsigned int index = 0; index < max;index++)
		{

			eigval.push_back(tempEigVal[index].real());
			unsigned int start = index*l_dim;
			Eigen::VectorXcd aState = Eigen_eigvecs.col(index);
			for (unsigned int dimIndex = 0; dimIndex<l_dim;dimIndex++)
			{
				it = eigvec.insert(it,aState(dimIndex).real());
				it = eigvec.insert(it,aState(dimIndex).imag());

			}
		}
		//eigvec = eigensolvers.eigenvectors();
/*

		SpMat.cols();
		SpMat.rows();
		SpMat.innerNonZeroPtr();
		//Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> > solver;
		Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> >, Eigen::COLAMDOrdering<int> > solver_lu;
		solver_lu.analyzePattern(SpMat);
		solver_lu.factorize(SpMat);
		solver_lu.compute(SpMat);
		std::cerr << solver_lu.info();
		Eigen::SparseQR<Eigen::SparseMatrix<std::complex<double> >, Eigen::COLAMDOrdering<int> > solver_qr;

		solver_qr.analyzePattern(SpMat);
		solver_qr.factorize(SpMat);
		solver_qr.compute(SpMat);
		std::cerr << solver_qr.info();
	*/
		//SpMat.
	/*
		Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<std::complex<double> > > es(SpMat,B,Eigen::EigenvaluesOnly);
		es.compute(SpMat, B);
		*/
}



/*! \fn void writeFirstGuessToFile(int eigsToFind,double *x)
 */
void Operator::writeFirstGuessToFile(int eigsToFind,int dim,double *x)
{  
	/* Format of the file: 1) sizeof(dim), dim
                         2) sizeof(eigsToFind), eigsToFind
			 3) data (i.e. eigsToFind*dim doubles 
			          = eigsToFind vectors)
	 */
	int tmp_size;
	glLogger.info ("Writing vectors for initial guess to file.\n");

	fseek(outf_firstGuess,0,SEEK_SET);  // write to the beginning of the file
	tmp_size=sizeof(dim);
	fwrite(&tmp_size,sizeof(int),1,outf_firstGuess);
	fwrite(&dim,sizeof(dim),1,outf_firstGuess);
	tmp_size=sizeof(eigsToFind);
	fwrite(&tmp_size,sizeof(int),1,outf_firstGuess);
	fwrite(&eigsToFind,sizeof(eigsToFind),1,outf_firstGuess);
	fwrite(x,sizeof(*x),dim*eigsToFind,outf_firstGuess);
}


/*! \fn int getFirstGuessFromFile(double *x)
  \brief Assuming format as from writeFirstGuessToFile
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
                    /* no diagnostic for this one */

int Operator::getFirstGuessFromFile(int dim,double *x)
{  
	int tmp_size,tmp_dim;
	int eigsToGet,eigsGot;
	fseek(inpf_firstGuess,0,SEEK_SET);  // write to the beginning of the file
	// read the dimension of vectors in the file
	fread(&tmp_size,sizeof(int),1,inpf_firstGuess);
	if(tmp_size!=sizeof(dim)) throw CxErrors(__FILE__, __LINE__);
	fread(&tmp_dim,tmp_size,1,inpf_firstGuess);
	if(tmp_dim!=dim) throw CxErrors(__FILE__, __LINE__);
	// read the Nr. of vectors in the file
	fread(&tmp_size,sizeof(int),1,inpf_firstGuess);
	if(tmp_size!=sizeof(eigsToGet)) throw CxErrors(__FILE__, __LINE__);
	fread(&eigsToGet,tmp_size,1,inpf_firstGuess);
	glLogger.info("Getting %d vectors as initial guess, ", eigsToGet);
	eigsGot=fread(x,sizeof(*x),dim*eigsToGet,inpf_firstGuess);
	glLogger.info("got %d vectors, ",eigsGot/dim);
	return eigsGot/dim;
}
#pragma GCC diagnostic pop

/*! \fn int Operator::diagonalizeFULLCOMPLEX(int eigsToFind,double *eigvec,double *eigval)
  \todo Proper returning of the eigvecs.
  \param eigvec [on exit only] Eigenvectors 
  \param eigval array of eigenvalues (same ordering as in eigvec)
  \param eigsToFind number of requested eigenvalues/vectors
  \brief mehtod cretaes a matrix_full_cpx object and uses the diagonalize method on it.
 */
int Operator::diagonalizeFULLCOMPLEX(int eigsToFind,double *eigvec,double *eigval)
{
	long int iltmp = 0;
	char stmp[10],stmp1[10];
	int dim = 0;
	time_t
	t1 = 0,
	t2 =0;
	double tdiff = 0.0;



	// Allocate the arrays to store matrix elements in
	dim=my_basis->dimension();
	if(dim>MAX_DIM_FOR_FULL) {
		throw CxErrors("Size of basis in FULL COMPLEX too large for being numbered by int/long int or whatever you have set",
				__FILE__, __LINE__);
		// see /usr/include/limits.h
	}
	int arraySize = dim*dim;
	try{
		mat_full_cplx= new dComplex [arraySize];
		matrix_allocated=1;
		x_cplx= new dComplex [dim*eigsToFind]; // array for eigenvectors
	}
	catch (...)
	{
		throw CxErrors("Could not allocate enough memory at",
				__FILE__,
				__LINE__);
	}

	// Compute the matrix elements
	std::complex<double> b; // temporary variable for keeping matrix element
	std::ofstream matrixFile("matrix.dat");

	glLogger.info ("Computing matrix elements in diagonalizeFULLCOMPLEX");

	t1=time(0);
	for(int i1=0; i1<dim; i1++)
	{
		for(int i2=0;i2<=i1;i2++)
		{
			b=matEl(i1,i2);
			mat_full_cplx[i1+i2*dim] = b;

			glLogger.debug ("%d, %d: %20.15f + i*%20.15f.\n", i1, i2, b.real(), b.imag());
			if ( m_writeMatrixOnly == true)
			{
				matrixFile << i1<< i2<< b.real()<< b.imag();

			}

		}
		if(!(i1%FREQ_REPORT_MATEL))
		{
			if(!i1) Npercent(-1.,stdout);
			Npercent(1.*i1/dim,stdout);
		}
	}
	t2=time(0);
	tdiff=difftime(t2,t1);
	Ntime(tdiff,stmp1);
	iltmp=dim*dim*sizeof(*mat_full_cplx);
	Nbytes(iltmp,stmp);
	glLogger.info ("...done (%s used, took %s).\n",stmp,stmp1);

	if (m_writeMatrixOnly == true)
	{
		matrixFile.close();
		ERROR("MATRIX WRITTEN, EXITING NOW");
		exit(0);
	}
	// Diagonalize it
	glLogger.info ("Entering diagonalization in FULL COMPLEX MODE...");

	t1=time(0);
	int retVal =  diag_full_cplx_matrix(dim,eigsToFind,mat_full_cplx,eigval,x_cplx);
	glLogger.info("Diagonalization returns (%d)",retVal);
	t2=time(0);
	tdiff=difftime(t2,t1);

	Ntime(tdiff,stmp1);
	//glLogger.error ("done with diagonalization, took %s.\n", stmp1);


	//test, if i-th eigenvector is eigenvector to eigenvalue i
	/*#ifdef DEBUG
  std::complex<double> * testV1 = new std::complex<double>[dim];
  std::complex<double> * testV2 = new std::complex<double>[dim];

  double absDiff = 0.0;
  for (int n=0; n < eigsToFind; n= n+1)
    {
		glLogger.info("Analyse eigenvector no %d, dim = (%d), energy = %f", n,dim,eigval[n] );
      for (int i=0; i<dim; i++)
		{
	  	testV1 [i] = std::complex<double>(0.0,0.0);
	  	for (int j=0; j<dim; j++) 
	    	{
	  		if (glLogger.getLogLevel() == INFO)
	  		{
	  			if ((n == 0 ) && (i==0))
	  			{
	  				// print once
	  				std::complex<double> testOnly = ((std::complex<double>*)x_cplx)[dim*n+j];

	  				glLogger.info("eigenvector [%d],[%d] = (%f, %f)",n,j,testOnly.real(),testOnly.imag());
	  			}
	  		}
	      		testV1[i] += matEl(i,j) * ((std::complex<double>*)x_cplx)[dim*n+j];
	    	}
	  	testV2[i] = ((std::complex<double>*)x_cplx)[dim*n+i] * eigval[n];

	  	absDiff += abs(testV1[i]-testV2[i]); 

		}
		glLogger.info("absdiff = %e", absDiff);
    }

  delete[] testV1;
  delete[] testV2;

  glLogger.error ("Test of eigenvectors, difference^2 = %f (should be nearly 0.0)", absDiff);
  if (absDiff > 1e-06)
    { 
    	ERROR("Difference of eigenvectors too large, leaving!");
      throw (CxErrors("Difference too large"));
    }

#endif
	 */
	// deallocate arrays which keep matEls
	delete[] mat_full_cplx;

	// translate x into eigvec

	//for(int i=0;i<eigsToFind;i++)
	//coefToState(my_basis,x+i*my_basis->size(),eigvec[i])
	/* // DREHEN DER EIGENVEKTOR_KOMPONENTEN

  dComplex *pointer = x_cplx;

  for (int i = 0; i<eigsToFind; i++)
  {
	  for (int j = 0; j <dim;j++)
	  {
		  pointer->real = pointer->real*-1;
		  pointer->imag = pointer->imag*-1;
		  pointer++;
	  }
  }
  // ENDE DREHEN --> FÜR NAG NICHT NÖTIG*/

	memcpy (eigvec, x_cplx, sizeof(double)*2*dim*eigsToFind);

	// deallocate array which keeps the comptd eigvecs

	delete[] x_cplx;
	matrix_allocated=0;
	return 0;
}


bool Operator::writeMatrixToFile(std::string& fileName, int maxDim)
{
	int upperLimit = 0;
	if (maxDim == -1)
	{
		upperLimit = my_basis->dimension();
	}
	else
	{
		upperLimit = maxDim;
	}
	std::complex<double> locMatEl;
	std::ofstream outFile(fileName.c_str());
	outFile << "# dimension of Matrix \n";
	outFile << upperLimit << "\n";
	for (int xIndex = 0; xIndex < upperLimit; ++xIndex) {
		for (int yIndex = 0; yIndex < upperLimit; ++yIndex) {
			outFile << xIndex <<"  "<< yIndex;

			if(!matEl_nonzero(xIndex,yIndex,&locMatEl))
			{
				outFile << "  " << 0.0 <<"\n";
			}
			else {
				outFile <<"   " << locMatEl.real() <<"  "<<locMatEl.imag() <<"\n";

			}
		}

	}

	return true;
}










