/*!
  \file yosEigenState.cpp
  \author Christian Mueller
  \brief Implementation of the class yosEigenState 
 */
#include <LLLlib/yosEigenState.hpp>
#include <vector>
#include <fstream>
#include <complex>



#include <complex>

#include <LLLlib/LLLlib.h>


#include <numerics/matrix_full_cplx.hpp>
#include <math.h>
#include <utils/CRandomizerGsl.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multimin.h>
#include <utils/CRandomizerGsl.h>

#include <utils/CRandomizerMkl.h>


using std::complex;

#define TINY 1e-07

void objfun(int *, int *, int *, int *, double *, double *, double *, int *, int *, double *);

void (*pObjfun) (int *, int *, int *, int *, double *, double *, double *, int *, int *, double *);
// For Intel minimizer
void jacobiObjFun(MKL_INT  *, MKL_INT *, double *, double *);
void (*pjacobiObjFun)(MKL_INT *, MKL_INT *, double *, double *);


yosEigenState * yosEigenState::m_instance = 0;   



void my_f (MKL_INT *m, MKL_INT *n, double *x, double *fvec)
{
	//ignore m and n

	glLogger.info("Called MY FUnction");

	std::complex <double> fctVal = yosEigenState::getInstance()->sumofSPWavefunctions(0,x[0],x[1]);
	fvec[0]=fctVal.real();
	fvec[1]=fctVal.imag();



}

double my_f(const gsl_vector *v, void *params)
{
	double x[2];

	x[0] = gsl_vector_get(v,0);
	x[1] = gsl_vector_get(v,1);
	if ( (x[0]<0.0) || (x[0]>1.0) || (x[1] < 0.0) || (x[1] > 1.0)  )
	{

		return (1.0e300);

	}
	else
	{

		std::complex <double> fctVal = yosEigenState::getInstance()->sumofSPWavefunctions(0,x[0],x[1]);
		//	   return (fctVal.real());


		gsl_complex z ;

		GSL_SET_COMPLEX(&z, fctVal.real(), fctVal.imag());

		double retVal = gsl_complex_abs (z);
		glLogger.info("Called MY FUnction at %f, %f, return %f ", x[0],x[1], retVal);
		return (retVal);
	}
}
void my_df(double *x, void *params, double * fjac)
{

	gsl_vector *in = gsl_vector_alloc(2);
	gsl_vector *result = gsl_vector_alloc(4);

	gsl_vector_set(in,0,x[0]);
	gsl_vector_set(in,1,x[1]);
	my_df_gsl(in,params,result);
	// copy back
	fjac[0] = gsl_vector_get(result,0);
	fjac[1] = gsl_vector_get(result,1);
	fjac[2] = gsl_vector_get(result,2);
	 fjac[3] = gsl_vector_get(result,3);
	 gsl_vector_free(result);
	 gsl_vector_free(in);


}
/*! calculate the jacobian matrix (2x2) with [
 * [0]= real df/dx
 * [1] = imag df/dx
 * [2] = real df/dy
 * [3]= imag df/dy
 * \param x : Position
 * \param params Any parameters needed (0 mostly)
 * \param g result vector of size 4
 */
void my_df_gsl (const gsl_vector * x, void * params, gsl_vector * g)
{
	INFO("called my_df");
	// check if position is inside rectangle, if so return function value, if not return huge value

	if ( (gsl_vector_get(x,0)<0.0) || (gsl_vector_get(x,0)>1.0) || (gsl_vector_get(x,1) < 0.0) || (gsl_vector_get(x,1) > 1.0)  )
	{
		if ((gsl_vector_get(x,0)<0.0 || (gsl_vector_get(x,1)<0.0)))
		{
			gsl_vector_set(g,0,-1000);
			gsl_vector_set(g,0,-1000);

		}
		if ((gsl_vector_get(x,0)>1.0 || (gsl_vector_get(x,1)>0.0)))
				{
					gsl_vector_set(g,0,1000);
					gsl_vector_set(g,0,1000);

				}
		return;


	}
	else
	{
		CxPeriodicPosition position(gsl_vector_get(x,0),gsl_vector_get(x,1));

		complex <double> fctdx = yosEigenState::getInstance()->sumofdxSPWavefunctions(0,position);
		complex <double> fctdy =yosEigenState::getInstance()->sumofdySPWavefunctions(0,position);
		// write back to vector g
		gsl_complex z ;

		GSL_SET_COMPLEX(&z, fctdx.real(), fctdx.imag());



		gsl_vector_set(g,0,fctdx.real());
		gsl_vector_set(g,1,fctdx.imag());
		gsl_vector_set(g,2,fctdy.real());
		gsl_vector_set(g,3,fctdy.imag());



	}

}

void my_fdf(const gsl_vector * x, void * params, double * f, gsl_vector * g)
{
	double fctVal = my_f(x,params);
	f = &fctVal;
	INFO ("called myfdf");
	my_df_gsl(x,params,g);
}

/*!
  \fn yosEigenState::yosEigenState(yosBasis *newBasis, complex <double> * coefficents) : yosState()
  \param newBasis The Basis of yoshioka states
  \param coefficents the coefficents (the so called eigenvector form diagonalization routines)
  \brief constructor
 */
yosEigenState::yosEigenState(yosBasis *newBasis, std::complex <double> * coefficents) : 
		  persist_eigSt((Basis *) newBasis, true),
		  fixedPart(0),
		  isSearchInitialized(false),
		  fvec(0),
		  fjac(0)
{
	if ((newBasis == 0) || (coefficents == 0))
	{
		throw (CxErrors(__FILE__, __LINE__));
	}

	init (newBasis, coefficents);
	m_instance = this;
}


yosEigenState::yosEigenState(yosBasis *newBasis, std::istream & in):
		  persist_eigSt(newBasis, true),
		  myBasis(newBasis),
		  fixedPart(0),
		  isSearchInitialized(false),
		  fvec(0),
		  fjac(0)
{

	if (!readState(in))
	{
		throw CxBadValueError(__FILE__,__LINE__,"Could not read a yosEigenState from file");
	}
}

/*!
  \fn yosEigenState::yosEigenState(yosBasis *newBasis, double * eigenvec, double * eigenVal, int stateNo)
  \param newBasis Pointer to basis
  \param eigenvec array of eigenvector coefficents
  \param eigenVal Array of eigenvalues (is always double)
  \param stateNo Number of tate to be used (defaults to 0 which is the groundstate)
  \brief Constructor which uses the output of LLLhamiltonian::diagonalize 


 */

yosEigenState::yosEigenState(yosBasis *newBasis, 
		double * eigenvec,
		double * eigenVal,
		int stateNo)
: persist_eigSt((Basis *) newBasis, true),
  fixedPart(0),
  isSearchInitialized(false),
  fvec(0),
  fjac(0)



{
	if ( (!eigenvec) || (!eigenVal) )
	{
		init (newBasis, 0);
		return;
	}


	complex <double> *coefficents = new complex <double> [newBasis->dimension()];
	/*
    Put together a complex array from eigenvec
	 */

	int cIndex = 0;
	/*
    Calculate the start and stop point from the array as received from Operator::diagonalize()
    Assumes the matrix to be organized in C like <Re(a11) Im(a11)...Im(a1,dimension, ... re(dim), Im(dim))
	 */
	int startIndex = stateNo*2*newBasis->dimension()+1;
	int stopIndex = startIndex + newBasis->dimension()*2;


	for (int index=startIndex; index < stopIndex; index=index+2)
	{

		coefficents[cIndex] = complex<double>(eigenvec[index],eigenvec[index+1]);

		cIndex++;

	}

	m_instance = this;

}

/*!

 */

yosEigenState::yosEigenState( const yosEigenState &rhs)
: persist_eigSt(rhs),
  fixedPart(0),
  isSearchInitialized(false),
  fvec(0),
  fjac(0)
{
	glLogger.debug("Copy ctor of yosEigenState");
	if (&rhs != this)
	{
		initialized = false;
		electronVector = rhs.getElectronPositions();
		calculateFixedPart();

	}
	m_instance = this;
}


/*!
  \brief ctor from persist_eigSt and basis
  \param newBasis the basis
  \param oldState the persist_eigst in the given basis

 */
yosEigenState::yosEigenState( yosBasis *newBasis, 
		persist_eigSt &oldState):
		  persist_eigSt(oldState),
		  fixedPart(0),
		  isSearchInitialized(false),
		  fvec(0),
		  fjac(0)
{
	glLogger.info("Ctor of yosEigenState with persist_eigSt");
	if (newBasis == 0)
	{
		throw (CxErrors("Null basis pointer given at ",
				__FILE__, __LINE__));
	}
	glLogger.info("Dimension of basis is now (%d)",
			newBasis->dimension());
	for (int index = 0 ; index < 2*dimension(); index=index+2)
	{
		glLogger.info("Coefficents in init (%f), (%f)",
				oldState.outCoef(index), oldState.outCoef(index+1));

	}
	fixedPart = 0;
	zeroPositions.clear();
	vortexPositions.clear();
	myBasis = newBasis;
	initialized = true;
	noOfElectrons = newBasis->getNe()-1;
	aspectRatio = myBasis->getAspectRatio();
	glLogger.info("Cell aspect ratio is (%f)", aspectRatio);
	glLogger.info("Number of electron sites to be sampled =(%d)",
			noOfElectrons);
	// electronSites = new CxPosition [ noOfElectrons ];
	readElectronPositions("positions.dat");

	if (!calculateFixedPart())
	{
		throw (CxErrors ("Problems with calculateFixedPart",
				__FILE__, __LINE__));
	}
	m_instance = this;
}

/*!

 */
yosEigenState::yosEigenState (yosBasis *newBasis,
		persist_eigSt &oldState,
		std::vector<CxPeriodicPosition>& elecPos):
		  persist_eigSt(oldState),
		  isSearchInitialized(false),
		  fvec(0),
		  fjac(0)
{
	glLogger.info("Ctor of yosEigenState with persist_eigSt and elecPos");
	if (newBasis == 0)
	{
		throw (CxErrors("Null basis pointer given at ",
				__FILE__, __LINE__));
	}
	glLogger.info("Dimension of basis is now (%d)",
			newBasis->dimension());

	noOfElectrons = newBasis->getNe()-1;

	if (elecPos.size() < (unsigned int) noOfElectrons)
	{
		throw ( CxErrors("Bad number of electrons in vector",
				__FILE__, __LINE__));
	}

	//  electronSites = new CxPosition [ elecPos.size()];

	electronVector.clear();
	for (unsigned int index = 0; index < elecPos.size(); index++)
	{
		glLogger.info("yosEigenState ctor receives electronPos at (%f), (%f) index (%f)",
				elecPos[index].getXPosition(),
				elecPos[index].getYPosition(),
				&index);

	}
	electronVector = elecPos;
	fixedPart = 0;
	zeroPositions.clear();
	vortexPositions.clear();
	myBasis = newBasis;
	initialized = true;
	m_instance = this;

	aspectRatio = myBasis->getAspectRatio();
	if (!calculateFixedPart())
	{
		throw (CxErrors ("Problems with calculateFixedPart",
				__FILE__, __LINE__));
	}
}

/*!
  \brief ctor from persist_eigSt and basis
  \param newBasis the basis
  \param oldState the persist_eigst in the given basis
  \param positionFile the name of the file holding the positions to be used as electron sites
 */
yosEigenState::yosEigenState( yosBasis *newBasis, 
		persist_eigSt &oldState,
		string positionFile):
		  persist_eigSt(oldState),
		  fixedPart(0),
		  isSearchInitialized(false),
		  fvec(0),
		  fjac(0)
{
	glLogger.info("Ctor of yosEigenState with persist_eigSt");
	if (newBasis == 0)
	{
		throw (CxErrors("Null basis pointer given at ",
				__FILE__, __LINE__));
	}
	glLogger.info("Dimension of basis is now (%d)",
			newBasis->dimension());
	for (int index = 0 ; index < 2*dimension(); index=index+2)
	{
		glLogger.info("Coefficents in init (%f), (%f)",
				oldState.outCoef(index), oldState.outCoef(index+1));

	}
	fixedPart = 0;
	zeroPositions.clear();
	vortexPositions.clear();
	myBasis = newBasis;
	initialized = true;
	noOfElectrons = newBasis->getNe()-1;
	aspectRatio = myBasis->getAspectRatio();
	glLogger.info("Cell aspect ratio is (%f)", aspectRatio);
	glLogger.info("Number of electron sites to be sampled =(%d)",
			noOfElectrons);
	//  electronSites = new CxPosition [ noOfElectrons ];
	//  cerr << electronSites[0];
	readElectronPositions(positionFile);
	m_instance = this;
	if (!calculateFixedPart())
	{
		throw (CxErrors ("Problems with calculateFixedPart",
				__FILE__, __LINE__));
	}
}
/*!
 * 
 */


yosEigenState::yosEigenState(yosBasis *n_Basis, std::string stateFile, int stateNo):
		  persist_eigSt(n_Basis,true),
		  myBasis(n_Basis),
		  noOfElectrons(0),
		  aspectRatio(1.0),
		  fixedPart(0),
		  initialized(false),
		  isSearchInitialized(false),
		  fvec(0),
		  fjac(0)


{
	std::ifstream in(stateFile.c_str());
	//in >> this;
	m_instance = this;
} 



/*!
  \fn yosEigenState::~yosEigenState()
  \brief the destructor
 */


yosEigenState::~yosEigenState()
{
	glLogger.info("Entering destructor of yosEigenState");

	if (fixedPart != 0)
	{
		delete [] fixedPart;
		fixedPart = 0;
	}

	electronVector.clear();
	zeroPositions.clear();
	vortexPositions.clear();
	delete [] fvec;
	fvec = 0;
	delete [] fjac;
	fjac = 0;
	glLogger.info("Leaving destructor of yosEigenState");
	m_instance= 0;
	return;

}

/*!
  \fn bool yosEigenState::init(yosBasis *newBasis, complex<double> *coeff)
  \brief A private helper function to initialize the class
  \param newBasis The Basis of yoshioka states
  \param coefficents the coefficents (the so called eigenvector form diagonalization routines)
  \return true if initialization was successfull
 */


bool yosEigenState::init(yosBasis *newBasis, complex<double> *coefficents)
{

	myBasis = newBasis;
	noOfElectrons = myBasis->getNe() -1;
	aspectRatio = myBasis->getAspectRatio();
	glLogger.info("Aspect ratio is (%f)", aspectRatio);
	glLogger.info(" init: dimension = (%d)", newBasis->dimension());
	glLogger.info("Number of electrons -1 =(%d)", noOfElectrons);
	if (coefficents)
	{
		inAllCoef (coefficents);
		zeroPositions.clear();
		vortexPositions.clear();
		initialized = false;

		if (!calculateFixedPart())
		{
			throw (CxErrors ("Problems with calculateFixedPart",
					__FILE__, __LINE__));
		}
	}
	return true;
}

/*!
  \fn void  yosEigenState::addElectronPositions(CxPeriodicPosition *newElectronSites)
  \author Christian Mueller
  \date  04 Jun 03
  \return nothing 
  \param newElectronSites an array of Ne-1 electron positions 
  \brief Sets the electrons positions explicitly. If there are already positons set they will be overwrtten
 */
/*
  Pre	: 
  Post	: 

 */

void yosEigenState::addElectronPositions(CxPeriodicPosition *newElectronSites)
{
	// Test of Parameter
	if (newElectronSites == 0)
	{
		throw CxNullPointerError(__FILE__, __LINE__);
	}


	// Allocation of Resources

	for (unsigned int index = 0; index < myBasis->getNe(); index++)
	{
		if (newElectronSites == 0)
		{
			throw CxNullPointerError(__FILE__, __LINE__);
		}

		electronVector[index].setXPosition(newElectronSites[index].getXPosition());
		electronVector[index].setYPosition(newElectronSites[index].getYPosition());
	}
	initialized = true;


}




/*!
  \fn bool yosEigenState::calculateFixedPart(void)
  \brief Calculates the fixed part of the wavefunction.
  \return true if everything went fine, false if no electronpostions have been specified
  Stores it in fixedPart. 
 */

bool yosEigenState::calculateFixedPart(void)
{
	glLogger.info("Entering bool yosEigenState::calculateFixedPart(void)");
	int Ne = myBasis->getNe();
	if ( !initialized)
	{
		glLogger.error("Could not start calculatefixedpart");
		return false;
	}

	if (Ne < 2)
	{
		glLogger.error("Only 1 electron provided");
		return false;
	}
	if (fixedPart)
	{
		delete [] fixedPart;
		fixedPart = 0;
	}
	fixedPart = new complex <double> [myBasis->getNm()];
	glLogger.info("Generated fixepart array of size (%d)",myBasis->getNm());
	glLogger.info("dimension of basis is (%d)",dimension());
	/*
    Work part starts here. Basically many determinant calculations
	 */
	complex <double> *slaterMat = new complex <double>[ noOfElectrons * myBasis->getNe() ];
	dComplex *subMatrix = new dComplex [(Ne-1)*(Ne-1)];
	for (int index = 0 ; index < dimension(); index=index+1)
	{
		/*
	Loop over eigenvectors yosStates
		 */
		complex <double> testEl (outCoef(2*index), outCoef((2*index)+1));
		glLogger.debug("Index is %d ", index);
		glLogger.debug(*(myBasis->getState(index)));

		if ( (testEl.real() == 0.0)
				&& ( testEl.imag()== 0.0 ) )
		{
			glLogger.info("Zero Element found, continuing");
			continue;
		}

		int slaterIndex = 0;
		for (int eIndex = 0 ;
				eIndex < Ne;
				eIndex ++)
		{
			/*
	    Loop over occNumbers in yosState 
	    (state is  retrieved via getState(index)
			 */
			glLogger.debug(*(myBasis->getState(index)));
			int wfIndex = myBasis->getState(index)->j(eIndex);
			for (int eiIndex = 0;
					eiIndex < Ne-1;
					eiIndex++)
			{
				/*
		Loop  over < Ne - 1
				 */

				double x = electronVector[eiIndex].getXPosition();
				double y = electronVector[eiIndex].getYPosition();
				slaterMat[slaterIndex] = myBasis->getState(index)->SPwaveFct(wfIndex, x,y);
				slaterIndex++;
			} //end eIndex


		} // end index
		glLogger.debug("Slater matrix  constructed");
		/*
	Now Constructing the quadratic matrix to be determinanted
	"Streichdeterminanten". Rowindex is the line to be cut
		 */
		int slaterSize = slaterIndex;
		//      std::complex<double> determinante;



		for (int rowIndex = 0; rowIndex <  Ne ; rowIndex++)
		{
			/*
	     Loop over rows to be cut out (each row shall be affected)
	     Start and stop of row to be cut out from slaterMat
			 */
			slaterIndex =0;
			int startIndex = rowIndex * (Ne-1);

			int stopIndex;
			stopIndex = startIndex + Ne -1;
			int subIndex = 0;
			while( slaterIndex < slaterSize)
			{
				/*
		 Copy a row
				 */

				if (slaterIndex == startIndex)
				{
					// Skip the row
					slaterIndex = stopIndex;
					//	       glLogger.debug( "Found the row to be cut");
					continue;
				}
				subMatrix[subIndex] = slaterMat[slaterIndex];
				//subMatrix[subIndex].imag = slaterMat[slaterIndex].imag();



				subIndex++;
				slaterIndex++;
			}

			/*
	     Now we have the cut down matrix in
	     subMatrix. Get determinant of it!
			 */
			std::complex<double> tempDeterminante;
			if (calculateComplexDeterminant(Ne-1, subMatrix, tempDeterminante))
			{
				delete [] subMatrix;
				delete [] slaterMat;
				throw CxErrors ("Error in diagonalization, called from",
						__FILE__, __LINE__);
			}


			/*
	    Put it into the element state[rowIndex] of fixedPart
			 */
			int fixedPosition = myBasis->getState(index)->j(rowIndex);

			if ( (floor((float)rowIndex/2)*2 )!= rowIndex)
			{
				tempDeterminante = tempDeterminante*-1.0;
			}

			fixedPart[fixedPosition] = tempDeterminante*testEl + fixedPart[fixedPosition];


		} // end rowIndex

	} // End of Loop over eigenStates
	/*
     Now multiplying with sqrt(NE) (factor from slaterderterminante)
	 */
	double vorFaktor = sqrt((float)factorial(Ne));
	for (unsigned int index = 0; index < myBasis->getNm(); index++)
	{
		fixedPart[index] = fixedPart[index]/vorFaktor;
		glLogger.info("FixedPart returns fixedPart[%d] = (%f) + i(%f)", index,fixedPart[index].real(),fixedPart[index].imag() );
	}
	delete [] subMatrix;
	delete [] slaterMat;
	glLogger.info("Successfully returning from calculateFixedPart");
	return true;
	/*


  This routine has been tested and shows the same results as fixedpart.f in the
  Big-Mac code. 03-07-2003 CM


	 */

}



/*!
  \brief sets the positions of the sample electrons(must be Ne -1)!
  \param numberOfElectrons
  \param newPositions Array of Electron positions, size must be numberofElectrons
  \return false if something failed

  throws errors for serious errors

 */
bool yosEigenState::setElectronPositions(int numberOfElectrons, CxPeriodicPosition *newPositions)
{
	if (!reset()) return (false);
	if (noOfElectrons-1 != numberOfElectrons)
	{
		//! \todo throw something here
		return (false);
	}
	if (0== newPositions)
	{
		//! \todo throw something here
		throw (CxNullPointerError(__FILE__,__LINE__));

	}
	for (int index = 0; index < numberOfElectrons; index++)
	{
		electronVector[index] = newPositions[index];
		glLogger.info("Electron no (%d) at position (%d), (%d)",
				newPositions[index].getXPosition(),
				newPositions[index].getYPosition());
	}
	initialized = true;
	return (!calculateFixedPart());

} 

/*!

 */
bool yosEigenState::setElectronPositions (std::vector<CxPeriodicPosition>& inPos)
{
	if ((int)inPos.size() != getNe() -1)
	{
		glLogger.warning("Bad number sample electrons passed. Expected (%d), received (%d)", (getNe()-1), inPos.size());
		return (false);
	}
	if (!reset())
	{
		glLogger.warning("Reset returned false in setElectrons <vector>");
		return (false);
	}

	electronVector.clear();
	electronVector = inPos;
	initialized = true;

	return (!calculateFixedPart());
}



/*!
  \fn void yosEigenState::readElectronPositions (const char* fileName)
  \param fileName guess what
  \brief Reads electron positions from a file in a suitable format. 
 */

void yosEigenState::readElectronPositions (string fileName)
{
	if (fileName == "")
	{
		throw (CxNullPointerError(__FILE__,__LINE__));
	}


	ifstream inFile(fileName.c_str());
	if (!inFile)
	{
		throw (CxFileNotFoundError(fileName.c_str(), __FILE__,__LINE__));
		return;
	}
	/*
    Read two double separated by C-style newline
	 */
	double
	tempXVal =0.0,
	tempYVal =0.0;
	if (!reset()) return;


	for (int index =0; index < noOfElectrons; index ++)
	{
		if (inFile)
		{
			inFile >> tempXVal;
		}
		else
		{
			return;
		}
		if (inFile)
		{
			inFile >> tempYVal;
		}
		else
		{
			return;
		}

		electronVector.push_back( CxPeriodicPosition(tempXVal, tempYVal));

		glLogger.info(" Electron Nr. (%d)  position at x=(%f), y=(%f)",
				index, tempXVal, tempYVal);

	}

	if (fixedPart != 0)
	{
		delete [] fixedPart;
	}
	if (! calculateFixedPart())
	{
		throw (CxErrors ("Problems with calculateFixedPart",
				__FILE__, __LINE__));
	}


	inFile.close();

}


/*!
  \fn  double yosEigenState::getWindingNumber(CxPeriodicPosition &point)
  \brief Calculates the windingnumber at a given point
  \return Winding number around a given point
  \param point [in] to be looked at
 */

double yosEigenState::getWindingNumber(CxPeriodicPosition &point ) const
{
	glLogger.debug(" yosEigenState::getWindingNumber(CxPeriodicPosition &point )");
	glLogger.debug(" Position is (%f), (%f)",
			point.getXPosition(),
			point.getYPosition());
	if (fixedPart == 0)
	{
		throw( CxErrors("yosEigenState::getWindingNumber not properly initialized"));
		return (0);
	}



	/*
    Check vorticity. Basicaly a line integral of the sp-function around the zero.
    Number of iterations as in Danielas Code
	 */
	complex <double> curvInt(0,0);
	for (int iIndex=0; iIndex < curvIntStepsMax; iIndex++)
	{
		double theta = iIndex*2*M_PI/curvIntStepsMax;
		double dx = point.getXPosition() + curvRadius*cos(theta);
		double dy = point.getYPosition() + curvRadius*sin(theta);
		//	   New point to be defined
		CxPeriodicPosition newPos(dx,dy, point.getAspectRatio());
		complex <double> tempVal(0,0);
		double dfyReal =  sumofdySPWavefunctions(0, newPos).real();
		double dfyImag = sumofdySPWavefunctions(0, newPos).imag();
		tempVal = sumofdxSPWavefunctions(0,newPos) - complex<double>(0,dfyReal);
		tempVal = tempVal + dfyImag;
		tempVal = tempVal / (sumofSPWavefunctions(0, newPos) - sumofSPWavefunctions(0, point));

		tempVal = tempVal *exp(complex<double>(0, theta));
		curvInt = curvInt + tempVal;

	}
	curvInt = curvInt*curvRadius/(std::complex <double>(2*curvIntStepsMax,0));
	glLogger.info("Curve integral yields (%f, %f)", curvInt.real(), curvInt.imag());

	/*
     Check Vortex
	 */
	if (curvInt.imag() < vortexTolerance)
	{
		if (fabs(rint(curvInt.real())-curvInt.real()) <  vortexTolerance)
		{
			glLogger.info("Vortex found with vorticity:(%f)", rint(curvInt.real()));
			return (rint(curvInt.real()));
		}
		else
		{
			return (rint(curvInt.real()));
		}
	}
	else
	{
		glLogger.info("Zero which is not a Vortex found");
	}


	return (0.0);
}






/*!
  \fn double yosEigenState::sampleElectronSite(CxPeriodicPosition &point ) const
  \param point Position to be looked up (usually the site where an electron resides)
  \return -1 if no zero is found, the abs of the reduced quasi sp wavefunction otherwise
 */

double yosEigenState::sampleElectronSite(CxPeriodicPosition &point ) const
{
	glLogger.debug(" Entering yosEigenState::sampleElectronSite");
	glLogger.debug("Electronposition at x=(%f), y=(%f)",
			point.getXPosition(),
			point.getYPosition());

	double retVal = 0;
	complex <double> sumWfct(0, 0);
	int groundState = 0;
//writeReducedWfct(groundState,"reducedFile2.txt",0.005);
	sumWfct = sumofSPWavefunctions(groundState, point);

	glLogger.info(", wavefct = (%f), (%f)",  sumWfct.real(), sumWfct.imag());
	if (abs(sumWfct) > 1e-09)
	{
		glLogger.info("Wavefunction is not zero. Returning");
		return (-1.0);
	}
	if (abs(sumWfct) < 1e-09)
	{

		glLogger.info("Found a promising zero");
		return (0.0);

	}
	glLogger.debug("Leaving sampleElectronSitres");
	return (retVal);
}
/*!

\fn complex <double>yosEigenState::sumofSPWavefunctions(const int  stateNo, CxPeriodicPosition &point) const
\brief sum the SPWavefct contributions
\param point where the function shall be evaluated
\param stateNo Number of state to be used
 */
complex <double> yosEigenState::sumofSPWavefunctions(const int  stateNo, CxPosition &point) const
{
	complex <double> retVal(0,0);

	//  glLogger.debug("Entering yosEigenState::sumofSPWavefunctions");
	for (unsigned  int jIndex = 0; jIndex < myBasis->getNm(); jIndex++)
	{
		complex <double> tempVal;
		tempVal = myBasis->getState(stateNo)->SPwaveFct(jIndex,
				point.getXPosition(),
				point.getYPosition());

	glLogger.debug("Loop over j's, j=(%d), value = ((%f), (%f))",
	jIndex, tempVal.real(), tempVal.imag());

		retVal = tempVal * fixedPart[jIndex]  + retVal;

	}

    glLogger.debug("Leaving sumofSPWavefunctions"
    "with (%f), (%f)",  retVal.real(), retVal.imag());
	return (retVal);

}


/*!

 */
complex <double> yosEigenState::sumofSPWavefunctions(const int  stateNo, double xPos, double yPos) const
{
	CxPosition position(xPos, yPos);
	return (sumofSPWavefunctions(stateNo, position));

}

/*!
  \fn complex <double> yosEigenState::sumofdySPWavefunctions(const int  stateNo, CxPosition &point) const
  \brief   y-Derivative of the summed up two-dimensional wavefunction
  \param point where the function shall be evaluated
  \param stateNo Number of state to be used
 */

complex <double> yosEigenState::sumofdySPWavefunctions(const int  stateNo, CxPosition &point) const
{
	complex <double> retVal(0,0);


	for (unsigned  int jIndex = 0; jIndex < myBasis->getNm(); jIndex++)
	{

		retVal = myBasis->getState(stateNo)->SPwaveFctdy(jIndex,
				point.getXPosition(),
				point.getYPosition())*
				fixedPart[jIndex]  + retVal;
	}
	return (retVal);

}



/*!
  \fn complex <double> yosEigenState::sumofdxSPWavefunctions(const int  stateNo, CxPosition &point) const
  \brief   x-Derivative of the summed-up Wavefunctions
  \param point where the function shal lbe evaluated
  \param stateNo Number of state to be used
 */

complex <double> yosEigenState::sumofdxSPWavefunctions(const int  stateNo, CxPosition &point) const
{
	complex <double> retVal(0,0);


	for (unsigned int jIndex = 0; jIndex < myBasis->getNm(); jIndex++)
	{

		retVal = myBasis->getState(stateNo)->SPwaveFctdx(jIndex,
				point.getXPosition(),
				point.getYPosition())*
				fixedPart[jIndex]  + retVal;
	}

	return (retVal);

}


/*!
  \fn void yosEigenState::findZeros(void)
  \brief Trys to find the zeros. 

  Needs the fixedPartcalculation and the electronSites


 */

void yosEigenState::findZeros(void)
{
	glLogger.info("Entering yosEigenState::findZeros");

	if (fixedPart == 0)
	{

		glLogger.warning("yosEigenState::findZeros receives null pointer, leaving here");
		return;
	}


	for (int elecIndex = 0; elecIndex < noOfElectrons; elecIndex++)
	{
		/*
	 Loop over electronSites, if wfct is zero there, 
	 calculate Windingnumber
		 */
		if ( sampleElectronSite(electronVector[elecIndex]) != -1.0)
		{
			double phase  = getWindingNumber(electronVector[elecIndex]);
			glLogger.info("Electron Pos sampling :Zero at (%f), (%f) with wind.number (%f)",
					electronVector[elecIndex].getXPosition(),
					electronVector[elecIndex].getYPosition(),
					phase);

			vortexPositions.addPosition( CxVortex(electronVector[elecIndex].getXPosition(),
					electronVector[elecIndex].getYPosition(),
					(int) rint(phase)));
			zeroPositions.addPosition( CxPeriodicPosition(electronVector[elecIndex].getXPosition(),
					electronVector[elecIndex].getYPosition(),
					1.0, 1.0));
		}
		else
		{
			glLogger.info("No zero at Electronposition (%d)",elecIndex);
		}


	}

	/*
    Call search over whole cell now
	 */

	//  searchAllZeros();
	freezeInstance();
#ifdef NAG
	glLogger.info("Calling searchNag");
	searchNag();
	glLogger.info("returning from searchNag");
#else
#ifdef MKL
	glLogger.info("Calling searchIntel");
	searchIntel();
	glLogger.info("returning from searchIntel");
#endif
#ifdef GSL
	glLogger.info("searchGSL called");
	searchGsl();
#endif
#ifdef LAPACK
	glLogger.info("searchGSL called");
	searchGsl();
#endif
#endif

}	


/*!

\fn int yosEigenState::writeReducedWfct(const char* fileName) const
\brief  Writes  the quasi-sp wavefunction into a file
\param fileName The name of the file to be written to
\return 0 if everything went fine, -1 otherwise
 */

int yosEigenState::writeReducedWfct(int stateNo, std::string fileName, double gridSpacing =0.0) const
{
	if (glLogger.getLogLevel()!= DEBUG)
	{
		return 1;
	}
	int retVal = 0;
	if (fileName.empty())
	{
		return (-1);
	}
	if (fixedPart == 0)
	{
		return (-1);
	}
	ofstream theFile(fileName.c_str());
	if (!theFile)
	{
		return (-1);
	}
	/*
    File is open and ready to write
    Loop over area 
	 */
	theFile << "# X \t Y \t RealPart \t Imaginary Part\n";
	for (double xPos = 0; xPos < 1.0 ; xPos = xPos+gridSpacing)
	{
		for (double yPos =0;
				yPos < myBasis->getState(0)->getAspect()*1.0;
				yPos = yPos + gridSpacing)
		{
			std::complex<double> tempVal = sumofSPWavefunctions(stateNo,xPos, yPos);
			theFile << xPos << " \t "<<  yPos << " \t ";
			theFile << tempVal.real()<< " \t";
			theFile << tempVal.imag()<< "\t";
			theFile << abs(tempVal)<<"\n";

		} // End of y-loop
		theFile <<std::endl;
	}// end of x-Loop





	return retVal;
}
#ifdef NAG
/*!

\brief does zero searching uses the NAG library (hopefully). Uses the same trick as Karel in matrix_sparse_real.cpp
\return -1 on error, number of zeros found otherwise
 */

int yosEigenState::searchNag(void)
{

	/*
    Set up the NAG variables, names are the same as in the NAG specification
    Integers first, the doubles, arrays of integers and arrays of doubles
	 */
	glLogger.info("Entering searchNAG");
	/*  int
      m, n, nclin, ncnln,lda, ldcj, ldfj, 
      ldr, iter, liwork, lwork, ifail;   
	 */

	zeroPositions.clear();
	int
	m = 2,
	n = 2,
	nclin =  0,
	ncnln = 0,
	lda = 1,
	ldcj = 1,
	ldfj = 2,
	ldr = 2,
	iter = 0,
	liwork =6,
	lwork = 50,
	ifail =0 ;

	int *istate = new int [(n)+(nclin)+(ncnln)];
	int *iwork = new  int [liwork];
	int *iuser;
	*iuser = 1;
	//  iuser = (int *) this;
	/*
    Double arrays next
	 */
	double *a = new double [(lda)];
	double *bl = new double [n + nclin + ncnln];
	setDoubleArray(bl, n+nclin+ncnln, 0.0);

	double *bu = new double [n+nclin+ncnln];

	setDoubleArray(bu, (n+nclin+ncnln), 0.0);
	bu[1] = 1.0;
	bu[0] = 1.0;

	double *y = new double[(m)];

	setDoubleArray(y, m, 0.0);


	double *c = new double [1];
	c[0] = 0.0;
	double *cjac = new double[(ldcj)];
	double *fjac = new double[(ldfj)*(n)];
	double *clamda = new double [(n)+(nclin)+(ncnln)];
	double * f = new double [(m)];
	//  double *fjac = new [(*ldfj)*(*n)];
	double *r = new double [(ldr)*(n)];
	double *x = new double[n];


	double *work = new double [lwork];
	//  double user[1];
	double *user = new double [1];
	double objf = 0;
	/*
    function pointers
	 */

	pObjfun = &objfun;
	//  pConfun = &confun;

	pConfun = &e04udm_;


	/*
     Set the positions accordingly
	 */


	int trialCount = 0;
	srand(1234);
	while ( (zeroPositions.size() <  (unsigned int)myBasis->getNm()) &&
			(trialCount < 100 ))
	{
		trialCount++;
		ifail = 1;
		x[0] = 0;
		x[1] = 0;
		while ((x[0]< curvRadius ) || x[0] > 1.0 -curvRadius)
		{

			x[0] = (double) rand() / (double) RAND_MAX;
		}

		while ((x[1]< curvRadius ) || x[1] > 1.0 -curvRadius)
		{
			x[1] = (double) rand() / (double) RAND_MAX;
		}
		//! \todo extend to atob
		/*
	Use with OLD NAG (pre MARK 18, use e04usf_ with newer versions
		 */

		freezeInstance();
		// suppress output from NAG
		int nerr = -1;
		int iflag1 = 1;
#ifdef NAG
		x04aaf_(&iflag1, &nerr);
#ifdef NAG_19 
		e04usf_(&m, &n, &nclin, &ncnln, &lda, &ldcj, &ldfj, &ldr, a, bl, bu, y,
				pConfun, pObjfun, &iter, istate, c, cjac, f, fjac,
				clamda, &objf, r, x, iwork, &liwork, work, &lwork,iuser,
				user, &ifail);
#else   
		// That could also be e04unf!!

		e04unf_(&m, &n, &nclin, &ncnln, &lda, &ldcj, &ldfj, &ldr, a, bl, bu, y,
				pConfun, pObjfun, &iter, istate, c, cjac, f, fjac,
				clamda, &objf, r, x, iwork, &liwork, work, &lwork,iuser,
				user, &ifail);
#endif
/*
	e04unf_(&m, &n, &nclin, &ncnln, &lda, &ldcj, &ldfj, &ldr, a, bl, bu, 
	pConfun, pObjfun, &iter, istate, c, cjac, f, fjac, 
	clamda, &objf, r, x, iwork, &liwork, work, &lwork,iuser, 
	user, &ifail);
 */
#endif
		glLogger.info("NAG returns (%d) after (%d) iterations", ifail, iter);
		double aToB = electronVector[0].getAspectRatio();
		if (ifail == 0)
		{
			glLogger.info("Position of zero is x = (%f), y=(%f)", x[0], x[1]);

		}
		else {
			continue;
		}

		CxPeriodicPosition tempPosition(x[0], x[1], aToB );
		tempPosition.setEpsilon(1e-05);



		double newPhase = getWindingNumber(tempPosition);
		if (newPhase != 0.0)
		{
			//	      CxPosition *myPos = new CxPosition(tempPosition);
			glLogger.info("Adding a position with x=(%f), y=(%f), xSize=(%f), ySize=(%f)",
					tempPosition.getXPosition(),tempPosition.getYPosition(),
					tempPosition.get_xCellSize(), tempPosition.get_yCellSize());
			zeroPositions.addPosition(tempPosition);;
			CxVortex tempVortex(tempPosition,(int) rint(newPhase));
			vortexPositions.addPosition( tempVortex);

		}
	}



	glLogger.warning("Found (%d) zeros in (%d) loops",
			zeroPositions.size(),trialCount);
	glLogger.info("Dumpinng vortices ");
	std::vector<CxVortex> vorPos = vortexPositions.getPositions();
	std::vector<CxPeriodicPosition> zeroPos = zeroPositions.getPositions();
	std::vector<CxPeriodicPosition>::iterator mZIt = zeroPos.begin();
	std::vector<CxVortex>::iterator mIt = vorPos.begin();
	int count = 0;
	while (vorPos.end() != mIt)
	{
		glLogger.info("Vortex Position No.%d, x=(%f), y=(%f), wind=(%f)", count,
				mIt->getXPosition(),
				mIt->getYPosition(),
				mIt->getWindingNumber());
		mIt++;
		count++;
	}
	count = 0;
	while (mZIt != zeroPos.end())
	{

		glLogger.info("Zero position No.%d, x=(%f), y=(%f)", count,
				mZIt->getXPosition(),
				mZIt->getYPosition());
		count++;
		mZIt++;
	}

	glLogger.info("End dump of vortices here ");

	/*
    Clean up 
	 */
	delete [] istate;
	istate = 0;
	delete [] iwork;
	iwork = 0;
	delete [] a;
	a = 0;
	delete [] bl;
	bl = 0;
	delete [] bu;
	bu = 0;

	delete [] y;
	y = 0;

	delete [] c;
	c = 0;
	delete [] cjac;
	cjac = 0;
	delete [] f;
	f = 0;
	delete [] fjac;
	fjac = 0;
	delete [] x;
	x = 0;
	delete [] r;
	r = 0;
	delete [] clamda;
	clamda = 0;
	delete [] work;
	work = 0;
	delete [] user;
	user = 0;
	return (42);
}

#endif


const int yosEigenState::curvIntStepsMax = 200;
const double yosEigenState::curvRadius = 1e-05;

const double yosEigenState::vortexTolerance = 1e-06;

/*

External defined functions needed for use of external C/Fortran type libraris

 */



extern "C" {

dComplex summedWavefunction(double x, double y, dComplex *prefactors, int Nm, int Ne)
{
	if (x && y && prefactors && Ne > 0)
	{
	}
	dComplex retVal(0.0,0.0);


	for (int jIndex =0; jIndex < Nm ; jIndex++)
	{
		/*
	  Sum up the sp-yoshioka type functions
		 */
	}

	return retVal;

}

}

bool complexLarger(  complex<double> lhs, complex<double> rhs)
{
	return (!complexSmaller(rhs, lhs));
}

bool complexSmaller( complex<double> lhs, complex<double> rhs)
{

	if (lhs.real() < rhs.real())
	{
		if (lhs.imag() < rhs.imag())
		{
			return (true);
		}
		else
		{
			return (false);
		}
	}
	else
	{
		return false;
	}

}


/*!
  \brief function to be minimized by NAG routine
  \param iuser holds the user defined passed values. In this case a pointer to class yosEigenState.
  \remark This is bad practice but it works...
 */
void objfun(int * mode, int * m, int *n, int *ldfj,  
		double * x, double * f, double * fjac, int * nstate, int * iuser, double * user)
{
	// Set precision of output stream

	//  c << "\t mode is "<< *mode << "\n";
	yosEigenState *myClass = yosEigenState::getInstance();
	//(yosEigenState *)iuser;
	//  cerr << "Position received from NAg are x= ("<<x[0] <<") y= ("<<x[1]<<"\n";
	if (m && n && ldfj && nstate && iuser && user  )
	{
		// Possibly throw something here, mainly to keep the compiler happy
	}
	if ( (x[0] < 0)
			|| x[1] < 0)
	{
		/*
	 Map it back to unit cell!
		 */
		if (x[0] < 0)
		{
			x[0] = x[0]+1.0;
		}
		if (x[1] < 0)
		{
			x[1] = x[1] + myClass->getAspectRatio();
		}
	}
	if (x[0] > 1)
	{
		x[0] = x[0]-1.0;
		glLogger.debug( "Mapped function into unit cell with x= (%d) y= (%d",x[0],x[1]);
	}
	if (x[1] > myClass->getAspectRatio())
	{
		x[1] = x[1] - myClass->getAspectRatio();
		glLogger.debug( "Mapped function into unit cell with x= (%d) y= (%d",x[0],x[1]);

	}




	CxPeriodicPosition position(x[0], x[1]);



	complex <double> fctVal = myClass->sumofSPWavefunctions(0,x[0],x[1]);
	f[0] = fctVal.real();
	f[1] = fctVal.imag();
	/*
    N.B. the array storage difference between C and Fortran
	 */

	if ((*mode) == 0)
	{

		return;
	}

	complex <double> fctdx = myClass->sumofdxSPWavefunctions(0,position);
	complex <double> fctdy =  myClass->sumofdySPWavefunctions(0,position);

	fjac[0] = fctdx.real();
	fjac[2] = fctdy.real();
	fjac[1] = fctdx.imag();
	fjac[3] = fctdy.imag();

}

void jacobiObjFun(int* m, int *n, double *x, double *f)
{

	yosEigenState *myClass = yosEigenState::getInstance();


	if ( (x[0] < 0)
			|| x[1] < 0)
	{
		/*
	 Map it back to unit cell!
		 */
		if (x[0] < 0)
		{
			x[0] = x[0]+1.0;
		}
		if (x[1] < 0)
		{
			x[1] = x[1] + myClass->getAspectRatio();
		}
	}
	if (x[0] > 1)
	{
		x[0] = x[0]-1.0;
		glLogger.error( "Mapped function into unit cell with x= (%d) y= (%d",x[0],x[1]);
	}
	if (x[1] > myClass->getAspectRatio())
	{
		x[1] = x[1] - myClass->getAspectRatio();
		glLogger.error( "Mapped function into unit cell with x= (%d) y= (%d",x[0],x[1]);

	}




	CxPeriodicPosition position(x[0], x[1]);



	complex <double> fctVal = myClass->sumofSPWavefunctions(0,x[0],x[1]);
	f[0] = fctVal.real();
	f[1] = fctVal.imag();
	/*
    N.B. the array storage difference between C and Fortran
	 */


	CMDCOMPLEX fctdx = myClass->sumofdxSPWavefunctions(0,position);
	CMDCOMPLEX fctdy =  myClass->sumofdySPWavefunctions(0,position);

	f[0] = fctdx.real();
	f[2] = fctdy.real();
	f[1] = fctdx.imag();
	f[3] = fctdy.imag();

}
/*!
  Just a dummy for linkage, later use the NAG dummy

 */
/*
  void confun(int *a, int *b, int *c, int *d, int *e, double *f, double *g, double *h)
  {
  return;

  }
 */
/*!

\brief convenience function to zero an double array
 */

void setDoubleArray(double *array, int length, double value)
{

	for (int index = 0; index < length ; index ++)
	{

		*array = value;
		array++;

	}

}
/*!
  \brief Shalloow accessor function. Be carefull, no copying is done, all changes wl affect the actual values in the class.

 */

std::vector<CxPeriodicPosition>  yosEigenState::getShallowZeros(void) const
{

	glLogger.info("entering  yosEigenState::getShallowZeros");
	return(zeroPositions.getPositions());

}

std::vector<CxVortex> yosEigenState::getShallowVortices(void) const
{

	return(vortexPositions.getPositions());
}
/*!
  \brief Access function for aspect ratio

 */
double yosEigenState::getAspectRatio(void) const
{
	return (aspectRatio);
}
/*!
  \fn  vector <CxPeriodicPosition>  yosEigenState::getVectorOfZeros(void )
  \author Christian Mueller
  \date  27 Aug 03
  \return 
  \param void 
  \brief 
 */
/*
  Pre	: 
  Post	: 

 */

std::vector <CxPeriodicPosition> yosEigenState::getVectorOfZeros(void ) const
{
	// Test of Parameter
	// Definition and Initialisation of Variables
	std::vector <CxPeriodicPosition> retVal =zeroPositions.getPositions() ;
	// Allocation of Resources


	// Work Part


	return (retVal);

}


/*!
  \brief retrieves the electronPositions
  \return A vector containing the electron Positions
 */
std::vector <CxPeriodicPosition> yosEigenState::getElectronPositions(void) const
{
	std::vector <CxPeriodicPosition> retVal(electronVector);;

	return (retVal);

}
/*!

\brief moves electrons in a star like way away/towards the impurity
 */

void yosEigenState::moveElectrons(const double newFac, const CxPeriodicPosition & impurity)
{
	std::vector<CxPeriodicPosition>::iterator elecIt = electronVector.begin();
	while (elecIt != electronVector.end())
	{

		glLogger.info("Positions before move %f),(%f)",
				elecIt->getXPosition(),
				elecIt->getYPosition());

		CxPeriodicPosition differenz;
		CxPeriodicPosition electron = *(elecIt);
		elecIt++;
	}
	if (!calculateFixedPart())
	{
		throw (CxErrors(__FILE__, __LINE__));
	}



}

/*!
  \brief reset function to be used when resetting electron positions etc.
  \return false if anythiung strange happened
 */

bool yosEigenState::reset(void)
{
	vortexPositions.clear() ;
	zeroPositions.clear();
	electronVector.clear();
	noOfElectrons = myBasis->getNe()-1;
	if (fixedPart) {
		delete [] fixedPart;
		fixedPart = 0;
	}
	return (true);
}
/*!
 * 
 */
int yosEigenState::getNe(void)const
{
	if (!myBasis)
	{
		return (0);
	}
	else
	{
		return (myBasis->getNe());
	}
}

/*! \fn float yosEigenState::getWeight(void) 
 * \brief return weight of the state
 * \return the weight
 */
float yosEigenState::getWeight(void) 
{
	float retVal = 0;
	if (!fixedPart)
	{
		if (!calculateFixedPart())
		{
			throw CxErrors("No fixed electrons available");
		}
	}
	for (unsigned int index = 0; index < myBasis->getNm(); index++)
	{
		retVal = retVal + abs(fixedPart[index]);
	}
	return (retVal);
}


/*!
  input / output operators
 */
std::ostream & operator<< (std::ostream &out, const yosEigenState &eigst) 
{

	// if (eigst.isComplex())
	glLogger.error ("Saving complex vector");

	int nCount = eigst.isComplex() ? 2*eigst.dimension() : eigst.dimension();

	if ( eigst.is_binary() )
	{
		// convert doubles to unformatted binary string
		double d = eigst.getEn();
		out.write( (char*) &d, sizeof(double)/sizeof(char) );
		out.write ( (char*) eigst.getCoefPtr(), nCount*sizeof(double)/sizeof(char) );

	}
	else
	{
		out << eigst.getEn()<<std::endl;
		for (int i=0; i< nCount; i++)
		{
			out << ' ' << eigst.outCoef(i);
			if (i == 0)
			{
				std::cerr << eigst.outCoef(i) << std::endl;
			}
		}

		out << '\n';
	}

	return out;
}


/*
 * read next state from stream
 * 
 */
bool yosEigenState::readState(std::istream &in)
{
	if ( !in.good() || in.eof() )
	{
		return (false);
	}
	if (is_binary())
	{
		//BINARY
		int nCount = isComplex() ? 2*dimension() : dimension();
		double *pCoefs = new double [nCount];
		in.read( (char*) &En, sizeof(double)/sizeof(char) );
		setEn(En);
		in.read ( (char*) pCoefs, sizeof(double)*nCount/sizeof(char) );
		if (isComplex())
		{
			inAllCoef( (std::complex<double>*)pCoefs );
		}
		else
		{
			inAllCoef(pCoefs);
		}
		delete [] pCoefs;



	}
	else
	{
		// ASCII
		in >> En;
		setEn(En);
		if (isComplex())
		{

			CMDCOMPLEX cCoeff;
			for (int i = 0; i<dimension(); i++)
			{
				in >>cCoeff;
				inCoef(2*i,cCoeff.real());
				inCoef(2*i+1, cCoeff.imag());
			}
		}
		else
		{
		int nCount =  dimension();
		double dCoefi;
		for (int i = 0; i < nCount; i++)
		{
			in >> dCoefi;
			inCoef(i, dCoefi);
		}
		}

	}

	return (true);

}


/*!
  \fn std::istream & operator>> (std::istream &in, persist_eigSt &eigst) 
  \brief Streams the content of state to stream
 */
std::istream & operator>> (std::istream &in, yosEigenState &eigst) 
{
	double En = 0.0;
	// reading complex values
	if (eigst.is_binary())
	{
		int nCount = eigst.isComplex() ? 2*eigst.dimension() : eigst.dimension();
		double *pCoefs = new double [nCount];



		in.read( (char*) &En, sizeof(double)/sizeof(char) );
		eigst.setEn(En);

		in.read ( (char*) pCoefs, sizeof(double)*nCount/sizeof(char) );

		if (eigst.isComplex())
			eigst.inAllCoef( (std::complex<double>*)pCoefs );
		else
			eigst.inAllCoef(pCoefs);

		delete [] pCoefs;
	}
	else
	{
		// handle first line in file
		in >> En;
		eigst.setEn(En);

		int nCount = eigst.isComplex() ? 2*eigst.dimension() : eigst.dimension();
		double dCoefi = 0.0;
		for (int i = 0; i < nCount; i++)
		{
			in >> dCoefi;
			eigst.inCoef(i, dCoefi);
		}
	}

	return in;
}

bool yosEigenState::freezeInstance()
{
	m_instance = this;
	return (true);

}

yosEigenState * yosEigenState::getInstance()
{

	return (m_instance);

}

void yosEigenState::initializeSearch()
{
	if (isSearchInitialized)
	{
		return;
	}
	/*
    Initialize NAG routine (set precisions)
	 */
#ifdef NAG
	const char * e04Text = "Function precision = 1.12d-12";
	const char * e04TolText = "Optimality Tolerance = 1.12d-16";

	e04urf_ ((char *) e04Text,29);
	e04urf_((char *)e04TolText, 32);
	string printLeveltext = 	"Print Level = 0";

	e04urf_( strdup(printLeveltext.c_str()), 15);
#else
#ifdef MKL
	// Initialize the Intel Search algorithm
	MKL_INT
	n = 2,
	m = 2,
	iter1 = 1000,
	iter2 = 100;
	double rs = 100.0;
	double x[2] = {CRandomizerGsl::getInstance()->fRand(0.0f,1.0f),CRandomizerGsl::getInstance()->fRand(0.0f,1.0f)};
	fvec = new double[m];

	for (int var = 0; var < m; ++var) {
		fvec[var] = 0.0;
	}
	fjac = new double[4];
	for (int index = 0; index < 4; ++index) {
		fjac[index] = 0.0;
	}
	double LW[2] = {0.0,0.0};
	double UP[2] = {1.0, 1.0};
	eps[0] =0.00001;
	eps[1] = 0.00001;
	eps[2] = 0.00001;
	eps[3] = 0.00001;
	eps[4] = 0.00001;
	eps[5] = 0.00001;
	int successflag  = 0;
	successflag = dtrnlspbc_init (&m_solverHandle, &n, &m, x, LW, UP, eps, &iter1, &iter2,&rs);
	/* set precisions for stop-criteria */
	if ( successflag!= TR_SUCCESS){
		glLogger.error("Could not initialize solver");
		return;
	}
	else
	{
		glLogger.info("Initialized Intel Solver");
	}
#endif
#endif
	isSearchInitialized = true;


}
#ifdef GSL
int yosEigenState::searchGsl(void)
{


	/*
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_vector *x;
  gsl_multimin_function_fdf my_func;
  double p[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 }; // do I need this??

  my_func.n = 2;
  my_func.f = my_f;

  my_func.df = my_df;
  my_func.fdf = my_fdf;

  my_func.params = 0;
    //parameter;
	 */
	// same here for f only minimax.
	const gsl_multimin_fminimizer_type *T;
	gsl_multimin_fminimizer *s;
	gsl_multimin_function my_func;
	gsl_vector *x;
	my_func.params = 0;
	my_func.n = 2;
	my_func.f = my_f;
	x = gsl_vector_alloc (2);
	int trialCount = 0;
	CRandomizerGsl::getInstance()->randomize();
	// main search loop starts here
	while (
			(zeroPositions.size() <  (unsigned int)myBasis->getNm() ) &&
			(trialCount < 500 ))
	{
		trialCount++;
		double m_x = CRandomizerGsl::getInstance()->dRand();
		double m_y = CRandomizerGsl::getInstance()->dRand();

		gsl_vector_set (x, 0, m_x);
		gsl_vector_set (x, 1, m_y);

		// choose here between the different minimax methods
		//T = gsl_multimin_fdfminimizer_conjugate_fr;
		T = gsl_multimin_fminimizer_nmsimplex2;
		//s = gsl_multimin_fdfminimizer_alloc (T, 2);
		s = gsl_multimin_fminimizer_alloc(T,2);
		// s = gsl_multimin_fminimizer_alloc
		// gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.0001, 1e-5);
		/*
		 *
		 * Only needed for simplex search
		 */
		gsl_vector *ss = gsl_vector_alloc (2);
		gsl_vector_set_all (ss, 0.05);
		gsl_multimin_fminimizer_set (s, &my_func, x, ss);

		//      gsl_multimin_fminimizer_set (s, &my_func, x, 0.0001, 1e-5);
		/*
		 * drive the minimizer , the loop is over several
		 * steps of minimum searching for one starting point
		 */
		int iter =  0;
		int status = 0;
		double size = 0.0;
		do
		{
			glLogger.info("In search loop ");
			glLogger.info("position in serach %5f,  %5f",gsl_vector_get(s->x,0),gsl_vector_get(s->x,1));
			iter++;
			// status = gsl_multimin_fdfminimizer_iterate (s);
			status = gsl_multimin_fminimizer_iterate (s);
			size = gsl_multimin_fminimizer_size (s);

			if (status)
			{
				glLogger.info("Status % i", status);
				break;
			}


			//status = gsl_multimin_test_gradient (s->gradient, 1e-3);
			double fctValTemp = s->fval;

			status = gsl_multimin_test_size (size, 1e-8);

			if ( (status == GSL_SUCCESS) || (fabs(fctValTemp) < 1e-08) )
			{
				std::complex<double> fctVal = sumofSPWavefunctions(0,gsl_vector_get (s->x, 0),gsl_vector_get (s->x, 1));

				glLogger.info("Minimum found at%5d %.5f %.5f %10.5f with value %f + i %f \n",
						iter,
						gsl_vector_get (s->x, 0),
						gsl_vector_get (s->x, 1),
						s->f,
						fctVal.real(),fctVal.imag());

				// Refactor later into single helper function
				double aToB = 1.0;
				CxPeriodicPosition tempPosition(gsl_vector_get (s->x, 0),gsl_vector_get (s->x, 1) , aToB );
				double testValue = sumofSPWavefunctions(0,tempPosition).real();
				tempPosition.setEpsilon(1e-05);
				//double fctValAbs = cabs(fctVal);
				//glLogger.info("Function value is %f", fctValAbs);


				double newPhase = getWindingNumber(tempPosition);
				if ( ((fabs(newPhase -1.0 ) < 1e-30 ) && (fctVal.real() < 1e-10 ) ) )
				{

					glLogger.info("Adding a position with x=(%f), y=(%f), xSize=(%f), ySize=(%f)",
							tempPosition.getXPosition(),tempPosition.getYPosition(),
							tempPosition.get_xCellSize(), tempPosition.get_yCellSize());
					zeroPositions.addPosition(tempPosition);
					CxVortex tempVortex(tempPosition,(int) rint(newPhase));
					vortexPositions.addPosition( tempVortex);

				}
				else
				{
					glLogger.info("Rejected zero! ");
				}


			}
			INFO("End of an iteration ");
		}
		while (status == GSL_CONTINUE && iter < 1000);
	}
	// end of main search loop
	return (-1);
}
#endif


int yosEigenState::searchIntel(void)
{	
#ifdef MKL
	// Initialization here for now
	MKL_INT
		n = 2,
		m = 2,
		iter1 = 1000,
		iter2 = 100;
		double rs = 100.0;
		double x[2] = {CRandomizerGsl::getInstance()->fRand(0.0f,1.0f),CRandomizerGsl::getInstance()->fRand(0.0f,1.0f)};
		fvec = new double[m];

		for (int var = 0; var < m; ++var) {
			fvec[var] = 0.0;
		}
		fjac = new double[4];
		for (int index = 0; index < 4; ++index) {
			fjac[index] = 0.0;
		}
		double LW[2] = {0.0,0.0};
		double UP[2] = {1.0, 1.0};

		eps[0] = 1e-08;
		eps[1] = 1e-08;
		eps[2] = 1e-08;
		eps[3] = 1e-08;
		eps[4] = 1e-08;
		eps[5] = 1e-08;


	glLogger.info("Entering yosEigenState::searchIntel!");
	/*int
	n = 2,
	m = 2;
	pjacobiObjFun = &my_f;
	double x[2] = {CRandomizer::getInstance()->fRand(0.0f,1.0f),CRandomizer::getInstance()->fRand(0.0f,1.0f)};


	int successflag  = 0;

			successflag = dtrnlspbc_init (&m_solverHandle, &n, &m, x, LW, UP, eps, &iter1, &iter2,&rs);
	*/
			/* set precisions for stop-criteria */

	/*
	if ( successflag!= TR_SUCCESS){
				glLogger.error("Could not initialize solver");
				return -1;
			}
			*/
	pjacobiObjFun = &my_f;
	// loop 2 cntrol
	bool outerLoopControl = true;
	//check handle
	if (!m_solverHandle)
	{
		throw CxErrors("No handle");
	}

	MKL_INT RCI_Request = 0;
	freezeInstance();
	int totalLoopCount = 0;
	while (outerLoopControl)
	{
		totalLoopCount++;
		bool loopControl = true;
		x[0] = CRandomizerGsl::getInstance()->fRand(0.0f,1.0f);
		x[1] = CRandomizerGsl::getInstance()->fRand(0.0f,1.0f);
		MKL_INT successflag  = 0;
				successflag = dtrnlspbc_init (&m_solverHandle, &n, &m, x, LW, UP, eps, &iter1, &iter2,&rs);
				/* set precisions for stop-criteria */
				if ( successflag!= TR_SUCCESS){
					glLogger.error("Could not initialize solver");
					return -1;
				}
	while (loopControl)
	{
		if (dtrnlspbc_solve (&m_solverHandle, fvec, fjac, & RCI_Request) != TR_SUCCESS)
		{
			throw CxErrors("Error while finding minimum ");
		}
		std::cerr << "solver returned " << RCI_Request<<std::endl;
		glLogger.info("solver returned (%d)", RCI_Request);
		if (RCI_Request == -1 || RCI_Request == -2 || RCI_Request == -3 ||
				RCI_Request == -4 || RCI_Request == -5 || RCI_Request == -6 || RCI_Request >10)
		{
			loopControl = 0;

			continue;
		}
		std::complex <double> fctVal;
		switch (RCI_Request) {
		case 1:
			fctVal = sumofSPWavefunctions(0,x[0],x[1]);
			//glLogger.info("Calculate function at (%f), (%f), = (%f) + i (%f)", x[0], x[1], fctVal.real(), fctVal.imag());
			std::cerr << fctVal <<std::endl;
			fvec[0]=fctVal.real();
			fvec[1]=fctVal.imag();
			break;
		case 2:
			//glLogger.info("calculate Jacobi Matrix");
			/* hier sollte die Jakobi MAtrix selber berechnet werden
			 * die Benutzung von djacobi tut nicht das was es soll
			 * Ergebnisse müssen in fjac landen
			 * fjac[0]=
			 * [1]=
			 * [2]=
			 * [3] =
			 */
			//my_df(x,0,fjac);
			if (djacobi(pjacobiObjFun,&n,&m,fjac, x,eps)!= TR_SUCCESS)
			{
				// Error handling here
				ERROR ("Error In calculation of Jacobi Matrix");
			}


			break;
		default:
			break;
		} // end of switch

	} //end of while
	MKL_INT
	stop_criterion = 0,
	no_iterations = 0,
	result = 0;
	double
	r1 = 0.0,
	r2 = 0.0;


	result = dtrnlsp_get(&m_solverHandle,&no_iterations,&stop_criterion,&r1, &r2);
	std::complex <double> fctVal = sumofSPWavefunctions(0,x[0],x[1]);
	if ((result == TR_SUCCESS) && std::abs(fctVal)<1e-08)
	{

		//glLogger.error("Success at (%f),(%f) with fct. Val (%f) + i (%f)", x[0],x[1],fctVal.real(), fctVal.imag());
		/*glLogger.error("After (%n) iterations with stop criterion (%n) and result (%n)",
				no_iterations,stop_criterion, result);
				*/
		double aToB = 1.0;
						CxPeriodicPosition tempPosition(x[0],x[1] , aToB );
						//double testValue = sumofSPWavefunctions(0,tempPosition).real();
						tempPosition.setEpsilon(1e-05);
						//double fctValAbs = cabs(fctVal);
						//glLogger.info("Function value is %f", fctValAbs);


						double newPhase = getWindingNumber(tempPosition);
						if ( ((fabs(newPhase -1.0 ) < 1e-10 ) && (fctVal.real() < 1e-9 ) ) )
						{

						/*	glLogger.info("Adding a position with x=(%f), y=(%f), xSize=(%f), ySize=(%f)",
									tempPosition.getXPosition(),tempPosition.getYPosition(),
									tempPosition.get_xCellSize(), tempPosition.get_yCellSize());
*/
							zeroPositions.addPosition(tempPosition);
							CxVortex tempVortex(tempPosition,(int) rint(newPhase));
							vortexPositions.addPosition( tempVortex);

						}
						else
						{
						//	glLogger.error("Rejected zero! ");
						}


	}
	else
	{
		//glLogger.error("could  not retrieve result");
	}
	if ((vortexPositions.size() == myBasis->getNm() ) || (totalLoopCount > 100) )
	{
		outerLoopControl = false;
	}
	dtrnlspbc_delete(&m_solverHandle);
	}
	return (-1);
#else
	glLogger.error("yosEigenState::searchIntel somehow called on the wrong plattform fix -D option");
	throw CxErrors("yosEigenState::searchIntel somehow called on the wrong plattform fix -D option");
#endif	
	// endif MKL defined
}
