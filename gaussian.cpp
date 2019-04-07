/*!
  \file gaussian.cpp
  \brief Implementation of class Gaussian
  \author Christian Mueller
 */

#include <stdio.h>
#include <iostream>
#include <cassert>
#include <string>
#include <fstream>
#include <sstream>  
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>



#include <gaussian.hpp>
#include <myaccessories.hpp>
#include <utils/logger.hpp>
#include <yosEigenState.hpp>
#include <utils/CFileParser.h>
#include <utils/CLine.h>
#include <geometry/CxPeriodicPosition.hpp>
#include <geometry/CImpurityArray.h>
#include <LLLhamGaussianArray.h>
#include <geometry/CGaussianImpurityArray.h>
#include <utils/CRandomizer.h>
#include <utils/CRandomizerGsl.h>
#include <LLLRandomLandauMatrix.h> 
#include <vector>

/*
  save in fortran's "unformatted" format
  real valued eigenvectors
 */
//extern "C" void saveunformatted_ (char *filename, double *eigvals, double *coefs, int *eigsFound, int *dim);

/*
  save in fortran's "unformatted" format
  complex valued eigenvectors ( coefs = {(real,imag), ...} )
 */
//extern "C" void savecplxfort_ (char *filename, double *eigvals, double *coefs, int *eigsFound, int *dim);

/*
  Constants
 */


/*
  constructor
  int
    m_iNe ,
    m_iNm0,
    m_iSpinYes,
    m_iReqSz;

  bool
    m_bComplex,
    isRandomArray,
    isRandomMatrix; //! Use a Random Landau Matrix (RLM)
 	double m_impStrength;
 	double m_correlationLength;
	int m_randomIMax;

  int
    m_iMatType,
    m_iPotentialType,
    m_iInteractionType,
    m_iXKmax,
    m_iYKmax;


  double
    m_dA,
    m_dB,
    m_dBli,
    m_dEnergyOffset,
    m_dAlpha1,
    m_dAlpha2;

    CImpurityArray m_impArray;
 */
gaussian::gaussian (bool n_potentialOnly):
m_potentialOnly(n_potentialOnly),
m_iNe(0) ,
m_iNm ( 0),
m_iSpinYes (0) ,
m_iReqSz (0)
{

	m_strFileBas.clear();
	m_strFileVec.clear();
	m_strFileDens.clear();
	m_strFileV.clear();
	m_strFileLand.clear();
	isRandomArray = false;
	isRandomMatrix = false;
	m_spectrumFilePrefix = "spectrum.dat";


}

/*
  read the parameters from 'gaussian.par'
 */
int gaussian::readParameters(std::string & name_of_param_file)
{

	// Filenames

	CFileParser *myParser = CFileParser::getInstance();
	myParser->readFile(name_of_param_file);
	baseFileprefix = myParser->getNextLine().getString();
	stateFilePrefix = myParser->getNextLine().getString();
	densFilePrefix   = myParser->getNextLine().getString();
	landauFilePrefix = myParser->getNextLine().getString();
	m_strFileV  = myParser->getNextLine().getString();

	// related to whole system
	m_iNe = (myParser->getNextLine().getIntegers())[0];
	m_iNm  = (myParser->getNextLine().getIntegers())[0];


	// spin-stuff
	m_iSpinYes =(myParser->getNextLine().getIntegers())[0];;
	m_iReqSz= (myParser->getNextLine().getIntegers())[0];;

	// matrix type
	m_iMatType  = (myParser->getNextLine().getIntegers())[0];;

	// type of vectors file to generate
	m_iVectorFileFormat  = (myParser->getNextLine().getIntegers())[0];

	m_iEigsToFind  = (myParser->getNextLine().getIntegers())[0];
	m_dA  = (myParser->getNextLine().getFloats())[0];
	m_dB =  (myParser->getNextLine().getFloats())[0];
	m_dBli=(myParser->getNextLine().getFloats())[0];
	// gaussian parameters
	m_iPotentialType = (myParser->getNextLine().getIntegers())[0];
	m_iInteractionType= (myParser->getNextLine().getIntegers())[0];

	m_dEnergyOffset= (myParser->getNextLine().getFloats())[0];

	// flux of solenoids
	m_dAlpha1= (myParser->getNextLine().getFloats())[0];
	m_dAlpha2= (myParser->getNextLine().getFloats())[0];

	m_iXKmax = (myParser->getNextLine().getIntegers())[0];
	m_iYKmax = (myParser->getNextLine().getIntegers())[0];
	std::string impurityFileName = myParser->getNextLine().getString();



	if (impurityFileName == "random.dat")
	{
		isRandomArray = true;
		CRandomizerGsl::getInstance()->randomize();
		CGaussianImpurityArray tArray(impurityFileName);
		m_impArray = tArray;
		tArray.writePotentialToFile(std::string("./PotentialArray.dat"));
		tArray.writePositionsToFile(std::string("./CheckGaussianArray.dat"));
		if (m_potentialOnly == true)
		{
			exit(0);
		}

	}
	else
	{
		if( impurityFileName == "RANDOM_MATRIX")
		{
			isRandomMatrix = true;
			m_correlationLength = (myParser->getNextLine().getFloats())[0];
			m_impStrength = (myParser->getNextLine().getFloats())[0];
			m_randomIMax = (myParser->getNextLine().getIntegers())[0];
			glLogger.error("Random Matrix with correlationLength (%f), strength (%f), iMax (%d)",
					m_correlationLength, m_impStrength, m_randomIMax);
		}
		else
		{
			m_impArray.readFromFile(impurityFileName);
			glLogger.error("Impurities read from file");
			glLogger.error(impurityFileName.c_str());
			isRandomMatrix = false;
			isRandomArray = false;
		}
	}
	// energy-offset to keep eigenvalues negative

	std::string appendix("");
	createStandardFilenames(appendix);

	return 0;
}


/*!
  \brief No idea (yet). prepares the system for diagonalization

 */
int gaussian::diagonalize (TRESULT & listEigVectors, yosBasis & basis)
{
	glLogger.error("entering gaussian::diagonalize");
	glLogger.error ("   * Type of potential is %s", ( (m_iPotentialType == 0) ? "gaussian" : "delta" ) );
	glLogger.info ("   * Type of e-e interaction is %s", ( (m_iInteractionType == 0) ? "Coulomb" : "HARDCORE" ) );
	glLogger.error(" Number of impurities %d", m_impArray.getImpurities().size());
	if (m_impArray.getImpurities().size() == 1)
	{
		glLogger.error ("   * Strength = %f", m_impArray.getImpurities()[0].getStrength());
		glLogger.error ("   * xWidth = %f", m_impArray.getImpurities()[0].getSigmaX());
		glLogger.error ("   * xPos = %f", m_impArray.getImpurities()[0].getYPosition());
		glLogger.error ("   * yWidth = %f", m_impArray.getImpurities()[0].getSigmaY());
		glLogger.error ("   * yPos = %f", m_impArray.getImpurities()[0].getYPosition());
	} else {
		if (isRandomArray)
		{
			glLogger.error("Random array");
		}
		else if (isRandomMatrix)
		{
			glLogger.error("Random Matrix");
		}
	}
	glLogger.error ("   * Offset = %f", m_dEnergyOffset);
	glLogger.error ("   * Type of e-e-interaction is %s", (m_iInteractionType == 0) ? "Coulomb" : "hardcore" );
	glLogger.error( "   * geometry a =(%f), b= (%f)", m_dA, m_dB);
	glLogger.info ("   * Solenoid 1 = %f", m_dAlpha1);
	glLogger.info ("   * Solenoid 2 = %f\n", m_dAlpha2);
	glLogger.info(" With (%d) electrons", basis.getNe());

	LLLhamiltonian *hamBarrier = 0;
	if (isRandomMatrix)
	{
		hamBarrier = new LLLRandomLandauMatrix(&basis, m_iMatType, m_dBli,m_correlationLength,
				m_impStrength, m_randomIMax, m_dA, m_dB);
	}
	else
	{
		if (m_impArray.getImpurities().size() == 1)
		{
			//  int laufIndex=3;
			hamBarrier = new LLLhamGaussian(&basis,
					m_iMatType,
					m_dBli,
					m_dA,
					m_dB,
					m_impArray.getImpurities()[0],
					m_dEnergyOffset,
					m_iXKmax,
					m_iYKmax);

		}
		else
		{
			hamBarrier = new LLLhamGaussianArray(&basis,
					m_iMatType,
					m_dBli,
					m_dA,
					m_dB,
					m_dEnergyOffset,
					m_impArray.getImpurities());


		}
	}



	/*
    calculation of interaction matrixelements  
	 */

	switch(m_iInteractionType)
	  {
	  case (0):
	    	hamBarrier->computeCoulomb();
		break;
	  case (1):
	    hamBarrier->computeHardcore();
	    break;
	    /*  case (2):
	  hamBarrier->computeCoulomb(false);
	  break;*/
	  }
	/*
	if (m_iInteractionType == 0)
		hamBarrier->computeCoulomb();
	else
		hamBarrier->computeHardcore();
	*/
	// diagonalization
	// m_iEigsToFind will be modified by diagonalize if less then
	// required eigenvectors are found
	glLogger.debug("Back from computeCoulomb/hardcore");
	int nCount = 0;
	switch (m_iMatType)
	{
	case Operator::FULL_REAL:
		m_bComplex = false;
		break;
	case Operator::SPARSE_REAL:
		m_bComplex = false;
		nCount = basis.dimension();
		break;
	case Operator::FULL_COMPLEX:
		//case Operator::SPARSE_COMPLEX:
		m_bComplex = true;
		nCount = 2 * basis.dimension();
		break;
	case Operator::SPARSE_COMPLEX:
		m_bComplex = true;
		nCount = basis.dimension() * 2;
		break;
	}
	double
	*pdEigVec = 0,
	*pdEigVal = 0;
	glLogger.info("ncount =(%d), m_iEigsToFind =(%d), dimension = (%d)",
			nCount, m_iEigsToFind, basis.dimension() );
	int trialSize = nCount*m_iEigsToFind;

	try
	{

		pdEigVec = new double [trialSize];
	}
	catch(...)
	{
		glLogger.error("Could not allocate memory");
	}

	trialSize = basis.dimension();
	try
	{
		pdEigVal = new double [ trialSize ]; // this is too much, n_eigsTo>Find should be enough here
	}
	catch(...)
	{
		glLogger.error("Could not allocate memory");
	}



	// adjust # of Eigenvalues to find to dimension
	if (m_iEigsToFind > basis.dimension()-1)
	{
		glLogger.warning ("trying to find more eigenvalues than dimension adjusting to %d", basis.dimension()-1);
		m_iEigsToFind = basis.dimension()-1;
	}
	glLogger.error("Of for diagonalization");
	hamBarrier->diagonalize ( m_iEigsToFind, pdEigVec, pdEigVal );
	glLogger.error("back from diagonalization");
	std::ofstream spectrumFile(m_spectrumFilePrefix.c_str());
	spectrumFile.precision(10);
	// DUMP THE EIGENVECTORS AND VALUES HERE
	for (int eigVIndex = 0; eigVIndex < m_iEigsToFind; ++eigVIndex)
	{
		glLogger.info("Eigenvalue (%i)= (%f)", eigVIndex,pdEigVal[eigVIndex]);
		spectrumFile << pdEigVal[eigVIndex] << "\n";
		int offSet = nCount * eigVIndex;
		for (int eigVecIndex = 0;eigVecIndex  < basis.dimension(); eigVecIndex= eigVecIndex+1)
		{
			double reAl = pdEigVec[(2*eigVecIndex)+offSet];
			double imag = pdEigVec[(2*eigVecIndex)+1+offSet];
			std::complex<double> temp(reAl, imag);
			glLogger.debug("Coeffic. (%i) = (%f)+ i * (%f)", (2*eigVecIndex+offSet), reAl, imag);
		}
	}
	spectrumFile.close();


	// store re sults in list<persist_eigSt>
	persist_eigSt EigSt ( &basis, m_bComplex );

	for (int i = 0; i < m_iEigsToFind; i++) {
		if (m_bComplex)
		{
			EigSt.inAllCoef ( &( ((std::complex<double>*) pdEigVec)
					[i*basis.dimension()] ) );

			//EigSt.inAllCoef ( &( pdEigVec [i*basis.dimension()] ) );

			EigSt.setEn ( pdEigVal[i] );
			glLogger.info("Eigenvector %d with norm %f", i,EigSt.normsq());
			listEigVectors.push_back ( EigSt );
		}
		else
		{
			EigSt.inAllCoef(&pdEigVec[i*basis.dimension()]);
			EigSt.setEn(pdEigVal[i]);
			listEigVectors.push_back( EigSt );
		}
	}
	listEigVectors.sort();     // sort by energy
	// Keep
	ERROR("Write vectors now to file");
	saveVectors(listEigVectors);
	yosEigenState myState(&basis,
			listEigVectors.front());

	//	  listEigVectors.clear();

	std::string posFileName("./positions.dat");
	try{
		myState.readElectronPositions(posFileName.c_str());
	}
	catch (CxErrors &e)
	{
		ERROR("Position File not found or not valid");
		delete hamBarrier;
		hamBarrier = 0;
		delete [] pdEigVal;
		delete [] pdEigVec;
		listEigVectors.clear();
		return m_iEigsToFind;
	}
	//	 myState.writeReducedWfct("out/reducedWfct.txt", 0.05);
	myState.initializeSearch();
	myState.findZeros();

	std::vector<CxPeriodicPosition> zeroValues = myState.getShallowZeros();
	/*!
    \todo Make the fileName configurable  
    Write to file
	 */
	int zeroCount = zeroValues.size();
	std::vector<CxVortex>vortexPosition = myState.getShallowVortices();
	int vortexCount = (int) vortexPosition.size();
	vector <CxPosition> zeros;
	vector <CxPeriodicPosition> electrons = myState.getElectronPositions();
	std::stringstream convStream;
	convStream << "./results/zeros.dat_position_strength_";
	if (m_impArray.getImpurities().size() == 1)
	{
		convStream << m_impArray.getImpurities()[0].getStrength()<<"_sigma_"<<m_impArray.getImpurities()[0].getSigmaX()<<".dat";

	}
	else
	{
		convStream << "multiple";
	}

	std::string  fileName;
	convStream >> fileName;
	convStream.clear();
	convStream << "./results/vortex.dat_position_strength_";
	if (m_impArray.getImpurities().size() == 1)
	{
		convStream << m_impArray.getImpurities()[0].getStrength()<<"_sigma_"<<m_impArray.getImpurities()[0].getSigmaX()<<".dat";

	}
	else
	{
		convStream << "multiple";
	}

	std::string vortexFileName;
	convStream >> vortexFileName;
	glLogger.info( "Number of zeros is (%d)",zeroCount);

	if (zeroValues.size() != 0)
	{
		glLogger.error("About to write zeros to file (%s)", fileName.c_str());
		ofstream of(fileName.c_str());
		for (unsigned int index = 0; index < electrons.size(); index++)
		{
			of << " #Electron Number "<< index;
			of <<"x= "<< electrons[index].getXPosition();
			of << " y = " << electrons[index].getYPosition()<<"\n";
		}
		//     of << "# Strength " << impStrength <<"\n";
		of << "# Xposition \t       YPosition \t Next Electron X \t Y \t Abstand \n";
		if (m_iInteractionType == 0)
		{
			of << "# Coulomb interaction\n ";
		}
		else
		{
			of << "#Hardcore interaction \n";
		}

		for (int index = 0; index < zeroCount; index++)
		{
			of << zeroValues[index].getXPosition() << " \t  " << zeroValues[index].getYPosition() <<"\t";
			CxPeriodicPosition tempElectron =  electrons[zeroValues[index].getNearestPosition(electrons)];
			of << tempElectron.getXPosition() <<"\t\t "<< tempElectron.getYPosition() << "\t\t ";
			of << (tempElectron-zeroValues[index]).vabs()<<"\n";
		}
		of.close();
	}
	else
	{
		glLogger.info("Did not get zeros back :-(");
	} // End of if != zeroValues
	/*
    Write vortices next
	 */

	if (vortexPosition.size() != 0)
	{
		glLogger.error("About to write vortices to file (%s)", vortexFileName.c_str());
		ofstream vof(vortexFileName.c_str());
		for (unsigned int index = 0; index < electrons.size(); index++)
		{
			vof << " #Electron Number "<< index;
			vof << "x= "<< electrons[index].getXPosition();
			vof << " y = " << electrons[index].getYPosition()<<"\n";


		}
		//      vof << "# Strength " << impStrength <<"\n";
		vof << "# Xposition \t       YPosition \t Vorticity \n";
		if (m_iInteractionType == 0)
		{
			vof << "# Coulomb interaction\n ";

		}
		else
		{
			vof << "#Hardcore interaction \n";
		}

		for (int index = 0; index < vortexCount; index++)
		{
			vof << vortexPosition[index].getXPosition() << " \t\t  ";
			vof << vortexPosition[index].getYPosition() <<"\t\t";
			vof << vortexPosition[index].getWindingNumber() <<"\n";
		}
		vof.close();
	}
	else
	{
		glLogger.info("Did not get zeros back :-(");
	} // End of if != zeroValues

	/*
    DONE
	 */
	glLogger.error("Done, cleaning up now");

	delete hamBarrier;
	hamBarrier = 0;
	delete [] pdEigVal;
	delete [] pdEigVec;
	listEigVectors.clear();
	return m_iEigsToFind;
}

int gaussian::calculate()
{

	// result of diagonalization
	TRESULT *plistEigVectors;

	// yoshioka basis
	yosBasis basis (m_iNe, m_iNm, 0);

	// for barrier
	int nEigValuesFound  = 0;

	// list holds results
	plistEigVectors = new TRESULT();


	// generate the basis (all j included, spin polarized)
	glLogger.info("Building basis.");
	basis.genBasis_allJ_noSpin();

	glLogger.info("...done. Dimension = %d", basis.dimension());


	// do the calculation
	nEigValuesFound = diagonalize ( *plistEigVectors, basis );
	glLogger.info ("Found %d eigenvalues\n", nEigValuesFound);
	// basis
	FILE *pfileBasis = testfopenW ((char *) m_strFileBas.c_str() );
	assert ( pfileBasis != NULL );
	basis.writeBasisHeaderToFile ( pfileBasis );
	basis.writeBasisToFile ( pfileBasis );
	fclose ( pfileBasis );

	// vectors
	//    saveVectors ( *plistEigVectors );

	delete plistEigVectors;

	return 0;
}


void gaussian::createStandardFilenames (std::string & strExtraAppendix)
{


	glLogger.debug("Entering createStandardfilenames");

	std::stringstream tempStream;

	tempStream << m_iNe <<"-" << m_iNm;
	if (m_iInteractionType == 0)
	{
		tempStream << 'C';
	}
	else
	{
		tempStream << 'H';
	}
	if (m_impArray.getImpurities().size() == 1)
	{
		tempStream << "-strength-" << m_impArray.getImpurities()[0].getStrength();

		switch (m_iPotentialType)
		{
		case 0:
			tempStream <<"-wx-"<<m_impArray.getImpurities()[0].getSigmaX() <<"-wy-"<<m_impArray.getImpurities()[0].getSigmaY()<<".dat";
			break;
		case 1:
			tempStream <<"-wy-"<<m_impArray.getImpurities()[0].getSigmaY()<<".dat";
			break;
		}
	} else {
		glLogger.info("multiple branch in createStandardFileName");
		tempStream << "multiple_"<< m_impArray.getImpurities().size()<<".dat";
	}
	std::string tempString;
	tempStream >> tempString;

	m_strFileDens =  m_strFileDens+tempString;
	m_strFileLand =  m_strFileLand+tempString;
	tempStream.clear();
	tempStream << m_iNe<<"-"<<m_iNm;

	tempStream >> tempString;
	m_strFileBas= (baseFileprefix+tempString).append(".dat");
	glLogger.debug("basisfileprefix (%s)",baseFileprefix.c_str());
	tempString.clear();
	tempStream.clear();
	if  (m_impArray.getImpurities().size() == 1)
	{

		tempStream << stateFilePrefix << "strength"<<m_impArray.getImpurities()[0].getStrength()<<"-wx-"<<m_impArray.getImpurities()[0].getSigmaX()<<"-wy-"<<m_impArray.getImpurities()[0].getSigmaY()<<".dat";;

	}
	else
	{
		tempStream <<stateFilePrefix<<"multiple "<< m_impArray.getImpurities().size();
	}
	tempStream >> tempString;
	//  m_strFileVec.clear();
	m_strFileVec =  tempString;
	glLogger.error ("saving Basis in %s", m_strFileBas.c_str());
	glLogger.error ("saving Vectors in %s", m_strFileVec.c_str());

}

/*!
  \fn int gaussian::saveVectors (TRESULT & EigVectors)
  saving of eigenvectors
  filename from m_strFileVec
 */
int gaussian::saveVectors (TRESULT & EigVectors)
{
	// save to disk

	int i;
	// int iDim = EigVectors.begin()->dimension();

	//int iEigsFound = EigVectors.size();
	int nEigCount = 0;

	double *pdCoefs = 0;
	double *pdEigvals = 0;
	std::string fileNameBasis = m_strFileVec;
	int stateCount = 0;
	switch (m_iVectorFileFormat)
	{
	case 0: // save as ascii "c-code standard"
	{
		glLogger.info ("saving in c-standard ascii format\n");




		for (TRESULT::iterator run = EigVectors.begin(); run != EigVectors.end(); run++)
		{
			std::string fileName = fileNameBasis + boost::lexical_cast<std::string>(stateCount);
			std::ofstream fsVectors(fileName.c_str());
			fsVectors.precision(12); // precision
			glLogger.error ("eigenvalue #%d = %f", nEigCount++, run->getEn());
			glLogger.error (fileName.c_str());

			run->set_binary(false);
			fsVectors << (*run) << '\n';
			fsVectors.close();
			stateCount++;
		}
		break;
		/*

		glLogger.info ("saving in c-standard ascii format\n");
		std::ofstream fsVectors(m_strFileVec.c_str());
		fsVectors.precision(12); // precision
		fsVectors << "# Eigenvectors for gaussian impurities\n";

		for (TRESULT::iterator run = EigVectors.begin(); run != EigVectors.end(); run++)
		{
			glLogger.info ("eigenvalue #%d = %f", nEigCount++, run->getEn());

			run->set_binary(false);
			fsVectors << (*run) << '\n';
		}

		fsVectors.close();
	}
	break;
		 */
	}
	case 1: // save as binary data, c-style
	{
		glLogger.info ("saving in c-standard binary format");
		for (TRESULT::iterator run = EigVectors.begin(); run != EigVectors.end(); run++)
			{
		std::string fileName = fileNameBasis + boost::lexical_cast<std::string>(stateCount);
		std::ofstream bfsVectors(fileName.c_str(), std::ios::binary | std::ios::out);

		glLogger.error ("eigenvalue #%d = %f", nEigCount++, run->getEn());
		glLogger.error (fileName.c_str());

		run->set_binary(true);
		bfsVectors << (*run) << '\n';
		bfsVectors.close();
		stateCount++;

			}






	}
	break;

	case 2: // save as fortran-unformatted data
		glLogger.info ("saving in fortran's unformatted format");

		throw (CxErrors("fortran is not supported anymore",__FILE__,__LINE__));

		break;

	} // of switch

	return 0;
	}


	/*
  create evalDensity.par
	 */
	int gaussian::writeEvalDensity ()
	{
		std::ofstream fsEvDens("evalDensity_better.par");
		if (!fsEvDens) return 1;

		fsEvDens  << m_strFileBas << "# basisFilename (in)\n";
		fsEvDens  << m_strFileVec  << "0# vectorFilename (in)\n";
		fsEvDens <<  m_strFileDens  << "0# densityFilename (out)\n";
		fsEvDens << m_strFileLand  << "# diagonal of landau-matrix (out)\n";
		fsEvDens << "1\t\t# which vector to investigate\n";
		fsEvDens << (m_bComplex ? '1': '0' ) << "\t\t# wether to use complex eigenvectors\n";
		fsEvDens << m_iVectorFileFormat << "\t\t# type of file to read 0-> normal c-code 2->fortran unformatted\n";
		fsEvDens << m_dA/m_dB << "\t\t# aspect ratio a/b\n";
		fsEvDens << m_dAlpha1 << "\t\t# flux in solenoid1 affecting x-momentum in units of h/e\n";
		fsEvDens << m_dAlpha2 << "\t\t# flux in solenoid2 affecting y-momentum in units of h/e\n";
		fsEvDens << "GRID\n";
		fsEvDens << "100 \t\t# nrXStep\n";
		fsEvDens << "100 \t\t# nrYStep\n";
		fsEvDens << "3 \t\t# nr. of colums to evaluate\n";
		fsEvDens << "1 \t\t# what to put into (see notes): 3rd column (density)\n";
		fsEvDens << "8 \t\t# current density in x-direction\n";
		fsEvDens << "9 \t\t# current density in y-direction\n";

		fsEvDens << "\nNotes:\n- 'what to put into': 1st and 2nd column are the x and y coords.\n";
		fsEvDens << "what comes into 3rd, 4th, ... column can be directed by the following codes\n";
		fsEvDens << " 0	constant (everytime = 1.00)\n";
		fsEvDens << "		1	particle density\n";
		fsEvDens << "		2	spin density\n";
		fsEvDens << "		3	spin-down density\n";
		fsEvDens << "		4	Sz density\n";
		fsEvDens << "		5	spin density/particle density\n";
		fsEvDens << "		6	spin-down density/particle density\n";
		fsEvDens << "		7	Sz density/particle density\n";

		fsEvDens.close();

		return 0;
	}

	/*
   create file containig potential shape
	 */
	int gaussian::writePotential ()
	{
		glLogger.debug("Entering  gaussian::writePotential %s",m_strFileV.c_str());
		/*
    int 
    XSTEP = (int) (30.0/impXWidth), 
    YSTEP = (int) (30.0/impYWidth);

    double *V_x = new double [XSTEP];
    double *V_y = new double [YSTEP];

    int k0 = (int)( impXPos*m_iNm + 0.5 );

    double VDeltaY[100];
    double VDeltaX[100];



    std::ofstream fsV( m_strFileV.c_str() );

    if (!fsV) return 1;

    switch (m_iPotentialType)
    {
    case 0: // gaussfoermige Barriere
    double d;
    for ( int ix = 0; ix < XSTEP; ix++ ) 
    {
    V_x[ix] = 0.0;
    for (int k=-10; k<=10; k++) 
    {
    d = ( ((double)ix)/XSTEP - impXPos + k) / impStrength;
    V_x[ix] += exp( -d*d );
    }
    }

    for ( int iy = 0; iy < YSTEP; iy++ ) 
    {
    V_y[iy] = 0.0;

    }

    for (int iy = 0; iy < YSTEP; iy++) 
    for (int ix = 0; ix < XSTEP; ix++) 
    {
    double V = impStrength * V_x[ix] + V_y[iy]*V_x[ix];
    if (V >= 0.001 * impStrength)
    fsV << ((double)ix)/XSTEP << ' ' << ((double)iy)/YSTEP << ' ' << V << '\n';
	  }

      // "corners"
      fsV << 0.0 << ' ' << 0.0 << ' ' << impStrength*V_x[0] + V_y[0]*V_x[0] << '\n';
      fsV << 1.0 << ' ' << 0.0 << ' ' << impStrength*V_x[XSTEP-1] + V_y[0]*V_x[XSTEP-1] << '\n';
      fsV << 0.0 << ' ' << 1.0 << ' ' << impStrength*V_x[0] + V_y[YSTEP-1]*V_x[0] << '\n';
      fsV << 1.0 << ' ' << 1.0 << ' ' << impStrength*V_x[XSTEP-1] + V_y[YSTEP-1]*V_x[XSTEP-1] << '\n';
      break;


    case 1: // delta-potential

      for (int iy=0; iy < 100; iy++)
	{
	  VDeltaY[iy] = 0.0;
	  for (int i=0; i<10; i++) 
	    VDeltaY[iy] += cos (4.0*M_PI*m_iNm*i * iy/100.0) * exp ( -2*M_PI*m_iNm * m_dA/m_dB * i*i );
	}

      for (int ix=0; ix < 100; ix++)
	{
	  VDeltaX[ix] = 0.0;
	  for (int i=-20; i<= 20; i++)
	    {
	      VDeltaX[ix] += exp ( -2*M_PI*m_iNm*m_dA/m_dB * ( k0/(1.0*m_iNm) - ix/100.0 - i)*( k0/(1.0*m_iNm) - ix/100.0 - i) );
	    }
	}

      for (int iy = 0; iy < 100; iy++) 
	for (int ix = 0; ix < 100; ix++) 
	  fsV << ((double)ix)/100.0 << ' ' << ((double)iy)/100.0 << ' ' << VDeltaX[ix] * VDeltaY[iy] << '\n';

      fsV.close();
      break;

    }// of switch
  cerr << "Leave writePotential \n";
  delete [] V_x;
  delete [] V_y;
		 */
		return 0;
	}











