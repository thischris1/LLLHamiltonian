#include <gaussian.hpp>
#include <exception>
#include <utils/logger.hpp>
#ifdef MKL
#include <mkl.h>
#endif
/*!
  \file barrier_main.cpp
  \brief
  Main file for executable 'barrier'.
  Program to calculate eigenstates of a fqh-system with a barrier.
  Implemented barrier-types are: gaussian-shaped and delta-shaped (in k-space).
  The behavior is controlled by 'barrier.par'.   
 */

void usage();
void setup();

/*!
  \fn int main (int argc, char **argv)
  \brief main function of the program
  \param At commandline you can pass -v 
  \return 0 if successful
 */
int main (int argc, char **argv)
{
  std::set_terminate(__gnu_cxx::__verbose_terminate_handler);
  glLogger.error("Entering main");
  bool potentialOnly = false;
  std::cerr << argc<<"\n";
  switch (argc)
    {
    case (1): 
      {
	glLogger.setLogLevel (ERROR);
	glLogger.setFileName("./out/log");

		break;
		  }
	case (2):
		{
		if (strcmp(argv[1],"-p") == 0)
		{
			glLogger.error("Only potentials requested");
			glLogger.setLogLevel (ERROR);
					glLogger.setFileName("./out/log");
			potentialOnly = true;

		}

		break;
		}
	case (3):
		  {
		if (strcmp(argv[1],"-v") == 0)
		{
			glLogger.setLogLevel(argv[2]);
		}
		else
		{
			if (strcmp(argv[1],"-p") == 0)
			{
				potentialOnly = true;

			}
			//
			usage();

		}

		break;
		  }
	case (5):
		  {
		if (strcmp(argv[1],"-v") == 0)
		{
			glLogger.setLogLevel(argv[2]);
		}
		if (strcmp(argv[3],"-f") != 0)
		{
			if (strcmp(argv[3],"-p") == 0)
			{
				potentialOnly = true;

			}
			else {
				usage();
			}

		}

		else
		{
			glLogger.setFileName(argv[4]);
			glLogger.disableConsole();
		}
		break;
		  }


	default:
	{
		usage();

	}


	}
	ERROR("BEfore calling setup");

		// set up the threaded execution
		ERROR("Calling setup now");
		setup();
	glLogger.debug("After command line parsing");

	try {
		gaussian myBarrier(potentialOnly);

		// reading configuration
		std::string parameterFile("./gaussian.par");
		myBarrier.readParameters (parameterFile);

		// posibility to append some extra info in filenames
		//  myBarrier.createStandardFilenames ("");

		// calculate and save eigenvectors / -values
		myBarrier.calculate();
		std::cerr << "Back from calculate gaussian \n";
		// write Potential to file
		myBarrier.writePotential();

		//   write evalDensity.par
		myBarrier.writeEvalDensity();
	}
	catch (CxErrors &e)
	{
		glLogger.error("Error thrown from gaussian program (%s)", e.getMessage());
		e.print();
		return (-1);
	}
	return 0;
}



void usage() 
{
	std::cerr << "\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cerr << "\tThis is Christians  wonderfull Hamilton Diagonalizer\n";
	std::cerr <<"\n\t usage:  -v [1..4] Verbosity 1 smallest 4 highest\n";
	std::cerr <<"\t        -f name of logfile\n";
	std::cerr <<"\t        -p calculate only potentials\n";
	std::cerr << "\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	std::cerr <<" NOTE: When you only want the matirx in file matrix.dat, create an (empty) file WRITEMATRIX.dat in exec dir, works for sparse complex matrix type only"<<std::endl;
	exit (-1);


}      

void setup() {
#ifdef MKL_PARALLEL
	std::cerr << "Setting number of threads to 3 \n";
	mkl_set_num_threads(3);
#endif
}

