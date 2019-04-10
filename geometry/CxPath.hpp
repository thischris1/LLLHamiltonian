#ifndef CXPATH_HPP
#define CXPATH_HPP


//
/*!
	\file        CxPath.hpp

	\brief The file CxPath.hpp contains the declaration of the classCxPath 
	\author      Christian Mueller
	\date        27 Aug 03

        
  */ 
 
//
//    System Include Files
//


//
//    Application Include Files
//

#include <geometry/CxPosition.hpp>
//
//    Preprocessor Defines
//



//
//    Function Declarations
//


//
//    Class Declarations
//
/*!
	\class CxPath
	\brief The class CxPath encapsulates a Path in a 2D system. A path consists of 1-N points (which are CxPosition objects.) The Geometry of the system is specified as well
	\author Christian Mueller
	\date 27 Aug 03


 */ 
 
  class CxPath
{
public:
/*! 
	\fn  CxPath()
	\brief Constructor
*/
    CxPath();


    /*  \fn     : CxPath(const CxPath& rOrig);
     \brief Copy Constructor
     \param rOrig the value to be copied
    */
    CxPath(const CxPath& rOrig);


    //  \fn     : operator = 
    //! \brief Assignment Operator

    CxPath& operator = (const CxPath& rOrig);


    /*!
      \fn     : ~CxPath
     \brief Destructor
    */
    ~CxPath();


protected:
private:
    

  };




#endif // CXPATH_HH



