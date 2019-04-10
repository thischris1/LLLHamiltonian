/*!
  \file CxVortexSystem.hpp
  \author Christian Mueller
*/

#ifndef CXVORTEXSYSTEM_HPP
#define CXVORTEXSYSTEM_HPP

#include <geometry/CxPosition.hpp>
#include <vector>
/*!
  \class CxVortexSystem
  \brief Describes a system with electrons and vortices.
  Analysis of vortex-electron relations is encapsulated here.
  

*/
using namespace std;
class CxVortexSystem
{
public:
  CxVortexSystem(double x, double y);
  
  CxVortexSystem();
  
  CxVortexSystem (std::vector<CxPosition>&, std::vector<CxPosition>&, double, double);
  bool addVortex(CxPosition&);

  bool addElectron(CxPosition&);

  bool addVortices (CxPosition *, int);
  
  bool addElectrons (CxPosition *, int);

  
protected:
  bool init(double, double);
  std::vector<unsigned int>  vorticesForElectron(int);

  std::vector<unsigned int> vorticesForElectron(CxPosition &);

  int electronForVortex(int);

  //  CxPosition electronForVortex(int);
  
  int electronForVortex(CxPosition &);

  
private:
  std::vector <CxPosition> *electrons;
  std::vector <CxPosition> *vortices;
  double xSize;
  double ySize;



};










#endif
