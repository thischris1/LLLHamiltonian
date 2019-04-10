/***************************************************************************
                          cxperiodicposition.h  -  description
                             -------------------
    begin                : Sat Oct 23 2004
    copyright            : (C) 2004 by Christian Mueller
    email                : cmueller@physnet.uni-hamburg.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CXPERIODICPOSITION_H
#define CXPERIODICPOSITION_H

#include <geometry/CxPosition.hpp>

/**A class which resembles positions in a periodic system
  *@author Christian Mueller
  */

class CxPeriodicPosition : public CxPosition  {
public: 
	CxPeriodicPosition();


	virtual ~CxPeriodicPosition();
  
  /** Full qualified ctor for a position in a periodic system */
   CxPeriodicPosition(double n_X, double n_Y, double n_XSize, double n_ySize);
  /** No descriptions */
   CxPeriodicPosition(double  n_X, double  n_Y,  double  aToB = 1.0);

   CxPeriodicPosition( const CxPosition &rhs, double  aspectRatio = 1.0);

/** Write property of double m_xCellSize. */
  void set_xCellSize( const double& _newVal);
  /** Read property of double m_xCellSize. */
  /**  */
 
  double get_xCellSize() const;

 /** Write property of double  m_yCellSize. */
  void set_yCellSize( const double & _newVal);

  /** Read property of double  m_YCellSize. */
   double  get_yCellSize() const;
  /** Constructor with default values */

  /** Addition operator */
  const CxPeriodicPosition  operator + (const CxPeriodicPosition &rhs);
  /** No descriptions */
  double vabs() const;
  float getAspectRatio(void) const { return ((float)m_xCellSize/ m_yCellSize);}

 
  /** Subtraktion  operator */
  const CxPeriodicPosition  operator -(const CxPeriodicPosition &rhs);
  /** Equal if cell and Position identical */
  bool operator ==(const CxPeriodicPosition & rhs);
  /** Copy constructor */
   CxPeriodicPosition ( const CxPeriodicPosition & rhs);
  /** assignement opertor**/
  CxPeriodicPosition & operator =( const CxPeriodicPosition & rhs);
  //! get nearest position from a given list 
   int getNearestPosition(std::vector<CxPeriodicPosition> vPositions); 
private: // Private attributes
  /**  */
  double m_xCellSize;
  /**  */
  double m_yCellSize;
	private:
	
  /** Maps the y Position into unit cell */
  double mapYpos(double yPos) const;
  /** Maps a x Position into the unit cell */
  double mapXPos(double xPos) const;
  /** No descriptions */
  CxPosition mapPosition(CxPosition & inPos);
	
};

#endif
