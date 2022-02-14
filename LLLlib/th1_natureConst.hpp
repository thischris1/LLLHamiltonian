//#######################################################
// File: th1_physConst.hpp    
// Description:
// Rev:
// Created: 04/10/2001
// Author: Michael Tews
// mail: mgtews@mtu.edu
//
// Copyright 
//           
//#######################################################
#ifndef TH1_PHYSCONST_HPP
#define TH1_PHYSCONST_HPP

// ANSI-Header
#include "math.h"


#define NC_ME  9.1093897e-31
#define NC_QE  1.60217733e-19
#define NC_H   6.6260755e-34
#define NC_HB  1.054588664e-34

// SI-Units
class th1_physConst
{
  public:
    th1_physConst()
    {
      me = 9.1093897e-31;
      ce = 1.60217733e-19;
      h  = 6.6260755e-34;
      hb = h/(2.0*M_PI);
    }
    double me;
    double ce;
    double h;
    double hb;
};

#endif 
