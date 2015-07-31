/***************************************************************************//**
\file Interp.cc
\brief Gives an example of using the interpolated torus cell

*                                                                              *
* Interp.cc                                                                 *
*                                                                              *
* Jason Sanders (2015)                                                       *
*******************************************************************************/


#include <iostream>
#include <iomanip>
#include <fstream>

#include "PJMebf.h"
#include "Toy_Isochrone.h"
#include "Torus.h"

#include "falPot.h"
#include "LogPot.h"
#include "MiyamotoNagaiPot.h"

#include "PJM_cline.h"

using std::cerr;

int main(int argc,char *argv[])
{
  Potential *Phi = new LogPotential(220.*Units::kms,0.8,0.,0.);
  Actions  J,deltaJ;
  J[0] = 0.4; J[1] = 0.4; J[2] = 1.8;
  deltaJ[0] = 0.1; deltaJ[1] = 0.1; deltaJ[2] = 0.1;
  double   dJ=0.0003;int N=400;
  TorusInterpCell TIC(J,Phi,deltaJ,dJ,N);
  Angles A;A[0]=0.4;A[1]=1.4;A[2]=5.2;
  J=J+deltaJ*.2;
  TIC.test(J,A);
}


