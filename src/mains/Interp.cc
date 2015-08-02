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

#include "Random.h"
using std::cerr;

int main(int argc,char *argv[])
{
  Random3 R3(6);
  Potential *Phi = new LogPotential(220.*Units::kms,0.9,0.,0.);
  Actions  J,J2,deltaJ;
  J[0] = 0.4; J[1] = 0.4; J[2] = 1.8;
  deltaJ[0] = 0.04; deltaJ[1] = 0.04; deltaJ[2] = 0.04;
  double   dJ=0.00003;int N=300;
  TorusInterpCell TIC(J,Phi,deltaJ,dJ,N);
  Angles A;
  for(int k=0;k<100;++k){
      for(int i=0;i<3;++i) A[i]=2.*Pi*R3();
      for(int i=0;i<3;++i) J2[i]=J[i]-deltaJ[i]/2.+R3()*deltaJ[i];
      TIC.test(J2,A);
  }
}


