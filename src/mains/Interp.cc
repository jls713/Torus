/***************************************************************************//**
\file AUTOTORUS.cc
\brief Fits a single torus from actions given at the command line

*                                                                              *
* AUTOTORUS.cc                                                                 *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2007                                      *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
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

class TorusInterpCell{
private:
  Actions J0,J,deltaJ,deltaJh;
  Potential *Phi;
  std::vector<Torus*> Tori;
  Torus* TorusToUse;
  double dJ,deltaprod;
  double flag;
  const int N;
  GenPar GP,GP1,GP2,GP3;
  AngPar AP;
public:
  TorusInterpCell(Actions J0, Potential* Phi, Actions deltaJ, double dJ,int N):J0(J0),Phi(Phi),deltaJ(deltaJ),dJ(dJ),N(N){
    deltaJh=.5*deltaJ;
    deltaprod=deltaJ[0]*deltaJ[1]*deltaJ[2];deltaprod=1./deltaprod;
    TorusToUse=new Torus;
    flag = TorusToUse->AutoFit(J0,Phi,dJ,N,300,15,5,24,200,24,0);
    IsoPar IM = TorusToUse->TP();
    std::vector<int> plus_minus(3,0);
    for(int i=0;i<8;++i){
      Tori.push_back(new Torus);
      plus_minus[0]=(i%2)*-2+1;
      plus_minus[1]=((i/2)%2)*-2+1;
      plus_minus[2]=((i/4)%2)*-2+1;
      for(int k=0;k<3;++k)J[k]=J0[k]+plus_minus[k]*deltaJh[k];
      Tori[i]->AutoFit(J,Phi,dJ,N,300,15,5,24,200,24,0,IM);
    }
    TorusToUse->SetTP(IM);
  }
  PSPT Map3D(Actions J,Angles A){
    Actions ddJ = J-J0+deltaJh,mddJ = deltaJ-ddJ;
    TorusToUse->SetActions(J);
    GP =  ddJ[2]*(ddJ[1]*(ddJ[0]*Tori[0]->SN()+mddJ[0]*Tori[1]->SN())
                +mddJ[1]*(ddJ[0]*Tori[2]->SN()+mddJ[0]*Tori[3]->SN()))
        +mddJ[2]*(ddJ[1]*(ddJ[0]*Tori[4]->SN()+mddJ[0]*Tori[5]->SN())
                +mddJ[1]*(ddJ[0]*Tori[6]->SN()+mddJ[0]*Tori[7]->SN()));
    GP *= deltaprod;
    GP1 =  ddJ[2]*(ddJ[1]*(ddJ[0]*Tori[0]->AP().dSdJ1()
                         +mddJ[0]*Tori[1]->AP().dSdJ1())
                +mddJ[1]*(ddJ[0]*Tori[2]->AP().dSdJ1()
                         +mddJ[0]*Tori[3]->AP().dSdJ1()))
        +mddJ[2]*(ddJ[1]*(ddJ[0]*Tori[4]->AP().dSdJ1()
                         +mddJ[0]*Tori[5]->AP().dSdJ1())
                +mddJ[1]*(ddJ[0]*Tori[6]->AP().dSdJ1()
                         +mddJ[0]*Tori[7]->AP().dSdJ1()));
    GP1 *= deltaprod;
    GP2 =  ddJ[2]*(ddJ[1]*(ddJ[0]*Tori[0]->AP().dSdJ2()
                         +mddJ[0]*Tori[1]->AP().dSdJ2())
                +mddJ[1]*(ddJ[0]*Tori[2]->AP().dSdJ2()
                         +mddJ[0]*Tori[3]->AP().dSdJ2()))
        +mddJ[2]*(ddJ[1]*(ddJ[0]*Tori[4]->AP().dSdJ2()
                         +mddJ[0]*Tori[5]->AP().dSdJ2())
                +mddJ[1]*(ddJ[0]*Tori[6]->AP().dSdJ2()
                         +mddJ[0]*Tori[7]->AP().dSdJ2()));
    GP2 *= deltaprod;
    GP3 =  ddJ[2]*(ddJ[1]*(ddJ[0]*Tori[0]->AP().dSdJ3()
                         +mddJ[0]*Tori[1]->AP().dSdJ3())
                +mddJ[1]*(ddJ[0]*Tori[2]->AP().dSdJ3()
                         +mddJ[0]*Tori[3]->AP().dSdJ3()))
        +mddJ[2]*(ddJ[1]*(ddJ[0]*Tori[4]->AP().dSdJ3()
                         +mddJ[0]*Tori[5]->AP().dSdJ3())
                +mddJ[1]*(ddJ[0]*Tori[6]->AP().dSdJ3()
                         +mddJ[0]*Tori[7]->AP().dSdJ3()));
    GP3 *= deltaprod;
    TorusToUse->SetSN(GP);
    TorusToUse->SetAP(AngPar(GP1,GP2,GP3));
    return TorusToUse->Map3D(A);
  }
  void test(Actions J,Angles A){
    std::cout<<Map3D(J,A)<<std::endl;
    Torus *T=new Torus;
    T->AutoFit(J,Phi,dJ,N,300,15,5,24,200,24,0);
    std::cout<<T->Map3D(A)<<std::endl;
    return;
  }
};

int main(int argc,char *argv[])
{
  Potential *Phi = new LogPotential(220.*Units::kms,0.8,0.,0.);
  Actions  J,deltaJ;
  J[0] = 0.4; J[1] = 0.4; J[2] = 1.8;
  deltaJ[0] = 0.1; deltaJ[1] = 0.1; deltaJ[2] = 0.1;
  double   dJ=0.003;int N=200;
  TorusInterpCell TIC(J,Phi,deltaJ,dJ,N);
  Angles A;A[0]=0.4;A[1]=1.4;A[2]=5.2;
  J=J+deltaJ*.3;
  TIC.test(J,A);
}


