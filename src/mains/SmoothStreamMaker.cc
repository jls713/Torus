//
// C++ code by Paul McMillan
//
// Using the Torus machinery to populate a very simple stream model
//

/**
\file SmoothStreamMaker.cc
\brief Make a stream, moderately simple prescription

Creates a stream. Parameters input are potential, progenitor actions & angles,
initial velocity dispersion & size of satellite,
time since first stripping (Myr) and number of stars wanted.

We assume that stars have been stripped at a constant rate since the
first one was stripped, and that the actions and initial angles are
independant of the time of stripping.

Following Eyre & Binney 2012 we take the spread in J to be

sig_JR = sig_v*(r_ap-r_peri)/pi
sig_Jz = 2*sig_v*Zmax/pi
sig_J  = sig_v*r_peri

And the initial spread in theta is (in terms of an effective radius)

sig_thR   = R_e * pi/(r_ap-r_peri)
sig_thz   = R_e * pi/(2*Zmax)
sig_thphi = R_e / r_peri

Which rather pleasingly means that we have
del^3 J del^3 th = del^3 x del^3 v
without any fine tuning.

Output is in code units and galactocentric cylindrical polar coordinates.

 */
#include <iostream>
#include <iomanip>
#include <fstream>



#include "Random.h"
#include "LogPot.h"
#include "falPot.h"
#include "Torus.h"
#include "PJM_utils.h"

using std::cout;
using std::cerr;




int main(int argc,char *argv[])
{
  Random3 R3(6), R3b(72);
  Gaussian Gau(&R3,&R3b);
  ifstream from;
  ofstream to;
  Torus T;
  bool interp;
  double stripwidth;

  if(argc<15) {
    cerr << "Creates a stream given potential, progenitor J & th, initial velocity dispersion (code units), progenitor radius, "
	 << "time since first stripping (Myr) and number of stars\n";
    cerr << "Input: Potential J_r J_z Jphi th_r th_z th_phi delv delx tmax nstars output_file interp fraction_of_period_strip_width (0. continuous)\n";
    cerr << "e.g. "<<argv[0]<<" pot/PJM11.Tpot 0.2 0.4 3.5 0.5 0.9 2 0.002 0.1 2000 400 tmp.st 1 0.1\n";
    exit(0);
  }

  interp = atoi(argv[13]);
  stripwidth = atof(argv[14]);

  int nstars;
  Actions J0, delJ, J;
  Angles A0, delA, A;
  Frequencies Om0, Om;
  double tmax, t, rp, ra, zmax, delv, delx;

  // Read in parameters
  my_open(to,argv[12]);
  Potential *Phi;
  if(std::string(argv[1])=="Log") Phi = new LogPotential(220.*Units::kms,0.8,0.,0.);
  else{
    my_open(from,argv[1]);
    Phi = new GalaxyPotential(from);
    from.close();
  }

  for(int i=0;i!=3;i++) {
    J0[i] = atof(argv[i+2]);
    A0[i] = atof(argv[i+5]);
  }
  T.AutoFit(J0,Phi);
  Om0 =  T.omegas();
  rp = T.minR();
  ra = T.maxR();
  zmax = T.maxz();

  delv = atof(argv[8]);
  delx = atof(argv[9]);

  delJ[0] = delv*(ra-rp)/Pi;
  delJ[1] = 2*delv*zmax/Pi;
  delJ[2] = delv*rp;

  for(int i=0;i!=3;i++) delA[i] = delx*delv/delJ[i];

  for(int i=0;i!=3;i++) cerr << delJ[i] << '\t';
  cerr << '\n';
  for(int i=0;i!=3;i++) cerr << delA[i] << '\t';
  cerr << '\n';

  tmax = atof(argv[10]);
  nstars = atoi(argv[11]);

  double dJ=1e99,tJ;
  for(int i=0;i<3;++i){ tJ = delJ[i]/J0[i]*0.01; if(tJ<dJ) dJ=tJ; }

  TorusInterpCell *TIC;
  if(interp)
    TIC = new TorusInterpCell(J0,Phi,delJ*5.,dJ,400);

  // Create stars
  for(int i=0;i!=nstars;i++) {
    for(int j=0;j!=3;j++) {
      do{
	     J[j] = J0[j] + delJ[j]*Gau();
      } while(J[j]<0 && j<2);
      A[j] = A0[j] + delA[j]*Gau();
    }
    t = tmax * R3();
    if(stripwidth){
      double tlastpericentre = A0[0]/Om0[0];
      int peri = (int)((t-tlastpericentre)*Om0[0]/TPi);
      t = (peri+stripwidth*Gau())*TPi/Om0[0]+tlastpericentre;
    }
    if(interp) Om=TIC->omegas(J);
    else{ T.AutoFit(J,Phi,dJ);
          Om = T.omegas(); }

    // Add the drift away in angle
    for(int j=0;j!=3;j++) {
      A[j] += (Om[j]-Om0[j])*t;
    }

    if(interp) to << TIC->Map3D(J,A);
    else to << T.Map3D(A);
    to << t << " " << J << " " << Om << '\n';

  }



}


