/***************************************************************************//**
\file DF.h
\brief File containing classes that give distribution functions of standard forms

*                                                                              *
* DF.h                                                                         *
*                                                                              *
* C++ code written by Paul McMillan, 2010-                                     *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*       Classes which output the value of a distribution function              *
*                                                                              *
*                                                                              *
*             DFs are of the quasi-isothermal form introduced by               *
*                        Binney & McMillan 2011 (BM11)                         *
*******************************************************************************/
#include "Potential.h"
#include "Random.h"
#include "Units.h"
#include "Pi.h"

#ifndef DFbase_h
#define DFbase_h


//*****************************************************************************
// Base classes

/** \brief Base for classes that describe basic DF

 */
class DF {
public:
  virtual double df(Potential*,Actions) =0;        // give df  
  virtual int    setup(istream&)              =0;
  virtual int    setup_full(istream&)         =0;
  virtual void   Parameters(double*)          =0;
  virtual int    NumberofParameters()         =0;
};

// Base for classes that can quickly find the df for a single J, potential 
// repeatedly for multiple parameter sets.
class quickDF_lowmem {
 public:
  virtual double df()                                              =0;
  virtual double df(double*,bool*,int)                             =0;
  virtual int    setup(Potential*,Actions,double*,bool*,int) =0;
  virtual int    NumberofParameters()                              =0;
};

//*****************************************************************************
// quasi-isothermal disk, BM11 paramters
class quasi_iso_DF_JJB : public DF {
public:
  double R0, sig2R, sig2z, Rd, q;
  int setup(istream&);
  int setup_full(istream&);
  quasi_iso_DF_JJB(double,double,double,double,double);
  quasi_iso_DF_JJB() {;}
  ~quasi_iso_DF_JJB() {;}
  int    NumberofParameters() {return 5;}
  void   Parameters(double*);
  double df(Potential*,Actions);
  //void describe(ostream&);
};

//*****************************************************************************
// quasi-isothermal disk, Ting et al 2013 decription of parameters
class quasi_iso_DF_YST : public DF {
public:
  double p_sigR, p_sigz, Rd, h_sig;
  void setup(double,double,double,double);
  int setup(istream&);
  int setup_full(istream&);
  quasi_iso_DF_YST(double,double,double,double);
  quasi_iso_DF_YST() {;}
  ~quasi_iso_DF_YST() {;}
  int    NumberofParameters() {return 4;}
  void   Parameters(double*);
  double df(Potential*,Actions);
  //void describe(ostream&);
};



//*****************************************************************************
// Multiple quasi-isothermal disk, BM11 parameters
class multidisk_DF : public DF {
  double *params; //R0, sig2Ra, sig2za, Rda, Rsa, L0a, frac_a, etc;
  int ndiscs;
  double isumfracs;
public:
  int setup(istream&);
  int setup_full(istream&);
  multidisk_DF()  {ndiscs=0;}
  ~multidisk_DF();
  int    NumberofParameters() {return 2+6*ndiscs;}
  void   Parameters(double*);
  double df(Potential*,Actions);
};


class multidisk_quickDF : public quickDF_lowmem {
  double* params; //ndiscs, R0, sigRa0, sigza0, Rda, Rsa, L0a, fraca, etc 
  bool *change_params;
  int ndiscs;
  float R, *Sigrat, *prefac;
  Vector <float,2> epiJ; // kappa*J_R, nu*J_z
public:
  int setup(Potential*,Actions,double*,bool*,int);
  multidisk_quickDF() {ndiscs=0;}
  ~multidisk_quickDF();
  double df();
  double df(double*,bool*,int);
  int    NumberofParameters() {return 2+6*ndiscs;}
};




//*****************************************************************************
// quasi-isothermal disk, BM11 parameters - the meat
inline quasi_iso_DF_JJB::quasi_iso_DF_JJB(double r0, double sigr,double sigz, 
					  double RD, double Q) 
{
  setup(r0,sigr,sigz,RD,Q);
}

inline void quasi_iso_DF_JJB::setup(double r0, double sigr,double sigz, 
				    double RD, double Q) 
{
  R0    = r0;
  sig2R = sigr*sigr;
  sig2z = sigz*sigz;
  Rd    = RD;
  q     = Q;
}

inline int quasi_iso_DF_JJB::setup(istream &from) {
  double r0, sigr, sigz, RD, Q;
  from >> r0 >> sigr >> sigz >> RD;
  if(from.eof()) 
    cerr << "Not enough parameters in file for quasi_iso_DF_JJB\n";
  from >> Q;
  sigr *= Units::kms;
  sigz *= Units::kms;
  setup(r0,sigr,sigz,RD,Q);
}

inline int quasi_iso_DF_JJB::setup_full(istream &from) {  
  char type1;
  int n;
  from >> type1;
  if(type1!='j' && type1!='J') {
    cerr << "df not understood\n";
    return 0;
  }
  from >> n;
  if(n!=1) {
    cerr << "wrong number of disks in quasi_iso_DF_JJB\n";
    return 0;
  }
  setup(from);
}

inline void quasi_iso_DF_JJB::Parameters(double* output) {
  output[0] = R0;
  output[1] = sqrt(sig2R);
  output[2] = sqrt(sig2z);
  output[3] = Rd;
  output[4] = q;
}

inline double quasi_iso_DF_JJB::df(Potential* Phi,Actions J) {

  double R = Phi->RfromLc(fabs(J(2)));
  double Sig = exp(-R/Rd);
  Frequencies epi = Phi->KapNuOm(R);
  double Sigrat = exp(2.*q*(R0-R)/Rd);
  double sig2Rl = Sigrat*sig2R, sig2zl = Sigrat*sig2z;
  double L0 = 10.*Units::kms*Units::kpc;
  return 0.5*(1-tanh(J(2)/L0)) * epi(1)*epi(2)*Sig
    *exp(-epi(0)*J(0)/sig2Rl-epi(1)*J(1)/sig2zl)/
    (4*Pi*Pi*Pi*epi(0)*sig2Rl*sig2zl*Rd*Rd);
}

//*****************************************************************************
// quasi-isothermal disk, Yuen-Sen's style - the meat
inline quasi_iso_DF_YST::quasi_iso_DF_YST(double sr,double sz, 
					  double hR, double hs) 
{
  setup(sr,sz,hR,hs);
}

inline void quasi_iso_DF_YST::setup(double sr,double sz, 
				    double hR, double hs) 
{
  p_sigR = sr;
  p_sigz = sz;
  Rd     = hR;
  h_sig  = hs;
}

inline int quasi_iso_DF_YST::setup(istream &from) {
  double sr, sz, hR, hs;
  from >> sr >> sz >> hR;
  if(from.eof()) 
    cerr << "Not enough parameters in file for quasi_iso_DF_YST\n";
  from >> hs;
  sr *= Units::kms*Units::kms;
  sz *= Units::kms*Units::kms;
  setup(sr,sz,hR,hs);
}

inline int quasi_iso_DF_YST::setup_full(istream &from) {  
  char type1;
  int n;
  from >> type1;
  if(type1!='y' && type1!='Y') {
    cerr << "df not understood\n";
    return 0;
  }
  from >> n;
  if(n!=1) {
    cerr << "wrong number of disks in quasi_iso_DF_YST\n";
    return 0;
  }
  setup(from);
}

inline void quasi_iso_DF_YST::Parameters(double* output) {
 
  output[0] = p_sigR;
  output[1] = p_sigz;
  output[2] = Rd;
  output[3] = h_sig;
}

inline double quasi_iso_DF_YST::df(Potential* Phi,Actions J) {

  double R = Phi->RfromLc(fabs(J(2)));
  double Sig = exp(-R/Rd);
  Frequencies epi = Phi->KapNuOm(R);
  double Sigrat = exp(-R/h_sig);
  double sig2Rl = Sigrat*p_sigR, sig2zl = Sigrat*p_sigz;
  return epi(1)*epi(2)*Sig * exp(-epi(0)*J(0)/sig2Rl-epi(1)*J(1)/sig2zl) /
    (4*Pi*Pi*Pi*epi(0)*sig2Rl*sig2zl*Rd*Rd);
}



//*****************************************************************************
// Many quasi-isothermal disks, BM11 parameters - the meat

inline int multidisk_DF::setup(istream &from) {
  if(ndiscs) delete[] params;
  from >> ndiscs;
  params = new double[1+6*ndiscs]; 
  // r0,ndiscs_parameter_sets(isig2R0,isig2z0,iRd,iRs,L0,frac)

  //double r0, sigr,sigz,RD,Q,L0,f;//,sigrb,sigzb,RDb,Qb,fb;
  for(int i=0;i!=1+6*ndiscs;i++) {
    if(from.eof()) { 
      cerr << "Not enough parameters in file for multidisk_DF\n";
      return 0;
    }
    from >> params[i];//
  }
  isumfracs=0.;
  for(int i=0;i!=ndiscs;i++) {
    params[1+6*i] = 
      1./(params[1+6*i]*params[1+6*i]*Units::kms*Units::kms);
    params[2+6*i] = 
      1./(params[2+6*i]*params[2+6*i]*Units::kms*Units::kms);
    params[3+6*i] = 1./params[3+6*i];
    params[4+6*i] = 1./params[4+6*i];
    params[5+6*i] = params[5+6*i]*Units::kms;
    isumfracs += params[6+6*i];
  }
  isumfracs = 1./isumfracs;
  return 1;
}

inline int multidisk_DF::setup_full(istream &from) {  
  char type1;
  int n;
  from >> type1;
  
  if(type1 !='m' && type1!='M') {
    cerr <<type1<< " df not understood\n";
    return 0;
  }
  return setup(from);
}

inline void multidisk_DF::Parameters(double* output) {
  output[0] = ndiscs;
  output[1] = params[0];
  for(int i=0;i!=ndiscs;i++) {
    output[2+6*i] = 1./sqrt(params[1+6*i]);
    output[3+6*i] = 1./sqrt(params[2+6*i]);
    output[4+6*i] = 1./params[3+6*i];
    output[5+6*i] = 1./params[4+6*i];
    output[6+6*i] = params[5+6*i];
    output[7+6*i] = params[6+6*i];
  }
}

inline double multidisk_DF::df(Potential *Phi, Actions J) 
{
  if(J(0)<0. || J(1) <0.) return 0.;
  double R = Phi->RfromLc(fabs(J(2)));
  double R0 = params[0];
  Frequencies epi = Phi->KapNuOm(R);
  double prefac = 0.5*epi(1)*epi(2)/(4.*Pi*Pi*Pi*epi(0));
  double out = 0., oldout=0.;
  for(int i=0;i!=ndiscs;i++) {
    double is2r = params[1+6*i], is2z = params[2+6*i], 
      iRd = params[3+6*i], iRs = params[4+6*i], 
      L0 = params[5+6*i], frac = params[6+6*i];
    double Sigrat = exp(2.*(R-R0)*iRs);
    double cut = (1-tanh(J(2)/L0));
    double isig2Rl = Sigrat*is2r, isig2zl = Sigrat*is2z;
    out += frac*prefac*cut * exp(-epi(0)*J(0)*isig2Rl-epi(1)*J(1)*isig2zl-R*iRd)
      *isig2Rl*isig2zl*iRd*iRd;
    oldout = out;
  }
  return out*isumfracs;
}

inline  multidisk_DF::~multidisk_DF() {
  if(ndiscs) delete[] params;
}


//------------------------------------------------------------------------------
// multidisk_quickDF                   -----------------------------------------
//
//

inline int multidisk_quickDF::setup(Potential* Phi, 
				    Actions j, double* p_in, 
				    bool* c_p_in, int np)
{
  if(np<1) { cerr << "really need positive number of parameters np\n"; exit(0);}
  ndiscs = int(p_in[0]);
  if(np!=NumberofParameters()) {
    cerr << "multidisk_quickDF requires " 
	 << NumberofParameters() 
	 << "parameters (given the claimed number of discs) as input\n";
    exit(0);
  }
  params = p_in;
  change_params = c_p_in;
  for(int i=0;i!=2;i++) epiJ[i] = j[i];
  R = Phi->RfromLc(fabs(j(2)));
  Frequencies epi_tmp = Phi->KapNuOm(R);
  for(int i=0;i!=2;i++) epiJ[i] *= epi_tmp[i];
  double R0  = params[1];

  Sigrat = new float[ndiscs];
  prefac = new float[ndiscs];

  for(int i=0;i!=ndiscs;i++) {
    float Rd  = params[4+6*i];
    float Rs  = params[5+6*i];
    float L0  = params[6+6*i];
    Sigrat[i] = exp(2.*(R0-R)/Rs);
    prefac[i] =  0.5*(1-tanh(j(2)/L0))*epi_tmp(1)*epi_tmp(2)/
      (4.*Pi*Pi*Pi*epi_tmp(0));
  }
}

inline double multidisk_quickDF::df(double* p_in, bool* c_p_in, 
					       int np)
{ 
  if(np!=NumberofParameters()) {
    cerr << "multidisk_quickDF requires " << NumberofParameters() 
	 << "parameters as input now that it's set up with "
	 << ndiscs <<" discs\n";
    exit(0);
  }
  params = p_in;
  change_params = c_p_in;
  
  return df();
}

inline double multidisk_quickDF::df()
{ 
  if(epiJ(0)<0. || epiJ(1)<0.) return 0.;
  double prefac_full, isig2Rl, isig2zl, weight = 0., sum_fracs=0.;
  double R0 = params[1];
  for(int i=0;i!=ndiscs;i++) {
    if(change_params[6+6*i]) { 
      cerr << "Can't change L0, as set up\n";
      exit(0);
    }
    double iRd = 1./params[4+6*i]; // N.B. save time by storing this?
    double Rs  = params[5+6*i];  
    if(change_params[5+6*i] || change_params[1]) 
      Sigrat[i] = exp(2.*(R0-R)/Rs);

    isig2Rl = 1./(Sigrat[i]*params[2+6*i]*params[2+6*i]);
    isig2zl = 1./(Sigrat[i]*params[3+6*i]*params[3+6*i]);

    prefac_full = prefac[i]*exp(-R*iRd)*iRd*iRd;
    weight += params[7+6*i]*prefac_full*
      exp(-epiJ(0)*isig2Rl-epiJ(1)*isig2zl)*isig2Rl*isig2zl;
    sum_fracs += params[7+6*i];  // N.B. save time (maybe) by storing 1/this
  }
  return weight/sum_fracs;
}

inline multidisk_quickDF::~multidisk_quickDF() {
  if(ndiscs) {
    delete[] Sigrat;
    delete[] prefac;
  }
}







//------------------------------------------------------------------------------
// function which sets up a DF from file
//

inline DF* set_DF(istream &from) { // return success?
  // check first two characters - tells me what it is.
  char type1;
  int n;
  DF *tmpdf;
  from >> type1;
  if(type1=='m' || type1=='M') {
      tmpdf = new multidisk_DF;
      tmpdf->setup(from);
      return tmpdf;
    return NULL;
  } else if(type1=='j' || type1=='J') {
    from >> n;
    if(n==1) {
      tmpdf = new quasi_iso_DF_JJB;
      tmpdf->setup(from);
      return tmpdf;
    } else {
      cerr << "asking for more disks ("<<n<<") than I can provide yet\n";
      return NULL;
    }
  } else if(type1=='y' || type1=='Y') {
    from >> n;
    if(n==1) {
      tmpdf = new quasi_iso_DF_YST;
      tmpdf->setup(from);
      return tmpdf;
    } else {
      cerr << "asking for more disks ("<<n<<") than I can provide yet\n";
      return NULL;
    }
  } else 
    cerr << "df not understood\n";
  return NULL;
}







#endif
