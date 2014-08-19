/*******************************************************************************
*                                                                              *
* noPT.h                                                                       *
*                                                                              *
* C++ code written by Walter Dehnen, 1994-97,                                  *
*                     Paul McMillan, 2007                                      *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/
//
// Header file for a no point transformation. 
// Input (r,TH,{phi},pr,pTH,{pphi}), Output (R,z,{phi},vR,vz,{v_phi})
//
// I.e. just converts from spherical polar coords to cylindrical.
//
//

#ifndef _NoTra_
#define _NoTra_ 1

#include "Maps.h"

////////////////////////////////////////////////////////////////////////////////
typedef Vector<double,6> PoiPar; // Vector.h
////////////////////////////////////////////////////////////////////////////////
class NoTra : public CanMap {

public:
    NoTra();
    NoTra(const NoTra& );
    ~NoTra() {}
    mutable bool derivs_ok;
    mutable double r,th,pr,pt,ir, ct,st,R,pR,z,pz;
    void    parameters	      (double *)		   const;
    int     NumberofParameters()                           const { return 0; }
    PSPD    Forward           (const PSPD&)                const;
    PSPD    Backward          (const PSPD&)                const;
    PSPT    Forward3D         (const PSPT&)                const;
    PSPT    Backward3D        (const PSPT&)                const;
    void    Derivatives       (double[4][4])               const;
};

inline void NoTra::parameters (double *tmp) const
{
  // Not applicable for this map. It has no parameters.
}

#endif
