/***************************************************************************//**
\file test_triaxial.cc
\brief Tests for the MultiPot class


*                                                                              *
*  test_triaxial.cc                                                            *
*                                                                              *
*  C++ code written by Paul McMillan, 2015,                                    *
* e-mail:  paul@astro.lu.se                                                    *
* github:  https://github.com/PaulMcMillan-Astro/Torus                         *
*                                                                              *
*******************************************************************************/
#include "ebf.hpp"
#include "Numerics.h"
#include <cmath>
#include "GeneratingFunction.h"
#include <fstream>

int main(int argc,char *argv[])
{
    GenPar G;
    G.MakeGeneric();
    std::ofstream out;
    out.open("tmp");
    G.write_log(out);
}
