/*******************************************************************************
*                                                                              *
* READLIST.cc                                                                  *
*                                                                              *
* C++ code written by Walter Dehnen, 1994/95,                                  *
*                     Paul McMillan, 2007                                      *
* Oxford University, Department of Physics, Theoretical Physics.               *
* address: 1 Keble Road, Oxford OX1 3NP, United Kingdom                        *
* e-mail:  p.mcmillan1@physics.ox.ac.uk                                        *
*                                                                              *
*******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>

#include "Torus.h"
#include "PJMebf.h"
#include "PJM_cline.h"


using std::cout;
using std::cin;

int main(int argc,char *argv[])
{
  int      temp;
  Torus    *T, T0;

//---------------------------------------------------------------------------
  if(argc<2){
    cerr << "Input: TorusList_file (Torus_name)\n"
	 << "READLIST gives full information about one Torus in a list\n";
    exit(1);
  }
//---------------------------------------------------------------------------

  //TorusList Tlist(argv[1],0);
  
  std::string Tlist = string(argv[1]);
  std::vector<std::string> tori;
  if(Ebf_GetToriNames(Tlist,tori)) {
    cerr << "Not a torus list\n"; exit(1);
  }

  T = new Torus;
  if(argc >2) {
    for(int i = 2;i!=argc;i++) { 
      //temp = atoi(argv[i]);
      if(!(T->read_ebf(Tlist,string(argv[i])))) {
	cerr << "Torus "+string(argv[i])+"is not in the list.\n";
	  } else { 
	//Tlist.ExTorus(temp, *T);
	T->show(cout);
      }    
    } 
  } else {
    for(;;){
      string tname;
      cout << "There are " << tori.size() 
	   << " tori in the list, probably t1 to t" << tori.size() 
	   << ", which one do you want to know about?\n";
      cin >> tname;
      if ( !(T->read_ebf(Tlist,tname))) break;
      T->show(cout);
    }
  }
}
