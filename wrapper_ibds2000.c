//#include <iostream>
//#include "TROOT.h"

void wrapper_ibds2000(int m_entry_num){
  gROOT->ProcessLine(".L findIBDs_2000.C+");

  all(m_entry_num,1500);
   
}
