//#include <iostream>
//#include "TROOT.h"

void wrapper_singles400(int m_entry_num){
  gROOT->ProcessLine(".L findSingles_new.C+");

  all(m_entry_num, 1500);
   
}
