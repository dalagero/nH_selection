#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TArrow.h"
#include "TFile.h"

#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TString.h"
#include "TMatrixD.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>
#include <TSQLiteServer.h>

void write_reactorInfo(){

	//Creating the new file
	  ifstream myfile ("/global/homes/d/dalagero/fromBeda/DybBerkFit-master/ReactorPowerCalculator/WeeklyAvg/WeeklyAvg_P17B_by_Beda_org.txt"); //File from Beda
	  ofstream newFile;
	  newFile.open("/global/homes/d/dalagero/fromBeda/DybBerkFit-master/ReactorPowerCalculator/WeeklyAvg/WeeklyAvg_P17B_by_Beda.txt"); //Altered file
	int week=0;
	int core=0;
	long start_s=0;
	long end_s=0;
	double power_frac=0;
	int zero=0;
	double frac_u235=0;
	double frac_u238=0;
	double frac_pu239=0;
	double frac_pu241=0;
	
	double new_power_frac=1.00;
	double new_frac_u235=0.564;
	double new_frac_u238=0.076;
	double new_frac_pu239=0.304;
	double new_frac_pu241=0.056;


	  if (myfile.is_open())
	  {
	    while (myfile >> week >> core >> start_s >> end_s >> power_frac >> zero >> frac_u235 >> frac_u238 >> frac_pu239 >> frac_pu241){ //Go line-by-line through the file
		//Write the code to write the new file
		if(core==1) newFile << week << "\t" << core << "\t" << start_s << "\t" << end_s << "\t" << new_power_frac << "\t" << zero << "\t" << frac_u235 << "\t" << frac_u238 << "\t" << frac_pu239 << "\t" << frac_pu241 << endl;
		else newFile << week << "\t" << core << "\t" << start_s << "\t" << end_s << "\t" << new_power_frac << "\t" << zero << "\t" << frac_u235 << "\t" << frac_u238 << "\t" << frac_pu239 << "\t" << frac_pu241 << endl;
	    }
	    myfile.close();
	  }
	  else cout << "Unable to open file" << endl; //End of the reactor information file

	newFile.close();
}


static int callback(void *data, int argc, char **argv, char **azColName){
  int i;
  fprintf(stderr, "%s: ", (const char*)data);
   
/*  for(i = 0; i<argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
  }
  	bkgd_spec_val = stod(argv[argc-1]);
    
  printf("\n");*/
  return 0;
}


void fill_reactorInfo(){
	char period_label[64];
	char reactor_label[64];
	sprintf(reactor_label,"Full Power Check 8/2/2022");

	double energy=0;
	int period[3]={6,8,7};
	double NuPerMeVPerSec[6];
	for(int icore=0; icore<6; icore++){
		NuPerMeVPerSec[icore]=0;
	}

	  sqlite3 *db_reactor;
	  char *zErrMsg_reactor = 0;
	  int rc_reactor;
	  char *sql_reactor;
	  const char* data = "Callback function called";

	  // Open database 
	  rc_reactor = sqlite3_open("parameters_THU.db", &db_reactor);

	  if( rc_reactor ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_reactor));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
//Need to loop over the 3 periods and go through all the files
	for(int iperiod=0; iperiod<3; iperiod++){
		sprintf(period_label,"%dad",period[iperiod]);

		//Reading the inputted file
		  ifstream myfile (Form("/global/homes/d/dalagero/fromBeda/DybBerkFit-master/ReactorPowerCalculator/isotope_spectra_by_Beda/reactor_%dAD_SNF_nonEq.txt", period[iperiod])); //File from toy

		  if (myfile.is_open())
		  {
		    while (myfile >> energy >> NuPerMeVPerSec[0] >> NuPerMeVPerSec[1] >> NuPerMeVPerSec[2] >> NuPerMeVPerSec[3] >> NuPerMeVPerSec[4] >> NuPerMeVPerSec[5]){ //Go line-by-line through the file
		    	for(int icore=0; icore<6; icore++){
			//CHANGE ALL OF THIS
			/*We want:
				CREATE TABLE IF NOT EXISTS "reactor_emitted_spectrum" (
				`Energy`	REAL NOT NULL,
				`Core`	INTEGER NOT NULL,
				`NuPerMeVPerSec`	REAL,
				`DataPeriod`	TEXT NOT NULL,
				`Source`	TEXT NOT NULL,
				PRIMARY KEY(`Energy`,`Core`,`Source`,`DataPeriod`)
			*/

			  // Create SQL statement
			  sql_reactor = Form("INSERT OR REPLACE INTO reactor_emitted_spectrum VALUES (%f, %d, %f, \"%s\", \"%s\");", energy, icore+1, NuPerMeVPerSec[icore], period_label, reactor_label);

			cout << sql_reactor << endl;

			   // Execute SQL statement
			   rc_reactor = sqlite3_exec(db_reactor, sql_reactor, callback, (void*)data, &zErrMsg_reactor);
			   }
		    }
		    myfile.close();
		  }
		  else cout << "Unable to open file" << endl; //End of the reactor information file

	}

	   if( rc_reactor != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_reactor);
	     sqlite3_free(zErrMsg_reactor);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_reactor);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;

}
