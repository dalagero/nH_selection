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

const int maxAD=8;
const int numStages=3;
int EH[maxAD]={1,1,2,2,3,3,3,3};
int AD[maxAD]={1,2,1,2,1,2,3,4};
const int stage[3]={6,8,7};
const int numBins=34;
const int numColumns=10;

char toyConfig[128];
char exConfig[128];
char fitConfig[128];

char parameters_db[64];
char background_db[64];
char prompt_db[64];

char toy_label[64];
char bkgd_counts_label[64];
char bkgd_spectra_label[64];
char acc_spectra_label[64];
char num_coincs_label[64];
char det_resp_label[128];
char muon_label[64];

char livetime_table[64];
char bkgd_table[64];
char acc_spec_table[64];
char li9_spec_table[64];
char fastn_spec_table[64];
char amc_spec_table[64];
char radn_spec_table[64];
char num_coincs_table[64];

char toy_file[64];
char toy_hist[64];
char acc_spec_file[128];
char acc_spec_hist[64];
char li9_spec_file[128];
char li9_spec_hist[64];
char fastn_spec_file[128];
char fastn_spec_hist[64];
char amc_spec_file[128];
char amc_spec_hist[64];
char radn_spec_file[128];
char radn_spec_hist[64];

//Numbers from toy config file
double DAQ[maxAD][numStages];
double totDAQ[maxAD];
double muonEff[maxAD];
double multEff[maxAD];
double accRate[maxAD];
double accRate_uncert[maxAD];
double li9Rate[maxAD];
double li9Rate_uncert[maxAD];
double fastnRate[maxAD];
double fastnRate_uncert[maxAD];
double amcRate[maxAD];
double amcRate_uncert[maxAD];
double alphanRate[maxAD];
double alphanRate_uncert[maxAD];


//Other
const double muon_rate[8] = {200.32, 200.32, 150.08, 149.80, 15.748, 15.748, 15.748, 15.747}; //rate in Hz


void init_names(){
	sprintf(toy_label, "toys_100");

	sprintf(toyConfig,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/input/Theta13-inputs_P17B_inclusive_");
	sprintf(exConfig,"fit_config_nH_resid_flash_example.json");
	sprintf(fitConfig,"fit_config_nH_resid_flash_%s.json",toy_label);
	
	sprintf(parameters_db,"parameters.db");
	sprintf(background_db,"nH_backgrounds.db");
	sprintf(prompt_db,"adtime_dt_eff_study.db");

	sprintf(bkgd_counts_label,"%s", toy_label);
	sprintf(bkgd_spectra_label,"nH backgrounds, %s", toy_label);
	sprintf(acc_spectra_label,"accidentals, %s", toy_label);
	sprintf(num_coincs_label,"%s", toy_label);
//	sprintf(det_resp_label,"Response based on Jinjing's matrices by Beda, nH modified 2 binning");
	sprintf(det_resp_label,"Matrix for fitter from Beda, 8/1/22, nH modified 2 binning");
	sprintf(muon_label,"%s", toy_label);

	sprintf(livetime_table,"muon_rates");
	sprintf(bkgd_table,"bg_counts");
	sprintf(acc_spec_table,"accidentals_spectrum");
	sprintf(li9_spec_table,"li9_spectrum");
	sprintf(fastn_spec_table,"fast_neutron_spectrum");
	sprintf(amc_spec_table,"amc_spectrum");
	sprintf(radn_spec_table,"rad_n_spectrum");
	sprintf(num_coincs_table,"num_coincidences");

	sprintf(toy_file,"./toy/%s.root", toy_label);
	sprintf(acc_spec_file,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/input/accidental_eprompt_shapes_8ad.root");
	sprintf(li9_spec_file,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/li9_spectrum/8he9li_nominal_spectrum.root");
	sprintf(fastn_spec_file,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/fn_spectrum/P15A_fn_spectrum.root");
	sprintf(amc_spec_file,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/amc_spectrum/amc_spectrum.root");
	sprintf(radn_spec_file,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/muon_decay_spectrum/MuonDecaySpec.root");

	sprintf(toy_hist,"h_nominal_stage");
	sprintf(acc_spec_hist,"h_accidental_eprompt_fine_inclusive");
	sprintf(li9_spec_hist,"h_nominal");
	sprintf(fastn_spec_hist,"h_2AD_fn_fine");
	sprintf(amc_spec_hist,"h_toy");
	sprintf(radn_spec_hist,"MdSpec_EH1");

}

void parse_config(){ //FIXME: Going to only look at all 3 files for DAQ time, and stay on 8ad period for the rest
	init_names();
	
	//Reading the config file
	  ifstream myfile (Form("%s%dad.txt",toyConfig,stage[1])); //Get all the values from the 8 AD period
	  string line;
	  int iline=0;
	  string tab="\t";
	  string temp[numColumns];
	  double value[numColumns];

	  if (myfile.is_open())
	  {
	    while (getline(myfile, line, '\n')){
	    	iline+=1;
	    	if(iline<40) continue;
	    	
	    	//Breaking down the line
	    	for(int i=0; i< numColumns; i++){
	    		temp[i] = line.substr(0,line.find(tab));
	    		line.erase(0,line.find(tab)+tab.length());
	    		value[i] = stod(temp[i]);
	    	}
	    	
	    	if(value[1]==2){//DAQ time
	    		for(int iad=0; iad<maxAD; iad++){
	    			DAQ[iad][1] = value[iad+2];
				totDAQ[iad] = DAQ[iad][1];
	    		}
	    	}
	    	if(value[1]==3){//Muon Veto Efficiency
	    		for(int iad=0; iad<maxAD; iad++){
	    			muonEff[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==4){//Multiplicity Efficiency
	    		for(int iad=0; iad<maxAD; iad++){
	    			multEff[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==11){//Accidentals Rate
	    		for(int iad=0; iad<maxAD; iad++){
	    			accRate[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==12){//Accidentals Rate Uncertainty
	    		for(int iad=0; iad<maxAD; iad++){
	    			accRate_uncert[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==13){//Li9 Rate
	    		for(int iad=0; iad<maxAD; iad++){
	    			li9Rate[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==14){//Li9 Rate Uncertainty
	    		for(int iad=0; iad<maxAD; iad++){
	    			li9Rate_uncert[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==15){//Fast Neutron Rate
	    		for(int iad=0; iad<maxAD; iad++){
	    			fastnRate[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==16){//Fast Neutron Rate Uncertainty
	    		for(int iad=0; iad<maxAD; iad++){
	    			fastnRate_uncert[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==17){//AmC Rate
	    		for(int iad=0; iad<maxAD; iad++){
	    			amcRate[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==18){//AmC Rate Uncertainty
	    		for(int iad=0; iad<maxAD; iad++){
	    			amcRate_uncert[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==19){//Alpha-n Rate
	    		for(int iad=0; iad<maxAD; iad++){
	    			alphanRate[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==20){//Alpha-n Rate Uncertainty
	    		for(int iad=0; iad<maxAD; iad++){
	    			alphanRate_uncert[iad] = value[iad+2];
	    		}
	    	}
	    }
	    myfile.close();
	  }
	  else cout << "Unable to open file" << endl;
	  
//	  cout << "8AD period done. DAQ time for AD1=" << DAQ[0] <<" days.\t On to 6AD period." << endl;

	//Reading the config file
	  ifstream myfile_6 (Form("%s%dad.txt",toyConfig,stage[0])); //6AD period
		iline=0;
	  if (myfile_6.is_open())
	  {
	    while (getline(myfile_6, line, '\n')){
	    	iline+=1;
	    	if(iline<40) continue;
	    	
	    	//Breaking down the line
	    	for(int i=0; i< numColumns; i++){
	    		temp[i] = line.substr(0,line.find(tab));
	    		line.erase(0,line.find(tab)+tab.length());
	    		value[i] = stod(temp[i]);
	    	}
	    	
	    	if(value[1]==2){//DAQ time
	    		for(int iad=0; iad<maxAD; iad++){
	    			DAQ[iad][0] = value[iad+2];
				totDAQ[iad] += DAQ[iad][0];
	    		}
	    	}

	    }
	    myfile_6.close();
	  }
	  else cout << "Unable to open file" << endl;

//	  cout << "6AD period done. DAQ time for AD1=" << DAQ[0] <<" days.\t On to 7AD period." << endl;

	//Reading the config file
	  ifstream myfile_7 (Form("%s%dad.txt",toyConfig,stage[2])); //7AD period
		iline=0;
	  if (myfile_7.is_open())
	  {
	    while (getline(myfile_7, line, '\n')){
	    	iline+=1;
	    	if(iline<40) continue;
	    	
	    	//Breaking down the line
	    	for(int i=0; i< numColumns; i++){
	    		temp[i] = line.substr(0,line.find(tab));
	    		line.erase(0,line.find(tab)+tab.length());
	    		value[i] = stod(temp[i]);
	    	}
	    	
	    	if(value[1]==2){//DAQ time
	    		for(int iad=0; iad<maxAD; iad++){
	    			DAQ[iad][2] = value[iad+2];
				totDAQ[iad] += DAQ[iad][2];
	    		}
	    	}

	    }
	    myfile_7.close();
	  }
	  else cout << "Unable to open file" << endl;

//	  cout << "Total DAQ time for AD1=" << DAQ[0] <<" days." << endl;
}


void write_fitConfig(){
	parse_config();

 	ofstream fit_file;
 	fit_file.open (fitConfig); //Making the config file for the fitter

	int iline=0;
	string line;
	string slash="/";
	string comma=",";
	string quotes="\"";
	string colon=":";
	string bracket="[";
	  ifstream ex_file (exConfig); //Reading in the example config file
	  if (ex_file.is_open())
	  {
	    while (getline(ex_file, line, '\n')){
		iline+=1;
		if(iline==2){ //parameters database
			while(line.find(slash) != string::npos){
				fit_file << line.substr(0,line.find(slash)+slash.length());
	    			line.erase(0, line.find(slash)+slash.length());
	    		}
	    		fit_file << parameters_db << "\",\n";
			continue;
		}
		
		else if(iline==4){ //backgrounds database
			while(line.find(slash) != string::npos){
				fit_file << line.substr(0,line.find(slash)+slash.length());
	    			line.erase(0, line.find(slash)+slash.length());
	    		}
	    		fit_file << background_db << "\",\n";
			continue;
		}

		else if(iline==5){ //backgrounds counts
			fit_file << line.substr(0,line.find(colon)+colon.length());
	    		fit_file << " \"" << bkgd_counts_label << "\",\n";
			continue;
		}
		
		else if(iline==8){ //accidentals spectra
			fit_file << line.substr(0,line.find(quotes)+quotes.length()) << acc_spectra_label << "\",\n";
			continue;
		}
		
		else if(iline==9){ //background spectra
			fit_file << line.substr(0,line.find(quotes)+quotes.length()) << bkgd_spectra_label << "\"\n";
			continue;
		}
		
		else if(iline==13){ //multiplicity efficiency
			fit_file << line.substr(0,line.find(bracket)+bracket.length());
			for(int iad=0; iad<maxAD; iad++){
				if(iad != maxAD-1) fit_file << multEff[iad] << ",";
				else fit_file << multEff[iad];
			}
			fit_file << "],\n";
			continue;
		}
		
		else if(iline==15){ //muon veto efficiency
			fit_file << line.substr(0,line.find(bracket)+bracket.length());
			for(int iad=0; iad<maxAD; iad++){
				if(iad != maxAD-1) fit_file << muonEff[iad] << ",";
				else fit_file << muonEff[iad];
			}
			fit_file << "],\n";
			continue;
		}

		if(iline==21){ //prompt database
			while(line.find(slash) != string::npos){
				fit_file << line.substr(0,line.find(slash)+slash.length());
	    			line.erase(0, line.find(slash)+slash.length());
	    		}
	    		fit_file << prompt_db << "\",\n";
			continue;
		}

		else if(iline==22){ //prompt source
			fit_file << line.substr(0,line.find(colon)+colon.length());
	    		fit_file << " \"" << num_coincs_label << "\",\n";
			continue;
		}

		else if(iline==25){ //prompt source
			fit_file << line.substr(0,line.find(colon)+colon.length());
	    		fit_file << " \"" << det_resp_label << "\",\n";
			continue;
		}

		else{
			fit_file << line << endl;
			continue;
		}
	    }
	    ex_file.close();
	  }
	  else cout << "Unable to open file" << endl; //End of looking at the example config file
	  
	  fit_file.close(); //Closing the fitter config file

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


void fill_livetimes(){

	parse_config();

	  sqlite3 *db_livetime;
	  char *zErrMsg_livetime = 0;
	  int rc_livetime;
	  char *sql_livetime;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_livetime = sqlite3_open(parameters_db, &db_livetime);

	  if( rc_livetime ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_livetime));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
	  cout << "Database:\t" << parameters_db << endl;

	  //First need to delete the previous entries!!!!
	  sql_livetime = Form("DELETE FROM %s;",livetime_table);
	  rc_livetime = sqlite3_exec(db_livetime, sql_livetime, callback, (void*)data, &zErrMsg_livetime);
	  
/*	for(int iad=0; iad<maxAD; iad++){
		  // Create SQL statement //
		  sql_livetime = Form("INSERT OR REPLACE INTO %s VALUES (%d,%d,\"%s\",%ld,%f,%f,%f);",livetime_table,EH[iad],AD[iad],muon_label,long(muon_rate[iad]*totDAQ[iad]*86400),totDAQ[iad]*muonEff[iad]*86400*1.e9,muon_rate[iad],muonEff[iad]);

		cout << sql_livetime << endl;
		
		   // Execute SQL statement //
		   rc_livetime = sqlite3_exec(db_livetime, sql_livetime, callback, (void*)data, &zErrMsg_livetime);
	}*/
	
	for(int istage=0; istage<numStages; istage++){
		for(int iad=0; iad<maxAD; iad++){
			  // Create SQL statement //
			  sql_livetime = Form("INSERT OR REPLACE INTO %s VALUES (%d%d,%d,\"%s\",%ld,%f,%f,%f);",livetime_table,EH[iad],stage[istage],AD[iad],muon_label,long(muon_rate[iad]*DAQ[iad][istage]*86400),DAQ[iad][istage]*muonEff[iad]*86400*1.e9,muon_rate[iad],muonEff[iad]);

			cout << sql_livetime << endl;
//			if(EH[iad]==1 && AD[iad]==1 && stage[istage]==7) continue; //Skip for 7AD period of EH1AD1
//			if(EH[iad]==2 && AD[iad]==2 && stage[istage]==6) continue; //Skip for 6AD period of EH2AD2
//			if(EH[iad]==3 && AD[iad]==4 && stage[istage]==6) continue; //Skip for 6AD period of EH3AD4
			   // Execute SQL statement //
			   rc_livetime = sqlite3_exec(db_livetime, sql_livetime, callback, (void*)data, &zErrMsg_livetime);
		}
	}
	   if( rc_livetime != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_livetime);
	     sqlite3_free(zErrMsg_livetime);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_livetime);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;



}


void fill_bkgd_counts(){

	parse_config();

	  sqlite3 *db_bkgd;
	  char *zErrMsg_bkgd = 0;
	  int rc_bkgd;
	  char *sql_bkgd;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_bkgd = sqlite3_open(background_db, &db_bkgd);

	  if( rc_bkgd ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_bkgd));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
	  cout << "Database:\t" << background_db << endl;
	
	for(int iad=0; iad<maxAD; iad++){ //ACCIDENTALS
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"accidental\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],accRate[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad],accRate_uncert[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad]); 
		//Check it
		cout << sql_bkgd << endl;
		   // Execute SQL statement
		   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
	}
	
	for(int iad=0; iad<maxAD; iad++){ //Li9
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"li9\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],li9Rate[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad],li9Rate_uncert[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad]); 
		//Check it
		cout << sql_bkgd << endl;
		   // Execute SQL statement
		   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
	}

	for(int iad=0; iad<maxAD; iad++){ //Fast-n
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"fast-neutron\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],fastnRate[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad],fastnRate_uncert[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad]); 
		//Check it
		cout << sql_bkgd << endl;
		   // Execute SQL statement
		   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
	}

	for(int iad=0; iad<maxAD; iad++){ //Amc
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"amc\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],amcRate[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad],amcRate_uncert[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad]); 
		//Check it
		cout << sql_bkgd << endl;
		   // Execute SQL statement
		   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
	}

	for(int iad=0; iad<maxAD; iad++){ //Radiogenic neutrons = alpha n .......
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"rad-n\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],alphanRate[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad],alphanRate_uncert[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad]); 
		//Check it
		cout << sql_bkgd << endl;
		   // Execute SQL statement
		   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
	}


	
	   if( rc_bkgd != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_bkgd);
	     sqlite3_free(zErrMsg_bkgd);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_bkgd);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;



}

void fill_prompt(){

	parse_config();
	
	TFile *in_file = new TFile(toy_file);

	  sqlite3 *db_prompt;
	  char *zErrMsg_prompt = 0;
	  int rc_prompt;
	  char *sql_prompt;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_prompt = sqlite3_open(prompt_db, &db_prompt);

	  if( rc_prompt ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_prompt));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
	  cout << "Database:\t" << prompt_db << endl;
	
	for(int iad=0; iad<maxAD; iad++){
		TH1F *in_hist = (TH1F*)in_file->Get(Form("%s1_ad%d",toy_hist,iad+1)); //Stage1 hist
		for(int istage=2; istage <=3; istage++){
			in_hist->Add((TH1F*)in_file->Get(Form("%s%d_ad%d",toy_hist,istage,iad+1))); //Adding Stage2 and Stage3 hists
		}
		// Create SQL statement
		sql_prompt = Form("INSERT OR REPLACE INTO %s VALUES (%d,%d,1,\"[",num_coincs_table,EH[iad],AD[iad]);
		for(int iBin=1; iBin<numBins+1; iBin++){
			if(iBin == numBins) sql_prompt = Form("%s%f]\",\"%s\");",sql_prompt,in_hist->GetBinContent(iBin),num_coincs_label);
			else sql_prompt = Form("%s%f,",sql_prompt,in_hist->GetBinContent(iBin));
		}
		//Check it
		cout << sql_prompt << endl;
		   // Execute SQL statement
		   rc_prompt = sqlite3_exec(db_prompt, sql_prompt, callback, (void*)data, &zErrMsg_prompt);
	}


	
	   if( rc_prompt != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_prompt);
	     sqlite3_free(zErrMsg_prompt);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_prompt);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;

}

void fill_acc_spec(){

	parse_config();
	
	TFile *in_file = new TFile(acc_spec_file);

	  sqlite3 *db_acc_spec;
	  char *zErrMsg_acc_spec = 0;
	  int rc_acc_spec;
	  char *sql_acc_spec;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_acc_spec = sqlite3_open(background_db, &db_acc_spec);

	  if( rc_acc_spec ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_acc_spec));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
	  cout << "Database:\t" << background_db << endl;
	
	for(int iad=0; iad<maxAD; iad++){
	
		TH1F *in_hist = (TH1F*)in_file->Get(Form("%s_eh%d_ad%d",acc_spec_hist,EH[iad],AD[iad]));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spectrum=new TH1F(Form("h_rebinned_%d",iad+1),Form("h_rebinned_%d",iad+1),numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){ //+1 not +2
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		}
		
		double spec_integral= h_rebinned_spectrum->Integral();
		
		for(int iBin=1; iBin<numBins+1; iBin++){
			// Create SQL statement
			sql_acc_spec = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,1,%d,%f);",acc_spec_table,acc_spectra_label,EH[iad],AD[iad],iBin-1,(h_rebinned_spectrum->GetBinContent(iBin))/spec_integral);
			//Check it
			cout << sql_acc_spec << endl;
			   // Execute SQL statement
			   rc_acc_spec = sqlite3_exec(db_acc_spec, sql_acc_spec, callback, (void*)data, &zErrMsg_acc_spec);
		}
	}


	
	   if( rc_acc_spec != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_acc_spec);
	     sqlite3_free(zErrMsg_acc_spec);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_acc_spec);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;

}


void fill_li9_spec(){

	parse_config();
	
	TFile *in_file = new TFile(li9_spec_file);
	
		TH1F *in_hist = (TH1F*)in_file->Get(Form("%s",li9_spec_hist));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spectrum=new TH1F("h_rebinned","h_rebinned",numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		}

		double spec_integral= h_rebinned_spectrum->Integral();

	  sqlite3 *db_li9_spec;
	  char *zErrMsg_li9_spec = 0;
	  int rc_li9_spec;
	  char *sql_li9_spec;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_li9_spec = sqlite3_open(background_db, &db_li9_spec);

	  if( rc_li9_spec ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_li9_spec));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
	  cout << "Database:\t" << background_db << endl;


		for(int iBin=1; iBin<numBins+1; iBin++){
			// Create SQL statement
			sql_li9_spec = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",1,%d,%f);",li9_spec_table,bkgd_spectra_label,iBin-1,(h_rebinned_spectrum->GetBinContent(iBin))/spec_integral);
			//Check it
			cout << sql_li9_spec << endl;

			   // Execute SQL statement
			   rc_li9_spec = sqlite3_exec(db_li9_spec, sql_li9_spec, callback, (void*)data, &zErrMsg_li9_spec);
		}


	
	   if( rc_li9_spec != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_li9_spec);
	     sqlite3_free(zErrMsg_li9_spec);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_li9_spec);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;

}


void fill_fastn_spec(){

	parse_config();
	
	TFile *in_file = new TFile(fastn_spec_file);
	
		TH1F *in_hist = (TH1F*)in_file->Get(Form("%s",fastn_spec_hist));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spectrum=new TH1F("h_rebinned","h_rebinned",numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		}

		double spec_integral= h_rebinned_spectrum->Integral();

	  sqlite3 *db_fastn_spec;
	  char *zErrMsg_fastn_spec = 0;
	  int rc_fastn_spec;
	  char *sql_fastn_spec;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_fastn_spec = sqlite3_open(background_db, &db_fastn_spec);

	  if( rc_fastn_spec ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_fastn_spec));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
	  cout << "Database:\t" << background_db << endl;


		for(int iBin=1; iBin<numBins+1; iBin++){
			// Create SQL statement
			sql_fastn_spec = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",1,%d,%f);",fastn_spec_table,bkgd_spectra_label,iBin-1,(h_rebinned_spectrum->GetBinContent(iBin))/spec_integral);
			//Check it
			cout << sql_fastn_spec << endl;
			
			   // Execute SQL statement
			   rc_fastn_spec = sqlite3_exec(db_fastn_spec, sql_fastn_spec, callback, (void*)data, &zErrMsg_fastn_spec);
		}



	
	   if( rc_fastn_spec != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_fastn_spec);
	     sqlite3_free(zErrMsg_fastn_spec);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_fastn_spec);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;

}

void fill_amc_spec(){

	parse_config();
	
	TFile *in_file = new TFile(amc_spec_file);
	
		TH1F *in_hist = (TH1F*)in_file->Get(Form("%s",amc_spec_hist));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spectrum=new TH1F("h_rebinned","h_rebinned",numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		}

		double spec_integral= h_rebinned_spectrum->Integral();

	  sqlite3 *db_amc_spec;
	  char *zErrMsg_amc_spec = 0;
	  int rc_amc_spec;
	  char *sql_amc_spec;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_amc_spec = sqlite3_open(background_db, &db_amc_spec);

	  if( rc_amc_spec ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_amc_spec));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
	  cout << "Database:\t" << background_db << endl;


		for(int iBin=1; iBin<numBins+1; iBin++){
			// Create SQL statement
			sql_amc_spec = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",1,%d,%f);",amc_spec_table,bkgd_spectra_label,iBin-1,(h_rebinned_spectrum->GetBinContent(iBin))/spec_integral);
			//Check it
			cout << sql_amc_spec << endl;

			   // Execute SQL statement
			   rc_amc_spec = sqlite3_exec(db_amc_spec, sql_amc_spec, callback, (void*)data, &zErrMsg_amc_spec);
		}



	
	   if( rc_amc_spec != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_amc_spec);
	     sqlite3_free(zErrMsg_amc_spec);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_amc_spec);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;

}


void fill_radn_spec(){

	parse_config();
	
	TFile *in_file = new TFile(radn_spec_file);
	
		TH1F *in_hist = (TH1F*)in_file->Get(Form("%s",radn_spec_hist));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spectrum=new TH1F("h_rebinned","h_rebinned",numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		}

		double spec_integral= h_rebinned_spectrum->Integral();

	  sqlite3 *db_radn_spec;
	  char *zErrMsg_radn_spec = 0;
	  int rc_radn_spec;
	  char *sql_radn_spec;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_radn_spec = sqlite3_open(background_db, &db_radn_spec);

	  if( rc_radn_spec ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_radn_spec));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }
	  
	  cout << "Database:\t" << background_db << endl;


		for(int iBin=1; iBin<numBins+1; iBin++){
			// Create SQL statement
			sql_radn_spec = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",1,%d,%f);",radn_spec_table,bkgd_spectra_label,iBin-1,(h_rebinned_spectrum->GetBinContent(iBin))/spec_integral);
			//Check it
			cout << sql_radn_spec << endl;

			   // Execute SQL statement
			   rc_radn_spec = sqlite3_exec(db_radn_spec, sql_radn_spec, callback, (void*)data, &zErrMsg_radn_spec);
		}



	
	   if( rc_radn_spec != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_radn_spec);
	     sqlite3_free(zErrMsg_radn_spec);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_radn_spec);

	  std::cout << " Test_sqlite3: ended " << std::endl;


	 ab:;

}


void plot_sub(int iad){

	parse_config();
	
	TFile *in_file_prompt = new TFile(toy_file);
		TH1F *in_hist_prompt = (TH1F*)in_file_prompt->Get(Form("%s1_ad%d",toy_hist,iad+1));
		for(int istage=2; istage <=3; istage++){
			in_hist_prompt->Add((TH1F*)in_file_prompt->Get(Form("%s%d_ad%d",toy_hist,istage,iad+1))); //Adding Stage2 and Stage3 hists
		}
		TH1F *h_sub;
		h_sub = (TH1F*)in_hist_prompt->Clone();
		cout << "Number of events in prompt spectrum:\t" << in_hist_prompt->Integral() << endl;
	
	TFile *in_file_acc = new TFile(acc_spec_file);
		TH1F *in_hist_acc = (TH1F*)in_file_acc->Get(Form("%s_eh%d_ad%d",acc_spec_hist,EH[iad],AD[iad]));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spectrum=new TH1F("h_rebinned","h_rebinned",numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist_acc->GetNbinsX())+2; iBin++){
			if((in_hist_acc->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist_acc->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spectrum->Fill(in_hist_acc->GetBinCenter(iBin),in_hist_acc->GetBinContent(iBin));
		}
		//Scale:
		h_rebinned_spectrum->Scale((accRate[iad]*multEff[iad]*muonEff[iad]*totDAQ[iad])/(h_rebinned_spectrum->Integral()));
		

	//Writing to the file
	TFile* outfile=new TFile(Form("./ckDB/inputSub_eh%dad%d.root",EH[iad],AD[iad]), "RECREATE");
	outfile->cd();
		in_hist_prompt->Write();

		h_rebinned_spectrum->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_rebinned_spectrum->GetYaxis()->SetTitle("Counts");
		h_rebinned_spectrum->Write();

		h_sub->Add(h_rebinned_spectrum,-1);
		h_sub->Write();

}




void all(){
//	write_fitConfig();
	fill_livetimes();
//	fill_bkgd_counts();
//	fill_prompt();
//	fill_acc_spec();
//	fill_li9_spec();
//	fill_fastn_spec();
//	fill_amc_spec();
//	fill_radn_spec();
}
