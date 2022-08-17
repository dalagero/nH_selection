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
int EH[maxAD]={1,1,2,2,3,3,3,3};
int AD[maxAD]={1,2,1,2,1,2,3,4};
const int numBins=34;
const int numColumns=10;
const int numPeriods=3;
const int period[numPeriods]={6,8,7};

char suffix[64];

char toyConfig[128];
char exConfig[128];
char fitConfig[128];

char parameters_db[64];
char background_db[64];
char prompt_db[64];

char bkgd_counts_label[64];
char bkgd_spectra_label[64];
char acc_spectra_label[64];
char num_coincs_label[64];
char det_resp_label[64];
char muon_label[64];

char livetime_table[64];
char bkgd_table[64];
char acc_spec_table[64];
char li9_spec_table[64];
char fastn_spec_table[64];
char amc_spec_table[64];
char radn_spec_table[64];
char num_coincs_table[64];

char data_file[64];
char data_hist[64];
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
double DAQ[maxAD];
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
//livetime in ns
const long livetime_NU_6[maxAD] = {8735696643816849,8706935160442762,10113181067367890,0,15441985656608942,15435958971995376,15435603078800936,0};
const long livetime_NU_8[maxAD] = {63508588292383088,63179191034858984,73224868299200272,73158035484463568,111480011143546768,111475373079604368,111457215448017744,111470651623536112}; 
const long livetime_NU_7[maxAD] = {0,9415080550633996,10866176610885526,10843931304682508,16518459379144330,16518755844606340,16515902872001016,16518546445582180}; 
const long livetime_noNU_6[maxAD] = {8740488756690403,8711184680746069,10118954001172654,0,15442900588282114,15436269946704184,15436039997326514,0};
const long livetime_noNU_8[maxAD] = {63542674439044656,72630117585539184,84123946367271600,84031181570224704,128005061996741616,127998306940008000,127977498531740224,127994617784696048};
const long livetime_noNU_7[maxAD] = {0,9419728361800096,10868317383974194,10843039451090238,16519117096602220,16519252000999172,16516429412179362,16518802730496970};
const int maxIdent=3;
const int noNU=0;
const int NU=1;
const int Sam=2;
double livetime[maxIdent][maxAD][numPeriods];
double tot_livetime[maxIdent][maxAD];
double Ed_eff[maxIdent] = {0.934,0.936};
double Ed_err[maxIdent] = {0.0021,0.0023};
double DT_eff[maxIdent] = {0.706,0.706};
double DT_err[maxIdent] = {0.0035,0.0034};
double tarp_relErr = 0.0037;
int identifier;
const double muon_rate[maxAD] = {200.32, 200.32, 150.08, 149.80, 15.748, 15.748, 15.748, 15.747}; //rate in Hz
const int DTcut = 800;
const double convertToDays = 1.e-9/(60*60*24);


void init_names(){ //FIXME: Organize better.... This is pretty ugly.
	for(int iad=0; iad<maxAD; iad++){
		livetime[noNU][iad][0]=livetime_noNU_6[iad];
		livetime[noNU][iad][1]=livetime_noNU_8[iad];
		livetime[noNU][iad][2]=livetime_noNU_7[iad];
		livetime[NU][iad][0]=livetime_NU_6[iad];
		livetime[NU][iad][1]=livetime_NU_8[iad];
		livetime[NU][iad][2]=livetime_NU_7[iad];
		livetime[Sam][iad][0]=0;
		livetime[Sam][iad][1]=0;
		livetime[Sam][iad][2]=0;
		for(int thisIdent=0; thisIdent<maxIdent; thisIdent++){
			tot_livetime[thisIdent][iad]=0;
		}
	}
	if(strcmp(suffix,"_NU")==0) identifier=NU;
	if(strcmp(suffix,"_Sam")==0) identifier=Sam;
	else identifier=noNU;

	if(strcmp(suffix,"_NU")==0) sprintf(toyConfig,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/input/Theta13-inputs_P17B_inclusive_8ad.txt");
	else if(strcmp(suffix,"_Sam")==0) sprintf(toyConfig,"/global/homes/d/dalagero/dyb-event-selection/dyb_analysis/fitter/data/Theta13-inputs_P17B_inclusive_8ad_Sam.txt");
	else sprintf(toyConfig,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/input/Theta13-inputs_P17B_inclusive_8ad_noNU.txt");
	sprintf(exConfig,"fit_config_nH_resid_flash_example.json");
	sprintf(fitConfig,"fit_config_nH_resid_flash%s.json",suffix);
	
	sprintf(parameters_db,"parameters.db");
		if(strcmp(suffix,"_Sam")==0) sprintf(parameters_db,"parameters_Sam.db");
	sprintf(background_db,"nH_backgrounds.db");
	sprintf(prompt_db,"adtime_dt_eff_study.db");

	sprintf(bkgd_counts_label,"Olivia%s 8/10/2022",suffix);
//	sprintf(bkgd_spectra_label,"nH backgrounds, Olivia%s 6/14/2022",suffix); //FIXME and the 2 following rows: Something is up here and doing a bad job.
		if(identifier==NU) sprintf(bkgd_spectra_label,"nH backgrounds, Olivia 3/15/2022 NU");
		else sprintf(bkgd_spectra_label,"nH backgrounds, Olivia 1/21/2022");
	sprintf(acc_spectra_label,"accidentals, Olivia%s 6/14/2022",suffix);
	sprintf(num_coincs_label,"Olivia%s 8/10/2022",suffix);
	sprintf(det_resp_label,"Matrix for fitter from Beda, 8/1/22, nH modified 2 binning");
	sprintf(muon_label,"Olivia%s 8/10/2022",suffix);

	sprintf(livetime_table,"muon_rates");
	sprintf(bkgd_table,"bg_counts");
	sprintf(acc_spec_table,"accidentals_spectrum");
	sprintf(li9_spec_table,"li9_spectrum");
	sprintf(fastn_spec_table,"fast_neutron_spectrum");
	sprintf(amc_spec_table,"amc_spectrum");
	sprintf(radn_spec_table,"rad_n_spectrum");
	sprintf(num_coincs_table,"num_coincidences");

	sprintf(data_file,"./data/TotaledPlots%s_EH", suffix);
	sprintf(acc_spec_file,"./data/TotaledSingles%s_1500_EH",suffix);
	sprintf(li9_spec_file,"/global/homes/d/dalagero/dyb-event-selection/dyb_analysis/fitter/data/Li9He8shape.root");
	sprintf(fastn_spec_file,"/global/homes/d/dalagero/dyb-event-selection/dyb_analysis/fitter/data/FastNshape.root");
	sprintf(amc_spec_file,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/amc_spectrum/amc_spectrum.root");
	sprintf(radn_spec_file,"/global/homes/d/dalagero/fromBeda/DybBerkFit-master/muon_decay_spectrum/MuonDecaySpec.root");

	sprintf(data_hist,"%s",Form("h_total_prompt_energy_DT%d_3sig_ad",DTcut));
	sprintf(acc_spec_hist,"%s",Form("h_total_prompt_energy_DT%d_3sig_scaled_ad",DTcut));
	sprintf(li9_spec_hist,"Li9 prompt energy");
	sprintf(fastn_spec_hist,"OWP_EH%d_nH_Ep1.5_DT0.8m_Ep_EpLt12MeV",EH[0]); //decide what to do with the fixed EH number
	sprintf(amc_spec_hist,"h_toy");
	sprintf(radn_spec_hist,"MdSpec_EH1");


/*	cout << "Toy Config File: " << toyConfig << endl;
	cout << "Fit Config File: " << fitConfig << endl;
	cout << "Data File: " << data_file << endl;
	cout << "Acc File: " << acc_spec_file << endl;*/

}

void parse_config(){
	init_names();
	
	//Reading the config file
	  ifstream myfile (toyConfig);
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
	    			DAQ[iad] = value[iad+2];
	    		}
	    	}
	    	if(value[1]==3){//Muon Veto Efficiency
	    		for(int iad=0; iad<maxAD; iad++){
	    			muonEff[iad] = value[iad+2];
				if(identifier==Sam) livetime[identifier][iad][1]=DAQ[iad]*muonEff[iad]/convertToDays;
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

		else if(iline==17){ //Detection efficiency (DT & Ed)
			fit_file << line.substr(0,line.find(bracket)+bracket.length());
			for(int iad=0; iad<maxAD; iad++){
				if(iad != maxAD-1) fit_file << Ed_eff[identifier]*DT_eff[identifier] << ",";
				else fit_file << Ed_eff[identifier]*DT_eff[identifier];
			}
			fit_file << "],\n";
			continue;
		}

		else if(iline==21){ //prompt database
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

		else if(iline==25){ //detector response matrix
			fit_file << line.substr(0,line.find(colon)+colon.length());
	    		fit_file << " \"" << det_resp_label << "\",\n";
			continue;
		}

		else if(iline==30){ //detection efficiency error
			fit_file << line.substr(0,line.find(colon)+colon.length());
	    		fit_file << " " << sqrt(pow((Ed_err[identifier]/Ed_eff[identifier]),2)+pow((DT_err[identifier]/DT_eff[identifier]),2)+pow(tarp_relErr,2)) << ",\n";
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
	
	identifier=NU;

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
	
	for(int iad=0; iad<maxAD; iad++){
		for(int iperiod=0; iperiod<numPeriods; iperiod++){
		  /* Create SQL statement */
		  sql_livetime = Form("INSERT OR REPLACE INTO %s VALUES (%d%d,%d,\"%s\",%ld,%ld,%f,%f);",livetime_table,EH[iad],period[iperiod],AD[iad],muon_label,long(muon_rate[iad]*DAQ[iad]*86400),long(livetime[identifier][iad][iperiod]),muon_rate[iad],muonEff[iad]);

		cout << sql_livetime << endl;

		   // Execute SQL statement
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

	for(int iad=0; iad<maxAD; iad++){
		cout << "***Livetime in ns for EH" << EH[iad] << "-AD" << AD[iad] << "*** 6AD: " << livetime[identifier][iad][0] << "\t8AD: " << livetime[identifier][iad][1] << "\t7AD: " << livetime[identifier][iad][2] << endl;
	}

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

	for(int iad=0; iad<maxAD; iad++){
		  for(int iperiod=0; iperiod<numPeriods; iperiod++){
		  	for(int thisIdent=0; thisIdent<maxIdent; thisIdent++){
		  		tot_livetime[thisIdent][iad]+=livetime[thisIdent][iad][iperiod];
		  	}
		  }
	}

	if(strcmp(suffix,"_Sam")==0){
		for(int iad=0; iad<maxAD; iad++){ //Accidentals for Sam's numbers
			  // Create SQL statement
			  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"accidental\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],accRate[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays,accRate_uncert[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays);  
			//Check it
			cout << sql_bkgd << endl;
			   // Execute SQL statement
			   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
		}
	}
	else{
		for(int iad=0; iad<maxAD; iad++){ //ACCIDENTALS
			TFile *in_file = new TFile(Form("%s%d.root",acc_spec_file,EH[iad]));
			TH1F *in_hist = (TH1F*)in_file->Get(Form("%s%d",acc_spec_hist,AD[iad]));
			double counts = 0;
			counts = in_hist->Integral();
			double scale = 0;
			TH1F *scaled_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_scaled_ad%d",AD[iad]));
			TH1F *before_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_before_ad%d",AD[iad]));
			scale = (scaled_hist->Integral())/(before_hist->Integral());

			  // Create SQL statement
			  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"accidental\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],counts,sqrt(counts/scale)*scale); 
			//Check it
			cout << sql_bkgd << endl;
			   // Execute SQL statement
			   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
			
			in_file->Close();
		}
	}
	
	for(int iad=0; iad<maxAD; iad++){ //Li9
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"li9\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],li9Rate[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays,li9Rate_uncert[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays); 
		//Check it
		cout << sql_bkgd << endl;
		   // Execute SQL statement
		   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
	}

	for(int iad=0; iad<maxAD; iad++){ //Fast-n
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"fast-neutron\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],fastnRate[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays,fastnRate_uncert[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays); 
		//Check it
		cout << sql_bkgd << endl;
		   // Execute SQL statement
		   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
	}

	for(int iad=0; iad<maxAD; iad++){ //Amc
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"amc\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],amcRate[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays,amcRate_uncert[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays); 
		//Check it
		cout << sql_bkgd << endl;
		   // Execute SQL statement
		   rc_bkgd = sqlite3_exec(db_bkgd, sql_bkgd, callback, (void*)data, &zErrMsg_bkgd);
	}

	for(int iad=0; iad<maxAD; iad++){ //Radiogenic neutrons = alpha n .......
		  // Create SQL statement
		  sql_bkgd = Form("INSERT OR REPLACE INTO %s VALUES (\"%s\",%d,%d,\"rad-n\",%f,%f);",bkgd_table,bkgd_counts_label,EH[iad],AD[iad],alphanRate[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays,alphanRate_uncert[iad]*multEff[iad]*tot_livetime[identifier][iad]*convertToDays); 
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
		TFile *in_file = new TFile(Form("%s%d_1500.root",data_file,EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get(Form("%s%d",data_hist,AD[iad]));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spectrum=new TH1F(Form("h_rebinned_%d",iad+1),Form("h_rebinned_%d",iad+1),numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		}
		
		// Create SQL statement
		sql_prompt = Form("INSERT OR REPLACE INTO %s VALUES (%d,%d,1,\"[",num_coincs_table,EH[iad],AD[iad]);
		for(int iBin=1; iBin<numBins+1; iBin++){
			if(iBin == numBins) sql_prompt = Form("%s%f]\",\"%s\");",sql_prompt,h_rebinned_spectrum->GetBinContent(iBin),num_coincs_label);
			else sql_prompt = Form("%s%f,",sql_prompt,h_rebinned_spectrum->GetBinContent(iBin));
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
		TFile *in_file = new TFile(Form("%s%d.root",acc_spec_file,EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get(Form("%s%d",acc_spec_hist,AD[iad]));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spectrum=new TH1F(Form("h_rebinned_%d",iad+1),Form("h_rebinned_%d",iad+1),numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
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
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
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
		}
		   // Execute SQL statement
		   rc_li9_spec = sqlite3_exec(db_li9_spec, sql_li9_spec, callback, (void*)data, &zErrMsg_li9_spec);


	
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
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
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
		}
		   // Execute SQL statement
		   rc_fastn_spec = sqlite3_exec(db_fastn_spec, sql_fastn_spec, callback, (void*)data, &zErrMsg_fastn_spec);


	
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
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
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
		}
		   // Execute SQL statement
		   rc_amc_spec = sqlite3_exec(db_amc_spec, sql_amc_spec, callback, (void*)data, &zErrMsg_amc_spec);


	
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
		for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
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
		}
		   // Execute SQL statement
		   rc_radn_spec = sqlite3_exec(db_radn_spec, sql_radn_spec, callback, (void*)data, &zErrMsg_radn_spec);


	
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
	
	TFile *in_file_prompt = new TFile(Form("%s%d_1500.root",data_file,EH[iad]));
		TH1F *in_hist_prompt = (TH1F*)in_file_prompt->Get(Form("%s%d",data_hist,AD[iad]));
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F* h_rebinned_spec_prompt=new TH1F("h_rebinned_prompt","h_rebinned_prompt",numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist_prompt->GetNbinsX())+1; iBin++){
			if((in_hist_prompt->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist_prompt->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spec_prompt->Fill(in_hist_prompt->GetBinCenter(iBin),in_hist_prompt->GetBinContent(iBin));
		}
	
	TFile *in_file_acc = new TFile(Form("%s%d.root",acc_spec_file,EH[iad]));
		TH1F *in_hist_acc = (TH1F*)in_file_acc->Get(Form("%s%d",acc_spec_hist,AD[iad]));
		TH1F* h_rebinned_spec_acc=new TH1F("h_rebinned_acc","h_rebinned_acc",numBins, binEdges);
		
		//Rebin:
		for(int iBin=1; iBin < (in_hist_acc->GetNbinsX())+1; iBin++){
			if((in_hist_acc->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((in_hist_acc->GetBinCenter(iBin))<binEdges[0]) continue;
			h_rebinned_spec_acc->Fill(in_hist_acc->GetBinCenter(iBin),in_hist_acc->GetBinContent(iBin));
		}

	TH1F *h_sub;
	h_sub = (TH1F*)h_rebinned_spec_prompt->Clone();
	cout << "Number of events in prompt spectrum:\t" << h_rebinned_spec_prompt->Integral() << endl;

	//Writing to the file
	TFile* outfile=new TFile(Form("./ckDB/inputSub_eh%dad%d.root",EH[iad],AD[iad]), "RECREATE");
	outfile->cd();
		h_rebinned_spec_prompt->Write();

		h_rebinned_spec_acc->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_rebinned_spec_acc->GetYaxis()->SetTitle("Counts");
		h_rebinned_spec_acc->Write();

		h_sub->Add(h_rebinned_spec_acc,-1);
		h_sub->Write();

}




void all(){
	write_fitConfig();
	fill_livetimes();
	fill_bkgd_counts();
	if(strcmp(suffix,"_Sam")!=0) fill_prompt();
	if(strcmp(suffix,"_Sam")!=0) fill_acc_spec(); //don't fill if it's Sam's inputs
/*	fill_li9_spec();
	fill_fastn_spec();
	fill_amc_spec();
	fill_radn_spec();*/
}

void process(string input){
	if(input == ""){
		sprintf(suffix, "");
		all();
	}
	else if(input == "NU"){
		sprintf(suffix,"_%s", input.c_str());
		cout << "Input: " << input << "\tSuffix: " << suffix << endl;
		all();
	}
	else if(input == "Sam"){
		sprintf(suffix,"_%s", input.c_str());
		cout << "Input: " << input << "\tSuffix: " << suffix << endl;
		all();
	}
/*	else if(input == "THU"){
		sprintf(suffix,"_%s", input.c_str());
		cout << "Input: " << input << "\tSuffix: " << suffix << endl; //change
		init_names();
	}*/
	else{
		cout << "Unknown input!! No information passed to the database!!" << endl;
	}



}
