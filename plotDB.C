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

const int maxAD=8;
int EH[maxAD]={1,1,2,2,3,3,3,3};
int AD[maxAD]={1,2,1,2,1,2,3,4};
const int numBins=34;

char prompt_db[64];
char bkgd_db[64];

char prompt_fileName[64];
char promptSource[64];
char prompt_output[128];

char bkgd_counts_fileName[64];
double bkgd_counts_vals[2];
char bkgd_spectra_fileName[64];
char BgName[64];
char BgLabel[64];
char BgLabel_counts[64];
char AccLabel[64];
int num_bkgds = 5;
double bkgd_spec_val;

char db_hists[64];

void init_config(string config_file){

	//The names of the output files
	sprintf(prompt_fileName,"./ckDB/prompt_spectra_data.txt");
	sprintf(bkgd_counts_fileName,"./ckDB/bkgd_counts_data.txt");
	sprintf(bkgd_spectra_fileName,"./ckDB/bkgd_spectra_data.txt");
	sprintf(db_hists,"./ckDB/db_spectra_data.root");

	//Reading the config file
	  ifstream myfile (config_file);
	  string line;
	  int iline=0;
	  string slash="/";
	  string quotes="\"";
	  string comma=",";
	  if (myfile.is_open())
	  {
	    while (getline(myfile, line, '\n')){
	    	iline+=1;
	    	if(iline == 4){//backgrounds database
	    		while(line.find(slash) != string::npos){
	    			line.erase(0, line.find(slash)+slash.length());
	    		}
	    		line.erase(line.find(quotes), line.find(quotes)+2);
	    		sprintf(bkgd_db,"%s", line.c_str());
	    		cout << bkgd_db << endl;
	    	}
	   	if(iline == 5){//background counts label
	   		line.erase(line.find(comma)-1, line.find(comma)+1);
	    		while(line.find(quotes) != string::npos){
	    			line.erase(0, line.find(quotes)+1);
	    		}
	    		sprintf(BgLabel_counts,"%s", line.c_str());
	    		cout << BgLabel_counts << endl;
	   	}
 		if(iline == 8){//accidentals spectra label
 			line.erase(0, line.find(comma)+2);
 			line.erase(line.find(comma)-1, line.find(comma)+1);
	    		sprintf(AccLabel,"%s", line.c_str());
			cout << AccLabel << endl;
		}
		if(iline == 9){//backgrounds spectra label
 			line.erase(0, line.find(comma)+2);
 			line.erase(line.find(quotes), line.find(quotes)+1);
	    		sprintf(BgLabel,"%s", line.c_str());
			cout << BgLabel << endl;
		}
		if(iline == 21){//prompt database
	    		while(line.find(slash) != string::npos){
	    			line.erase(0, line.find(slash)+slash.length());
	    		}
	    		line.erase(line.find(quotes), line.find(quotes)+2);
	    		sprintf(prompt_db,"%s", line.c_str());
	    		cout << prompt_db << endl;
		}
		if(iline == 22){//prompt label
	   		line.erase(line.find(comma)-1, line.find(comma)+1);
	    		while(line.find(quotes) != string::npos){
	    			line.erase(0, line.find(quotes)+1);
	    		}
	    		sprintf(promptSource,"%s", line.c_str());
	    		cout << promptSource << endl;
		}
	    }
	    myfile.close();
	  }
	  else cout << "Unable to open file" << endl;

}


static int callback_prompt(void *data, int argc, char **argv, char **azColName){
  int i;
  fprintf(stderr, "%s: ", (const char*)data);
//  ofstream myfile;
//  myfile.open (prompt_fileName);
   
  for(i = 0; i<argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
	if(i == 3){
//		myfile << argv[i] << endl;
		sprintf(prompt_output,"%s",argv[i]);
	}
  }
  
//    myfile.close();
    
  printf("\n");
  return 0;
}

static int callback_bkgd_counts(void *data, int argc, char **argv, char **azColName){
  int i;
  fprintf(stderr, "%s: ", (const char*)data);
   
  for(i = 0; i<argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
  }
  	bkgd_counts_vals[0] = stod(argv[4]);
  	bkgd_counts_vals[1] = stod(argv[5]);
    
  printf("\n");
  return 0;
}

static int callback_bkgd_spectra(void *data, int argc, char **argv, char **azColName){
  int i;
  fprintf(stderr, "%s: ", (const char*)data);
   
  for(i = 0; i<argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
  }
  	bkgd_spec_val = stod(argv[argc-1]);
    
  printf("\n");
  return 0;
}


void prompt_spectra(string config_file){
 	 ofstream myfile;
 	 myfile.open (prompt_fileName);

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


	for(int iad=0; iad<maxAD; iad++){
		init_config(config_file);
	
		  /* Create SQL statement */
		  sql_prompt = Form("select * from  num_coincidences where Source = \"%s\" and Hall=%d and DetNo=%d",promptSource,EH[iad],AD[iad]); 

		   /* Execute SQL statement */
		   rc_prompt = sqlite3_exec(db_prompt, sql_prompt, callback_prompt, (void*)data, &zErrMsg_prompt);
		   
	//	   cout << sql_prompt << endl;
		   myfile << EH[iad] << "\t" << AD[iad] << "\t" << prompt_output << endl;
	}   
	   if( rc_prompt != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_prompt);
	     sqlite3_free(zErrMsg_prompt);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_prompt);

	  std::cout << " Test_sqlite3: ended " << std::endl;

    myfile.close();

	 ab:;
}



void bkgd_counts(string config_file){

	init_config(config_file);

 	 ofstream myfile;
 	 myfile.open (bkgd_counts_fileName);

	  sqlite3 *db_bkgd_counts;
	  char *zErrMsg_bkgd_counts = 0;
	  int rc_bkgd_counts;
	  char *sql_bkgd_counts;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_bkgd_counts = sqlite3_open(bkgd_db, &db_bkgd_counts);

	  if( rc_bkgd_counts ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_bkgd_counts));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }

	for(int iad=0; iad<maxAD; iad++){
		for(int i=1; i<=num_bkgds; i++){
			if(i==1) sprintf(BgName,"accidental");
			if(i==2) sprintf(BgName,"amc");
			if(i==3) sprintf(BgName,"fast-neutron");
			if(i==4) sprintf(BgName,"li9");
			if(i==5) sprintf(BgName,"rad-n");
		  // Create SQL statement
		  sql_bkgd_counts = Form("select * from bg_counts where Label = \"%s\" and Hall=%d and DetNo=%d and BgName=\"%s\"",BgLabel_counts,EH[iad],AD[iad],BgName); 
		   // Execute SQL statement
		   rc_bkgd_counts = sqlite3_exec(db_bkgd_counts, sql_bkgd_counts, callback_bkgd_counts, (void*)data, &zErrMsg_bkgd_counts);
		   myfile << BgName << "\t" << EH[iad] << "\t" << AD[iad] << "\t" << bkgd_counts_vals[0] << "\t" << bkgd_counts_vals[1] << endl;
		}
	}
	
	   if( rc_bkgd_counts != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_bkgd_counts);
	     sqlite3_free(zErrMsg_bkgd_counts);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_bkgd_counts);

	  std::cout << " Test_sqlite3: ended " << std::endl;

    myfile.close();

	 ab:;
}

void bkgd_spectra(string config_file){

	init_config(config_file);
	
 	 ofstream myfile;
 	 myfile.open (bkgd_spectra_fileName);

	  sqlite3 *db_bkgd_spectra;
	  char *zErrMsg_bkgd_spectra = 0;
	  int rc_bkgd_spectra;
	  char *sql_bkgd_spectra;
	  const char* data = "Callback function called";

	  /* Open database */
	  rc_bkgd_spectra = sqlite3_open(bkgd_db, &db_bkgd_spectra);

	  if( rc_bkgd_spectra ) {
	    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db_bkgd_spectra));
	    goto ab;
	  } else {
	    fprintf(stderr, "Opened database successfully\n");
	  }


	for(int i=1; i<=num_bkgds; i++){
		if(i==1) sprintf(BgName,"accidentals_spectrum");
		if(i==2) sprintf(BgName,"amc_spectrum");
		if(i==3) sprintf(BgName,"fast_neutron_spectrum");
		if(i==4) sprintf(BgName,"li9_spectrum");
		if(i==5) sprintf(BgName,"rad_n_spectrum");

		if(i==1){
			for(int iad=0; iad<maxAD; iad++){
				myfile << BgName << "\t" << EH[iad] << "\t" << AD[iad];
				for(int iBin=0; iBin<numBins; iBin++){
					  // Create SQL statement
					  sql_bkgd_spectra = Form("select * from %s where Label = \"accidentals, %s\" and Hall=%d and DetNo=%d and BinIndex=%d",BgName,AccLabel,EH[iad],AD[iad],iBin);
					   // Execute SQL statement
					   rc_bkgd_spectra = sqlite3_exec(db_bkgd_spectra, sql_bkgd_spectra, callback_bkgd_spectra, (void*)data, &zErrMsg_bkgd_spectra);
					   myfile << " " << bkgd_spec_val;
				}
				myfile << endl;
			}
		}
		else{
			myfile << BgName << "\t" << EH[0] << "\t" << AD[0];
			for(int iBin=0; iBin<numBins; iBin++){
				  // Create SQL statement
				  sql_bkgd_spectra = Form("select * from %s where Label = \"nH backgrounds, %s\" and BinIndex=%d",BgName,BgLabel,iBin);
				   // Execute SQL statement
				   rc_bkgd_spectra = sqlite3_exec(db_bkgd_spectra, sql_bkgd_spectra, callback_bkgd_spectra, (void*)data, &zErrMsg_bkgd_spectra);
				   myfile << " " << bkgd_spec_val;
			}
			myfile << endl;
		}
	}
	
	   if( rc_bkgd_spectra != SQLITE_OK ) {
	     fprintf(stderr, "SQL error: %s\n", zErrMsg_bkgd_spectra);
	     sqlite3_free(zErrMsg_bkgd_spectra);
	   } else {
	     fprintf(stdout, "Operation done successfully\n");
	   }  

	  sqlite3_close(db_bkgd_spectra);

	  std::cout << " Test_sqlite3: ended " << std::endl;

    myfile.close();

	 ab:;
}

void readOutput(string config_file){

	init_config(config_file);

	//Building the histograms
	char name[64];
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	TH1F* h_prompt[maxAD];
	TH1F* h_acc[maxAD];
	TH1F* h_amc[maxAD];
	TH1F* h_fastn[maxAD];
	TH1F* h_li9[maxAD];
	TH1F* h_radn[maxAD];
	TH1F* h_ibds[maxAD];

	for(int iad=0; iad<maxAD; ++iad){
		sprintf(name, "h_numCoincs_eh%dad%d", EH[iad],AD[iad]); //Prompt
		h_prompt[iad]=new TH1F(name,name,numBins, binEdges);
		
		sprintf(name, "h_acc_eh%dad%d", EH[iad],AD[iad]); //Accidentals
		h_acc[iad]=new TH1F(name,name,numBins, binEdges);

		sprintf(name, "h_amc_eh%dad%d", EH[iad],AD[iad]); //Amc
		h_amc[iad]=new TH1F(name,name,numBins, binEdges);

		sprintf(name, "h_fastn_eh%dad%d", EH[iad],AD[iad]); //Fast-n
		h_fastn[iad]=new TH1F(name,name,numBins, binEdges);

		sprintf(name, "h_li9_eh%dad%d", EH[iad],AD[iad]); //li9
		h_li9[iad]=new TH1F(name,name,numBins, binEdges);

		sprintf(name, "h_radn_eh%dad%d", EH[iad],AD[iad]); //rad-n
		h_radn[iad]=new TH1F(name,name,numBins, binEdges);
		
		sprintf(name, "h_ibds_eh%dad%d", EH[iad],AD[iad]); //ibds: all backgrounds subtracted from the prompt spectra
		h_ibds[iad]=new TH1F(name,name,numBins, binEdges);
	}
	
	  int hall=0;
	  int det=0;
	  string line;
	  string delimiter = ",";
	  string values[numBins];

	  double prompt_vals[maxAD][numBins]; //Start of Prompt
	  ifstream myfile (prompt_fileName);
	  if (myfile.is_open())
	  {
	    while (myfile >> hall >> det >> line){
		for(int i=0; i<numBins; i++){
		    values[i] = line.substr(0, line.find(delimiter));
		    if(i==0) values[i].erase(0,1);
		    if(i==33) values[i].erase(line.find("]")-1,line.find("]"));
		    line.erase(0, line.find(delimiter) + delimiter.length());
		}
		for(int i=0; i<numBins; i++){
	    		prompt_vals[2*(hall-1)+det-1][i] = stod(values[i]);
			h_prompt[2*(hall-1)+det-1]->SetBinContent(i+1, prompt_vals[2*(hall-1)+det-1][i]);
		}
	    }
	    myfile.close();
	  }
	  else cout << "Unable to open file" << endl; //End of Prompt


//Start of Background Counts
		double counts=0;
		double error=0;
		double BgCounts[maxAD][num_bkgds];
		double BgError[maxAD][num_bkgds];
		int iline_counts=0;
	  ifstream myfile_bkgd_counts (bkgd_counts_fileName);
	  if (myfile_bkgd_counts.is_open()){
	    while (myfile_bkgd_counts >> line >> hall >> det >> counts >> error){
		BgCounts[2*(hall-1)+det-1][iline_counts%5]=counts;
		BgError[2*(hall-1)+det-1][iline_counts%5]=error;
		iline_counts+=1;
	    }
	    myfile_bkgd_counts.close();
	  }
	  else cout << "Unable to open file" << endl; //End of Background Counts

	//Start of Backgrounds Spectra
	  int iline=0;
	  string bg_line;
	  double bg_values[numBins];
	  ifstream bkgd_file (bkgd_spectra_fileName);
	  if (bkgd_file.is_open())
	  {
	    while (bkgd_file >> bg_line >> hall >> det >> bg_values[0] >> bg_values[1] >> bg_values[2] >> bg_values[3] >> bg_values[4] >> bg_values[5] >> bg_values[6] >> bg_values[7] >> bg_values[8] >> bg_values[9] >> bg_values[10] >> bg_values[11] >> bg_values[12] >> bg_values[13] >> bg_values[14] >> bg_values[15] >> bg_values[16] >> bg_values[17] >> bg_values[18] >> bg_values[19] >> bg_values[20] >> bg_values[21] >> bg_values[22] >> bg_values[23] >> bg_values[24] >> bg_values[25] >> bg_values[26] >> bg_values[27] >> bg_values[28] >> bg_values[29] >> bg_values[30] >> bg_values[31] >> bg_values[32] >> bg_values[33]){
		iline+=1;
		if(iline <= 8){ //Accidentals
			for(int iBin=0; iBin<numBins; iBin++){
				h_acc[2*(hall-1)+det-1]->SetBinContent(iBin+1, bg_values[iBin]*BgCounts[2*(hall-1)+det-1][0]);
				h_acc[2*(hall-1)+det-1]->SetBinError(iBin+1, bg_values[iBin]*BgError[2*(hall-1)+det-1][0]);
			}
			if(h_acc[2*(hall-1)+det-1]->Integral() != BgCounts[2*(hall-1)+det-1][0]) cout << "ACCIDENTALS FOR EH" << hall << "AD" << det << " ARE OFF BY " << (h_acc[2*(hall-1)+det-1]->Integral()-BgCounts[2*(hall-1)+det-1][0])/(h_acc[2*(hall-1)+det-1]->Integral()) << " %" << endl;
		}
		else if(iline == 9){ //Amc
			for(int iad=0; iad<maxAD; iad++){
				for(int iBin=0; iBin<numBins; iBin++){
					h_amc[iad]->SetBinContent(iBin+1, bg_values[iBin]*BgCounts[iad][1]);
					h_amc[iad]->SetBinError(iBin+1, bg_values[iBin]*BgError[iad][1]);
				}
				if(h_amc[iad]->Integral() != BgCounts[iad][1]) cout << "AMC FOR EH" << EH[iad] << "AD" << AD[iad] << " ARE OFF BY " << (h_amc[iad]->Integral()-BgCounts[iad][1])/(h_amc[iad]->Integral()) << " %" << endl;
			}
		}
		else if(iline == 10){ //Fast-n
			for(int iad=0; iad<maxAD; iad++){
				for(int iBin=0; iBin<numBins; iBin++){
					h_fastn[iad]->SetBinContent(iBin+1, bg_values[iBin]*BgCounts[iad][2]);
					h_fastn[iad]->SetBinError(iBin+1, bg_values[iBin]*BgError[iad][2]);
				}
				if(h_fastn[iad]->Integral() != BgCounts[iad][2]) cout << "FAST N FOR EH" << EH[iad] << "AD" << AD[iad] << " ARE OFF BY " << (h_fastn[iad]->Integral()-BgCounts[iad][2])/(h_fastn[iad]->Integral()) << " %" << endl;
			}
		}
		else if(iline == 11){ //Li9
			for(int iad=0; iad<maxAD; iad++){
				for(int iBin=0; iBin<numBins; iBin++){
					h_li9[iad]->SetBinContent(iBin+1, bg_values[iBin]*BgCounts[iad][3]);
					h_li9[iad]->SetBinError(iBin+1, bg_values[iBin]*BgError[iad][3]);
				}
				if(h_li9[iad]->Integral() != BgCounts[iad][3]) cout << "LI9 FOR EH" << EH[iad] << "AD" << AD[iad] << " ARE OFF BY " << (h_li9[iad]->Integral()-BgCounts[iad][3])/(h_li9[iad]->Integral()) << " %" << endl;
			}
		}
		else if(iline == 12){ //Rad-n
			for(int iad=0; iad<maxAD; iad++){
				for(int iBin=0; iBin<numBins; iBin++){
					h_radn[iad]->SetBinContent(iBin+1, bg_values[iBin]*BgCounts[iad][4]);
					h_radn[iad]->SetBinError(iBin+1, bg_values[iBin]*BgError[iad][4]);
				}
				if(h_radn[iad]->Integral() != BgCounts[iad][4]) cout << "RAD N FOR EH" << EH[iad] << "AD" << AD[iad] << " ARE OFF BY " << (h_radn[iad]->Integral()-BgCounts[iad][4])/(h_radn[iad]->Integral()) << " %" << endl;
			}
		}
	    }
	    bkgd_file.close();
	  }
	  else cout << "Unable to open file" << endl; //End of Background Spectra
	


	//Writing to the file
	TFile* outfile=new TFile(db_hists, "RECREATE");
	outfile->cd();
	for(int iad=0; iad<maxAD; iad++){
		h_prompt[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_prompt[iad]->GetYaxis()->SetTitle("Counts");
		h_prompt[iad]->Write();

		h_acc[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_acc[iad]->GetYaxis()->SetTitle("Counts");
		h_acc[iad]->Write();

		h_amc[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_amc[iad]->GetYaxis()->SetTitle("Counts");
		h_amc[iad]->Write();

		h_fastn[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_fastn[iad]->GetYaxis()->SetTitle("Counts");
		h_fastn[iad]->Write();

		h_li9[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_li9[iad]->GetYaxis()->SetTitle("Counts");
		h_li9[iad]->Write();

		h_radn[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_radn[iad]->GetYaxis()->SetTitle("Counts");
		h_radn[iad]->Write();
}

	for(int iad=0; iad<maxAD; iad++){
		h_ibds[iad]->Add(h_prompt[iad]);
		h_ibds[iad]->Add(h_acc[iad],-1);
		h_ibds[iad]->Add(h_amc[iad],-1);
		h_ibds[iad]->Add(h_fastn[iad],-1);
		h_ibds[iad]->Add(h_li9[iad],-1);
		h_ibds[iad]->Add(h_radn[iad],-1);
		h_ibds[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ibds[iad]->GetYaxis()->SetTitle("Counts");
		h_ibds[iad]->Write();
	}

}


void allDB_check(string config_file){
	
	init_config(config_file);
	
	prompt_spectra(config_file);
	bkgd_counts(config_file);
	bkgd_spectra(config_file);

	readOutput(config_file);
}


