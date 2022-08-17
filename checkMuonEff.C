#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TMinuit.h"
#include <typeinfo>

using namespace std;

void checkMuonEff(int hall_num, int period){

//	char hallList[64];
//	sprintf(hallList, "run_list_good_sync.txt",hall_num);
	FILE* runfile=fopen("./run_list_good_sync.txt","r");

	int maxAD = 0;
	int nRuns = 0;
//	int numPeriods=3;
	if(hall_num == 1){
		maxAD = 2;
		nRuns = 1134;
	}
	if(hall_num == 2){
		maxAD = 2;
		nRuns = 1125;
	}
	if(hall_num == 3){
		maxAD = 4;
		nRuns = 1215;
	}

	int run_num = 0;
	int EH = 0;

	//Go through the list and count the number of entries:
/*	while(1){
		fscanf(runfile,"%d %d",&run_num,&EH);
		if(feof(runfile)) break; //If it's the end of the file, break.
		if(EH == hall_num){
			nRuns += 1;
		}
		else if(EH != hall_num){
			continue;
		}
	}
	cout << "Number of runs for EH" << hall_num << " is: " << nRuns << endl;
*/
	//Making Histograms:

		//Efficiency Plot
/*		TH1F* h_muon_efficiency[maxAD]; //efficiency histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_muon_efficiency_ad%d", iad+1);
			h_muon_efficiency[iad]=new TH1F(name,name,nRuns,0.5,nRuns+0.5);
		}
*/

	int run_order = 0;
	double total_DAQ_ad1 = 0;
	double total_DAQ_ad2 = 0;
	double total_DAQ_ad3 = 0;
	double total_DAQ_ad4 = 0;
	double tot_live_ad1 = 0;
	double tot_veto_ad1 = 0;
	double tot_live_ad2 = 0;
	double tot_veto_ad2 = 0;
	double tot_live_ad3 = 0;
	double tot_veto_ad3 = 0;
	double tot_live_ad4 = 0;
	double tot_veto_ad4 = 0;
	double efficiency[maxAD];
	double TOTAL_DAQ[maxAD];
	double TOTAL_LIVE[maxAD];
//	double TOTAL_DAQ_PERIOD[maxAD][numPeriods];
//	double TOTAL_LIVE_PERIOD[maxAD][numPeriods];
	double AVG_EFF[maxAD];
	int ad_on[maxAD];
	int ad_start[maxAD];
	int ad_end[maxAD];
	const int start_8ad = 34523;
	const int end_8ad = 67012;

	for(int iad=0; iad<maxAD; iad++){
		efficiency[iad] = 0;
		AVG_EFF[iad] = 0;
		TOTAL_DAQ[iad] = 0;
		TOTAL_LIVE[iad] = 0;
		ad_on[iad]=0;
		ad_start[iad]=0;
		ad_end[iad]=0;
//		for(int iperiod=0; iperiod<numPeriods; iperiod++){
//			TOTAL_DAQ_PERIOD[iad][iperiod]=0;
//			TOTAL_LIVE_PERIOD[iad][iperiod]=0;
//		}
	}

	while(1){
		fscanf(runfile,"%d %d",&run_num,&EH);
		if(feof(runfile)) break; //If it's the end of the file, break.
		if(EH != hall_num){
			continue;
		}
		
		if(period==6 && run_num >= start_8ad) break;
		else if(period==7 && run_num <= end_8ad) continue;
		else if(period==8){
			if(run_num < start_8ad) continue;
			if(run_num > end_8ad) break;
		}
		
		char runFileName[64];
		sprintf(runFileName, "./accResults/EH%d/AccidentalsPlots_NU_1500_%d.root",hall_num, run_num);
		TFile *runFile = new TFile(runFileName);
		
		//TTree tr_sum;
		TTree *tr_efficiency = (TTree*)runFile->Get("tr_accidentals_time;1");
		tr_efficiency->SetBranchStatus("*",1);
		tr_efficiency->SetBranchAddress("p_DAQ_ad1",&total_DAQ_ad1);
		tr_efficiency->SetBranchAddress("p_DAQ_ad2",&total_DAQ_ad2);
		tr_efficiency->SetBranchAddress("p_DAQ_ad3",&total_DAQ_ad3);
		tr_efficiency->SetBranchAddress("p_DAQ_ad4",&total_DAQ_ad4);
		tr_efficiency->SetBranchAddress("p_live_ad1",&tot_live_ad1);
		tr_efficiency->SetBranchAddress("p_veto_ad1",&tot_veto_ad1);
		tr_efficiency->SetBranchAddress("p_live_ad2",&tot_live_ad2);
		tr_efficiency->SetBranchAddress("p_veto_ad2",&tot_veto_ad2);
			tr_efficiency->SetBranchAddress("p_live_ad3",&tot_live_ad3);
			tr_efficiency->SetBranchAddress("p_veto_ad3",&tot_veto_ad3);
			tr_efficiency->SetBranchAddress("p_live_ad4",&tot_live_ad4);
			tr_efficiency->SetBranchAddress("p_veto_ad4",&tot_veto_ad4);

		tr_efficiency->GetEntry(0);

		cout << "DAQ time for run" << run_num << " is: " << total_DAQ_ad1 << "s" << endl;

		if(total_DAQ_ad1 == 0 && total_DAQ_ad2 == 0){
			cout << "DAQ time = 0...?!?!?!" << endl;
			run_order += 1;
			runFile->Close();
			continue;
		}

		if(ad_on[0]==0 && total_DAQ_ad1 !=0){ //For AD1 to turn on
			ad_on[0]=1;
			ad_start[0]=run_num;
		}

		if(ad_on[1]==0 && total_DAQ_ad2 !=0){ //For AD2 to turn on
			ad_on[1]=1;
			ad_start[1]=run_num;
		}

		if(ad_on[2]==0 && total_DAQ_ad3 !=0){ //For AD3 to turn on
			ad_on[2]=1;
			ad_start[2]=run_num;
		}

		if(ad_on[3]==0 && total_DAQ_ad4 !=0){ //For AD4 to turn on
			ad_on[3]=1;
			ad_start[3]=run_num;
		}

			if(ad_on[0]==1 && total_DAQ_ad1 ==0){ //For AD1 to turn off
				ad_on[0]=0;
			}

			if(ad_on[1]==1 && total_DAQ_ad2 ==0){ //For AD2 to turn off
				ad_on[1]=0;
			}

			if(ad_on[2]==1 && total_DAQ_ad3 ==0){ //For AD3 to turn off
				ad_on[2]=0;
			}

			if(ad_on[3]==1 && total_DAQ_ad4 ==0){ //For AD4 to turn off
				ad_on[3]=0;
			}
		
		if(ad_on[0]==1){ //For AD1 to update the end run number
			ad_end[0]=run_num;
		}

		if(ad_on[1]==1){ //For AD2 to update the end run number
			ad_end[1]=run_num;
		}

		if(ad_on[2]==1){ //For AD3 to update the end run number
			ad_end[2]=run_num;
		}

		if(ad_on[3]==1){ //For AD4 to update the end run number
			ad_end[3]=run_num;
		}

		if(total_DAQ_ad1 != 0) efficiency[0] = tot_live_ad1/total_DAQ_ad1;
		if(total_DAQ_ad2 != 0) efficiency[1] = tot_live_ad2/total_DAQ_ad2;
		if(total_DAQ_ad3 != 0) efficiency[2] = tot_live_ad3/total_DAQ_ad3;
		if(total_DAQ_ad4 != 0) efficiency[3] = tot_live_ad4/total_DAQ_ad4;

		for(int iad=0; iad<maxAD; iad++){
			if(hall_num == 2 && efficiency[1] == efficiency[2] && iad == 1) continue;
//			h_muon_efficiency[iad]->Fill(run_order+1, efficiency[iad]);
			AVG_EFF[iad] = (AVG_EFF[iad]*(run_order+1) + efficiency[iad])/(run_order+2);
		}

		TOTAL_DAQ[0] += total_DAQ_ad1;
		TOTAL_DAQ[1] += total_DAQ_ad2;
		TOTAL_DAQ[2] += total_DAQ_ad3;
		TOTAL_DAQ[3] += total_DAQ_ad4;

		TOTAL_LIVE[0] += tot_live_ad1;
		TOTAL_LIVE[1] += tot_live_ad2;
		TOTAL_LIVE[2] += tot_live_ad3;
		TOTAL_LIVE[3] += tot_live_ad4;

		cout << "Done with run_order #" << run_order << " out of " << nRuns << "runs." << endl;


		run_order += 1;
		runFile->Close();
	}


  /*      char outputname[64];
        sprintf(outputname,"./MuonEff_sync_EH%d.root",hall_num);
	TFile* outfile=new TFile(outputname, "RECREATE");
		outfile->cd();
		for(int iad=0; iad<maxAD; ++iad){
			h_muon_efficiency[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_muon_efficiency[iad]->GetYaxis()->SetTitle("Muon Efficiency");
			h_muon_efficiency[iad]->Write();
		}
*/

cout << fixed<< setprecision(10);

	for(int iad=0; iad<maxAD; iad++){
		cout << "Total DAQ time (days) for EH" << hall_num << " AD" << iad+1 << " is: " << TOTAL_DAQ[iad]/(24.*60.*60.) << endl;
		cout << "Total Livetime (days) for EH" << hall_num << " AD" << iad+1 << " is: " << TOTAL_LIVE[iad]/(24.*60.*60.) << endl;
		cout << "Total Vetotime (days) for EH" << hall_num << " AD" << iad+1 << " is: " << (TOTAL_DAQ[iad]-TOTAL_LIVE[iad])/(24.*60.*60.) << endl;
	//	cout << "Average (over runs) Efficiency for EH" << hall_num << " AD" << iad+1 << " is: " << AVG_EFF[iad] << endl << endl;
		AVG_EFF[iad] = TOTAL_LIVE[iad]/TOTAL_DAQ[iad];
		cout << "Average (over time) Efficiency for EH" << hall_num << " AD" << iad+1 << " is: " << TOTAL_LIVE[iad]/TOTAL_DAQ[iad] << endl << endl;

		cout << "Total DAQ from livetime/efficency (days) for EH" << hall_num << " AD" << iad+1 << " is: " << TOTAL_LIVE[iad]/(AVG_EFF[iad]*24.*60.*60.) << endl;
		cout << "Total Livetime from livetime (days) for EH" << hall_num << " AD" << iad+1 << " is: " << TOTAL_LIVE[iad]/(24.*60.*60.) << endl;
		cout << "Total Livetime from DAQ*efficency (days) for EH" << hall_num << " AD" << iad+1 << " is: " << TOTAL_DAQ[iad]*AVG_EFF[iad]/(24.*60.*60.) << endl;
		cout << "Total Livetime (ns) for EH" << hall_num << " AD" << iad+1 << " is: " << TOTAL_LIVE[iad]*1.e9 << endl << endl << endl;
	}
	
	for(int iad=0; iad<maxAD; iad++){
		cout << "START RUN FOR AD " << iad+1 << ": " << ad_start[iad] << "\tEND RUN FOR AD " << iad+1 << ": " << ad_end[iad] << endl;
	}

//	outfile->Close();

}
