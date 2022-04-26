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
#include "TH2F.h"
#include <typeinfo>
#include <TRandom3.h>
#include <iomanip>

using namespace std;

char outputname_findPrompt[64];
char outputname_findDelayed[64];
char outputname_summarize[64];
char outputname_summarize_long1[64];
char outputname_summarize_long2[64];
char outputname_summarize_long3[64];
char outputname_plot[64];

//varibales for the phi-z(-r) and three r2-z correction maps
TH2F* h_phiz_nu_correction_map[8][2]; //8ADs and 2 time periods
TH2F* h_r2z_nu_correction_map_alphas_only[8][2]; //8 ADs and 2 time periods
TH2F* h_r2z_nu_correction_map_alphas_ngd[8][2]; //8 ADs and 2 time periods //a bit dull coding but easier to understand the code
TH2F* h_r2z_nu_correction_map_nh[8][2]; //8 ADs and 2 time periods //a bit dull coding but easier to understand the code

//points of division between pre and post dead PMT period
float period_division_sec[8] = {1490889600., 1519747200., 1427731200., 1483113600., 1459353600., 1519747200., 1522425600., 1443542400.};
string period_division_date[8]={"20161220","20180228","20150531","20161231","20170331","20180228","20180331","20150730"};

//constants to remove the "original" phi-correction
float amplitude[8] = {-0.00348262, -0.00657182, 0.00674897, 0.00423918, 0.00427006, 0.0111522, 0.00670342, 0.00997318};
float phase[8] = {0.727572, 0.671381, 0.0868127, -0.129778, -1.18759, -1.05132, -0.610239, -0.829069};

void LoadMaps(){ //load the correction maps from the file
    int EH[8]={1,1,2,2,3,3,3,3};
    int AD[8]={1,2,1,2,1,2,3,4};
    TFile infile("./non_uniformity_correction_maps.root","READ");
    for(int iad=0; iad<8; ++iad){
        h_phiz_nu_correction_map[iad][0]=(TH2F*)infile.Get(Form("h_nu_corr_map_alpha_only_phiz_wholeAD_eh%d_ad%d_time0",EH[iad],AD[iad]));
        h_phiz_nu_correction_map[iad][1]=(TH2F*)infile.Get(Form("h_nu_corr_map_alpha_only_phiz_wholeAD_eh%d_ad%d_time1",EH[iad],AD[iad]));
        h_phiz_nu_correction_map[iad][0]->SetDirectory(0);
        h_phiz_nu_correction_map[iad][1]->SetDirectory(0);
        
        h_r2z_nu_correction_map_alphas_only[iad][0]=(TH2F*)infile.Get(Form("h_nu_corr_map_alpha_only_r2z_eh%d_ad%d_before%s",EH[iad],AD[iad],period_division_date[iad].c_str()));
        h_r2z_nu_correction_map_alphas_only[iad][1]=(TH2F*)infile.Get(Form("h_nu_corr_map_alpha_only_r2z_eh%d_ad%d_after%s",EH[iad],AD[iad],period_division_date[iad].c_str()));
        h_r2z_nu_correction_map_alphas_only[iad][0]->SetDirectory(0);
        h_r2z_nu_correction_map_alphas_only[iad][1]->SetDirectory(0);
        
        h_r2z_nu_correction_map_alphas_ngd[iad][0]=(TH2F*)infile.Get(Form("h_nu_corr_map_alpha_ngd_r2z_eh%d_ad%d_before%s",EH[iad],AD[iad],period_division_date[iad].c_str()));
        h_r2z_nu_correction_map_alphas_ngd[iad][1]=(TH2F*)infile.Get(Form("h_nu_corr_map_alpha_ngd_r2z_eh%d_ad%d_after%s",EH[iad],AD[iad],period_division_date[iad].c_str()));
        h_r2z_nu_correction_map_alphas_ngd[iad][0]->SetDirectory(0);
        h_r2z_nu_correction_map_alphas_ngd[iad][1]->SetDirectory(0);
        
        h_r2z_nu_correction_map_nh[iad][0]=(TH2F*)infile.Get(Form("h_nu_corr_map_nh_r2z_eh%d_ad%d_before%s",EH[iad],AD[iad],period_division_date[iad].c_str()));
        h_r2z_nu_correction_map_nh[iad][1]=(TH2F*)infile.Get(Form("h_nu_corr_map_nh_r2z_eh%d_ad%d_after%s",EH[iad],AD[iad],period_division_date[iad].c_str()));
        h_r2z_nu_correction_map_nh[iad][0]->SetDirectory(0);
        h_r2z_nu_correction_map_nh[iad][1]->SetDirectory(0);
    }
    infile.Close();
}

float CorrectedEnergy(int eh, int detector, float energy, float x, float y, float z, int time_sec, bool use_alphas_only=false, bool use_nh=false){
    //eh number 1,2 or 3
    //ad number 1,2,3 or 4
    //energy - energy in AdSimpleNL TTree
    // x, y, z - in mm
    //time_sec - time of the trigger in sec (enough to determine the period)
    
    float corrected_energy = energy; //no correction yet
    
    int iad=2*(eh-1)+detector-1; //give you the proper iad number, EH1-AD1 iad=0, EH1-AD2 iad=1, EH2-AD1 iad=2, EH2-AD2 iad=3, EH3-AD1 iad=4, EH3-AD2 iad=5, EH3-AD3 iad=6, EH3-AD4 iad=7
    
    //detemination of the pre or post dead PMT period
    int period=-1;
    if(time_sec<=period_division_sec[iad]) period=0;
    else period=1;
    
    float phi_rad = atan2(y, x);
    float original_anticorrection_factor = 1.+ amplitude[iad]*sin(phi_rad + phase[iad]);
    
    //removing "original" phi-correction
    corrected_energy*=original_anticorrection_factor;
    
    //eventually loading the input non-uniformity maps
    if(h_phiz_nu_correction_map[iad][period]==NULL || h_r2z_nu_correction_map_nh[iad][period]==NULL || h_r2z_nu_correction_map_alphas_only[iad][period]==NULL || h_r2z_nu_correction_map_alphas_ngd[iad][period]==NULL) LoadMaps();
    
    float phi_deg = 180.*phi_rad/3.1415926;  //fyi, original correction used radians and phi in (-pi, pi), new correction uses degress and phi in (0,360)
    if (phi_deg < 0) phi_deg += 360.;
    
    int phi_bin = int(phi_deg/30.) + 1;
    if (phi_bin > 12) phi_bin = 12;
    
    int z_bin_for_phiz = int((z+2000.)/500.) + 1;
    if(z_bin_for_phiz<1) z_bin_for_phiz=1; //events reconstructed below the map are corrected as being at the edge in z, i.e. -2 m,
    if(z_bin_for_phiz>8) z_bin_for_phiz=8; //events reconstructed above the map are corrected as being at the edge in z, i.e. 2 m
    
    float r = sqrt(x*x + y*y)/1.E3;//in m, NOT mm
    
    //application of the new phi-z(-r) correction
    float correction_factor_phiz = 1. - r*h_phiz_nu_correction_map[iad][period]->GetBinContent(phi_bin,z_bin_for_phiz); //phi-z correction based on map and linearly on radius
    corrected_energy*=correction_factor_phiz;
    
    int r2_bin=int((x*x+y*y)/0.4E6) + 1; //hard-coded nu correction map dimensions - r^2 has 0.4 m^2 = 0.4E6 mm^2 bins, + 1 is due to ROOT bin numbering convention
    if(r2_bin>10) r2_bin=10; //events outside the map are corrected as being at the edge in r^2, i.e. 4 m^2
    
    int z_bin=int((z+2000.)/400.) + 1; //hard-coded nu correction map dimensions - z has 0.4 m = 400 mm bins, + 1 is due to ROOT bin numbering convention
    if(z_bin<1) z_bin=1; //events reconstucted below the map are corrected as being at the edge in z, i.e. -2 m,
    if(z_bin>10) z_bin=10; //events reconstucted above the map are corrected as being at the edge in z, i.e. 2 m
    
    if(use_alphas_only == true && use_nh==true) return -1; //if this happens, you did not select the r2-z map correctly
    
    //retriving the r2-z correction factor from a proper map - by defaults, nGd+alpha map is used
    float correction_factor_r2z=-1.;
    if(use_alphas_only) correction_factor_r2z=h_r2z_nu_correction_map_alphas_only[iad][period]->GetBinContent(r2_bin,z_bin);
    else if(use_nh) correction_factor_r2z=h_r2z_nu_correction_map_nh[iad][period]->GetBinContent(r2_bin,z_bin);
    else correction_factor_r2z=h_r2z_nu_correction_map_alphas_ngd[iad][period]->GetBinContent(r2_bin,z_bin);
    
    //applying r2-z correction
    corrected_energy*=correction_factor_r2z;
    
    //returns final corrected energy
    return corrected_energy;
}


bool IsFlasher(float m_MaxQ, float m_Quadrant, float m_MaxQ_2inchPMT, float time_PSD, float time_PSD1, float Q1, float Q2, float x, float y, float z){ //function that decides if event is light from flashing PMT (This is mostly for low energy events, not muons)
    if(m_MaxQ_2inchPMT>100.) return true;   //2-inch PMT cut
    if((pow(m_Quadrant,2)+pow(m_MaxQ/0.45,2))>1) return true; //ellipse cut
	if((4.*pow(1.-time_PSD,2)+1.8*pow(1.-time_PSD1,2))>1){
//		cout << "PSD cut" << endl;
		 return true; //PSD cut
	}
	if((Q1/Q2 > 0.6) && (log10(pow(m_Quadrant,2)+pow(m_MaxQ/0.45,2)) > 0.5*Q1/Q2 - 0.8) && (log10(pow(m_Quadrant,2)+pow(m_MaxQ/0.45,2)) > -0.3)){
//		cout << "cluster cut" << endl;
		 return true; //residual "cluster" flasher cut
	}
	if((z>2.4 *1.e3) && (pow(x,2)+pow(y,2) > 0.5*1.e3)){
//		cout << "top cut" << endl;
		 return true; //top ring flasher cut
	}
	if(sqrt(pow(x,2)+pow(y,2))>2.2*1.e3){
//		cout << "large r cut" << endl;
		 return true; //large r flasher cut
	}

    else return false;
}


void find(int run_order, float singleElow, float singleEhigh, int timePreSingle, int timePostSingle, int pd_window_microsec){
	TChain* reconT_AdSimple=new TChain("/Event/Rec/AdSimpleNL");  //TChain - is like TTree but multiple files can be added  //reconT_AdSimple- containes reconstucted data
	TChain* calibStatsT=new TChain("/Event/Data/CalibStats"); //contains some other event characteristics
    
	int run_num=-1;
	int EH=-1;

	FILE* runfile=fopen("./run_list_good_sync.txt","r");  //This file contains all good runs numbers and experimental hall they belong to
	int iline=0;
	while(1){ //go though the file and retrieve run_num (run number) which corresponds to run_order provided as a parameter of GetMuons function
		fscanf(runfile,"%d %d",&run_num,&EH);
		if(feof(runfile)) break;
		if(iline==run_order) break;
		iline++;
	}
	fclose(runfile);
    
	cout<<"Run identified as: run#"<<run_num<<" EH"<<EH<<endl;

	ifstream inputfile("/dybfs/rec/P17B/GoodRun_v3/paths.physics.good.p17b.v3.sync.txt"); //this file contains paths to all good run files - each run (run_num) is devided to several .root files with TTrees
	std::string line;
	int countfiles = 0;
	while(1){  //reading the data from file
		std::getline(inputfile,line);
		if(inputfile.eof()) break;
		std::size_t found = line.find("recon.Neutrino.00");
		std::size_t found2 = line.find("recon.NoTag.00");
        
		if(found>600) found=found2-3;
        
		std::string st_runnum = line.substr(found+17,found+22);
		const char* tempchar=st_runnum.c_str();
		int temp_runnum=atoi(tempchar);
        
		if(temp_runnum!=run_num) continue;
		if(temp_runnum>run_num) break;
        
		const char* file_path=line.c_str();
        
		cout<<file_path<<endl;
        
		//above some tedious procedure to select proper lines with paths corresponding to our run_num
		reconT_AdSimple->AddFile(file_path); //finally adding those to TTrees
		calibStatsT->AddFile(file_path);

//		countfiles += 1;
//		if(countfiles == 5) break;
	}
	int nentries=reconT_AdSimple->GetEntries();
	cout<<"There are "<<nentries<<" entries"<<endl;
    
	reconT_AdSimple->SetBranchStatus("*",0);  //turning off reading all TBranches of TTrees to save data transfer
	calibStatsT->SetBranchStatus("*",0);
    
    //turning on just ones we will use
	reconT_AdSimple->SetBranchStatus("detector",1); //detector number 1,2,3,4 - AD1,AD2,AD3,AD4 in given experimental hall (EH) to which run belong to. AD3 and AD4 only for EH3. detector number 5,6 - IWS, OWS. detector number 7 - RPC
	reconT_AdSimple->SetBranchStatus("energy",1); //reconstucted energy - non-zero only for ADs since we do not reconstruct energy in other detectors //in MeV
	reconT_AdSimple->SetBranchStatus("triggerTimeSec",1);  //time is seconds from 1.1.1970 (I think)
	reconT_AdSimple->SetBranchStatus("triggerTimeNanoSec",1); //time is nanoseconds from beginning of last second
	reconT_AdSimple->SetBranchStatus("x",1);
	reconT_AdSimple->SetBranchStatus("y",1);
	reconT_AdSimple->SetBranchStatus("z",1);

	reconT_AdSimple->SetBranchStatus("triggerType",1);
    
	calibStatsT->SetBranchStatus("nHit",1); //number of PMTs hit - useful for IWS and OWS muon veto
	calibStatsT->SetBranchStatus("Quadrant",1);  //some characteristics for flashing PMT veto
	calibStatsT->SetBranchStatus("MaxQ",1);
	calibStatsT->SetBranchStatus("MaxQ_2inchPMT",1);
	calibStatsT->SetBranchStatus("time_PSD",1);
	calibStatsT->SetBranchStatus("time_PSD1",1);
	calibStatsT->SetBranchStatus("Q1",1);
	calibStatsT->SetBranchStatus("Q2",1);
    
	reconT_AdSimple->SetMakeClass(1);
	calibStatsT->SetMakeClass(1);
    
	int detector=0;
	int time_sec=0;
	int time_nanosec=0;
	int nhit=0;  //detector, times and nhit are integers
	float energy=0;
	float x=0;
	float y=0;
	float z=0;
	int trigType=0;
    
	reconT_AdSimple->SetBranchAddress("detector",&detector);  //setting where to save them if getting entry from TTree (TChain)
	reconT_AdSimple->SetBranchAddress("energy",&energy);
	reconT_AdSimple->SetBranchAddress("triggerTimeSec",&time_sec);
	reconT_AdSimple->SetBranchAddress("triggerTimeNanoSec",&time_nanosec);
	reconT_AdSimple->SetBranchAddress("x",&x);
	reconT_AdSimple->SetBranchAddress("y",&y);
	reconT_AdSimple->SetBranchAddress("z",&z);
    	reconT_AdSimple->SetBranchAddress("triggerType",&trigType);
    
	float Quadrant=0;
	float MaxQ=0;
	float MaxQ_2inchPMT=0;
	float time_PSD=0;
	float time_PSD1=0;
	float Q1=0;
	float Q2=0;
    
	calibStatsT->SetBranchAddress("Quadrant",&Quadrant);
	calibStatsT->SetBranchAddress("MaxQ",&MaxQ);
	calibStatsT->SetBranchAddress("MaxQ_2inchPMT",&MaxQ_2inchPMT);
	calibStatsT->SetBranchAddress("nHit",&nhit);
	calibStatsT->SetBranchAddress("time_PSD",&time_PSD);
	calibStatsT->SetBranchAddress("time_PSD1",&time_PSD1);
	calibStatsT->SetBranchAddress("Q1",&Q1);
	calibStatsT->SetBranchAddress("Q2",&Q2);

	//TTree Variables
	int check_ad[4];
	float s_energy=0;
	int s_sec = 0;
	int s_nanosec = 0;
	float s_x = 0;
	float s_y = 0;
	float s_z = 0;
	float distance = 0;
	int hall_num = 0;
	int det_num = 0;
	double time_between = 0;
	
	int window_start_sec = 0;
	int window_start_nanosec = 0;
	int window_end_sec = 0;
	int window_end_nanosec = 0;
	float window_length = 0;
	
	double total_DAQ_ad1 = 0;
	double total_DAQ_ad2 = 0;
	double total_DAQ_ad3 = 0;
	double total_DAQ_ad4 = 0;
	double tot_veto_ad1 = 0;
	double tot_live_ad1 = 0;
	double tot_veto_ad2 = 0;
	double tot_live_ad2 = 0;
	double tot_veto_ad3 = 0;
	double tot_live_ad3 = 0;
	double tot_veto_ad4 = 0;
	double tot_live_ad4 = 0;
	
	

	//Variables - DAQ calculation
	int veto_time_ad_muon_microsec=800;
	int veto_time_shower_muon_sec=1;
	int veto_time_ws_muon_microsec=400;

	int prev_trig_sec = 0; //For total DAQ
	int prev_trig_nanosec = 0; //For total DAQ

	int window_time_sec = 0.;
	int window_time_nanosec = 0;
	double total_DAQ_time = 0.;
	int trig_gap_sec = 0;
	int trig_gap_nanosec = 0;
	int muon_start_sec[4];
	int muon_start_nanosec[4];
	int muon_end_sec[4];
	int muon_end_nanosec[4];
	int muon_veto_sec[4];
	int muon_veto_nanosec[4];
	double total_veto_time[4];
	int proj_end_sec[4];
	int proj_end_nanosec[4];
	int pre_veto_nanosec = timePostSingle*1000; //additional vetoed time because delayed event can't be that close
	int post_veto_nanosec = timePreSingle*1000; //additional vetoed time because delayed event can't be that close
	int muon_veto_start_sec = 0;
	int muon_veto_start_nanosec = 0;
	int muon_veto_end_sec = 0;
	int muon_veto_end_nanosec = 0;
	int first_trig_sec = 0;
	int first_trig_nanosec = 0;
	int last_trig_sec = 0;
	int last_trig_nanosec = 0;
	
	//Initializing:
	for(int i=0; i<4; i++){
		muon_start_sec[i]=0;
		muon_start_nanosec[i]=0;
		muon_end_sec[i]=0;
		muon_end_nanosec[i]=0;	
		muon_veto_sec[i]=0;
		muon_veto_nanosec[i]=0;
		total_veto_time[i]=0;
		proj_end_sec[i]=0;
		proj_end_nanosec[i]=0;
		check_ad[i] = 0;
	}
	
	//Variables for IBD seletion
	int last_ws_muon_sec = 0;
	int last_ws_muon_nanosec = 0;
	int last_shower_muon_sec[4];
	int last_shower_muon_nanosec[4];
	int last_ad_muon_sec[4];
	int last_ad_muon_nanosec[4];
//	float dElow = 1.75;
//	float dEhigh = 2.89;
	int s_entry[4];
	int single_det = 0;
	int single_sec[4];
	int single_nanosec[4];
	float single_energy[4];
	float single_x[4];
	float single_y[4];
	float single_z[4];
	int pre_counter = 0;
	float sElow = 1.5; //For checking that no other events are in the isolation times for the singles
	float sEhigh = 12.; //For checking that no other events are in the isolation times for the singles
	int in_pre = 0;
	int post_counter = 0;
	int incomplete_post_window = 0;
	int temp_prev_trig_sec = 0;
	int temp_prev_trig_nanosec = 0;
	
	//Initializing
	for(int i=0; i<4; i++){
		last_shower_muon_sec[i] = 0;
		last_shower_muon_nanosec[i] = 0;
		last_ad_muon_sec[i] = 0;
		last_ad_muon_nanosec[i] = 0;
		s_entry[i] = 0;
		single_sec[i] = 0;
		single_nanosec[i] = 0;
		single_energy[i] = 0;
		single_x[i] = 0;
		single_y[i] = 0;
		single_z[i] = 0;
	}

	
    
    //----Here you create some object you want to save things to

	//IBD TTree:
	TTree* tr_sing=new TTree("tr_sing", "My Singles Tree");
	tr_sing->Branch("s_energy", &s_energy, "Energy of single event/Float");
	tr_sing->Branch("s_sec", &s_sec, "Time (sec) of single event/Int");
	tr_sing->Branch("s_nanosec", &s_nanosec, "Time (sec) of single event/Int");
	tr_sing->Branch("s_x", &s_x, "X position of single event/Float");
	tr_sing->Branch("s_y", &s_y, "Y position of single event/Float");
	tr_sing->Branch("s_z", &s_z, "Z position of single event/Float");
	tr_sing->Branch("hall_num", &hall_num, "Hall Number/Int");
	tr_sing->Branch("det_num", &det_num, "Detector Number/Int");


	//DAQ TTree:
	TTree* tr_DAQ = new TTree("tr_DAQ", "My DAQ Tree");
	tr_DAQ->Branch("window_start_sec", &window_start_sec, "Sec of window starting/Int");
	tr_DAQ->Branch("window_start_nanosec", &window_start_nanosec, "Nanosec of window starting/Int");
	tr_DAQ->Branch("window_end_sec", &window_end_sec, "Sec of window ending/Int");
	tr_DAQ->Branch("window_end_nanosec", &window_end_nanosec, "Nanosec of window ending/Int");
	tr_DAQ->Branch("window_length", &window_length, "Length of window in seconds/Float");

	//Summary TTree
	TTree* tr_summary = new TTree("tr_summary", "Summary Tree");
	tr_summary->Branch("total_DAQ_ad1", &total_DAQ_ad1, "Total Data Acquisition Time for AD1 in s/Double");
	tr_summary->Branch("total_DAQ_ad2", &total_DAQ_ad2, "Total Data Acquisition Time for AD2 in s/Double");
	tr_summary->Branch("total_DAQ_ad3", &total_DAQ_ad3, "Total Data Acquisition Time for AD3 in s/Double");
	tr_summary->Branch("total_DAQ_ad4", &total_DAQ_ad4, "Total Data Acquisition Time for AD4 in s/Double");
	tr_summary->Branch("tot_veto_ad1", &tot_veto_ad1, "Total Veto Time for AD1 in s/Double");
	tr_summary->Branch("tot_live_ad1", &tot_live_ad1, "Total Live Time for AD1 in s/Double");
	tr_summary->Branch("tot_veto_ad2", &tot_veto_ad2, "Total Veto Time for AD2 in s/Double");
	tr_summary->Branch("tot_live_ad2", &tot_live_ad2, "Total Live Time for AD2 in s/Double");
	tr_summary->Branch("tot_veto_ad3", &tot_veto_ad3, "Total Veto Time for AD3 in s/Double");
	tr_summary->Branch("tot_live_ad3", &tot_live_ad3, "Total Live Time for AD3 in s/Double");
	tr_summary->Branch("tot_veto_ad4", &tot_veto_ad4, "Total Veto Time for AD4 in s/Double");
	tr_summary->Branch("tot_live_ad4", &tot_live_ad4, "Total Live Time for AD4 in s/Double");

	int first_entry = 0;

	cout << fixed;
	cout << setprecision(9);

	while(1){ //determining the first non-RPC entry
		reconT_AdSimple->GetEntry(first_entry);
		calibStatsT->GetEntry(first_entry);
		if(detector == 7) first_entry = first_entry+1;
		if(detector != 7){
			first_trig_sec = time_sec;
			first_trig_nanosec = time_nanosec;
			break;
		}
	}

	while(1){ //determining the last non-RPC entry
		reconT_AdSimple->GetEntry(nentries-1);
		calibStatsT->GetEntry(nentries-1);
		if(detector == 7) nentries = nentries-1;
		if(detector != 7){
			last_trig_sec = time_sec;
			last_trig_nanosec = time_nanosec;
			break;
		}
	}

	for(int ientry=first_entry; ientry<nentries; ientry++){//Going through the data
		//Tracking the progress:
		if(ientry%1000000==0) cout<<"Running "<<ientry/1.e6<<" M out of "<<nentries/1.e6<<" M"<<endl;

		reconT_AdSimple->GetEntry(ientry); //Get the event information
		calibStatsT->GetEntry(ientry);

		if(detector == 7) continue; //Ignore RPC events

		//For livetime:
		det_num = detector-1;
		if(det_num < 4 && check_ad[det_num] == 0){
			check_ad[det_num] = 1; //Seeing if the AD is on or not
			cout << "Event in AD: " << det_num << endl;
		}

		//First entry - Total DAQ calc
		if(ientry == first_entry){
			window_start_sec = time_sec;
			window_start_nanosec = time_nanosec;
			for(int iad = 0; iad<4; iad++){ //First event counts as the end of a muon b/c we don't know the rest of the info
				muon_start_sec[iad] = time_sec;
				muon_start_nanosec[iad] = time_nanosec;
				muon_end_sec[iad] = time_sec;
				muon_end_nanosec[iad] = time_nanosec + post_veto_nanosec;
			}
		}
		
		if(ientry != first_entry){
			trig_gap_sec = time_sec - prev_trig_sec;
			trig_gap_nanosec = time_nanosec - prev_trig_nanosec;
		}
		
		//gap larger than 1 sec-->end window, start new window
		if(trig_gap_sec + 1.e-9*trig_gap_nanosec > 1){
			window_end_sec = prev_trig_sec;
			window_end_nanosec = prev_trig_nanosec;
			window_time_sec = window_end_sec - window_start_sec;
			window_time_nanosec = window_end_nanosec - window_start_nanosec;
			window_length = window_time_sec + 1.e-9*window_time_nanosec;
			total_DAQ_time = total_DAQ_time + window_length;
			tr_DAQ->Fill();

			for(int iad = 0; iad<4; iad++){ //closing out muon vetoes
				if(muon_end_sec[iad]-window_end_sec > 1.e-9*(window_end_nanosec - muon_end_nanosec[iad]-pre_veto_nanosec)){ //if there is already a muon window open
					muon_end_sec[iad] = window_end_sec;
					muon_end_nanosec[iad] = window_end_nanosec;
				}
				else{ //if there is not a muon window open during the closing of the window
					//Save the previous muon veto
					muon_veto_sec[iad] = muon_end_sec[iad]-muon_start_sec[iad];
					muon_veto_nanosec[iad] = muon_end_nanosec[iad]-muon_start_nanosec[iad];
					total_veto_time[iad] = total_veto_time[iad] + muon_veto_sec[iad] + 1.e-9*muon_veto_nanosec[iad];
					muon_veto_start_sec = muon_start_sec[iad];
					muon_veto_start_nanosec = muon_start_nanosec[iad];
					muon_veto_end_sec = muon_end_sec[iad];
					muon_veto_end_nanosec = muon_end_nanosec[iad];

					//Save the veto of the end of the window
					muon_start_sec[iad] = window_end_sec;
					muon_start_nanosec[iad] = window_end_nanosec - pre_veto_nanosec;
					muon_end_sec[iad] = window_end_sec;
					muon_end_nanosec[iad] = window_end_nanosec;
				}
				muon_veto_sec[iad] = muon_end_sec[iad]-muon_start_sec[iad];
				muon_veto_nanosec[iad] = muon_end_nanosec[iad]-muon_start_nanosec[iad];
				total_veto_time[iad] = total_veto_time[iad] + muon_veto_sec[iad] + 1.e-9*muon_veto_nanosec[iad];
				muon_veto_start_sec = muon_start_sec[iad];
				muon_veto_start_nanosec = muon_start_nanosec[iad];
				muon_veto_end_sec = muon_end_sec[iad];
				muon_veto_end_nanosec = muon_end_nanosec[iad];
				
			    //New window start up:
				muon_start_sec[iad] = time_sec;
				muon_start_nanosec[iad] = time_nanosec;
				muon_end_sec[iad] = time_sec;
				muon_end_nanosec[iad] = time_nanosec + post_veto_nanosec;
			}

			window_start_sec = time_sec;
			window_start_nanosec = time_nanosec;

		}//End of time gap stuff

		//Last entry - Total DAQ calc
		if(ientry == nentries-1){ //need to close out windows and muon vetoes
			window_end_sec = time_sec;
			window_end_nanosec = time_nanosec;
			window_time_sec = window_end_sec - window_start_sec;
			window_time_nanosec = window_end_nanosec - window_start_nanosec;
			window_length = window_time_sec + window_time_nanosec*1.e-9;
			total_DAQ_time = total_DAQ_time + window_length;
			tr_DAQ->Fill();

			for(int iad = 0; iad<4; iad++){
				//if there is already a muon window open
				if(muon_end_sec[iad]-window_end_sec > 1.e-9*(window_end_nanosec - muon_end_nanosec[iad]-pre_veto_nanosec)){
					muon_end_sec[iad] = window_end_sec;
					muon_end_nanosec[iad] = window_end_nanosec;
				}
				else{ //if there is not a muon window open during the closing of the window
					//Save the previous muon veto
					muon_veto_sec[iad] = muon_end_sec[iad]-muon_start_sec[iad];
					muon_veto_nanosec[iad] = muon_end_nanosec[iad]-muon_start_nanosec[iad];
					total_veto_time[iad] = total_veto_time[iad] + muon_veto_sec[iad] + 1.e-9*muon_veto_nanosec[iad];
					muon_veto_start_sec = muon_start_sec[iad];
					muon_veto_start_nanosec = muon_start_nanosec[iad];
					muon_veto_end_sec = muon_end_sec[iad];
					muon_veto_end_nanosec = muon_end_nanosec[iad];

					//Starting new vetoed period right before the end.
					muon_start_sec[iad] = window_end_sec;
					muon_start_nanosec[iad] = window_end_nanosec - pre_veto_nanosec;
					muon_end_sec[iad] = window_end_sec;
					muon_end_nanosec[iad] = window_end_nanosec;			
				}

				muon_veto_sec[iad] = muon_end_sec[iad]-muon_start_sec[iad];
				muon_veto_nanosec[iad] = muon_end_nanosec[iad]-muon_start_nanosec[iad];
				total_veto_time[iad] = total_veto_time[iad] + muon_veto_sec[iad] + 1.e-9*muon_veto_nanosec[iad];
				muon_veto_start_sec = muon_start_sec[iad];
				muon_veto_start_nanosec = muon_start_nanosec[iad];
				muon_veto_end_sec = muon_end_sec[iad];
				muon_veto_end_nanosec = muon_end_nanosec[iad];

			}

		}//end of last entry conditions
		
		//WS muon - affects all ADs
		if(((detector==5 && nhit>12) || (detector==6 && nhit>15))&&(trigType==268439552 || trigType==268435712 || trigType==268439808)){ //Trig type specified to avoid cross trigger issues
			for(int iad=0; iad<4; iad++){
				if(time_sec-muon_end_sec[iad]>1.e-9*(muon_end_nanosec[iad]-time_nanosec+pre_veto_nanosec)){ //if after the last muon time has ended
					//Finishing last muon period
					muon_veto_start_sec = muon_start_sec[iad];
					muon_veto_start_nanosec = muon_start_nanosec[iad];
					muon_veto_end_sec = muon_end_sec[iad];
					muon_veto_end_nanosec = muon_end_nanosec[iad];
					muon_veto_sec[iad] = muon_end_sec[iad] - muon_start_sec[iad];
					muon_veto_nanosec[iad] = muon_end_nanosec[iad] - muon_start_nanosec[iad];
					total_veto_time[iad] = total_veto_time[iad] + muon_veto_sec[iad] + 1.e-9*muon_veto_nanosec[iad];
				
					//Starting new muon period
					muon_start_sec[iad] = time_sec;
					muon_start_nanosec[iad] = time_nanosec - pre_veto_nanosec;
					muon_end_sec[iad] = time_sec;
					muon_end_nanosec[iad] = time_nanosec + veto_time_ws_muon_microsec*1.e3 + post_veto_nanosec;
				}
				else{ //if mingling with previous muon veto - determine which one goes for longer
					proj_end_sec[iad] = time_sec;
					proj_end_nanosec[iad] = time_nanosec + veto_time_ws_muon_microsec*1.e3 + post_veto_nanosec;
					if(muon_end_sec[iad]-proj_end_sec[iad] < 1.e-9*(proj_end_nanosec[iad]-muon_end_nanosec[iad])){
						muon_end_sec[iad] = proj_end_sec[iad];
						muon_end_nanosec[iad] = proj_end_nanosec[iad];
					}
				}
			}
			
		} //End of WS muon

		if(detector < 5) energy = CorrectedEnergy(EH, detector, energy, x, y, z, time_sec, false, false);

		//Shower muon
		if(energy>2.5e3 && detector < 5){
			if(time_sec-muon_end_sec[det_num]>1.e-9*(muon_end_nanosec[det_num]-time_nanosec+pre_veto_nanosec)){ //if after the last muon time has ended
			    //Finishing last muon period
				muon_veto_start_sec = muon_start_sec[det_num];
				muon_veto_start_nanosec = muon_start_nanosec[det_num];
				muon_veto_end_sec = muon_end_sec[det_num];
				muon_veto_end_nanosec = muon_end_nanosec[det_num];
				muon_veto_sec[det_num] = muon_end_sec[det_num] - muon_start_sec[det_num];
				muon_veto_nanosec[det_num] = muon_end_nanosec[det_num] - muon_start_nanosec[det_num];
				total_veto_time[det_num] = total_veto_time[det_num] + muon_veto_sec[det_num] + 1.e-9*muon_veto_nanosec[det_num];
			    //Starting new muon period
				muon_start_sec[det_num] = time_sec;
				muon_start_nanosec[det_num] = time_nanosec - pre_veto_nanosec;
				muon_end_sec[det_num] = time_sec + veto_time_shower_muon_sec;
				muon_end_nanosec[det_num] = time_nanosec + post_veto_nanosec;
			}
			else{ //if mingling with previous muon veto - determine which one goes for longer
				proj_end_sec[det_num] = time_sec + veto_time_shower_muon_sec;
				proj_end_nanosec[det_num] = time_nanosec + post_veto_nanosec;
				if(muon_end_sec[det_num]-proj_end_sec[det_num] < 1.e-9*(proj_end_nanosec[det_num]-muon_end_nanosec[det_num])){
					muon_end_sec[det_num] = proj_end_sec[det_num];
					muon_end_nanosec[det_num] = proj_end_nanosec[det_num];
				}
			}

		}//End of Shower muon
		
		//AD muon
		if(energy>20. && energy <2.5e3 && detector <5){
			if(time_sec-muon_end_sec[det_num]>1.e-9*(muon_end_nanosec[det_num]-time_nanosec+pre_veto_nanosec)){ //if after the last muon time has ended
				//Finishing last muon period
				muon_veto_start_sec = muon_start_sec[det_num];
				muon_veto_start_nanosec = muon_start_nanosec[det_num];
				muon_veto_end_sec = muon_end_sec[det_num];
				muon_veto_end_nanosec = muon_end_nanosec[det_num];
				muon_veto_sec[det_num] = muon_end_sec[det_num] - muon_start_sec[det_num];
				muon_veto_nanosec[det_num] = muon_end_nanosec[det_num] - muon_start_nanosec[det_num];
				total_veto_time[det_num] = total_veto_time[det_num] + muon_veto_sec[det_num] + 1.e-9*muon_veto_nanosec[det_num];
			    //Starting new muon period
				muon_start_sec[det_num] = time_sec;
				muon_start_nanosec[det_num] = time_nanosec - pre_veto_nanosec;
				muon_end_sec[det_num] = time_sec;
				muon_end_nanosec[det_num] = time_nanosec + veto_time_ad_muon_microsec*1.e3 + post_veto_nanosec;
			}
			else{ //if mingling with previous muon veto - determine which one goes for longer
				proj_end_sec[det_num] = time_sec;
				proj_end_nanosec[det_num] = time_nanosec + veto_time_ad_muon_microsec*1.e3 + post_veto_nanosec;

				if(muon_end_sec[det_num]-proj_end_sec[det_num] < 1.e-9*(proj_end_nanosec[det_num]-muon_end_nanosec[det_num])){
					muon_end_sec[det_num] = proj_end_sec[det_num];
					muon_end_nanosec[det_num] = proj_end_nanosec[det_num];
				}
			}
			
		} //End of AD muon
		
		prev_trig_sec = time_sec;
		prev_trig_nanosec = time_nanosec;

		
		//DONE GOING THROUGH THE DAQ CALCULATION
		
		//TIME TO LOOK FOR DELAYED EVENTS
		if(EH == 1 && detector==1 && (run_num == 67633 || run_num == 67749 || run_num == 67755 || run_num == 67768)) continue; //EH1AD1 was run by JUNO for these runs
		
		if(((detector==5 && nhit>12) || (detector==6 && nhit>15))&&(trigType==268439552 || trigType==268435712 || trigType==268439808)){ //checking for ws muon event
			last_ws_muon_sec = time_sec;
			last_ws_muon_nanosec = time_nanosec;
			continue;
		}
		
		if(detector>4) continue; //Only want events in ADs now
		
		if(energy>2.5e3){//checking for shower muon event
			last_shower_muon_sec[detector-1] = time_sec;
			last_shower_muon_nanosec[detector-1] = time_nanosec;
			continue;
		}
		
		if(energy>20.){//checking for ad muon event
			last_ad_muon_sec[detector-1] = time_sec;
			last_ad_muon_nanosec[detector-1] = time_nanosec;
			continue;
		}
		
		if(IsFlasher(MaxQ, Quadrant, MaxQ_2inchPMT, time_PSD, time_PSD1, Q1, Q2, x, y, z)==1) continue; //Ignore flashers
		
		if(energy < singleElow || energy > singleEhigh) continue; //Only want events inside the specified energy range
		
		//Outside ws muon veto?
		if(time_sec - last_ws_muon_sec < 1.e-9*(last_ws_muon_nanosec - time_nanosec + post_veto_nanosec) + 1.e-6*veto_time_ws_muon_microsec){
			continue;
		}
		
		//Outside shower muon veto?
		if(time_sec - last_shower_muon_sec[detector-1] - veto_time_shower_muon_sec < 1.e-9*(last_shower_muon_nanosec[detector-1] - time_nanosec + post_veto_nanosec)){
			continue;
		}
		
		//Outside ad muon veto?
		if(time_sec - last_ad_muon_sec[detector-1] < 1.e-9*(last_ad_muon_nanosec[detector-1] - time_nanosec + post_veto_nanosec) + 1.e-6*veto_time_ad_muon_microsec){
			continue;
		}
		
		//Far enough away from window start time?
		if(time_sec - window_start_sec < 1.e-9*(window_start_nanosec - time_nanosec + post_veto_nanosec)){
			continue;//no? then it's not a delayed event
		}
		
		//Possible delayed event! Saving the information
		single_det = detector-1;
		s_entry[single_det] = ientry;
		single_sec[single_det] = time_sec;
		single_nanosec[single_det] = time_nanosec;
		single_energy[single_det] = energy;
		single_x[single_det] = x;
		single_y[single_det] = y;
		single_z[single_det] = z;
		
		pre_counter = 0;
		in_pre = 0;
		//Moving on to look for prompt events or other events in the window
		while(1){
			pre_counter += 1;
			//Going from the delayed event back in time
			reconT_AdSimple->GetEntry(s_entry[single_det]-pre_counter);
			calibStatsT->GetEntry(s_entry[single_det]-pre_counter);
			
			if(detector == 7) continue; //Ignore RPCs
			
			if(time_sec - single_sec[single_det] < 1.e-9*(single_nanosec[single_det] - time_nanosec) - timePreSingle*1.e-6) break; //it's outside the pre_delayed window, stop looking for events
			
			if(detector-1 != single_det) continue; //Want to look for events in the same detector

			energy = CorrectedEnergy(EH, detector, energy, x, y, z, time_sec, false, false);
			
			if(IsFlasher(MaxQ, Quadrant, MaxQ_2inchPMT, time_PSD, time_PSD1, Q1, Q2, x, y, z)==1) continue; //Ignore flashers
			
			if(energy < sElow || energy > sEhigh) continue;
			
			in_pre += 1; //Prompt-like event in 800us window before the delayed-like event
			
		} //end of looking for prompt events

		if(in_pre != 0) continue; //Don't want any events around the single
		
		post_counter = 0;
		incomplete_post_window = 0;
		temp_prev_trig_sec = single_sec[single_det];
		temp_prev_trig_nanosec = single_nanosec[single_det];
		//Now to look for post-single events
		while(1){
			post_counter += 1;
			
			if(s_entry[single_det]+post_counter > nentries-1){ //if no more data for run
				incomplete_post_window = 1;
				break;
			}
			
			//Going from delayed to after
			reconT_AdSimple->GetEntry(s_entry[single_det]+post_counter);
			calibStatsT->GetEntry(s_entry[single_det]+post_counter);
			
			if(detector == 7) continue; //Ignore RPCs
			
			//Check for time gaps:
			if(time_sec - temp_prev_trig_sec - 1. > 1.e-9*(temp_prev_trig_nanosec - time_nanosec)){
				incomplete_post_window = 1;
				break;
			}
			
			temp_prev_trig_sec = time_sec;
			temp_prev_trig_nanosec = time_nanosec;
			
			//If it is outside the post_delayed window, then it's good & we have a single event
			if(time_sec - single_sec[single_det] > 1.e-9*(single_nanosec[single_det] - time_nanosec) + 1.e-6*timePostSingle){
				break; //if this is the case, we found one!
			}
			
			//Check for ws muon
			if(((detector==5 && nhit>12) || (detector==6 && nhit>15))&&(trigType==268439552 || trigType==268435712 || trigType==268439808)){
				incomplete_post_window = 1;
				break;
			}

			if(detector<5) energy = CorrectedEnergy(EH, detector, energy, x, y, z, time_sec, false, false);
			
			//Check for shower or ad muon
			if(energy > 20.){
				incomplete_post_window = 1;
				break;
			}
			
			if(detector-1 != single_det) continue; //Want only events in the same detector
			
			if(IsFlasher(MaxQ, Quadrant, MaxQ_2inchPMT, time_PSD, time_PSD1, Q1, Q2, x, y, z)==1) continue; //Ignore flashers
			
			//if there's a prompt-like event in the post-delayed window, not counted
			if(energy > sElow && energy < sEhigh){
				incomplete_post_window = 1;
				break;
			}
			
		} //End of looking for post-delayed events
		
		//If the post-delayed window was broken, go on looking for other IBD-like events
		if(incomplete_post_window == 1) continue;
		
		
		//----------------SAVING TO TTREE----------------
		hall_num = EH;
		det_num = single_det;
		s_energy = single_energy[single_det];
		s_sec = single_sec[single_det];
		s_nanosec = single_nanosec[single_det];
		s_x = single_x[single_det];
		s_y = single_y[single_det];
		s_z = single_z[single_det];

		
		tr_sing->Fill();
		
		
	}//Done going through the data
	
	//Finish the time calculations and print information
	cout << "Time between last and first triggers: " << last_trig_sec - first_trig_sec + 1.e-9*(last_trig_nanosec - first_trig_nanosec) << endl;
	cout << "Total DAQ time: " << total_DAQ_time << "s" << endl;
	for(int iad=0; iad<4; iad++){
		if(check_ad[iad] == 0){
			if(iad == 0){
				tot_veto_ad1 = 0;
				tot_live_ad1 = 0;
				total_DAQ_ad1 = 0;
				continue;
			}
			if(iad == 1){
				tot_veto_ad2 = 0;
				tot_live_ad2 = 0;
				total_DAQ_ad2 = 0;
				continue;
			}
			if(iad == 2){
				tot_veto_ad3 = 0;
				tot_live_ad3 = 0;
				total_DAQ_ad3 = 0;
				continue;
			}
			if(iad == 3){
				tot_veto_ad4 = 0;
				tot_live_ad4 = 0;
				total_DAQ_ad4 = 0;
				continue;
			}
		}
		else{
			if(iad == 0){
				tot_veto_ad1 = total_veto_time[iad];
				tot_live_ad1 = total_DAQ_time - total_veto_time[iad];
				total_DAQ_ad1 = total_DAQ_time;
			}
			if(iad == 1){
				tot_veto_ad2 = total_veto_time[iad];
				tot_live_ad2 = total_DAQ_time - total_veto_time[iad];
				total_DAQ_ad2 = total_DAQ_time;
			}
			if(iad == 2){
				tot_veto_ad3 = total_veto_time[iad];
				tot_live_ad3 = total_DAQ_time - total_veto_time[iad];
				total_DAQ_ad3 = total_DAQ_time;
			}
			if(iad == 3){
				tot_veto_ad4 = total_veto_time[iad];
				tot_live_ad4 = total_DAQ_time - total_veto_time[iad];
				total_DAQ_ad4 = total_DAQ_time;
			}
		}
		if(run_num == 67633 || run_num == 67749 || run_num == 67755 || run_num == 67768){
			tot_veto_ad1 = 0;
			tot_live_ad1 = 0;
			total_DAQ_ad1 = 0;
		}
		cout << "AD#" << iad << "\t Total vetoed time: " << total_veto_time[iad] << "s\t Percentage live: " << (total_DAQ_time - total_veto_time[iad])*100/total_DAQ_time << "%" << endl;
	}
	
	tr_summary->Fill();
	
	
	//-----------SAVE TO OUTPUT FILE------------
        char* outputname;
	if(timePreSingle == pd_window_microsec+400){
//		sprintf(outputname,"./delayedSingles/EH%d/foundDelayedSingles_%d_%d.root",EH,pd_window_microsec,run_num);
		outputname = outputname_findDelayed;
	}
	else{
		outputname = outputname_findPrompt;
	}
/*	if(timePreSingle == pd_window_microsec+300){
		sprintf(outputname,"./promptSingles/EH%d/foundPromptSingles_%d_%d.root",EH,pd_window_microsec,run_num);
	}
	if(timePostSingle == pd_window_microsec+300){
		sprintf(outputname,"./promptSingles/EH%d/foundUNPromptSingles_%d_%d.root",EH,pd_window_microsec,run_num);
	}*/
/*	if(timePreSingle == 2400){
		sprintf(outputname,"./delayedSingles/EH%d/foundDelayedSingles_TcLong_%d.root",EH,run_num);
	}
	if(timePreSingle == 2300 || timePreSingle == 500){
		sprintf(outputname,"./promptSingles/EH%d/foundPromptSingles_TcLong_%d.root",EH,run_num);
	}*/

	cout << "Writing to: " << outputname << endl;

	TFile* outfile=new TFile(outputname, "RECREATE");
		outfile->cd();
		tr_sing->Write();
		tr_DAQ->Write();
		tr_summary->Write();
	outfile->Close();

}





void summarize(int run_num, float pElow, float pEhigh, float dElow, float dEhigh, int pd_window_microsec){ //Taking the data from above and pairing the events

	int EH = 0;
	int runNum = 0;
	int maxAD = 0;
    FILE* runfile=fopen("./run_list_good_sync.txt","r");  //This file contains all good runs numbers and experimental hall they belong to
    while(1){ //go though the file and retrieve EH which corresponds to run_num provided as a parameter
        fscanf(runfile,"%d %d",&runNum,&EH);
        if(feof(runfile)) break;
        if(run_num == runNum) break;
    }
    fclose(runfile);
    
//	const char* OutputFileName_summarize = outputname_summarize;

    cout<<"Run identified as: run#"<<run_num<<" EH"<<EH<<endl;
	cout << "n-4 " << outputname_summarize << endl;
	if(EH == 1 || EH == 2) maxAD = 2;
	if(EH == 3) maxAD = 4;

	//Making Histograms and TTree:
		TH1F* h_p_singles_energy_raw[maxAD]; //energy histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_raw_ad%d",iad+1);
			h_p_singles_energy_raw[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_d_singles_energy_raw[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_raw_ad%d",iad+1);
		//	h_d_singles_energy_raw[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_raw[iad]=new TH1F(name,name,230,0.7,3.);
		}


	double p_DAQ_ad1 = 0;
	double p_veto_ad1 = 0;
	double p_live_ad1 = 0;
	double p_DAQ_ad2 = 0;
	double p_veto_ad2 = 0;
	double p_live_ad2 = 0;
	double p_DAQ_ad3 = 0;
	double p_veto_ad3 = 0;
	double p_live_ad3 = 0;
	double p_DAQ_ad4 = 0;
	double p_veto_ad4 = 0;
	double p_live_ad4 = 0;
	double p_counts = 0;
	double d_DAQ_ad1 = 0;
	double d_veto_ad1 = 0;
	double d_live_ad1 = 0;
	double d_DAQ_ad2 = 0;
	double d_veto_ad2 = 0;
	double d_live_ad2 = 0;
	double d_DAQ_ad3 = 0;
	double d_veto_ad3 = 0;
	double d_live_ad3 = 0;
	double d_DAQ_ad4 = 0;
	double d_veto_ad4 = 0;
	double d_live_ad4 = 0;
	double d_counts = 0;


		TTree* tr_live=new TTree("tr_singles_summary", "Summary of Singles Search");
			tr_live->Branch("p_DAQ_ad1", &p_DAQ_ad1, "DAQ for AD1/Double");
			tr_live->Branch("p_veto_ad1", &p_veto_ad1, "Veto for AD1/Double");
			tr_live->Branch("p_live_ad1", &p_live_ad1, "Live for AD1/Double");
			tr_live->Branch("p_DAQ_ad2", &p_DAQ_ad2, "DAQ for AD2/Double");
			tr_live->Branch("p_veto_ad2", &p_veto_ad2, "Veto for AD2/Double");
			tr_live->Branch("p_live_ad2", &p_live_ad2, "Live for AD2/Double");
			tr_live->Branch("p_DAQ_ad3", &p_DAQ_ad3, "DAQ for AD3/Double");
			tr_live->Branch("p_veto_ad3", &p_veto_ad3, "Veto for AD3/Double");
			tr_live->Branch("p_live_ad3", &p_live_ad3, "Live for AD3/Double");
			tr_live->Branch("p_DAQ_ad4", &p_DAQ_ad4, "DAQ for AD4/Double");
			tr_live->Branch("p_veto_ad4", &p_veto_ad4, "Veto for AD4/Double");
			tr_live->Branch("p_live_ad4", &p_live_ad4, "Live for AD4/Double");
			tr_live->Branch("d_DAQ_ad1", &d_DAQ_ad1, "DAQ for AD1/Double");
			tr_live->Branch("d_veto_ad1", &d_veto_ad1, "Veto for AD1/Double");
			tr_live->Branch("d_live_ad1", &d_live_ad1, "Live for AD1/Double");
			tr_live->Branch("d_DAQ_ad2", &d_DAQ_ad2, "DAQ for AD2/Double");
			tr_live->Branch("d_veto_ad2", &d_veto_ad2, "Veto for AD2/Double");
			tr_live->Branch("d_live_ad2", &d_live_ad2, "Live for AD2/Double");
			tr_live->Branch("d_DAQ_ad3", &d_DAQ_ad3, "DAQ for AD3/Double");
			tr_live->Branch("d_veto_ad3", &d_veto_ad3, "Veto for AD3/Double");
			tr_live->Branch("d_live_ad3", &d_live_ad3, "Live for AD3/Double");
			tr_live->Branch("d_DAQ_ad4", &d_DAQ_ad4, "DAQ for AD4/Double");
			tr_live->Branch("d_veto_ad4", &d_veto_ad4, "Veto for AD4/Double");
			tr_live->Branch("d_live_ad4", &d_live_ad4, "Live for AD4/Double");


	float p_energy = 0;
	float d_energy = 0;
	float p_x = 0;
	float p_y = 0;
	float p_z = 0;
	float d_x = 0;
	float d_y = 0;
	float d_z = 0;
	float distance = 0;
	int det_num = 0;

	int d_det_num = 0;
	int p_det_num = 0;
	int p_entries = 0;
	int d_entries = 0;
	int p_rand_entry = 0;
	int d_rand_entry = 0;
	int d_hall_num = 0;
	int p_hall_num = 0;

//Here is where we load in the singles:

	//Reading in the prompt file
//        char inputname[64];
//        sprintf(inputname,"./promptSingles/EH%d/foundUNPromptSingles_%d_%d.root",EH,pd_window_microsec,run_num);
//        sprintf(inputname,"./promptSingles/EH%d/foundPromptSingles_TcLong_%d.root",EH,run_num);
//        sprintf(inputname,"./foundPromptSingles_%d.root",run_num);
	TFile *promptInfile = new TFile(outputname_findPrompt,"READ");


		//TTree tr_DAQ:
		TTree *tr_DAQ = (TTree*)promptInfile->Get("tr_DAQ");
		tr_DAQ->SetBranchStatus("*",1);

			int window_start_sec = 0;
			int window_start_nanosec = 0;
			int window_end_sec = 0;
			int window_end_nanosec = 0;
			int Dentries = 0;

			tr_DAQ->SetBranchAddress("window_start_sec",&window_start_sec);
			tr_DAQ->SetBranchAddress("window_start_nanosec",&window_start_nanosec);
			tr_DAQ->SetBranchAddress("window_end_sec",&window_end_sec);
			tr_DAQ->SetBranchAddress("window_end_nanosec",&window_end_nanosec);


	//Reading in the delayed file
 //       sprintf(inputname,"./delayedSingles/EH%d/foundDelayedSingles_%d_%d.root",EH,pd_window_microsec,run_num);
//        sprintf(inputname,"./delayedSingles/EH%d/foundDelayedSingles_TcLong_%d.root",EH,run_num);
//        sprintf(inputname,"./foundDelayedSingles_%d.root",run_num);
	TFile *delayedInfile = new TFile(outputname_findDelayed,"READ");
	cout << "n-3 " << outputname_summarize << endl;
		//Establishing the new TTree of randomly paired singles
		TTree* tr_pair=new TTree("tr_pair", "Random pair of Prompt and Delayed Singles");
		tr_pair->Branch("p_energy", &p_energy, "Prompt Energy MeV/Float");
		tr_pair->Branch("d_energy", &d_energy, "Delayed Energy MeV/Float");
		tr_pair->Branch("p_x", &p_x, "X position of prompt in cm/Float");
		tr_pair->Branch("d_x", &d_x, "X position of delayed in cm/Float");
		tr_pair->Branch("p_y", &p_y, "Y position of prompt in cm/Float");
		tr_pair->Branch("d_y", &d_y, "Y position of delayed in cm/Float");
		tr_pair->Branch("p_z", &p_z, "Z position of prompt in cm/Float");
		tr_pair->Branch("d_z", &d_z, "z position of delayed in cm/Float");
		tr_pair->Branch("det_num", &det_num, "Detector number/Int");
		tr_pair->Branch("distance", &distance, "Distance between Prompt and Delayed in cm/Float");


		//TTree tr_sum
		TTree *tr_p_sum = (TTree*)promptInfile->Get("tr_summary");
		tr_p_sum->SetBranchStatus("*",1);

		TTree *tr_d_sum = (TTree*)delayedInfile->Get("tr_summary");
		tr_d_sum->SetBranchStatus("*",1);

	double total_DAQ_ad1 = 0;
	double tot_veto_ad1 = 0;
	double tot_live_ad1 = 0;
	double total_DAQ_ad2 = 0;
	double tot_veto_ad2 = 0;
	double tot_live_ad2 = 0;
	double total_DAQ_ad3 = 0;
	double tot_veto_ad3 = 0;
	double tot_live_ad3 = 0;
	double total_DAQ_ad4 = 0;
	double tot_veto_ad4 = 0;
	double tot_live_ad4 = 0;

	tr_p_sum->SetBranchAddress("total_DAQ_ad1",&p_DAQ_ad1);
	tr_p_sum->SetBranchAddress("tot_veto_ad1",&p_veto_ad1);
	tr_p_sum->SetBranchAddress("tot_live_ad1",&p_live_ad1);
	tr_p_sum->SetBranchAddress("total_DAQ_ad2",&p_DAQ_ad2);
	tr_p_sum->SetBranchAddress("tot_veto_ad2",&p_veto_ad2);
	tr_p_sum->SetBranchAddress("tot_live_ad2",&p_live_ad2);
	tr_p_sum->SetBranchAddress("total_DAQ_ad3",&p_DAQ_ad3);
	tr_p_sum->SetBranchAddress("tot_veto_ad3",&p_veto_ad3);
	tr_p_sum->SetBranchAddress("tot_live_ad3",&p_live_ad3);
	tr_p_sum->SetBranchAddress("total_DAQ_ad4",&p_DAQ_ad4);
	tr_p_sum->SetBranchAddress("tot_veto_ad4",&p_veto_ad4);
	tr_p_sum->SetBranchAddress("tot_live_ad4",&p_live_ad4);

	tr_d_sum->SetBranchAddress("total_DAQ_ad1",&d_DAQ_ad1);
	tr_d_sum->SetBranchAddress("tot_veto_ad1",&d_veto_ad1);
	tr_d_sum->SetBranchAddress("tot_live_ad1",&d_live_ad1);
	tr_d_sum->SetBranchAddress("total_DAQ_ad2",&d_DAQ_ad2);
	tr_d_sum->SetBranchAddress("tot_veto_ad2",&d_veto_ad2);
	tr_d_sum->SetBranchAddress("tot_live_ad2",&d_live_ad2);
	tr_d_sum->SetBranchAddress("total_DAQ_ad3",&d_DAQ_ad3);
	tr_d_sum->SetBranchAddress("tot_veto_ad3",&d_veto_ad3);
	tr_d_sum->SetBranchAddress("tot_live_ad3",&d_live_ad3);
	tr_d_sum->SetBranchAddress("total_DAQ_ad4",&d_DAQ_ad4);
	tr_d_sum->SetBranchAddress("tot_veto_ad4",&d_veto_ad4);
	tr_d_sum->SetBranchAddress("tot_live_ad4",&d_live_ad4);
		
		//TTree p_tree
		TTree *p_tree = (TTree*)promptInfile->Get("tr_sing");
		p_tree->SetBranchStatus("*",1);

		//TTree d_tree
		TTree *d_tree = (TTree*)delayedInfile->Get("tr_sing");
		d_tree->SetBranchStatus("*",1);
		d_tree->SetBranchStatus("s_sec",0);
		d_tree->SetBranchStatus("s_nanosec",0);
//		d_tree->SetBranchStatus("run",0);
//		d_tree->SetBranchStatus("s_radius",0);



	float s_energy=0;
	int s_sec=0;
	int s_nanosec=0;
	float s_x=0;
	float s_y=0;
	float s_z=0;
	float s_radius=0;
	int hall_num=0;
	int run=0;


	p_tree->SetBranchAddress("s_energy",&s_energy);
	p_tree->SetBranchAddress("s_sec",&s_sec);
	p_tree->SetBranchAddress("s_nanosec",&s_nanosec);
	p_tree->SetBranchAddress("s_x",&s_x);
	p_tree->SetBranchAddress("s_y",&s_y);
	p_tree->SetBranchAddress("s_z",&s_z);
//	p_tree->SetBranchAddress("s_radius", &s_radius);
	p_tree->SetBranchAddress("det_num",&det_num);
	p_tree->SetBranchAddress("hall_num",&hall_num);
//	p_tree->SetBranchAddress("run",&run);


	d_tree->SetBranchAddress("s_energy",&d_energy);
	d_tree->SetBranchAddress("s_x",&d_x);
	d_tree->SetBranchAddress("s_y",&d_y);
	d_tree->SetBranchAddress("s_z",&d_z);
	d_tree->SetBranchAddress("hall_num",&d_hall_num);
	d_tree->SetBranchAddress("det_num",&d_det_num);

	int nentries = tr_p_sum->GetEntries();

	cout << "nentries = " << nentries << endl;
	cout << "n-2 " << outputname_summarize << endl;
		tr_p_sum->GetEntry(0);
		tr_d_sum->GetEntry(0);
		tr_live->Fill();

	cout << fixed;
	cout << "Delayed live: " << setprecision(9) << d_live_ad1 << endl;
	cout << "Prompt live: " << setprecision(9) << p_live_ad1 << endl;

	p_entries = p_tree->GetEntries();
	d_entries = d_tree->GetEntries();

	cout << "Delayed entries: " << d_entries << endl;
	cout << "Prompt entries: " << p_entries << endl;

//	d_entries = 15000;
//	p_entries = 10000;

	//DEALING WITH DATA
	for(int ientry=0; ientry<d_entries; ientry++){
		d_tree->GetEntry(ientry);

		if(d_det_num > (maxAD-1)) continue;

		h_d_singles_energy_raw[d_det_num]->Fill(d_energy);
	}

	for(int ientry=0; ientry<p_entries; ientry++){
		p_tree->GetEntry(ientry);

		if(det_num > (maxAD-1)) continue;

		h_p_singles_energy_raw[det_num]->Fill(s_energy);

	} //End of filling histograms

	cout << "n-1 " << outputname_summarize << endl;

	//Random Funtion:
	TRandom3* random = new TRandom3();

	//Making tr_pair:
//For pairing prompt to each delayed:
	if(d_entries <= p_entries){
		cout << "Pairing a prompt to each delayed event" << endl;
		for(int ientry=0; ientry<d_entries; ientry++){
	//	for(int ientry=0; ientry<(d_entries/2); ientry++){
			d_x = 0;
			d_y = 0;
			d_z = 0;
			if(ientry%10000 ==0) cout<< ientry << "/" << d_entries << endl;
			d_tree->GetEntry(ientry);
			if(d_energy < dElow || d_energy > dEhigh) continue;

			p_rand_entry =(int) random->Uniform(0,p_entries);
	//		p_rand_entry = ientry;
			p_x = 0;
			p_y = 0;
			p_z = 0;
			p_tree->GetEntry(p_rand_entry);

			//Making sure it's in the same detector
			while(d_hall_num != hall_num || d_det_num != det_num || s_energy < pElow || s_energy > pEhigh){
				p_rand_entry =(int) random->Uniform(0,p_entries);
		//		p_rand_entry+=1;
				p_tree->GetEntry(p_rand_entry);
			}
			
			//Renaming prompt variables to save to new TTree
			p_x = s_x;
			p_y = s_y;
			p_z = s_z;
			p_energy = s_energy;
			distance = sqrt(pow(p_x-d_x,2)+pow(p_y-d_y,2)+pow(p_z-d_z,2));
			tr_pair->Fill();
		}
	}

//For pairing each delayed to a prompt:
	else{
		cout << "Pairing a delayed to each prompt event" << endl;
		for(int ientry=0; ientry<p_entries; ientry++){
			p_x = 0;
			p_y = 0;
			p_z = 0;
			if(ientry%10000 ==0) cout<< ientry << "/" << p_entries << endl;
			p_tree->GetEntry(ientry);
			p_x = s_x;
			p_y = s_y;
			p_z = s_z;
			p_energy = s_energy;
			p_hall_num = hall_num;
			p_det_num = det_num;
			if(p_energy < pElow || p_energy > pEhigh) continue;

			d_rand_entry =(int) random->Uniform(0,d_entries);
	//		d_rand_entry = ientry;
			d_x = 0;
			d_y = 0;
			d_z = 0;
			d_tree->GetEntry(d_rand_entry);

			//Making sure it's in the same detector
			while(p_hall_num != d_hall_num || p_det_num != d_det_num || d_energy < dElow || d_energy > dEhigh){
				d_rand_entry =(int) random->Uniform(0,d_entries);
		//		d_rand_entry+=1;
				d_tree->GetEntry(d_rand_entry);
			}

			distance = sqrt(pow(p_x-d_x,2)+pow(p_y-d_y,2)+pow(p_z-d_z,2));
			tr_pair->Fill();
		}
	}
	cout << "n " << outputname_summarize << endl;
	cout << "Number of pairs: " << tr_pair->GetEntries() << endl << endl << endl << endl;

	//SAVING TO OUTPUT FILE
//        char outputname[64];
//        sprintf(outputname,"./accResults/EH%d/AccidentalsSummary_%d_%d.root",EH,pd_window_microsec,run_num);
//	if(dElow > 1.5) sprintf(outputname,"./accResults/EH%d/3sig/AccidentalsSummary_%d_%d.root",EH,pd_window_microsec,run_num);
//        sprintf(outputname,"./accResults/EH%d/AccidentalsSummary_TcLong_%d.root",EH,run_num);
//        sprintf(outputname,"./accResults/EH%d/practiceSummary_%d.root",EH,run_num);
//        sprintf(outputname,"./AccidentalsSummary_%d.root",run_num);

	cout << "Saving to: " << outputname_summarize << endl;
	TFile* outfile=new TFile(outputname_summarize, "RECREATE");
		outfile->cd();
		tr_live->Write();
		tr_pair->Write();

		for(int iad = 0; iad < maxAD; ++iad){
		//	h_p_singles_energy_raw[iad]->SetStats(0);
			h_p_singles_energy_raw[iad]->GetXaxis()->SetTitle("Prompt Singles Energy [MeV]");
			h_p_singles_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_raw[iad]->Write();

		//	h_d_singles_energy_raw[iad]->SetStats(0);
			h_d_singles_energy_raw[iad]->GetXaxis()->SetTitle("Delayed Singles Energy [MeV]");
			h_d_singles_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_raw[iad]->Write();

		}

	outfile->Close();

}

void summarize_long(int run_num, int part_num, float pElow, float pEhigh, float dElow, float dEhigh, int pd_window_microsec){ //Taking the data from above and making plots of it

	int EH = 0;
	int runNum = 0;
	int maxAD = 0;
    FILE* runfile=fopen("./run_list_good_sync.txt","r");  //This file contains all good runs numbers and experimental hall they belong to
    while(1){ //go though the file and retrieve EH which corresponds to run_num provided as a parameter
        fscanf(runfile,"%d %d",&runNum,&EH);
        if(feof(runfile)) break;
        if(run_num == runNum) break;
    }
    fclose(runfile);
    
    cout<<"Run identified as: run#"<<run_num<<" EH"<<EH<<endl;

	if(EH == 1 || EH == 2) maxAD = 2;
	if(EH == 3) maxAD = 4;

	//Making Histograms and TTree:
		TH1F* h_p_singles_energy_raw[maxAD]; //energy histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_raw_ad%d",iad+1);
			h_p_singles_energy_raw[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_d_singles_energy_raw[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_raw_ad%d",iad+1);
		//	h_d_singles_energy_raw[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_raw[iad]=new TH1F(name,name,230,0.7,3.);
		}


	double p_DAQ_ad1 = 0;
	double p_veto_ad1 = 0;
	double p_live_ad1 = 0;
	double p_DAQ_ad2 = 0;
	double p_veto_ad2 = 0;
	double p_live_ad2 = 0;
	double p_DAQ_ad3 = 0;
	double p_veto_ad3 = 0;
	double p_live_ad3 = 0;
	double p_DAQ_ad4 = 0;
	double p_veto_ad4 = 0;
	double p_live_ad4 = 0;
	double p_counts = 0;
	double d_DAQ_ad1 = 0;
	double d_veto_ad1 = 0;
	double d_live_ad1 = 0;
	double d_DAQ_ad2 = 0;
	double d_veto_ad2 = 0;
	double d_live_ad2 = 0;
	double d_DAQ_ad3 = 0;
	double d_veto_ad3 = 0;
	double d_live_ad3 = 0;
	double d_DAQ_ad4 = 0;
	double d_veto_ad4 = 0;
	double d_live_ad4 = 0;
	double d_counts = 0;


		TTree* tr_live=new TTree("tr_singles_summary", "Summary of Singles Search");
			tr_live->Branch("p_DAQ_ad1", &p_DAQ_ad1, "DAQ for AD1/Double");
			tr_live->Branch("p_veto_ad1", &p_veto_ad1, "Veto for AD1/Double");
			tr_live->Branch("p_live_ad1", &p_live_ad1, "Live for AD1/Double");
			tr_live->Branch("p_DAQ_ad2", &p_DAQ_ad2, "DAQ for AD2/Double");
			tr_live->Branch("p_veto_ad2", &p_veto_ad2, "Veto for AD2/Double");
			tr_live->Branch("p_live_ad2", &p_live_ad2, "Live for AD2/Double");
			tr_live->Branch("p_DAQ_ad3", &p_DAQ_ad3, "DAQ for AD3/Double");
			tr_live->Branch("p_veto_ad3", &p_veto_ad3, "Veto for AD3/Double");
			tr_live->Branch("p_live_ad3", &p_live_ad3, "Live for AD3/Double");
			tr_live->Branch("p_DAQ_ad4", &p_DAQ_ad4, "DAQ for AD4/Double");
			tr_live->Branch("p_veto_ad4", &p_veto_ad4, "Veto for AD4/Double");
			tr_live->Branch("p_live_ad4", &p_live_ad4, "Live for AD4/Double");
			tr_live->Branch("d_DAQ_ad1", &d_DAQ_ad1, "DAQ for AD1/Double");
			tr_live->Branch("d_veto_ad1", &d_veto_ad1, "Veto for AD1/Double");
			tr_live->Branch("d_live_ad1", &d_live_ad1, "Live for AD1/Double");
			tr_live->Branch("d_DAQ_ad2", &d_DAQ_ad2, "DAQ for AD2/Double");
			tr_live->Branch("d_veto_ad2", &d_veto_ad2, "Veto for AD2/Double");
			tr_live->Branch("d_live_ad2", &d_live_ad2, "Live for AD2/Double");
			tr_live->Branch("d_DAQ_ad3", &d_DAQ_ad3, "DAQ for AD3/Double");
			tr_live->Branch("d_veto_ad3", &d_veto_ad3, "Veto for AD3/Double");
			tr_live->Branch("d_live_ad3", &d_live_ad3, "Live for AD3/Double");
			tr_live->Branch("d_DAQ_ad4", &d_DAQ_ad4, "DAQ for AD4/Double");
			tr_live->Branch("d_veto_ad4", &d_veto_ad4, "Veto for AD4/Double");
			tr_live->Branch("d_live_ad4", &d_live_ad4, "Live for AD4/Double");


	float p_energy = 0;
	float d_energy = 0;
	float p_x = 0;
	float p_y = 0;
	float p_z = 0;
	float d_x = 0;
	float d_y = 0;
	float d_z = 0;
	float distance = 0;
	int det_num = 0;

	int d_det_num = 0;
	int p_det_num = 0;
	int p_entries = 0;
	int d_entries = 0;
	int p_rand_entry = 0;
	int d_rand_entry = 0;
	int d_hall_num = 0;
	int p_hall_num = 0;

	//Reading in the prompt file
//        char inputname[64];
//        sprintf(inputname,"./promptSingles/EH%d/foundUNPromptSingles_%d_%d.root",EH,pd_window_microsec,run_num);
//        sprintf(inputname,"./promptSingles/EH%d/foundPromptSingles_TcLong_%d.root",EH,run_num);
//        sprintf(inputname,"./foundPromptSingles_%d.root",run_num);
	TFile *promptInfile = new TFile(outputname_findPrompt,"READ");


		//TTree tr_DAQ:
		TTree *tr_DAQ = (TTree*)promptInfile->Get("tr_DAQ");
		tr_DAQ->SetBranchStatus("*",1);

			int window_start_sec = 0;
			int window_start_nanosec = 0;
			int window_end_sec = 0;
			int window_end_nanosec = 0;
			int Dentries = 0;

			tr_DAQ->SetBranchAddress("window_start_sec",&window_start_sec);
			tr_DAQ->SetBranchAddress("window_start_nanosec",&window_start_nanosec);
			tr_DAQ->SetBranchAddress("window_end_sec",&window_end_sec);
			tr_DAQ->SetBranchAddress("window_end_nanosec",&window_end_nanosec);


	//Reading in the delayed file
//        sprintf(inputname,"./delayedSingles/EH%d/foundDelayedSingles_%d_%d.root",EH,pd_window_microsec,run_num);
//        sprintf(inputname,"./delayedSingles/EH%d/foundDelayedSingles_TcLong_%d.root",EH,run_num);
//        sprintf(inputname,"./foundDelayedSingles_%d.root",run_num);
	TFile *delayedInfile = new TFile(outputname_findDelayed,"READ");

		//Establishing the new TTree of randomly paired singles
		TTree* tr_pair=new TTree("tr_pair", "Random pair of Prompt and Delayed Singles");
		tr_pair->Branch("p_energy", &p_energy, "Prompt Energy MeV/Float");
		tr_pair->Branch("d_energy", &d_energy, "Delayed Energy MeV/Float");
		tr_pair->Branch("p_x", &p_x, "X position of prompt in cm/Float");
		tr_pair->Branch("d_x", &d_x, "X position of delayed in cm/Float");
		tr_pair->Branch("p_y", &p_y, "Y position of prompt in cm/Float");
		tr_pair->Branch("d_y", &d_y, "Y position of delayed in cm/Float");
		tr_pair->Branch("p_z", &p_z, "Z position of prompt in cm/Float");
		tr_pair->Branch("d_z", &d_z, "z position of delayed in cm/Float");
		tr_pair->Branch("det_num", &det_num, "Detector number/Int");
		tr_pair->Branch("distance", &distance, "Distance between Prompt and Delayed in cm/Float");


		//TTree tr_sum
		TTree *tr_p_sum = (TTree*)promptInfile->Get("tr_summary");
		tr_p_sum->SetBranchStatus("*",1);

		TTree *tr_d_sum = (TTree*)delayedInfile->Get("tr_summary");
		tr_d_sum->SetBranchStatus("*",1);

	double total_DAQ_ad1 = 0;
	double tot_veto_ad1 = 0;
	double tot_live_ad1 = 0;
	double total_DAQ_ad2 = 0;
	double tot_veto_ad2 = 0;
	double tot_live_ad2 = 0;
	double total_DAQ_ad3 = 0;
	double tot_veto_ad3 = 0;
	double tot_live_ad3 = 0;
	double total_DAQ_ad4 = 0;
	double tot_veto_ad4 = 0;
	double tot_live_ad4 = 0;

	tr_p_sum->SetBranchAddress("total_DAQ_ad1",&p_DAQ_ad1);
	tr_p_sum->SetBranchAddress("tot_veto_ad1",&p_veto_ad1);
	tr_p_sum->SetBranchAddress("tot_live_ad1",&p_live_ad1);
	tr_p_sum->SetBranchAddress("total_DAQ_ad2",&p_DAQ_ad2);
	tr_p_sum->SetBranchAddress("tot_veto_ad2",&p_veto_ad2);
	tr_p_sum->SetBranchAddress("tot_live_ad2",&p_live_ad2);
	tr_p_sum->SetBranchAddress("total_DAQ_ad3",&p_DAQ_ad3);
	tr_p_sum->SetBranchAddress("tot_veto_ad3",&p_veto_ad3);
	tr_p_sum->SetBranchAddress("tot_live_ad3",&p_live_ad3);
	tr_p_sum->SetBranchAddress("total_DAQ_ad4",&p_DAQ_ad4);
	tr_p_sum->SetBranchAddress("tot_veto_ad4",&p_veto_ad4);
	tr_p_sum->SetBranchAddress("tot_live_ad4",&p_live_ad4);

	tr_d_sum->SetBranchAddress("total_DAQ_ad1",&d_DAQ_ad1);
	tr_d_sum->SetBranchAddress("tot_veto_ad1",&d_veto_ad1);
	tr_d_sum->SetBranchAddress("tot_live_ad1",&d_live_ad1);
	tr_d_sum->SetBranchAddress("total_DAQ_ad2",&d_DAQ_ad2);
	tr_d_sum->SetBranchAddress("tot_veto_ad2",&d_veto_ad2);
	tr_d_sum->SetBranchAddress("tot_live_ad2",&d_live_ad2);
	tr_d_sum->SetBranchAddress("total_DAQ_ad3",&d_DAQ_ad3);
	tr_d_sum->SetBranchAddress("tot_veto_ad3",&d_veto_ad3);
	tr_d_sum->SetBranchAddress("tot_live_ad3",&d_live_ad3);
	tr_d_sum->SetBranchAddress("total_DAQ_ad4",&d_DAQ_ad4);
	tr_d_sum->SetBranchAddress("tot_veto_ad4",&d_veto_ad4);
	tr_d_sum->SetBranchAddress("tot_live_ad4",&d_live_ad4);
		
		//TTree p_tree
		TTree *p_tree = (TTree*)promptInfile->Get("tr_sing");
		p_tree->SetBranchStatus("*",1);

		//TTree d_tree
		TTree *d_tree = (TTree*)delayedInfile->Get("tr_sing");
		d_tree->SetBranchStatus("*",1);
		d_tree->SetBranchStatus("s_sec",0);
		d_tree->SetBranchStatus("s_nanosec",0);
//		d_tree->SetBranchStatus("run",0);
//		d_tree->SetBranchStatus("s_radius",0);



	float s_energy=0;
	int s_sec=0;
	int s_nanosec=0;
	float s_x=0;
	float s_y=0;
	float s_z=0;
	float s_radius=0;
	int hall_num=0;
	int run=0;


	p_tree->SetBranchAddress("s_energy",&s_energy);
	p_tree->SetBranchAddress("s_sec",&s_sec);
	p_tree->SetBranchAddress("s_nanosec",&s_nanosec);
	p_tree->SetBranchAddress("s_x",&s_x);
	p_tree->SetBranchAddress("s_y",&s_y);
	p_tree->SetBranchAddress("s_z",&s_z);
//	p_tree->SetBranchAddress("s_radius", &s_radius);
	p_tree->SetBranchAddress("det_num",&det_num);
	p_tree->SetBranchAddress("hall_num",&hall_num);
//	p_tree->SetBranchAddress("run",&run);


	d_tree->SetBranchAddress("s_energy",&d_energy);
	d_tree->SetBranchAddress("s_x",&d_x);
	d_tree->SetBranchAddress("s_y",&d_y);
	d_tree->SetBranchAddress("s_z",&d_z);
	d_tree->SetBranchAddress("hall_num",&d_hall_num);
	d_tree->SetBranchAddress("det_num",&d_det_num);

	int nentries = tr_p_sum->GetEntries();

	cout << "nentries = " << nentries << endl;

		tr_p_sum->GetEntry(0);
		tr_d_sum->GetEntry(0);
		tr_live->Fill();

	cout << fixed;
	cout << "Delayed live: " << setprecision(9) << d_live_ad1 << endl;
	cout << "Prompt live: " << setprecision(9) << p_live_ad1 << endl;

	p_entries = p_tree->GetEntries();
	d_entries = d_tree->GetEntries();

	cout << "Delayed entries: " << d_entries << endl;
	cout << "Prompt entries: " << p_entries << endl;

//	d_entries = 15000;
//	p_entries = 10000;

	//DEALING WITH DATA
	if(part_num == 1){
		for(int ientry=0; ientry<d_entries; ientry++){
			d_tree->GetEntry(ientry);
			if(d_det_num > (maxAD-1)) continue;
			h_d_singles_energy_raw[d_det_num]->Fill(d_energy);
		}

		for(int ientry=0; ientry<p_entries; ientry++){
			p_tree->GetEntry(ientry);
			if(det_num > (maxAD-1)) continue;
			h_p_singles_energy_raw[det_num]->Fill(s_energy);
		}
	}



	//Random Funtion:
	TRandom3* random = new TRandom3();

	//Making tr_pair:
//For pairing prompt to each delayed:
	if(d_entries <= p_entries){
		cout << "Pairing a prompt to each delayed event" << endl;
		for(int ientry=int(float((part_num-1.)/3.)*d_entries); ientry<int(float((part_num)/3.)*d_entries); ientry++){
			d_x = 0;
			d_y = 0;
			d_z = 0;
			if(ientry%10000 ==0) cout<< ientry << "/" << d_entries << endl;
			d_tree->GetEntry(ientry);
			if(d_energy < dElow || d_energy > dEhigh) continue;

			p_rand_entry =(int) random->Uniform(0,p_entries);
	//		p_rand_entry = ientry;
			p_x = 0;
			p_y = 0;
			p_z = 0;
			p_tree->GetEntry(p_rand_entry);

			//Making sure it's in the same detector
			while(d_hall_num != hall_num || d_det_num != det_num || s_energy < pElow || s_energy > pEhigh){
				p_rand_entry =(int) random->Uniform(0,p_entries);
		//		p_rand_entry+=1;
				p_tree->GetEntry(p_rand_entry);
			}
			
			//Renaming prompt variables to save to new TTree
			p_x = s_x;
			p_y = s_y;
			p_z = s_z;
			p_energy = s_energy;
			distance = sqrt(pow(p_x-d_x,2)+pow(p_y-d_y,2)+pow(p_z-d_z,2));
			tr_pair->Fill();
		}
	}

//For pairing each delayed to a prompt:
	else{
		cout << "Pairing a delayed to each prompt event" << endl;
		for(int ientry=int(float((part_num-1.)/3.)*p_entries); ientry<int(float((part_num)/3.)*p_entries); ientry++){
			p_x = 0;
			p_y = 0;
			p_z = 0;
			if(ientry%10000 ==0) cout<< ientry << "/" << p_entries << endl;
			p_tree->GetEntry(ientry);
			p_x = s_x;
			p_y = s_y;
			p_z = s_z;
			p_energy = s_energy;
			p_hall_num = hall_num;
			p_det_num = det_num;
			if(p_energy < pElow || p_energy > pEhigh) continue;

			d_rand_entry =(int) random->Uniform(0,d_entries);
	//		d_rand_entry = ientry;
			d_x = 0;
			d_y = 0;
			d_z = 0;
			d_tree->GetEntry(d_rand_entry);

			//Making sure it's in the same detector
			while(p_hall_num != d_hall_num || p_det_num != d_det_num || d_energy < dElow || d_energy > dEhigh){
				d_rand_entry =(int) random->Uniform(0,d_entries);
		//		d_rand_entry+=1;
				d_tree->GetEntry(d_rand_entry);
			}

			distance = sqrt(pow(p_x-d_x,2)+pow(p_y-d_y,2)+pow(p_z-d_z,2));
			tr_pair->Fill();
		}
	}

	cout << "Number of pairs: " << tr_pair->GetEntries() << endl << endl << endl << endl;

	//SAVING TO OUTPUT FILE
        char* outputname;
//        sprintf(outputname,"./accResults/EH%d/AccidentalsSummary_%d_%d_part%d.root",EH,pd_window_microsec,run_num,part_num);
//	if(dElow >1.5) sprintf(outputname,"./accResults/EH%d/3sig/AccidentalsSummary_%d_%d_part%d.root",EH,pd_window_microsec,run_num,part_num);
//        sprintf(outputname,"./accResults/EH%d/AccidentalsSummary_TcLong_%d_part%d.root",EH,run_num,part_num);
//        sprintf(outputname,"./accResults/EH%d/practiceSummary_%d.root",EH,run_num);
//        sprintf(outputname,"./AccidentalsSummary_%d.root",run_num);

	if(part_num ==1) outputname = outputname_summarize_long1;
	if(part_num ==2) outputname = outputname_summarize_long2;
	if(part_num ==3) outputname = outputname_summarize_long3;

	cout << "Saving to: " << outputname << endl;
	TFile* outfile=new TFile(outputname, "RECREATE");
		outfile->cd();
		tr_live->Write();
		tr_pair->Write();

		if(part_num == 1){
			for(int iad = 0; iad < maxAD; ++iad){
			//	h_p_singles_energy_raw[iad]->SetStats(0);
				h_p_singles_energy_raw[iad]->GetXaxis()->SetTitle("Prompt Singles Energy [MeV]");
				h_p_singles_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_raw[iad]->Write();

			//	h_d_singles_energy_raw[iad]->SetStats(0);
				h_d_singles_energy_raw[iad]->GetXaxis()->SetTitle("Delayed Singles Energy [MeV]");
				h_d_singles_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_raw[iad]->Write();
			}
		}

	outfile->Close();

}







void plot(int run_num, float pElow, float pEhigh, float dElow, float dEhigh, int pd_window_microsec){ //Taking the data from above and making plots of it

	int EH = 0;
	int runNum = 0;
	int maxAD = 0;
    FILE* runfile=fopen("./run_list_good_sync.txt","r");  //This file contains all good runs numbers and experimental hall they belong to
    while(1){ //go though the file and retrieve EH which corresponds to run_num provided as a parameter
        fscanf(runfile,"%d %d",&runNum,&EH);
        if(feof(runfile)) break;
        if(run_num == runNum) break;
    }
    fclose(runfile);
    
    cout<<"Run identified as: run#"<<run_num<<" EH"<<EH<<endl;

	if(EH == 1 || EH == 2) maxAD = 2;
	if(EH == 3) maxAD = 4;

//Hard-coded: Delayed Energy 3 sigma cut for Prompt Energy Plots
//	double peak_Ed[8] = {2.25265, 2.25556, 2.25893, 2.25995, 2.26061, 2.26219, 2.25559, 2.26959};
//	double sigma_Ed[8] = {0.136083, 0.137839, 0.135009, 0.134626, 0.135502, 0.134837, 0.136808, 0.134884};

//For NU Results
	double peak_Ed[8] = {2.26019, 2.26252, 2.26865, 2.26805, 2.26442, 2.26818, 2.26226, 2.27394};
	double sigma_Ed[8] = {0.136494, 0.137613, 0.135735, 0.135387, 0.136701, 0.134367, 0.135997, 0.135052};

	const int NzBins = 20;
	double nZbins = NzBins;

	const int Nr2Bins = 20;
	double nR2bins = Nr2Bins;

	//Making Histograms and TTree:
		TH2F* h_acc_energy_before[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_acc_energy_before_ad%d",iad+1);
		//	h_acc_energy_before[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_acc_energy_before[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH2F* h_acc_energy_1m[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_acc_energy_1m_ad%d",iad+1);
		//	h_acc_energy_1m[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_acc_energy_1m[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH1F* h_p_singles_energy[maxAD]; //energy histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_ad%d",iad+1);
			h_p_singles_energy[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_p_singles_energy_DT800[maxAD]; //energy histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_DT800_ad%d",iad+1);
			h_p_singles_energy_DT800[iad]=new TH1F(name,name,113,0.7,12.);
		}

			TH1F* h_p_singles_energy_DT800_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT800_3sig_ad%d", iad+1);
				h_p_singles_energy_DT800_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT300_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT300_3sig_ad%d", iad+1);
				h_p_singles_energy_DT300_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT500_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT500_3sig_ad%d", iad+1);
				h_p_singles_energy_DT500_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT1000_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT1000_3sig_ad%d", iad+1);
				h_p_singles_energy_DT1000_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT1500_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT1500_3sig_ad%d", iad+1);
				h_p_singles_energy_DT1500_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT2000_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT2000_3sig_ad%d", iad+1);
				h_p_singles_energy_DT2000_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

		TH1F* h_p_singles_energy_raw[maxAD]; //energy histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_raw_ad%d",iad+1);
			h_p_singles_energy_raw[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_p_singles_energy_largeDist[maxAD]; //energy histogram for prompt with less than 1.5 meters separation
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_largeDist_ad%d",iad+1);
			h_p_singles_energy_largeDist[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_p_singles_energy_shortDist[maxAD]; //energy histogram for prompt with less than half a meter separation
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_shortDist_ad%d",iad+1);
			h_p_singles_energy_shortDist[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH2F* h_p_singles_locations[maxAD]; //where the events are happening in the AD histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_locations_ad%d",iad+1);
			h_p_singles_locations[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH1F* h_d_singles_energy_raw[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_raw_ad%d",iad+1);
		//	h_d_singles_energy_raw[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_raw[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_d_singles_energy[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_ad%d",iad+1);
		//	h_d_singles_energy[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_d_singles_energy_DT800[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_DT800_ad%d",iad+1);
			h_d_singles_energy_DT800[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_d_singles_energy_fine[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_fine_ad%d",iad+1);
		//	h_d_singles_energy[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_fine[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_d_singles_energy_fine_DT800[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_fine_DT800_ad%d",iad+1);
			h_d_singles_energy_fine_DT800[iad]=new TH1F(name,name,23000,0.7,3.);
		}

			TH1F* h_d_singles_energy_fine_DT300[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT300_ad%d",iad+1);
				h_d_singles_energy_fine_DT300[iad]=new TH1F(name,name,23000,0.7,3.);
			}

			TH1F* h_d_singles_energy_fine_DT500[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT500_ad%d",iad+1);
				h_d_singles_energy_fine_DT500[iad]=new TH1F(name,name,23000,0.7,3.);
			}

			TH1F* h_d_singles_energy_fine_DT1000[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT1000_ad%d",iad+1);
				h_d_singles_energy_fine_DT1000[iad]=new TH1F(name,name,23000,0.7,3.);
			}

			TH1F* h_d_singles_energy_fine_DT1500[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT1500_ad%d",iad+1);
				h_d_singles_energy_fine_DT1500[iad]=new TH1F(name,name,23000,0.7,3.);
			}

			TH1F* h_d_singles_energy_fine_DT2000[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT2000_ad%d",iad+1);
				h_d_singles_energy_fine_DT2000[iad]=new TH1F(name,name,23000,0.7,3.);
			}

		TH1F* h_delayed_energy_scaled_Ep35[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_delayed_energy_scaled_Ep35_ad%d",iad+1);
		//	h_d_singles_energy[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_Ep35[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_DT800_Ep35[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_delayed_energy_scaled_DT800_Ep35_ad%d",iad+1);
			h_delayed_energy_scaled_DT800_Ep35[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_d_singles_energy_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				char name[64];
				sprintf(name,"h_d_singles_energy_z_ad%d_iz%d", iad+1, iz+1);
				h_d_singles_energy_z[iad][iz]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_d_singles_energy_DT800_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				char name[64];
				sprintf(name,"h_d_singles_energy_DT800_z_ad%d_iz%d", iad+1, iz+1);
				h_d_singles_energy_DT800_z[iad][iz]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_d_singles_energy_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				char name[64];
				sprintf(name,"h_d_singles_energy_r2_ad%d_ir2%d", iad+1, ir2+1);
				h_d_singles_energy_r2[iad][ir2]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_d_singles_energy_DT800_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				char name[64];
				sprintf(name,"h_d_singles_energy_DT800_r2_ad%d_ir2%d", iad+1, ir2+1);
				h_d_singles_energy_DT800_r2[iad][ir2]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_d_singles_energy_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					char name[64];
					sprintf(name,"h_d_singles_energy_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1,iz+1);
					h_d_singles_energy_zVSr2[iad][ir2][iz]=new TH1F(name,name,230,0.7,3.);
				}
			}
		}

		TH1F* h_d_singles_energy_DT800_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					char name[64];
					sprintf(name,"h_d_singles_energy_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1);
					h_d_singles_energy_DT800_zVSr2[iad][ir2][iz]=new TH1F(name,name,230,0.7,3.);
				}
			}
		}

			TH2F* h_acc_delayedVsDT[maxAD]; //delayed vs. DT histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_delayedVsDT_ad%d",iad+1);
				h_acc_delayedVsDT[iad]=new TH2F(name,name,500,0,10,2300,0.7,3.);
			}

			TH2F* h_acc_delayedVsDT_Ep35[maxAD]; //delayed vs. DT histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_delayedVsDT_Ep35_ad%d",iad+1);
				h_acc_delayedVsDT_Ep35[iad]=new TH2F(name,name,500,0,10,2300,0.7,3.);
			}


		TH1F* h_d_singles_energy_largeDist[maxAD]; //energy histogram for delayed with less than 1.5 meters in separation
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_largeDist_ad%d",iad+1);
		//	h_d_singles_energy_largeDist[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_largeDist[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_d_singles_energy_shortDist[maxAD]; //energy histogram for delayed with less than half a meter in separation
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_shortDist_ad%d",iad+1);
		//	h_d_singles_energy_shortDist[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_shortDist[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH2F* h_d_singles_locations[maxAD]; //where the events are happening in the AD histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_locations_ad%d",iad+1);
			h_d_singles_locations[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH1F* h_acc_distance[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_distance_ad%d", iad+1);
			h_acc_distance[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_acc_distance_3sig[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_distance_3sig_ad%d", iad+1);
			h_acc_distance_3sig[iad]=new TH1F(name,name,700,0,7.);
		}

//Delayed energy plots for various distance ranges
		TH1F* h_delayed_energy_scaled_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_delayed_energy_scaled_dist%d_ad%d",dist, iad+1);
			//	h_delayed_energy_scaled_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_delayed_energy_scaled_dist[iad][dist]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_prompt_energy_scaled_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_prompt_energy_scaled_dist%d_ad%d",dist, iad+1);
			//	h_prompt_energy_scaled_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_prompt_energy_scaled_dist[iad][dist]=new TH1F(name,name,113,0.7,12.);
			}
		}

		TH1F* h_delayed_energy_scaled_p35[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_p35_ad%d", iad+1);
		//	h_delayed_energy_scaled_p35[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_p35[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_15_22[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_15_22_ad%d", iad+1);
		//	h_delayed_energy_scaled_15_22[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_15_22[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_22_29[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_22_29_ad%d", iad+1);
		//	h_delayed_energy_scaled_22_29[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_22_29[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_29_36[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_29_36_ad%d", iad+1);
		//	h_delayed_energy_scaled_29_36[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_29_36[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_36_43[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_36_43_ad%d", iad+1);
		//	h_delayed_energy_scaled_36_43[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_36_43[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_43_50[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_43_50_ad%d", iad+1);
		//	h_delayed_energy_scaled_43_50[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_43_50[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH2F* h_delayed_energy_vs_distance[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_vs_distance_ad%d", iad+1);
		//	h_delayed_energy_vs_distance[iad]=new TH2F(name,name,150,1.5,3.,10,0,5.);
			h_delayed_energy_vs_distance[iad]=new TH2F(name,name,230,0.7,3.,11,0,5.5);
		}

		TH2F* h_prompt_energy_vs_distance[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_prompt_energy_vs_distance_ad%d", iad+1);
		//	h_prompt_energy_vs_distance[iad]=new TH2F(name,name,150,1.5,3.,10,0,5.);
			h_prompt_energy_vs_distance[iad]=new TH2F(name,name,113,0.7,12.,11,0,5.5);
		}

			TH2F* h_acc_distVStime[maxAD];
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_distVStime_ad%d",iad+1);
				h_acc_distVStime[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

			TH2F* h_acc_distVStime_Ep35[maxAD];
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_distVStime_Ep35_ad%d",iad+1);
				h_acc_distVStime_Ep35[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

			TH2F* h_acc_promptVStime[maxAD]; //prompt vs. time histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_promptVStime_ad%d",iad+1);
				h_acc_promptVStime[iad]=new TH2F(name,name,1999,1,2000,113,0.7,12.);
			}

			TH2F* h_acc_promptVStime_DT800[maxAD]; //prompt vs. time histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_promptVStime_DT800_ad%d",iad+1);
				h_acc_promptVStime_DT800[iad]=new TH2F(name,name,1999,1,2000,113,0.7,12.);
			}

		TH1F* h_acc_distance_Ep35[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_distance_Ep35_ad%d", iad+1);
			h_acc_distance_Ep35[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1D* h_acc_DT[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_DT_ad%d", iad+1);
			h_acc_DT[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_acc_DT_3sig[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_DT_3sig_ad%d", iad+1);
			h_acc_DT_3sig[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_acc_DT_Ep35[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_DT_Ep35_ad%d", iad+1);
			h_acc_DT_Ep35[iad]=new TH1D(name,name,500,0,10);
		}
//End of energy plots for different distance requirements


	double p_DAQ_ad1 = 0;
	double p_veto_ad1 = 0;
	double p_live_ad1 = 0;
	double p_DAQ_ad2 = 0;
	double p_veto_ad2 = 0;
	double p_live_ad2 = 0;
	double p_DAQ_ad3 = 0;
	double p_veto_ad3 = 0;
	double p_live_ad3 = 0;
	double p_DAQ_ad4 = 0;
	double p_veto_ad4 = 0;
	double p_live_ad4 = 0;
	double p_counts = 0;
	double d_DAQ_ad1 = 0;
	double d_veto_ad1 = 0;
	double d_live_ad1 = 0;
	double d_DAQ_ad2 = 0;
	double d_veto_ad2 = 0;
	double d_live_ad2 = 0;
	double d_DAQ_ad3 = 0;
	double d_veto_ad3 = 0;
	double d_live_ad3 = 0;
	double d_DAQ_ad4 = 0;
	double d_veto_ad4 = 0;
	double d_live_ad4 = 0;
	double d_counts = 0;


		TTree* tr_live=new TTree("tr_accidentals_time", "Summary of Singles Search");
			tr_live->Branch("p_DAQ_ad1", &p_DAQ_ad1, "DAQ for AD1/Double");
			tr_live->Branch("p_veto_ad1", &p_veto_ad1, "Veto for AD1/Double");
			tr_live->Branch("p_live_ad1", &p_live_ad1, "Live for AD1/Double");
			tr_live->Branch("p_DAQ_ad2", &p_DAQ_ad2, "DAQ for AD2/Double");
			tr_live->Branch("p_veto_ad2", &p_veto_ad2, "Veto for AD2/Double");
			tr_live->Branch("p_live_ad2", &p_live_ad2, "Live for AD2/Double");
			tr_live->Branch("p_DAQ_ad3", &p_DAQ_ad3, "DAQ for AD3/Double");
			tr_live->Branch("p_veto_ad3", &p_veto_ad3, "Veto for AD3/Double");
			tr_live->Branch("p_live_ad3", &p_live_ad3, "Live for AD3/Double");
			tr_live->Branch("p_DAQ_ad4", &p_DAQ_ad4, "DAQ for AD4/Double");
			tr_live->Branch("p_veto_ad4", &p_veto_ad4, "Veto for AD4/Double");
			tr_live->Branch("p_live_ad4", &p_live_ad4, "Live for AD4/Double");
			tr_live->Branch("d_DAQ_ad1", &d_DAQ_ad1, "DAQ for AD1/Double");
			tr_live->Branch("d_veto_ad1", &d_veto_ad1, "Veto for AD1/Double");
			tr_live->Branch("d_live_ad1", &d_live_ad1, "Live for AD1/Double");
			tr_live->Branch("d_DAQ_ad2", &d_DAQ_ad2, "DAQ for AD2/Double");
			tr_live->Branch("d_veto_ad2", &d_veto_ad2, "Veto for AD2/Double");
			tr_live->Branch("d_live_ad2", &d_live_ad2, "Live for AD2/Double");
			tr_live->Branch("d_DAQ_ad3", &d_DAQ_ad3, "DAQ for AD3/Double");
			tr_live->Branch("d_veto_ad3", &d_veto_ad3, "Veto for AD3/Double");
			tr_live->Branch("d_live_ad3", &d_live_ad3, "Live for AD3/Double");
			tr_live->Branch("d_DAQ_ad4", &d_DAQ_ad4, "DAQ for AD4/Double");
			tr_live->Branch("d_veto_ad4", &d_veto_ad4, "Veto for AD4/Double");
			tr_live->Branch("d_live_ad4", &d_live_ad4, "Live for AD4/Double");


	float p_energy = 0;
	float d_energy = 0;
	float p_x = 0;
	float p_y = 0;
	float p_z = 0;
	float d_x = 0;
	float d_y = 0;
	float d_z = 0;
	float distance = 0;
	int det_num = 0;
	double ranTime = 0;
	double DTvalue = 0;

	int pair_entries = 0;

	//Reading in the accidentals file
//        char inputname[64];
//        sprintf(inputname,"./accResults/EH%d/AccidentalsSummary_%d_%d.root",EH,pd_window_microsec,run_num);
//		if(dElow > 1.5) sprintf(inputname,"./accResults/EH%d/3sig/AccidentalsSummary_%d_%d.root",EH,pd_window_microsec,run_num);
//        sprintf(inputname,"./accResults/EH%d/AccidentalsSummary_TcLong_%d.root",EH,run_num);
	cout << "Retrieving: " << outputname_summarize << endl;
	TFile *accInfile = new TFile(outputname_summarize,"READ");

	for(int iad = 0; iad < maxAD; ++iad){
		char name[64];
		sprintf(name,"h_d_singles_energy_raw_ad%d",iad+1);
		TH1F *h_d_raw = (TH1F*)accInfile->Get(name);
		h_d_singles_energy_raw[iad]->Add(h_d_raw);

		sprintf(name,"h_p_singles_energy_raw_ad%d",iad+1);
		TH1F *h_p_raw = (TH1F*)accInfile->Get(name);
		h_p_singles_energy_raw[iad]->Add(h_p_raw);
	} 

	//Getting the pair tree
	TTree *tr_pair = (TTree*)accInfile->Get("tr_pair");
	tr_pair->SetBranchStatus("*",1);

		//Establishing the new TTree of randomly paired singles
		tr_pair->SetBranchAddress("p_energy", &p_energy);
		tr_pair->SetBranchAddress("d_energy", &d_energy);
		tr_pair->SetBranchAddress("p_x", &p_x);
		tr_pair->SetBranchAddress("d_x", &d_x);
		tr_pair->SetBranchAddress("p_y", &p_y);
		tr_pair->SetBranchAddress("d_y", &d_y);
		tr_pair->SetBranchAddress("p_z", &p_z);
		tr_pair->SetBranchAddress("d_z", &d_z);
		tr_pair->SetBranchAddress("det_num", &det_num);
		tr_pair->SetBranchAddress("distance", &distance);


		//TTree tr_sum
		TTree *tr_time = (TTree*)accInfile->Get("tr_singles_summary");
		tr_time->SetBranchStatus("*",1);


			tr_time->SetBranchAddress("p_DAQ_ad1", &p_DAQ_ad1);
			tr_time->SetBranchAddress("p_veto_ad1", &p_veto_ad1);
			tr_time->SetBranchAddress("p_live_ad1", &p_live_ad1);
			tr_time->SetBranchAddress("p_DAQ_ad2", &p_DAQ_ad2);
			tr_time->SetBranchAddress("p_veto_ad2", &p_veto_ad2);
			tr_time->SetBranchAddress("p_live_ad2", &p_live_ad2);
			tr_time->SetBranchAddress("p_DAQ_ad3", &p_DAQ_ad3);
			tr_time->SetBranchAddress("p_veto_ad3", &p_veto_ad3);
			tr_time->SetBranchAddress("p_live_ad3", &p_live_ad3);
			tr_time->SetBranchAddress("p_DAQ_ad4", &p_DAQ_ad4);
			tr_time->SetBranchAddress("p_veto_ad4", &p_veto_ad4);
			tr_time->SetBranchAddress("p_live_ad4", &p_live_ad4);
			tr_time->SetBranchAddress("d_DAQ_ad1", &d_DAQ_ad1);
			tr_time->SetBranchAddress("d_veto_ad1", &d_veto_ad1);
			tr_time->SetBranchAddress("d_live_ad1", &d_live_ad1);
			tr_time->SetBranchAddress("d_DAQ_ad2", &d_DAQ_ad2);
			tr_time->SetBranchAddress("d_veto_ad2", &d_veto_ad2);
			tr_time->SetBranchAddress("d_live_ad2", &d_live_ad2);
			tr_time->SetBranchAddress("d_DAQ_ad3", &d_DAQ_ad3);
			tr_time->SetBranchAddress("d_veto_ad3", &d_veto_ad3);
			tr_time->SetBranchAddress("d_live_ad3", &d_live_ad3);
			tr_time->SetBranchAddress("d_DAQ_ad4", &d_DAQ_ad4);
			tr_time->SetBranchAddress("d_veto_ad4", &d_veto_ad4);
			tr_time->SetBranchAddress("d_live_ad4", &d_live_ad4);




	int nentries = tr_time->GetEntries();

	cout << "nentries = " << nentries << endl;

	for(int ientry=0; ientry<nentries; ientry++){
		tr_time->GetEntry(0);
		tr_live->Fill();
	}

	cout << fixed;
	cout << "Delayed live: " << setprecision(9) << d_live_ad1 << endl;
	cout << "Prompt live: " << setprecision(9) << p_live_ad1 << endl;

	pair_entries = tr_pair->GetEntries();

	//Random Funtion:
	TRandom3* random = new TRandom3();

	cout << "Num Pair_entries: " << pair_entries << endl;

	//Plotting tr_pair:
	for(int ientry=0; ientry<pair_entries; ientry++){
		tr_pair->GetEntry(ientry);
		if(p_energy < pElow || p_energy > pEhigh) continue;
		if(d_energy < dElow || d_energy > dEhigh) continue;

		ranTime = (double) random->Uniform(1,pd_window_microsec);
		DTvalue = distance*1.e-3+(ranTime)/600;

	//	if(DTvalue >= 0.8) continue;

			h_acc_distVStime[det_num]->Fill(ranTime,distance*1.e-3);
			if(p_energy >=3.5){
				h_acc_distVStime_Ep35[det_num]->Fill(ranTime,distance*1.e-3);
				h_acc_distance_Ep35[det_num]->Fill(distance*1.e-3);
				h_acc_DT_Ep35[det_num]->Fill(DTvalue);
				h_acc_delayedVsDT_Ep35[det_num]->Fill(DTvalue, d_energy);
			}

			h_acc_DT[det_num]->Fill(DTvalue);

			h_acc_promptVStime[det_num]->Fill(ranTime, p_energy);

			int temp_iZ = int((d_z*1.e-3+2.5)/(5./nZbins));
			int temp_iR2 = int(((pow(p_x*1.e-3,2.)+pow(p_y*1.e-3,2.)))/(5./nR2bins));


		if(DTvalue < .8){
			h_d_singles_energy_DT800[det_num]->Fill(d_energy);
			if(temp_iZ < nZbins && temp_iZ >= 0) h_d_singles_energy_DT800_z[det_num][temp_iZ]->Fill(d_energy);
			if(temp_iR2 < nR2bins && temp_iR2 >= 0) h_d_singles_energy_DT800_r2[det_num][temp_iR2]->Fill(d_energy);
			if(temp_iR2 < nR2bins && temp_iR2 >= 0 && temp_iZ < nZbins && temp_iZ >= 0 && p_energy >= 3.5) h_d_singles_energy_DT800_zVSr2[det_num][temp_iR2][temp_iZ]->Fill(d_energy);
			h_d_singles_energy_fine_DT800[det_num]->Fill(d_energy);
			if(p_energy >= 3.5) h_delayed_energy_scaled_DT800_Ep35[det_num]->Fill(d_energy);
			h_p_singles_energy_DT800[det_num]->Fill(p_energy);
			h_acc_promptVStime_DT800[det_num]->Fill(ranTime, p_energy);
			if(d_energy >= (peak_Ed[2*(EH-1)+det_num] - 3*sigma_Ed[2*(EH-1)+det_num]) && d_energy <= (peak_Ed[2*(EH-1)+det_num] + 3*sigma_Ed[2*(EH-1)+det_num])) h_p_singles_energy_DT800_3sig[det_num]->Fill(p_energy);
		}

		if(d_energy >= (peak_Ed[2*(EH-1)+det_num] - 3*sigma_Ed[2*(EH-1)+det_num]) && d_energy <= (peak_Ed[2*(EH-1)+det_num] + 3*sigma_Ed[2*(EH-1)+det_num])){
			if(DTvalue < 0.3) h_p_singles_energy_DT300_3sig[det_num]->Fill(p_energy);
			if(DTvalue < 0.5) h_p_singles_energy_DT500_3sig[det_num]->Fill(p_energy);
			if(DTvalue < 1.0) h_p_singles_energy_DT1000_3sig[det_num]->Fill(p_energy);
			if(DTvalue < 1.5) h_p_singles_energy_DT1500_3sig[det_num]->Fill(p_energy);
			if(DTvalue < 2.0) h_p_singles_energy_DT2000_3sig[det_num]->Fill(p_energy);
		}

		if(DTvalue < 0.3) h_d_singles_energy_fine_DT300[det_num]->Fill(d_energy);
		if(DTvalue < 0.5) h_d_singles_energy_fine_DT500[det_num]->Fill(d_energy);
		if(DTvalue < 1.0) h_d_singles_energy_fine_DT1000[det_num]->Fill(d_energy);
		if(DTvalue < 1.5) h_d_singles_energy_fine_DT1500[det_num]->Fill(d_energy);
		if(DTvalue < 2.0) h_d_singles_energy_fine_DT2000[det_num]->Fill(d_energy);
		
		if(d_energy >= (peak_Ed[2*(EH-1)+det_num] - 3*sigma_Ed[2*(EH-1)+det_num]) && d_energy <= (peak_Ed[2*(EH-1)+det_num] + 3*sigma_Ed[2*(EH-1)+det_num])){
			h_acc_distance_3sig[det_num]->Fill(distance*1.e-3);
			h_acc_DT_3sig[det_num]->Fill(DTvalue);
		}

		h_d_singles_energy[det_num]->Fill(d_energy);
		h_d_singles_energy_fine[det_num]->Fill(d_energy);
		h_acc_delayedVsDT[det_num]->Fill(DTvalue, d_energy);
		if(temp_iZ < nZbins && temp_iZ >= 0) h_d_singles_energy_z[det_num][temp_iZ]->Fill(d_energy);
		if(temp_iR2 < nR2bins && temp_iR2 >= 0) h_d_singles_energy_r2[det_num][temp_iR2]->Fill(d_energy);
		if(temp_iR2 < nR2bins && temp_iR2 >= 0 && temp_iZ < nZbins && temp_iZ >= 0 && p_energy >= 3.5) h_d_singles_energy_zVSr2[det_num][temp_iR2][temp_iZ]->Fill(d_energy);
		h_d_singles_locations[det_num]->Fill((pow(d_x*1.e-3,2.)+pow(d_y*1.e-3,2.)),d_z*1.e-3);

		if(distance*1.e-3 > 5.0)h_delayed_energy_vs_distance[det_num]->Fill(d_energy, 5.25);
		else h_delayed_energy_vs_distance[det_num]->Fill(d_energy, distance*1.e-3);

		if(distance*1.e-3 > 5.0)h_prompt_energy_vs_distance[det_num]->Fill(p_energy, 5.25);
		else h_prompt_energy_vs_distance[det_num]->Fill(p_energy, distance*1.e-3);

		h_acc_energy_before[det_num]->Fill(p_energy, d_energy);

		h_p_singles_energy[det_num]->Fill(p_energy);
		h_p_singles_locations[det_num]->Fill((pow(p_x*1.e-3,2.)+pow(p_y*1.e-3,2.)),p_z*1.e-3);
		h_acc_distance[det_num]->Fill(distance*1.e-3);

		if(p_energy >= 3.5){
			h_delayed_energy_scaled_p35[det_num]->Fill(d_energy);
			h_delayed_energy_scaled_Ep35[det_num]->Fill(d_energy);
		}

		if(distance*1.e-3 > 1.5 && distance*1.e-3 <= 2.2) h_delayed_energy_scaled_15_22[det_num]->Fill(d_energy);
		if(distance*1.e-3 > 2.2 && distance*1.e-3 <= 2.9) h_delayed_energy_scaled_22_29[det_num]->Fill(d_energy);
		if(distance*1.e-3 > 2.9 && distance*1.e-3 <= 3.6) h_delayed_energy_scaled_29_36[det_num]->Fill(d_energy);
		if(distance*1.e-3 > 3.6 && distance*1.e-3 <= 4.3) h_delayed_energy_scaled_36_43[det_num]->Fill(d_energy);
		if(distance*1.e-3 > 4.3) h_delayed_energy_scaled_43_50[det_num]->Fill(d_energy);

		if(distance*1.e-3 <= 1.0){
			h_delayed_energy_scaled_dist[det_num][0]->Fill(d_energy);
			h_prompt_energy_scaled_dist[det_num][0]->Fill(p_energy);
			h_acc_energy_1m[det_num]->Fill(p_energy, d_energy);
		}
		if(distance*1.e-3 > 1.0 && distance*1.e-3 <= 2.0){
			h_delayed_energy_scaled_dist[det_num][1]->Fill(d_energy);
			h_prompt_energy_scaled_dist[det_num][1]->Fill(p_energy);
		}
		if(distance*1.e-3 > 2.0 && distance*1.e-3 <= 3.0){
			h_delayed_energy_scaled_dist[det_num][2]->Fill(d_energy);
			h_prompt_energy_scaled_dist[det_num][2]->Fill(p_energy);
		}
		if(distance*1.e-3 > 3.0 && distance*1.e-3 <= 4.0){
			h_delayed_energy_scaled_dist[det_num][3]->Fill(d_energy);
			h_prompt_energy_scaled_dist[det_num][3]->Fill(p_energy);
		}
		if(distance*1.e-3 > 4.0 && distance*1.e-3 <= 5.0){
			h_delayed_energy_scaled_dist[det_num][4]->Fill(d_energy);
			h_prompt_energy_scaled_dist[det_num][4]->Fill(p_energy);
		}

		if(distance*1.e-3 > 5.){
			h_delayed_energy_scaled_dist[det_num][5]->Fill(d_energy);
			h_prompt_energy_scaled_dist[det_num][5]->Fill(p_energy);
		}

		if(distance*1.e-3 > 1.5) continue;

		h_p_singles_energy_largeDist[det_num]->Fill(p_energy);
		h_d_singles_energy_largeDist[det_num]->Fill(d_energy);

		if(distance*1.e-3 > 0.5) continue;

		h_p_singles_energy_shortDist[det_num]->Fill(p_energy);
		h_d_singles_energy_shortDist[det_num]->Fill(d_energy);

	}



	//SAVING TO OUTPUT FILE
//        char outputname[64];
//        sprintf(outputname,"./accResults/EH%d/AccidentalsPlots_%d_%d.root",EH,pd_window_microsec,run_num);
//	if(dElow >1.5) sprintf(outputname,"./accResults/EH%d/3sig/AccidentalsPlots_%d_%d.root",EH,pd_window_microsec,run_num);
//        sprintf(outputname,"./accResults/EH%d/AccidentalsPlots_TcLong_%d.root",EH,run_num);

	cout << "Saving to: " << outputname_plot << endl;
//        sprintf(outputname,"./AccidentalsSummary_%d.root",run_num);
	TFile* outfile=new TFile(outputname_plot, "RECREATE");
		outfile->cd();
		tr_live->Write();
		for(int iad=0; iad<maxAD; ++iad){
			h_acc_energy_before[iad]->SetOption("COLZ");
			h_acc_energy_before[iad]->SetStats(0);
			h_acc_energy_before[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_acc_energy_before[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_acc_energy_before[iad]->Write();

			h_acc_energy_1m[iad]->SetOption("COLZ");
			h_acc_energy_1m[iad]->SetStats(0);
			h_acc_energy_1m[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_acc_energy_1m[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_acc_energy_1m[iad]->Write();

			h_acc_distVStime[iad]->SetStats(0);
			h_acc_distVStime[iad]->SetOption("COLZ");
			h_acc_distVStime[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_acc_distVStime[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_acc_distVStime[iad]->Write();

			h_acc_distVStime_Ep35[iad]->SetStats(0);
			h_acc_distVStime_Ep35[iad]->SetOption("COLZ");
			h_acc_distVStime_Ep35[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_acc_distVStime_Ep35[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_acc_distVStime_Ep35[iad]->Write();

		//		h_acc_DT[iad]->SetStats(0);
				h_acc_DT[iad]->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT[iad]->GetYaxis()->SetTitle("Counts");
				h_acc_DT[iad]->Write();

		//		h_acc_DT_3sig[iad]->SetStats(0);
				h_acc_DT_3sig[iad]->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_acc_DT_3sig[iad]->Write();

		//		h_acc_DT_Ep35[iad]->SetStats(0);
				h_acc_DT_Ep35[iad]->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_Ep35[iad]->GetYaxis()->SetTitle("Counts");
				h_acc_DT_Ep35[iad]->Write();

			h_p_singles_locations[iad]->SetOption("COLZ");
			h_p_singles_locations[iad]->SetStats(0);
			h_p_singles_locations[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_p_singles_locations[iad]->GetYaxis()->SetTitle("z [m]");
			h_p_singles_locations[iad]->Write();

		//	h_p_singles_energy_raw[iad]->SetStats(0);
			h_p_singles_energy_raw[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_raw[iad]->Write();

		//	h_p_singles_energy[iad]->SetStats(0);
			h_p_singles_energy[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy[iad]->Write();

		//	h_p_singles_energy_DT800[iad]->SetStats(0);
			h_p_singles_energy_DT800[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_DT800[iad]->Write();

		//	h_p_singles_energy_DT800_3sig[iad]->SetStats(0);
			h_p_singles_energy_DT800_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_DT800_3sig[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_DT800_3sig[iad]->Write();

			//	h_p_singles_energy_DT300_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT300_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT300_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT300_3sig[iad]->Write();

			//	h_p_singles_energy_DT500_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT500_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT500_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT500_3sig[iad]->Write();

			//	h_p_singles_energy_DT1000_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT1000_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT1000_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT1000_3sig[iad]->Write();

			//	h_p_singles_energy_DT1500_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT1500_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT1500_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT1500_3sig[iad]->Write();

			//	h_p_singles_energy_DT2000_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT2000_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT2000_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT2000_3sig[iad]->Write();

		//	h_p_singles_energy_largeDist[iad]->SetStats(0);
			h_p_singles_energy_largeDist[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_largeDist[iad]->Write();

		//	h_p_singles_energy_shortDist[iad]->SetStats(0);
			h_p_singles_energy_shortDist[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_shortDist[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_shortDist[iad]->Write();

			h_d_singles_locations[iad]->SetOption("COLZ");
			h_d_singles_locations[iad]->SetStats(0);
			h_d_singles_locations[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_d_singles_locations[iad]->GetYaxis()->SetTitle("z [m]");
			h_d_singles_locations[iad]->Write();

		//	h_d_singles_energy_raw[iad]->SetStats(0);
			h_d_singles_energy_raw[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_raw[iad]->Write();

		//	h_d_singles_energy[iad]->SetStats(0);
			h_d_singles_energy[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy[iad]->Write();


		//	h_d_singles_energy_DT800[iad]->SetStats(0);
			h_d_singles_energy_DT800[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_DT800[iad]->Write();

		//	h_d_singles_energy_fine[iad]->SetStats(0);
			h_d_singles_energy_fine[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_fine[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_fine[iad]->Write();


		//	h_d_singles_energy_fine_DT800[iad]->SetStats(0);
			h_d_singles_energy_fine_DT800[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_fine_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_fine_DT800[iad]->Write();

			//	h_d_singles_energy_fine_DT300[iad]->SetStats(0);
				h_d_singles_energy_fine_DT300[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT300[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT300[iad]->Write();

			//	h_d_singles_energy_fine_DT500[iad]->SetStats(0);
				h_d_singles_energy_fine_DT500[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT500[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT500[iad]->Write();

			//	h_d_singles_energy_fine_DT1000[iad]->SetStats(0);
				h_d_singles_energy_fine_DT1000[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT1000[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT1000[iad]->Write();

			//	h_d_singles_energy_fine_DT1500[iad]->SetStats(0);
				h_d_singles_energy_fine_DT1500[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT1500[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT1500[iad]->Write();

			//	h_d_singles_energy_fine_DT2000[iad]->SetStats(0);
				h_d_singles_energy_fine_DT2000[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT2000[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT2000[iad]->Write();

		//	h_delayed_energy_scaled_Ep35[iad]->SetStats(0);
			h_delayed_energy_scaled_Ep35[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_delayed_energy_scaled_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_Ep35[iad]->Write();

		//	h_delayed_energy_scaled_DT800_Ep35[iad]->SetStats(0);
			h_delayed_energy_scaled_DT800_Ep35[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_delayed_energy_scaled_DT800_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_DT800_Ep35[iad]->Write();

				h_acc_delayedVsDT[iad]->SetStats(0);
				h_acc_delayedVsDT[iad]->SetOption("COLZ");
				h_acc_delayedVsDT[iad]->GetXaxis()->SetTitle("DT Value [m]");
				h_acc_delayedVsDT[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
				h_acc_delayedVsDT[iad]->Write();

				h_acc_delayedVsDT_Ep35[iad]->SetStats(0);
				h_acc_delayedVsDT_Ep35[iad]->SetOption("COLZ");
				h_acc_delayedVsDT_Ep35[iad]->GetXaxis()->SetTitle("DT Value [m]");
				h_acc_delayedVsDT_Ep35[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
				h_acc_delayedVsDT_Ep35[iad]->Write();

			for(int iz = 0; iz<NzBins; iz++){
				h_d_singles_energy_z[iad][iz]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_z[iad][iz]->Write();

				h_d_singles_energy_DT800_z[iad][iz]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_DT800_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_DT800_z[iad][iz]->Write();
			}

			for(int ir2 = 0; ir2<Nr2Bins; ir2++){
				h_d_singles_energy_r2[iad][ir2]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_r2[iad][ir2]->Write();

				h_d_singles_energy_DT800_r2[iad][ir2]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_DT800_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_DT800_r2[iad][ir2]->Write();
			}

			for(int iz = 0; iz<NzBins; iz++){
				for(int ir2 = 0; ir2<Nr2Bins; ir2++){
					h_d_singles_energy_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Singles Energy [MeV]");
					h_d_singles_energy_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_d_singles_energy_zVSr2[iad][ir2][iz]->Write();

					h_d_singles_energy_DT800_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Singles Energy [MeV]");
					h_d_singles_energy_DT800_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_d_singles_energy_DT800_zVSr2[iad][ir2][iz]->Write();
				}
			}

		//	h_d_singles_energy_largeDist[iad]->SetStats(0);
			h_d_singles_energy_largeDist[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_largeDist[iad]->Write();

		//	h_d_singles_energy_shortDist[iad]->SetStats(0);
			h_d_singles_energy_shortDist[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_shortDist[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_shortDist[iad]->Write();

			h_acc_distance[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_acc_distance[iad]->GetYaxis()->SetTitle("Counts");
			h_acc_distance[iad]->Write();

			h_acc_distance_3sig[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_acc_distance_3sig[iad]->GetYaxis()->SetTitle("Counts");
			h_acc_distance_3sig[iad]->Write();

			h_acc_distance_Ep35[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_acc_distance_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_acc_distance_Ep35[iad]->Write();

			h_delayed_energy_scaled_15_22[iad]->SetStats(0);
			h_delayed_energy_scaled_15_22[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_15_22[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_15_22[iad]->Write();

			h_delayed_energy_scaled_22_29[iad]->SetStats(0);
			h_delayed_energy_scaled_22_29[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_22_29[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_22_29[iad]->Write();

			h_delayed_energy_scaled_29_36[iad]->SetStats(0);
			h_delayed_energy_scaled_29_36[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_29_36[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_29_36[iad]->Write();

			h_delayed_energy_scaled_36_43[iad]->SetStats(0);
			h_delayed_energy_scaled_36_43[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_36_43[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_36_43[iad]->Write();

			h_delayed_energy_scaled_43_50[iad]->SetStats(0);
			h_delayed_energy_scaled_43_50[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_43_50[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_43_50[iad]->Write();

			h_delayed_energy_scaled_p35[iad]->SetStats(0);
			h_delayed_energy_scaled_p35[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_p35[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_p35[iad]->Write();

			h_delayed_energy_vs_distance[iad]->SetStats(0);
			h_delayed_energy_vs_distance[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_vs_distance[iad]->GetYaxis()->SetTitle("Distance between prompt and delayed");
			h_delayed_energy_vs_distance[iad]->SetOption("COLZ");
			h_delayed_energy_vs_distance[iad]->Write();

			h_prompt_energy_vs_distance[iad]->SetStats(0);
			h_prompt_energy_vs_distance[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_prompt_energy_vs_distance[iad]->GetYaxis()->SetTitle("Distance between prompt and delayed");
			h_prompt_energy_vs_distance[iad]->SetOption("COLZ");
			h_prompt_energy_vs_distance[iad]->Write();

			h_acc_promptVStime[iad]->SetStats(0);
			h_acc_promptVStime[iad]->SetOption("COLZ");
			h_acc_promptVStime[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_acc_promptVStime[iad]->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_acc_promptVStime[iad]->Write();

			h_acc_promptVStime_DT800[iad]->SetStats(0);
			h_acc_promptVStime_DT800[iad]->SetOption("COLZ");
			h_acc_promptVStime_DT800[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_acc_promptVStime_DT800[iad]->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_acc_promptVStime_DT800[iad]->Write();

			for(int dist=0; dist<6; dist++){
				h_delayed_energy_scaled_dist[iad][dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_delayed_energy_scaled_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_delayed_energy_scaled_dist[iad][dist]->Write();
				h_prompt_energy_scaled_dist[iad][dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_prompt_energy_scaled_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_prompt_energy_scaled_dist[iad][dist]->Write();
			}
		}


	outfile->Close();

}





void plot_long(int run_num, float pElow, float pEhigh, float dElow, float dEhigh, int pd_window_microsec){ //Taking the data from above and making plots of it

	int EH = 0;
	int runNum = 0;
	int maxAD = 0;
    FILE* runfile=fopen("./run_list_good_sync.txt","r");  //This file contains all good runs numbers and experimental hall they belong to
    while(1){ //go though the file and retrieve EH which corresponds to run_num provided as a parameter
        fscanf(runfile,"%d %d",&runNum,&EH);
        if(feof(runfile)) break;
        if(run_num == runNum) break;
    }
    fclose(runfile);
    
    cout<<"Run identified as: run#"<<run_num<<" EH"<<EH<<endl;

	if(EH == 1 || EH == 2) maxAD = 2;
	if(EH == 3) maxAD = 4;

//Hard-coded: Delayed Energy 3 sigma cut for Prompt Energy Plots
//	double peak_Ed[8] = {2.25265, 2.25556, 2.25893, 2.25995, 2.26061, 2.26219, 2.25559, 2.26959};
//	double sigma_Ed[8] = {0.136083, 0.137839, 0.135009, 0.134626, 0.135502, 0.134837, 0.136808, 0.134884};

//For NU Results
	double peak_Ed[8] = {2.26019, 2.26252, 2.26865, 2.26805, 2.26442, 2.26818, 2.26226, 2.27394};
	double sigma_Ed[8] = {0.136494, 0.137613, 0.135735, 0.135387, 0.136701, 0.134367, 0.135997, 0.135052};

	const int NzBins = 20;
	double nZbins = NzBins;

	const int Nr2Bins = 20;
	double nR2bins = Nr2Bins;

	//Making Histograms and TTree:
		TH2F* h_acc_energy_before[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_acc_energy_before_ad%d",iad+1);
		//	h_acc_energy_before[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_acc_energy_before[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH2F* h_acc_energy_1m[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_acc_energy_1m_ad%d",iad+1);
		//	h_acc_energy_1m[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_acc_energy_1m[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH1F* h_p_singles_energy[maxAD]; //energy histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_ad%d",iad+1);
			h_p_singles_energy[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_p_singles_energy_DT800[maxAD]; //energy histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_DT800_ad%d",iad+1);
			h_p_singles_energy_DT800[iad]=new TH1F(name,name,113,0.7,12.);
		}

			TH1F* h_p_singles_energy_DT800_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT800_3sig_ad%d", iad+1);
				h_p_singles_energy_DT800_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT300_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT300_3sig_ad%d", iad+1);
				h_p_singles_energy_DT300_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT500_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT500_3sig_ad%d", iad+1);
				h_p_singles_energy_DT500_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT1000_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT1000_3sig_ad%d", iad+1);
				h_p_singles_energy_DT1000_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT1500_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT1500_3sig_ad%d", iad+1);
				h_p_singles_energy_DT1500_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_p_singles_energy_DT2000_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_p_singles_energy_DT2000_3sig_ad%d", iad+1);
				h_p_singles_energy_DT2000_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

		TH1F* h_p_singles_energy_raw[maxAD]; //energy histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_raw_ad%d",iad+1);
			h_p_singles_energy_raw[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_p_singles_energy_largeDist[maxAD]; //energy histogram for prompt with less than 1.5 meters separation
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_largeDist_ad%d",iad+1);
			h_p_singles_energy_largeDist[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_p_singles_energy_shortDist[maxAD]; //energy histogram for prompt with less than half a meter separation
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_energy_shortDist_ad%d",iad+1);
			h_p_singles_energy_shortDist[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH2F* h_p_singles_locations[maxAD]; //where the events are happening in the AD histogram for prompt
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_p_singles_locations_ad%d",iad+1);
			h_p_singles_locations[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH1F* h_d_singles_energy_raw[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_raw_ad%d",iad+1);
		//	h_d_singles_energy_raw[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_raw[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_d_singles_energy[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_ad%d",iad+1);
		//	h_d_singles_energy[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_d_singles_energy_DT800[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_DT800_ad%d",iad+1);
			h_d_singles_energy_DT800[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_d_singles_energy_fine[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_fine_ad%d",iad+1);
		//	h_d_singles_energy[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_fine[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_d_singles_energy_fine_DT800[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_fine_DT800_ad%d",iad+1);
			h_d_singles_energy_fine_DT800[iad]=new TH1F(name,name,23000,0.7,3.);
		}

			TH1F* h_d_singles_energy_fine_DT300[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT300_ad%d",iad+1);
				h_d_singles_energy_fine_DT300[iad]=new TH1F(name,name,23000,0.7,3.);
			}

			TH1F* h_d_singles_energy_fine_DT500[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT500_ad%d",iad+1);
				h_d_singles_energy_fine_DT500[iad]=new TH1F(name,name,23000,0.7,3.);
			}

			TH1F* h_d_singles_energy_fine_DT1000[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT1000_ad%d",iad+1);
				h_d_singles_energy_fine_DT1000[iad]=new TH1F(name,name,23000,0.7,3.);
			}

			TH1F* h_d_singles_energy_fine_DT1500[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT1500_ad%d",iad+1);
				h_d_singles_energy_fine_DT1500[iad]=new TH1F(name,name,23000,0.7,3.);
			}

			TH1F* h_d_singles_energy_fine_DT2000[maxAD]; //energy histogram for delayed
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_fine_DT2000_ad%d",iad+1);
				h_d_singles_energy_fine_DT2000[iad]=new TH1F(name,name,23000,0.7,3.);
			}

		TH1F* h_delayed_energy_scaled_Ep35[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_delayed_energy_scaled_Ep35_ad%d",iad+1);
		//	h_d_singles_energy[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_Ep35[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_DT800_Ep35[maxAD]; //energy histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_delayed_energy_scaled_DT800_Ep35_ad%d",iad+1);
			h_delayed_energy_scaled_DT800_Ep35[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_d_singles_energy_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				char name[64];
				sprintf(name,"h_d_singles_energy_z_ad%d_iz%d", iad+1, iz+1);
				h_d_singles_energy_z[iad][iz]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_d_singles_energy_DT800_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				char name[64];
				sprintf(name,"h_d_singles_energy_DT800_z_ad%d_iz%d", iad+1, iz+1);
				h_d_singles_energy_DT800_z[iad][iz]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_d_singles_energy_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				char name[64];
				sprintf(name,"h_d_singles_energy_r2_ad%d_ir2%d", iad+1, ir2+1);
				h_d_singles_energy_r2[iad][ir2]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_d_singles_energy_DT800_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				char name[64];
				sprintf(name,"h_d_singles_energy_DT800_r2_ad%d_ir2%d", iad+1, ir2+1);
				h_d_singles_energy_DT800_r2[iad][ir2]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_d_singles_energy_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					char name[64];
					sprintf(name,"h_d_singles_energy_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1,iz+1);
					h_d_singles_energy_zVSr2[iad][ir2][iz]=new TH1F(name,name,230,0.7,3.);
				}
			}
		}

		TH1F* h_d_singles_energy_DT800_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					char name[64];
					sprintf(name,"h_d_singles_energy_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1);
					h_d_singles_energy_DT800_zVSr2[iad][ir2][iz]=new TH1F(name,name,230,0.7,3.);
				}
			}
		}

			TH2F* h_acc_delayedVsDT[maxAD]; //delayed vs. DT histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_delayedVsDT_ad%d",iad+1);
				h_acc_delayedVsDT[iad]=new TH2F(name,name,500,0,10,2300,0.7,3.);
			}

			TH2F* h_acc_delayedVsDT_Ep35[maxAD]; //delayed vs. DT histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_delayedVsDT_Ep35_ad%d",iad+1);
				h_acc_delayedVsDT_Ep35[iad]=new TH2F(name,name,500,0,10,2300,0.7,3.);
			}

		TH1F* h_d_singles_energy_largeDist[maxAD]; //energy histogram for delayed with less than 1.5 meters in separation
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_largeDist_ad%d",iad+1);
		//	h_d_singles_energy_largeDist[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_largeDist[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_d_singles_energy_shortDist[maxAD]; //energy histogram for delayed with less than half a meter in separation
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_energy_shortDist_ad%d",iad+1);
		//	h_d_singles_energy_shortDist[iad]=new TH1F(name,name,150,1.5,3.);
			h_d_singles_energy_shortDist[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH2F* h_d_singles_locations[maxAD]; //where the events are happening in the AD histogram for delayed
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_d_singles_locations_ad%d",iad+1);
			h_d_singles_locations[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH1F* h_acc_distance[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_distance_ad%d", iad+1);
			h_acc_distance[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_acc_distance_3sig[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_distance_3sig_ad%d", iad+1);
			h_acc_distance_3sig[iad]=new TH1F(name,name,700,0,7.);
		}

//Delayed energy plots for various distance ranges
		TH1F* h_delayed_energy_scaled_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_delayed_energy_scaled_dist%d_ad%d",dist, iad+1);
			//	h_delayed_energy_scaled_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_delayed_energy_scaled_dist[iad][dist]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_prompt_energy_scaled_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_prompt_energy_scaled_dist%d_ad%d",dist, iad+1);
			//	h_prompt_energy_scaled_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_prompt_energy_scaled_dist[iad][dist]=new TH1F(name,name,113,0.7,12.);
			}
		}

		TH1F* h_delayed_energy_scaled_p35[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_p35_ad%d", iad+1);
		//	h_delayed_energy_scaled_p35[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_p35[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_15_22[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_15_22_ad%d", iad+1);
		//	h_delayed_energy_scaled_15_22[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_15_22[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_22_29[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_22_29_ad%d", iad+1);
		//	h_delayed_energy_scaled_22_29[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_22_29[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_29_36[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_29_36_ad%d", iad+1);
		//	h_delayed_energy_scaled_29_36[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_29_36[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_36_43[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_36_43_ad%d", iad+1);
		//	h_delayed_energy_scaled_36_43[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_36_43[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_delayed_energy_scaled_43_50[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_scaled_43_50_ad%d", iad+1);
		//	h_delayed_energy_scaled_43_50[iad]=new TH1F(name,name,150,1.5,3.);
			h_delayed_energy_scaled_43_50[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH2F* h_delayed_energy_vs_distance[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_delayed_energy_vs_distance_ad%d", iad+1);
		//	h_delayed_energy_vs_distance[iad]=new TH2F(name,name,150,1.5,3.,10,0,5.);
			h_delayed_energy_vs_distance[iad]=new TH2F(name,name,230,0.7,3.,11,0,5.5);
		}

		TH2F* h_prompt_energy_vs_distance[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_prompt_energy_vs_distance_ad%d", iad+1);
		//	h_prompt_energy_vs_distance[iad]=new TH2F(name,name,150,1.5,3.,10,0,5.);
			h_prompt_energy_vs_distance[iad]=new TH2F(name,name,113,0.7,12.,11,0,5.5);
		}

			TH2F* h_acc_distVStime[maxAD];
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_distVStime_ad%d",iad+1);
				h_acc_distVStime[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

			TH2F* h_acc_distVStime_Ep35[maxAD];
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_distVStime_Ep35_ad%d",iad+1);
				h_acc_distVStime_Ep35[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

			TH2F* h_acc_promptVStime[maxAD]; //prompt vs. time histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_promptVStime_ad%d",iad+1);
				h_acc_promptVStime[iad]=new TH2F(name,name,1999,1,2000,113,0.7,12.);
			}

			TH2F* h_acc_promptVStime_DT800[maxAD]; //prompt vs. time histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_acc_promptVStime_DT800_ad%d",iad+1);
				h_acc_promptVStime_DT800[iad]=new TH2F(name,name,1999,1,2000,113,0.7,12.);
			}

		TH1F* h_acc_distance_Ep35[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_distance_Ep35_ad%d", iad+1);
			h_acc_distance_Ep35[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1D* h_acc_DT[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_DT_ad%d", iad+1);
			h_acc_DT[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_acc_DT_3sig[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_DT_3sig_ad%d", iad+1);
			h_acc_DT_3sig[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_acc_DT_Ep35[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_DT_Ep35_ad%d", iad+1);
			h_acc_DT_Ep35[iad]=new TH1D(name,name,500,0,10);
		}
//End of energy plots for different distance requirements


	double p_DAQ_ad1 = 0;
	double p_veto_ad1 = 0;
	double p_live_ad1 = 0;
	double p_DAQ_ad2 = 0;
	double p_veto_ad2 = 0;
	double p_live_ad2 = 0;
	double p_DAQ_ad3 = 0;
	double p_veto_ad3 = 0;
	double p_live_ad3 = 0;
	double p_DAQ_ad4 = 0;
	double p_veto_ad4 = 0;
	double p_live_ad4 = 0;
	double p_counts = 0;
	double d_DAQ_ad1 = 0;
	double d_veto_ad1 = 0;
	double d_live_ad1 = 0;
	double d_DAQ_ad2 = 0;
	double d_veto_ad2 = 0;
	double d_live_ad2 = 0;
	double d_DAQ_ad3 = 0;
	double d_veto_ad3 = 0;
	double d_live_ad3 = 0;
	double d_DAQ_ad4 = 0;
	double d_veto_ad4 = 0;
	double d_live_ad4 = 0;
	double d_counts = 0;


		TTree* tr_live=new TTree("tr_accidentals_time", "Summary of Singles Search");
			tr_live->Branch("p_DAQ_ad1", &p_DAQ_ad1, "DAQ for AD1/Double");
			tr_live->Branch("p_veto_ad1", &p_veto_ad1, "Veto for AD1/Double");
			tr_live->Branch("p_live_ad1", &p_live_ad1, "Live for AD1/Double");
			tr_live->Branch("p_DAQ_ad2", &p_DAQ_ad2, "DAQ for AD2/Double");
			tr_live->Branch("p_veto_ad2", &p_veto_ad2, "Veto for AD2/Double");
			tr_live->Branch("p_live_ad2", &p_live_ad2, "Live for AD2/Double");
			tr_live->Branch("p_DAQ_ad3", &p_DAQ_ad3, "DAQ for AD3/Double");
			tr_live->Branch("p_veto_ad3", &p_veto_ad3, "Veto for AD3/Double");
			tr_live->Branch("p_live_ad3", &p_live_ad3, "Live for AD3/Double");
			tr_live->Branch("p_DAQ_ad4", &p_DAQ_ad4, "DAQ for AD4/Double");
			tr_live->Branch("p_veto_ad4", &p_veto_ad4, "Veto for AD4/Double");
			tr_live->Branch("p_live_ad4", &p_live_ad4, "Live for AD4/Double");
			tr_live->Branch("d_DAQ_ad1", &d_DAQ_ad1, "DAQ for AD1/Double");
			tr_live->Branch("d_veto_ad1", &d_veto_ad1, "Veto for AD1/Double");
			tr_live->Branch("d_live_ad1", &d_live_ad1, "Live for AD1/Double");
			tr_live->Branch("d_DAQ_ad2", &d_DAQ_ad2, "DAQ for AD2/Double");
			tr_live->Branch("d_veto_ad2", &d_veto_ad2, "Veto for AD2/Double");
			tr_live->Branch("d_live_ad2", &d_live_ad2, "Live for AD2/Double");
			tr_live->Branch("d_DAQ_ad3", &d_DAQ_ad3, "DAQ for AD3/Double");
			tr_live->Branch("d_veto_ad3", &d_veto_ad3, "Veto for AD3/Double");
			tr_live->Branch("d_live_ad3", &d_live_ad3, "Live for AD3/Double");
			tr_live->Branch("d_DAQ_ad4", &d_DAQ_ad4, "DAQ for AD4/Double");
			tr_live->Branch("d_veto_ad4", &d_veto_ad4, "Veto for AD4/Double");
			tr_live->Branch("d_live_ad4", &d_live_ad4, "Live for AD4/Double");


	float p_energy = 0;
	float d_energy = 0;
	float p_x = 0;
	float p_y = 0;
	float p_z = 0;
	float d_x = 0;
	float d_y = 0;
	float d_z = 0;
	float distance = 0;
	int det_num = 0;
	double ranTime = 0;
	double DTvalue = 0;

	int pair_entries = 0;


	for(int part_num = 1; part_num <4; part_num++){
		//Reading in the accidentals file
		char* inputname;
//		sprintf(inputname,"./accResults/EH%d/AccidentalsSummary_%d_%d_part%d.root",EH,pd_window_microsec,run_num,part_num);
//		if(dElow > 1.5) sprintf(inputname,"./accResults/EH%d/3sig/AccidentalsSummary_%d_%d_part%d.root",EH,pd_window_microsec,run_num,part_num);
//		sprintf(inputname,"./accResults/EH%d/AccidentalsSummary_TcLong_%d_part%d.root",EH,run_num,part_num);
		if(part_num == 1) inputname = outputname_summarize_long1;
		if(part_num == 2) inputname = outputname_summarize_long2;
		if(part_num == 3) inputname = outputname_summarize_long3;
		cout << "Retrieving: " << inputname << endl;
		TFile *accInfile = new TFile(inputname,"READ");

		if(part_num == 1){
			for(int iad = 0; iad < maxAD; ++iad){
				char name[64];
				sprintf(name,"h_d_singles_energy_raw_ad%d",iad+1);
				TH1F *h_d_raw = (TH1F*)accInfile->Get(name);
				h_d_singles_energy_raw[iad]->Add(h_d_raw);

				sprintf(name,"h_p_singles_energy_raw_ad%d",iad+1);
				TH1F *h_p_raw = (TH1F*)accInfile->Get(name);
				h_p_singles_energy_raw[iad]->Add(h_p_raw);
			}
		}

		//Getting the pair tree
		TTree *tr_pair = (TTree*)accInfile->Get("tr_pair");
		tr_pair->SetBranchStatus("*",1);

			//Establishing the new TTree of randomly paired singles
			tr_pair->SetBranchAddress("p_energy", &p_energy);
			tr_pair->SetBranchAddress("d_energy", &d_energy);
			tr_pair->SetBranchAddress("p_x", &p_x);
			tr_pair->SetBranchAddress("d_x", &d_x);
			tr_pair->SetBranchAddress("p_y", &p_y);
			tr_pair->SetBranchAddress("d_y", &d_y);
			tr_pair->SetBranchAddress("p_z", &p_z);
			tr_pair->SetBranchAddress("d_z", &d_z);
			tr_pair->SetBranchAddress("det_num", &det_num);
			tr_pair->SetBranchAddress("distance", &distance);


			//TTree tr_sum
			TTree *tr_time = (TTree*)accInfile->Get("tr_singles_summary");
			tr_time->SetBranchStatus("*",1);


				tr_time->SetBranchAddress("p_DAQ_ad1", &p_DAQ_ad1);
				tr_time->SetBranchAddress("p_veto_ad1", &p_veto_ad1);
				tr_time->SetBranchAddress("p_live_ad1", &p_live_ad1);
				tr_time->SetBranchAddress("p_DAQ_ad2", &p_DAQ_ad2);
				tr_time->SetBranchAddress("p_veto_ad2", &p_veto_ad2);
				tr_time->SetBranchAddress("p_live_ad2", &p_live_ad2);
				tr_time->SetBranchAddress("p_DAQ_ad3", &p_DAQ_ad3);
				tr_time->SetBranchAddress("p_veto_ad3", &p_veto_ad3);
				tr_time->SetBranchAddress("p_live_ad3", &p_live_ad3);
				tr_time->SetBranchAddress("p_DAQ_ad4", &p_DAQ_ad4);
				tr_time->SetBranchAddress("p_veto_ad4", &p_veto_ad4);
				tr_time->SetBranchAddress("p_live_ad4", &p_live_ad4);
				tr_time->SetBranchAddress("d_DAQ_ad1", &d_DAQ_ad1);
				tr_time->SetBranchAddress("d_veto_ad1", &d_veto_ad1);
				tr_time->SetBranchAddress("d_live_ad1", &d_live_ad1);
				tr_time->SetBranchAddress("d_DAQ_ad2", &d_DAQ_ad2);
				tr_time->SetBranchAddress("d_veto_ad2", &d_veto_ad2);
				tr_time->SetBranchAddress("d_live_ad2", &d_live_ad2);
				tr_time->SetBranchAddress("d_DAQ_ad3", &d_DAQ_ad3);
				tr_time->SetBranchAddress("d_veto_ad3", &d_veto_ad3);
				tr_time->SetBranchAddress("d_live_ad3", &d_live_ad3);
				tr_time->SetBranchAddress("d_DAQ_ad4", &d_DAQ_ad4);
				tr_time->SetBranchAddress("d_veto_ad4", &d_veto_ad4);
				tr_time->SetBranchAddress("d_live_ad4", &d_live_ad4);




		int nentries = tr_time->GetEntries();

		cout << "nentries = " << nentries << endl;

		if(part_num == 1){
			for(int ientry=0; ientry<nentries; ientry++){
				tr_time->GetEntry(0);
				tr_live->Fill();
			}
		}

		cout << fixed;
		cout << "Delayed live: " << setprecision(9) << d_live_ad1 << endl;
		cout << "Prompt live: " << setprecision(9) << p_live_ad1 << endl;

		pair_entries = tr_pair->GetEntries();

		//Random Funtion:
		TRandom3* random = new TRandom3();

		cout << "Num Pair_entries: " << pair_entries << endl;

		//Plotting tr_pair:
		for(int ientry=0; ientry<pair_entries; ientry++){
			tr_pair->GetEntry(ientry);
			if(p_energy < pElow || p_energy > pEhigh) continue;
			if(d_energy < dElow || d_energy > dEhigh) continue;

			ranTime = (double) random->Uniform(1,pd_window_microsec);
			DTvalue = distance*1.e-3+(ranTime)/600;

	//		if(DTvalue >= .8) continue;

			int temp_iZ = int((d_z*1.e-3+2.5)/(5./nZbins));
			int temp_iR2 = int(((pow(p_x*1.e-3,2.)+pow(p_y*1.e-3,2.)))/(5./nR2bins));

			h_d_singles_energy[det_num]->Fill(d_energy);
			h_d_singles_energy_fine[det_num]->Fill(d_energy);
			h_acc_delayedVsDT[det_num]->Fill(DTvalue, d_energy);
			if(temp_iZ < nZbins && temp_iZ >= 0) h_d_singles_energy_z[det_num][temp_iZ]->Fill(d_energy);
			if(temp_iR2 < nR2bins && temp_iR2 >= 0) h_d_singles_energy_r2[det_num][temp_iR2]->Fill(d_energy);
			if(temp_iR2 < nR2bins && temp_iR2 >= 0 && temp_iZ < nZbins && temp_iZ >= 0 && p_energy >= 3.5) h_d_singles_energy_zVSr2[det_num][temp_iR2][temp_iZ]->Fill(d_energy);
			h_d_singles_locations[det_num]->Fill((pow(d_x*1.e-3,2.)+pow(d_y*1.e-3,2.)),d_z*1.e-3);

			h_acc_distVStime[det_num]->Fill(ranTime,distance*1.e-3);
			if(p_energy >=3.5){
				h_acc_distVStime_Ep35[det_num]->Fill(ranTime,distance*1.e-3);
				h_acc_distance_Ep35[det_num]->Fill(distance*1.e-3);
				h_acc_DT_Ep35[det_num]->Fill(DTvalue);
				h_acc_delayedVsDT_Ep35[det_num]->Fill(DTvalue, d_energy);
			}

			h_acc_DT[det_num]->Fill(DTvalue);

			h_acc_promptVStime[det_num]->Fill(ranTime, p_energy);

			if(DTvalue < .8){
				h_d_singles_energy_DT800[det_num]->Fill(d_energy);
				h_d_singles_energy_fine_DT800[det_num]->Fill(d_energy);
				if(p_energy >= 3.5) h_delayed_energy_scaled_DT800_Ep35[det_num]->Fill(d_energy);
				if(temp_iZ < nZbins && temp_iZ >= 0) h_d_singles_energy_DT800_z[det_num][temp_iZ]->Fill(d_energy);
				if(temp_iR2 < nR2bins && temp_iR2 >= 0) h_d_singles_energy_DT800_r2[det_num][temp_iR2]->Fill(d_energy);
				if(temp_iR2 < nR2bins && temp_iR2 >= 0 && temp_iZ < nZbins && temp_iZ >= 0 && p_energy >= 3.5) h_d_singles_energy_DT800_zVSr2[det_num][temp_iR2][temp_iZ]->Fill(d_energy);
				h_p_singles_energy_DT800[det_num]->Fill(p_energy);
				h_acc_promptVStime_DT800[det_num]->Fill(ranTime, p_energy);
				if(d_energy >= (peak_Ed[2*(EH-1)+det_num] - 3*sigma_Ed[2*(EH-1)+det_num]) && d_energy <= (peak_Ed[2*(EH-1)+det_num] + 3*sigma_Ed[2*(EH-1)+det_num])) h_p_singles_energy_DT800_3sig[det_num]->Fill(p_energy);
			}
			
		if(d_energy >= (peak_Ed[2*(EH-1)+det_num] - 3*sigma_Ed[2*(EH-1)+det_num]) && d_energy <= (peak_Ed[2*(EH-1)+det_num] + 3*sigma_Ed[2*(EH-1)+det_num])){
			if(DTvalue < 0.3) h_p_singles_energy_DT300_3sig[det_num]->Fill(p_energy);
			if(DTvalue < 0.5) h_p_singles_energy_DT500_3sig[det_num]->Fill(p_energy);
			if(DTvalue < 1.0) h_p_singles_energy_DT1000_3sig[det_num]->Fill(p_energy);
			if(DTvalue < 1.5) h_p_singles_energy_DT1500_3sig[det_num]->Fill(p_energy);
			if(DTvalue < 2.0) h_p_singles_energy_DT2000_3sig[det_num]->Fill(p_energy);
		}

			if(DTvalue < 0.3) h_d_singles_energy_fine_DT300[det_num]->Fill(d_energy);
			if(DTvalue < 0.5) h_d_singles_energy_fine_DT500[det_num]->Fill(d_energy);
			if(DTvalue < 1.0) h_d_singles_energy_fine_DT1000[det_num]->Fill(d_energy);
			if(DTvalue < 1.5) h_d_singles_energy_fine_DT1500[det_num]->Fill(d_energy);
			if(DTvalue < 2.0) h_d_singles_energy_fine_DT2000[det_num]->Fill(d_energy);

			if(d_energy >= (peak_Ed[2*(EH-1)+det_num] - 3*sigma_Ed[2*(EH-1)+det_num]) && d_energy <= (peak_Ed[2*(EH-1)+det_num] + 3*sigma_Ed[2*(EH-1)+det_num])){
				h_acc_distance_3sig[det_num]->Fill(distance*1.e-3);
				h_acc_DT_3sig[det_num]->Fill(DTvalue);
			}

			if(distance*1.e-3 > 5.0)h_delayed_energy_vs_distance[det_num]->Fill(d_energy, 5.25);
			else h_delayed_energy_vs_distance[det_num]->Fill(d_energy, distance*1.e-3);

			if(distance*1.e-3 > 5.0)h_prompt_energy_vs_distance[det_num]->Fill(p_energy, 5.25);
			else h_prompt_energy_vs_distance[det_num]->Fill(p_energy, distance*1.e-3);

			h_acc_energy_before[det_num]->Fill(p_energy, d_energy);

			h_p_singles_energy[det_num]->Fill(p_energy);
			h_p_singles_locations[det_num]->Fill((pow(p_x*1.e-3,2.)+pow(p_y*1.e-3,2.)),p_z*1.e-3);
			h_acc_distance[det_num]->Fill(distance*1.e-3);

			if(p_energy >= 3.5){
				h_delayed_energy_scaled_p35[det_num]->Fill(d_energy);
				h_delayed_energy_scaled_Ep35[det_num]->Fill(d_energy);
			}

			if(distance*1.e-3 > 1.5 && distance*1.e-3 <= 2.2) h_delayed_energy_scaled_15_22[det_num]->Fill(d_energy);
			if(distance*1.e-3 > 2.2 && distance*1.e-3 <= 2.9) h_delayed_energy_scaled_22_29[det_num]->Fill(d_energy);
			if(distance*1.e-3 > 2.9 && distance*1.e-3 <= 3.6) h_delayed_energy_scaled_29_36[det_num]->Fill(d_energy);
			if(distance*1.e-3 > 3.6 && distance*1.e-3 <= 4.3) h_delayed_energy_scaled_36_43[det_num]->Fill(d_energy);
			if(distance*1.e-3 > 4.3) h_delayed_energy_scaled_43_50[det_num]->Fill(d_energy);

			if(distance*1.e-3 <= 1.0){
				h_delayed_energy_scaled_dist[det_num][0]->Fill(d_energy);
				h_prompt_energy_scaled_dist[det_num][0]->Fill(p_energy);
				h_acc_energy_1m[det_num]->Fill(p_energy, d_energy);
			}
			if(distance*1.e-3 > 1.0 && distance*1.e-3 <= 2.0){
				h_delayed_energy_scaled_dist[det_num][1]->Fill(d_energy);
				h_prompt_energy_scaled_dist[det_num][1]->Fill(p_energy);
			}
			if(distance*1.e-3 > 2.0 && distance*1.e-3 <= 3.0){
				h_delayed_energy_scaled_dist[det_num][2]->Fill(d_energy);
				h_prompt_energy_scaled_dist[det_num][2]->Fill(p_energy);
			}
			if(distance*1.e-3 > 3.0 && distance*1.e-3 <= 4.0){
				h_delayed_energy_scaled_dist[det_num][3]->Fill(d_energy);
				h_prompt_energy_scaled_dist[det_num][3]->Fill(p_energy);
			}
			if(distance*1.e-3 > 4.0 && distance*1.e-3 <= 5.0){
				h_delayed_energy_scaled_dist[det_num][4]->Fill(d_energy);
				h_prompt_energy_scaled_dist[det_num][4]->Fill(p_energy);
			}

			if(distance*1.e-3 > 5.){
				h_delayed_energy_scaled_dist[det_num][5]->Fill(d_energy);
				h_prompt_energy_scaled_dist[det_num][5]->Fill(p_energy);
			}

			if(distance*1.e-3 > 1.5) continue;

			h_p_singles_energy_largeDist[det_num]->Fill(p_energy);
			h_d_singles_energy_largeDist[det_num]->Fill(d_energy);

			if(distance*1.e-3 > 0.5) continue;

			h_p_singles_energy_shortDist[det_num]->Fill(p_energy);
			h_d_singles_energy_shortDist[det_num]->Fill(d_energy);

		}
	}


	//SAVING TO OUTPUT FILE
//        char outputname[64];
//        sprintf(outputname,"./accResults/EH%d/AccidentalsPlots_%d_%d.root",EH,pd_window_microsec,run_num);
//		if(dElow > 1.5) sprintf(outputname,"./accResults/EH%d/3sig/AccidentalsPlots_%d_%d.root",EH,pd_window_microsec,run_num);
//        sprintf(outputname,"./accResults/EH%d/AccidentalsPlots_TcLong_%d.root",EH,run_num);

	cout << "Saving to: " << outputname_plot << endl;
//        sprintf(outputname,"./AccidentalsSummary_%d.root",run_num);
	TFile* outfile=new TFile(outputname_plot, "RECREATE");
		outfile->cd();
		tr_live->Write();
		for(int iad=0; iad<maxAD; ++iad){
			h_acc_energy_before[iad]->SetOption("COLZ");
			h_acc_energy_before[iad]->SetStats(0);
			h_acc_energy_before[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_acc_energy_before[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_acc_energy_before[iad]->Write();

			h_acc_energy_1m[iad]->SetOption("COLZ");
			h_acc_energy_1m[iad]->SetStats(0);
			h_acc_energy_1m[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_acc_energy_1m[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_acc_energy_1m[iad]->Write();

			h_acc_distVStime[iad]->SetStats(0);
			h_acc_distVStime[iad]->SetOption("COLZ");
			h_acc_distVStime[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_acc_distVStime[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_acc_distVStime[iad]->Write();

			h_acc_distVStime_Ep35[iad]->SetStats(0);
			h_acc_distVStime_Ep35[iad]->SetOption("COLZ");
			h_acc_distVStime_Ep35[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_acc_distVStime_Ep35[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_acc_distVStime_Ep35[iad]->Write();

		//		h_acc_DT[iad]->SetStats(0);
				h_acc_DT[iad]->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT[iad]->GetYaxis()->SetTitle("Counts");
				h_acc_DT[iad]->Write();

		//		h_acc_DT_3sig[iad]->SetStats(0);
				h_acc_DT_3sig[iad]->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_acc_DT_3sig[iad]->Write();

		//		h_acc_DT_Ep35[iad]->SetStats(0);
				h_acc_DT_Ep35[iad]->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_Ep35[iad]->GetYaxis()->SetTitle("Counts");
				h_acc_DT_Ep35[iad]->Write();

			h_p_singles_locations[iad]->SetOption("COLZ");
			h_p_singles_locations[iad]->SetStats(0);
			h_p_singles_locations[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_p_singles_locations[iad]->GetYaxis()->SetTitle("z [m]");
			h_p_singles_locations[iad]->Write();

		//	h_p_singles_energy_raw[iad]->SetStats(0);
			h_p_singles_energy_raw[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_raw[iad]->Write();

		//	h_p_singles_energy[iad]->SetStats(0);
			h_p_singles_energy[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy[iad]->Write();

		//	h_p_singles_energy_DT800[iad]->SetStats(0);
			h_p_singles_energy_DT800[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_DT800[iad]->Write();

		//	h_p_singles_energy_DT800_3sig[iad]->SetStats(0);
			h_p_singles_energy_DT800_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_DT800_3sig[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_DT800_3sig[iad]->Write();

			//	h_p_singles_energy_DT300_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT300_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT300_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT300_3sig[iad]->Write();

			//	h_p_singles_energy_DT500_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT500_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT500_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT500_3sig[iad]->Write();

			//	h_p_singles_energy_DT1000_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT1000_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT1000_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT1000_3sig[iad]->Write();

			//	h_p_singles_energy_DT1500_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT1500_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT1500_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT1500_3sig[iad]->Write();

			//	h_p_singles_energy_DT2000_3sig[iad]->SetStats(0);
				h_p_singles_energy_DT2000_3sig[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_p_singles_energy_DT2000_3sig[iad]->GetYaxis()->SetTitle("Counts");
				h_p_singles_energy_DT2000_3sig[iad]->Write();

		//	h_p_singles_energy_largeDist[iad]->SetStats(0);
			h_p_singles_energy_largeDist[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_largeDist[iad]->Write();

		//	h_p_singles_energy_shortDist[iad]->SetStats(0);
			h_p_singles_energy_shortDist[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_p_singles_energy_shortDist[iad]->GetYaxis()->SetTitle("Counts");
			h_p_singles_energy_shortDist[iad]->Write();

			h_d_singles_locations[iad]->SetOption("COLZ");
			h_d_singles_locations[iad]->SetStats(0);
			h_d_singles_locations[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_d_singles_locations[iad]->GetYaxis()->SetTitle("z [m]");
			h_d_singles_locations[iad]->Write();

		//	h_d_singles_energy_raw[iad]->SetStats(0);
			h_d_singles_energy_raw[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_raw[iad]->Write();

		//	h_d_singles_energy[iad]->SetStats(0);
			h_d_singles_energy[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy[iad]->Write();

		//	h_d_singles_energy_DT800[iad]->SetStats(0);
			h_d_singles_energy_DT800[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_DT800[iad]->Write();

		//	h_d_singles_energy_fine[iad]->SetStats(0);
			h_d_singles_energy_fine[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_fine[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_fine[iad]->Write();

		//	h_d_singles_energy_fine_DT800[iad]->SetStats(0);
			h_d_singles_energy_fine_DT800[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_fine_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_fine_DT800[iad]->Write();

			//	h_d_singles_energy_fine_DT300[iad]->SetStats(0);
				h_d_singles_energy_fine_DT300[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT300[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT300[iad]->Write();

			//	h_d_singles_energy_fine_DT500[iad]->SetStats(0);
				h_d_singles_energy_fine_DT500[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT500[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT500[iad]->Write();

			//	h_d_singles_energy_fine_DT1000[iad]->SetStats(0);
				h_d_singles_energy_fine_DT1000[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT1000[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT1000[iad]->Write();

			//	h_d_singles_energy_fine_DT1500[iad]->SetStats(0);
				h_d_singles_energy_fine_DT1500[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT1500[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT1500[iad]->Write();

			//	h_d_singles_energy_fine_DT2000[iad]->SetStats(0);
				h_d_singles_energy_fine_DT2000[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_fine_DT2000[iad]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_fine_DT2000[iad]->Write();

		//	h_delayed_energy_scaled_Ep35[iad]->SetStats(0);
			h_delayed_energy_scaled_Ep35[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_delayed_energy_scaled_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_Ep35[iad]->Write();

		//	h_delayed_energy_scaled_DT800_Ep35[iad]->SetStats(0);
			h_delayed_energy_scaled_DT800_Ep35[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_delayed_energy_scaled_DT800_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_DT800_Ep35[iad]->Write();

		/*//	h_delayed_energy_fine_DT800_Ep35[iad]->SetStats(0);
			h_delayed_energy_fine_DT800_Ep35[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_delayed_energy_fine_DT800_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_fine_DT800_Ep35[iad]->Write();*/

				h_acc_delayedVsDT[iad]->SetStats(0);
				h_acc_delayedVsDT[iad]->SetOption("COLZ");
				h_acc_delayedVsDT[iad]->GetXaxis()->SetTitle("DT Value [m]");
				h_acc_delayedVsDT[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
				h_acc_delayedVsDT[iad]->Write();

				h_acc_delayedVsDT_Ep35[iad]->SetStats(0);
				h_acc_delayedVsDT_Ep35[iad]->SetOption("COLZ");
				h_acc_delayedVsDT_Ep35[iad]->GetXaxis()->SetTitle("DT Value [m]");
				h_acc_delayedVsDT_Ep35[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
				h_acc_delayedVsDT_Ep35[iad]->Write();

			for(int iz = 0; iz<NzBins; iz++){
				h_d_singles_energy_z[iad][iz]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_z[iad][iz]->Write();

				h_d_singles_energy_DT800_z[iad][iz]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_DT800_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_DT800_z[iad][iz]->Write();
			}

			for(int ir2 = 0; ir2<Nr2Bins; ir2++){
				h_d_singles_energy_r2[iad][ir2]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_r2[iad][ir2]->Write();

				h_d_singles_energy_DT800_r2[iad][ir2]->GetXaxis()->SetTitle("Singles Energy [MeV]");
				h_d_singles_energy_DT800_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_d_singles_energy_DT800_r2[iad][ir2]->Write();
			}

			for(int iz = 0; iz<NzBins; iz++){
				for(int ir2 = 0; ir2<Nr2Bins; ir2++){
					h_d_singles_energy_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Singles Energy [MeV]");
					h_d_singles_energy_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_d_singles_energy_zVSr2[iad][ir2][iz]->Write();

					h_d_singles_energy_DT800_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Singles Energy [MeV]");
					h_d_singles_energy_DT800_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_d_singles_energy_DT800_zVSr2[iad][ir2][iz]->Write();
				}
			}

		//	h_d_singles_energy_largeDist[iad]->SetStats(0);
			h_d_singles_energy_largeDist[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_largeDist[iad]->Write();

		//	h_d_singles_energy_shortDist[iad]->SetStats(0);
			h_d_singles_energy_shortDist[iad]->GetXaxis()->SetTitle("Singles Energy [MeV]");
			h_d_singles_energy_shortDist[iad]->GetYaxis()->SetTitle("Counts");
			h_d_singles_energy_shortDist[iad]->Write();

			h_acc_distance[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_acc_distance[iad]->GetYaxis()->SetTitle("Counts");
			h_acc_distance[iad]->Write();

			h_acc_distance_3sig[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_acc_distance_3sig[iad]->GetYaxis()->SetTitle("Counts");
			h_acc_distance_3sig[iad]->Write();

			h_acc_distance_Ep35[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_acc_distance_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_acc_distance_Ep35[iad]->Write();

			h_delayed_energy_scaled_15_22[iad]->SetStats(0);
			h_delayed_energy_scaled_15_22[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_15_22[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_15_22[iad]->Write();

			h_delayed_energy_scaled_22_29[iad]->SetStats(0);
			h_delayed_energy_scaled_22_29[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_22_29[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_22_29[iad]->Write();

			h_delayed_energy_scaled_29_36[iad]->SetStats(0);
			h_delayed_energy_scaled_29_36[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_29_36[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_29_36[iad]->Write();

			h_delayed_energy_scaled_36_43[iad]->SetStats(0);
			h_delayed_energy_scaled_36_43[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_36_43[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_36_43[iad]->Write();

			h_delayed_energy_scaled_43_50[iad]->SetStats(0);
			h_delayed_energy_scaled_43_50[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_43_50[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_43_50[iad]->Write();

			h_delayed_energy_scaled_p35[iad]->SetStats(0);
			h_delayed_energy_scaled_p35[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_scaled_p35[iad]->GetYaxis()->SetTitle("Counts");
			h_delayed_energy_scaled_p35[iad]->Write();

			h_delayed_energy_vs_distance[iad]->SetStats(0);
			h_delayed_energy_vs_distance[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_vs_distance[iad]->GetYaxis()->SetTitle("Distance between prompt and delayed");
			h_delayed_energy_vs_distance[iad]->SetOption("COLZ");
			h_delayed_energy_vs_distance[iad]->Write();

			h_prompt_energy_vs_distance[iad]->SetStats(0);
			h_prompt_energy_vs_distance[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_prompt_energy_vs_distance[iad]->GetYaxis()->SetTitle("Distance between prompt and delayed");
			h_prompt_energy_vs_distance[iad]->SetOption("COLZ");
			h_prompt_energy_vs_distance[iad]->Write();

			h_acc_promptVStime[iad]->SetStats(0);
			h_acc_promptVStime[iad]->SetOption("COLZ");
			h_acc_promptVStime[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_acc_promptVStime[iad]->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_acc_promptVStime[iad]->Write();

			h_acc_promptVStime_DT800[iad]->SetStats(0);
			h_acc_promptVStime_DT800[iad]->SetOption("COLZ");
			h_acc_promptVStime_DT800[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_acc_promptVStime_DT800[iad]->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_acc_promptVStime_DT800[iad]->Write();

			for(int dist=0; dist<6; dist++){
				h_delayed_energy_scaled_dist[iad][dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_delayed_energy_scaled_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_delayed_energy_scaled_dist[iad][dist]->Write();
				h_prompt_energy_scaled_dist[iad][dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_prompt_energy_scaled_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_prompt_energy_scaled_dist[iad][dist]->Write();
			}
		}


	outfile->Close();

}








void all(int run_order, int pd_window_microsec){
	int run_num=-1;
	int EH=-1;

	FILE* runfile=fopen("./run_list_good_sync.txt","r");  //This file contains all good runs numbers and experimental hall they belong to
	int iline=0;
	while(1){ //go though the file and retrieve run_num (run number) which corresponds to run_order provided as a parameter of GetMuons function
		fscanf(runfile,"%d %d",&run_num,&EH);
		if(feof(runfile)) break;
		if(iline==run_order) break;
		iline++;
	}
	fclose(runfile);

        sprintf(outputname_findPrompt,"./promptSingles/EH%d/foundPromptSingles_NU_%d_%d.root",EH,pd_window_microsec,run_num);
        sprintf(outputname_findDelayed,"./delayedSingles/EH%d/foundDelayedSingles_NU_%d_%d.root",EH,pd_window_microsec,run_num);
        sprintf(outputname_summarize,"./accResults/EH%d/AccidentalsSummary_NU_%d_%d.root",EH,pd_window_microsec,run_num);
        sprintf(outputname_summarize_long1,"./accResults/EH%d/AccidentalsSummary_NU_%d_%d_part1.root",EH,pd_window_microsec,run_num);
        sprintf(outputname_summarize_long2,"./accResults/EH%d/AccidentalsSummary_NU_%d_%d_part2.root",EH,pd_window_microsec,run_num);
        sprintf(outputname_summarize_long3,"./accResults/EH%d/AccidentalsSummary_NU_%d_%d_part3.root",EH,pd_window_microsec,run_num);
        sprintf(outputname_plot,"./accResults/EH%d/AccidentalsPlots_NU_%d_%d.root",EH,pd_window_microsec,run_num);

//	find(run_order, 1.5, 3, pd_window_microsec+400, 400, pd_window_microsec); //For delayed singles -- usually 1.75 to 2.89 MeV
//	find(run_order, 1.5, 12., pd_window_microsec+300, 500, pd_window_microsec); //For prompt singles
//	find(run_order, 1.5, 12., 500, pd_window_microsec+300, pd_window_microsec); //For "inverted" prompt singles

	if(run_num != 24614 && run_num != 37322 && run_num != 37645 && run_num != 63825 && run_num != 63941 && run_order != 672 && run_order != 2636 && run_order != 2641 && run_order != 3004 && run_order != 3123 && run_order != 625 && run_order != 2786){
//		summarize(run_num, 1.5, 12., 1.5, 3,pd_window_microsec); //1.5-3 or 1.9-2.74
		plot(run_num, 1.5, 12., 1.5, 3, pd_window_microsec);
	}

	else{
//		summarize_long(run_num, 1, 1.5, 12., 1.5, 3,pd_window_microsec); //Part 1
//		summarize_long(run_num, 2, 1.5, 12., 1.5, 3,pd_window_microsec); //Part 2
//		summarize_long(run_num, 3, 1.5, 12., 1.5, 3,pd_window_microsec); //Part 3
		plot_long(run_num, 1.5, 12., 1.5, 3, pd_window_microsec);
	}

}
