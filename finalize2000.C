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

void finalize(int hall_num, int pd_window_microsec){

	char hallList[64];
	sprintf(hallList, "EH%druns_sync.txt",hall_num);
	FILE* runfile=fopen(hallList,"r");

/*	char accList[64];
	sprintf(accList, "./IBDs/2mCounts_EH%d.txt",hall_num);
	ofstream accCountsFile;
	accCountsFile.open(accList);*/


	int maxAD = 0;
	int nRuns = 0;
	if(hall_num == 1){
		maxAD = 2;
//		nRuns = 1134;
		nRuns = 1035;
	}
	if(hall_num == 2){
		maxAD = 2;
//		nRuns = 1125;
		nRuns = 1034;
	}
	if(hall_num == 3){
		maxAD = 4;
//		nRuns = 1215;
		nRuns = 1099;
	}

	const int NzBins = 20;
	const int Nr2Bins = 20;

	//Making Histograms:
		//BEFORE DISTANCE CUT HISTS
		TH2F* h_total_ibd_energy_before[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_ibd_energy_before_ad%d",iad+1);
		//	h_total_ibd_energy_before[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_total_ibd_energy_before[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH2F* h_total_ibd_energy_1m[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_ibd_energy_1m_ad%d",iad+1);
		//	h_total_ibd_energy_1m[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_total_ibd_energy_1m[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH2F* h_total_locations_before[maxAD]; //where the events are happening in the AD histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_locations_before_ad%d",iad+1);
			h_total_locations_before[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH2F* h_total_plocations_before[maxAD]; //where the prompt events are happening in the AD histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_plocations_before_ad%d",iad+1);
			h_total_plocations_before[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH2F* h_total_dlocations_before[maxAD]; //where the delay events are happening in the AD histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_dlocations_before_ad%d",iad+1);
			h_total_dlocations_before[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH1F* h_total_delayed_energy_before[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_before_ad%d", iad+1);
			h_total_delayed_energy_before[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_before_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_before_z[iad][iz]=new TH1F(Form("h_total_delayed_energy_before_z_ad%d_iz%d", iad+1, iz+1),Form("h_total_delayed_energy_before_z_ad%d_iz%d", iad+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_DT800_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_DT800_z[iad][iz]=new TH1F(Form("h_total_delayed_energy_DT800_z_ad%d_iz%d", iad+1, iz+1),Form("h_total_delayed_energy_DT800_z_ad%d_iz%d", iad+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_before_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_before_r2[iad][ir2]=new TH1F(Form("h_total_delayed_energy_before_r2_ad%d_ir2%d", iad+1, ir2+1),Form("h_total_delayed_energy_before_r2_ad%d_ir2%d", iad+1, ir2+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_DT800_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_DT800_r2[iad][ir2]=new TH1F(Form("h_total_delayed_energy_DT800_r2_ad%d_ir2%d", iad+1, ir2+1),Form("h_total_delayed_energy_DT800_r2_ad%d_ir2%d", iad+1, ir2+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_before_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_before_zVSr2[iad][ir2][iz]=new TH1F(Form("h_total_delayed_energy_before_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),Form("h_total_delayed_energy_before_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),230,0.7,3.);
				}
			}
		}

		TH1F* h_total_delayed_energy_DT800_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_DT800_zVSr2[iad][ir2][iz]=new TH1F(Form("h_total_delayed_energy_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),Form("h_total_delayed_energy_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),230,0.7,3.);
				}
			}
		}

		TH1F* h_total_delayed_energy_DT800[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_DT800_ad%d", iad+1);
			h_total_delayed_energy_DT800[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_before[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_before_ad%d", iad+1);
			h_total_delayed_energy_fine_before[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DT800[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_ad%d", iad+1);
			h_total_delayed_energy_fine_DT800[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_Ep35[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_Ep35_ad%d", iad+1);
			h_total_delayed_energy_fine_Ep35[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DT800_Ep35[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_Ep35_ad%d", iad+1);
			h_total_delayed_energy_fine_DT800_Ep35[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_prompt_energy_before[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_before_ad%d", iad+1);
			h_total_prompt_energy_before[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_DT800[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_DT800_ad%d", iad+1);
			h_total_prompt_energy_DT800[iad]=new TH1F(name,name,113,0.7,12.);
		}

			TH1F* h_total_prompt_energy_DT800_3sig[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_prompt_energy_DT800_3sig_ad%d", iad+1);
				h_total_prompt_energy_DT800_3sig[iad]=new TH1F(name,name,113,0.7,12.);
			}

		TH1D* h_total_delta_time_before[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_before_ad%d", iad+1);
			h_total_delta_time_before[iad]=new TH1D(name,name,2000,0,2000);
		}

		TH1D* h_total_delta_time_finer[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_finer_ad%d", iad+1);
			h_total_delta_time_finer[iad]=new TH1D(name,name,20000,0,2000);
		}

		TH1D* h_total_delta_time_1us[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_ad%d", iad+1);
			h_total_delta_time_1us[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_2m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_2m_ad%d", iad+1);
			h_total_delta_time_1us_2m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_15m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_15m_ad%d", iad+1);
			h_total_delta_time_1us_15m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_1m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_1m_ad%d", iad+1);
			h_total_delta_time_1us_1m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_075m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_075m_ad%d", iad+1);
			h_total_delta_time_1us_075m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_05m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_05m_ad%d", iad+1);
			h_total_delta_time_1us_05m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_gdls[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_gdls_ad%d", iad+1);
			h_total_delta_time_1us_gdls[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_gdls_1m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_gdls_1m_ad%d", iad+1);
			h_total_delta_time_1us_gdls_1m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_gdls_075m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_gdls_075m_ad%d", iad+1);
			h_total_delta_time_1us_gdls_075m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_gdls_05m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_gdls_05m_ad%d", iad+1);
			h_total_delta_time_1us_gdls_05m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_ls[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_ls_ad%d", iad+1);
			h_total_delta_time_1us_ls[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_ls_1m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_ls_1m_ad%d", iad+1);
			h_total_delta_time_1us_ls_1m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_ls_075m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_ls_075m_ad%d", iad+1);
			h_total_delta_time_1us_ls_075m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_ls_05m[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_ls_05m_ad%d", iad+1);
			h_total_delta_time_1us_ls_05m[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1D* h_total_delta_time_1us_largeDist[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_1us_largeDist_ad%d", iad+1);
			h_total_delta_time_1us_largeDist[iad]=new TH1D(name,name,1999,1,2000);
		}

		TH1F* h_total_distance_before[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_distance_before_ad%d", iad+1);
			h_total_distance_before[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_distance_before_Ep35[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_distance_before_Ep35_ad%d", iad+1);
			h_total_distance_before_Ep35[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_distance_before_400[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_distance_before_400_ad%d", iad+1);
			h_total_distance_before_400[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_distance_before_600[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_distance_before_600_ad%d", iad+1);
			h_total_distance_before_600[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_distance_before_800[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_distance_before_800_ad%d", iad+1);
			h_total_distance_before_800[iad]=new TH1F(name,name,700,0,7.);
		}

			TH2F* h_total_ibd_distVStime[maxAD]; //prompt vs. delayed energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_ibd_distVStime_ad%d",iad+1);
				h_total_ibd_distVStime[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

			TH2F* h_total_ibd_distVStime_Ep35[maxAD]; //prompt vs. delayed energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_ibd_distVStime_Ep35_ad%d",iad+1);
				h_total_ibd_distVStime_Ep35[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

			TH2F* h_total_ibd_promptVStime[maxAD]; //prompt vs. time histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_ibd_promptVStime_ad%d",iad+1);
				h_total_ibd_promptVStime[iad]=new TH2F(name,name,1999,1,2000,113,0.7,12.);
			}

			TH2F* h_total_ibd_promptVStime_DT800[maxAD]; //prompt vs. time histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_ibd_promptVStime_DT800_ad%d",iad+1);
				h_total_ibd_promptVStime_DT800[iad]=new TH2F(name,name,1999,1,2000,113,0.7,12.);
			}

		TH1D* h_total_ibd_DT[maxAD]; //DT histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_ibd_DT_ad%d", iad+1);
			h_total_ibd_DT[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_ibd_DT_Ep35[maxAD]; //DT histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_ibd_DT_Ep35_ad%d", iad+1);
			h_total_ibd_DT_Ep35[iad]=new TH1D(name,name,500,0,10);
		}

//Start of section for looking at distance ranges of delayed energy scaled plots
		TH1F* h_total_delayed_energy_scaled_15_22[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_ibd_delayed_energy_scaled_15_22_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_15_22[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_15_22[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_22_29[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_ibd_delayed_energy_scaled_22_29_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_22_29[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_22_29[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_29_36[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_ibd_delayed_energy_scaled_29_36_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_29_36[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_29_36[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_36_43[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_ibd_delayed_energy_scaled_36_43_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_36_43[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_36_43[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_43_50[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_ibd_delayed_energy_scaled_43_50_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_43_50[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_43_50[iad]=new TH1F(name,name,230,0.7,3.);
		}
//End of section for looking at distance ranges of delayed energy scaled plots

		//AFTER DISTANCE CUT HISTOGRAMS
		TH2F* h_total_ibd_energy_after[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_ibd_energy_after_ad%d",iad+1);
		//	h_total_ibd_energy_after[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_total_ibd_energy_after[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH2F* h_total_locations_after[maxAD]; //where the events are happening in the AD histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_locations_after_ad%d",iad+1);
			h_total_locations_after[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH2F* h_total_plocations_after[maxAD]; //where the prompt events are happening in the AD histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_plocations_after_ad%d",iad+1);
			h_total_plocations_after[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH2F* h_total_dlocations_after[maxAD]; //where the delay events are happening in the AD histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_dlocations_after_ad%d",iad+1);
			h_total_dlocations_after[iad]=new TH2F(name,name,100,0.,5.,100,-3.,3.);
		}

		TH1F* h_total_delayed_energy_largeDist[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_largeDist_ad%d", iad+1);
		//	h_total_delayed_energy_largeDist[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_largeDist[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_after[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_after_ad%d", iad+1);
		//	h_total_delayed_energy_after[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_after[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_prompt_energy_largeDist[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_largeDist_ad%d", iad+1);
			h_total_prompt_energy_largeDist[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_after[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_after_ad%d", iad+1);
			h_total_prompt_energy_after[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1D* h_total_delta_time_after[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delta_time_after_ad%d", iad+1);
			h_total_delta_time_after[iad]=new TH1D(name,name,2000,0,2000);
		}

		TH1F* h_total_distance_after[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_distance_after_ad%d", iad+1);
			h_total_distance_after[iad]=new TH1F(name,name,700,0,7.);
		}

		TH2F* h_total_delayed_energy_vs_distance[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_ibd_delayed_energy_vs_distance_ad%d", iad+1);
		//	h_total_delayed_energy_vs_distance[iad]=new TH2F(name,name,150,1.5,3.,10,0,5.);
			h_total_delayed_energy_vs_distance[iad]=new TH2F(name,name,230,0.7,3.,11,0,5.5);
		}

		TH2F* h_total_prompt_energy_vs_distance[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_vs_distance_ad%d", iad+1);
			h_total_prompt_energy_vs_distance[iad]=new TH2F(name,name,113,0.7,12.,11,0,5.5);
		}

/*		TH1F* h_total_GdLSdistance_before[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_GdLSdistance_before_ad%d", iad+1);
			h_total_GdLSdistance_before[iad]=new TH1F(name,name,500,0,5.);
		}

		TH2F* h_total_GdLS_ibd_energy_before[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_GdLS_ibd_energy_before_ad%d",iad+1);
			h_total_GdLS_ibd_energy_before[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
		}

		TH1F* h_total_GdLSdistance_after[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_GdLSdistance_after_ad%d", iad+1);
			h_total_GdLSdistance_after[iad]=new TH1F(name,name,500,0,5.);
		}

		TH2F* h_total_GdLS_ibd_energy_after[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_GdLS_ibd_energy_after_ad%d",iad+1);
			h_total_GdLS_ibd_energy_after[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
		}*/


		//Efficiency Plot
		TH1F* h_muon_efficiency[maxAD]; //efficiency histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_muon_efficiency_ad%d", iad+1);
			h_muon_efficiency[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_IBDperTime[maxAD]; //efficiency histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_IBDperTime_ad%d", iad+1);
			h_IBDperTime[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_total_delayed_energy_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_dist%d_ad%d",dist, iad+1);
			//	h_total_delayed_energy_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_total_delayed_energy_dist[iad][dist]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_total_prompt_energy_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_prompt_energy_dist%d_ad%d",dist, iad+1);
			//	h_total_prompt_energy_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_total_prompt_energy_dist[iad][dist]=new TH1F(name,name,113,0.7,12.);
			}
		}

		TH1F* h_total_delayed_energy_p35[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_p35_ad%d", iad+1);
		//	h_total_delayed_energy_p35[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_p35[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_dist_2plus[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_dist_2plus_ad%d", iad+1);
		//	h_total_delayed_energy_dist_2plus[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_dist_2plus[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_prompt_energy_dist_2plus[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_dist_2plus_ad%d", iad+1);
		//	h_total_prompt_energy_dist_2plus[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_prompt_energy_dist_2plus[iad]=new TH1F(name,name,113,0.7,12.);
		}

	int run_num = 0;
/*	double accCounts_ad1 = 0;
	double accCounts_ad2 = 0;
	double accCounts_ad3 = 0;
	double accCounts_ad4 = 0;

	TTree* tr_2m=new TTree("tr_2m", "Tree for counting events > 2m");
	tr_2m->Branch("run_num", &run_num, "Run Number/Int");
	tr_2m->Branch("accCounts_ad1", &accCounts_ad1, "Number of accidentals > 2m in AD1 /Double");
	tr_2m->Branch("accCounts_ad2", &accCounts_ad2, "Number of accidentals > 2m in AD2 /Double");
	tr_2m->Branch("accCounts_ad3", &accCounts_ad3, "Number of accidentals > 2m in AD3 /Double");
	tr_2m->Branch("accCounts_ad4", &accCounts_ad4, "Number of accidentals > 2m in AD4 /Double");*/

	int run_order = 0;
	int EH = 0;
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
	int ibdEntries[maxAD];
//	double accidentals[4];

	for(int iad=0; iad<maxAD; iad++){
		efficiency[iad] = 0;
		ibdEntries[iad] = 0;
	}

	ofstream resub;
	char resubName[64];
	sprintf(resubName, "./resub_EH%d.sh",hall_num);
	resub.open(resubName);
	resub << "#!/bin/bash" << endl;


	while(1){
		fscanf(runfile,"%d %d",&run_num,&EH);
		if(feof(runfile)) break; //If it's the end of the file, break.
		if(EH != hall_num){
			cout << " WARNING: HALL NUMBERS DO NOT MATCH!!! Run #" << run_num << endl;
			continue;
		}

		//For failed runs:
/*		if(run_num == 24614 || run_num == 37322 || run_num == 37645 || run_num == 63825 || run_num == 63941){
			run_order += 1;
			continue;
		}
		//For failed runs:
		if(run_num == 63825){
			run_order += 1;
			continue;
		}*/

		//For failed runs:
/*		if(run_num == 22914){
			run_order += 1;
			continue;
		}*/

	/*	if(run_order == 192 || run_order == 205 ||  run_order == 852 || run_order == 853){
			run_order+=1;
			continue;
		} 


		if(run_order >= 207){
			run_order+=1;
			continue;
		}

		if(run_order < 207 || run_order >= 414 ){
			run_order+=1;
			continue;
		}

		if(run_order < 414 || run_order >= 621 ){
			run_order+=1;
			continue;
		}

		if(run_order < 621 || run_order >= 828 ){
			run_order+=1;
			continue;
		}

		if(run_order < 828){
			run_order+=1;
			continue;
		}

		if(run_order == 94 || run_order == 192 || run_order == 205 ||  run_order == 852 || run_order == 853){
			run_order+=1;
			continue;
		}*/

/*		if(run_order == 94 || run_order == 192 || run_order == 205 ||  run_order == 852 || run_order == 853|| run_order == 70 || run_order == 71 || run_num == 86){
			run_order+=1;
			continue;
		}
*/

/*		if(run_num == 24614 || run_num == 37322 || run_num == 37645 || run_num == 63825 || run_num == 63941){
			run_order += 1;
			continue;
		}
*/
		/*if(run_order > 100){
			run_order+=1;
			continue;
		}*/

/*		if(run_num == 22008 || run_num == 22278 || run_num == 25211 || run_num == 25971 || run_num == 34548){
			run_order += 1;
			continue;
		}*/



		char runFileName[64];
//		sprintf(runFileName, "./IBDs/EH%d/summary_TcLong_2_%d.root",EH,run_num);
		sprintf(runFileName, "./IBDs/EH%d/summary_NU_%d_%d.root",EH,pd_window_microsec,run_num);
//		sprintf(runFileName, "./IBDs/EH%d/round1/summary_%d.root",EH,run_num);

		if(FILE *file = fopen(runFileName, "r")) fclose(file);
		else{
			cout << "Skipping run: " << run_num << endl;

				FILE* allRunsFile=fopen("./run_list_good_sync.txt","r");  //This file contains all good runs numbers and experimental hall they belong to
				int iline=0;
				int temp_run_num = 0;
				int temp_EH = 0;
				while(1){ //go though the file and retrieve run_order given the run number
					fscanf(allRunsFile,"%d %d",&temp_run_num,&temp_EH);
					if(feof(allRunsFile)) break;
  					if(run_num == temp_run_num) break;
					iline++;
				}
				fclose(allRunsFile);

			resub << "hep_sub -os SL6 run_ibds2000.sh -argu " << iline << " -mem 2000" << endl << "sleep 1" << endl;
			continue;
		}

//		run_order += 1;
//		continue;
		

		TFile *runFile = new TFile(runFileName);
		
		//TTree tr_sum;
		TTree *tr_efficiency = (TTree*)runFile->Get("tr_summary");
		tr_efficiency->SetBranchStatus("*",1);
		tr_efficiency->SetBranchAddress("total_DAQ_ad1",&total_DAQ_ad1);
		tr_efficiency->SetBranchAddress("total_DAQ_ad2",&total_DAQ_ad2);
		tr_efficiency->SetBranchAddress("tot_live_ad1",&tot_live_ad1);
		tr_efficiency->SetBranchAddress("tot_veto_ad1",&tot_veto_ad1);
		tr_efficiency->SetBranchAddress("tot_live_ad2",&tot_live_ad2);
		tr_efficiency->SetBranchAddress("tot_veto_ad2",&tot_veto_ad2);
		tr_efficiency->SetBranchAddress("total_DAQ_ad3",&total_DAQ_ad3);
		tr_efficiency->SetBranchAddress("total_DAQ_ad4",&total_DAQ_ad4);
		tr_efficiency->SetBranchAddress("tot_live_ad3",&tot_live_ad3);
		tr_efficiency->SetBranchAddress("tot_veto_ad3",&tot_veto_ad3);
		tr_efficiency->SetBranchAddress("tot_live_ad4",&tot_live_ad4);
		tr_efficiency->SetBranchAddress("tot_veto_ad4",&tot_veto_ad4);

		tr_efficiency->GetEntry(0);

		cout << "DAQ time for run" << run_num << " is: " << total_DAQ_ad2 << "s" << endl;

		if(total_DAQ_ad2 == 0 && total_DAQ_ad1 == 0){
			cout << "DAQ time = 0...?!?!?!" << endl;
			run_order += 1;
			runFile->Close();
			continue;
		}

	/*	for(int iad=0; iad<4; iad++){
			accidentals[iad] = 0;
		}*/


		for(int iad=0; iad<maxAD; iad++){
			if(iad == 0){
				if(total_DAQ_ad1 == 0) efficiency[iad] = 0;
				else efficiency[iad] = tot_live_ad1/total_DAQ_ad1;
			}
			if(iad == 1){
				if(total_DAQ_ad2 == 0) efficiency[iad] = 0;
				else efficiency[iad] = tot_live_ad2/total_DAQ_ad2;
			}
			if(hall_num == 3 && iad == 2){
				if(total_DAQ_ad3 == 0) efficiency[iad] = 0;
				else efficiency[iad] = tot_live_ad3/total_DAQ_ad3;
			}
			if(hall_num == 3 && iad == 3){
				if(total_DAQ_ad4 == 0) efficiency[iad] = 0;
				else efficiency[iad] = tot_live_ad4/total_DAQ_ad4;
			}
			h_muon_efficiency[iad]->Fill(run_order, efficiency[iad]);

			char name[64];
			sprintf(name,"h_ibd_energy_before_ad%d",iad+1);
			TH2F *h_run_ibd_energy_before = (TH2F*)runFile->Get(name);
			h_total_ibd_energy_before[iad]->Add(h_run_ibd_energy_before);

			sprintf(name,"h_ibd_energy_1m_ad%d",iad+1);
			TH2F *h_run_ibd_energy_1m = (TH2F*)runFile->Get(name);
			h_total_ibd_energy_1m[iad]->Add(h_run_ibd_energy_1m);

			sprintf(name,"h_ibd_distVStime_ad%d",iad+1);
			TH2F *h_run_ibd_distVStime = (TH2F*)runFile->Get(name);
			h_total_ibd_distVStime[iad]->Add(h_run_ibd_distVStime);

			sprintf(name,"h_ibd_distVStime_Ep35_ad%d",iad+1);
			TH2F *h_run_ibd_distVStime_Ep35 = (TH2F*)runFile->Get(name);
			h_total_ibd_distVStime_Ep35[iad]->Add(h_run_ibd_distVStime_Ep35);

			sprintf(name,"h_ibd_promptVStime_ad%d",iad+1);
			TH2F *h_run_ibd_promptVStime = (TH2F*)runFile->Get(name);
			h_total_ibd_promptVStime[iad]->Add(h_run_ibd_promptVStime);

			sprintf(name,"h_ibd_promptVStime_DT800_ad%d",iad+1);
			TH2F *h_run_ibd_promptVStime_DT800 = (TH2F*)runFile->Get(name);
			h_total_ibd_promptVStime_DT800[iad]->Add(h_run_ibd_promptVStime_DT800);

			sprintf(name,"h_ibd_DT_ad%d",iad+1);
			TH1D *h_run_ibd_DT = (TH1D*)runFile->Get(name);
			h_total_ibd_DT[iad]->Add(h_run_ibd_DT);

			sprintf(name,"h_ibd_DT_Ep35_ad%d",iad+1);
			TH1D *h_run_ibd_DT_Ep35 = (TH1D*)runFile->Get(name);
			h_total_ibd_DT_Ep35[iad]->Add(h_run_ibd_DT_Ep35);

			ibdEntries[iad] = h_run_ibd_energy_before->GetEntries();
			if(iad == 0){
				if(tot_live_ad1 != 0) h_IBDperTime[iad]->Fill(run_order, ibdEntries[iad]/tot_live_ad1);
			}
			if(iad == 1){
				if(tot_live_ad2 !=0) h_IBDperTime[iad]->Fill(run_order, ibdEntries[iad]/tot_live_ad2);
			}
			if(hall_num == 3 && iad == 2){
				if(tot_live_ad3 != 0) h_IBDperTime[iad]->Fill(run_order, ibdEntries[iad]/tot_live_ad3);
			}
			if(hall_num == 3 && iad == 3){
				if(tot_live_ad4 != 0) h_IBDperTime[iad]->Fill(run_order, ibdEntries[iad]/tot_live_ad4);
			}

			sprintf(name,"h_ibd_energy_after_ad%d",iad+1);
			TH2F *h_run_ibd_energy_after = (TH2F*)runFile->Get(name);
			h_total_ibd_energy_after[iad]->Add(h_run_ibd_energy_after);

			sprintf(name,"h_locations_before_ad%d",iad+1);
			TH2F *h_run_locations_before = (TH2F*)runFile->Get(name);
			h_total_locations_before[iad]->Add(h_run_locations_before);

			sprintf(name,"h_locations_after_ad%d",iad+1);
			TH2F *h_run_locations_after = (TH2F*)runFile->Get(name);
			h_total_locations_after[iad]->Add(h_run_locations_after);

			sprintf(name,"h_plocations_before_ad%d",iad+1);
			TH2F *h_run_plocations_before = (TH2F*)runFile->Get(name);
			h_total_plocations_before[iad]->Add(h_run_plocations_before);

			sprintf(name,"h_plocations_after_ad%d",iad+1);
			TH2F *h_run_plocations_after = (TH2F*)runFile->Get(name);
			h_total_plocations_after[iad]->Add(h_run_plocations_after);

			sprintf(name,"h_dlocations_before_ad%d",iad+1);
			TH2F *h_run_dlocations_before = (TH2F*)runFile->Get(name);
			h_total_dlocations_before[iad]->Add(h_run_dlocations_before);

			sprintf(name,"h_dlocations_after_ad%d",iad+1);
			TH2F *h_run_dlocations_after = (TH2F*)runFile->Get(name);
			h_total_dlocations_after[iad]->Add(h_run_dlocations_after);

			sprintf(name,"h_delayed_energy_before_ad%d",iad+1);
			TH1F *h_run_delayed_energy_before = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_before[iad]->Add(h_run_delayed_energy_before);

			for(int iz = 0; iz < NzBins; iz++){
				TH1F *h_run_delayed_energy_before_z = (TH1F*)runFile->Get(Form("h_delayed_energy_before_z_ad%d_iz%d", iad+1, iz+1));
				h_total_delayed_energy_before_z[iad][iz]->Add(h_run_delayed_energy_before_z);

				TH1F *h_run_delayed_energy_DT800_z = (TH1F*)runFile->Get(Form("h_delayed_energy_DT800_z_ad%d_iz%d", iad+1, iz+1));
				h_total_delayed_energy_DT800_z[iad][iz]->Add(h_run_delayed_energy_DT800_z);
			}

			for(int ir2 = 0; ir2 < NzBins; ir2++){
				TH1F *h_run_delayed_energy_before_r2 = (TH1F*)runFile->Get(Form("h_delayed_energy_before_r2_ad%d_ir2%d", iad+1, ir2+1));
				h_total_delayed_energy_before_r2[iad][ir2]->Add(h_run_delayed_energy_before_r2);

				TH1F *h_run_delayed_energy_DT800_r2 = (TH1F*)runFile->Get(Form("h_delayed_energy_DT800_r2_ad%d_ir2%d", iad+1, ir2+1));
				h_total_delayed_energy_DT800_r2[iad][ir2]->Add(h_run_delayed_energy_DT800_r2);
			}

			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < NzBins; ir2++){
					TH1F *h_run_delayed_energy_before_zVSr2 = (TH1F*)runFile->Get(Form("h_delayed_energy_before_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1,iz+1));
					h_total_delayed_energy_before_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_before_zVSr2);

					TH1F *h_run_delayed_energy_DT800_zVSr2 = (TH1F*)runFile->Get(Form("h_delayed_energy_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1,iz+1));
					h_total_delayed_energy_DT800_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_DT800_zVSr2);
				}
			}

			sprintf(name,"h_delayed_energy_DT800_ad%d",iad+1);
			TH1F *h_run_delayed_energy_DT800 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_DT800[iad]->Add(h_run_delayed_energy_DT800);

			sprintf(name,"h_delayed_energy_fine_before_ad%d",iad+1);
			TH1F *h_run_delayed_energy_fine_before = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_fine_before[iad]->Add(h_run_delayed_energy_fine_before);

			sprintf(name,"h_delayed_energy_fine_DT800_ad%d",iad+1);
			TH1F *h_run_delayed_energy_fine_DT800 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_fine_DT800[iad]->Add(h_run_delayed_energy_fine_DT800);

				sprintf(name,"h_delayed_energy_fine_Ep35_ad%d",iad+1);
				TH1F *h_run_delayed_energy_fine_Ep35 = (TH1F*)runFile->Get(name);
				h_total_delayed_energy_fine_Ep35[iad]->Add(h_run_delayed_energy_fine_Ep35);

				sprintf(name,"h_delayed_energy_fine_DT800_Ep35_ad%d",iad+1);
				TH1F *h_run_delayed_energy_fine_DT800_Ep35 = (TH1F*)runFile->Get(name);
				h_total_delayed_energy_fine_DT800_Ep35[iad]->Add(h_run_delayed_energy_fine_DT800_Ep35);

			sprintf(name,"h_delayed_energy_largeDist_ad%d",iad+1);
			TH1F *h_run_delayed_energy_largeDist = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_largeDist[iad]->Add(h_run_delayed_energy_largeDist);

			sprintf(name,"h_delayed_energy_after_ad%d",iad+1);
			TH1F *h_run_delayed_energy_after = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_after[iad]->Add(h_run_delayed_energy_after);

			sprintf(name,"h_prompt_energy_before_ad%d",iad+1);
			TH1F *h_run_prompt_energy_before = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_before[iad]->Add(h_run_prompt_energy_before);

			sprintf(name,"h_prompt_energy_DT800_ad%d",iad+1);
			TH1F *h_run_prompt_energy_DT800 = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_DT800[iad]->Add(h_run_prompt_energy_DT800);

			sprintf(name,"h_prompt_energy_DT800_3sig_ad%d",iad+1);
			TH1F *h_run_prompt_energy_DT800_3sig = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_DT800_3sig[iad]->Add(h_run_prompt_energy_DT800_3sig);

			sprintf(name,"h_prompt_energy_largeDist_ad%d",iad+1);
			TH1F *h_run_prompt_energy_largeDist = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_largeDist[iad]->Add(h_run_prompt_energy_largeDist);

			sprintf(name,"h_prompt_energy_after_ad%d",iad+1);
			TH1F *h_run_prompt_energy_after = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_after[iad]->Add(h_run_prompt_energy_after);

			sprintf(name,"h_delta_time_before_ad%d",iad+1);
			TH1D *h_run_delta_time_before = (TH1D*)runFile->Get(name);
			h_total_delta_time_before[iad]->Add(h_run_delta_time_before);

			sprintf(name,"h_delta_time_finer_ad%d",iad+1);
			TH1D *h_run_delta_time_finer = (TH1D*)runFile->Get(name);
			h_total_delta_time_finer[iad]->Add(h_run_delta_time_finer);

			sprintf(name,"h_delta_time_1us_ad%d",iad+1);
			TH1D *h_run_delta_time_1us = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us[iad]->Add(h_run_delta_time_1us);

			sprintf(name,"h_delta_time_1us_2m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_2m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_2m[iad]->Add(h_run_delta_time_1us_2m);

			sprintf(name,"h_delta_time_1us_15m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_15m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_15m[iad]->Add(h_run_delta_time_1us_15m);

			sprintf(name,"h_delta_time_1us_1m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_1m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_1m[iad]->Add(h_run_delta_time_1us_1m);

			sprintf(name,"h_delta_time_1us_075m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_075m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_075m[iad]->Add(h_run_delta_time_1us_075m);

			sprintf(name,"h_delta_time_1us_05m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_05m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_05m[iad]->Add(h_run_delta_time_1us_05m);


			sprintf(name,"h_delta_time_1us_gdls_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_gdls = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_gdls[iad]->Add(h_run_delta_time_1us_gdls);

			sprintf(name,"h_delta_time_1us_gdls_1m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_gdls_1m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_gdls_1m[iad]->Add(h_run_delta_time_1us_gdls_1m);

			sprintf(name,"h_delta_time_1us_gdls_075m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_gdls_075m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_gdls_075m[iad]->Add(h_run_delta_time_1us_gdls_075m);

			sprintf(name,"h_delta_time_1us_gdls_05m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_gdls_05m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_gdls_05m[iad]->Add(h_run_delta_time_1us_gdls_05m);


			sprintf(name,"h_delta_time_1us_ls_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_ls = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_ls[iad]->Add(h_run_delta_time_1us_ls);

			sprintf(name,"h_delta_time_1us_ls_1m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_ls_1m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_ls_1m[iad]->Add(h_run_delta_time_1us_ls_1m);

			sprintf(name,"h_delta_time_1us_ls_075m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_ls_075m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_ls_075m[iad]->Add(h_run_delta_time_1us_ls_075m);

			sprintf(name,"h_delta_time_1us_ls_05m_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_ls_05m = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_ls_05m[iad]->Add(h_run_delta_time_1us_ls_05m);

			sprintf(name,"h_delta_time_1us_largeDist_ad%d",iad+1);
			TH1D *h_run_delta_time_1us_largeDist = (TH1D*)runFile->Get(name);
			h_total_delta_time_1us_largeDist[iad]->Add(h_run_delta_time_1us_largeDist);


			sprintf(name,"h_delta_time_after_ad%d",iad+1);
			TH1D *h_run_delta_time_after = (TH1D*)runFile->Get(name);
			h_total_delta_time_after[iad]->Add(h_run_delta_time_after);

			sprintf(name,"h_distance_before_ad%d",iad+1);
			TH1F *h_run_distance_before = (TH1F*)runFile->Get(name);
			h_total_distance_before[iad]->Add(h_run_distance_before);

			sprintf(name,"h_distance_before_Ep35_ad%d",iad+1);
			TH1F *h_run_distance_before_Ep35 = (TH1F*)runFile->Get(name);
			h_total_distance_before_Ep35[iad]->Add(h_run_distance_before_Ep35);

			sprintf(name,"h_distance_before_400_ad%d",iad+1);
			TH1F *h_run_distance_before_400 = (TH1F*)runFile->Get(name);
			h_total_distance_before_400[iad]->Add(h_run_distance_before_400);

			sprintf(name,"h_distance_before_600_ad%d",iad+1);
			TH1F *h_run_distance_before_600 = (TH1F*)runFile->Get(name);
			h_total_distance_before_600[iad]->Add(h_run_distance_before_600);

			sprintf(name,"h_distance_before_800_ad%d",iad+1);
			TH1F *h_run_distance_before_800 = (TH1F*)runFile->Get(name);
			h_total_distance_before_800[iad]->Add(h_run_distance_before_800);

/*				int startBin = 0;
				for(int ibin=0; ibin<702; ibin++){
					if(h_total_distance_before[iad]->GetBinCenter(ibin)>2){
						startBin = ibin;
						break;
					}
				}*/
				
	//	accidentals[iad] = (h_total_distance_before[iad]->Integral(startBin,701));

			sprintf(name,"h_distance_after_ad%d",iad+1);
			TH1F *h_run_distance_after = (TH1F*)runFile->Get(name);
			h_total_distance_after[iad]->Add(h_run_distance_after);

			sprintf(name, "h_delayed_energy_scaled_15_22_ad%d", iad+1);
			TH1F *h_run_delayed_15_22 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_15_22[iad]->Add(h_run_delayed_15_22);

			sprintf(name, "h_delayed_energy_scaled_22_29_ad%d", iad+1);
			TH1F *h_run_delayed_22_29 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_22_29[iad]->Add(h_run_delayed_22_29);

			sprintf(name, "h_delayed_energy_scaled_29_36_ad%d", iad+1);
			TH1F *h_run_delayed_29_36 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_29_36[iad]->Add(h_run_delayed_29_36);

			sprintf(name, "h_delayed_energy_scaled_36_43_ad%d", iad+1);
			TH1F *h_run_delayed_36_43 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_36_43[iad]->Add(h_run_delayed_36_43);

			sprintf(name, "h_delayed_energy_scaled_43_50_ad%d", iad+1);
			TH1F *h_run_delayed_43_50 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_43_50[iad]->Add(h_run_delayed_43_50);

			sprintf(name, "h_delayed_energy_vs_distance_ad%d", iad+1);
			TH1F *h_run_delayed_vs_dist = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_vs_distance[iad]->Add(h_run_delayed_vs_dist);

			sprintf(name, "h_prompt_energy_vs_distance_ad%d", iad+1);
			TH1F *h_run_prompt_vs_dist = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_vs_distance[iad]->Add(h_run_prompt_vs_dist);

			sprintf(name, "h_delayed_energy_scaled_p35_ad%d", iad+1);
			TH1F *h_run_delayed_p35 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_p35[iad]->Add(h_run_delayed_p35);

			for(int dist = 0; dist < 6; dist++){
				sprintf(name, "h_delayed_energy_scaled_dist%d_ad%d",dist, iad+1);
				TH1F *h_run_delayed_dist = (TH1F*)runFile->Get(name);
				h_total_delayed_energy_dist[iad][dist]->Add(h_run_delayed_dist);
				if(dist >= 2) h_total_delayed_energy_dist_2plus[iad]->Add(h_run_delayed_dist);

				sprintf(name, "h_prompt_energy_scaled_dist%d_ad%d",dist, iad+1);
				TH1F *h_run_prompt_dist = (TH1F*)runFile->Get(name);
				h_total_prompt_energy_dist[iad][dist]->Add(h_run_prompt_dist);
				if(dist >= 2) h_total_prompt_energy_dist_2plus[iad]->Add(h_run_prompt_dist);
			}


/*			sprintf(name,"h_GdLSdistance_before_ad%d",iad+1);
			TH1F *h_run_GdLSdistance_before = (TH1F*)runFile->Get(name);
			h_total_GdLSdistance_before[iad]->Add(h_run_GdLSdistance_before);

			sprintf(name,"h_GdLSdistance_after_ad%d",iad+1);
			TH1F *h_run_GdLSdistance_after = (TH1F*)runFile->Get(name);
			h_total_GdLSdistance_after[iad]->Add(h_run_GdLSdistance_after);

			sprintf(name,"h_GdLS_ibd_energy_before_ad%d",iad+1);
			TH1F *h_run_GdLS_ibd_energy_before = (TH1F*)runFile->Get(name);
			h_total_GdLS_ibd_energy_before[iad]->Add(h_run_GdLS_ibd_energy_before);

			sprintf(name,"h_GdLS_ibd_energy_after_ad%d",iad+1);
			TH1F *h_run_GdLS_ibd_energy_after = (TH1F*)runFile->Get(name);
			h_total_GdLS_ibd_energy_after[iad]->Add(h_run_GdLS_ibd_energy_after);*/

		}

//		accCountsFile << run_num << "\t" << accidentals[0] << "\t" << accidentals[1] << "\t" << accidentals[2] << "\t" << accidentals[3] << endl;

		/*accCounts_ad1 = accidentals[0];
		accCounts_ad2 = accidentals[1];
		accCounts_ad3 = accidentals[2];
		accCounts_ad4 = accidentals[3];
		tr_2m->Fill();*/

		cout << "Done with run_order #" << run_order << " out of " << nRuns << "runs." << endl;

		runFile->Close();
		//if(run_order >= 99) break;
		run_order += 1;
	}

//	accCountsFile.close();

//Rebin and add error bars to 1us time plot
	double timeCounts = 0;
	for(int iad = 0; iad<maxAD; iad++){
		for(int ibin=0; ibin <= 402; ibin++){
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us[iad]->GetBinContent(ibin);
			h_total_delta_time_1us[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_2m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_2m[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_15m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_15m[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_1m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_1m[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_075m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_075m[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_05m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_05m[iad]->SetBinError(ibin, sqrt(timeCounts));

			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_gdls[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_gdls[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_gdls_1m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_gdls_1m[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_gdls_075m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_gdls_075m[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_gdls_05m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_gdls_05m[iad]->SetBinError(ibin, sqrt(timeCounts));

			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_ls[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_ls[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_ls_1m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_ls_1m[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_ls_075m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_ls_075m[iad]->SetBinError(ibin, sqrt(timeCounts));
			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_ls_05m[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_ls_05m[iad]->SetBinError(ibin, sqrt(timeCounts));

			timeCounts = 0;
			timeCounts = h_total_delta_time_1us_largeDist[iad]->GetBinContent(ibin);
			h_total_delta_time_1us_largeDist[iad]->SetBinError(ibin, sqrt(timeCounts));
		}
	}

/*	TF1* time_fit=new TF1("time_fit","[0]+[1]*TMath::Exp(-x/[2])+[3]*TMath::Exp(-x/[4])",1,400);
		time_fit->SetParameter(1,333);
		time_fit->SetParameter(2,210);
		time_fit->SetParameter(0,268);
		time_fit->SetParameter(3,378);
		time_fit->SetParameter(4,29);

	TF1* time_fit1=new TF1("time_fit1","[0]+[1]*TMath::Exp(-x/[2])+[3]*TMath::Exp(-x/[4])",1,400);
		time_fit1->SetParameter(1,333);
		time_fit1->SetParameter(2,210);
		time_fit1->SetParameter(0,268);
		time_fit1->SetParameter(3,378);
		time_fit1->SetParameter(4,29);

	TF1* time_fit075=new TF1("time_fit075","[0]+[1]*TMath::Exp(-x/[2])+[3]*TMath::Exp(-x/[4])",1,400);
		time_fit075->SetParameter(1,333);
		time_fit075->SetParameter(2,210);
		time_fit075->SetParameter(0,268);
		time_fit075->SetParameter(3,378);
		time_fit075->SetParameter(4,29);

	TF1* time_fit05=new TF1("time_fit05","[0]+[1]*TMath::Exp(-x/[2])+[3]*TMath::Exp(-x/[4])",1,400);
		time_fit05->SetParameter(1,333);
		time_fit05->SetParameter(2,210);
		time_fit05->SetParameter(0,268);
		time_fit05->SetParameter(3,378);
		time_fit05->SetParameter(4,29);

	double offset = 0;*/


        char outputname[64];
//	sprintf(outputname,"./IBDs/TotaledPlots_TcLong_EH%d_Ep2.root",hall_num);
	sprintf(outputname,"./IBDs/TotaledPlots_NU_EH%d_%d.root",hall_num,pd_window_microsec);
//	sprintf(outputname,"./IBDs/TotaledPlots_4sigma_EH%d.root",hall_num);
	TFile* outfile=new TFile(outputname, "RECREATE");
		outfile->cd();
	//	tr_2m->Write();
		for(int iad=0; iad<maxAD; ++iad){
			h_total_ibd_energy_before[iad]->SetOption("COLZ");
			h_total_ibd_energy_before[iad]->SetStats(0);
			h_total_ibd_energy_before[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_ibd_energy_before[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_ibd_energy_before[iad]->Write();

			h_total_ibd_energy_1m[iad]->SetOption("COLZ");
			h_total_ibd_energy_1m[iad]->SetStats(0);
			h_total_ibd_energy_1m[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_ibd_energy_1m[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_ibd_energy_1m[iad]->Write();

			h_total_ibd_distVStime[iad]->SetStats(0);
			h_total_ibd_distVStime[iad]->SetOption("COLZ");
			h_total_ibd_distVStime[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_ibd_distVStime[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_ibd_distVStime[iad]->Write();

			h_total_ibd_distVStime_Ep35[iad]->SetStats(0);
			h_total_ibd_distVStime_Ep35[iad]->SetOption("COLZ");
			h_total_ibd_distVStime_Ep35[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_ibd_distVStime_Ep35[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_ibd_distVStime_Ep35[iad]->Write();

		//		h_total_ibd_DT[iad]->SetStats(0);
				h_total_ibd_DT[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_ibd_DT[iad]->GetYaxis()->SetTitle("Counts");
				h_total_ibd_DT[iad]->Write();

		//		h_total_ibd_DT_Ep35[iad]->SetStats(0);
				h_total_ibd_DT_Ep35[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_ibd_DT_Ep35[iad]->GetYaxis()->SetTitle("Counts");
				h_total_ibd_DT_Ep35[iad]->Write();

			h_total_ibd_energy_after[iad]->SetOption("COLZ");
			h_total_ibd_energy_after[iad]->SetStats(0);
			h_total_ibd_energy_after[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_ibd_energy_after[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_ibd_energy_after[iad]->Write();

			h_total_locations_before[iad]->SetOption("COLZ");
			h_total_locations_before[iad]->SetStats(0);
			h_total_locations_before[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_total_locations_before[iad]->GetYaxis()->SetTitle("z [m]");
			h_total_locations_before[iad]->Write();

			h_total_locations_after[iad]->SetOption("COLZ");
			h_total_locations_after[iad]->SetStats(0);
			h_total_locations_after[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_total_locations_after[iad]->GetYaxis()->SetTitle("z [m]");
			h_total_locations_after[iad]->Write();

			h_total_plocations_before[iad]->SetOption("COLZ");
			h_total_plocations_before[iad]->SetStats(0);
			h_total_plocations_before[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_total_plocations_before[iad]->GetYaxis()->SetTitle("z [m]");
			h_total_plocations_before[iad]->Write();

			h_total_plocations_after[iad]->SetOption("COLZ");
			h_total_plocations_after[iad]->SetStats(0);
			h_total_plocations_after[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_total_plocations_after[iad]->GetYaxis()->SetTitle("z [m]");
			h_total_plocations_after[iad]->Write();

			h_total_dlocations_before[iad]->SetOption("COLZ");
			h_total_dlocations_before[iad]->SetStats(0);
			h_total_dlocations_before[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_total_dlocations_before[iad]->GetYaxis()->SetTitle("z [m]");
			h_total_dlocations_before[iad]->Write();

			h_total_dlocations_after[iad]->SetOption("COLZ");
			h_total_dlocations_after[iad]->SetStats(0);
			h_total_dlocations_after[iad]->GetXaxis()->SetTitle("r^{2} [m^{2}]");
			h_total_dlocations_after[iad]->GetYaxis()->SetTitle("z [m]");
			h_total_dlocations_after[iad]->Write();

			h_total_delayed_energy_before[iad]->SetStats(0);
			h_total_delayed_energy_before[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_before[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_before[iad]->Write();

			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_before_z[iad][iz]->SetStats(0);
				h_total_delayed_energy_before_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_before_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_before_z[iad][iz]->Write();

				h_total_delayed_energy_DT800_z[iad][iz]->SetStats(0);
				h_total_delayed_energy_DT800_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_DT800_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_DT800_z[iad][iz]->Write();
			}

			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_before_r2[iad][ir2]->SetStats(0);
				h_total_delayed_energy_before_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_before_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_before_r2[iad][ir2]->Write();

				h_total_delayed_energy_DT800_r2[iad][ir2]->SetStats(0);
				h_total_delayed_energy_DT800_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_DT800_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_DT800_r2[iad][ir2]->Write();
			}

			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_before_zVSr2[iad][ir2][iz]->SetStats(0);
					h_total_delayed_energy_before_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
					h_total_delayed_energy_before_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_total_delayed_energy_before_zVSr2[iad][ir2][iz]->Write();

					h_total_delayed_energy_DT800_zVSr2[iad][ir2][iz]->SetStats(0);
					h_total_delayed_energy_DT800_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
					h_total_delayed_energy_DT800_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_total_delayed_energy_DT800_zVSr2[iad][ir2][iz]->Write();
				}
			}

			h_total_delayed_energy_DT800[iad]->SetStats(0);
			h_total_delayed_energy_DT800[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_DT800[iad]->Write();

			h_total_delayed_energy_fine_before[iad]->SetStats(0);
			h_total_delayed_energy_fine_before[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_before[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_before[iad]->Write();

			h_total_delayed_energy_fine_DT800[iad]->SetStats(0);
			h_total_delayed_energy_fine_DT800[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DT800[iad]->Write();

			h_total_delayed_energy_fine_Ep35[iad]->SetStats(0);
			h_total_delayed_energy_fine_Ep35[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_Ep35[iad]->Write();

			h_total_delayed_energy_fine_DT800_Ep35[iad]->SetStats(0);
			h_total_delayed_energy_fine_DT800_Ep35[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DT800_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DT800_Ep35[iad]->Write();

			h_total_delayed_energy_largeDist[iad]->SetStats(0);
			h_total_delayed_energy_largeDist[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_largeDist[iad]->Write();

			h_total_delayed_energy_after[iad]->SetStats(0);
			h_total_delayed_energy_after[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_after[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_after[iad]->Write();

			h_total_prompt_energy_before[iad]->SetStats(0);
			h_total_prompt_energy_before[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_before[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_before[iad]->Write();

			h_total_prompt_energy_DT800[iad]->SetStats(0);
			h_total_prompt_energy_DT800[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800[iad]->Write();

			h_total_prompt_energy_DT800_3sig[iad]->SetStats(0);
			h_total_prompt_energy_DT800_3sig[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800_3sig[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800_3sig[iad]->Write();

			h_total_prompt_energy_largeDist[iad]->SetStats(0);
			h_total_prompt_energy_largeDist[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_largeDist[iad]->Write();

			h_total_prompt_energy_after[iad]->SetStats(0);
			h_total_prompt_energy_after[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_after[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_after[iad]->Write();

			h_total_delta_time_before[iad]->SetStats(0);
			h_total_delta_time_before[iad]->SetMinimum(0);
			h_total_delta_time_before[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_before[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delta_time_before[iad]->Write();

			h_total_delta_time_finer[iad]->SetStats(0);
			h_total_delta_time_finer[iad]->SetMinimum(0);
			h_total_delta_time_finer[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_finer[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delta_time_finer[iad]->Write();


		//	cout << "no dist cut:" << endl;
			//h_total_delta_time_1us[iad]->SetStats(0);
			h_total_delta_time_1us[iad]->SetMinimum(0);
			h_total_delta_time_1us[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us[iad]->Fit("time_fit","R");
			h_total_delta_time_1us[iad]->Write();

		//	cout << "2m cut:" << endl;
			//h_total_delta_time_1us_2m[iad]->SetStats(0);
			h_total_delta_time_1us_2m[iad]->SetMinimum(0);
			h_total_delta_time_1us_2m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_2m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_2m[iad]->Fit("time_fit1","R");
			h_total_delta_time_1us_2m[iad]->Write();

		//	cout << "1.5m cut:" << endl;
			//h_total_delta_time_1us_15m[iad]->SetStats(0);
			h_total_delta_time_1us_15m[iad]->SetMinimum(0);
			h_total_delta_time_1us_15m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_15m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_15m[iad]->Fit("time_fit1","R");
			h_total_delta_time_1us_15m[iad]->Write();

		//	cout << "1m cut:" << endl;
			//h_total_delta_time_1us_1m[iad]->SetStats(0);
			h_total_delta_time_1us_1m[iad]->SetMinimum(0);
			h_total_delta_time_1us_1m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_1m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_1m[iad]->Fit("time_fit1","R");
			h_total_delta_time_1us_1m[iad]->Write();

		//	cout << "0.75m cut:" << endl;
			//h_total_delta_time_1us_075m[iad]->SetStats(0);
			h_total_delta_time_1us_075m[iad]->SetMinimum(0);
			h_total_delta_time_1us_075m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_075m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_075m[iad]->Fit("time_fit075","R");
			h_total_delta_time_1us_075m[iad]->Write();

		//	cout << "0.5m cut:" << endl;
			//h_total_delta_time_1us_05m[iad]->SetStats(0);
			h_total_delta_time_1us_05m[iad]->SetMinimum(0);
			h_total_delta_time_1us_05m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_05m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_05m[iad]->Fit("time_fit05","R");
			h_total_delta_time_1us_05m[iad]->Write();

		//	cout << "no dist cut:" << endl;
			//h_total_delta_time_1us_gdls[iad]->SetStats(0);
			h_total_delta_time_1us_gdls[iad]->SetMinimum(0);
			h_total_delta_time_1us_gdls[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_gdls[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_gdls[iad]->Fit("time_fit","R");
			h_total_delta_time_1us_gdls[iad]->Write();

		//	cout << "1m cut:" << endl;
			//h_total_delta_time_1us_gdls_1m[iad]->SetStats(0);
			h_total_delta_time_1us_gdls_1m[iad]->SetMinimum(0);
			h_total_delta_time_1us_gdls_1m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_gdls_1m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_gdls_1m[iad]->Fit("time_fit1","R");
			h_total_delta_time_1us_gdls_1m[iad]->Write();

		//	cout << "0.75m cut:" << endl;
			//h_total_delta_time_1us_gdls_075m[iad]->SetStats(0);
			h_total_delta_time_1us_gdls_075m[iad]->SetMinimum(0);
			h_total_delta_time_1us_gdls_075m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_gdls_075m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_gdls_075m[iad]->Fit("time_fit075","R");
			h_total_delta_time_1us_gdls_075m[iad]->Write();

		//	cout << "0.5m cut:" << endl;
			//h_total_delta_time_1us_gdls_05m[iad]->SetStats(0);
			h_total_delta_time_1us_gdls_05m[iad]->SetMinimum(0);
			h_total_delta_time_1us_gdls_05m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_gdls_05m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_gdls_05m[iad]->Fit("time_fit05","R");
			h_total_delta_time_1us_gdls_05m[iad]->Write();


		//	cout << "no dist cut:" << endl;
			//h_total_delta_time_1us_ls[iad]->SetStats(0);
			h_total_delta_time_1us_ls[iad]->SetMinimum(0);
			h_total_delta_time_1us_ls[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_ls[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_ls[iad]->Fit("time_fit","R");
			h_total_delta_time_1us_ls[iad]->Write();

		//	cout << "1m cut:" << endl;
			//h_total_delta_time_1us_ls_1m[iad]->SetStats(0);
			h_total_delta_time_1us_ls_1m[iad]->SetMinimum(0);
			h_total_delta_time_1us_ls_1m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_ls_1m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_ls_1m[iad]->Fit("time_fit1","R");
			h_total_delta_time_1us_ls_1m[iad]->Write();

		//	cout << "0.75m cut:" << endl;
			//h_total_delta_time_1us_ls_075m[iad]->SetStats(0);
			h_total_delta_time_1us_ls_075m[iad]->SetMinimum(0);
			h_total_delta_time_1us_ls_075m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_ls_075m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_ls_075m[iad]->Fit("time_fit075","R");
			h_total_delta_time_1us_ls_075m[iad]->Write();

		//	cout << "0.5m cut:" << endl;
			//h_total_delta_time_1us_ls_05m[iad]->SetStats(0);
			h_total_delta_time_1us_ls_05m[iad]->SetMinimum(0);
			h_total_delta_time_1us_ls_05m[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_ls_05m[iad]->GetYaxis()->SetTitle("Counts");
		//	h_total_delta_time_1us_ls_05m[iad]->Fit("time_fit05","R");
			h_total_delta_time_1us_ls_05m[iad]->Write();

			h_total_delta_time_1us_largeDist[iad]->SetMinimum(0);
			h_total_delta_time_1us_largeDist[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_1us_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delta_time_1us_largeDist[iad]->Write();

			h_total_delta_time_after[iad]->SetStats(0);
			h_total_delta_time_after[iad]->SetMinimum(0);
			h_total_delta_time_after[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_delta_time_after[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delta_time_after[iad]->Write();

			//h_total_distance_before[iad]->SetStats(0);
			h_total_distance_before[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_distance_before[iad]->GetYaxis()->SetTitle("Counts");
			h_total_distance_before[iad]->Write();

			//h_total_distance_before_Ep35[iad]->SetStats(0);
			h_total_distance_before_Ep35[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_distance_before_Ep35[iad]->GetYaxis()->SetTitle("Counts");
			h_total_distance_before_Ep35[iad]->Write();

			//h_total_distance_before_400[iad]->SetStats(0);
			h_total_distance_before_400[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_distance_before_400[iad]->GetYaxis()->SetTitle("Counts");
			h_total_distance_before_400[iad]->Write();

			//h_total_distance_before_600[iad]->SetStats(0);
			h_total_distance_before_600[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_distance_before_600[iad]->GetYaxis()->SetTitle("Counts");
			h_total_distance_before_600[iad]->Write();

			//h_total_distance_before_800[iad]->SetStats(0);
			h_total_distance_before_800[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_distance_before_800[iad]->GetYaxis()->SetTitle("Counts");
			h_total_distance_before_800[iad]->Write();

			//h_total_distance_after[iad]->SetStats(0);
			h_total_distance_after[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_distance_after[iad]->GetYaxis()->SetTitle("Counts");
			h_total_distance_after[iad]->Write();

			h_muon_efficiency[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_muon_efficiency[iad]->GetYaxis()->SetTitle("Muon Efficiency");
			h_muon_efficiency[iad]->Write();

			h_IBDperTime[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_IBDperTime[iad]->GetYaxis()->SetTitle("IBDs/livetime");
			h_IBDperTime[iad]->Write();

			h_total_delayed_energy_scaled_15_22[iad]->SetStats(0);
			h_total_delayed_energy_scaled_15_22[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_scaled_15_22[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_scaled_15_22[iad]->Write();

			h_total_delayed_energy_scaled_22_29[iad]->SetStats(0);
			h_total_delayed_energy_scaled_22_29[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_scaled_22_29[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_scaled_22_29[iad]->Write();

			h_total_delayed_energy_scaled_29_36[iad]->SetStats(0);
			h_total_delayed_energy_scaled_29_36[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_scaled_29_36[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_scaled_29_36[iad]->Write();

			h_total_delayed_energy_scaled_36_43[iad]->SetStats(0);
			h_total_delayed_energy_scaled_36_43[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_scaled_36_43[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_scaled_36_43[iad]->Write();

			h_total_delayed_energy_scaled_43_50[iad]->SetStats(0);
			h_total_delayed_energy_scaled_43_50[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_scaled_43_50[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_scaled_43_50[iad]->Write();

			h_total_delayed_energy_p35[iad]->SetStats(0);
			h_total_delayed_energy_p35[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_p35[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_p35[iad]->Write();

			h_total_delayed_energy_vs_distance[iad]->SetStats(0);
			h_total_delayed_energy_vs_distance[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_vs_distance[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_delayed_energy_vs_distance[iad]->SetOption("COLZ");
			h_total_delayed_energy_vs_distance[iad]->Write();

			h_total_prompt_energy_vs_distance[iad]->SetStats(0);
			h_total_prompt_energy_vs_distance[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_vs_distance[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_prompt_energy_vs_distance[iad]->SetOption("COLZ");
			h_total_prompt_energy_vs_distance[iad]->Write();

			h_total_ibd_promptVStime[iad]->SetStats(0);
			h_total_ibd_promptVStime[iad]->SetOption("COLZ");
			h_total_ibd_promptVStime[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_ibd_promptVStime[iad]->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_ibd_promptVStime[iad]->Write();

			h_total_ibd_promptVStime_DT800[iad]->SetStats(0);
			h_total_ibd_promptVStime_DT800[iad]->SetOption("COLZ");
			h_total_ibd_promptVStime_DT800[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_ibd_promptVStime_DT800[iad]->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_ibd_promptVStime_DT800[iad]->Write();

			for(int dist=0; dist<6; dist++){
				h_total_delayed_energy_dist[iad][dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_dist[iad][dist]->Write();
				h_total_prompt_energy_dist[iad][dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_dist[iad][dist]->Write();
			}

			h_total_delayed_energy_dist_2plus[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_dist_2plus[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_dist_2plus[iad]->Write();

			h_total_prompt_energy_dist_2plus[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_dist_2plus[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_dist_2plus[iad]->Write();

		}

	for(int iad = 0; iad < maxAD; iad++){
		cout << endl << "Total N_DC for EH" << EH << " AD" << iad+1 << ":\t" << h_total_prompt_energy_DT800_3sig[iad]->Integral() << endl;
	}

	outfile->Close();

	resub.close();

}
