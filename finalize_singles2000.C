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

double correct(double uncorrRate, double timing){
	double croot = 0;
	double rate = 0;
//	float timing = 1200.e-6;
//	float timing = 2800.e-6;
	croot = pow(27*pow(timing,4.)*uncorrRate-10*pow(timing,3)+3*sqrt(3)*sqrt(27*pow(timing,8)*pow(uncorrRate,2)-20*pow(timing,7)*uncorrRate+4*pow(timing,6)),1/3.);

	rate = -2/(3*croot) + croot/(3*pow(timing,2))+2/(3*timing);

//	cout << "Croot is: " << croot << endl;
//	rate = uncorrRate;
//	rate = (1-sqrt(1-4*timing*uncorrRate))/(2*timing);
//	cout << "Uncorrected Rate is: " << uncorrRate << endl;
//	cout << "Corrected Rate is: " << rate << endl;

	return rate;
}





void finalize(int hall_num, int pd_window_microsec){

	char hallList[64];
	sprintf(hallList, "./EH%druns_sync.txt",hall_num);
	FILE* runfile=fopen(hallList,"r");

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

	//Variables:
	int num_prompt[8];
	int num_delayed[8];
	double prompt_div[8];
	double delayed_div[8];
	float t_sec = 1200*1.e-6;
	double prompt_rate[8];
	double delayed_rate[8];
	double acc_rate[8];
	double acc_rate_400[8];
	double acc_rate_600[8];
	double acc_rate_800[8];
	double scale[8];
	double scale_400[8];
	double scale_600[8];
	double scale_800[8];
	double normScale[8];
	double normScale_400[8];
	double normScale_600[8];
	double normScale_800[8];
	double DTscale[8];
	double DTscale_Ep35[8];
	double acc_counts[8];
	double prompt_live[8];
	double prompt_DAQ[8];
	double delayed_live[8];
	double delayed_DAQ[8];
	double d_rate_error[8];
	double p_rate_error[8];
	int run[3];
	int run_num = 0;
	int EH = 0;


	for(int i=0; i<3; i++){
		run[i] = 0;
	}

	//Making Histograms:
		//BEFORE DISTANCE CUT HISTS
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

		TH1F* h_total_delayed_energy_fine_before[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_before_ad%d", iad+1);
			h_total_delayed_energy_fine_before[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_prompt_energy_before[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_before_ad%d", iad+1);
			h_total_prompt_energy_before[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_delayed_energy_DT800[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_DT800_ad%d", iad+1);
			h_total_delayed_energy_DT800[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DT800[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_ad%d", iad+1);
			h_total_delayed_energy_fine_DT800[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_prompt_energy_DT800[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_DT800_ad%d", iad+1);
			h_total_prompt_energy_DT800[iad]=new TH1F(name,name,113,0.7,12.);
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

		TH1F* h_total_delayed_energy_scaled_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_scaled_z[iad][iz]=new TH1F(Form("h_total_delayed_energy_scaled_z_ad%d_iz%d", iad+1, iz+1),Form("h_total_delayed_energy_scaled_z_ad%d_iz%d", iad+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_scaled_DT800_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_scaled_DT800_z[iad][iz]=new TH1F(Form("h_total_delayed_energy_scaled_DT800_z_ad%d_iz%d", iad+1, iz+1),Form("h_total_delayed_energy_scaled_DT800_z_ad%d_iz%d", iad+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_norm_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_norm_z[iad][iz]=new TH1F(Form("h_total_delayed_energy_norm_z_ad%d_iz%d", iad+1, iz+1),Form("h_total_delayed_energy_norm_z_ad%d_iz%d", iad+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_norm_DT800_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_norm_DT800_z[iad][iz]=new TH1F(Form("h_total_delayed_energy_norm_DT800_z_ad%d_iz%d", iad+1, iz+1),Form("h_total_delayed_energy_norm_DT800_z_ad%d_iz%d", iad+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_DTnorm_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_DTnorm_z[iad][iz]=new TH1F(Form("h_total_delayed_energy_DTnorm_z_ad%d_iz%d", iad+1, iz+1),Form("h_total_delayed_energy_DTnorm_z_ad%d_iz%d", iad+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_DTnorm_DT800_z[maxAD][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				h_total_delayed_energy_DTnorm_DT800_z[iad][iz]=new TH1F(Form("h_total_delayed_energy_DTnorm_DT800_z_ad%d_iz%d", iad+1, iz+1),Form("h_total_delayed_energy_DTnorm_DT800_z_ad%d_iz%d", iad+1, iz+1),230,0.7,3.);
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

		TH1F* h_total_delayed_energy_scaled_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_scaled_r2[iad][ir2]=new TH1F(Form("h_total_delayed_energy_scaled_r2_ad%d_ir2%d", iad+1, ir2+1),Form("h_total_delayed_energy_scaled_r2_ad%d_ir2%d", iad+1, ir2+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_scaled_DT800_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_scaled_DT800_r2[iad][ir2]=new TH1F(Form("h_total_delayed_energy_scaled_DT800_r2_ad%d_ir2%d", iad+1, ir2+1),Form("h_total_delayed_energy_scaled_DT800_r2_ad%d_ir2%d", iad+1, ir2+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_norm_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_norm_r2[iad][ir2]=new TH1F(Form("h_total_delayed_energy_norm_r2_ad%d_ir2%d", iad+1, ir2+1),Form("h_total_delayed_energy_norm_r2_ad%d_ir2%d", iad+1, ir2+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_norm_DT800_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_norm_DT800_r2[iad][ir2]=new TH1F(Form("h_total_delayed_energy_norm_DT800_r2_ad%d_ir2%d", iad+1, ir2+1),Form("h_total_delayed_energy_norm_DT800_r2_ad%d_ir2%d", iad+1, ir2+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_DTnorm_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_DTnorm_r2[iad][ir2]=new TH1F(Form("h_total_delayed_energy_DTnorm_r2_ad%d_ir2%d", iad+1, ir2+1),Form("h_total_delayed_energy_DTnorm_r2_ad%d_ir2%d", iad+1, ir2+1),230,0.7,3.);
			}
		}

		TH1F* h_total_delayed_energy_DTnorm_DT800_r2[maxAD][Nr2Bins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_total_delayed_energy_DTnorm_DT800_r2[iad][ir2]=new TH1F(Form("h_total_delayed_energy_DTnorm_DT800_r2_ad%d_ir2%d", iad+1, ir2+1),Form("h_total_delayed_energy_DTnorm_DT800_r2_ad%d_ir2%d", iad+1, ir2+1),230,0.7,3.);
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

		TH1F* h_total_delayed_energy_scaled_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_scaled_zVSr2[iad][ir2][iz]=new TH1F(Form("h_total_delayed_energy_scaled_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),Form("h_total_delayed_energy_scaled_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),230,0.7,3.);
				}
			}
		}

		TH1F* h_total_delayed_energy_scaled_DT800_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_scaled_DT800_zVSr2[iad][ir2][iz]=new TH1F(Form("h_total_delayed_energy_scaled_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),Form("h_total_delayed_energy_scaled_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),230,0.7,3.);
				}
			}
		}

		TH1F* h_total_delayed_energy_norm_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_norm_zVSr2[iad][ir2][iz]=new TH1F(Form("h_total_delayed_energy_norm_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),Form("h_total_delayed_energy_norm_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),230,0.7,3.);
				}
			}
		}

		TH1F* h_total_delayed_energy_norm_DT800_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_norm_DT800_zVSr2[iad][ir2][iz]=new TH1F(Form("h_total_delayed_energy_norm_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),Form("h_total_delayed_energy_norm_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),230,0.7,3.);
				}
			}
		}

		TH1F* h_total_delayed_energy_DTnorm_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_DTnorm_zVSr2[iad][ir2][iz]=new TH1F(Form("h_total_delayed_energy_DTnorm_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),Form("h_total_delayed_energy_DTnorm_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),230,0.7,3.);
				}
			}
		}

		TH1F* h_total_delayed_energy_DTnorm_DT800_zVSr2[maxAD][Nr2Bins][NzBins]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_total_delayed_energy_DTnorm_DT800_zVSr2[iad][ir2][iz]=new TH1F(Form("h_total_delayed_energy_DTnorm_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),Form("h_total_delayed_energy_DTnorm_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1),230,0.7,3.);
				}
			}
		}


			TH2F* h_total_acc_promptVStime_DTnorm[maxAD]; //prompt vs. time histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_acc_promptVStime_DTnorm_ad%d",iad+1);
				h_total_acc_promptVStime_DTnorm[iad]=new TH2F(name,name,1999,1,2000,113,0.7,12.);
			}

			TH2F* h_total_acc_promptVStime_DT800_DTnorm[maxAD]; //prompt vs. time histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_acc_promptVStime_DT800_DTnorm_ad%d",iad+1);
				h_total_acc_promptVStime_DT800_DTnorm[iad]=new TH2F(name,name,1999,1,2000,113,0.7,12.);
			}

		TH1F* h_total_prompt_energy_before_raw[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_before_raw_ad%d", iad+1);
			h_total_prompt_energy_before_raw[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_acc_distance_before[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_before_ad%d", iad+1);
			h_total_acc_distance_before[iad]=new TH1F(name,name,700,0,7.);
		}

		//AFTER DISTANCE CUT HISTOGRAMS
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

		TH1F* h_total_delayed_energy_after[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_after_ad%d", iad+1);
		//	h_total_delayed_energy_after[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_after[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_largeDist[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_largeDist_ad%d", iad+1);
		//	h_total_delayed_energy_largeDist[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_largeDist[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_shortDist[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_shortDist_ad%d", iad+1);
		//	h_total_delayed_energy_shortDist[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_shortDist[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_prompt_energy_after[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_after_ad%d", iad+1);
			h_total_prompt_energy_after[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_largeDist[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_largeDist_ad%d", iad+1);
			h_total_prompt_energy_largeDist[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_shortDist[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_shortDist_ad%d", iad+1);
			h_total_prompt_energy_shortDist[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_acc_distance_after[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_after_ad%d", iad+1);
			h_total_acc_distance_after[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_run_acc_distance_scaled[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_run_acc_distance_scaled_ad%d", iad+1);
			h_run_acc_distance_scaled[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_scaled[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_scaled_ad%d", iad+1);
			h_total_acc_distance_scaled[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_norm[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_norm_ad%d", iad+1);
			h_total_acc_distance_norm[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_3sig_scaled[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_3sig_scaled_ad%d", iad+1);
			h_total_acc_distance_3sig_scaled[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_3sig_norm[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_3sig_norm_ad%d", iad+1);
			h_total_acc_distance_3sig_norm[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_3sig_DTnorm[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_3sig_DTnorm_ad%d", iad+1);
			h_total_acc_distance_3sig_DTnorm[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_Ep35_scaled[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_Ep35_scaled_ad%d", iad+1);
			h_total_acc_distance_Ep35_scaled[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_Ep35_norm[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_Ep35_norm_ad%d", iad+1);
			h_total_acc_distance_Ep35_norm[iad]=new TH1F(name,name,700,0,7.);
		}

		//400 us distance plots
		TH1F* h_total_acc_distance_400_scaled[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_400_scaled_ad%d", iad+1);
			h_total_acc_distance_400_scaled[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_400_norm[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_400_norm_ad%d", iad+1);
			h_total_acc_distance_400_norm[iad]=new TH1F(name,name,700,0,7.);
		}

		//600 us distance plots
		TH1F* h_total_acc_distance_600_scaled[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_600_scaled_ad%d", iad+1);
			h_total_acc_distance_600_scaled[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_600_norm[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_600_norm_ad%d", iad+1);
			h_total_acc_distance_600_norm[iad]=new TH1F(name,name,700,0,7.);
		}

		//800 us distance plots
		TH1F* h_total_acc_distance_800_scaled[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_800_scaled_ad%d", iad+1);
			h_total_acc_distance_800_scaled[iad]=new TH1F(name,name,700,0,7.);
		}

		TH1F* h_total_acc_distance_800_norm[maxAD]; //distance histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_distance_800_norm_ad%d", iad+1);
			h_total_acc_distance_800_norm[iad]=new TH1F(name,name,700,0,7.);
		}




		TH1F* h_run_delayed_energy_scaled[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_run_delayed_energy_scaled_ad%d", iad+1);
		//	h_run_delayed_energy_scaled[iad]=new TH1F(name,name,150,1.5,3.);
			h_run_delayed_energy_scaled[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_run_delayed_energy_fine_scaled[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_run_delayed_energy_fine_scaled_ad%d", iad+1);
			h_run_delayed_energy_fine_scaled[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_run_prompt_energy_scaled[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_run_prompt_energy_scaled_ad%d", iad+1);
			h_run_prompt_energy_scaled[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_run_delayed_energy_DT800_scaled[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_run_delayed_energy_DT800_scaled_ad%d", iad+1);
			h_run_delayed_energy_DT800_scaled[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_run_delayed_energy_fine_DT800_scaled[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_run_delayed_energy_fine_DT800_scaled_ad%d", iad+1);
			h_run_delayed_energy_fine_DT800_scaled[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_run_prompt_energy_DT800_scaled[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_run_prompt_energy_DT800_scaled_ad%d", iad+1);
			h_run_prompt_energy_DT800_scaled[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_delayed_energy_raw[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_raw_ad%d", iad+1);
		//	h_total_delayed_energy_raw[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_raw[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_scaled_ad%d", iad+1);
		//	h_total_delayed_energy_scaled[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_norm[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_norm_ad%d", iad+1);
		//	h_total_delayed_energy_norm[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_norm[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_DTnorm[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_DTnorm_ad%d", iad+1);
			h_total_delayed_energy_DTnorm[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_DT800_scaled[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_DT800_scaled_ad%d", iad+1);
		//	h_total_delayed_energy_DT800_scaled[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_DT800_scaled[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_DT800_norm[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_DT800_norm_ad%d", iad+1);
		//	h_total_delayed_energy_DT800_norm[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_DT800_norm[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_DT800_DTnorm[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_DT800_DTnorm_ad%d", iad+1);
			h_total_delayed_energy_DT800_DTnorm[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_scaled[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_scaled_ad%d", iad+1);
		//	h_total_delayed_energy_fine_scaled[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_fine_scaled[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_norm[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_norm_ad%d", iad+1);
		//	h_total_delayed_energy_fine_norm[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_fine_norm[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DTnorm[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DTnorm_ad%d", iad+1);
			h_total_delayed_energy_fine_DTnorm[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DT800_scaled_1[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_scaled_1_ad%d", iad+1);
		//	h_total_delayed_energy_fine_DT800_scaled[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_fine_DT800_scaled_1[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DT800_norm_1[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_norm_1_ad%d", iad+1);
		//	h_total_delayed_energy_fine_DT800_norm[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_fine_DT800_norm_1[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DT800_DTnorm_1[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_DTnorm_1_ad%d", iad+1);
			h_total_delayed_energy_fine_DT800_DTnorm_1[iad]=new TH1F(name,name,23000,0.7,3.);
		}
		
		TH1F* h_total_delayed_energy_fine_DT800_scaled_2[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_scaled_2_ad%d", iad+1);
		//	h_total_delayed_energy_fine_DT800_scaled[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_fine_DT800_scaled_2[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DT800_norm_2[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_norm_2_ad%d", iad+1);
		//	h_total_delayed_energy_fine_DT800_norm[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_fine_DT800_norm_2[iad]=new TH1F(name,name,23000,0.7,3.);
		}

		TH1F* h_total_delayed_energy_fine_DT800_DTnorm_2[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_fine_DT800_DTnorm_2_ad%d", iad+1);
			h_total_delayed_energy_fine_DT800_DTnorm_2[iad]=new TH1F(name,name,23000,0.7,3.);
		}

			TH1F* h_total_delayed_energy_fine_Ep35_scaled[maxAD]; //delayed energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_fine_Ep35_scaled_ad%d", iad+1);
				h_total_delayed_energy_fine_Ep35_scaled[iad]=new TH1F(name,name,2300,0.7,3.);
			}

			TH1F* h_total_delayed_energy_fine_Ep35_norm[maxAD]; //delayed energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_fine_Ep35_norm_ad%d", iad+1);
				h_total_delayed_energy_fine_Ep35_norm[iad]=new TH1F(name,name,2300,0.7,3.);
			}

			TH1F* h_total_delayed_energy_fine_Ep35_DTnorm[maxAD]; //delayed energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_fine_Ep35_DTnorm_ad%d", iad+1);
				h_total_delayed_energy_fine_Ep35_DTnorm[iad]=new TH1F(name,name,2300,0.7,3.);
			}

			TH1F* h_total_delayed_energy_fine_DT800_Ep35_scaled[maxAD]; //delayed energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_fine_DT800_Ep35_scaled_ad%d", iad+1);
				h_total_delayed_energy_fine_DT800_Ep35_scaled[iad]=new TH1F(name,name,2300,0.7,3.);
			}

			TH1F* h_total_delayed_energy_fine_DT800_Ep35_norm[maxAD]; //delayed energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_fine_DT800_Ep35_norm_ad%d", iad+1);
				h_total_delayed_energy_fine_DT800_Ep35_norm[iad]=new TH1F(name,name,2300,0.7,3.);
			}

			TH1F* h_total_delayed_energy_fine_DT800_Ep35_DTnorm[maxAD]; //delayed energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_fine_DT800_Ep35_DTnorm_ad%d", iad+1);
				h_total_delayed_energy_fine_DT800_Ep35_DTnorm[iad]=new TH1F(name,name,2300,0.7,3.);
			}


//Start of section for looking at distance ranges of delayed energy scaled plots
		TH1F* h_total_delayed_energy_scaled_15_22[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_scaled_15_22_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_15_22[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_15_22[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_22_29[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_scaled_22_29_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_22_29[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_22_29[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_29_36[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_scaled_29_36_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_29_36[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_29_36[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_36_43[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_scaled_36_43_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_36_43[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_36_43[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_43_50[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_scaled_43_50_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_43_50[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_43_50[iad]=new TH1F(name,name,230,0.7,3.);
		}
//End of section for looking at distance ranges of delayed energy scaled plots

		TH1F* h_total_prompt_energy_scaled[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_scaled_ad%d", iad+1);
			h_total_prompt_energy_scaled[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_norm[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_norm_ad%d", iad+1);
			h_total_prompt_energy_norm[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_DTnorm[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_DTnorm_ad%d", iad+1);
			h_total_prompt_energy_DTnorm[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_DT800_scaled[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_DT800_scaled_ad%d", iad+1);
			h_total_prompt_energy_DT800_scaled[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_DT800_norm[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_DT800_norm_ad%d", iad+1);
			h_total_prompt_energy_DT800_norm[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_prompt_energy_DT800_DTnorm[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_DT800_DTnorm_ad%d", iad+1);
			h_total_prompt_energy_DT800_DTnorm[iad]=new TH1F(name,name,113,0.7,12.);
		}

			TH1F* h_total_prompt_energy_DT800_3sig_scaled[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_prompt_energy_DT800_3sig_scaled_ad%d", iad+1);
				h_total_prompt_energy_DT800_3sig_scaled[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_total_prompt_energy_DT800_3sig_norm[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_prompt_energy_DT800_3sig_norm_ad%d", iad+1);
				h_total_prompt_energy_DT800_3sig_norm[iad]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_total_prompt_energy_DT800_3sig_DTnorm[maxAD]; //prompt energy histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_prompt_energy_DT800_3sig_DTnorm_ad%d", iad+1);
				h_total_prompt_energy_DT800_3sig_DTnorm[iad]=new TH1F(name,name,113,0.7,12.);
			}

				TH1F* h_total_prompt_energy_DT300_3sig_scaled[maxAD]; //prompt energy histogram
				for(int iad=0; iad<maxAD; ++iad){
					char name[64];
					sprintf(name, "h_total_prompt_energy_DT300_3sig_scaled_ad%d", iad+1);
					h_total_prompt_energy_DT300_3sig_scaled[iad]=new TH1F(name,name,113,0.7,12.);
				}

				TH1F* h_total_prompt_energy_DT500_3sig_scaled[maxAD]; //prompt energy histogram
				for(int iad=0; iad<maxAD; ++iad){
					char name[64];
					sprintf(name, "h_total_prompt_energy_DT500_3sig_scaled_ad%d", iad+1);
					h_total_prompt_energy_DT500_3sig_scaled[iad]=new TH1F(name,name,113,0.7,12.);
				}

				TH1F* h_total_prompt_energy_DT1000_3sig_scaled[maxAD]; //prompt energy histogram
				for(int iad=0; iad<maxAD; ++iad){
					char name[64];
					sprintf(name, "h_total_prompt_energy_DT1000_3sig_scaled_ad%d", iad+1);
					h_total_prompt_energy_DT1000_3sig_scaled[iad]=new TH1F(name,name,113,0.7,12.);
				}

				TH1F* h_total_prompt_energy_DT1500_3sig_scaled[maxAD]; //prompt energy histogram
				for(int iad=0; iad<maxAD; ++iad){
					char name[64];
					sprintf(name, "h_total_prompt_energy_DT1500_3sig_scaled_ad%d", iad+1);
					h_total_prompt_energy_DT1500_3sig_scaled[iad]=new TH1F(name,name,113,0.7,12.);
				}

				TH1F* h_total_prompt_energy_DT2000_3sig_scaled[maxAD]; //prompt energy histogram
				for(int iad=0; iad<maxAD; ++iad){
					char name[64];
					sprintf(name, "h_total_prompt_energy_DT2000_3sig_scaled_ad%d", iad+1);
					h_total_prompt_energy_DT2000_3sig_scaled[iad]=new TH1F(name,name,113,0.7,12.);
				}


	//Making the efficiency hists:
		TH1F* h_p_muon_efficiency[maxAD]; //efficiency histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_p_muon_efficiency_ad%d", iad+1);
			h_p_muon_efficiency[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}
		TH1F* h_d_muon_efficiency[maxAD]; //efficiency histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_d_muon_efficiency_ad%d", iad+1);
			h_d_muon_efficiency[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_p_rate_before[maxAD]; //prompt rate histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_p_rate_before_ad%d", iad+1);
			h_p_rate_before[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}
		TH1F* h_d_rate_before[maxAD]; //delayed rate histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_d_rate_before_ad%d", iad+1);
			h_d_rate_before[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_p_rate[maxAD]; //prompt rate histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_p_rate_ad%d", iad+1);
			h_p_rate[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}
		TH1F* h_d_rate[maxAD]; //delayed rate histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_d_rate_ad%d", iad+1);
			h_d_rate[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_pd_ratio[maxAD]; //ratio of prompt to delayed histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_pd_ratio_ad%d", iad+1);
			h_pd_ratio[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_acc_rate[maxAD]; //efficiency histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_acc_rate_ad%d", iad+1);
			h_acc_rate[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_DAQ[maxAD]; //DAQ histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_DAQ_ad%d", iad+1);
			h_DAQ[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_scale[maxAD]; //scale histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_rate_scale_ad%d", iad+1);
			h_scale[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_normScale[maxAD]; //scale histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_norm_scale_ad%d", iad+1);
			h_normScale[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_diffScales_time[maxAD]; //scale histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_diffScales_time_ad%d", iad+1);
			h_diffScales_time[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

		TH1F* h_diffScales_sigma_time[maxAD]; //scale histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_diffScales_sigma_time_ad%d", iad+1);
			h_diffScales_sigma_time[iad]=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);
		}

			TH1F* h_diffScales[maxAD]; //scale histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_diffScales_ad%d", iad+1);
				h_diffScales[iad]=new TH1F(name,name,80,-0.01,0.01);
			}

			TH1F* h_diffScales_sigma[maxAD]; //scale histogram
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_diffScales_sigma_ad%d", iad+1);
				h_diffScales_sigma[iad]=new TH1F(name,name,100,-5.0,5.0);
			}

		TH2F* h_total_delayed_energy_vs_distance_rate[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_vs_distance_rate_ad%d", iad+1);
		//	h_total_delayed_energy_vs_distance_rate[iad]=new TH2F(name,name,150,1.5,3.,10,0,5.);
			h_total_delayed_energy_vs_distance_rate[iad]=new TH2F(name,name,230,0.7,3.,11,0,5.5);
		}

		TH2F* h_total_delayed_energy_vs_distance_norm[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_vs_distance_norm_ad%d", iad+1);
		//	h_total_delayed_energy_vs_distance_norm[iad]=new TH2F(name,name,150,1.5,3.,10,0,5.);
			h_total_delayed_energy_vs_distance_norm[iad]=new TH2F(name,name,230,0.7,3.,11,0,5.5);
		}

		TH2F* h_total_prompt_energy_vs_distance_rate[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_vs_distance_rate_ad%d", iad+1);
			h_total_prompt_energy_vs_distance_rate[iad]=new TH2F(name,name,113,0.7,12.,11,0,5.5);
		}

		TH2F* h_total_prompt_energy_vs_distance_norm[maxAD]; //prompt energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_vs_distance_norm_ad%d", iad+1);
			h_total_prompt_energy_vs_distance_norm[iad]=new TH2F(name,name,113,0.7,12.,11,0,5.5);
		}

		TH2F* h_total_acc_energy_before[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_acc_energy_before_ad%d",iad+1);
		//	h_total_acc_energy_before[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_total_acc_energy_before[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH2F* h_total_acc_energy_1m[maxAD]; //prompt vs. delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name,"h_total_acc_energy_1m_ad%d",iad+1);
		//	h_total_acc_energy_1m[iad]=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
			h_total_acc_energy_1m[iad]=new TH2F(name,name,175,0.7,12.,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_scaled_dist%d_ad%d",dist, iad+1);
			//	h_total_delayed_energy_scaled_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_total_delayed_energy_scaled_dist[iad][dist]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_total_prompt_energy_scaled_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_prompt_energy_scaled_dist%d_ad%d",dist, iad+1);
			//	h_total_prompt_energy_scaled_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_total_prompt_energy_scaled_dist[iad][dist]=new TH1F(name,name,113,0.7,12.);
			}
		}

		TH1F* h_total_delayed_energy_norm_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_delayed_energy_norm_dist%d_ad%d",dist, iad+1);
			//	h_total_delayed_energy_norm_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_total_delayed_energy_norm_dist[iad][dist]=new TH1F(name,name,230,0.7,3.);
			}
		}

		TH1F* h_total_prompt_energy_norm_dist[maxAD][6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name, "h_total_prompt_energy_norm_dist%d_ad%d",dist, iad+1);
			//	h_total_prompt_energy_norm_dist[iad][dist]=new TH1F(name,name,150,1.5,3.);
				h_total_prompt_energy_norm_dist[iad][dist]=new TH1F(name,name,113,0.7,12.);
			}
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

		TH1F* h_total_delayed_energy_scaled_p35[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_scaled_p35_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_p35[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_p35[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_delayed_energy_scaled_dist_2plus[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_scaled_dist_2plus_ad%d", iad+1);
		//	h_total_delayed_energy_scaled_dist_2plus[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_scaled_dist_2plus[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_prompt_energy_scaled_dist_2plus[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_scaled_dist_2plus_ad%d", iad+1);
		//	h_total_prompt_energy_scaled_dist_2plus[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_prompt_energy_scaled_dist_2plus[iad]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_total_delayed_energy_norm_dist_2plus[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_delayed_energy_norm_dist_2plus_ad%d", iad+1);
		//	h_total_delayed_energy_norm_dist_2plus[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_delayed_energy_norm_dist_2plus[iad]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_total_prompt_energy_norm_dist_2plus[maxAD]; //delayed energy histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_prompt_energy_norm_dist_2plus_ad%d", iad+1);
		//	h_total_prompt_energy_norm_dist_2plus[iad]=new TH1F(name,name,150,1.5,3.);
			h_total_prompt_energy_norm_dist_2plus[iad]=new TH1F(name,name,113,0.7,12.);
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

			TH2F* h_total_acc_distVStime[maxAD];
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_acc_distVStime_ad%d",iad+1);
				h_total_acc_distVStime[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

			TH2F* h_total_acc_distVStime_norm[maxAD];
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_acc_distVStime_norm_ad%d",iad+1);
				h_total_acc_distVStime_norm[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

			TH2F* h_total_acc_distVStime_rate[maxAD];
			for(int iad=0; iad<maxAD; ++iad){
				char name[64];
				sprintf(name,"h_total_acc_distVStime_rate_ad%d",iad+1);
				h_total_acc_distVStime_rate[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
			}

				TH2F* h_total_acc_distVStime_Ep35[maxAD];
				for(int iad=0; iad<maxAD; ++iad){
					char name[64];
					sprintf(name,"h_total_acc_distVStime_Ep35_ad%d",iad+1);
					h_total_acc_distVStime_Ep35[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
				}

				TH2F* h_total_acc_distVStime_Ep35_norm[maxAD];
				for(int iad=0; iad<maxAD; ++iad){
					char name[64];
					sprintf(name,"h_total_acc_distVStime_Ep35_norm_ad%d",iad+1);
					h_total_acc_distVStime_Ep35_norm[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
				}

				TH2F* h_total_acc_distVStime_Ep35_rate[maxAD];
				for(int iad=0; iad<maxAD; ++iad){
					char name[64];
					sprintf(name,"h_total_acc_distVStime_Ep35_rate_ad%d",iad+1);
					h_total_acc_distVStime_Ep35_rate[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);
				}

		TH1D* h_total_acc_DT[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_ad%d", iad+1);
			h_total_acc_DT[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_rate[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_rate_ad%d", iad+1);
			h_total_acc_DT_rate[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_norm[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_norm_ad%d", iad+1);
			h_total_acc_DT_norm[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_DTnorm[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_DTnorm_ad%d", iad+1);
			h_total_acc_DT_DTnorm[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_3sig_rate[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_3sig_rate_ad%d", iad+1);
			h_total_acc_DT_3sig_rate[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_3sig_norm[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_3sig_norm_ad%d", iad+1);
			h_total_acc_DT_3sig_norm[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_3sig_DTnorm[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_3sig_DTnorm_ad%d", iad+1);
			h_total_acc_DT_3sig_DTnorm[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_Ep35[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_Ep35_ad%d", iad+1);
			h_total_acc_DT_Ep35[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_Ep35_norm[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_Ep35_norm_ad%d", iad+1);
			h_total_acc_DT_Ep35_norm[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_Ep35_rate[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_Ep35_rate_ad%d", iad+1);
			h_total_acc_DT_Ep35_rate[iad]=new TH1D(name,name,500,0,10);
		}

		TH1D* h_total_acc_DT_Ep35_DTnorm[maxAD]; //delta time histogram
		for(int iad=0; iad<maxAD; ++iad){
			char name[64];
			sprintf(name, "h_total_acc_DT_Ep35_DTnorm_ad%d", iad+1);
			h_total_acc_DT_Ep35_DTnorm[iad]=new TH1D(name,name,500,0,10);
		}

	int run_order = 0;
	int binCounts = 0;
	double p_DAQ_ad1 = 0;
	double p_DAQ_ad2 = 0;
	double p_DAQ_ad3 = 0;
	double p_DAQ_ad4 = 0;
	double d_DAQ_ad1 = 0;
	double d_DAQ_ad2 = 0;
	double d_DAQ_ad3 = 0;
	double d_DAQ_ad4 = 0;
	double p_live_ad1 = 0;
	double p_veto_ad1 = 0;
	double p_live_ad2 = 0;
	double p_veto_ad2 = 0;
	double p_live_ad3 = 0;
	double p_veto_ad3 = 0;
	double p_live_ad4 = 0;
	double p_veto_ad4 = 0;
	double d_live_ad1 = 0;
	double d_veto_ad1 = 0;
	double d_live_ad2 = 0;
	double d_veto_ad2 = 0;
	double d_live_ad3 = 0;
	double d_veto_ad3 = 0;
	double d_live_ad4 = 0;
	double d_veto_ad4 = 0;
	double p_efficiency[maxAD];
	double d_efficiency[maxAD];
	double p_counts[maxAD];
	double d_counts[maxAD];
	double p_rate[maxAD];
	double d_rate[maxAD];
	float p_energy = 0;
	float d_energy = 0;
	float p_x = 0;
	float p_y = 0;
	float p_z = 0;
	float d_x = 0;
	float d_y = 0;
	float d_z = 0;
	int det_num = 0;
	float distance = 0;
//	int pair_entries = 0;
//	float timing = 400.e-6;
//	float timing = 2000.e-6;
	double timing = (pd_window_microsec+800)*1.e-6;

	int runNum = 0;
	double ibd2m[4];
	double ibd2m_400[4];
	double ibd2m_600[4];
	double ibd2m_800[4];
	double sigma[4];
	for(int iad=0; iad<4; iad++){
		ibd2m[iad] = 0;
		ibd2m_400[iad] = 0;
		ibd2m_600[iad] = 0;
		ibd2m_800[iad] = 0;
		sigma[iad] = 0;
	}

	double TOTAL_DAQ[8];
	double TOTAL_LIVE[8];
	double run_Eff_mult[8];
	double TOTAL_Eff_mult[8];
	double TOTAL_prompt_rate[8];
	double TOTAL_delayed_rate[8];
	double AVERAGE_SCALE[8];

		for(int i=0; i<8; i++){ //Initialized the vaules
			prompt_live[i] = 0;
			prompt_DAQ[i] = 0;
			delayed_live[i] = 0;
			delayed_DAQ[i] = 0;
			prompt_rate[i] = 0;
			delayed_rate[i] = 0;
			num_prompt[i] = 0;
			num_delayed[i] = 0;
			prompt_div[i] = 0;
			delayed_div[i] = 0;
			acc_rate[i] = 0;
			acc_rate_400[i] = 0;
			acc_rate_600[i] = 0;
			acc_rate_800[i] = 0;
			scale[i] = 0;
			scale_400[i] = 0;
			scale_600[i] = 0;
			scale_800[i] = 0;
			normScale[i] = 0;
			normScale_400[i] = 0;
			normScale_600[i] = 0;
			normScale_800[i] = 0;
				DTscale[i] = 0;
				DTscale_Ep35[i] = 0;
			acc_counts[i]=0;
			p_rate_error[i] = 0;
			d_rate_error[i] = 0;
			TOTAL_DAQ[i] = 0;
			TOTAL_LIVE[i] = 0;
			run_Eff_mult[i] = 0;
			TOTAL_Eff_mult[i] = 0;
			TOTAL_prompt_rate[i] = 0;
			TOTAL_delayed_rate[i] = 0;
			AVERAGE_SCALE[i] = 0;
		}

	ofstream resub;
	char resubName[64];
	sprintf(resubName, "./resub_singles_EH%d.sh",hall_num);
	resub.open(resubName);
	resub << "#!/bin/bash" << endl;

	while(1){
		fscanf(runfile,"%d %d",&run_num,&EH);
		if(feof(runfile)) break; //If it's the end of the file, break.
		if(EH != hall_num){
			cout << " WARNING: HALL NUMBERS DO NOT MATCH!!! Run #" << run_num << endl;
			continue;
		}


		char ibdFileName[64];
//		sprintf(ibdFileName, "./IBDs/EH%d/summary_TcLong_%d.root",EH,run_num);
		sprintf(ibdFileName, "./IBDs/EH%d/summary_ktrain_%d_%d.root",EH,pd_window_microsec,run_num);
		TFile *ibdFile = new TFile(ibdFileName);


		//For failed runs:
	/*	if(run_num == 24614 || run_num == 37322 || run_num == 37645 || run_num == 63825 || run_num == 63941){
			run[EH-1] +=1;
			run_order += 1;
			continue;
		}*/
		//For failed runs:
	/*	if(run_num == 65844){
			run[EH-1] +=1;
			run_order += 1;
			continue;
		}*/

		/*if(run_order > 100){
			run[EH-1] +=1;
			run_order+=1;
			continue;

		}*/


/*		if(run_order == 192 || run_order == 205 ||  run_order == 852 || run_order == 853){
			run[EH-1] +=1;
			run_order+=1;
			continue;
		}
*/		

	/*	if(run_order >= 207){
			run[EH-1] +=1;
			run_order+=1;
			continue;

		}

		if(run_order < 207 || run_order >= 414 ){
			run[EH-1] +=1;
			run_order+=1;
			continue;
		}

		if(run_order < 414 || run_order >= 621 ){
			run[EH-1] +=1;
			run_order+=1;
			continue;
		}

		if(run_order < 621 || run_order >= 828 ){
			run[EH-1] +=1;
			run_order+=1;
			continue;
		}
*/
/*		if(run_num < 65800){
			run[EH-1] +=1;
			run_order+=1;
			continue;
		}*/

/*		if(run_order == 94 || run_order == 192 || run_order == 205 ||  run_order == 852 || run_order == 853 || run_order == 70 || run_order == 71 || run_num == 86){
			run[EH-1] +=1;
			run_order+=1;
			continue;
		}
*/

/*		if(run_order < 94){
			run_order += 1;
			ibdFile->Close();
			continue;
		}*/

//cout << "Made it here" << endl;

		char runFileName[64];
//		sprintf(runFileName, "./accResults/EH%d/AccidentalsSummary_%d.root",hall_num, run_num);
//		sprintf(runFileName, "./accResults/EH%d/AccidentalsPlots_TcLong_%d.root",hall_num, run_num);
		sprintf(runFileName, "./accResults/EH%d/AccidentalsPlots_ktrain_%d_%d.root",hall_num, pd_window_microsec, run_num);
//		sprintf(runFileName, "./accResults/EH%d/round1/AccidentalsPlots_%d.root",hall_num, run_num);

//		sprintf(runFileName, "./delayedSingles/EH%d/foundDelayedSingles_%d_%d.root",hall_num, pd_window_microsec, run_num);
cout << "On run# " << run_order+1 << " out of " << nRuns << endl;
		if(FILE *file = fopen(runFileName, "r")){
			fclose(file);
//			continue;
		}
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

			resub << "hep_sub -os SL6 run_singles400.sh -argu " << iline << " -mem 2000" << endl << "sleep 1" << endl;
			continue;
		}

/*		run_order += 1;
		ibdFile->Close();
		continue;*/

		TFile *runFile = new TFile(runFileName);

//cout << "Made it here" << endl;
		
		//TTree tr_sum;
//		TTree *tr_efficiency = (TTree*)runFile->Get("tr_singles_summary");
		TTree *tr_efficiency = (TTree*)runFile->Get("tr_accidentals_time");
		tr_efficiency->SetBranchStatus("*",1);
		tr_efficiency->SetBranchAddress("p_DAQ_ad1",&p_DAQ_ad1);
		tr_efficiency->SetBranchAddress("p_DAQ_ad2",&p_DAQ_ad2);
		tr_efficiency->SetBranchAddress("p_DAQ_ad3",&p_DAQ_ad3);
		tr_efficiency->SetBranchAddress("p_DAQ_ad4",&p_DAQ_ad4);
		tr_efficiency->SetBranchAddress("p_live_ad1",&p_live_ad1);
		tr_efficiency->SetBranchAddress("p_veto_ad1",&p_veto_ad1);
		tr_efficiency->SetBranchAddress("p_live_ad2",&p_live_ad2);
		tr_efficiency->SetBranchAddress("p_veto_ad2",&p_veto_ad2);
		tr_efficiency->SetBranchAddress("d_DAQ_ad1",&d_DAQ_ad1);
		tr_efficiency->SetBranchAddress("d_DAQ_ad2",&d_DAQ_ad2);
		tr_efficiency->SetBranchAddress("d_DAQ_ad3",&d_DAQ_ad3);
		tr_efficiency->SetBranchAddress("d_DAQ_ad4",&d_DAQ_ad4);
		tr_efficiency->SetBranchAddress("d_live_ad1",&d_live_ad1);
		tr_efficiency->SetBranchAddress("d_veto_ad1",&d_veto_ad1);
		tr_efficiency->SetBranchAddress("d_live_ad2",&d_live_ad2);
		tr_efficiency->SetBranchAddress("d_veto_ad2",&d_veto_ad2);
		if(hall_num ==3){
			tr_efficiency->SetBranchAddress("p_live_ad3",&p_live_ad3);
			tr_efficiency->SetBranchAddress("p_veto_ad3",&p_veto_ad3);
			tr_efficiency->SetBranchAddress("p_live_ad4",&p_live_ad4);
			tr_efficiency->SetBranchAddress("p_veto_ad4",&p_veto_ad4);
			tr_efficiency->SetBranchAddress("d_live_ad3",&d_live_ad3);
			tr_efficiency->SetBranchAddress("d_veto_ad3",&d_veto_ad3);
			tr_efficiency->SetBranchAddress("d_live_ad4",&d_live_ad4);
			tr_efficiency->SetBranchAddress("d_veto_ad4",&d_veto_ad4);
		}

		tr_efficiency->GetEntry(0);

		cout << "DAQ time for run" << run_num << " for prompt is: " << p_DAQ_ad2 << "s" << endl;
		cout << "DAQ time for run" << run_num << " for delayed is: " << d_DAQ_ad2 << "s" << endl;

		if(hall_num ==3){
			cout << "DAQ time for run" << run_num << " for prompt is: " << p_DAQ_ad4 << "s" << endl;
			cout << "DAQ time for run" << run_num << " for delayed is: " << d_DAQ_ad4 << "s" << endl;
		}

		//TTree tr_pair;
		/*TTree *tr_data = (TTree*)runFile->Get("tr_pair");
		tr_data->SetBranchStatus("*",1);
		tr_data->SetBranchAddress("p_energy",&p_energy);
		tr_data->SetBranchAddress("d_energy",&d_energy);
		tr_data->SetBranchAddress("p_x",&p_x);
		tr_data->SetBranchAddress("p_y",&p_y);
		tr_data->SetBranchAddress("p_z",&p_z);
		tr_data->SetBranchAddress("d_x",&d_x);
		tr_data->SetBranchAddress("d_y",&d_y);
		tr_data->SetBranchAddress("d_z",&d_z);
		tr_data->SetBranchAddress("det_num",&det_num);
		tr_data->SetBranchAddress("distance",&distance);*/
	
	//cout << "Got the branches..." << endl;

		//pair_entries = tr_data->GetEntries();
	//cout << "Got the entries..... " << pair_entries << endl;
		if(p_DAQ_ad2 == 0 && p_DAQ_ad1 == 0){
			cout << "DAQ time = 0...?!?!?!" << endl;
			run_order += 1;
			runFile->Close();
			continue;
		}

	//cout << "Going into the loop for ADs" << endl;

		for(int iad=0; iad<maxAD; iad++){
			char name[64];

			sprintf(name,"h_distance_before_ad%d",iad+1);
			TH1F *h_ibd_distance = (TH1F*)ibdFile->Get(name);

				int startBin = 0;
				for(int ibin=0; ibin<702; ibin++){
					if(h_ibd_distance->GetBinCenter(ibin)>2){
						startBin = ibin;
						break;
					}
				}
				
		ibd2m[iad] = (h_ibd_distance->Integral(startBin,500));
		

			sprintf(name,"h_p_singles_locations_ad%d",iad+1);
			TH2F *h_run_plocations_before = (TH2F*)runFile->Get(name);
			h_total_plocations_before[iad]->Add(h_run_plocations_before);

			sprintf(name,"h_d_singles_locations_ad%d",iad+1);
			TH2F *h_run_dlocations_before = (TH2F*)runFile->Get(name);
			h_total_dlocations_before[iad]->Add(h_run_dlocations_before);

			sprintf(name,"h_p_singles_energy_raw_ad%d",iad+1);
			TH1F *h_run_prompt_energy_before_raw = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_before_raw[iad]->Add(h_run_prompt_energy_before_raw);

			sprintf(name,"h_d_singles_energy_raw_ad%d",iad+1);
			TH1F *h_run_delayed_energy_before_raw = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_raw[iad]->Add(h_run_delayed_energy_before_raw);

			sprintf(name,"h_p_singles_energy_ad%d",iad+1);
			TH1F *h_run_prompt_energy_before = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_before[iad]->Add(h_run_prompt_energy_before);

			sprintf(name,"h_d_singles_energy_ad%d",iad+1);
			TH1F *h_run_delayed_energy_before = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_before[iad]->Add(h_run_delayed_energy_before);


			for(int iz = 0; iz < NzBins; iz++){
				TH1F *h_run_delayed_energy_before_z = (TH1F*)runFile->Get(Form("h_d_singles_energy_z_ad%d_iz%d", iad+1, iz+1));
				h_total_delayed_energy_before_z[iad][iz]->Add(h_run_delayed_energy_before_z);

				TH1F *h_run_delayed_energy_DT800_z = (TH1F*)runFile->Get(Form("h_d_singles_energy_DT800_z_ad%d_iz%d", iad+1, iz+1));
				h_total_delayed_energy_DT800_z[iad][iz]->Add(h_run_delayed_energy_DT800_z);
			}

			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				TH1F *h_run_delayed_energy_before_r2 = (TH1F*)runFile->Get(Form("h_d_singles_energy_r2_ad%d_ir2%d", iad+1, ir2+1));
				h_total_delayed_energy_before_r2[iad][ir2]->Add(h_run_delayed_energy_before_r2);

				TH1F *h_run_delayed_energy_DT800_r2 = (TH1F*)runFile->Get(Form("h_d_singles_energy_DT800_r2_ad%d_ir2%d", iad+1, ir2+1));
				h_total_delayed_energy_DT800_r2[iad][ir2]->Add(h_run_delayed_energy_DT800_r2);
			}

			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					TH1F *h_run_delayed_energy_before_zVSr2 = (TH1F*)runFile->Get(Form("h_d_singles_energy_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1));
					h_total_delayed_energy_before_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_before_zVSr2);

					TH1F *h_run_delayed_energy_DT800_zVSr2 = (TH1F*)runFile->Get(Form("h_d_singles_energy_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1));
					h_total_delayed_energy_DT800_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_DT800_zVSr2);
				}
			}

			sprintf(name,"h_d_singles_energy_fine_ad%d",iad+1);
			TH1F *h_run_delayed_energy_fine_before = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_fine_before[iad]->Add(h_run_delayed_energy_fine_before);

			sprintf(name,"h_p_singles_energy_DT800_ad%d",iad+1);
			TH1F *h_run_prompt_energy_DT800 = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_DT800[iad]->Add(h_run_prompt_energy_DT800);

			sprintf(name,"h_d_singles_energy_DT800_ad%d",iad+1);
			TH1F *h_run_delayed_energy_DT800 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_DT800[iad]->Add(h_run_delayed_energy_DT800);

			sprintf(name,"h_d_singles_energy_fine_DT800_ad%d",iad+1);
			h_run_delayed_energy_fine_DT800_scaled[iad] = (TH1F*)runFile->Get(name);
	//		h_total_delayed_energy_fine_DT800[iad]->Add(h_run_delayed_energy_fine_DT800);

				sprintf(name,"h_acc_distVStime_ad%d",iad+1);
				TH2F *h_run_acc_distVStime = (TH2F*)runFile->Get(name);
				h_total_acc_distVStime[iad]->Add(h_run_acc_distVStime);

				sprintf(name,"h_acc_distVStime_Ep35_ad%d",iad+1);
				TH2F *h_run_acc_distVStime_Ep35 = (TH2F*)runFile->Get(name);
				h_total_acc_distVStime_Ep35[iad]->Add(h_run_acc_distVStime_Ep35);

				sprintf(name,"h_acc_DT_ad%d",iad+1);
				TH1D *h_run_acc_DT = (TH1D*)runFile->Get(name);
				h_total_acc_DT[iad]->Add(h_run_acc_DT);

				sprintf(name,"h_acc_DT_Ep35_ad%d",iad+1);
				TH1D *h_run_acc_DT_Ep35 = (TH1D*)runFile->Get(name);
				h_total_acc_DT_Ep35[iad]->Add(h_run_acc_DT_Ep35);

				sprintf(name,"h_acc_DT_3sig_ad%d",iad+1);
				TH1D *h_run_acc_DT_3sig = (TH1D*)runFile->Get(name);

			//For largeDist and shortDist plots:
			sprintf(name,"h_p_singles_energy_largeDist_ad%d",iad+1);
			TH1F *h_run_prompt_energy_largeDist = (TH1F*)runFile->Get(name);

			sprintf(name,"h_p_singles_energy_shortDist_ad%d",iad+1);
			TH1F *h_run_prompt_energy_shortDist = (TH1F*)runFile->Get(name);

			sprintf(name,"h_d_singles_energy_largeDist_ad%d",iad+1);
			TH1F *h_run_delayed_energy_largeDist = (TH1F*)runFile->Get(name);

			sprintf(name,"h_d_singles_energy_shortDist_ad%d",iad+1);
			TH1F *h_run_delayed_energy_shortDist = (TH1F*)runFile->Get(name);

				//Calculating rates and efficiencies for Prompt
				num_prompt[2*(EH-1) + iad]=h_run_prompt_energy_before_raw->GetEntries();


				if(num_prompt[2*(EH-1) +iad] == 0) continue;
				prompt_rate[2*(EH-1)+iad]=0;
				if(iad==0){
					prompt_live[2*(EH-1) + iad] = p_live_ad1;
					prompt_DAQ[2*(EH-1) + iad] = p_DAQ_ad1;
					prompt_div[2*(EH-1) + iad] = num_prompt[2*(EH-1) + iad]/p_live_ad1;
					prompt_rate[2*(EH-1) + iad] = correct(prompt_div[2*(EH-1) + iad],timing);
				}
				if(iad==1){
					prompt_live[2*(EH-1) + iad] = p_live_ad2;
					prompt_DAQ[2*(EH-1) + iad] = p_DAQ_ad2;
					prompt_div[2*(EH-1) + iad] = num_prompt[2*(EH-1) + iad]/p_live_ad2;
					prompt_rate[2*(EH-1) + iad] = correct(prompt_div[2*(EH-1) + iad],timing);
				}
				if(iad==2){
					prompt_live[2*(EH-1) + iad] = p_live_ad3;
					prompt_DAQ[2*(EH-1) + iad] = p_DAQ_ad3;
					prompt_div[2*(EH-1) + iad] = num_prompt[2*(EH-1) + iad]/p_live_ad3;
					prompt_rate[2*(EH-1) + iad] = correct(prompt_div[2*(EH-1) + iad],timing);
				}
				if(iad==3){
					prompt_live[2*(EH-1) + iad] = p_live_ad4;
					prompt_DAQ[2*(EH-1) + iad] = p_DAQ_ad4;
					prompt_div[2*(EH-1) + iad] = num_prompt[2*(EH-1) + iad]/p_live_ad4;
					prompt_rate[2*(EH-1) + iad] = correct(prompt_div[2*(EH-1) + iad],timing);
				}

			//	cout << "Calculated the prompt rate" << endl;

			h_p_rate_before[iad]->Fill(run[EH-1],prompt_div[2*(EH-1) + iad]);
			h_p_rate[iad]->Fill(run[EH-1],prompt_rate[2*(EH-1) + iad]);
			h_DAQ[iad]->Fill(run[EH-1],prompt_DAQ[2*(EH-1) + iad]);
			p_rate_error[iad] = sqrt(num_prompt[2*(EH-1)+iad])/prompt_live[2*(EH-1)+iad];
			h_p_muon_efficiency[iad]->Fill(run[EH-1],(prompt_live[2*(EH-1) + iad]/prompt_DAQ[2*(EH-1)+iad])*100);
		//		cout << "Added to the graph..." << endl;
			 //Done doing prompt efficiencies and rates

//			sprintf(name,"h_d_singles_energy_ad%d",iad+1);
//			TH1F *h_run_delayed_energy_before = (TH1F*)runFile->Get(name);
//			h_total_delayed_energy_before[iad]->Add(h_run_delayed_energy_before);

				//Calculating delayed rates and efficiencies
			//	num_delayed[2*(EH-1) + iad]=h_run_delayed_energy_before_raw->GetEntries();
			//	if(num_delayed[2*(EH-1)+iad] == 0) num_delayed[2*(EH-1) + iad]=h_run_delayed_energy_before->GetEntries();

		num_delayed[2*(EH-1) + iad]=h_run_delayed_energy_before->GetEntries(); //For 1.75 MeV cut

				if(num_delayed[2*(EH-1) +iad] == 0) continue;
				delayed_rate[2*(EH-1)+iad]=0;
				if(iad==0){
					delayed_live[2*(EH-1) + iad] = d_live_ad1;
					delayed_DAQ[2*(EH-1) + iad] = d_DAQ_ad1;
					delayed_div[2*(EH-1) + iad] = num_delayed[2*(EH-1) + iad]/d_live_ad1;
//					delayed_rate[2*(EH-1) + iad] = delayed_div[2*(EH-1) + iad]/(exp(-1200.e-6*prompt_rate[2*(EH-1)+iad]));
					delayed_rate[2*(EH-1) + iad] = delayed_div[2*(EH-1) + iad]/(exp(-(timing)*prompt_rate[2*(EH-1)+iad]));
				}
				if(iad==1){
					delayed_live[2*(EH-1) + iad] = d_live_ad2;
					delayed_DAQ[2*(EH-1) + iad] = d_DAQ_ad2;
					delayed_div[2*(EH-1) + iad] = num_delayed[2*(EH-1) + iad]/d_live_ad2;
//					delayed_rate[2*(EH-1) + iad] = delayed_div[2*(EH-1) + iad]/(exp(-1200.e-6*prompt_rate[2*(EH-1)+iad]));
					delayed_rate[2*(EH-1) + iad] = delayed_div[2*(EH-1) + iad]/(exp(-(timing)*prompt_rate[2*(EH-1)+iad]));
				}
				if(iad==2){
					delayed_live[2*(EH-1) + iad] = d_live_ad3;
					delayed_DAQ[2*(EH-1) + iad] = d_DAQ_ad3;
					delayed_div[2*(EH-1) + iad] = num_delayed[2*(EH-1) + iad]/d_live_ad3;
//					delayed_rate[2*(EH-1) + iad] = delayed_div[2*(EH-1) + iad]/(exp(-1200.e-6*prompt_rate[2*(EH-1)+iad]));
					delayed_rate[2*(EH-1) + iad] = delayed_div[2*(EH-1) + iad]/(exp(-(timing)*prompt_rate[2*(EH-1)+iad]));
				}
				if(iad==3){
					delayed_live[2*(EH-1) + iad] = d_live_ad4;
					delayed_DAQ[2*(EH-1) + iad] = d_DAQ_ad4;
					delayed_div[2*(EH-1) + iad] = num_delayed[2*(EH-1) + iad]/d_live_ad4;
//					delayed_rate[2*(EH-1) + iad] = delayed_div[2*(EH-1) + iad]/(exp(-1200.e-6*prompt_rate[2*(EH-1)+iad]));
					delayed_rate[2*(EH-1) + iad] = delayed_div[2*(EH-1) + iad]/(exp(-(timing)*prompt_rate[2*(EH-1)+iad]));
				}
		//		cout << "Calculated the delayed rate" << endl;
			h_d_rate_before[iad]->Fill(run[EH-1],delayed_div[2*(EH-1) + iad]);
			h_d_rate[iad]->Fill(run[EH-1],delayed_rate[2*(EH-1) + iad]);
		//		cout << "Filled the rate..." << endl;
			if(delayed_rate[2*(EH-1)+iad] != 0) h_pd_ratio[iad]->Fill(run[EH-1],prompt_rate[2*(EH-1)+iad]/delayed_rate[2*(EH-1)+iad]);
			d_rate_error[iad] = sqrt(num_delayed[2*(EH-1)+iad])/delayed_live[2*(EH-1)+iad];
		//		cout << "Calculated the error..." << endl;
		//		cout << "Set the error... I think the next line is the issue" << endl;
		//	cout << "Prompt_rate is:\t" << prompt_rate[2*(EH-1)+iad] << endl;
		//	cout << "Delayed_rate is:\t" << delayed_rate[2*(EH-1)+iad] << endl;
		//	cout << "Run order is:\t" << run[EH-1] << endl;
		//		cout << "Was I right? ..." << endl;
		//
		//		cout << "Or maybe the one after...?" << endl;
			h_d_muon_efficiency[iad]->Fill(run[EH-1],(delayed_live[2*(EH-1) + iad]/delayed_DAQ[2*(EH-1)+iad])*100);
	//			cout << "Added to the graph..." << endl;

			h_d_rate[iad]->SetBinError(run[EH-1]+1, d_rate_error[iad]);
			h_p_rate[iad]->SetBinError(run[EH-1]+1, p_rate_error[iad]);
		//	if(delayed_rate[2*(EH-1)+iad] != 0) h_pd_ratio[iad]->SetBinError(run[EH-1]+1, prompt_rate[2*(EH-1)+iad]/d_rate_error[iad]*sqrt(pow(d_rate_error[iad]/delayed_rate[2*(EH-1)+iad],2)+pow(p_rate_error[iad]/prompt_rate[2*(EH-1)+iad],2)));

			//Data stuff!!!!
/*			for(int ientry=0; ientry < pair_entries; ientry++){
				tr_data->GetEntry(ientry);
	//			h_total_delayed_energy_before[iad]->Fill(d_energy);
	//			cout << "Det_num is: " << det_num << endl;

				if(distance*1.e-3 > 0.5) continue;
				h_total_delayed_energy_after[det_num]->Fill(d_energy);
				h_total_prompt_energy_after[det_num]->Fill(p_energy);
				h_total_acc_distance_after[det_num]->Fill(distance*1.e-3);
				h_total_plocations_after[det_num]->Fill((pow(p_x*1.e-3,2.)+pow(p_y*1.e-3,2.)),p_z*1.e-3);
				h_total_dlocations_after[det_num]->Fill((pow(d_x*1.e-3,2.)+pow(d_y*1.e-3,2.)),d_z*1.e-3);
	
			}*/
			acc_rate[2*(EH-1)+iad] = prompt_rate[2*(EH-1)+iad]*delayed_rate[2*(EH-1)+iad]*((pd_window_microsec-1.)*1.e-6)*exp(-prompt_rate[2*(EH-1)+iad]*((pd_window_microsec+800)*1.e-6)); //399
//			acc_rate[2*(EH-1)+iad] = prompt_rate[2*(EH-1)+iad]*delayed_rate[2*(EH-1)+iad]*(timing)*exp(-prompt_rate[2*(EH-1)+iad]*3*timing); //400
		//	cout << "Prompt_rate is:\t" << prompt_rate[2*(EH-1)+iad] << endl;
		//	cout << "Delayed_rate is:\t" << delayed_rate[2*(EH-1)+iad] << endl;
		//	cout << "timing is:\t" << timing << endl;

			cout << "Accidental rate for EH " << EH << " AD " << iad+1 << " is: \t" << acc_rate[2*(EH-1)+iad] << endl;
			h_acc_rate[iad]->Fill(run[EH-1],acc_rate[2*(EH-1)+iad]);

			sprintf(name,"h_acc_distance_ad%d",iad+1);
			TH1F *h_run_distance_before = (TH1F*)runFile->Get(name);
			acc_counts[2*(EH-1)+iad] = h_run_distance_before->GetEntries();
			h_total_acc_distance_before[iad]->Add(h_run_distance_before);

			for(int ibin=0; ibin<702; ibin++){
				if((h_run_distance_before->GetBinCenter(ibin))>2){
					startBin = ibin;
					break;
				}
			}
				
		normScale[2*(EH-1)+iad] = ibd2m[iad]/(h_run_distance_before->Integral(startBin,500));
		sigma[iad] = normScale[2*(EH-1)+iad]*sqrt(1/(ibd2m[iad])+1/(h_run_distance_before->Integral(startBin,500)));


			scale[2*(EH-1)+iad] = acc_rate[2*(EH-1)+iad]*prompt_live[2*(EH-1)+iad]/acc_counts[2*(EH-1)+iad];
			h_scale[iad]->Fill(run[EH-1],scale[2*(EH-1)+iad]);
			h_normScale[iad]->Fill(run[EH-1],normScale[2*(EH-1)+iad]);

//			h_diffScales_time[iad]->Fill(run[EH-1],scale[2*(EH-1)+iad]-normScale[2*(EH-1)+iad]);
			h_diffScales[iad]->Fill(scale[2*(EH-1)+iad]-normScale[2*(EH-1)+iad]);
			if(sigma[iad] != 0){
				h_diffScales_sigma_time[iad]->Fill(run[EH-1],(scale[2*(EH-1)+iad]-normScale[2*(EH-1)+iad])/sigma[iad]);
				h_diffScales_sigma[iad]->Fill((scale[2*(EH-1)+iad]-normScale[2*(EH-1)+iad])/sigma[iad]);
			}

			AVERAGE_SCALE[2*(EH-1)+iad] = (AVERAGE_SCALE[2*(EH-1)+iad]*TOTAL_DAQ[2*(EH-1)+iad] + scale[2*(EH-1)+iad])/(TOTAL_DAQ[2*(EH-1)+iad]+prompt_DAQ[2*(EH-1)+iad]);
			run_Eff_mult[2*(EH-1)+iad] = exp(-(timing)*prompt_rate[2*(EH-1)+iad]);
			TOTAL_Eff_mult[2*(EH-1)+iad] = (TOTAL_Eff_mult[2*(EH-1)+iad]*TOTAL_LIVE[2*(EH-1)+iad] + run_Eff_mult[2*(EH-1)+iad]*prompt_live[2*(EH-1)+iad])/(TOTAL_LIVE[2*(EH-1)+iad]+prompt_live[2*(EH-1)+iad]);
			TOTAL_prompt_rate[2*(EH-1)+iad] = (TOTAL_prompt_rate[2*(EH-1)+iad]*TOTAL_LIVE[2*(EH-1)+iad] + prompt_rate[2*(EH-1)+iad]*prompt_live[2*(EH-1)+iad])/(TOTAL_LIVE[2*(EH-1)+iad]+prompt_live[2*(EH-1)+iad]);
			TOTAL_delayed_rate[2*(EH-1)+iad] = (TOTAL_delayed_rate[2*(EH-1)+iad]*TOTAL_LIVE[2*(EH-1)+iad] + delayed_rate[2*(EH-1)+iad]*delayed_live[2*(EH-1)+iad])/(TOTAL_LIVE[2*(EH-1)+iad]+delayed_live[2*(EH-1)+iad]);

			TOTAL_DAQ[2*(EH-1)+iad] += prompt_DAQ[2*(EH-1)+iad];
			TOTAL_LIVE[2*(EH-1)+iad] += prompt_live[2*(EH-1)+iad];

			//Go bin by bin to scale down the accidentals
			for(int ibin=0; ibin<702; ibin++){ // for distance hists
				binCounts = 0;
				h_run_acc_distance_scaled[iad]->SetBinContent(ibin, 0);
				binCounts = h_run_distance_before->GetBinContent(ibin);
				h_run_acc_distance_scaled[iad]->SetBinContent(ibin, binCounts);
				h_run_acc_distance_scaled[iad]->SetBinError(ibin, sqrt(binCounts));
			}

		//	for(int ibin=0; ibin <152; ibin++){ //for delayed energy hists
			for(int ibin=0; ibin <232; ibin++){ //for delayed energy hists
				binCounts = 0;
				h_run_delayed_energy_scaled[iad]->SetBinContent(ibin, 0);
				binCounts = h_run_delayed_energy_before->GetBinContent(ibin);
				h_run_delayed_energy_scaled[iad]->SetBinContent(ibin, binCounts);
				h_run_delayed_energy_scaled[iad]->SetBinError(ibin,sqrt(binCounts));

				binCounts = 0;
				binCounts = h_run_delayed_energy_DT800->GetBinContent(ibin);
				h_run_delayed_energy_DT800_scaled[iad]->SetBinContent(ibin, binCounts);
				h_run_delayed_energy_DT800_scaled[iad]->SetBinError(ibin,sqrt(binCounts));

				binCounts = 0;
				binCounts = h_run_delayed_energy_largeDist->GetBinContent(ibin);
				h_run_delayed_energy_largeDist->SetBinError(ibin,sqrt(binCounts));

				binCounts = 0;
				binCounts = h_run_delayed_energy_shortDist->GetBinContent(ibin);
				h_run_delayed_energy_shortDist->SetBinError(ibin,sqrt(binCounts));
			}
			for(int ibin=0; ibin < 115; ibin++){ //for prompt energy hists
				binCounts = 0;
				h_run_prompt_energy_scaled[iad]->SetBinContent(ibin, 0);
				binCounts = h_run_prompt_energy_before->GetBinContent(ibin);
				h_run_prompt_energy_scaled[iad]->SetBinContent(ibin, binCounts);
				h_run_prompt_energy_scaled[iad]->SetBinError(ibin,sqrt(binCounts));

				binCounts = 0;
				binCounts = h_run_prompt_energy_DT800->GetBinContent(ibin);
				h_run_prompt_energy_DT800_scaled[iad]->SetBinContent(ibin, binCounts);
				h_run_prompt_energy_DT800_scaled[iad]->SetBinError(ibin,sqrt(binCounts));

				binCounts = 0;
				binCounts = h_run_prompt_energy_largeDist->GetBinContent(ibin);
				h_run_prompt_energy_largeDist->SetBinError(ibin,sqrt(binCounts));

				binCounts = 0;
				binCounts = h_run_prompt_energy_shortDist->GetBinContent(ibin);
				h_run_prompt_energy_shortDist->SetBinError(ibin,sqrt(binCounts));

			}

			for(int ibin = 0; ibin< 2302; ibin++){
				binCounts = 0;
				h_run_delayed_energy_fine_scaled[iad]->SetBinContent(ibin, 0);
				binCounts = h_run_delayed_energy_fine_before->GetBinContent(ibin);
				h_run_delayed_energy_fine_scaled[iad]->SetBinContent(ibin, binCounts);
				h_run_delayed_energy_fine_scaled[iad]->SetBinError(ibin,sqrt(binCounts));

	//			binCounts = 0;
	//			binCounts = h_run_delayed_energy_fine_DT800->GetBinContent(ibin);
	//			h_run_delayed_energy_fine_DT800_scaled[iad]->SetBinContent(ibin, binCounts);
	//			h_run_delayed_energy_fine_DT800_scaled[iad]->SetBinError(ibin,sqrt(binCounts));
			}

	/*		h_run_acc_distance_scaled[iad]->Scale(scale[2*(EH-1)+iad]);
			h_run_delayed_energy_scaled[iad]->Scale(scale[2*(EH-1)+iad]);
			h_run_prompt_energy_scaled[iad]->Scale(scale[2*(EH-1)+iad]);
*/

			sprintf(name,"h_acc_energy_before_ad%d",iad+1);
			TH2F *h_run_acc_energy_before = (TH2F*)runFile->Get(name);
			h_total_acc_energy_before[iad]->Add(h_run_acc_energy_before,scale[2*(EH-1)+iad]);

			sprintf(name,"h_acc_energy_1m_ad%d",iad+1);
			TH2F *h_run_acc_energy_1m = (TH2F*)runFile->Get(name);
			h_total_acc_energy_1m[iad]->Add(h_run_acc_energy_1m,scale[2*(EH-1)+iad]);

			sprintf(name,"h_acc_distance_Ep35_ad%d",iad+1);
			TH1F *h_run_acc_distance_Ep35_scaled = (TH1F*)runFile->Get(name);


			sprintf(name,"h_p_singles_energy_DT800_3sig_ad%d",iad+1);
			TH1F *h_run_prompt_energy_DT800_3sig = (TH1F*)runFile->Get(name);

				sprintf(name,"h_p_singles_energy_DT300_3sig_ad%d",iad+1);
				TH1F *h_run_prompt_energy_DT300_3sig = (TH1F*)runFile->Get(name);

				sprintf(name,"h_p_singles_energy_DT500_3sig_ad%d",iad+1);
				TH1F *h_run_prompt_energy_DT500_3sig = (TH1F*)runFile->Get(name);

				sprintf(name,"h_p_singles_energy_DT1000_3sig_ad%d",iad+1);
				TH1F *h_run_prompt_energy_DT1000_3sig = (TH1F*)runFile->Get(name);

				sprintf(name,"h_p_singles_energy_DT1500_3sig_ad%d",iad+1);
				TH1F *h_run_prompt_energy_DT1500_3sig = (TH1F*)runFile->Get(name);

				sprintf(name,"h_p_singles_energy_DT2000_3sig_ad%d",iad+1);
				TH1F *h_run_prompt_energy_DT2000_3sig = (TH1F*)runFile->Get(name);

			sprintf(name,"h_acc_distance_3sig_ad%d",iad+1);
			TH1F *h_run_acc_distance_3sig = (TH1F*)runFile->Get(name);

			h_total_acc_distance_scaled[iad]->Add(h_run_acc_distance_scaled[iad],scale[2*(EH-1)+iad]);
				h_total_acc_distance_3sig_scaled[iad]->Add(h_run_acc_distance_3sig,scale[2*(EH-1)+iad]);
			h_total_acc_distance_Ep35_scaled[iad]->Add(h_run_acc_distance_Ep35_scaled,scale[2*(EH-1)+iad]);
			h_total_delayed_energy_scaled[iad]->Add(h_run_delayed_energy_scaled[iad],scale[2*(EH-1)+iad]);
				h_total_delayed_energy_fine_scaled[iad]->Add(h_run_delayed_energy_fine_scaled[iad],scale[2*(EH-1)+iad]);
			h_total_prompt_energy_scaled[iad]->Add(h_run_prompt_energy_scaled[iad],scale[2*(EH-1)+iad]);
			h_total_delayed_energy_DT800_scaled[iad]->Add(h_run_delayed_energy_DT800_scaled[iad],scale[2*(EH-1)+iad]);
				if(run_order < 550) h_total_delayed_energy_fine_DT800_scaled_1[iad]->Add(h_run_delayed_energy_fine_DT800_scaled[iad],scale[2*(EH-1)+iad]);
				else h_total_delayed_energy_fine_DT800_scaled_2[iad]->Add(h_run_delayed_energy_fine_DT800_scaled[iad],scale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT800_scaled[iad]->Add(h_run_prompt_energy_DT800_scaled[iad],scale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT800_3sig_scaled[iad]->Add(h_run_prompt_energy_DT800_3sig,scale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT300_3sig_scaled[iad]->Add(h_run_prompt_energy_DT300_3sig,scale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT500_3sig_scaled[iad]->Add(h_run_prompt_energy_DT500_3sig,scale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT1000_3sig_scaled[iad]->Add(h_run_prompt_energy_DT1000_3sig,scale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT1500_3sig_scaled[iad]->Add(h_run_prompt_energy_DT1500_3sig,scale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT2000_3sig_scaled[iad]->Add(h_run_prompt_energy_DT2000_3sig,scale[2*(EH-1)+iad]);


//	cout << "Before Ep35" << endl;
				sprintf(name,"h_delayed_energy_scaled_Ep35_ad%d",iad+1);
				TH1F *h_run_delayed_energy_fine_Ep35 = (TH1F*)runFile->Get(name);
					h_run_delayed_energy_fine_Ep35->Rebin(10);
					h_total_delayed_energy_fine_Ep35_scaled[iad]->Add(h_run_delayed_energy_fine_Ep35,scale[2*(EH-1)+iad]);
//	cout << "After Ep35, Before DT800_Ep35" << endl;
				sprintf(name,"h_delayed_energy_scaled_DT800_Ep35_ad%d",iad+1);
				TH1F *h_run_delayed_energy_fine_DT800_Ep35 = (TH1F*)runFile->Get(name);
					h_run_delayed_energy_fine_DT800_Ep35->Rebin(10);
					h_total_delayed_energy_fine_DT800_Ep35_scaled[iad]->Add(h_run_delayed_energy_fine_DT800_Ep35,scale[2*(EH-1)+iad]);
//	cout << "After Ep35 and DT800_Ep35" << endl;

			h_total_acc_distance_norm[iad]->Add(h_run_acc_distance_scaled[iad],normScale[2*(EH-1)+iad]);
				h_total_acc_distance_3sig_norm[iad]->Add(h_run_acc_distance_3sig,normScale[2*(EH-1)+iad]);
			h_total_acc_distance_Ep35_norm[iad]->Add(h_run_acc_distance_Ep35_scaled,normScale[2*(EH-1)+iad]);
			h_total_delayed_energy_norm[iad]->Add(h_run_delayed_energy_scaled[iad],normScale[2*(EH-1)+iad]);
				h_total_delayed_energy_fine_norm[iad]->Add(h_run_delayed_energy_fine_scaled[iad],normScale[2*(EH-1)+iad]);
			h_total_prompt_energy_norm[iad]->Add(h_run_prompt_energy_scaled[iad],normScale[2*(EH-1)+iad]);
			h_total_delayed_energy_DT800_norm[iad]->Add(h_run_delayed_energy_DT800_scaled[iad],normScale[2*(EH-1)+iad]);
				if(run_order < 550) h_total_delayed_energy_fine_DT800_norm_1[iad]->Add(h_run_delayed_energy_fine_DT800_scaled[iad],normScale[2*(EH-1)+iad]);
				else h_total_delayed_energy_fine_DT800_norm_2[iad]->Add(h_run_delayed_energy_fine_DT800_scaled[iad],normScale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT800_norm[iad]->Add(h_run_prompt_energy_DT800_scaled[iad],normScale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT800_3sig_norm[iad]->Add(h_run_prompt_energy_DT800_3sig,normScale[2*(EH-1)+iad]);
					h_total_delayed_energy_fine_Ep35_norm[iad]->Add(h_run_delayed_energy_fine_Ep35,normScale[2*(EH-1)+iad]);
					h_total_delayed_energy_fine_DT800_Ep35_norm[iad]->Add(h_run_delayed_energy_fine_DT800_Ep35,normScale[2*(EH-1)+iad]);

		/*	h_total_prompt_energy_largeDist[iad]->Add(h_run_prompt_energy_largeDist,scale[2*(EH-1)+iad]);

			h_total_prompt_energy_shortDist[iad]->Add(h_run_prompt_energy_shortDist,scale[2*(EH-1)+iad]);

			h_total_delayed_energy_largeDist[iad]->Add(h_run_delayed_energy_largeDist,scale[2*(EH-1)+iad]);

			h_total_delayed_energy_shortDist[iad]->Add(h_run_delayed_energy_shortDist,scale[2*(EH-1)+iad]);*/

/*			sprintf(name,"h_acc_distVStime_ad%d",iad+1);
			TH2F *h_run_acc_distVStime_norm = (TH2F*)runFile->Get(name);
			h_run_acc_distVStime_norm->Scale(normScale[2*(EH-1)+iad]);*/

			h_total_acc_distVStime_norm[iad]->Add(h_run_acc_distVStime,normScale[2*(EH-1)+iad]);
			h_total_acc_distVStime_Ep35_norm[iad]->Add(h_run_acc_distVStime_Ep35,normScale[2*(EH-1)+iad]);
			h_total_acc_DT_norm[iad]->Add(h_run_acc_DT,normScale[2*(EH-1)+iad]);
			h_total_acc_DT_3sig_norm[iad]->Add(h_run_acc_DT_3sig,normScale[2*(EH-1)+iad]);
			h_total_acc_DT_Ep35_norm[iad]->Add(h_run_acc_DT_Ep35,normScale[2*(EH-1)+iad]);

			sprintf(name, "h_delayed_energy_scaled_p35_ad%d", iad+1);
			TH1F *h_run_delayed_p35 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_p35[iad]->Add(h_run_delayed_p35,scale[2*(EH-1)+iad]);

/*			for(int dist = 0; dist < 6; dist++){
				sprintf(name, "h_delayed_energy_scaled_dist%d_ad%d",dist, iad+1);
				TH1F *h_run_delayed_dist = (TH1F*)runFile->Get(name);
				h_total_delayed_energy_dist[iad][dist]->Add(h_run_delayed_dist);
				if(dist >= 2) h_total_delayed_energy_dist_2plus[iad]->Add(h_run_delayed_dist);
				h_total_delayed_energy_scaled_dist[iad][dist]->Add(h_run_delayed_dist,scale[2*(EH-1)+iad]);
				if(dist >= 2) h_total_delayed_energy_scaled_dist_2plus[iad]->Add(h_run_delayed_dist,scale[2*(EH-1)+iad]);
				h_total_delayed_energy_norm_dist[iad][dist]->Add(h_run_delayed_dist,normScale[2*(EH-1)+iad]);
				if(dist >= 2) h_total_delayed_energy_norm_dist_2plus[iad]->Add(h_run_delayed_dist,normScale[2*(EH-1)+iad]);

				sprintf(name, "h_prompt_energy_scaled_dist%d_ad%d",dist, iad+1);
				TH1F *h_run_prompt_dist = (TH1F*)runFile->Get(name);
				h_total_prompt_energy_dist[iad][dist]->Add(h_run_prompt_dist);
				if(dist >= 2) h_total_prompt_energy_dist_2plus[iad]->Add(h_run_prompt_dist);
				h_total_prompt_energy_scaled_dist[iad][dist]->Add(h_run_prompt_dist,scale[2*(EH-1)+iad]);
				if(dist >= 2) h_total_prompt_energy_scaled_dist_2plus[iad]->Add(h_run_prompt_dist,scale[2*(EH-1)+iad]);
				h_total_prompt_energy_norm_dist[iad][dist]->Add(h_run_prompt_dist,normScale[2*(EH-1)+iad]);
				if(dist >= 2) h_total_prompt_energy_norm_dist_2plus[iad]->Add(h_run_prompt_dist,normScale[2*(EH-1)+iad]);
			}*/

	/*		sprintf(name, "h_delayed_energy_scaled_15_22_ad%d", iad+1);
			TH1F *h_run_delayed_15_22 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_15_22[iad]->Add(h_run_delayed_15_22,scale[2*(EH-1)+iad]);

			sprintf(name, "h_delayed_energy_scaled_22_29_ad%d", iad+1);
			TH1F *h_run_delayed_22_29 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_22_29[iad]->Add(h_run_delayed_22_29,scale[2*(EH-1)+iad]);

			sprintf(name, "h_delayed_energy_scaled_29_36_ad%d", iad+1);
			TH1F *h_run_delayed_29_36 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_29_36[iad]->Add(h_run_delayed_29_36,scale[2*(EH-1)+iad]);

			sprintf(name, "h_delayed_energy_scaled_36_43_ad%d", iad+1);
			TH1F *h_run_delayed_36_43 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_36_43[iad]->Add(h_run_delayed_36_43,scale[2*(EH-1)+iad]);

			sprintf(name, "h_delayed_energy_scaled_43_50_ad%d", iad+1);
			TH1F *h_run_delayed_43_50 = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_scaled_43_50[iad]->Add(h_run_delayed_43_50,scale[2*(EH-1)+iad]);

			sprintf(name, "h_delayed_energy_vs_distance_ad%d", iad+1);
			TH1F *h_run_delayed_vs_dist = (TH1F*)runFile->Get(name);
			h_total_delayed_energy_vs_distance_rate[iad]->Add(h_run_delayed_vs_dist,scale[2*(EH-1)+iad]);
			h_total_delayed_energy_vs_distance_norm[iad]->Add(h_run_delayed_vs_dist);

			sprintf(name, "h_prompt_energy_vs_distance_ad%d", iad+1);
			TH1F *h_run_prompt_vs_dist = (TH1F*)runFile->Get(name);
			h_total_prompt_energy_vs_distance_rate[iad]->Add(h_run_prompt_vs_dist,scale[2*(EH-1)+iad]);
			h_total_prompt_energy_vs_distance_norm[iad]->Add(h_run_prompt_vs_dist);*/

/*			sprintf(name,"h_acc_distVStime_ad%d",iad+1);
			TH2F *h_run_acc_distVStime_rate = (TH2F*)runFile->Get(name);
			h_run_acc_distVStime_norm->Scale(scale[2*(EH-1)+iad]);*/
			h_total_acc_distVStime_rate[iad]->Add(h_run_acc_distVStime,scale[2*(EH-1)+iad]);
			h_total_acc_distVStime_Ep35_rate[iad]->Add(h_run_acc_distVStime_Ep35,scale[2*(EH-1)+iad]);
			h_total_acc_DT_rate[iad]->Add(h_run_acc_DT,scale[2*(EH-1)+iad]);
			h_total_acc_DT_3sig_rate[iad]->Add(h_run_acc_DT_3sig,scale[2*(EH-1)+iad]);
			h_total_acc_DT_Ep35_rate[iad]->Add(h_run_acc_DT_Ep35,scale[2*(EH-1)+iad]);


			//DT normalized DT plots
			sprintf(name,"h_ibd_DT_ad%d",iad+1);
			TH1D *h_ibd_DT = (TH1D*)ibdFile->Get(name);

			sprintf(name,"h_ibd_DT_Ep35_ad%d",iad+1);
			TH1D *h_ibd_DT_Ep35 = (TH1D*)ibdFile->Get(name);



				int startBinDT = 151;
				
/*				for(int ibin=0; ibin<502; ibin++){
					if(h_ibd_DT_Ep35->GetBinCenter(ibin)>3){
						startBinDT = ibin;
						cout << startBinDT << endl;
						break;
					}
				}*/


/*			if(EH==3 && run_num==65844) DTscale[2*(EH-1)+iad] = scale[2*(EH-1)+iad];	
			else DTscale[2*(EH-1)+iad] = (h_ibd_DT->Integral(startBinDT,500))/(h_run_acc_DT->Integral(startBinDT,500));*/

			DTscale[2*(EH-1)+iad] = (h_ibd_DT->Integral(startBinDT,500))/(h_run_acc_DT->Integral(startBinDT,500));

			h_diffScales_time[iad]->Fill(run[EH-1],100*(DTscale[2*(EH-1)+iad]-scale[2*(EH-1)+iad])/DTscale[2*(EH-1)+iad]);

			h_total_acc_DT_DTnorm[iad]->Add(h_run_acc_DT,DTscale[2*(EH-1)+iad]);
			h_total_acc_DT_3sig_DTnorm[iad]->Add(h_run_acc_DT_3sig,DTscale[2*(EH-1)+iad]);
			
/*			if(EH==3 && run_num==65844) DTscale_Ep35[2*(EH-1)+iad] = scale[2*(EH-1)+iad];	
			else DTscale_Ep35[2*(EH-1)+iad] = (h_ibd_DT_Ep35->Integral(startBinDT,500))/(h_run_acc_DT_Ep35->Integral(startBinDT,500));*/
			DTscale_Ep35[2*(EH-1)+iad] = (h_ibd_DT_Ep35->Integral(startBinDT,500))/(h_run_acc_DT_Ep35->Integral(startBinDT,500));
			
			h_total_acc_DT_Ep35_DTnorm[iad]->Add(h_run_acc_DT_Ep35,DTscale_Ep35[2*(EH-1)+iad]);
			h_total_acc_DT_Ep35_DTnorm[iad]->Add(h_run_acc_DT_Ep35,DTscale[2*(EH-1)+iad]);

			h_total_prompt_energy_DTnorm[iad]->Add(h_run_prompt_energy_scaled[iad],DTscale[2*(EH-1)+iad]);
			h_total_delayed_energy_DTnorm[iad]->Add(h_run_delayed_energy_scaled[iad],DTscale[2*(EH-1)+iad]);
				h_total_delayed_energy_fine_DTnorm[iad]->Add(h_run_delayed_energy_fine_scaled[iad],DTscale[2*(EH-1)+iad]);

				h_total_acc_distance_3sig_DTnorm[iad]->Add(h_run_acc_distance_3sig,DTscale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT800_DTnorm[iad]->Add(h_run_prompt_energy_DT800_scaled[iad],DTscale[2*(EH-1)+iad]);
			h_total_prompt_energy_DT800_3sig_DTnorm[iad]->Add(h_run_prompt_energy_DT800_3sig,DTscale[2*(EH-1)+iad]);
			h_total_delayed_energy_DT800_DTnorm[iad]->Add(h_run_delayed_energy_DT800_scaled[iad],DTscale[2*(EH-1)+iad]);
				if(run_order < 550) h_total_delayed_energy_fine_DT800_DTnorm_1[iad]->Add(h_run_delayed_energy_fine_DT800_scaled[iad],DTscale[2*(EH-1)+iad]);
				else h_total_delayed_energy_fine_DT800_DTnorm_2[iad]->Add(h_run_delayed_energy_fine_DT800_scaled[iad],DTscale[2*(EH-1)+iad]);
					h_total_delayed_energy_fine_Ep35_DTnorm[iad]->Add(h_run_delayed_energy_fine_Ep35,DTscale[2*(EH-1)+iad]);
					h_total_delayed_energy_fine_DT800_Ep35_DTnorm[iad]->Add(h_run_delayed_energy_fine_DT800_Ep35,DTscale[2*(EH-1)+iad]);

			for(int iz = 0; iz < NzBins; iz++){
				TH1F *h_run_delayed_energy_before_z = (TH1F*)runFile->Get(Form("h_d_singles_energy_z_ad%d_iz%d", iad+1, iz+1));
				h_total_delayed_energy_scaled_z[iad][iz]->Add(h_run_delayed_energy_before_z,scale[2*(EH-1)+iad]);
				h_total_delayed_energy_norm_z[iad][iz]->Add(h_run_delayed_energy_before_z,normScale[2*(EH-1)+iad]);
				h_total_delayed_energy_DTnorm_z[iad][iz]->Add(h_run_delayed_energy_before_z,DTscale[2*(EH-1)+iad]);

				TH1F *h_run_delayed_energy_DT800_z = (TH1F*)runFile->Get(Form("h_d_singles_energy_DT800_z_ad%d_iz%d", iad+1, iz+1));
				h_total_delayed_energy_scaled_DT800_z[iad][iz]->Add(h_run_delayed_energy_DT800_z,scale[2*(EH-1)+iad]);
				h_total_delayed_energy_norm_DT800_z[iad][iz]->Add(h_run_delayed_energy_DT800_z,normScale[2*(EH-1)+iad]);
				h_total_delayed_energy_DTnorm_DT800_z[iad][iz]->Add(h_run_delayed_energy_DT800_z,DTscale[2*(EH-1)+iad]);
			}

			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				TH1F *h_run_delayed_energy_before_r2 = (TH1F*)runFile->Get(Form("h_d_singles_energy_r2_ad%d_ir2%d", iad+1, ir2+1));
				h_total_delayed_energy_scaled_r2[iad][ir2]->Add(h_run_delayed_energy_before_r2,scale[2*(EH-1)+iad]);
				h_total_delayed_energy_norm_r2[iad][ir2]->Add(h_run_delayed_energy_before_r2,normScale[2*(EH-1)+iad]);
				h_total_delayed_energy_DTnorm_r2[iad][ir2]->Add(h_run_delayed_energy_before_r2,DTscale[2*(EH-1)+iad]);

				TH1F *h_run_delayed_energy_DT800_r2 = (TH1F*)runFile->Get(Form("h_d_singles_energy_DT800_r2_ad%d_ir2%d", iad+1, ir2+1));
				h_total_delayed_energy_scaled_DT800_r2[iad][ir2]->Add(h_run_delayed_energy_DT800_r2,scale[2*(EH-1)+iad]);
				h_total_delayed_energy_norm_DT800_r2[iad][ir2]->Add(h_run_delayed_energy_DT800_r2,normScale[2*(EH-1)+iad]);
				h_total_delayed_energy_DTnorm_DT800_r2[iad][ir2]->Add(h_run_delayed_energy_DT800_r2,DTscale[2*(EH-1)+iad]);
			}

			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					TH1F *h_run_delayed_energy_before_zVSr2 = (TH1F*)runFile->Get(Form("h_d_singles_energy_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1));
					h_total_delayed_energy_scaled_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_before_zVSr2,scale[2*(EH-1)+iad]);
					h_total_delayed_energy_norm_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_before_zVSr2,normScale[2*(EH-1)+iad]);
					h_total_delayed_energy_DTnorm_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_before_zVSr2,DTscale[2*(EH-1)+iad]);

					TH1F *h_run_delayed_energy_DT800_zVSr2 = (TH1F*)runFile->Get(Form("h_d_singles_energy_DT800_zVSr2_ad%d_ir2%d_iz%d", iad+1, ir2+1, iz+1));
					h_total_delayed_energy_scaled_DT800_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_DT800_zVSr2,scale[2*(EH-1)+iad]);
					h_total_delayed_energy_norm_DT800_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_DT800_zVSr2,normScale[2*(EH-1)+iad]);
					h_total_delayed_energy_DTnorm_DT800_zVSr2[iad][ir2][iz]->Add(h_run_delayed_energy_DT800_zVSr2,DTscale[2*(EH-1)+iad]);
				}
			}

			sprintf(name,"h_acc_promptVStime_DT800_ad%d",iad+1);
			TH1F *h_run_promptVStime_DT800_DTnorm = (TH1F*)runFile->Get(name);
			h_total_acc_promptVStime_DT800_DTnorm[iad]->Add(h_run_promptVStime_DT800_DTnorm,DTscale[2*(EH-1)+iad]);

			sprintf(name,"h_acc_promptVStime_ad%d",iad+1);
			TH1F *h_run_promptVStime_DTnorm = (TH1F*)runFile->Get(name);
			h_total_acc_promptVStime_DTnorm[iad]->Add(h_run_promptVStime_DTnorm,DTscale[2*(EH-1)+iad]);

			if(EH==3 && run_num==65844) continue;
			//400 us distance plots:
			sprintf(name,"h_distance_before_400_ad%d",iad+1);
			TH1F *h_ibd_distance_400 = (TH1F*)ibdFile->Get(name);
			ibd2m_400[iad] = (h_ibd_distance_400->Integral(startBin,500));
			normScale_400[2*(EH-1)+iad] = ibd2m_400[iad]/(h_run_distance_before->Integral(startBin,500));
			h_total_acc_distance_400_norm[iad]->Add(h_run_acc_distance_scaled[iad],normScale_400[2*(EH-1)+iad]);

			acc_rate_400[2*(EH-1)+iad] = prompt_rate[2*(EH-1)+iad]*delayed_rate[2*(EH-1)+iad]*((400-1.)*1.e-6)*exp(-prompt_rate[2*(EH-1)+iad]*((pd_window_microsec+800)*1.e-6));
			scale_400[2*(EH-1)+iad] = acc_rate_400[2*(EH-1)+iad]*prompt_live[2*(EH-1)+iad]/acc_counts[2*(EH-1)+iad];
			h_total_acc_distance_400_scaled[iad]->Add(h_run_acc_distance_scaled[iad],scale_400[2*(EH-1)+iad]);

			//600 us distance plots:
			sprintf(name,"h_distance_before_600_ad%d",iad+1);
			TH1F *h_ibd_distance_600 = (TH1F*)ibdFile->Get(name);
			ibd2m_600[iad] = (h_ibd_distance_600->Integral(startBin,500));
			normScale_600[2*(EH-1)+iad] = ibd2m_600[iad]/(h_run_distance_before->Integral(startBin,500));
			h_total_acc_distance_600_norm[iad]->Add(h_run_acc_distance_scaled[iad],normScale_600[2*(EH-1)+iad]);

			acc_rate_600[2*(EH-1)+iad] = prompt_rate[2*(EH-1)+iad]*delayed_rate[2*(EH-1)+iad]*((600-1.)*1.e-6)*exp(-prompt_rate[2*(EH-1)+iad]*((pd_window_microsec+800)*1.e-6));
			scale_600[2*(EH-1)+iad] = acc_rate_600[2*(EH-1)+iad]*prompt_live[2*(EH-1)+iad]/acc_counts[2*(EH-1)+iad];
			h_total_acc_distance_600_scaled[iad]->Add(h_run_acc_distance_scaled[iad],scale_600[2*(EH-1)+iad]);

			//800 us distance plots:
			sprintf(name,"h_distance_before_800_ad%d",iad+1);
			TH1F *h_ibd_distance_800 = (TH1F*)ibdFile->Get(name);
			ibd2m_800[iad] = (h_ibd_distance_800->Integral(startBin,500));
			normScale_800[2*(EH-1)+iad] = ibd2m_800[iad]/(h_run_distance_before->Integral(startBin,500));
			h_total_acc_distance_800_norm[iad]->Add(h_run_acc_distance_scaled[iad],normScale_800[2*(EH-1)+iad]);

			acc_rate_800[2*(EH-1)+iad] = prompt_rate[2*(EH-1)+iad]*delayed_rate[2*(EH-1)+iad]*((800-1.)*1.e-6)*exp(-prompt_rate[2*(EH-1)+iad]*((pd_window_microsec+800)*1.e-6));
			scale_800[2*(EH-1)+iad] = acc_rate_800[2*(EH-1)+iad]*prompt_live[2*(EH-1)+iad]/acc_counts[2*(EH-1)+iad];
			h_total_acc_distance_800_scaled[iad]->Add(h_run_acc_distance_scaled[iad],scale_800[2*(EH-1)+iad]);

		}
		cout << "Done with run_order #" << run_order << " out of " << nRuns << "runs." << endl;

		runFile->Close();
		ibdFile->Close();

	//	if(run_order >= 99) break;
		run[EH-1] +=1;
		run_order += 1;
		
	//	if(run_order == 550) break;
	}

        char outputname[64];
//        sprintf(outputname,"./accResults/TotaledSingles_TcLong_EH%d.root",hall_num);
        sprintf(outputname,"./accResults/TotaledSingles_ktrain_%d_EH%d.root",pd_window_microsec,hall_num);
//        sprintf(outputname,"./accResults/TotaledSingles_4sigma_EH%d.root",hall_num);
	TFile* outfile=new TFile(outputname, "RECREATE");
		outfile->cd();
		for(int iad=0; iad<maxAD; ++iad){
			cout << "Total DAQ time for AD" << iad+1 << ": " << TOTAL_DAQ[2*(EH-1)+iad] << endl;
			cout << "Average scale factor for AD" << iad+1 << ": " << AVERAGE_SCALE[2*(EH-1)+iad] << endl << endl << endl;

			h_total_acc_energy_before[iad]->SetOption("COLZ");
			h_total_acc_energy_before[iad]->SetStats(0);
			h_total_acc_energy_before[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_acc_energy_before[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_acc_energy_before[iad]->Write();

			h_total_acc_energy_1m[iad]->SetOption("COLZ");
			h_total_acc_energy_1m[iad]->SetStats(0);
			h_total_acc_energy_1m[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_acc_energy_1m[iad]->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_acc_energy_1m[iad]->Write();

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

			//h_total_delayed_energy_raw[iad]->SetStats(0);
			h_total_delayed_energy_raw[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_delayed_energy_raw[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_raw[iad]->Write();

			//h_total_delayed_energy_before[iad]->SetStats(0);
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

				h_total_delayed_energy_scaled_z[iad][iz]->SetStats(0);
				h_total_delayed_energy_scaled_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_scaled_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_scaled_z[iad][iz]->Write();

				h_total_delayed_energy_scaled_DT800_z[iad][iz]->SetStats(0);
				h_total_delayed_energy_scaled_DT800_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_scaled_DT800_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_scaled_DT800_z[iad][iz]->Write();

				h_total_delayed_energy_norm_z[iad][iz]->SetStats(0);
				h_total_delayed_energy_norm_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_norm_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_norm_z[iad][iz]->Write();

				h_total_delayed_energy_norm_DT800_z[iad][iz]->SetStats(0);
				h_total_delayed_energy_norm_DT800_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_norm_DT800_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_norm_DT800_z[iad][iz]->Write();

				h_total_delayed_energy_DTnorm_z[iad][iz]->SetStats(0);
				h_total_delayed_energy_DTnorm_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_DTnorm_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_DTnorm_z[iad][iz]->Write();

				h_total_delayed_energy_DTnorm_DT800_z[iad][iz]->SetStats(0);
				h_total_delayed_energy_DTnorm_DT800_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_DTnorm_DT800_z[iad][iz]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_DTnorm_DT800_z[iad][iz]->Write();
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

				h_total_delayed_energy_scaled_r2[iad][ir2]->SetStats(0);
				h_total_delayed_energy_scaled_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_scaled_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_scaled_r2[iad][ir2]->Write();

				h_total_delayed_energy_scaled_DT800_r2[iad][ir2]->SetStats(0);
				h_total_delayed_energy_scaled_DT800_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_scaled_DT800_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_scaled_DT800_r2[iad][ir2]->Write();

				h_total_delayed_energy_norm_r2[iad][ir2]->SetStats(0);
				h_total_delayed_energy_norm_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_norm_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_norm_r2[iad][ir2]->Write();

				h_total_delayed_energy_norm_DT800_r2[iad][ir2]->SetStats(0);
				h_total_delayed_energy_norm_DT800_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_norm_DT800_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_norm_DT800_r2[iad][ir2]->Write();

				h_total_delayed_energy_DTnorm_r2[iad][ir2]->SetStats(0);
				h_total_delayed_energy_DTnorm_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_DTnorm_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_DTnorm_r2[iad][ir2]->Write();

				h_total_delayed_energy_DTnorm_DT800_r2[iad][ir2]->SetStats(0);
				h_total_delayed_energy_DTnorm_DT800_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_DTnorm_DT800_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_DTnorm_DT800_r2[iad][ir2]->Write();
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

					h_total_delayed_energy_scaled_zVSr2[iad][ir2][iz]->SetStats(0);
					h_total_delayed_energy_scaled_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
					h_total_delayed_energy_scaled_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_total_delayed_energy_scaled_zVSr2[iad][ir2][iz]->Write();

					h_total_delayed_energy_scaled_DT800_zVSr2[iad][ir2][iz]->SetStats(0);
					h_total_delayed_energy_scaled_DT800_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
					h_total_delayed_energy_scaled_DT800_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_total_delayed_energy_scaled_DT800_zVSr2[iad][ir2][iz]->Write();

					h_total_delayed_energy_norm_zVSr2[iad][ir2][iz]->SetStats(0);
					h_total_delayed_energy_norm_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
					h_total_delayed_energy_norm_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_total_delayed_energy_norm_zVSr2[iad][ir2][iz]->Write();

					h_total_delayed_energy_norm_DT800_zVSr2[iad][ir2][iz]->SetStats(0);
					h_total_delayed_energy_norm_DT800_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
					h_total_delayed_energy_norm_DT800_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_total_delayed_energy_norm_DT800_zVSr2[iad][ir2][iz]->Write();

					h_total_delayed_energy_DTnorm_zVSr2[iad][ir2][iz]->SetStats(0);
					h_total_delayed_energy_DTnorm_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
					h_total_delayed_energy_DTnorm_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_total_delayed_energy_DTnorm_zVSr2[iad][ir2][iz]->Write();

					h_total_delayed_energy_DTnorm_DT800_zVSr2[iad][ir2][iz]->SetStats(0);
					h_total_delayed_energy_DTnorm_DT800_zVSr2[iad][ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
					h_total_delayed_energy_DTnorm_DT800_zVSr2[iad][ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_total_delayed_energy_DTnorm_DT800_zVSr2[iad][ir2][iz]->Write();
				}
			}

			//h_total_delayed_energy_DT800[iad]->SetStats(0);
			h_total_delayed_energy_DT800[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_DT800[iad]->Write();

			h_total_delayed_energy_after[iad]->SetStats(0);
			h_total_delayed_energy_after[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_after[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_after[iad]->Write();

			//h_total_prompt_energy_before_raw[iad]->SetStats(0);
			h_total_prompt_energy_before_raw[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_before_raw[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_before_raw[iad]->Write();

			//h_total_prompt_energy_before[iad]->SetStats(0);
			h_total_prompt_energy_before[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_before[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_before[iad]->Write();

			//h_total_prompt_energy_DT800[iad]->SetStats(0);
			h_total_prompt_energy_DT800[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800[iad]->Write();

			h_total_prompt_energy_after[iad]->SetStats(0);
			h_total_prompt_energy_after[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_after[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_after[iad]->Write();

			//h_total_acc_distance_before[iad]->SetStats(0);
			h_total_acc_distance_before[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_before[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_before[iad]->Write();

			//h_total_acc_distance_after[iad]->SetStats(0);
			h_total_acc_distance_after[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_after[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_after[iad]->Write();

			//h_total_acc_distance_scaled[iad]->SetStats(0);
			h_total_acc_distance_scaled[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_scaled[iad]->Write();

			//h_total_acc_distance_norm[iad]->SetStats(0);
			h_total_acc_distance_norm[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_norm[iad]->Write();

			//h_total_acc_distance_3sig_scaled[iad]->SetStats(0);
			h_total_acc_distance_3sig_scaled[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_3sig_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_3sig_scaled[iad]->Write();

			//h_total_acc_distance_3sig_norm[iad]->SetStats(0);
			h_total_acc_distance_3sig_norm[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_3sig_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_3sig_norm[iad]->Write();

			//h_total_acc_distance_3sig_DTnorm[iad]->SetStats(0);
			h_total_acc_distance_3sig_DTnorm[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_3sig_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_3sig_DTnorm[iad]->Write();

			//h_total_acc_distance_Ep35_scaled[iad]->SetStats(0);
			h_total_acc_distance_Ep35_scaled[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_Ep35_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_Ep35_scaled[iad]->Write();

			//h_total_acc_distance_Ep35_norm[iad]->SetStats(0);
			h_total_acc_distance_Ep35_norm[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_Ep35_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_Ep35_norm[iad]->Write();



			//h_total_acc_distance_400_scaled[iad]->SetStats(0);
			h_total_acc_distance_400_scaled[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_400_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_400_scaled[iad]->Write();

			//h_total_acc_distance_400_scaled[iad]->SetStats(0);
			h_total_acc_distance_400_norm[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_400_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_400_norm[iad]->Write();

			//h_total_acc_distance_600_scaled[iad]->SetStats(0);
			h_total_acc_distance_600_scaled[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_600_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_600_scaled[iad]->Write();

			//h_total_acc_distance_600_scaled[iad]->SetStats(0);
			h_total_acc_distance_600_norm[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_600_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_600_norm[iad]->Write();

			//h_total_acc_distance_800_scaled[iad]->SetStats(0);
			h_total_acc_distance_800_scaled[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_800_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_800_scaled[iad]->Write();

			//h_total_acc_distance_800_scaled[iad]->SetStats(0);
			h_total_acc_distance_800_norm[iad]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_total_acc_distance_800_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_acc_distance_800_norm[iad]->Write();



			h_total_delayed_energy_scaled[iad]->SetStats(0);
			h_total_delayed_energy_scaled[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_scaled[iad]->Write();

			h_total_prompt_energy_scaled[iad]->SetStats(0);
			h_total_prompt_energy_scaled[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_scaled[iad]->Write();

			h_total_delayed_energy_norm[iad]->SetStats(0);
			h_total_delayed_energy_norm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_norm[iad]->Write();

			h_total_prompt_energy_norm[iad]->SetStats(0);
			h_total_prompt_energy_norm[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_norm[iad]->Write();

			h_total_delayed_energy_DTnorm[iad]->SetStats(0);
			h_total_delayed_energy_DTnorm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_DTnorm[iad]->Write();

			h_total_prompt_energy_DTnorm[iad]->SetStats(0);
			h_total_prompt_energy_DTnorm[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DTnorm[iad]->Write();

			h_total_delayed_energy_DT800_scaled[iad]->SetStats(0);
			h_total_delayed_energy_DT800_scaled[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_DT800_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_DT800_scaled[iad]->Write();

			h_total_prompt_energy_DT800_scaled[iad]->SetStats(0);
			h_total_prompt_energy_DT800_scaled[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800_scaled[iad]->Write();

			h_total_delayed_energy_DT800_norm[iad]->SetStats(0);
			h_total_delayed_energy_DT800_norm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_DT800_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_DT800_norm[iad]->Write();

			h_total_prompt_energy_DT800_norm[iad]->SetStats(0);
			h_total_prompt_energy_DT800_norm[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800_norm[iad]->Write();

			h_total_delayed_energy_DT800_DTnorm[iad]->SetStats(0);
			h_total_delayed_energy_DT800_DTnorm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_DT800_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_DT800_DTnorm[iad]->Write();

			h_total_prompt_energy_DT800_DTnorm[iad]->SetStats(0);
			h_total_prompt_energy_DT800_DTnorm[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800_DTnorm[iad]->Write();

			h_total_delayed_energy_fine_scaled[iad]->SetStats(0);
			h_total_delayed_energy_fine_scaled[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_scaled[iad]->Write();

			h_total_prompt_energy_DT800_3sig_scaled[iad]->SetStats(0);
			h_total_prompt_energy_DT800_3sig_scaled[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800_3sig_scaled[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800_3sig_scaled[iad]->Write();

				h_total_prompt_energy_DT300_3sig_scaled[iad]->SetStats(0);
				h_total_prompt_energy_DT300_3sig_scaled[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_DT300_3sig_scaled[iad]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_DT300_3sig_scaled[iad]->Write();

				h_total_prompt_energy_DT500_3sig_scaled[iad]->SetStats(0);
				h_total_prompt_energy_DT500_3sig_scaled[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_DT500_3sig_scaled[iad]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_DT500_3sig_scaled[iad]->Write();

				h_total_prompt_energy_DT1000_3sig_scaled[iad]->SetStats(0);
				h_total_prompt_energy_DT1000_3sig_scaled[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_DT1000_3sig_scaled[iad]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_DT1000_3sig_scaled[iad]->Write();

				h_total_prompt_energy_DT1500_3sig_scaled[iad]->SetStats(0);
				h_total_prompt_energy_DT1500_3sig_scaled[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_DT1500_3sig_scaled[iad]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_DT1500_3sig_scaled[iad]->Write();

				h_total_prompt_energy_DT2000_3sig_scaled[iad]->SetStats(0);
				h_total_prompt_energy_DT2000_3sig_scaled[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_DT2000_3sig_scaled[iad]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_DT2000_3sig_scaled[iad]->Write();

			h_total_prompt_energy_DT800_3sig_norm[iad]->SetStats(0);
			h_total_prompt_energy_DT800_3sig_norm[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800_3sig_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800_3sig_norm[iad]->Write();

			h_total_prompt_energy_DT800_3sig_DTnorm[iad]->SetStats(0);
			h_total_prompt_energy_DT800_3sig_DTnorm[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_DT800_3sig_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_DT800_3sig_DTnorm[iad]->Write();

			h_total_delayed_energy_fine_norm[iad]->SetStats(0);
			h_total_delayed_energy_fine_norm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_norm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_norm[iad]->Write();

			h_total_delayed_energy_fine_DTnorm[iad]->SetStats(0);
			h_total_delayed_energy_fine_DTnorm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DTnorm[iad]->Write();

			h_total_delayed_energy_fine_DT800_scaled_1[iad]->SetStats(0);
			h_total_delayed_energy_fine_DT800_scaled_1[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DT800_scaled_1[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DT800_scaled_1[iad]->Write();

			h_total_delayed_energy_fine_DT800_norm_1[iad]->SetStats(0);
			h_total_delayed_energy_fine_DT800_norm_1[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DT800_norm_1[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DT800_norm_1[iad]->Write();

			h_total_delayed_energy_fine_DT800_DTnorm_1[iad]->SetStats(0);
			h_total_delayed_energy_fine_DT800_DTnorm_1[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DT800_DTnorm_1[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DT800_DTnorm_1[iad]->Write();
			
			h_total_delayed_energy_fine_DT800_scaled_2[iad]->SetStats(0);
			h_total_delayed_energy_fine_DT800_scaled_2[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DT800_scaled_2[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DT800_scaled_2[iad]->Write();

			h_total_delayed_energy_fine_DT800_norm_2[iad]->SetStats(0);
			h_total_delayed_energy_fine_DT800_norm_2[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DT800_norm_2[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DT800_norm_2[iad]->Write();

			h_total_delayed_energy_fine_DT800_DTnorm_2[iad]->SetStats(0);
			h_total_delayed_energy_fine_DT800_DTnorm_2[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_fine_DT800_DTnorm_2[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_fine_DT800_DTnorm_2[iad]->Write();

				h_total_delayed_energy_fine_Ep35_scaled[iad]->SetStats(0);
				h_total_delayed_energy_fine_Ep35_scaled[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_fine_Ep35_scaled[iad]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_fine_Ep35_scaled[iad]->Write();

				h_total_delayed_energy_fine_Ep35_norm[iad]->SetStats(0);
				h_total_delayed_energy_fine_Ep35_norm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_fine_Ep35_norm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_fine_Ep35_norm[iad]->Write();

				h_total_delayed_energy_fine_Ep35_DTnorm[iad]->SetStats(0);
				h_total_delayed_energy_fine_Ep35_DTnorm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_fine_Ep35_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_fine_Ep35_DTnorm[iad]->Write();

				h_total_delayed_energy_fine_DT800_Ep35_scaled[iad]->SetStats(0);
				h_total_delayed_energy_fine_DT800_Ep35_scaled[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_fine_DT800_Ep35_scaled[iad]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_fine_DT800_Ep35_scaled[iad]->Write();

				h_total_delayed_energy_fine_DT800_Ep35_norm[iad]->SetStats(0);
				h_total_delayed_energy_fine_DT800_Ep35_norm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_fine_DT800_Ep35_norm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_fine_DT800_Ep35_norm[iad]->Write();

				h_total_delayed_energy_fine_DT800_Ep35_DTnorm[iad]->SetStats(0);
				h_total_delayed_energy_fine_DT800_Ep35_DTnorm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_fine_DT800_Ep35_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_fine_DT800_Ep35_DTnorm[iad]->Write();

			h_total_delayed_energy_largeDist[iad]->SetStats(0);
			h_total_delayed_energy_largeDist[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_largeDist[iad]->Write();

			h_total_prompt_energy_largeDist[iad]->SetStats(0);
			h_total_prompt_energy_largeDist[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_largeDist[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_largeDist[iad]->Write();

			h_total_delayed_energy_shortDist[iad]->SetStats(0);
			h_total_delayed_energy_shortDist[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_shortDist[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_shortDist[iad]->Write();

			h_total_prompt_energy_shortDist[iad]->SetStats(0);
			h_total_prompt_energy_shortDist[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_shortDist[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_shortDist[iad]->Write();


			h_p_muon_efficiency[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_p_muon_efficiency[iad]->GetYaxis()->SetTitle("Muon Efficiency for Prompt Singles(%)");
			h_p_muon_efficiency[iad]->Write();

			h_d_muon_efficiency[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_d_muon_efficiency[iad]->GetYaxis()->SetTitle("Muon Efficiency for Delayed Singles(%)");
			h_d_muon_efficiency[iad]->Write();

			h_p_rate_before[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_p_rate_before[iad]->GetYaxis()->SetTitle("Prompt Rate [Hz]");
			h_p_rate_before[iad]->Write();

			h_d_rate_before[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_d_rate_before[iad]->GetYaxis()->SetTitle("Delayed Rate [Hz]");
			h_d_rate_before[iad]->Write();

			h_p_rate[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_p_rate[iad]->GetYaxis()->SetTitle("Prompt Rate [Hz]");
			h_p_rate[iad]->Write();

			h_d_rate[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_d_rate[iad]->GetYaxis()->SetTitle("Delayed Rate [Hz]");
			h_d_rate[iad]->Write();

			h_pd_ratio[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_pd_ratio[iad]->GetYaxis()->SetTitle("Ratio of Rates (Prompt/Delayed)");
			h_pd_ratio[iad]->Write();

			h_acc_rate[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_acc_rate[iad]->GetYaxis()->SetTitle("Accidental Rate [Hz]");
			h_acc_rate[iad]->Write();

			h_DAQ[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_DAQ[iad]->GetYaxis()->SetTitle("DAQ time [s]");
			h_DAQ[iad]->Write();

			h_scale[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_scale[iad]->GetYaxis()->SetTitle("Rate Corrected Scale");
			h_scale[iad]->Write();

			h_normScale[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_normScale[iad]->GetYaxis()->SetTitle("Normalized Scale");
			h_normScale[iad]->Write();

			h_diffScales_time[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_diffScales_time[iad]->GetYaxis()->SetTitle("Difference in scales 100% * (DTnorm-Rate)/DTnorm");
			h_diffScales_time[iad]->Write();

			h_diffScales_sigma_time[iad]->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_diffScales_sigma_time[iad]->GetYaxis()->SetTitle("Difference in scales (Rate-Norm)/sigma");
			h_diffScales_sigma_time[iad]->Write();

			h_diffScales[iad]->GetXaxis()->SetTitle("Difference in scales (Rate-Norm)");
			h_diffScales[iad]->GetYaxis()->SetTitle("Counts");
			h_diffScales[iad]->Write();

			h_diffScales_sigma[iad]->GetXaxis()->SetTitle("Difference in scales (Rate-Norm)/sigma");
			h_diffScales_sigma[iad]->GetYaxis()->SetTitle("Counts");
			h_diffScales_sigma[iad]->Write();

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

			h_total_delayed_energy_scaled_p35[iad]->SetStats(0);
			h_total_delayed_energy_scaled_p35[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_scaled_p35[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_scaled_p35[iad]->Write();

			h_total_delayed_energy_vs_distance_rate[iad]->SetStats(0);
			h_total_delayed_energy_vs_distance_rate[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_vs_distance_rate[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_delayed_energy_vs_distance_rate[iad]->SetOption("COLZ");
			h_total_delayed_energy_vs_distance_rate[iad]->Write();

			h_total_delayed_energy_vs_distance_norm[iad]->SetStats(0);
			h_total_delayed_energy_vs_distance_norm[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_vs_distance_norm[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_delayed_energy_vs_distance_norm[iad]->SetOption("COLZ");
			h_total_delayed_energy_vs_distance_norm[iad]->Write();

			h_total_prompt_energy_vs_distance_rate[iad]->SetStats(0);
			h_total_prompt_energy_vs_distance_rate[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_vs_distance_rate[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_prompt_energy_vs_distance_rate[iad]->SetOption("COLZ");
			h_total_prompt_energy_vs_distance_rate[iad]->Write();

			h_total_prompt_energy_vs_distance_norm[iad]->SetStats(0);
			h_total_prompt_energy_vs_distance_norm[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_vs_distance_norm[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_prompt_energy_vs_distance_norm[iad]->SetOption("COLZ");
			h_total_prompt_energy_vs_distance_norm[iad]->Write();

			for(int dist=0; dist<6; dist++){
				h_total_delayed_energy_dist[iad][dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_dist[iad][dist]->Write();
				h_total_prompt_energy_dist[iad][dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_dist[iad][dist]->Write();

				h_total_delayed_energy_scaled_dist[iad][dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_scaled_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_scaled_dist[iad][dist]->Write();
				h_total_prompt_energy_scaled_dist[iad][dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_scaled_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_scaled_dist[iad][dist]->Write();

				h_total_delayed_energy_norm_dist[iad][dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_total_delayed_energy_norm_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_total_delayed_energy_norm_dist[iad][dist]->Write();
				h_total_prompt_energy_norm_dist[iad][dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_total_prompt_energy_norm_dist[iad][dist]->GetYaxis()->SetTitle("Counts");
				h_total_prompt_energy_norm_dist[iad][dist]->Write();
			}
/*			h_total_delayed_energy_dist_2plus[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_dist_2plus[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_dist_2plus[iad]->Write();

			h_total_prompt_energy_dist_2plus[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_dist_2plus[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_dist_2plus[iad]->Write();

			h_total_delayed_energy_scaled_dist_2plus[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_scaled_dist_2plus[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_scaled_dist_2plus[iad]->Write();

			h_total_prompt_energy_scaled_dist_2plus[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_scaled_dist_2plus[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_scaled_dist_2plus[iad]->Write();

			h_total_delayed_energy_norm_dist_2plus[iad]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_total_delayed_energy_norm_dist_2plus[iad]->GetYaxis()->SetTitle("Counts");
			h_total_delayed_energy_norm_dist_2plus[iad]->Write();

			h_total_prompt_energy_norm_dist_2plus[iad]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_prompt_energy_norm_dist_2plus[iad]->GetYaxis()->SetTitle("Counts");
			h_total_prompt_energy_norm_dist_2plus[iad]->Write();*/

			h_total_acc_distVStime[iad]->SetStats(0);
			h_total_acc_distVStime[iad]->SetOption("COLZ");
			h_total_acc_distVStime[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_acc_distVStime[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_acc_distVStime[iad]->Write();

			h_total_acc_distVStime_norm[iad]->SetStats(0);
			h_total_acc_distVStime_norm[iad]->SetOption("COLZ");
			h_total_acc_distVStime_norm[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_acc_distVStime_norm[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_acc_distVStime_norm[iad]->Write();

			h_total_acc_distVStime_rate[iad]->SetStats(0);
			h_total_acc_distVStime_rate[iad]->SetOption("COLZ");
			h_total_acc_distVStime_rate[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_acc_distVStime_rate[iad]->GetYaxis()->SetTitle("Distance [m]");
			h_total_acc_distVStime_rate[iad]->Write();

				h_total_acc_distVStime_Ep35[iad]->SetStats(0);
				h_total_acc_distVStime_Ep35[iad]->SetOption("COLZ");
				h_total_acc_distVStime_Ep35[iad]->GetXaxis()->SetTitle("Delta Time [us]");
				h_total_acc_distVStime_Ep35[iad]->GetYaxis()->SetTitle("Distance [m]");
				h_total_acc_distVStime_Ep35[iad]->Write();

				h_total_acc_distVStime_Ep35_norm[iad]->SetStats(0);
				h_total_acc_distVStime_Ep35_norm[iad]->SetOption("COLZ");
				h_total_acc_distVStime_Ep35_norm[iad]->GetXaxis()->SetTitle("Delta Time [us]");
				h_total_acc_distVStime_Ep35_norm[iad]->GetYaxis()->SetTitle("Distance [m]");
				h_total_acc_distVStime_Ep35_norm[iad]->Write();

				h_total_acc_distVStime_Ep35_rate[iad]->SetStats(0);
				h_total_acc_distVStime_Ep35_rate[iad]->SetOption("COLZ");
				h_total_acc_distVStime_Ep35_rate[iad]->GetXaxis()->SetTitle("Delta Time [us]");
				h_total_acc_distVStime_Ep35_rate[iad]->GetYaxis()->SetTitle("Distance [m]");
				h_total_acc_distVStime_Ep35_rate[iad]->Write();

		//		h_total_acc_DT[iad]->SetStats(0);
				h_total_acc_DT[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT[iad]->Write();

		//		h_total_acc_DT_rate[iad]->SetStats(0);
				h_total_acc_DT_rate[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_rate[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_rate[iad]->Write();

		//		h_total_acc_DT_norm[iad]->SetStats(0);
				h_total_acc_DT_norm[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_norm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_norm[iad]->Write();

		//		h_total_acc_DT_DTnorm[iad]->SetStats(0);
				h_total_acc_DT_DTnorm[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_DTnorm[iad]->Write();

		//		h_total_acc_DT_3sig_rate[iad]->SetStats(0);
				h_total_acc_DT_3sig_rate[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_3sig_rate[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_3sig_rate[iad]->Write();

		//		h_total_acc_DT_3sig_norm[iad]->SetStats(0);
				h_total_acc_DT_3sig_norm[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_3sig_norm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_3sig_norm[iad]->Write();

		//		h_total_acc_DT_3sig_DTnorm[iad]->SetStats(0);
				h_total_acc_DT_3sig_DTnorm[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_3sig_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_3sig_DTnorm[iad]->Write();

		//		h_total_acc_DT_Ep35[iad]->SetStats(0);
				h_total_acc_DT_Ep35[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_Ep35[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_Ep35[iad]->Write();

		//		h_total_acc_DT_Ep35_norm[iad]->SetStats(0);
				h_total_acc_DT_Ep35_norm[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_Ep35_norm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_Ep35_norm[iad]->Write();

		//		h_total_acc_DT_Ep35_rate[iad]->SetStats(0);
				h_total_acc_DT_Ep35_rate[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_Ep35_rate[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_Ep35_rate[iad]->Write();

		//		h_total_acc_DT_Ep35_DTnorm[iad]->SetStats(0);
				h_total_acc_DT_Ep35_DTnorm[iad]->GetXaxis()->SetTitle("DT [m]");
				h_total_acc_DT_Ep35_DTnorm[iad]->GetYaxis()->SetTitle("Counts");
				h_total_acc_DT_Ep35_DTnorm[iad]->Write();

			h_total_acc_promptVStime_DTnorm[iad]->SetStats(0);
			h_total_acc_promptVStime_DTnorm[iad]->SetOption("COLZ");
			h_total_acc_promptVStime_DTnorm[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_acc_promptVStime_DTnorm[iad]->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_acc_promptVStime_DTnorm[iad]->Write();

			h_total_acc_promptVStime_DT800_DTnorm[iad]->SetStats(0);
			h_total_acc_promptVStime_DT800_DTnorm[iad]->SetOption("COLZ");
			h_total_acc_promptVStime_DT800_DTnorm[iad]->GetXaxis()->SetTitle("Delta Time [us]");
			h_total_acc_promptVStime_DT800_DTnorm[iad]->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_total_acc_promptVStime_DT800_DTnorm[iad]->Write();
			

		}

	for(int iad=0; iad < maxAD; iad++){
		cout << endl << "Total eff_m for EH" << EH << " AD" << iad+1 << ":\t" << TOTAL_Eff_mult[2*(EH-1)+iad] << endl;
		cout << "Total R_s,prompt [Hz] for EH" << EH << " AD" << iad+1 << ":\t" << TOTAL_prompt_rate[2*(EH-1)+iad] << endl;
		cout << "Total R_s,delayed [Hz] (1.5-3MeV) for EH" << EH << " AD" << iad+1 << ":\t" << TOTAL_delayed_rate[2*(EH-1)+iad] << endl;
		cout << "Total N_Acc for EH" << EH << " AD" << iad+1 << ":\t" << h_total_prompt_energy_DT800_3sig_scaled[iad]->Integral() << endl;
		cout << "Total R_Acc [1/d] for EH" << EH << " AD" << iad+1 << ":\t" << ((h_total_prompt_energy_DT800_3sig_scaled[iad]->Integral())/(TOTAL_LIVE[2*(EH-1)+iad]/(60*60*24)))/(TOTAL_Eff_mult[2*(EH-1)+iad]) << endl;
		cout << "Total Percent Error on Prompt Singles Rate (uncorrected) for EH" << EH << " AD" << iad+1 << ":\t" << sqrt(h_total_prompt_energy_before_raw[iad]->Integral())/(h_total_prompt_energy_before_raw[iad]->Integral())*100 << endl << endl;
	}

	outfile->Close();

	resub.close();

}
