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
#include "TMath.h"
#include <typeinfo>

using namespace std;

/*Double_t efficiency(Double_t x, double p1, double tau_ls, double p3, double tau_gdls){
		TF1* fit=new TF1("fit","[0]*TMath::Exp(-x/[1])+[2]*TMath::Exp(-x/[3])",0,1.e7);
		fit->FixParameter(0,p1);
		fit->FixParameter(1,tau_ls);
		fit->FixParameter(2,p3);
		fit->FixParameter(3,tau_gdls);
   Double_t num = fit->Integral(0,x);
   Double_t denom = fit->Integral(0,1.e6);
   return num/denom;
}*/

void subtract(int hall_num, int ad_num, int pd_window_microsec){

	int nRuns = 0;
	if(hall_num == 1){
//		nRuns = 1134;
		nRuns = 1035;
	}
	if(hall_num == 2){
//		nRuns = 1125;
		nRuns = 1034;
	}
	if(hall_num == 3){
//		nRuns = 1215;
		nRuns = 1099;
	}

	const int NzBins = 20;
	const int Nr2Bins = 20;

	char name[64];

//	sprintf(name, "h_scales_ad%d", ad_num);
//	TH1F* h_scales=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);

	sprintf(name, "h_rate_scale_ad%d", ad_num);
	TH1F* h_rateScale=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);

	sprintf(name, "h_norm_scale_ad%d", ad_num);
	TH1F* h_normScale=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);

	sprintf(name, "h_ratio_scale_ad%d", ad_num);
	TH1F* h_ratioScale=new TH1F(name,name,nRuns,-0.5,nRuns-0.5);

	sprintf(name, "h_scale_diff_ad%d", ad_num);
	TH1F* h_scale_diff=new TH1F(name,name,51,-0.00931,0.00931);

	//Making the Histograms
	sprintf(name, "h_IBDDistance_ad%d", ad_num);
	TH1F* h_IBDDistance=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_accDistance_rate_ad%d", ad_num);
	TH1F* h_accDistance_rate=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_accDistance_norm_ad%d", ad_num);
	TH1F* h_accDistance_norm=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_normSub_distance_ad%d", ad_num);
	TH1F* h_normSub_distance=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_rateSub_distance_ad%d", ad_num);
	TH1F* h_rateSub_distance=new TH1F(name,name,700,0,7.); //distance histogram

		//400 subset
		sprintf(name, "h_normSub_distance_400_ad%d", ad_num);
		TH1F* h_normSub_distance_400=new TH1F(name,name,700,0,7.); //distance histogram

		sprintf(name, "h_rateSub_distance_400_ad%d", ad_num);
		TH1F* h_rateSub_distance_400=new TH1F(name,name,700,0,7.); //distance histogram

		//600 subset
		sprintf(name, "h_normSub_distance_600_ad%d", ad_num);
		TH1F* h_normSub_distance_600=new TH1F(name,name,700,0,7.); //distance histogram

		sprintf(name, "h_rateSub_distance_600_ad%d", ad_num);
		TH1F* h_rateSub_distance_600=new TH1F(name,name,700,0,7.); //distance histogram

		//800 subset
		sprintf(name, "h_normSub_distance_800_ad%d", ad_num);
		TH1F* h_normSub_distance_800=new TH1F(name,name,700,0,7.); //distance histogram

		sprintf(name, "h_rateSub_distance_800_ad%d", ad_num);
		TH1F* h_rateSub_distance_800=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_IBDDistance_3sig_ad%d", ad_num);
	TH1F* h_IBDDistance_3sig=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_accDistance_3sig_rate_ad%d", ad_num);
	TH1F* h_accDistance_3sig_rate=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_accDistance_3sig_norm_ad%d", ad_num);
	TH1F* h_accDistance_3sig_norm=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_accDistance_3sig_DTnorm_ad%d", ad_num);
	TH1F* h_accDistance_3sig_DTnorm=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_normSub_distance_3sig_ad%d", ad_num);
	TH1F* h_normSub_distance_3sig=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_rateSub_distance_3sig_ad%d", ad_num);
	TH1F* h_rateSub_distance_3sig=new TH1F(name,name,700,0,7.); //distance histogram

	sprintf(name, "h_DTnormSub_distance_3sig_ad%d", ad_num);
	TH1F* h_DTnormSub_distance_3sig=new TH1F(name,name,700,0,7.); //distance histogram

	int time = 0;
			TH1F* h_rate_modified[4][5]; //staged modification of rate subtraction
			for(int tim=0; tim<4; tim++){
				for(int stage=0; stage<5; stage++){
					if(tim == 0) time = 2000;
					if(tim == 1) time = 400;
					if(tim == 2) time = 600;
					if(tim == 3) time = 800;
					sprintf(name, "h_rate_modified_ad%d_%d_stage%d", ad_num, time, stage+1);
					h_rate_modified[tim][stage]=new TH1F(name,name,700,0,7.);
				}
			}



	sprintf(name, "h_dist_ratio_rateTOnorm_ad%d", ad_num);
	TH1F* h_dist_ratio_rateTOnorm=new TH1F(name,name,700,0,7.); //ratio histogram

	sprintf(name, "h_dist_sub_rateTOnorm_ad%d", ad_num);
	TH1F* h_dist_sub_rateTOnorm=new TH1F(name,name,700,0,7.); //acc sub histogram

	sprintf(name, "h_Eprompt_subtract_ad%d", ad_num);
	TH1F* h_Eprompt_sub=new TH1F(name,name,113,0.7,12.);

	sprintf(name, "h_Edelayed_IBD_ad%d", ad_num);
	TH1F* h_Edelayed_IBD=new TH1F(name,name,230,0.7,3.);

	sprintf(name, "h_Edelayed_IBD_DT800_ad%d", ad_num);
	TH1F* h_Edelayed_IBD_DT800=new TH1F(name,name,230,0.7,3.);

	sprintf(name, "h_Edelayed_IBD_fine_ad%d", ad_num);
	TH1F* h_Edelayed_IBD_fine=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_IBD_fine_DT800_ad%d", ad_num);
	TH1F* h_Edelayed_IBD_fine_DT800=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_IBD_fine_Ep35_ad%d", ad_num);
	TH1F* h_Edelayed_IBD_fine_Ep35=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_IBD_fine_DT800_Ep35_ad%d", ad_num);
	TH1F* h_Edelayed_IBD_fine_DT800_Ep35=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_subtract_ad%d", ad_num);
//	TH1F* h_Edelayed_sub=new TH1F(name,name,150,1.5,3.);
	TH1F* h_Edelayed_sub=new TH1F(name,name,230,0.7,3.);

	sprintf(name, "h_Eprompt_subtract_norm_ad%d", ad_num);
	TH1F* h_Eprompt_sub_norm=new TH1F(name,name,113,0.7,12.);

	sprintf(name, "h_EpromptToDelayed_norm_ad%d", ad_num);
	TH1F* h_EpromptToDelayed_norm=new TH1F(name,name,113,0.7,12.);

//	sprintf(name, "h_Eprompt_sub_noDTminusDT800_norm_ad%d", ad_num);
//	TH1F* h_Eprompt_sub_noDTminusDT800_norm=new TH1F(name,name,113,0.7,12.);

	sprintf(name, "h_Edelayed_subtract_norm_ad%d", ad_num);
//	TH1F* h_Edelayed_sub_norm=new TH1F(name,name,150,1.5,3.);
	TH1F* h_Edelayed_sub_norm=new TH1F(name,name,230,0.7,3.);

	sprintf(name, "h_Eprompt_subtract_DTnorm_ad%d", ad_num);
	TH1F* h_Eprompt_sub_DTnorm=new TH1F(name,name,113,0.7,12.);

	sprintf(name, "h_Edelayed_subtract_DTnorm_ad%d", ad_num);
	TH1F* h_Edelayed_sub_DTnorm=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_Eprompt_subtract_DT800_ad%d", ad_num);
		TH1F* h_Eprompt_sub_DT800=new TH1F(name,name,113,0.7,12.);

			sprintf(name, "h_Eprompt_subtract_DT800_3sig_ad%d", ad_num);
			TH1F* h_Eprompt_sub_DT800_3sig=new TH1F(name,name,113,0.7,12.);

		sprintf(name, "h_Edelayed_subtract_DT800_ad%d", ad_num);
		TH1F* h_Edelayed_sub_DT800=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_Eprompt_subtract_DT800_norm_ad%d", ad_num);
		TH1F* h_Eprompt_sub_DT800_norm=new TH1F(name,name,113,0.7,12.);

		sprintf(name, "h_Edelayed_subtract_DT800_norm_ad%d", ad_num);
		TH1F* h_Edelayed_sub_DT800_norm=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_Eprompt_subtract_DT800_DTnorm_ad%d", ad_num);
		TH1F* h_Eprompt_sub_DT800_DTnorm=new TH1F(name,name,113,0.7,12.);

		sprintf(name, "h_Edelayed_subtract_DT800_DTnorm_ad%d", ad_num);
		TH1F* h_Edelayed_sub_DT800_DTnorm=new TH1F(name,name,230,0.7,3.);

		TH1F* h_Edelayed_ibd_z[NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			h_Edelayed_ibd_z[iz]=new TH1F(Form("h_Edelayed_ibd_z_ad%d_iz%d", ad_num, iz+1),Form("h_Edelayed_ibd_z_ad%d_iz%d", ad_num, iz+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_ibd_DT800_z[NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			h_Edelayed_ibd_DT800_z[iz]=new TH1F(Form("h_Edelayed_ibd_DT800_z_ad%d_iz%d", ad_num, iz+1),Form("h_Edelayed_ibd_DT800_z_ad%d_iz%d", ad_num, iz+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_rate_z[NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			h_Edelayed_sub_rate_z[iz]=new TH1F(Form("h_Edelayed_subtract_rate_z_ad%d_iz%d", ad_num, iz+1),Form("h_Edelayed_subtract_rate_z_ad%d_iz%d", ad_num, iz+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_norm_z[NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			h_Edelayed_sub_norm_z[iz]=new TH1F(Form("h_Edelayed_subtract_norm_z_ad%d_iz%d", ad_num, iz+1),Form("h_Edelayed_subtract_norm_z_ad%d_iz%d", ad_num, iz+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_DTnorm_z[NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			h_Edelayed_sub_DTnorm_z[iz]=new TH1F(Form("h_Edelayed_subtract_DTnorm_z_ad%d_iz%d", ad_num, iz+1),Form("h_Edelayed_subtract_DTnorm_z_ad%d_iz%d", ad_num, iz+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_rate_DT800_z[NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			h_Edelayed_sub_rate_DT800_z[iz]=new TH1F(Form("h_Edelayed_subtract_rate_DT800_z_ad%d_iz%d", ad_num, iz+1),Form("h_Edelayed_subtract_rate_DT800_z_ad%d_iz%d", ad_num, iz+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_norm_DT800_z[NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			h_Edelayed_sub_norm_DT800_z[iz]=new TH1F(Form("h_Edelayed_subtract_norm_DT800_z_ad%d_iz%d", ad_num, iz+1),Form("h_Edelayed_subtract_norm_DT800_z_ad%d_iz%d", ad_num, iz+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_DTnorm_DT800_z[NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			h_Edelayed_sub_DTnorm_DT800_z[iz]=new TH1F(Form("h_Edelayed_subtract_DTnorm_DT800_z_ad%d_iz%d", ad_num, iz+1),Form("h_Edelayed_subtract_DTnorm_DT800_z_ad%d_iz%d", ad_num, iz+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_ibd_r2[Nr2Bins]; //delayed energy histogram
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			h_Edelayed_ibd_r2[ir2]=new TH1F(Form("h_Edelayed_ibd_r2_ad%d_ir2%d", ad_num, ir2+1),Form("h_Edelayed_ibd_r2_ad%d_ir2%d", ad_num, ir2+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_ibd_DT800_r2[Nr2Bins]; //delayed energy histogram
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			h_Edelayed_ibd_DT800_r2[ir2]=new TH1F(Form("h_Edelayed_ibd_DT800_r2_ad%d_ir2%d", ad_num, ir2+1),Form("h_Edelayed_ibd_DT800_r2_ad%d_ir2%d", ad_num, ir2+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_rate_r2[Nr2Bins]; //delayed energy histogram
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			h_Edelayed_sub_rate_r2[ir2]=new TH1F(Form("h_Edelayed_subtract_rate_r2_ad%d_ir2%d", ad_num, ir2+1),Form("h_Edelayed_subtract_rate_r2_ad%d_ir2%d", ad_num, ir2+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_norm_r2[Nr2Bins]; //delayed energy histogram
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			h_Edelayed_sub_norm_r2[ir2]=new TH1F(Form("h_Edelayed_subtract_norm_r2_ad%d_ir2%d", ad_num, ir2+1),Form("h_Edelayed_subtract_norm_r2_ad%d_ir2%d", ad_num, ir2+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_DTnorm_r2[Nr2Bins]; //delayed energy histogram
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			h_Edelayed_sub_DTnorm_r2[ir2]=new TH1F(Form("h_Edelayed_subtract_DTnorm_r2_ad%d_ir2%d", ad_num, ir2+1),Form("h_Edelayed_subtract_DTnorm_r2_ad%d_ir2%d", ad_num, ir2+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_rate_DT800_r2[Nr2Bins]; //delayed energy histogram
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			h_Edelayed_sub_rate_DT800_r2[ir2]=new TH1F(Form("h_Edelayed_subtract_rate_DT800_r2_ad%d_ir2%d", ad_num, ir2+1),Form("h_Edelayed_subtract_rate_DT800_r2_ad%d_ir2%d", ad_num, ir2+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_norm_DT800_r2[Nr2Bins]; //delayed energy histogram
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			h_Edelayed_sub_norm_DT800_r2[ir2]=new TH1F(Form("h_Edelayed_subtract_norm_DT800_r2_ad%d_ir2%d", ad_num, ir2+1),Form("h_Edelayed_subtract_norm_DT800_r2_ad%d_ir2%d", ad_num, ir2+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_sub_DTnorm_DT800_r2[Nr2Bins]; //delayed energy histogram
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			h_Edelayed_sub_DTnorm_DT800_r2[ir2]=new TH1F(Form("h_Edelayed_subtract_DTnorm_DT800_r2_ad%d_ir2%d", ad_num, ir2+1),Form("h_Edelayed_subtract_DTnorm_DT800_r2_ad%d_ir2%d", ad_num, ir2+1),230,0.7,3.);
		}

		TH1F* h_Edelayed_ibd_zVSr2[Nr2Bins][NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_ibd_zVSr2[ir2][iz]=new TH1F(Form("h_Edelayed_ibd_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),Form("h_Edelayed_ibd_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_Edelayed_ibd_DT800_zVSr2[Nr2Bins][NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_ibd_DT800_zVSr2[ir2][iz]=new TH1F(Form("h_Edelayed_ibd_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),Form("h_Edelayed_ibd_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_Edelayed_sub_rate_zVSr2[Nr2Bins][NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_sub_rate_zVSr2[ir2][iz]=new TH1F(Form("h_Edelayed_subtract_rate_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),Form("h_Edelayed_subtract_rate_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_Edelayed_sub_norm_zVSr2[Nr2Bins][NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_sub_norm_zVSr2[ir2][iz]=new TH1F(Form("h_Edelayed_subtract_norm_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),Form("h_Edelayed_subtract_norm_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_Edelayed_sub_DTnorm_zVSr2[Nr2Bins][NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]=new TH1F(Form("h_Edelayed_subtract_DTnorm_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),Form("h_Edelayed_subtract_DTnorm_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_Edelayed_sub_rate_DT800_zVSr2[Nr2Bins][NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_sub_rate_DT800_zVSr2[ir2][iz]=new TH1F(Form("h_Edelayed_subtract_rate_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),Form("h_Edelayed_subtract_rate_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_Edelayed_sub_norm_DT800_zVSr2[Nr2Bins][NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_sub_norm_DT800_zVSr2[ir2][iz]=new TH1F(Form("h_Edelayed_subtract_norm_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),Form("h_Edelayed_subtract_norm_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),230,0.7,3.);
			}
		}

		TH1F* h_Edelayed_sub_DTnorm_DT800_zVSr2[Nr2Bins][NzBins]; //delayed energy histogram
		for(int iz = 0; iz < NzBins; iz++){
			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_sub_DTnorm_DT800_zVSr2[ir2][iz]=new TH1F(Form("h_Edelayed_subtract_DTnorm_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),Form("h_Edelayed_subtract_DTnorm_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1),230,0.7,3.);
			}
		}



	sprintf(name, "h_Edelayed_subtract_fine_ad%d", ad_num);
	TH1F* h_Edelayed_sub_fine=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_subtract_fine_norm_ad%d", ad_num);
	TH1F* h_Edelayed_sub_fine_norm=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_subtract_fine_DTnorm_ad%d", ad_num);
	TH1F* h_Edelayed_sub_fine_DTnorm=new TH1F(name,name,2300,0.7,3.);

		sprintf(name, "h_Edelayed_subtract_fine_DT800_ad%d", ad_num);
		TH1F* h_Edelayed_sub_fine_DT800=new TH1F(name,name,2300,0.7,3.);

		sprintf(name, "h_Edelayed_subtract_fine_DT800_norm_ad%d", ad_num);
		TH1F* h_Edelayed_sub_fine_DT800_norm=new TH1F(name,name,2300,0.7,3.);

		sprintf(name, "h_Edelayed_subtract_fine_DT800_DTnorm_ad%d", ad_num);
		TH1F* h_Edelayed_sub_fine_DT800_DTnorm=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_subtract_fine_Ep35_ad%d", ad_num);
	TH1F* h_Edelayed_sub_fine_Ep35=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_subtract_fine_Ep35_norm_ad%d", ad_num);
	TH1F* h_Edelayed_sub_fine_Ep35_norm=new TH1F(name,name,2300,0.7,3.);

	sprintf(name, "h_Edelayed_subtract_fine_Ep35_DTnorm_ad%d", ad_num);
	TH1F* h_Edelayed_sub_fine_Ep35_DTnorm=new TH1F(name,name,2300,0.7,3.);

		sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_ad%d", ad_num);
		TH1F* h_Edelayed_sub_fine_DT800_Ep35=new TH1F(name,name,2300,0.7,3.);

		sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_norm_ad%d", ad_num);
		TH1F* h_Edelayed_sub_fine_DT800_Ep35_norm=new TH1F(name,name,2300,0.7,3.);

		sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_DTnorm_ad%d", ad_num);
		TH1F* h_Edelayed_sub_fine_DT800_Ep35_DTnorm=new TH1F(name,name,2300,0.7,3.);



	sprintf(name, "h_delayed_energy_vs_distance_rate_ratio_ad%d", ad_num);
//	TH2F* h_delayed_energy_vs_distance_rate_ratio=new TH2F(name,name,150,1.5,3.);
	TH2F* h_delayed_energy_vs_distance_rate_ratio=new TH2F(name,name,230,0.7,3.,11,0.,5.5);

	sprintf(name, "h_delayed_energy_vs_distance_rate_sub_ad%d", ad_num);
//	TH2F* h_delayed_energy_vs_distance_rate_sub=new TH2F(name,name,150,1.5,3.);
	TH2F* h_delayed_energy_vs_distance_rate_sub=new TH2F(name,name,230,0.7,3.,11,0.,5.5);

	sprintf(name, "h_delayed_energy_vs_distance_norm_ratio_ad%d", ad_num);
//	TH2F* h_delayed_energy_vs_distance_norm_ratio=new TH2F(name,name,150,1.5,3.);
	TH2F* h_delayed_energy_vs_distance_norm_ratio=new TH2F(name,name,230,0.7,3.,11,0.,5.5);

	sprintf(name, "h_delayed_energy_vs_distance_norm_sub_ad%d", ad_num);
//	TH2F* h_delayed_energy_vs_distance_norm_sub=new TH2F(name,name,150,1.5,3.);
	TH2F* h_delayed_energy_vs_distance_norm_sub=new TH2F(name,name,230,0.7,3.,11,0.,5.5);


	sprintf(name, "h_prompt_energy_vs_distance_rate_ratio_ad%d", ad_num);
	TH2F* h_prompt_energy_vs_distance_rate_ratio=new TH2F(name,name,113,0.7,12.,11,0.,5.5);

	sprintf(name, "h_prompt_energy_vs_distance_rate_sub_ad%d", ad_num);
	TH2F* h_prompt_energy_vs_distance_rate_sub=new TH2F(name,name,113,0.7,12.,11,0.,5.5);

	sprintf(name, "h_prompt_energy_vs_distance_norm_ratio_ad%d", ad_num);
	TH2F* h_prompt_energy_vs_distance_norm_ratio=new TH2F(name,name,113,0.7,12.,11,0.,5.5);

	sprintf(name, "h_prompt_energy_vs_distance_norm_sub_ad%d", ad_num);
	TH2F* h_prompt_energy_vs_distance_norm_sub=new TH2F(name,name,113,0.7,12.,11,0.,5.5);

/*	sprintf(name, "h_prompt_energy_vs_time_DTnorm_sub_ad%d", ad_num);
	TH2F* h_prompt_energy_vs_time_DTnorm_sub=new TH2F(name,name,1999,1,2000,113,0.7,12.);

	sprintf(name, "h_prompt_energy_vs_time_DTnorm_DT800_sub_ad%d", ad_num);
	TH2F* h_prompt_energy_vs_time_DTnorm_DT800_sub=new TH2F(name,name,1999,1,2000,113,0.7,12.);*/


			TH1F* h_acc_Denergy_dist_norm[6]; //delayed energy histogram
			for(int dist=0; dist<6; dist++){
				sprintf(name, "h_acc_Denergy_dist_norm%d_ad%d", dist, ad_num);
			//	h_acc_Denergy_dist_norm[dist]=new TH1F(name,name,150,1.5,3.);
				h_acc_Denergy_dist_norm[dist]=new TH1F(name,name,230,0.7,3.);
			}

			TH1F* h_acc_Penergy_dist_norm[6]; //delayed energy histogram
			for(int dist=0; dist<6; dist++){
				sprintf(name, "h_acc_Penergy_dist_norm%d_ad%d",dist, ad_num);
			//	h_acc_Penergy_dist_norm[dist]=new TH1F(name,name,150,1.5,3.);
				h_acc_Penergy_dist_norm[dist]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_acc_Denergy_dist[6]; //delayed energy histogram
			for(int dist=0; dist<6; dist++){
				sprintf(name, "h_acc_Denergy_dist%d_ad%d", dist, ad_num);
			//	h_acc_Denergy_dist[dist]=new TH1F(name,name,150,1.5,3.);
				h_acc_Denergy_dist[dist]=new TH1F(name,name,230,0.7,3.);
			}

			TH1F* h_acc_Penergy_dist[6]; //delayed energy histogram
			for(int dist=0; dist<6; dist++){
				sprintf(name, "h_acc_Penergy_dist%d_ad%d",dist, ad_num);
			//	h_acc_Penergy_dist[dist]=new TH1F(name,name,150,1.5,3.);
				h_acc_Penergy_dist[dist]=new TH1F(name,name,113,0.7,12.);
			}

			TH1F* h_ibd_Denergy_dist[6]; //delayed energy histogram
			for(int dist=0; dist<6; dist++){
				sprintf(name, "h_ibd_Denergy_dist%d_ad%d", dist, ad_num);
			//	h_ibd_Denergy_dist[dist]=new TH1F(name,name,150,1.5,3.);
				h_ibd_Denergy_dist[dist]=new TH1F(name,name,230,0.7,3.);
			}

			TH1F* h_ibd_Penergy_dist[6]; //delayed energy histogram
			for(int dist=0; dist<6; dist++){
				sprintf(name, "h_ibd_Penergy_dist%d_ad%d",dist, ad_num);
			//	h_ibd_Penergy_dist[dist]=new TH1F(name,name,150,1.5,3.);
				h_ibd_Penergy_dist[dist]=new TH1F(name,name,113,0.7,12.);
			}

		TH1F* h_sub_delayed_energy_norm_dist[6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_sub_delayed_energy_norm_dist%d_ad%d", dist, ad_num);
		//	h_sub_delayed_energy_norm_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_sub_delayed_energy_norm_dist[dist]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_sub_prompt_energy_norm_dist[6]; //prompt energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_sub_prompt_energy_norm_dist%d_ad%d",dist, ad_num);
		//	h_sub_prompt_energy_norm_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_sub_prompt_energy_norm_dist[dist]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_sub_delayed_energy_scaled_dist[6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_sub_delayed_energy_scaled_dist%d_ad%d", dist, ad_num);
		//	h_sub_delayed_energy_scaled_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_sub_delayed_energy_scaled_dist[dist]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_sub_prompt_energy_scaled_dist[6]; //prompt energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_sub_prompt_energy_scaled_dist%d_ad%d",dist, ad_num);
		//	h_sub_prompt_energy_scaled_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_sub_prompt_energy_scaled_dist[dist]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_ratio_delayed_energy_norm_dist[6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_ratio_delayed_energy_norm_dist%d_ad%d", dist, ad_num);
		//	h_ratio_delayed_energy_norm_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_ratio_delayed_energy_norm_dist[dist]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_ratio_prompt_energy_norm_dist[6]; //prompt energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_ratio_prompt_energy_norm_dist%d_ad%d",dist, ad_num);
		//	h_ratio_prompt_energy_norm_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_ratio_prompt_energy_norm_dist[dist]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_ratio_delayed_energy_scaled_dist[6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_ratio_delayed_energy_scaled_dist%d_ad%d", dist, ad_num);
		//	h_ratio_delayed_energy_scaled_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_ratio_delayed_energy_scaled_dist[dist]=new TH1F(name,name,230,0.7,3.);
		}

		TH1F* h_ratio_prompt_energy_scaled_dist[6]; //prompt energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_ratio_prompt_energy_scaled_dist%d_ad%d",dist, ad_num);
		//	h_ratio_prompt_energy_scaled_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_ratio_prompt_energy_scaled_dist[dist]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_ratio_prompt_energy_scaledTOnorm_dist[6]; //prompt energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_ratio_prompt_energy_scaledTOnorm_dist%d_ad%d",dist, ad_num);
		//	h_ratio_prompt_energy_scaledTOnorm_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_ratio_prompt_energy_scaledTOnorm_dist[dist]=new TH1F(name,name,113,0.7,12.);
		}

		TH1F* h_ratio_delayed_energy_scaledTOnorm_dist[6]; //delayed energy histogram
		for(int dist=0; dist<6; dist++){
			sprintf(name, "h_ratio_delayed_energy_scaledTOnorm_dist%d_ad%d", dist, ad_num);
		//	h_ratio_delayed_energy_scaledTOnorm_dist[dist]=new TH1F(name,name,150,1.5,3.);
			h_ratio_delayed_energy_scaledTOnorm_dist[dist]=new TH1F(name,name,230,0.7,3.);
		}


		sprintf(name, "h_sub_delayed_energy_scaled_p35_ad%d", ad_num);
	//	h_sub_delayed_energy_scaled_p35=new TH1F(name,name,150,1.5,3.);
		TH1F* h_sub_delayed_energy_scaled_p35=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_sub_delayed_energy_norm_dist_2plus_ad%d", ad_num);
	//	h_sub_delayed_energy_norm_dist_2plus=new TH1F(name,name,150,1.5,3.);
		TH1F* h_sub_delayed_energy_norm_dist_2plus=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_sub_prompt_energy_norm_dist_2plus_ad%d", ad_num);
	//	h_sub_prompt_energy_norm_dist_2plus=new TH1F(name,name,150,1.5,3.);
		TH1F* h_sub_prompt_energy_norm_dist_2plus=new TH1F(name,name,113,0.7,12.);

		sprintf(name, "h_sub_delayed_energy_scaled_dist_2plus_ad%d", ad_num);
	//	h_sub_delayed_energy_scaled_dist_2plus=new TH1F(name,name,150,1.5,3.);
		TH1F* h_sub_delayed_energy_scaled_dist_2plus=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_sub_prompt_energy_scaled_dist_2plus_ad%d", ad_num);
	//	h_sub_prompt_energy_scaled_dist_2plus=new TH1F(name,name,150,1.5,3.);
		TH1F* h_sub_prompt_energy_scaled_dist_2plus=new TH1F(name,name,113,0.7,12.);

		sprintf(name, "h_ratio_delayed_energy_scaled_dist_2plus_ad%d", ad_num);
		TH1F* h_ratio_delayed_energy_scaled_dist_2plus=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_ratio_prompt_energy_scaled_dist_2plus_ad%d", ad_num);
		TH1F* h_ratio_prompt_energy_scaled_dist_2plus=new TH1F(name,name,113,0.7,12.);

		sprintf(name, "h_ratio_delayed_energy_norm_dist_2plus_ad%d", ad_num);
		TH1F* h_ratio_delayed_energy_norm_dist_2plus=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_ratio_prompt_energy_norm_dist_2plus_ad%d", ad_num);
		TH1F* h_ratio_prompt_energy_norm_dist_2plus=new TH1F(name,name,113,0.7,12.);

		sprintf(name, "h_ratio_delayed_energy_scaledTOnorm_dist_2plus_ad%d", ad_num);
		TH1F* h_ratio_delayed_energy_scaledTOnorm_dist_2plus=new TH1F(name,name,230,0.7,3.);

		sprintf(name, "h_ratio_prompt_energy_scaledTOnorm_dist_2plus_ad%d", ad_num);
		TH1F* h_ratio_prompt_energy_scaledTOnorm_dist_2plus=new TH1F(name,name,113,0.7,12.);

		sprintf(name,"h_sub_PvsD_energy_rate_ad%d",ad_num); //prompt vs. delayed energy histogram
	//	h_sub_PvsD_energy_rate=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
		TH2F* h_sub_PvsD_energy_rate=new TH2F(name,name,175,0.7,12.,230,0.7,3.);

		sprintf(name,"h_sub_PvsD_energy_rate_1m_ad%d",ad_num); //prompt vs. delayed energy histogram
	//	h_sub_PvsD_energy_rate_1m=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
		TH2F* h_sub_PvsD_energy_rate_1m=new TH2F(name,name,175,0.7,12.,230,0.7,3.);

		sprintf(name,"h_sub_PvsD_energy_norm_ad%d",ad_num); //prompt vs. delayed energy histogram
	//	h_sub_PvsD_energy_norm=new TH2F(name,name,175,0.7,12.,150,1.5,3.);
		TH2F* h_sub_PvsD_energy_norm=new TH2F(name,name,175,0.7,12.,230,0.7,3.);


	sprintf(name, "h_ratio_ad%d", ad_num);
	TH1F* h_ratio=new TH1F(name,name,700,0,7.); //ratio histogram

	sprintf(name, "h_coincidence_time_1us_ad%d", ad_num);
	TH1D* h_coincidence_time_1us=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_1m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_1m=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_075m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_075m=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_05m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_05m=new TH1D(name,name,1999,1,2000);

		sprintf(name, "h_time_efficiency_1us_ad%d", ad_num);
		TH1D* h_time_efficiency_1us=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_1m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_1m=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_075m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_075m=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_05m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_05m=new TH1D(name,name,100,1,2000);


	sprintf(name, "h_coincidence_time_1us_ls_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_ls=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_ls_1m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_ls_1m=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_ls_075m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_ls_075m=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_ls_05m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_ls_05m=new TH1D(name,name,1999,1,2000);


	sprintf(name, "h_coincidence_time_1us_gdls_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_gdls=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_gdls_1m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_gdls_1m=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_gdls_075m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_gdls_075m=new TH1D(name,name,1999,1,2000);

	sprintf(name, "h_coincidence_time_1us_gdls_05m_ad%d", ad_num);
	TH1D* h_coincidence_time_1us_gdls_05m=new TH1D(name,name,1999,1,2000);

		sprintf(name, "h_time_efficiency_1us_ls_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_ls=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_ls_1m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_ls_1m=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_ls_075m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_ls_075m=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_ls_05m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_ls_05m=new TH1D(name,name,100,1,2000);


		sprintf(name, "h_time_efficiency_1us_gdls_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_gdls=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_gdls_1m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_gdls_1m=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_gdls_075m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_gdls_075m=new TH1D(name,name,100,1,2000);

		sprintf(name, "h_time_efficiency_1us_gdls_05m_ad%d", ad_num);
		TH1D* h_time_efficiency_1us_gdls_05m=new TH1D(name,name,100,1,2000);


		sprintf(name,"h_total_sub_distVStime_norm_ad%d",ad_num);
		TH2F* h_total_sub_distVStime_norm=new TH2F(name,name,1999,1,2000,700,0,7.);

		sprintf(name,"h_total_sub_distVStime_rate_ad%d",ad_num);
		TH2F* h_total_sub_distVStime_rate=new TH2F(name,name,1999,1,2000,700,0,7.);

		sprintf(name,"h_total_sub_distVStime_Ep35_norm_ad%d",ad_num);
		TH2F* h_total_sub_distVStime_Ep35_norm=new TH2F(name,name,1999,1,2000,700,0,7.);

		sprintf(name,"h_total_sub_distVStime_Ep35_rate_ad%d",ad_num);
		TH2F* h_total_sub_distVStime_Ep35_rate=new TH2F(name,name,1999,1,2000,700,0,7.);

		sprintf(name,"h_sub_DT_norm_ad%d",ad_num);
		TH1D* h_sub_DT_norm=new TH1D(name,name,500,0,10);

		sprintf(name,"h_sub_DT_rate_ad%d",ad_num);
		TH1D* h_sub_DT_rate=new TH1D(name,name,500,0,10);

		sprintf(name,"h_sub_DT_DTnorm_ad%d",ad_num);
		TH1D* h_sub_DT_DTnorm=new TH1D(name,name,500,0,10);

			sprintf(name,"h_sub_DT_3sig_norm_ad%d",ad_num);
			TH1D* h_sub_DT_3sig_norm=new TH1D(name,name,500,0,10);

			sprintf(name,"h_sub_DT_3sig_rate_ad%d",ad_num);
			TH1D* h_sub_DT_3sig_rate=new TH1D(name,name,500,0,10);

			sprintf(name,"h_sub_DT_3sig_DTnorm_ad%d",ad_num);
			TH1D* h_sub_DT_3sig_DTnorm=new TH1D(name,name,500,0,10);

		sprintf(name,"h_sub_DT_Ep35_norm_ad%d",ad_num);
		TH1D* h_sub_DT_Ep35_norm=new TH1D(name,name,500,0,10);

		sprintf(name,"h_sub_DT_Ep35_rate_ad%d",ad_num);
		TH1D* h_sub_DT_Ep35_rate=new TH1D(name,name,500,0,10);

		sprintf(name,"h_sub_DT_Ep35_DTnorm_ad%d",ad_num);
		TH1D* h_sub_DT_Ep35_DTnorm=new TH1D(name,name,500,0,10);

		sprintf(name,"h_sub_DT_rate_mod_ad%d",ad_num);
		TH1D* h_sub_DT_rate_mod=new TH1D(name,name,500,0,10);

		sprintf(name,"h_sub_DT_Ep35_rate_mod_ad%d",ad_num);
		TH1D* h_sub_DT_Ep35_rate_mod=new TH1D(name,name,500,0,10);

	//Loading the files
	char IBDFileName[64];
	sprintf(IBDFileName, "../nH_files/TotaledPlots_EH%d_%d.root",hall_num,pd_window_microsec);
	TFile *IBDFile = new TFile(IBDFileName);

	char accFileName[64];
	sprintf(accFileName, "../nH_files/TotaledSingles_%d_EH%d.root",pd_window_microsec,hall_num);
	TFile *accFile = new TFile(accFileName);

	double accCounts=0;
	double IBDCounts=0;
	double icounts = 0;
	double acounts = 0;
	double rcounts = 0;
	double rateCounts = 0;
	double scale = 0;
	double rate_scale = 0;
	float IBDbinCounts = 0;
	float accBinCounts = 0;
		double IBD_Edelayed_counts = 0;
		double Acc_Edelayed_counts = 0;
		float dCounts_ibd = 0;
		float dCounts_ibd_total = 0;
		float dCounts_acc = 0;
		float dCounts_acc_total = 0;
		double IBD_Eprompt_counts = 0;
		double Acc_Eprompt_counts = 0;
		int startBin = 0;
		int endBin = 0;
	double ibdCOUNTS = 0;
	double accCOUNTS = 0;
	double ibdERROR = 0;
	double accERROR = 0;
	float ibd_EvsD_counts = 0;
	float acc_EvsD_rate_counts = 0;
	float acc_EvsD_norm_counts = 0;
	double scaledCOUNTS = 0;
	double rateCOUNTS = 0;
	double scaledERROR = 0;
	double rateERROR = 0;
	double accPdist_counts = 0;
	double accPdist_norm_counts = 0;
	double ibdPdist_counts = 0;
	double accDdist_counts = 0;
	double accDdist_norm_counts = 0;
	double ibdDdist_counts = 0;
	double pDistError = 0;
	double dDistError = 0;

	double modScale = 0;
	double modScale_Ep = 0;

	//Scale factors:
		sprintf(name, "h_rate_scale_ad%d", ad_num);
		TH1F *h_rScales = (TH1F*)accFile->Get(name);
		h_rateScale->Add(h_rScales);

		sprintf(name, "h_norm_scale_ad%d", ad_num);
		TH1F *h_nScales = (TH1F*)accFile->Get(name);
		h_normScale->Add(h_nScales);

		h_ratioScale->Divide(h_normScale, h_rateScale,1,1);

		double nScale = 0;
		double rScale = 0;

		for(int iBin = 1; iBin < nRuns+1; iBin++){
			rScale = h_rateScale->GetBinContent(iBin);
			nScale = h_normScale->GetBinContent(iBin);
			h_scale_diff->Fill(rScale-nScale);

		}



	//Distance plots:
		sprintf(name, "h_total_distance_before_ad%d",ad_num);
		TH1F *h_IBD_distance = (TH1F*)IBDFile->Get(name);
		h_IBDDistance->Add(h_IBD_distance);

		sprintf(name, "h_total_acc_distance_scaled_ad%d", ad_num);
		TH1F *h_acc_scaled = (TH1F*)accFile->Get(name);
		h_accDistance_rate->Add(h_acc_scaled);

		sprintf(name, "h_total_acc_distance_norm_ad%d", ad_num);
		TH1F *h_acc_norm = (TH1F*)accFile->Get(name);
		h_accDistance_norm->Add(h_acc_norm);

			sprintf(name, "h_total_distance_3sig_ad%d",ad_num);
			TH1F *h_IBD_distance_3sig = (TH1F*)IBDFile->Get(name);
			h_IBDDistance_3sig->Add(h_IBD_distance_3sig);

			sprintf(name, "h_total_acc_distance_3sig_scaled_ad%d", ad_num);
			TH1F *h_acc_scaled_3sig = (TH1F*)accFile->Get(name);
			h_accDistance_3sig_rate->Add(h_acc_scaled_3sig);

			sprintf(name, "h_total_acc_distance_3sig_norm_ad%d", ad_num);
			TH1F *h_acc_norm_3sig = (TH1F*)accFile->Get(name);
			h_accDistance_3sig_norm->Add(h_acc_norm_3sig);

			sprintf(name, "h_total_acc_distance_3sig_DTnorm_ad%d", ad_num);
			TH1F *h_acc_DTnorm_3sig = (TH1F*)accFile->Get(name);
			h_accDistance_3sig_DTnorm->Add(h_acc_DTnorm_3sig);

		for(int ibin=0; ibin<702; ibin++){
			icounts = 0;
			icounts = h_IBDDistance->GetBinContent(ibin);
			h_IBDDistance->SetBinError(ibin, sqrt(icounts));
		}

		h_rateSub_distance->Add(h_IBD_distance);
		h_rateSub_distance->Add(h_accDistance_rate,-1);
		h_normSub_distance->Add(h_IBD_distance);
		h_normSub_distance->Add(h_accDistance_norm,-1);

		h_rateSub_distance_3sig->Add(h_IBD_distance_3sig);
		h_rateSub_distance_3sig->Add(h_accDistance_3sig_rate,-1);
		h_normSub_distance_3sig->Add(h_IBD_distance_3sig);
		h_normSub_distance_3sig->Add(h_accDistance_3sig_norm,-1);
		h_DTnormSub_distance_3sig->Add(h_IBD_distance_3sig);
		h_DTnormSub_distance_3sig->Add(h_accDistance_3sig_norm,-1);

		h_ratio->Divide(h_rateSub_distance,h_normSub_distance,1,1);

		h_dist_ratio_rateTOnorm->Divide(h_acc_scaled, h_accDistance_norm,1,1);
		h_dist_sub_rateTOnorm->Add(h_acc_scaled,1);
		h_dist_sub_rateTOnorm->Add(h_accDistance_norm, -1);

		cout << "Rate corrected!!! Number of IBDs: \t" << h_rateSub_distance->Integral(0,101) << endl;
		cout << "Norm corrected!!! Number of IBDs: \t" << h_normSub_distance->Integral(0,101) << endl;


			//400 subset
			sprintf(name, "h_total_distance_before_400_ad%d",ad_num);
			TH1F *h_IBD_distance_400 = (TH1F*)IBDFile->Get(name);

			sprintf(name, "h_total_acc_distance_400_scaled_ad%d", ad_num);
			TH1F *h_acc_scaled_400 = (TH1F*)accFile->Get(name);

			sprintf(name, "h_total_acc_distance_400_norm_ad%d", ad_num);
			TH1F *h_acc_norm_400 = (TH1F*)accFile->Get(name);

			h_rateSub_distance_400->Add(h_IBD_distance_400);
			h_rateSub_distance_400->Add(h_acc_scaled_400,-1);
			h_normSub_distance_400->Add(h_IBD_distance_400);
			h_normSub_distance_400->Add(h_acc_norm_400,-1);

			//600 subset
			sprintf(name, "h_total_distance_before_600_ad%d",ad_num);
			TH1F *h_IBD_distance_600 = (TH1F*)IBDFile->Get(name);

			sprintf(name, "h_total_acc_distance_600_scaled_ad%d", ad_num);
			TH1F *h_acc_scaled_600 = (TH1F*)accFile->Get(name);

			sprintf(name, "h_total_acc_distance_600_norm_ad%d", ad_num);
			TH1F *h_acc_norm_600 = (TH1F*)accFile->Get(name);

			h_rateSub_distance_600->Add(h_IBD_distance_600);
			h_rateSub_distance_600->Add(h_acc_scaled_600,-1);
			h_normSub_distance_600->Add(h_IBD_distance_600);
			h_normSub_distance_600->Add(h_acc_norm_600,-1);

			//800 subset
			sprintf(name, "h_total_distance_before_800_ad%d",ad_num);
			TH1F *h_IBD_distance_800 = (TH1F*)IBDFile->Get(name);

			sprintf(name, "h_total_acc_distance_800_scaled_ad%d", ad_num);
			TH1F *h_acc_scaled_800 = (TH1F*)accFile->Get(name);

			sprintf(name, "h_total_acc_distance_800_norm_ad%d", ad_num);
			TH1F *h_acc_norm_800 = (TH1F*)accFile->Get(name);

			h_rateSub_distance_800->Add(h_IBD_distance_800);
			h_rateSub_distance_800->Add(h_acc_scaled_800,-1);
			h_normSub_distance_800->Add(h_IBD_distance_800);
			h_normSub_distance_800->Add(h_acc_norm_800,-1);

		//stages
		double factor = (h_acc_norm->Integral())/(h_acc_scaled->Integral())-1;
		double factor_400 = (h_acc_norm_400->Integral())/(h_acc_scaled_400->Integral())-1;
		double factor_600 = (h_acc_norm_600->Integral())/(h_acc_scaled_600->Integral())-1;
		double factor_800 = (h_acc_norm_800->Integral())/(h_acc_scaled_800->Integral())-1;

		cout << "FACTOR = " << factor << endl;
		cout << "FACTOR_400 = " << factor_400 << endl;
		cout << "FACTOR_600 = " << factor_600 << endl;
		cout << "FACTOR_800 = " << factor_800 << endl;

		for(int stage=0; stage<5; stage++){
			h_rate_modified[0][stage]->Add(h_IBD_distance);
			h_rate_modified[0][stage]->Add(h_acc_scaled,-1+(-1*stage*factor/4));
			h_rate_modified[1][stage]->Add(h_IBD_distance_400);
			h_rate_modified[1][stage]->Add(h_acc_scaled_400,-1+(-1*stage*factor_400/4));
			h_rate_modified[2][stage]->Add(h_IBD_distance_600);
			h_rate_modified[2][stage]->Add(h_acc_scaled_600,-1+(-1*stage*factor_600/4));
			h_rate_modified[3][stage]->Add(h_IBD_distance_800);
			h_rate_modified[3][stage]->Add(h_acc_scaled_800,-1+(-1*stage*factor_800/4));
		}


	//Dist vs time:
		sprintf(name, "h_total_ibd_distVStime_ad%d",ad_num);
		TH2F *h_IBD_distVStime = (TH2F*)IBDFile->Get(name);

		sprintf(name, "h_total_acc_distVStime_norm_ad%d",ad_num);
		TH2F *h_acc_norm_distVStime = (TH2F*)accFile->Get(name);

		sprintf(name, "h_total_acc_distVStime_rate_ad%d",ad_num);
		TH2F *h_acc_rate_distVStime = (TH2F*)accFile->Get(name);

		h_total_sub_distVStime_norm->Add(h_IBD_distVStime);
		h_total_sub_distVStime_norm->Add(h_acc_norm_distVStime,-1);

		h_total_sub_distVStime_rate->Add(h_IBD_distVStime);
		h_total_sub_distVStime_rate->Add(h_acc_rate_distVStime,-1);

			sprintf(name, "h_total_ibd_distVStime_Ep35_ad%d",ad_num);
			TH2F *h_IBD_distVStime_Ep35 = (TH2F*)IBDFile->Get(name);

			sprintf(name, "h_total_acc_distVStime_Ep35_norm_ad%d",ad_num);
			TH2F *h_acc_norm_distVStime_Ep35 = (TH2F*)accFile->Get(name);

			sprintf(name, "h_total_acc_distVStime_Ep35_rate_ad%d",ad_num);
			TH2F *h_acc_rate_distVStime_Ep35 = (TH2F*)accFile->Get(name);

			h_total_sub_distVStime_Ep35_norm->Add(h_IBD_distVStime_Ep35);
			h_total_sub_distVStime_Ep35_norm->Add(h_acc_norm_distVStime_Ep35,-1);

			h_total_sub_distVStime_Ep35_rate->Add(h_IBD_distVStime_Ep35);
			h_total_sub_distVStime_Ep35_rate->Add(h_acc_rate_distVStime_Ep35,-1);


	//DT plots:
		sprintf(name, "h_total_ibd_DT_ad%d",ad_num);
		TH1D *h_ibd_DT = (TH1D*)IBDFile->Get(name);

		sprintf(name, "h_total_acc_DT_norm_ad%d",ad_num);
		TH1D *h_acc_DT_norm = (TH1D*)accFile->Get(name);

		sprintf(name, "h_total_acc_DT_rate_ad%d",ad_num);
		TH1D *h_acc_DT_rate = (TH1D*)accFile->Get(name);

		sprintf(name, "h_total_acc_DT_DTnorm_ad%d",ad_num);
		TH1D *h_acc_DT_DTnorm = (TH1D*)accFile->Get(name);

			sprintf(name, "h_total_ibd_DT_3sig_ad%d",ad_num);
			TH1D *h_ibd_DT_3sig = (TH1D*)IBDFile->Get(name);

			sprintf(name, "h_total_acc_DT_3sig_rate_ad%d",ad_num);
			TH1D *h_acc_DT_3sig_rate = (TH1D*)accFile->Get(name);

			sprintf(name, "h_total_acc_DT_3sig_norm_ad%d",ad_num);
			TH1D *h_acc_DT_3sig_norm = (TH1D*)accFile->Get(name);

			sprintf(name, "h_total_acc_DT_3sig_DTnorm_ad%d",ad_num);
			TH1D *h_acc_DT_3sig_DTnorm = (TH1D*)accFile->Get(name);

		sprintf(name, "h_total_ibd_DT_Ep35_ad%d",ad_num);
		TH1D *h_ibd_DT_Ep35 = (TH1D*)IBDFile->Get(name);

		sprintf(name, "h_total_acc_DT_Ep35_norm_ad%d",ad_num);
		TH1D *h_acc_DT_Ep35_norm = (TH1D*)accFile->Get(name);

		sprintf(name, "h_total_acc_DT_Ep35_rate_ad%d",ad_num);
		TH1D *h_acc_DT_Ep35_rate = (TH1D*)accFile->Get(name);

		sprintf(name, "h_total_acc_DT_Ep35_DTnorm_ad%d",ad_num);
		TH1D *h_acc_DT_Ep35_DTnorm = (TH1D*)accFile->Get(name);

		h_sub_DT_norm->Add(h_ibd_DT);
		h_sub_DT_rate->Add(h_ibd_DT);
		h_sub_DT_DTnorm->Add(h_ibd_DT);

		h_sub_DT_norm->Add(h_acc_DT_norm,-1);
		h_sub_DT_rate->Add(h_acc_DT_rate,-1);
		h_sub_DT_DTnorm->Add(h_acc_DT_DTnorm,-1);

			h_sub_DT_3sig_norm->Add(h_ibd_DT_3sig);
			h_sub_DT_3sig_rate->Add(h_ibd_DT_3sig);
			h_sub_DT_3sig_DTnorm->Add(h_ibd_DT_3sig);

			h_sub_DT_3sig_norm->Add(h_acc_DT_3sig_norm,-1);
			h_sub_DT_3sig_rate->Add(h_acc_DT_3sig_rate,-1);
			h_sub_DT_3sig_DTnorm->Add(h_acc_DT_3sig_DTnorm,-1);

		h_sub_DT_Ep35_norm->Add(h_ibd_DT_Ep35);
		h_sub_DT_Ep35_rate->Add(h_ibd_DT_Ep35);
		h_sub_DT_Ep35_DTnorm->Add(h_ibd_DT_Ep35);

		h_sub_DT_Ep35_norm->Add(h_acc_DT_Ep35_norm,-1);
		h_sub_DT_Ep35_rate->Add(h_acc_DT_Ep35_rate,-1);
		h_sub_DT_Ep35_DTnorm->Add(h_acc_DT_Ep35_DTnorm,-1);

		cout << "N_DC:\t" << h_ibd_DT->Integral(h_ibd_DT->FindBin(3),h_ibd_DT->FindBin(7)) << "\t+/-\t" << sqrt(h_ibd_DT->Integral(h_ibd_DT->FindBin(3),h_ibd_DT->FindBin(7))) << endl;
		cout << "N_cor:\t" << h_sub_DT_rate->Integral(h_sub_DT_rate->FindBin(3),h_sub_DT_rate->FindBin(7)) << endl;

				modScale = (h_acc_DT_DTnorm->Integral())/(h_acc_DT_rate->Integral());
				modScale_Ep = (h_acc_DT_Ep35_DTnorm->Integral())/(h_acc_DT_Ep35_rate->Integral());
				h_sub_DT_rate_mod->Add(h_ibd_DT);
				h_sub_DT_rate_mod->Add(h_acc_DT_rate,-1.*modScale);
				h_sub_DT_Ep35_rate_mod->Add(h_ibd_DT_Ep35);
				h_sub_DT_Ep35_rate_mod->Add(h_acc_DT_Ep35_rate,-1.*modScale_Ep);
				cout << "Modifying scale for DT plots (rate->DTnorm):"<< endl;
				cout << "\tScale:\t" << modScale << endl << "\tEp Scale:\t" << modScale_Ep << endl;
				

		for(int iBin = 0; iBin <502; iBin++){
			h_ibd_DT->SetBinError(iBin,sqrt((h_ibd_DT->GetBinContent(iBin))));
			h_acc_DT_norm->SetBinError(iBin,1/10.*sqrt((h_acc_DT_norm->GetBinContent(iBin))));
			h_acc_DT_rate->SetBinError(iBin,1/10.*sqrt((h_acc_DT_rate->GetBinContent(iBin))));
			h_acc_DT_DTnorm->SetBinError(iBin,1/10.*sqrt((h_acc_DT_DTnorm->GetBinContent(iBin))));
			h_sub_DT_norm->SetBinError(iBin,sqrt((h_ibd_DT->GetBinContent(iBin))));
			h_sub_DT_rate->SetBinError(iBin,sqrt((h_ibd_DT->GetBinContent(iBin))));
			h_sub_DT_DTnorm->SetBinError(iBin,sqrt((h_ibd_DT->GetBinContent(iBin))));
			h_sub_DT_rate_mod->SetBinError(iBin,sqrt((h_ibd_DT->GetBinContent(iBin))));

			h_ibd_DT_Ep35->SetBinError(iBin,sqrt((h_ibd_DT_Ep35->GetBinContent(iBin))));
			h_acc_DT_Ep35_norm->SetBinError(iBin,1/10.*sqrt((h_acc_DT_Ep35_norm->GetBinContent(iBin))));
			h_acc_DT_Ep35_rate->SetBinError(iBin,1/10.*sqrt((h_acc_DT_Ep35_rate->GetBinContent(iBin))));
			h_acc_DT_Ep35_DTnorm->SetBinError(iBin,1/10.*sqrt((h_acc_DT_Ep35_DTnorm->GetBinContent(iBin))));
			h_sub_DT_Ep35_norm->SetBinError(iBin,sqrt((h_ibd_DT_Ep35->GetBinContent(iBin))));
			h_sub_DT_Ep35_rate->SetBinError(iBin,sqrt((h_ibd_DT_Ep35->GetBinContent(iBin))));
			h_sub_DT_Ep35_DTnorm->SetBinError(iBin,sqrt((h_ibd_DT_Ep35->GetBinContent(iBin))));
			h_sub_DT_Ep35_rate_mod->SetBinError(iBin,sqrt((h_ibd_DT_Ep35->GetBinContent(iBin))));
		}




	//Delayed Energy:
		sprintf(name, "h_total_delayed_energy_before_ad%d",ad_num);
		TH1F *h_IBD_Edelayed = (TH1F*)IBDFile->Get(name);
		sprintf(name, "h_total_delayed_energy_scaled_ad%d", ad_num);
		TH1F *h_acc_Edelayed_scaled = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_delayed_energy_norm_ad%d", ad_num);
		TH1F *h_acc_Edelayed_norm = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_delayed_energy_DTnorm_ad%d", ad_num);
		TH1F *h_acc_Edelayed_DTnorm = (TH1F*)accFile->Get(name);

		h_Edelayed_IBD->Add(h_IBD_Edelayed);
		h_Edelayed_sub->Add(h_IBD_Edelayed);
		h_Edelayed_sub->Add(h_acc_Edelayed_scaled, -1);
		h_Edelayed_sub_norm->Add(h_IBD_Edelayed);
		h_Edelayed_sub_norm->Add(h_acc_Edelayed_norm, -1);
		h_Edelayed_sub_DTnorm->Add(h_IBD_Edelayed);
		h_Edelayed_sub_DTnorm->Add(h_acc_Edelayed_DTnorm, -1);

		sprintf(name, "h_total_delayed_energy_DT800_ad%d",ad_num);
		TH1F *h_IBD_Edelayed_DT800 = (TH1F*)IBDFile->Get(name);
		sprintf(name, "h_total_delayed_energy_DT800_scaled_ad%d", ad_num);
		TH1F *h_acc_Edelayed_DT800_scaled = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_delayed_energy_DT800_norm_ad%d", ad_num);
		TH1F *h_acc_Edelayed_DT800_norm = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_delayed_energy_DT800_DTnorm_ad%d", ad_num);
		TH1F *h_acc_Edelayed_DT800_DTnorm = (TH1F*)accFile->Get(name);

		h_Edelayed_IBD_DT800->Add(h_IBD_Edelayed_DT800);
		h_Edelayed_sub_DT800->Add(h_IBD_Edelayed_DT800);
		h_Edelayed_sub_DT800->Add(h_acc_Edelayed_DT800_scaled, -(1.0 +0.00));
		h_Edelayed_sub_DT800_norm->Add(h_IBD_Edelayed_DT800);
		h_Edelayed_sub_DT800_norm->Add(h_acc_Edelayed_DT800_norm, -(1.0 +0.00));
		h_Edelayed_sub_DT800_DTnorm->Add(h_IBD_Edelayed_DT800);
		h_Edelayed_sub_DT800_DTnorm->Add(h_acc_Edelayed_DT800_DTnorm, -1);
cout << "1" << endl;
/*			sprintf(name, "h_total_delayed_energy_fine_before_ad%d",ad_num);
			TH1F *h_IBD_Edelayed_fine = (TH1F*)IBDFile->Get(name);
				h_IBD_Edelayed_fine->Rebin(10);
cout << "1.2" << endl;
			sprintf(name, "h_total_delayed_energy_fine_scaled_1_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_scaled = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_scaled->Rebin(10);
cout << "1.4" << endl;
			sprintf(name, "h_total_delayed_energy_fine_norm_1_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_norm = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_norm->Rebin(10);
cout << "1.6" << endl;
			sprintf(name, "h_total_delayed_energy_fine_DTnorm_1_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_DTnorm = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_DTnorm->Rebin(10);
cout << "1.8" << endl;

			h_Edelayed_IBD_fine->Add(h_IBD_Edelayed_fine);
			h_Edelayed_sub_fine->Add(h_IBD_Edelayed_fine);
			h_Edelayed_sub_fine->Add(h_acc_Edelayed_fine_scaled, -(1.0));
			h_Edelayed_sub_fine_norm->Add(h_IBD_Edelayed_fine);
			h_Edelayed_sub_fine_norm->Add(h_acc_Edelayed_fine_norm, -(1.0 +0.00));
			h_Edelayed_sub_fine_DTnorm->Add(h_IBD_Edelayed_fine);
			h_Edelayed_sub_fine_DTnorm->Add(h_acc_Edelayed_fine_DTnorm, -1);
cout << "2" << endl;
		
			sprintf(name, "h_total_delayed_energy_fine_DT800_ad%d",ad_num);
			TH1F *h_IBD_Edelayed_fine_DT800 = (TH1F*)IBDFile->Get(name);
				h_IBD_Edelayed_fine_DT800->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_DT800_scaled_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_DT800_scaled = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_DT800_scaled->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_DT800_norm_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_DT800_norm = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_DT800_norm->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_DT800_DTnorm_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_DT800_DTnorm = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_DT800_DTnorm->Rebin(10);

			h_Edelayed_IBD_fine_DT800->Add(h_IBD_Edelayed_fine_DT800);
			h_Edelayed_sub_fine_DT800->Add(h_IBD_Edelayed_fine_DT800);
			h_Edelayed_sub_fine_DT800->Add(h_acc_Edelayed_fine_DT800_scaled, -(1.0));
			h_Edelayed_sub_fine_DT800_norm->Add(h_IBD_Edelayed_fine_DT800);
			h_Edelayed_sub_fine_DT800_norm->Add(h_acc_Edelayed_fine_DT800_norm, -(1.0 +0.00));
			h_Edelayed_sub_fine_DT800_DTnorm->Add(h_IBD_Edelayed_fine_DT800);
			h_Edelayed_sub_fine_DT800_DTnorm->Add(h_acc_Edelayed_fine_DT800_DTnorm, -1);

cout << "3" << endl;

			sprintf(name, "h_total_delayed_energy_fine_Ep35_ad%d",ad_num);
			TH1F *h_IBD_Edelayed_fine_Ep35 = (TH1F*)IBDFile->Get(name);
				h_IBD_Edelayed_fine_Ep35->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_Ep35_scaled_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_Ep35_scaled = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_Ep35_scaled->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_Ep35_norm_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_Ep35_norm = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_Ep35_norm->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_Ep35_DTnorm_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_Ep35_DTnorm = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_Ep35_DTnorm->Rebin(10);

			h_Edelayed_IBD_fine_Ep35->Add(h_IBD_Edelayed_fine_DT800);
			h_Edelayed_sub_fine_Ep35->Add(h_IBD_Edelayed_fine_Ep35);
			h_Edelayed_sub_fine_Ep35->Add(h_acc_Edelayed_fine_Ep35_scaled, -(1.0 +0.00));
			h_Edelayed_sub_fine_Ep35_norm->Add(h_IBD_Edelayed_fine_Ep35);
			h_Edelayed_sub_fine_Ep35_norm->Add(h_acc_Edelayed_fine_Ep35_norm, -(1.0 +0.00));
			h_Edelayed_sub_fine_Ep35_DTnorm->Add(h_IBD_Edelayed_fine_Ep35);
			h_Edelayed_sub_fine_Ep35_DTnorm->Add(h_acc_Edelayed_fine_Ep35_DTnorm, -1);
cout << "4" << endl;
			sprintf(name, "h_total_delayed_energy_fine_DT800_Ep35_ad%d",ad_num);
			TH1F *h_IBD_Edelayed_fine_DT800_Ep35 = (TH1F*)IBDFile->Get(name);
				h_IBD_Edelayed_fine_DT800_Ep35->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_DT800_Ep35_scaled_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_DT800_Ep35_scaled = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_DT800_Ep35_scaled->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_DT800_Ep35_norm_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_DT800_Ep35_norm = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_DT800_Ep35_norm->Rebin(10);
			sprintf(name, "h_total_delayed_energy_fine_DT800_Ep35_DTnorm_ad%d", ad_num);
			TH1F *h_acc_Edelayed_fine_DT800_Ep35_DTnorm = (TH1F*)accFile->Get(name);
				h_acc_Edelayed_fine_DT800_Ep35_DTnorm->Rebin(10);

			h_Edelayed_IBD_fine_DT800_Ep35->Add(h_IBD_Edelayed_fine_DT800);
			h_Edelayed_sub_fine_DT800_Ep35->Add(h_IBD_Edelayed_fine_DT800_Ep35);
			h_Edelayed_sub_fine_DT800_Ep35->Add(h_acc_Edelayed_fine_DT800_Ep35_scaled, -(1.0 +0.00));
			h_Edelayed_sub_fine_DT800_Ep35_norm->Add(h_IBD_Edelayed_fine_DT800_Ep35);
			h_Edelayed_sub_fine_DT800_Ep35_norm->Add(h_acc_Edelayed_fine_DT800_Ep35_norm, -(1.0 +0.00));
			h_Edelayed_sub_fine_DT800_Ep35_DTnorm->Add(h_IBD_Edelayed_fine_DT800_Ep35);
			h_Edelayed_sub_fine_DT800_Ep35_DTnorm->Add(h_acc_Edelayed_fine_DT800_Ep35_DTnorm, -1);*/
cout << "5" << endl;

			for(int iz = 0; iz < NzBins; iz++){
				TH1F *h_IBD_Edelayed_z = (TH1F*)IBDFile->Get(Form("h_total_delayed_energy_before_z_ad%d_iz%d", ad_num, iz+1));
				TH1F *h_acc_Edelayed_scaled_z = (TH1F*)accFile->Get(Form("h_total_delayed_energy_scaled_z_ad%d_iz%d", ad_num, iz+1));
				TH1F *h_acc_Edelayed_norm_z = (TH1F*)accFile->Get(Form("h_total_delayed_energy_norm_z_ad%d_iz%d", ad_num, iz+1));
				TH1F *h_acc_Edelayed_DTnorm_z = (TH1F*)accFile->Get(Form("h_total_delayed_energy_DTnorm_z_ad%d_iz%d", ad_num, iz+1));

				h_Edelayed_ibd_z[iz]->Add(h_IBD_Edelayed_z);
				h_Edelayed_sub_rate_z[iz]->Add(h_IBD_Edelayed_z);
				h_Edelayed_sub_rate_z[iz]->Add(h_acc_Edelayed_scaled_z,-1);
				h_Edelayed_sub_norm_z[iz]->Add(h_IBD_Edelayed_z);
				h_Edelayed_sub_norm_z[iz]->Add(h_acc_Edelayed_norm_z,-1);
				h_Edelayed_sub_DTnorm_z[iz]->Add(h_IBD_Edelayed_z);
				h_Edelayed_sub_DTnorm_z[iz]->Add(h_acc_Edelayed_DTnorm_z,-1);


				TH1F *h_IBD_Edelayed_DT800_z = (TH1F*)IBDFile->Get(Form("h_total_delayed_energy_DT800_z_ad%d_iz%d", ad_num, iz+1));
				TH1F *h_acc_Edelayed_scaled_DT800_z = (TH1F*)accFile->Get(Form("h_total_delayed_energy_scaled_DT800_z_ad%d_iz%d", ad_num, iz+1));
				TH1F *h_acc_Edelayed_norm_DT800_z = (TH1F*)accFile->Get(Form("h_total_delayed_energy_norm_DT800_z_ad%d_iz%d", ad_num, iz+1));
				TH1F *h_acc_Edelayed_DTnorm_DT800_z = (TH1F*)accFile->Get(Form("h_total_delayed_energy_DTnorm_DT800_z_ad%d_iz%d", ad_num, iz+1));

				h_Edelayed_ibd_DT800_z[iz]->Add(h_IBD_Edelayed_DT800_z);
				h_Edelayed_sub_rate_DT800_z[iz]->Add(h_IBD_Edelayed_DT800_z);
				h_Edelayed_sub_rate_DT800_z[iz]->Add(h_acc_Edelayed_scaled_DT800_z,-1);
				h_Edelayed_sub_norm_DT800_z[iz]->Add(h_IBD_Edelayed_DT800_z);
				h_Edelayed_sub_norm_DT800_z[iz]->Add(h_acc_Edelayed_norm_DT800_z,-1);
				h_Edelayed_sub_DTnorm_DT800_z[iz]->Add(h_IBD_Edelayed_DT800_z);
				h_Edelayed_sub_DTnorm_DT800_z[iz]->Add(h_acc_Edelayed_DTnorm_DT800_z,-1);
			}

			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				TH1F *h_IBD_Edelayed_r2 = (TH1F*)IBDFile->Get(Form("h_total_delayed_energy_before_r2_ad%d_ir2%d", ad_num, ir2+1));
				TH1F *h_acc_Edelayed_scaled_r2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_scaled_r2_ad%d_ir2%d", ad_num, ir2+1));
				TH1F *h_acc_Edelayed_norm_r2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_norm_r2_ad%d_ir2%d", ad_num, ir2+1));
				TH1F *h_acc_Edelayed_DTnorm_r2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_DTnorm_r2_ad%d_ir2%d", ad_num, ir2+1));

				h_Edelayed_ibd_r2[ir2]->Add(h_IBD_Edelayed_r2);
				h_Edelayed_sub_rate_r2[ir2]->Add(h_IBD_Edelayed_r2);
				h_Edelayed_sub_rate_r2[ir2]->Add(h_acc_Edelayed_scaled_r2,-(1.0 +0.00));
				h_Edelayed_sub_norm_r2[ir2]->Add(h_IBD_Edelayed_r2);
				h_Edelayed_sub_norm_r2[ir2]->Add(h_acc_Edelayed_norm_r2,-(1.0 +0.00));
				h_Edelayed_sub_DTnorm_r2[ir2]->Add(h_IBD_Edelayed_r2);
				h_Edelayed_sub_DTnorm_r2[ir2]->Add(h_acc_Edelayed_DTnorm_r2,-(1.0 +0.00));


				TH1F *h_IBD_Edelayed_DT800_r2 = (TH1F*)IBDFile->Get(Form("h_total_delayed_energy_DT800_r2_ad%d_ir2%d", ad_num, ir2+1));
				TH1F *h_acc_Edelayed_scaled_DT800_r2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_scaled_DT800_r2_ad%d_ir2%d", ad_num, ir2+1));
				TH1F *h_acc_Edelayed_norm_DT800_r2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_norm_DT800_r2_ad%d_ir2%d", ad_num, ir2+1));
				TH1F *h_acc_Edelayed_DTnorm_DT800_r2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_DTnorm_DT800_r2_ad%d_ir2%d", ad_num, ir2+1));

				h_Edelayed_ibd_DT800_r2[ir2]->Add(h_IBD_Edelayed_DT800_r2);
				h_Edelayed_sub_rate_DT800_r2[ir2]->Add(h_IBD_Edelayed_DT800_r2);
				h_Edelayed_sub_rate_DT800_r2[ir2]->Add(h_acc_Edelayed_scaled_DT800_r2,-(1.0 +0.00));
				h_Edelayed_sub_norm_DT800_r2[ir2]->Add(h_IBD_Edelayed_DT800_r2);
				h_Edelayed_sub_norm_DT800_r2[ir2]->Add(h_acc_Edelayed_norm_DT800_r2,-(1.0 +0.00));
				h_Edelayed_sub_DTnorm_DT800_r2[ir2]->Add(h_IBD_Edelayed_DT800_r2);
				h_Edelayed_sub_DTnorm_DT800_r2[ir2]->Add(h_acc_Edelayed_DTnorm_DT800_r2,-(1.0 +0.00));
			}

			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					TH1F *h_IBD_Edelayed_zVSr2 = (TH1F*)IBDFile->Get(Form("h_total_delayed_energy_before_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1));
					TH1F *h_acc_Edelayed_scaled_zVSr2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_scaled_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1));
					TH1F *h_acc_Edelayed_norm_zVSr2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_norm_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1));
					TH1F *h_acc_Edelayed_DTnorm_zVSr2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_DTnorm_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1));

					h_Edelayed_ibd_zVSr2[ir2][iz]->Add(h_IBD_Edelayed_zVSr2);
					h_Edelayed_sub_rate_zVSr2[ir2][iz]->Add(h_IBD_Edelayed_zVSr2);
					h_Edelayed_sub_rate_zVSr2[ir2][iz]->Add(h_acc_Edelayed_scaled_zVSr2,-(1.0 +0.00));
					h_Edelayed_sub_norm_zVSr2[ir2][iz]->Add(h_IBD_Edelayed_zVSr2);
					h_Edelayed_sub_norm_zVSr2[ir2][iz]->Add(h_acc_Edelayed_norm_zVSr2,-(1.0 +0.00));
					h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]->Add(h_IBD_Edelayed_zVSr2);
					h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]->Add(h_acc_Edelayed_DTnorm_zVSr2,-(1.0 +0.00));


					TH1F *h_IBD_Edelayed_DT800_zVSr2 = (TH1F*)IBDFile->Get(Form("h_total_delayed_energy_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1));
					TH1F *h_acc_Edelayed_scaled_DT800_zVSr2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_scaled_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1));
					TH1F *h_acc_Edelayed_norm_DT800_zVSr2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_norm_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1));
					TH1F *h_acc_Edelayed_DTnorm_DT800_zVSr2 = (TH1F*)accFile->Get(Form("h_total_delayed_energy_DTnorm_DT800_zVSr2_ad%d_ir2%d_iz%d", ad_num, ir2+1, iz+1));

					h_Edelayed_ibd_DT800_zVSr2[ir2][iz]->Add(h_IBD_Edelayed_DT800_zVSr2);
					h_Edelayed_sub_rate_DT800_zVSr2[ir2][iz]->Add(h_IBD_Edelayed_DT800_zVSr2);
					h_Edelayed_sub_rate_DT800_zVSr2[ir2][iz]->Add(h_acc_Edelayed_scaled_DT800_zVSr2,-(1.0 +0.00));
					h_Edelayed_sub_norm_DT800_zVSr2[ir2][iz]->Add(h_IBD_Edelayed_DT800_zVSr2);
					h_Edelayed_sub_norm_DT800_zVSr2[ir2][iz]->Add(h_acc_Edelayed_norm_DT800_zVSr2,-(1.0 +0.00));
					h_Edelayed_sub_DTnorm_DT800_zVSr2[ir2][iz]->Add(h_IBD_Edelayed_DT800_zVSr2);
					h_Edelayed_sub_DTnorm_DT800_zVSr2[ir2][iz]->Add(h_acc_Edelayed_DTnorm_DT800_zVSr2,-(1.0 +0.00));
				}
			}


	//Prompt Energy:
		sprintf(name, "h_total_prompt_energy_before_ad%d",ad_num);
		TH1F *h_IBD_Eprompt = (TH1F*)IBDFile->Get(name);
		sprintf(name, "h_total_prompt_energy_scaled_ad%d", ad_num);
		TH1F *h_acc_Eprompt_scaled = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_prompt_energy_norm_ad%d", ad_num);
		TH1F *h_acc_Eprompt_norm = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_prompt_energy_DTnorm_ad%d", ad_num);
		TH1F *h_acc_Eprompt_DTnorm = (TH1F*)accFile->Get(name);

		h_Eprompt_sub->Add(h_IBD_Eprompt);
		h_Eprompt_sub->Add(h_acc_Eprompt_scaled, -1);
		h_Eprompt_sub_norm->Add(h_IBD_Eprompt);
		h_Eprompt_sub_norm->Add(h_acc_Eprompt_norm, -1);
		h_Eprompt_sub_DTnorm->Add(h_IBD_Eprompt);
		h_Eprompt_sub_DTnorm->Add(h_acc_Eprompt_DTnorm, -1);

		sprintf(name, "h_total_prompt_energy_DT800_ad%d",ad_num);
		TH1F *h_IBD_Eprompt_DT800 = (TH1F*)IBDFile->Get(name);
		sprintf(name, "h_total_prompt_energy_DT800_scaled_ad%d", ad_num);
		TH1F *h_acc_Eprompt_DT800_scaled = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_prompt_energy_DT800_norm_ad%d", ad_num);
		TH1F *h_acc_Eprompt_DT800_norm = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_prompt_energy_DT800_DTnorm_ad%d", ad_num);
		TH1F *h_acc_Eprompt_DT800_DTnorm = (TH1F*)accFile->Get(name);

		h_Eprompt_sub_DT800->Add(h_IBD_Eprompt_DT800);
		h_Eprompt_sub_DT800->Add(h_acc_Eprompt_DT800_scaled, -1);
		h_Eprompt_sub_DT800_norm->Add(h_IBD_Eprompt_DT800);
		h_Eprompt_sub_DT800_norm->Add(h_acc_Eprompt_DT800_norm, -1);
		h_Eprompt_sub_DT800_DTnorm->Add(h_IBD_Eprompt_DT800);
		h_Eprompt_sub_DT800_DTnorm->Add(h_acc_Eprompt_DT800_DTnorm, -1);

			sprintf(name, "h_total_prompt_energy_DT800_3sig_ad%d",ad_num);
			TH1F *h_IBD_Eprompt_DT800_3sig = (TH1F*)IBDFile->Get(name);
			sprintf(name, "h_total_prompt_energy_DT800_3sig_scaled_ad%d", ad_num);
			TH1F *h_acc_Eprompt_DT800_3sig_scaled = (TH1F*)accFile->Get(name);
			
			h_Eprompt_sub_DT800_3sig->Add(h_IBD_Eprompt_DT800_3sig);
			h_Eprompt_sub_DT800_3sig->Add(h_acc_Eprompt_DT800_3sig_scaled, -1);

cout << "delayed and prompt done" << endl;

	//Distance Slices:
		for(int dist=0; dist < 6; dist++){
			sprintf(name, "h_total_prompt_energy_dist%d_ad%d",dist,ad_num);
			h_ibd_Penergy_dist[dist] = (TH1F*)IBDFile->Get(name);
			sprintf(name, "h_total_prompt_energy_scaled_dist%d_ad%d",dist,ad_num);
			h_acc_Penergy_dist[dist] = (TH1F*)accFile->Get(name);
			sprintf(name, "h_total_prompt_energy_norm_dist%d_ad%d",dist,ad_num);
			h_acc_Penergy_dist_norm[dist] = (TH1F*)accFile->Get(name);


			h_sub_prompt_energy_scaled_dist[dist]->Add(h_ibd_Penergy_dist[dist]);
			h_sub_prompt_energy_scaled_dist[dist]->Add(h_acc_Penergy_dist[dist],-1);
			h_sub_prompt_energy_norm_dist[dist]->Add(h_ibd_Penergy_dist[dist]);
			h_sub_prompt_energy_norm_dist[dist]->Add(h_acc_Penergy_dist_norm[dist],-1);

			for(int ibin = 0; ibin < 115; ibin++){
				ibdPdist_counts = h_ibd_Penergy_dist[dist]->GetBinContent(ibin);
				accPdist_counts = h_acc_Penergy_dist[dist]->GetBinContent(ibin);
				pDistError = sqrt(ibdPdist_counts + accPdist_counts);
				h_sub_prompt_energy_scaled_dist[dist]->SetBinError(ibin, pDistError);
				h_sub_prompt_energy_norm_dist[dist]->SetBinError(ibin, pDistError);
			}

			h_ratio_prompt_energy_scaled_dist[dist]->Divide(h_ibd_Penergy_dist[dist],h_acc_Penergy_dist[dist],1,1);
			h_ratio_prompt_energy_norm_dist[dist]->Divide(h_ibd_Penergy_dist[dist],h_acc_Penergy_dist_norm[dist],1,1);
			h_ratio_prompt_energy_scaledTOnorm_dist[dist]->Divide(h_acc_Penergy_dist[dist],h_acc_Penergy_dist_norm[dist],1,1);

			sprintf(name, "h_total_delayed_energy_dist%d_ad%d",dist,ad_num);
			h_ibd_Denergy_dist[dist] = (TH1F*)IBDFile->Get(name);
			sprintf(name, "h_total_delayed_energy_scaled_dist%d_ad%d",dist,ad_num);
			h_acc_Denergy_dist[dist] = (TH1F*)accFile->Get(name);
			sprintf(name, "h_total_delayed_energy_norm_dist%d_ad%d",dist,ad_num);
			h_acc_Denergy_dist_norm[dist] = (TH1F*)accFile->Get(name);

			h_sub_delayed_energy_scaled_dist[dist]->Add(h_ibd_Denergy_dist[dist]);
			h_sub_delayed_energy_scaled_dist[dist]->Add(h_acc_Denergy_dist[dist],-1);
			h_sub_delayed_energy_norm_dist[dist]->Add(h_ibd_Denergy_dist[dist]);
			h_sub_delayed_energy_norm_dist[dist]->Add(h_acc_Denergy_dist_norm[dist],-1);

			for(int ibin = 0; ibin < 232; ibin++){
				ibdDdist_counts = h_ibd_Denergy_dist[dist]->GetBinContent(ibin);
				accDdist_counts = h_acc_Denergy_dist[dist]->GetBinContent(ibin);
				dDistError = sqrt(ibdDdist_counts + accDdist_counts);
				h_sub_delayed_energy_scaled_dist[dist]->SetBinError(ibin, dDistError);
				h_sub_delayed_energy_norm_dist[dist]->SetBinError(ibin, dDistError);
			}

			h_ratio_delayed_energy_scaled_dist[dist]->Divide(h_ibd_Denergy_dist[dist],h_acc_Denergy_dist[dist],1,1);
			h_ratio_delayed_energy_norm_dist[dist]->Divide(h_ibd_Denergy_dist[dist],h_acc_Denergy_dist_norm[dist],1,1);
			h_ratio_delayed_energy_scaledTOnorm_dist[dist]->Divide(h_acc_Denergy_dist[dist],h_acc_Denergy_dist_norm[dist],1,1);

		}

cout << "delayed slices done" << endl;

		//Delayed Energy where distance > 2m
			sprintf(name, "h_total_delayed_energy_dist_2plus_ad%d",ad_num);
			TH1F *h_ibd_Denergy_2plus = (TH1F*)IBDFile->Get(name);
			sprintf(name, "h_total_delayed_energy_scaled_dist_2plus_ad%d", ad_num);
			TH1F *h_acc_Denergy_2plus = (TH1F*)accFile->Get(name);
			sprintf(name, "h_total_delayed_energy_norm_dist_2plus_ad%d", ad_num);
			TH1F *h_acc_Denergy_2plus_norm = (TH1F*)accFile->Get(name);

			h_sub_delayed_energy_scaled_dist_2plus->Add(h_ibd_Denergy_2plus);
			h_sub_delayed_energy_scaled_dist_2plus->Add(h_acc_Denergy_2plus,-1);
			h_sub_delayed_energy_norm_dist_2plus->Add(h_ibd_Denergy_2plus);
			h_sub_delayed_energy_norm_dist_2plus->Add(h_acc_Denergy_2plus_norm,-1);

			h_ratio_delayed_energy_scaled_dist_2plus->Divide(h_ibd_Denergy_2plus,h_acc_Denergy_2plus,1,1);
			h_ratio_delayed_energy_norm_dist_2plus->Divide(h_ibd_Denergy_2plus,h_acc_Denergy_2plus_norm,1,1);
			h_ratio_delayed_energy_scaledTOnorm_dist_2plus->Divide(h_acc_Denergy_2plus,h_acc_Denergy_2plus_norm,1,1);

cout << "delayed >2m done" << endl;

		//Delayed energy where prompt energy > 3.5 MeV
			sprintf(name, "h_total_delayed_energy_p35_ad%d",ad_num);
			TH1F *h_ibd_Denergy_p35 = (TH1F*)IBDFile->Get(name);
			sprintf(name, "h_total_delayed_energy_scaled_p35_ad%d", ad_num);
			TH1F *h_acc_Denergy_p35 = (TH1F*)accFile->Get(name);
			h_sub_delayed_energy_scaled_p35->Add(h_ibd_Denergy_p35);
			h_sub_delayed_energy_scaled_p35->Add(h_acc_Denergy_p35,-1);

cout << "special delayed done" << endl;

	//Prompt Energy where distance > 2m
		sprintf(name, "h_total_prompt_energy_dist_2plus_ad%d",ad_num);
		TH1F *h_ibd_Penergy_2plus = (TH1F*)IBDFile->Get(name);
		sprintf(name, "h_total_prompt_energy_scaled_dist_2plus_ad%d", ad_num);
		TH1F *h_acc_Penergy_2plus = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_prompt_energy_norm_dist_2plus_ad%d", ad_num);
		TH1F *h_acc_Penergy_2plus_norm = (TH1F*)accFile->Get(name);
		//For scaled sub
		h_sub_prompt_energy_scaled_dist_2plus->Add(h_ibd_Penergy_2plus);
		h_sub_prompt_energy_scaled_dist_2plus->Add(h_acc_Penergy_2plus,-1);
		//For normalied sub
		h_sub_prompt_energy_norm_dist_2plus->Add(h_ibd_Penergy_2plus);
		h_sub_prompt_energy_norm_dist_2plus->Add(h_acc_Penergy_2plus_norm,-1);

		h_ratio_prompt_energy_scaled_dist_2plus->Divide(h_ibd_Penergy_2plus,h_acc_Penergy_2plus,1,1);
		h_ratio_prompt_energy_norm_dist_2plus->Divide(h_ibd_Penergy_2plus,h_acc_Penergy_2plus_norm,1,1);
		h_ratio_prompt_energy_scaledTOnorm_dist_2plus->Divide(h_acc_Penergy_2plus,h_acc_Penergy_2plus_norm,1,1);



	//Vs Distance:
		sprintf(name, "h_ibd_delayed_energy_vs_distance_ad%d", ad_num);
		TH1F *h_IBD_Edelayed_vs_dist = (TH1F*)IBDFile->Get(name);
		sprintf(name, "h_total_delayed_energy_vs_distance_rate_ad%d", ad_num);
		TH1F *h_acc_Edelayed_vs_dist_rate = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_delayed_energy_vs_distance_norm_ad%d", ad_num);
		TH1F *h_acc_Edelayed_vs_dist_norm = (TH1F*)accFile->Get(name);

		sprintf(name, "h_total_prompt_energy_vs_distance_ad%d", ad_num);
		TH1F *h_IBD_Eprompt_vs_dist = (TH1F*)IBDFile->Get(name);
		sprintf(name, "h_total_prompt_energy_vs_distance_rate_ad%d", ad_num);
		TH1F *h_acc_Eprompt_vs_dist_rate = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_prompt_energy_vs_distance_norm_ad%d", ad_num);
		TH1F *h_acc_Eprompt_vs_dist_norm = (TH1F*)accFile->Get(name);

		h_delayed_energy_vs_distance_rate_ratio->Divide(h_IBD_Edelayed_vs_dist,h_acc_Edelayed_vs_dist_rate,1,1);
		h_delayed_energy_vs_distance_rate_sub->Add(h_IBD_Edelayed_vs_dist);
		h_delayed_energy_vs_distance_rate_sub->Add(h_acc_Edelayed_vs_dist_rate,-1);

		h_prompt_energy_vs_distance_rate_ratio->Divide(h_IBD_Eprompt_vs_dist,h_acc_Eprompt_vs_dist_rate,1,1);
		h_prompt_energy_vs_distance_rate_sub->Add(h_IBD_Eprompt_vs_dist);
		h_prompt_energy_vs_distance_rate_sub->Add(h_acc_Eprompt_vs_dist_rate,-1);

		h_delayed_energy_vs_distance_norm_ratio->Divide(h_IBD_Edelayed_vs_dist,h_acc_Edelayed_vs_dist_norm,1,1);
		h_delayed_energy_vs_distance_norm_sub->Add(h_IBD_Edelayed_vs_dist);
		h_delayed_energy_vs_distance_norm_sub->Add(h_acc_Edelayed_vs_dist_norm,-1);

		h_prompt_energy_vs_distance_norm_ratio->Divide(h_IBD_Eprompt_vs_dist,h_acc_Eprompt_vs_dist_norm,1,1);
		h_prompt_energy_vs_distance_norm_sub->Add(h_IBD_Eprompt_vs_dist);
		h_prompt_energy_vs_distance_norm_sub->Add(h_acc_Eprompt_vs_dist_norm,-1);


		sprintf(name, "h_total_ibd_promptVStime_ad%d", ad_num);
		TH2F *h_IBD_promptVStime = (TH2F*)IBDFile->Get(name);
		sprintf(name, "h_total_ibd_promptVStime_DT800_ad%d", ad_num);
		TH2F *h_IBD_promptVStime_DT800 = (TH2F*)IBDFile->Get(name);

		sprintf(name, "h_total_acc_promptVStime_DTnorm_ad%d", ad_num);
		TH1F *h_acc_promptVStime_DTnorm = (TH1F*)accFile->Get(name);
		sprintf(name, "h_total_acc_promptVStime_DT800_DTnorm_ad%d", ad_num);
		TH1F *h_acc_promptVStime_DTnorm_DT800 = (TH1F*)accFile->Get(name);

/*h_prompt_energy_vs_time_DTnorm_sub->Add(h_IBD_promptVStime);
h_prompt_energy_vs_time_DTnorm_sub->Add(h_acc_promptVStime_DTnorm,-1);

h_prompt_energy_vs_time_DTnorm_DT800_sub->Add(h_IBD_promptVStime_DT800);
h_prompt_energy_vs_time_DTnorm_DT800_sub->Add(h_acc_promptVStime_DTnorm_DT800,-1);*/



	//PvsD 2D plots
		sprintf(name, "h_total_ibd_energy_before_ad%d",ad_num);
		TH2F *h_IBD_energy_2D = (TH2F*)IBDFile->Get(name);
		sprintf(name, "h_total_acc_energy_before_ad%d", ad_num);
		TH2F *h_acc_energy_2D = (TH2F*)accFile->Get(name);
		h_sub_PvsD_energy_rate->Add(h_IBD_energy_2D);
		h_sub_PvsD_energy_rate->Add(h_acc_energy_2D,-1);

		sprintf(name, "h_total_ibd_energy_1m_ad%d",ad_num);
		TH2F *h_IBD_energy_2D_1m = (TH2F*)IBDFile->Get(name);
		sprintf(name, "h_total_acc_energy_1m_ad%d", ad_num);
		TH2F *h_acc_energy_2D_1m = (TH2F*)accFile->Get(name);
		h_sub_PvsD_energy_rate_1m->Add(h_IBD_energy_2D_1m);
		h_sub_PvsD_energy_rate_1m->Add(h_acc_energy_2D_1m,-1);

	double binVal = 0;
	double binErr = 0;

	//Coincidence Time
		sprintf(name, "h_total_delta_time_1us_ad%d",ad_num);
		TH1D *h_coin_time = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us->Add(h_coin_time);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time->GetBinContent(ibin);
			binErr = h_coin_time->GetBinError(ibin);
			h_coincidence_time_1us->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_1m_ad%d",ad_num);
		TH1D *h_coin_time1 = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_1m->Add(h_coin_time1);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time1->GetBinContent(ibin);
			binErr = h_coin_time1->GetBinError(ibin);
			h_coincidence_time_1us_1m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_1m->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_075m_ad%d",ad_num);
		TH1D *h_coin_time075 = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_075m->Add(h_coin_time075);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time075->GetBinContent(ibin);
			binErr = h_coin_time075->GetBinError(ibin);
			h_coincidence_time_1us_075m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_075m->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_05m_ad%d",ad_num);
		TH1D *h_coin_time05 = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_05m->Add(h_coin_time05);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time05->GetBinContent(ibin);
			binErr = h_coin_time05->GetBinError(ibin);
			h_coincidence_time_1us_05m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_05m->SetBinError(ibin,binErr);
		}


		sprintf(name, "h_total_delta_time_1us_ls_ad%d",ad_num);
		TH1D *h_coin_time_ls = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_ls->Add(h_coin_time_ls);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time_ls->GetBinContent(ibin);
			binErr = h_coin_time_ls->GetBinError(ibin);
			h_coincidence_time_1us_ls->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_ls->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_ls_1m_ad%d",ad_num);
		TH1D *h_coin_time1_ls = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_ls_1m->Add(h_coin_time1_ls);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time1_ls->GetBinContent(ibin);
			binErr = h_coin_time1_ls->GetBinError(ibin);
			h_coincidence_time_1us_ls_1m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_ls_1m->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_ls_075m_ad%d",ad_num);
		TH1D *h_coin_time075_ls = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_ls_075m->Add(h_coin_time075_ls);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time075_ls->GetBinContent(ibin);
			binErr = h_coin_time075_ls->GetBinError(ibin);
			h_coincidence_time_1us_ls_075m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_ls_075m->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_ls_05m_ad%d",ad_num);
		TH1D *h_coin_time05_ls = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_ls_05m->Add(h_coin_time05_ls);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time05_ls->GetBinContent(ibin);
			binErr = h_coin_time05_ls->GetBinError(ibin);
			h_coincidence_time_1us_ls_05m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_ls_05m->SetBinError(ibin,binErr);
		}


		sprintf(name, "h_total_delta_time_1us_gdls_ad%d",ad_num);
		TH1D *h_coin_time_gdls = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_gdls->Add(h_coin_time_gdls);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time_gdls->GetBinContent(ibin);
			binErr = h_coin_time_gdls->GetBinError(ibin);
			h_coincidence_time_1us_gdls->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_gdls->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_gdls_1m_ad%d",ad_num);
		TH1D *h_coin_time1_gdls = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_gdls_1m->Add(h_coin_time1_gdls);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time1_gdls->GetBinContent(ibin);
			binErr = h_coin_time1_gdls->GetBinError(ibin);
			h_coincidence_time_1us_gdls_1m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_gdls_1m->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_gdls_075m_ad%d",ad_num);
		TH1D *h_coin_time075_gdls = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_gdls_075m->Add(h_coin_time075_gdls);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time075_gdls->GetBinContent(ibin);
			binErr = h_coin_time075_gdls->GetBinError(ibin);
			h_coincidence_time_1us_gdls_075m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_gdls_075m->SetBinError(ibin,binErr);
		}

		sprintf(name, "h_total_delta_time_1us_gdls_05m_ad%d",ad_num);
		TH1D *h_coin_time05_gdls = (TH1D*)IBDFile->Get(name);
//		h_coincidence_time_1us_gdls_05m->Add(h_coin_time05_gdls);
		for(int ibin = 0; ibin < 2002; ibin++){
			binVal = 0;
			binErr = 0;
			binVal = h_coin_time05_gdls->GetBinContent(ibin);
			binErr = h_coin_time05_gdls->GetBinError(ibin);
			h_coincidence_time_1us_gdls_05m->SetBinContent(ibin,binVal);
			if(binErr == 0) continue;
			h_coincidence_time_1us_gdls_05m->SetBinError(ibin,binErr);
		}


		h_coincidence_time_1us->Rebin(20);
		h_coincidence_time_1us_1m->Rebin(20);
		h_coincidence_time_1us_075m->Rebin(20);
		h_coincidence_time_1us_05m->Rebin(20);

		h_coincidence_time_1us_ls->Rebin(20);
		h_coincidence_time_1us_ls_1m->Rebin(20);
		h_coincidence_time_1us_ls_075m->Rebin(20);
		h_coincidence_time_1us_ls_05m->Rebin(20);

		h_coincidence_time_1us_gdls->Rebin(20);
		h_coincidence_time_1us_gdls_1m->Rebin(20);
		h_coincidence_time_1us_gdls_075m->Rebin(20);
		h_coincidence_time_1us_gdls_05m->Rebin(20);

/*		double offset = 0;
		double offset_1 = 0;
		double offset_075 = 0;
		double offset_05 = 0;

		double offset_ls = 0;
		double offset_1_ls = 0;
		double offset_075_ls = 0;
		double offset_05_ls = 0;

		double offset_gdls = 0;
		double offset_1_gdls = 0;
		double offset_075_gdls = 0;
		double offset_05_gdls = 0;

		TF1* time_fit=new TF1("time_fit","[0]+[1]*TMath::Exp(-x/[2])+[3]*TMath::Exp(-x/[4])",10,2000);
			offset = h_acc_scaled->Integral(0,700)/700.;
	//		time_fit->FixParameter(0,offset);
			time_fit->SetParameter(0,offset);
			time_fit->SetParName(0,"constant");
			time_fit->SetParameter(1,3330);
			time_fit->SetParameter(2,210);
			time_fit->SetParName(2,"tau_ls");
			time_fit->SetParameter(3,3780);
			time_fit->SetParameter(4,29);
			time_fit->SetParName(4,"tau_gdls");

			time_fit->SetParLimits(1,0.,2.*h_coincidence_time_1us->GetBinContent(h_coincidence_time_1us->GetMaximumBin()));
			time_fit->SetParLimits(2,0.,400.);
			time_fit->SetParLimits(3,0.,2.*h_coincidence_time_1us->GetBinContent(h_coincidence_time_1us->GetMaximumBin()));
			time_fit->SetParLimits(4,0.,60.);

	//	cout << "Integral of time_fit: " << time_fit->Integral(0,100) << endl;

		TF1* time_fit1=new TF1("time_fit1","[0]+[1]*TMath::Exp(-x/[2])+[3]*TMath::Exp(-x/[4])",10,2000);
			for(int ibin=0; ibin < 700; ibin++){
				if(h_acc_scaled->GetXaxis()->GetBinCenter(ibin) > 1){
					endBin = ibin-1;
					break;
				}
			}
			offset_1 = h_acc_scaled->Integral(0,endBin)/endBin;
	//		time_fit1->FixParameter(0,offset_1);
	//		time_fit1->SetParameter(0,2680);
			time_fit1->SetParameter(0,offset_1);
			time_fit1->SetParName(0,"constant");
			time_fit1->SetParameter(1,3330);
			time_fit1->SetParameter(2,150);
			time_fit1->SetParName(2,"tau_ls");
			time_fit1->SetParameter(3,3780);
			time_fit1->SetParameter(4,29);
			time_fit1->SetParName(4,"tau_gdls");

			time_fit1->SetParLimits(1,0.,2.*h_coincidence_time_1us_1m->GetBinContent(h_coincidence_time_1us_1m->GetMaximumBin()));
			time_fit1->SetParLimits(2,0.,400.);
			time_fit1->SetParLimits(3,0.,2.*h_coincidence_time_1us_1m->GetBinContent(h_coincidence_time_1us_1m->GetMaximumBin()));
			time_fit1->SetParLimits(4,0.,60.);

		TF1* time_fit075=new TF1("time_fit075","[0]+[1]*TMath::Exp(-x/[2])+[3]*TMath::Exp(-x/[4])",10,2000);
			for(int ibin=0; ibin < 700; ibin++){
				if(h_acc_scaled->GetXaxis()->GetBinCenter(ibin) > 0.75){
					endBin = ibin-1;
					break;
				}
			}
			offset_075 = h_acc_scaled->Integral(0,endBin)/endBin;
	//		time_fit075->FixParameter(0,offset_075);
	//		time_fit075->SetParameter(0,2680);
			time_fit075->SetParameter(0,offset_075);
			time_fit075->SetParName(0,"constant");
			time_fit075->SetParameter(1,h_coincidence_time_1us_075m->GetBinContent(h_coincidence_time_1us_075m->GetMaximumBin()));
			time_fit075->SetParameter(2,210);
			time_fit075->SetParName(2,"tau_ls");
			time_fit075->SetParameter(3,h_coincidence_time_1us_075m->GetBinContent(h_coincidence_time_1us_075m->GetMaximumBin()));
			time_fit075->SetParameter(4,29);
			time_fit075->SetParName(4,"tau_gdls");

			time_fit075->SetParLimits(1,0.,1.2*h_coincidence_time_1us_075m->GetBinContent(h_coincidence_time_1us_075m->GetMaximumBin()));
			time_fit075->SetParLimits(2,0.,400.);
			time_fit075->SetParLimits(3,0.,1.2*h_coincidence_time_1us_075m->GetBinContent(h_coincidence_time_1us_075m->GetMaximumBin()));
			time_fit075->SetParLimits(4,0.,60.);

		TF1* time_fit05=new TF1("time_fit05","[0]+[1]*TMath::Exp(-x/[2])+[3]*TMath::Exp(-x/[4])",10,2000);
			for(int ibin=0; ibin < 700; ibin++){
				if(h_acc_scaled->GetXaxis()->GetBinCenter(ibin) > 0.5){
					endBin = ibin-1;
					break;
				}
			}
			offset_05 = h_acc_scaled->Integral(0,endBin)/endBin;
			cout << "Offset for 0.5 m: " << offset_05 << endl;
	//		time_fit05->FixParameter(0,offset_05);
	//		time_fit05->SetParameter(0,2680);
			time_fit05->SetParameter(0,offset_05);
			time_fit05->SetParName(0,"constant");
			time_fit05->SetParameter(1,h_coincidence_time_1us_05m->GetBinContent(h_coincidence_time_1us_05m->GetMaximumBin()));
			time_fit05->SetParameter(2,210);
			time_fit05->SetParName(2,"tau_ls");
			time_fit05->SetParameter(3,h_coincidence_time_1us_05m->GetBinContent(h_coincidence_time_1us_05m->GetMaximumBin()));
			time_fit05->SetParameter(4,29);
			time_fit05->SetParName(4,"tau_gdls");

			time_fit05->SetParLimits(1,0.,1.2*h_coincidence_time_1us_05m->GetBinContent(h_coincidence_time_1us_05m->GetMaximumBin()));
			time_fit05->SetParLimits(2,0.,400.);
			time_fit05->SetParLimits(3,0.,1.2*h_coincidence_time_1us_05m->GetBinContent(h_coincidence_time_1us_05m->GetMaximumBin()));
			time_fit05->SetParLimits(4,0.,60.);

*/
/*			TF1* time_fit_ls=new TF1("time_fit_ls","[0]+[1]*TMath::Exp(-x/[2])",10,2000);
				offset_ls = h_coincidence_time_1us_ls->GetBinContent(2000);
				time_fit_ls->SetParameter(0,offset_ls);
				time_fit_ls->SetParName(0,"constant");
				time_fit_ls->SetParameter(1,3330);
				time_fit_ls->SetParameter(2,210);
				time_fit_ls->SetParName(2,"tau_ls");

			TF1* time_fit1_ls=new TF1("time_fit1_ls","[0]+[1]*TMath::Exp(-x/[2])",10,2000);
				offset_1_ls = h_coincidence_time_1us_ls_1m->GetBinContent(2000);
				time_fit_ls->SetParameter(0,offset_1_ls);
				time_fit1_ls->SetParName(0,"constant");
				time_fit1_ls->SetParameter(1,3330);
				time_fit1_ls->SetParameter(2,210);
				time_fit1_ls->SetParName(2,"tau_ls");

			TF1* time_fit075_ls=new TF1("time_fit075_ls","[0]+[1]*TMath::Exp(-x/[2])",10,2000);
				offset_075_ls = h_coincidence_time_1us_ls_075m->GetBinContent(2000);
				time_fit_ls->SetParameter(0,offset_075_ls);
				time_fit075_ls->SetParName(0,"constant");
				time_fit075_ls->SetParameter(1,3330);
				time_fit075_ls->SetParameter(2,210);
				time_fit075_ls->SetParName(2,"tau_ls");

			TF1* time_fit05_ls=new TF1("time_fit05_ls","[0]+[1]*TMath::Exp(-x/[2])",10,2000);
				offset_05_ls = h_coincidence_time_1us_ls_05m->GetBinContent(2000);
				time_fit05_ls->SetParameter(0,offset_05_ls);
				time_fit05_ls->SetParName(0,"constant");
				time_fit05_ls->SetParameter(1,3330);
				time_fit05_ls->SetParameter(2,210);
				time_fit05_ls->SetParName(2,"tau_ls");


		TF1* time_fit_gdls=new TF1("time_fit_gdls","[0]+[1]*TMath::Exp(-x/[2])",10,2000);
			offset_gdls = h_coincidence_time_1us_gdls->GetBinContent(2000);
			time_fit_gdls->SetParameter(0,offset_gdls);
			time_fit_gdls->SetParName(0,"constant");
			time_fit_gdls->SetParameter(1,3330);
			time_fit_gdls->SetParameter(2,29);
			time_fit_gdls->SetParName(2,"tau_gdls");

		TF1* time_fit1_gdls=new TF1("time_fit1_gdls","[0]+[1]*TMath::Exp(-x/[2])",10,2000);
			offset_1_gdls = h_coincidence_time_1us_gdls_1m->GetBinContent(2000);
			time_fit1_gdls->SetParameter(0,offset_1_gdls);
			time_fit1_gdls->SetParName(0,"constant");
			time_fit1_gdls->SetParameter(1,3330);
			time_fit1_gdls->SetParameter(2,29);
			time_fit1_gdls->SetParName(2,"tau_gdls");

		TF1* time_fit075_gdls=new TF1("time_fit075_gdls","[0]+[1]*TMath::Exp(-x/[2])",10,2000);
			offset_075_gdls = h_coincidence_time_1us_gdls_075m->GetBinContent(2000);
			time_fit075_gdls->SetParameter(0,offset_075_gdls);
			time_fit075_gdls->SetParName(0,"constant");
			time_fit075_gdls->SetParameter(1,3330);
			time_fit075_gdls->SetParameter(2,29);
			time_fit075_gdls->SetParName(2,"tau_gdls");

		TF1* time_fit05_gdls=new TF1("time_fit05_gdls","[0]+[1]*TMath::Exp(-x/[2])",10,2000);
			offset_05_gdls = h_coincidence_time_1us_gdls_05m->GetBinContent(2000);
			time_fit05_gdls->SetParameter(0,offset_05_gdls);
			time_fit05_gdls->SetParName(0,"constant");
			time_fit05_gdls->SetParameter(1,3330);
			time_fit05_gdls->SetParameter(2,29);
			time_fit05_gdls->SetParName(2,"tau_gdls");*/


//	h_normSub_distance->Rebin(10);
//	h_rateSub_distance->Rebin(10);

		cout << endl << endl << "Norm: " << 100*(h_normSub_distance->Integral(201,500))/(h_acc_norm->Integral(201,500)) << "\t400: " << 100*(h_normSub_distance_400->Integral(201,500))/(h_acc_norm_400->Integral(201,500)) << "\t600: " << 100*(h_normSub_distance_600->Integral(201,500))/(h_acc_norm_600->Integral(201,500)) << "\t800: " << 100*(h_normSub_distance_800->Integral(201,500))/(h_acc_norm_800->Integral(201,500)) << endl;
		cout << "Rate: " << 100*(h_rateSub_distance->Integral(201,500))/(h_acc_scaled->Integral(201,500)) << "\t400: " << 100*(h_rateSub_distance_400->Integral(201,500))/(h_acc_scaled_400->Integral(201,500)) << "\t600: " << 100*(h_rateSub_distance_600->Integral(201,500))/(h_acc_scaled_600->Integral(201,500)) << "\t800: " << 100*(h_rateSub_distance_800->Integral(201,500))/(h_acc_scaled_800->Integral(201,500)) << endl;

	TF1* subFit=new TF1("subFit","[0]",2,5);

			gStyle->SetOptFit(1111);
        char outputname[64];
        sprintf(outputname,"../nH_files/SubtractedAccidentals_%d_EH%dAD%d.root",pd_window_microsec,hall_num,ad_num);
	TFile* outfile=new TFile(outputname, "RECREATE");
		outfile->cd();
			h_normSub_distance->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_normSub_distance->GetYaxis()->SetTitle("Counts");
			cout << endl << endl << endl << "Subtracted Distance (Normalized) Plot" << endl;
			h_normSub_distance->Fit("subFit","R");
			h_normSub_distance->Write();

		//	h_rateSub_distance->Rebin(20);
			h_rateSub_distance->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_rateSub_distance->GetYaxis()->SetTitle("Counts");
			cout << endl << endl << endl << "Rate Corrected Distance (Rate) Plot" << endl;
			h_rateSub_distance->Fit("subFit","R");
			h_rateSub_distance->Write();

				//400 subset
				h_normSub_distance_400->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
				h_normSub_distance_400->GetYaxis()->SetTitle("Counts");
				cout << endl << endl << endl << "Subtracted Distance (Normalized) Plot [400]" << endl;
				h_normSub_distance_400->Fit("subFit","R");
				h_normSub_distance_400->Write();

				h_rateSub_distance_400->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
				h_rateSub_distance_400->GetYaxis()->SetTitle("Counts");
				cout << endl << endl << endl << "Rate Corrected Distance (Rate) Plot [400]" << endl;
				h_rateSub_distance_400->Fit("subFit","R");
				h_rateSub_distance_400->Write();

				//600 subset
				h_normSub_distance_600->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
				h_normSub_distance_600->GetYaxis()->SetTitle("Counts");
				cout << endl << endl << endl << "Subtracted Distance (Normalized) Plot [600]" << endl;
				h_normSub_distance_600->Fit("subFit","R");
				h_normSub_distance_600->Write();

				h_rateSub_distance_600->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
				h_rateSub_distance_600->GetYaxis()->SetTitle("Counts");
				cout << endl << endl << endl << "Rate Corrected Distance (Rate) Plot [600]" << endl;
				h_rateSub_distance_600->Fit("subFit","R");
				h_rateSub_distance_600->Write();

				//800 subset
				h_normSub_distance_800->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
				h_normSub_distance_800->GetYaxis()->SetTitle("Counts");
				cout << endl << endl << endl << "Subtracted Distance (Normalized) Plot [800]" << endl;
				h_normSub_distance_800->Fit("subFit","R");
				h_normSub_distance_800->Write();

				h_rateSub_distance_800->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
				h_rateSub_distance_800->GetYaxis()->SetTitle("Counts");
				cout << endl << endl << endl << "Rate Corrected Distance (Rate) Plot [800]" << endl;
				h_rateSub_distance_800->Fit("subFit","R");
				h_rateSub_distance_800->Write();

				for(int tim = 0; tim <4; tim++){
					for(int stage = 0; stage <5; stage++){
						h_rate_modified[tim][stage]->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
						h_rate_modified[tim][stage]->GetYaxis()->SetTitle("Counts");
						h_rate_modified[tim][stage]->Fit("subFit","R");
						h_rate_modified[tim][stage]->Write();
					}
				}

			h_rateSub_distance_3sig->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_rateSub_distance_3sig->GetYaxis()->SetTitle("Counts");
			h_rateSub_distance_3sig->Write();

			h_normSub_distance_3sig->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_normSub_distance_3sig->GetYaxis()->SetTitle("Counts");
			h_normSub_distance_3sig->Write();

			h_DTnormSub_distance_3sig->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_DTnormSub_distance_3sig->GetYaxis()->SetTitle("Counts");
			h_DTnormSub_distance_3sig->Write();

			h_dist_ratio_rateTOnorm->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_dist_ratio_rateTOnorm->GetYaxis()->SetTitle("Ratio of Rate/Normalized Distance");
			h_dist_ratio_rateTOnorm->Write();

			h_dist_sub_rateTOnorm->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_dist_sub_rateTOnorm->GetYaxis()->SetTitle("Ratio of Rate/Normalized Distance");
			h_dist_sub_rateTOnorm->Write();

			h_IBD_distVStime->SetStats(0);
			h_IBD_distVStime->GetXaxis()->SetTitle("Coincidence Time [us]");
			h_IBD_distVStime->GetYaxis()->SetTitle("Distance [m]");
			h_IBD_distVStime->SetOption("COLZ");
			h_IBD_distVStime->Write();

			h_total_sub_distVStime_norm->SetStats(0);
			h_total_sub_distVStime_norm->GetXaxis()->SetTitle("Coincidence Time [us]");
			h_total_sub_distVStime_norm->GetYaxis()->SetTitle("Distance [m]");
			h_total_sub_distVStime_norm->SetOption("COLZ");
			h_total_sub_distVStime_norm->Write();

			h_total_sub_distVStime_rate->SetStats(0);
			h_total_sub_distVStime_rate->GetXaxis()->SetTitle("Coincidence Time [us]");
			h_total_sub_distVStime_rate->GetYaxis()->SetTitle("Distance [m]");
			h_total_sub_distVStime_rate->SetOption("COLZ");
			h_total_sub_distVStime_rate->Write();

				h_IBD_distVStime_Ep35->SetStats(0);
				h_IBD_distVStime_Ep35->GetXaxis()->SetTitle("Coincidence Time [us]");
				h_IBD_distVStime_Ep35->GetYaxis()->SetTitle("Distance [m]");
				h_IBD_distVStime_Ep35->SetOption("COLZ");
				h_IBD_distVStime_Ep35->Write();

				h_total_sub_distVStime_Ep35_norm->SetStats(0);
				h_total_sub_distVStime_Ep35_norm->GetXaxis()->SetTitle("Coincidence Time [us]");
				h_total_sub_distVStime_Ep35_norm->GetYaxis()->SetTitle("Distance [m]");
				h_total_sub_distVStime_Ep35_norm->SetOption("COLZ");
				h_total_sub_distVStime_Ep35_norm->Write();

				h_total_sub_distVStime_Ep35_rate->SetStats(0);
				h_total_sub_distVStime_Ep35_rate->GetXaxis()->SetTitle("Coincidence Time [us]");
				h_total_sub_distVStime_Ep35_rate->GetYaxis()->SetTitle("Distance [m]");
				h_total_sub_distVStime_Ep35_rate->SetOption("COLZ");
				h_total_sub_distVStime_Ep35_rate->Write();

		//		h_ibd_DT->SetStats(0);
				h_ibd_DT->GetXaxis()->SetTitle("DT [m]");
				h_ibd_DT->GetYaxis()->SetTitle("Counts");
				h_ibd_DT->Write();

		//		h_acc_DT_norm->SetStats(0);
				h_acc_DT_norm->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_norm->GetYaxis()->SetTitle("Counts");
				h_acc_DT_norm->Write();

		//		h_acc_DT_rate->SetStats(0);
				h_acc_DT_rate->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_rate->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_acc_DT_rate->Write();

		//		h_acc_DT_DTnorm->SetStats(0);
				h_acc_DT_DTnorm->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_DTnorm->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_acc_DT_DTnorm->Write();

		//		h_sub_DT_norm->SetStats(0);
				h_sub_DT_norm->GetXaxis()->SetTitle("DT [m]");
				h_sub_DT_norm->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_sub_DT_norm->Write();

		//		h_sub_DT_rate->SetStats(0);
				h_sub_DT_rate->GetXaxis()->SetTitle("DT [m]");
				h_sub_DT_rate->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_sub_DT_rate->Write();

		//		h_sub_DT_DTnorm->SetStats(0);
				h_sub_DT_DTnorm->GetXaxis()->SetTitle("DT [m]");
				h_sub_DT_DTnorm->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_sub_DT_DTnorm->Write();

		//		h_sub_DT_rate_mod->SetStats(0);
				h_sub_DT_rate_mod->GetXaxis()->SetTitle("DT [m]");
				h_sub_DT_rate_mod->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_sub_DT_rate_mod->Write();

			//		h_ibd_DT_Ep35->SetStats(0);
					h_ibd_DT_Ep35->GetXaxis()->SetTitle("DT [m]");
					h_ibd_DT_Ep35->GetYaxis()->SetTitle("Counts/0.02 MeV");
					h_ibd_DT_Ep35->Write();

			//		h_acc_DT_Ep35_norm->SetStats(0);
					h_acc_DT_Ep35_norm->GetXaxis()->SetTitle("DT [m]");
					h_acc_DT_Ep35_norm->GetYaxis()->SetTitle("Counts/0.02 MeV");
					h_acc_DT_Ep35_norm->Write();

			//		h_acc_DT_Ep35_rate->SetStats(0);
					h_acc_DT_Ep35_rate->GetXaxis()->SetTitle("DT [m]");
					h_acc_DT_Ep35_rate->GetYaxis()->SetTitle("Counts/0.02 MeV");
					h_acc_DT_Ep35_rate->Write();

			//		h_acc_DT_Ep35_DTnorm->SetStats(0);
					h_acc_DT_Ep35_DTnorm->GetXaxis()->SetTitle("DT [m]");
					h_acc_DT_Ep35_DTnorm->GetYaxis()->SetTitle("Counts/0.02 MeV");
					h_acc_DT_Ep35_DTnorm->Write();

			//		h_sub_DT_Ep35_norm->SetStats(0);
					h_sub_DT_Ep35_norm->GetXaxis()->SetTitle("DT [m]");
					h_sub_DT_Ep35_norm->GetYaxis()->SetTitle("Counts/0.02 MeV");
					h_sub_DT_Ep35_norm->Write();

			//		h_sub_DT_Ep35_rate->SetStats(0);
					h_sub_DT_Ep35_rate->GetXaxis()->SetTitle("DT [m]");
					h_sub_DT_Ep35_rate->GetYaxis()->SetTitle("Counts/0.02 MeV");
					h_sub_DT_Ep35_rate->Write();

			//		h_sub_DT_Ep35_DTnorm->SetStats(0);
					h_sub_DT_Ep35_DTnorm->GetXaxis()->SetTitle("DT [m]");
					h_sub_DT_Ep35_DTnorm->GetYaxis()->SetTitle("Counts/0.02 MeV");
					h_sub_DT_Ep35_DTnorm->Write();

			//		h_sub_DT_Ep35_rate_mod->SetStats(0);
					h_sub_DT_Ep35_rate_mod->GetXaxis()->SetTitle("DT [m]");
					h_sub_DT_Ep35_rate_mod->GetYaxis()->SetTitle("Counts/0.02 MeV");
					h_sub_DT_Ep35_rate_mod->Write();

		//		h_ibd_DT_3sig->SetStats(0);
				h_ibd_DT_3sig->GetXaxis()->SetTitle("DT [m]");
				h_ibd_DT_3sig->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_ibd_DT_3sig->Write();

		//		h_acc_DT_3sig_norm->SetStats(0);
				h_acc_DT_3sig_norm->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_3sig_norm->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_acc_DT_3sig_norm->Write();

		//		h_acc_DT_3sig_rate->SetStats(0);
				h_acc_DT_3sig_rate->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_3sig_rate->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_acc_DT_3sig_rate->Write();

		//		h_acc_DT_3sig_DTnorm->SetStats(0);
				h_acc_DT_3sig_DTnorm->GetXaxis()->SetTitle("DT [m]");
				h_acc_DT_3sig_DTnorm->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_acc_DT_3sig_DTnorm->Write();

		//		h_sub_DT_3sig_norm->SetStats(0);
				h_sub_DT_3sig_norm->GetXaxis()->SetTitle("DT [m]");
				h_sub_DT_3sig_norm->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_sub_DT_3sig_norm->Write();

		//		h_sub_DT_3sig_rate->SetStats(0);
				h_sub_DT_3sig_rate->GetXaxis()->SetTitle("DT [m]");
				h_sub_DT_3sig_rate->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_sub_DT_3sig_rate->Write();

		//		h_sub_DT_3sig_DTnorm->SetStats(0);
				h_sub_DT_3sig_DTnorm->GetXaxis()->SetTitle("DT [m]");
				h_sub_DT_3sig_DTnorm->GetYaxis()->SetTitle("Counts/0.02 MeV");
				h_sub_DT_3sig_DTnorm->Write();

			h_Edelayed_IBD->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_IBD->GetYaxis()->SetTitle("Counts");
			h_Edelayed_IBD->Write();

			h_Edelayed_IBD_DT800->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_IBD_DT800->GetYaxis()->SetTitle("Counts");
			h_Edelayed_IBD_DT800->Write();

			h_Edelayed_IBD_fine->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_IBD_fine->GetYaxis()->SetTitle("Counts");
			h_Edelayed_IBD_fine->Write();

			h_Edelayed_IBD_fine_DT800->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_IBD_fine_DT800->GetYaxis()->SetTitle("Counts");
			h_Edelayed_IBD_fine_DT800->Write();

			h_Edelayed_IBD_fine_Ep35->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_IBD_fine_Ep35->GetYaxis()->SetTitle("Counts");
			h_Edelayed_IBD_fine_Ep35->Write();

			h_Edelayed_IBD_fine_DT800_Ep35->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_IBD_fine_DT800_Ep35->GetYaxis()->SetTitle("Counts");
			h_Edelayed_IBD_fine_DT800_Ep35->Write();

/*			for(int iz = 0; iz < NzBins; iz++){
				h_Edelayed_ibd_z[iz]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_ibd_z[iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_ibd_z[iz]->Write();

				h_Edelayed_sub_rate_z[iz]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_rate_z[iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_rate_z[iz]->Write();

				h_Edelayed_sub_norm_z[iz]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_norm_z[iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_norm_z[iz]->Write();

				h_Edelayed_sub_DTnorm_z[iz]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_DTnorm_z[iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_DTnorm_z[iz]->Write();

				h_Edelayed_ibd_DT800_z[iz]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_ibd_DT800_z[iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_ibd_DT800_z[iz]->Write();

				h_Edelayed_sub_rate_DT800_z[iz]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_rate_DT800_z[iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_rate_DT800_z[iz]->Write();

				h_Edelayed_sub_norm_DT800_z[iz]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_norm_DT800_z[iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_norm_DT800_z[iz]->Write();

				h_Edelayed_sub_DTnorm_DT800_z[iz]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_DTnorm_DT800_z[iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_DTnorm_DT800_z[iz]->Write();

			}

			for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
				h_Edelayed_ibd_r2[ir2]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_ibd_r2[ir2]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_ibd_r2[ir2]->Write();

				h_Edelayed_sub_rate_r2[ir2]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_rate_r2[ir2]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_rate_r2[ir2]->Write();

				h_Edelayed_sub_norm_r2[ir2]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_norm_r2[ir2]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_norm_r2[ir2]->Write();

				h_Edelayed_sub_DTnorm_r2[ir2]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_DTnorm_r2[ir2]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_DTnorm_r2[ir2]->Write();

				h_Edelayed_ibd_DT800_r2[ir2]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_ibd_DT800_r2[ir2]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_ibd_DT800_r2[ir2]->Write();

				h_Edelayed_sub_rate_DT800_r2[ir2]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_rate_DT800_r2[ir2]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_rate_DT800_r2[ir2]->Write();

				h_Edelayed_sub_norm_DT800_r2[ir2]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_norm_DT800_r2[ir2]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_norm_DT800_r2[ir2]->Write();

				h_Edelayed_sub_DTnorm_DT800_r2[ir2]->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_DTnorm_DT800_r2[ir2]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_DTnorm_DT800_r2[ir2]->Write();

			}

			for(int iz = 0; iz < NzBins; iz++){
				for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
					h_Edelayed_ibd_zVSr2[ir2][iz]->GetXaxis()->SetTitle("Energy [MeV]");
					h_Edelayed_ibd_zVSr2[ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_Edelayed_ibd_zVSr2[ir2][iz]->Write();

					h_Edelayed_sub_rate_zVSr2[ir2][iz]->GetXaxis()->SetTitle("Energy [MeV]");
					h_Edelayed_sub_rate_zVSr2[ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_Edelayed_sub_rate_zVSr2[ir2][iz]->Write();

					h_Edelayed_sub_norm_zVSr2[ir2][iz]->GetXaxis()->SetTitle("Energy [MeV]");
					h_Edelayed_sub_norm_zVSr2[ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_Edelayed_sub_norm_zVSr2[ir2][iz]->Write();

					h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]->GetXaxis()->SetTitle("Energy [MeV]");
					h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]->Write();

					h_Edelayed_ibd_DT800_zVSr2[ir2][iz]->GetXaxis()->SetTitle("Energy [MeV]");
					h_Edelayed_ibd_DT800_zVSr2[ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_Edelayed_ibd_DT800_zVSr2[ir2][iz]->Write();

					h_Edelayed_sub_rate_DT800_zVSr2[ir2][iz]->GetXaxis()->SetTitle("Energy [MeV]");
					h_Edelayed_sub_rate_DT800_zVSr2[ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_Edelayed_sub_rate_DT800_zVSr2[ir2][iz]->Write();

					h_Edelayed_sub_norm_DT800_zVSr2[ir2][iz]->GetXaxis()->SetTitle("Energy [MeV]");
					h_Edelayed_sub_norm_DT800_zVSr2[ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_Edelayed_sub_norm_DT800_zVSr2[ir2][iz]->Write();

					h_Edelayed_sub_DTnorm_DT800_zVSr2[ir2][iz]->GetXaxis()->SetTitle("Energy [MeV]");
					h_Edelayed_sub_DTnorm_DT800_zVSr2[ir2][iz]->GetYaxis()->SetTitle("Counts");
					h_Edelayed_sub_DTnorm_DT800_zVSr2[ir2][iz]->Write();

				}
			}*/

			h_Eprompt_sub->GetXaxis()->SetTitle("Energy [MeV]");
			h_Eprompt_sub->GetYaxis()->SetTitle("Counts");
			h_Eprompt_sub->Write();

			h_Edelayed_sub->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub->GetXaxis()->SetRangeUser(1.5,3);
/*			TF1* delayedFit = new TF1("delayedFit", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
				delayedFit->SetParameter(0,(h_Edelayed_sub->GetBinContent(h_Edelayed_sub->FindBin(2.3)))/3.); //normalization
				delayedFit->SetParameter(1,0.8); //alpha
				delayedFit->SetParameter(2,2.3); //mu
				delayedFit->SetParameter(3,0.135); //sigma
				delayedFit->SetParameter(4,3); //lambda

				delayedFit->SetParName(0,"N"); //normalization
				delayedFit->SetParName(1,"alpha"); //alpha
				delayedFit->SetParName(2,"mu"); //mu
				delayedFit->SetParName(3,"sigma"); //sigma
				delayedFit->SetParName(4,"lambda"); //lambda
			h_Edelayed_sub->Fit("delayedFit", "R");*/
			h_Edelayed_sub->Write();

			h_Eprompt_sub_norm->GetXaxis()->SetTitle("Energy [MeV]");
			h_Eprompt_sub_norm->GetYaxis()->SetTitle("Counts");
			h_Eprompt_sub_norm->Write();

			h_Edelayed_sub_norm->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub_norm->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_norm->GetXaxis()->SetRangeUser(1.5,3);
/*			TF1* delayedFit_norm = new TF1("delayedFit_norm", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
				delayedFit_norm->SetParameter(0,(h_Edelayed_sub_norm->GetBinContent(h_Edelayed_sub_norm->FindBin(2.3)))/3.); //normalization
				delayedFit_norm->SetParameter(1,0.8); //alpha
				delayedFit_norm->SetParameter(2,2.3); //mu
				delayedFit_norm->SetParameter(3,0.135); //sigma
				delayedFit_norm->SetParameter(4,3); //lambda

				delayedFit_norm->SetParName(0,"N"); //normalization
				delayedFit_norm->SetParName(1,"alpha"); //alpha
				delayedFit_norm->SetParName(2,"mu"); //mu
				delayedFit_norm->SetParName(3,"sigma"); //sigma
				delayedFit_norm->SetParName(4,"lambda"); //lambda
			h_Edelayed_sub_norm->Fit("delayedFit_norm", "R");*/
			h_Edelayed_sub_norm->Write();

			h_Eprompt_sub_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
			h_Eprompt_sub_DTnorm->GetYaxis()->SetTitle("Counts");
			h_Eprompt_sub_DTnorm->Write();

			h_Edelayed_sub_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub_DTnorm->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_DTnorm->GetXaxis()->SetRangeUser(1.5,3);
/*			TF1* delayedFit_DTnorm = new TF1("delayedFit_DTnorm", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
				delayedFit_DTnorm->SetParameter(0,(h_Edelayed_sub_DTnorm->GetBinContent(h_Edelayed_sub_DTnorm->FindBin(2.3)))/3.); //normalization
				delayedFit_DTnorm->SetParameter(1,0.8); //alpha
				delayedFit_DTnorm->SetParameter(2,2.3); //mu
				delayedFit_DTnorm->SetParameter(3,0.135); //sigma
				delayedFit_DTnorm->SetParameter(4,3); //lambda

				delayedFit_DTnorm->SetParName(0,"N"); //normalization
				delayedFit_DTnorm->SetParName(1,"alpha"); //alpha
				delayedFit_DTnorm->SetParName(2,"mu"); //mu
				delayedFit_DTnorm->SetParName(3,"sigma"); //sigma
				delayedFit_DTnorm->SetParName(4,"lambda"); //lambda
			h_Edelayed_sub_DTnorm->Fit("delayedFit_DTnorm", "R");*/
			h_Edelayed_sub_DTnorm->Write();

				h_Eprompt_sub_DT800->GetXaxis()->SetTitle("Energy [MeV]");
				h_Eprompt_sub_DT800->GetYaxis()->SetTitle("Counts");
				h_Eprompt_sub_DT800->Write();

				h_Edelayed_sub_DT800->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_DT800->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_DT800->GetXaxis()->SetRangeUser(1.5,3);
				
				h_Eprompt_sub_DT800_3sig->GetXaxis()->SetTitle("Energy [MeV]");
				h_Eprompt_sub_DT800_3sig->GetYaxis()->SetTitle("Counts");
				h_Eprompt_sub_DT800_3sig->Write();
/*		TF1* delayedFit_DT800 = new TF1("delayedFit_DT800", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
			delayedFit_DT800->SetParameter(0,(h_Edelayed_sub_DT800->GetBinContent(h_Edelayed_sub_DT800->FindBin(2.3)))/3.); //normalization
			delayedFit_DT800->SetParameter(1,0.8); //alpha
			delayedFit_DT800->SetParameter(2,2.3); //mu
			delayedFit_DT800->SetParameter(3,0.135); //sigma
			delayedFit_DT800->SetParameter(4,3); //lambda

			delayedFit_DT800->SetParName(0,"N"); //normalization
			delayedFit_DT800->SetParName(1,"alpha"); //alpha
			delayedFit_DT800->SetParName(2,"mu"); //mu
			delayedFit_DT800->SetParName(3,"sigma"); //sigma
			delayedFit_DT800->SetParName(4,"lambda"); //lambda
		h_Edelayed_sub_DT800->Fit("delayedFit_DT800", "R");*/
				h_Edelayed_sub_DT800->Write();

				h_Eprompt_sub_DT800_norm->GetXaxis()->SetTitle("Energy [MeV]");
				h_Eprompt_sub_DT800_norm->GetYaxis()->SetTitle("Counts");
				h_Eprompt_sub_DT800_norm->Write();

				h_Edelayed_sub_DT800_norm->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_DT800_norm->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_DT800_norm->GetXaxis()->SetRangeUser(1.5,3);
/*		TF1* delayedFit_DT800_norm = new TF1("delayedFit_DT800_norm", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
			delayedFit_DT800_norm->SetParameter(0,(h_Edelayed_sub_DT800_norm->GetBinContent(h_Edelayed_sub_DT800_norm->FindBin(2.3)))/3.); //normalization
			delayedFit_DT800_norm->SetParameter(1,0.8); //alpha
			delayedFit_DT800_norm->SetParameter(2,2.3); //mu
			delayedFit_DT800_norm->SetParameter(3,0.135); //sigma
			delayedFit_DT800_norm->SetParameter(4,3); //lambda

			delayedFit_DT800_norm->SetParName(0,"N"); //normalization
			delayedFit_DT800_norm->SetParName(1,"alpha"); //alpha
			delayedFit_DT800_norm->SetParName(2,"mu"); //mu
			delayedFit_DT800_norm->SetParName(3,"sigma"); //sigma
			delayedFit_DT800_norm->SetParName(4,"lambda"); //lambda
		h_Edelayed_sub_DT800_norm->Fit("delayedFit_DT800_norm", "R");*/
				h_Edelayed_sub_DT800_norm->Write();

				h_Eprompt_sub_DT800_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
				h_Eprompt_sub_DT800_DTnorm->GetYaxis()->SetTitle("Counts");
				h_Eprompt_sub_DT800_DTnorm->Write();

				h_Edelayed_sub_DT800_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_DT800_DTnorm->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_DT800_DTnorm->GetXaxis()->SetRangeUser(1.5,3);
/*		TF1* delayedFit_DT800_DTnorm = new TF1("delayedFit_DT800_DTnorm", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
			delayedFit_DT800_DTnorm->SetParameter(0,(h_Edelayed_sub_DT800_DTnorm->GetBinContent(h_Edelayed_sub_DT800_DTnorm->FindBin(2.3)))/3.); //normalization
			delayedFit_DT800_DTnorm->SetParameter(1,0.8); //alpha
			delayedFit_DT800_DTnorm->SetParameter(2,2.3); //mu
			delayedFit_DT800_DTnorm->SetParameter(3,0.135); //sigma
			delayedFit_DT800_DTnorm->SetParameter(4,3); //lambda

			delayedFit_DT800_DTnorm->SetParName(0,"N"); //normalization
			delayedFit_DT800_DTnorm->SetParName(1,"alpha"); //alpha
			delayedFit_DT800_DTnorm->SetParName(2,"mu"); //mu
			delayedFit_DT800_DTnorm->SetParName(3,"sigma"); //sigma
			delayedFit_DT800_DTnorm->SetParName(4,"lambda"); //lambda
		h_Edelayed_sub_DT800_DTnorm->Fit("delayedFit_DT800_DTnorm", "R");*/

		//Sam's version
/*		TF1* delayedFit_DT800_DTnorm_Sam = new TF1("delayedFit_DT800_DTnorm_Sam", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(exp([4]*[2])-1) * exp([3]*[3]*[4]*[4]) * exp(2*[4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
			delayedFit_DT800_DTnorm_Sam->SetParameter(0,(h_Edelayed_sub_DT800_DTnorm->GetBinContent(h_Edelayed_sub_DT800_DTnorm->FindBin(2.3)))/3.); //normalization
			delayedFit_DT800_DTnorm_Sam->SetParameter(1,0.8); //alpha
			delayedFit_DT800_DTnorm_Sam->SetParameter(2,2.3); //mu
			delayedFit_DT800_DTnorm_Sam->SetParameter(3,0.135); //sigma
			delayedFit_DT800_DTnorm_Sam->SetParameter(4,3); //lambda
			delayedFit_DT800_DTnorm_Sam->SetLineStyle(7);
			delayedFit_DT800_DTnorm_Sam->SetLineColor(kGreen);
		h_Edelayed_sub_DT800_DTnorm->Fit("delayedFit_DT800_DTnorm_Sam", "R");*/
				h_Edelayed_sub_DT800_DTnorm->Write();

			h_Edelayed_sub_fine->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub_fine->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_fine->GetXaxis()->SetRangeUser(1.5,3);
			h_Edelayed_sub_fine->Write();

			h_Edelayed_sub_fine_norm->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub_fine_norm->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_fine_norm->GetXaxis()->SetRangeUser(1.5,3);
			h_Edelayed_sub_fine_norm->Write();

			h_Edelayed_sub_fine_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub_fine_DTnorm->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_fine_DTnorm->GetXaxis()->SetRangeUser(1.5,3);
			h_Edelayed_sub_fine_DTnorm->Write();

				h_Edelayed_sub_fine_DT800->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_fine_DT800->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_fine_DT800->GetXaxis()->SetRangeUser(1.5,3);
				h_Edelayed_sub_fine_DT800->Write();

				h_Edelayed_sub_fine_DT800_norm->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_fine_DT800_norm->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_fine_DT800_norm->GetXaxis()->SetRangeUser(1.5,3);
				h_Edelayed_sub_fine_DT800_norm->Write();

				h_Edelayed_sub_fine_DT800_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_fine_DT800_DTnorm->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_fine_DT800_DTnorm->GetXaxis()->SetRangeUser(1.5,3);
				h_Edelayed_sub_fine_DT800_DTnorm->Write();


			h_Edelayed_sub_fine_Ep35->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub_fine_Ep35->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_fine_Ep35->GetXaxis()->SetRangeUser(1.5,3);
			h_Edelayed_sub_fine_Ep35->Write();

			h_Edelayed_sub_fine_Ep35_norm->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub_fine_Ep35_norm->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_fine_Ep35_norm->GetXaxis()->SetRangeUser(1.5,3);
			h_Edelayed_sub_fine_Ep35_norm->Write();

			h_Edelayed_sub_fine_Ep35_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
			h_Edelayed_sub_fine_Ep35_DTnorm->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_fine_Ep35_DTnorm->GetXaxis()->SetRangeUser(1.5,3);
			h_Edelayed_sub_fine_Ep35_DTnorm->Write();

				h_Edelayed_sub_fine_DT800_Ep35->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_fine_DT800_Ep35->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_fine_DT800_Ep35->GetXaxis()->SetRangeUser(1.5,3);
				h_Edelayed_sub_fine_DT800_Ep35->Write();

				h_Edelayed_sub_fine_DT800_Ep35_norm->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_fine_DT800_Ep35_norm->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_fine_DT800_Ep35_norm->GetXaxis()->SetRangeUser(1.5,3);
				h_Edelayed_sub_fine_DT800_Ep35_norm->Write();

				h_Edelayed_sub_fine_DT800_Ep35_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
				h_Edelayed_sub_fine_DT800_Ep35_DTnorm->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_fine_DT800_Ep35_DTnorm->GetXaxis()->SetRangeUser(1.5,3);
				h_Edelayed_sub_fine_DT800_Ep35_DTnorm->Write();


			h_ratio->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
			h_ratio->GetYaxis()->SetTitle("Ratio of rate/normalized");
			h_ratio->Write();


			h_delayed_energy_vs_distance_rate_ratio->SetStats(0);
			h_delayed_energy_vs_distance_rate_ratio->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_vs_distance_rate_ratio->GetYaxis()->SetTitle("Distance [m]");
			h_delayed_energy_vs_distance_rate_ratio->SetOption("COLZ");
			h_delayed_energy_vs_distance_rate_ratio->GetXaxis()->SetRangeUser(1.75,3);
			h_delayed_energy_vs_distance_rate_ratio->Write();

			h_delayed_energy_vs_distance_rate_sub->SetStats(0);
			h_delayed_energy_vs_distance_rate_sub->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_vs_distance_rate_sub->GetYaxis()->SetTitle("Distance [m]");
			h_delayed_energy_vs_distance_rate_sub->SetOption("COLZ");
			h_delayed_energy_vs_distance_rate_sub->GetXaxis()->SetRangeUser(1.75,3);
			h_delayed_energy_vs_distance_rate_sub->Write();

			h_delayed_energy_vs_distance_norm_ratio->SetStats(0);
			h_delayed_energy_vs_distance_norm_ratio->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_vs_distance_norm_ratio->GetYaxis()->SetTitle("Distance [m]");
			h_delayed_energy_vs_distance_norm_ratio->SetOption("COLZ");
			h_delayed_energy_vs_distance_norm_ratio->GetXaxis()->SetRangeUser(1.75,3);
			h_delayed_energy_vs_distance_norm_ratio->Write();

			h_delayed_energy_vs_distance_norm_sub->SetStats(0);
			h_delayed_energy_vs_distance_norm_sub->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_delayed_energy_vs_distance_norm_sub->GetYaxis()->SetTitle("Distance [m]");
			h_delayed_energy_vs_distance_norm_sub->SetOption("COLZ");
			h_delayed_energy_vs_distance_norm_sub->GetXaxis()->SetRangeUser(1.75,3);
			h_delayed_energy_vs_distance_norm_sub->Write();

			h_prompt_energy_vs_distance_rate_ratio->SetStats(0);
			h_prompt_energy_vs_distance_rate_ratio->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_prompt_energy_vs_distance_rate_ratio->GetYaxis()->SetTitle("Distance [m]");
			h_prompt_energy_vs_distance_rate_ratio->SetOption("COLZ");
			h_prompt_energy_vs_distance_rate_ratio->Write();

			h_prompt_energy_vs_distance_rate_sub->SetStats(0);
			h_prompt_energy_vs_distance_rate_sub->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_prompt_energy_vs_distance_rate_sub->GetYaxis()->SetTitle("Distance [m]");
			h_prompt_energy_vs_distance_rate_sub->SetOption("COLZ");
			h_prompt_energy_vs_distance_rate_sub->Write();

			h_prompt_energy_vs_distance_norm_ratio->SetStats(0);
			h_prompt_energy_vs_distance_norm_ratio->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_prompt_energy_vs_distance_norm_ratio->GetYaxis()->SetTitle("Distance [m]");
			h_prompt_energy_vs_distance_norm_ratio->SetOption("COLZ");
			h_prompt_energy_vs_distance_norm_ratio->Write();

			h_prompt_energy_vs_distance_norm_sub->SetStats(0);
			h_prompt_energy_vs_distance_norm_sub->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_prompt_energy_vs_distance_norm_sub->GetYaxis()->SetTitle("Distance [m]");
			h_prompt_energy_vs_distance_norm_sub->SetOption("COLZ");
			h_prompt_energy_vs_distance_norm_sub->Write();

			for(int dist=0; dist<6; dist++){


				/*h_ratio_delayed_energy_norm_dist[dist]->Rebin(4);
				h_ratio_prompt_energy_norm_dist[dist]->Rebin(4);
				h_ratio_delayed_energy_scaled_dist[dist]->Rebin(4);
				h_ratio_prompt_energy_scaled_dist[dist]->Rebin(4);*/


				//For scaled
				h_sub_delayed_energy_scaled_dist[dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_sub_delayed_energy_scaled_dist[dist]->GetYaxis()->SetTitle("Counts");
				h_sub_delayed_energy_scaled_dist[dist]->GetXaxis()->SetRangeUser(1.75,3);
				h_sub_delayed_energy_scaled_dist[dist]->Write();
				h_sub_prompt_energy_scaled_dist[dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_sub_prompt_energy_scaled_dist[dist]->GetYaxis()->SetTitle("Counts");
				h_sub_prompt_energy_scaled_dist[dist]->Write();

				h_ratio_delayed_energy_scaled_dist[dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_ratio_delayed_energy_scaled_dist[dist]->GetYaxis()->SetTitle("Ratio");
				h_ratio_delayed_energy_scaled_dist[dist]->GetXaxis()->SetRangeUser(1.75,3);
				h_ratio_delayed_energy_scaled_dist[dist]->Write();
				h_ratio_prompt_energy_scaled_dist[dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_ratio_prompt_energy_scaled_dist[dist]->GetYaxis()->SetTitle("Ratio");
				h_ratio_prompt_energy_scaled_dist[dist]->Write();

				//For normalized
				h_sub_delayed_energy_norm_dist[dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_sub_delayed_energy_norm_dist[dist]->GetYaxis()->SetTitle("Counts");
				h_sub_delayed_energy_norm_dist[dist]->GetXaxis()->SetRangeUser(1.75,3);
				h_sub_delayed_energy_norm_dist[dist]->Write();
				h_sub_prompt_energy_norm_dist[dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_sub_prompt_energy_norm_dist[dist]->GetYaxis()->SetTitle("Counts");
				h_sub_prompt_energy_norm_dist[dist]->Write();

				h_ratio_delayed_energy_norm_dist[dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_ratio_delayed_energy_norm_dist[dist]->GetYaxis()->SetTitle("Ratio");
				h_ratio_delayed_energy_norm_dist[dist]->GetXaxis()->SetRangeUser(1.75,3);
				h_ratio_delayed_energy_norm_dist[dist]->Write();
				h_ratio_prompt_energy_norm_dist[dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_ratio_prompt_energy_norm_dist[dist]->GetYaxis()->SetTitle("Ratio");
				h_ratio_prompt_energy_norm_dist[dist]->Write();

				h_ratio_delayed_energy_scaledTOnorm_dist[dist]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_ratio_delayed_energy_scaledTOnorm_dist[dist]->GetYaxis()->SetTitle("Ratio");
				h_ratio_delayed_energy_scaledTOnorm_dist[dist]->GetXaxis()->SetRangeUser(1.75,3);
				h_ratio_delayed_energy_scaledTOnorm_dist[dist]->Write();
				h_ratio_prompt_energy_scaledTOnorm_dist[dist]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
				h_ratio_prompt_energy_scaledTOnorm_dist[dist]->GetYaxis()->SetTitle("Ratio");
				h_ratio_prompt_energy_scaledTOnorm_dist[dist]->Write();

			}

			h_sub_delayed_energy_norm_dist_2plus->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_sub_delayed_energy_norm_dist_2plus->GetYaxis()->SetTitle("Counts");
			h_sub_delayed_energy_norm_dist_2plus->GetXaxis()->SetRangeUser(1.75,3);
			h_sub_delayed_energy_norm_dist_2plus->Write();

			h_sub_prompt_energy_norm_dist_2plus->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_sub_prompt_energy_norm_dist_2plus->GetYaxis()->SetTitle("Counts");
			h_sub_prompt_energy_norm_dist_2plus->Write();

			h_sub_delayed_energy_scaled_dist_2plus->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_sub_delayed_energy_scaled_dist_2plus->GetYaxis()->SetTitle("Counts");
			h_sub_delayed_energy_scaled_dist_2plus->GetXaxis()->SetRangeUser(1.75,3);
			h_sub_delayed_energy_scaled_dist_2plus->Write();

			h_sub_prompt_energy_scaled_dist_2plus->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_sub_prompt_energy_scaled_dist_2plus->GetYaxis()->SetTitle("Counts");
			h_sub_prompt_energy_scaled_dist_2plus->Write();

			h_ratio_delayed_energy_scaled_dist_2plus->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_ratio_delayed_energy_scaled_dist_2plus->GetYaxis()->SetTitle("IBD/Acc");
			h_ratio_delayed_energy_scaled_dist_2plus->GetXaxis()->SetRangeUser(1.75,3);
			h_ratio_delayed_energy_scaled_dist_2plus->Write();

			h_ratio_prompt_energy_scaled_dist_2plus->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_ratio_prompt_energy_scaled_dist_2plus->GetYaxis()->SetTitle("IBD/Acc");
			h_ratio_prompt_energy_scaled_dist_2plus->Write();

			h_ratio_delayed_energy_norm_dist_2plus->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_ratio_delayed_energy_norm_dist_2plus->GetYaxis()->SetTitle("IBD/Acc");
			h_ratio_delayed_energy_norm_dist_2plus->GetXaxis()->SetRangeUser(1.75,3);
			h_ratio_delayed_energy_norm_dist_2plus->Write();

			h_ratio_prompt_energy_norm_dist_2plus->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_ratio_prompt_energy_norm_dist_2plus->GetYaxis()->SetTitle("IBD/Acc");
			h_ratio_prompt_energy_norm_dist_2plus->Write();

			h_ratio_delayed_energy_scaledTOnorm_dist_2plus->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_ratio_delayed_energy_scaledTOnorm_dist_2plus->GetYaxis()->SetTitle("Scaled/Norm");
			h_ratio_delayed_energy_scaledTOnorm_dist_2plus->GetXaxis()->SetRangeUser(1.75,3);
			h_ratio_delayed_energy_scaledTOnorm_dist_2plus->Write();

			h_ratio_prompt_energy_scaledTOnorm_dist_2plus->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_ratio_prompt_energy_scaledTOnorm_dist_2plus->GetYaxis()->SetTitle("Scaled/Norm");
			h_ratio_prompt_energy_scaledTOnorm_dist_2plus->Write();

			h_sub_delayed_energy_scaled_p35->SetStats(0);
			h_sub_delayed_energy_scaled_p35->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_sub_delayed_energy_scaled_p35->GetYaxis()->SetTitle("Counts");
			h_sub_delayed_energy_scaled_p35->Write();

			h_sub_PvsD_energy_rate->SetOption("COLZ");
			h_sub_PvsD_energy_rate->SetStats(0);
			h_sub_PvsD_energy_rate->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_sub_PvsD_energy_rate->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_sub_PvsD_energy_rate->GetYaxis()->SetRangeUser(1.75,3);
			h_sub_PvsD_energy_rate->Write();

			h_sub_PvsD_energy_rate_1m->SetOption("COLZ");
			h_sub_PvsD_energy_rate_1m->SetStats(0);
			h_sub_PvsD_energy_rate_1m->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_sub_PvsD_energy_rate_1m->GetYaxis()->SetTitle("Delayed Energy [MeV]");
			h_sub_PvsD_energy_rate_1m->GetYaxis()->SetRangeUser(1.75,3);
			h_sub_PvsD_energy_rate_1m->Write();

					h_sub_PvsD_energy_norm->SetOption("COLZ");
					h_sub_PvsD_energy_norm->SetStats(0);
					h_sub_PvsD_energy_norm->GetXaxis()->SetTitle("Prompt Energy [MeV]");
					h_sub_PvsD_energy_norm->GetYaxis()->SetTitle("Delayed Energy [MeV]");
					h_sub_PvsD_energy_norm->GetYaxis()->SetRangeUser(1.75,3);
					h_sub_PvsD_energy_norm->Write();


			h_normScale->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_normScale->GetYaxis()->SetTitle("Scale");
			h_normScale->Write();

			h_rateScale->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_rateScale->GetYaxis()->SetTitle("Scale");
			h_rateScale->Write();

			h_ratioScale->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
			h_ratioScale->GetYaxis()->SetTitle("Norm/Rate");
			h_ratioScale->Write();

			h_scale_diff->GetXaxis()->SetTitle("Difference in scale (Rate - Norm)");
			h_scale_diff->GetYaxis()->SetTitle("Counts");
			h_scale_diff->Write();

/*			h_prompt_energy_vs_time_DTnorm_sub->SetStats(0);
			h_prompt_energy_vs_time_DTnorm_sub->GetXaxis()->SetTitle("Time [us]");
			h_prompt_energy_vs_time_DTnorm_sub->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_prompt_energy_vs_time_DTnorm_sub->SetOption("COLZ");
			h_prompt_energy_vs_time_DTnorm_sub->Write();

			h_prompt_energy_vs_time_DTnorm_DT800_sub->SetStats(0);
			h_prompt_energy_vs_time_DTnorm_DT800_sub->GetXaxis()->SetTitle("Time [us]");
			h_prompt_energy_vs_time_DTnorm_DT800_sub->GetYaxis()->SetTitle("Prompt Energy [MeV]");
			h_prompt_energy_vs_time_DTnorm_DT800_sub->SetOption("COLZ");
			h_prompt_energy_vs_time_DTnorm_DT800_sub->Write();*/

	//After checks are done, remove this part:
	


	//Check part done

	double xBins = 0;
	double eff_numer = 0;
	double eff_denom = 0;
	double acc_const = 0;

/*			cout << endl << endl << "no dist cut:" << endl;
			//h_coincidence_time_1us->SetStats(0);
			h_coincidence_time_1us->SetMinimum(0);
			h_coincidence_time_1us->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us->Fit("time_fit","R");
			h_coincidence_time_1us->Write();
					xBins = h_coincidence_time_1us->GetNbinsX();
					acc_const = time_fit->GetParameter("constant");
					eff_denom = h_coincidence_time_1us->Integral(0,xBins);
					eff_denom = eff_denom - acc_const*xBins;
					for(int iBin = 0; iBin < xBins+1; iBin++){
						eff_numer = 0;
						eff_numer = h_coincidence_time_1us->Integral(0,iBin);
						eff_numer = eff_numer - acc_const*iBin;
						h_time_efficiency_1us->SetBinContent(iBin, eff_numer/eff_denom);
					}
				eff_numer = h_coincidence_time_1us->Integral(0,30);
				eff_numer = eff_numer - acc_const*20;
				cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;
				cout << "Efficiency at " << h_time_efficiency_1us->GetBinCenter(40) << " us: " << h_time_efficiency_1us->GetBinContent(40) << endl;

			cout << endl << endl << "1m cut:" << endl;
			//h_coincidence_time_1us_1m->SetStats(0);
			h_coincidence_time_1us_1m->SetMinimum(0);
			h_coincidence_time_1us_1m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_1m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_1m->Fit("time_fit1","RI");
			h_coincidence_time_1us_1m->Write();
					xBins = h_coincidence_time_1us_1m->GetNbinsX();
					acc_const = time_fit1->GetParameter("constant");
					eff_denom = h_coincidence_time_1us_1m->Integral(0,xBins);
					eff_denom = eff_denom - acc_const*xBins;
					for(int iBin = 0; iBin < xBins+1; iBin++){
						eff_numer = 0;
						eff_numer = h_coincidence_time_1us_1m->Integral(0,iBin);
						eff_numer = eff_numer - acc_const*iBin;
						h_time_efficiency_1us_1m->SetBinContent(iBin, eff_numer/eff_denom);
					}
				eff_numer = h_coincidence_time_1us_1m->Integral(0,30);
				eff_numer = eff_numer - acc_const*20;
				cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;
				cout << "Efficiency at " << h_time_efficiency_1us_1m->GetBinCenter(40) << " us: " << h_time_efficiency_1us_1m->GetBinContent(40) << endl;


			cout << endl << endl << "0.75m cut:" << endl;
			//h_coincidence_time_1us_075m->SetStats(0);
			h_coincidence_time_1us_075m->SetMinimum(0);
			h_coincidence_time_1us_075m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_075m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_075m->Fit("time_fit075","RI");
			h_coincidence_time_1us_075m->Write();
					xBins = h_coincidence_time_1us_075m->GetNbinsX();
					acc_const = time_fit075->GetParameter("constant");
					eff_denom = h_coincidence_time_1us_075m->Integral(0,xBins);
					eff_denom = eff_denom - acc_const*xBins;
					for(int iBin = 0; iBin < xBins+1; iBin++){
						eff_numer = 0;
						eff_numer = h_coincidence_time_1us_075m->Integral(0,iBin);
						eff_numer = eff_numer - acc_const*iBin;
						h_time_efficiency_1us_075m->SetBinContent(iBin, eff_numer/eff_denom);
					}
				eff_numer = h_coincidence_time_1us_075m->Integral(0,30);
				eff_numer = eff_numer - acc_const*20;
				cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;
				cout << "Efficiency at " << h_time_efficiency_1us_075m->GetBinCenter(40) << " us: " << h_time_efficiency_1us_075m->GetBinContent(40) << endl;

			cout << endl << endl << "0.5m cut:" << endl;
			//h_coincidence_time_1us_05m->SetStats(0);
			h_coincidence_time_1us_05m->SetMinimum(0);
			h_coincidence_time_1us_05m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_05m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_05m->Fit("time_fit05","RI");
			h_coincidence_time_1us_05m->Write();
					xBins = h_coincidence_time_1us_05m->GetNbinsX();
					acc_const = time_fit05->GetParameter("constant");
					eff_denom = h_coincidence_time_1us_05m->Integral(0,xBins);
					eff_denom = eff_denom - acc_const*xBins;
					for(int iBin = 0; iBin < xBins+1; iBin++){
						eff_numer = 0;
						eff_numer = h_coincidence_time_1us_05m->Integral(0,iBin);
						eff_numer = eff_numer - acc_const*iBin;
						h_time_efficiency_1us_05m->SetBinContent(iBin, eff_numer/eff_denom);
					}
				eff_numer = h_coincidence_time_1us_05m->Integral(0,30);
				eff_numer = eff_numer - acc_const*20;
				cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;
				cout << "Efficiency at " << h_time_efficiency_1us_05m->GetBinCenter(40) << " us: " << h_time_efficiency_1us_05m->GetBinContent(40) << endl;

			h_time_efficiency_1us->SetStats(0);
			h_time_efficiency_1us->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us->GetYaxis()->SetTitle("Efficiency");
				TF1* eff=new TF1("eff","efficiency(x,[0],[1],[2],[3])",0,2000);
				eff->FixParameter(0,time_fit->GetParameter("p1"));
				eff->FixParameter(1,time_fit->GetParameter("tau_ls"));
				eff->FixParameter(2,time_fit->GetParameter("p3"));
				eff->FixParameter(3,time_fit->GetParameter("tau_gdls"));
			h_time_efficiency_1us->Fit("eff","R");
			h_time_efficiency_1us->Write();
			cout << "For no distance cut... " << endl << "\t400us: " << efficiency(400,time_fit->GetParameter("p1"),time_fit->GetParameter("tau_ls"),time_fit->GetParameter("p3"),time_fit->GetParameter("tau_gdls")) << endl << "\t600us: " << efficiency(600,time_fit->GetParameter("p1"),time_fit->GetParameter("tau_ls"),time_fit->GetParameter("p3"),time_fit->GetParameter("tau_gdls")) << endl <<  "\t800us: " << efficiency(800,time_fit->GetParameter("p1"),time_fit->GetParameter("tau_ls"),time_fit->GetParameter("p3"),time_fit->GetParameter("tau_gdls")) << endl;

			h_time_efficiency_1us_1m->SetStats(0);
			h_time_efficiency_1us_1m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_1m->GetYaxis()->SetTitle("Efficiency");
				TF1* eff_1=new TF1("eff_1","efficiency(x,[0],[1],[2],[3])",0,2000);
				eff_1->FixParameter(0,time_fit1->GetParameter("p1"));
				eff_1->FixParameter(1,time_fit1->GetParameter("tau_ls"));
				eff_1->FixParameter(2,time_fit1->GetParameter("p3"));
				eff_1->FixParameter(3,time_fit1->GetParameter("tau_gdls"));
		//	h_time_efficiency_1us_1m->Fit("eff_1","R");
			h_time_efficiency_1us_1m->Write();
			cout << "For 1 m... " << endl << "\t400us: " << efficiency(400,time_fit1->GetParameter("p1"),time_fit1->GetParameter("tau_ls"),time_fit1->GetParameter("p3"),time_fit1->GetParameter("tau_gdls")) << endl << "\t600us: " << efficiency(600,time_fit1->GetParameter("p1"),time_fit1->GetParameter("tau_ls"),time_fit1->GetParameter("p3"),time_fit1->GetParameter("tau_gdls")) << endl <<  "\t800us: " << efficiency(800,time_fit1->GetParameter("p1"),time_fit1->GetParameter("tau_ls"),time_fit1->GetParameter("p3"),time_fit1->GetParameter("tau_gdls")) << endl;

			h_time_efficiency_1us_075m->SetStats(0);
			h_time_efficiency_1us_075m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_075m->GetYaxis()->SetTitle("Efficiency");
				TF1* eff_075=new TF1("eff_075","efficiency(x,[0],[1],[2],[3])",0,2000);
				eff_075->FixParameter(0,time_fit075->GetParameter("p1"));
				eff_075->FixParameter(1,time_fit075->GetParameter("tau_ls"));
				eff_075->FixParameter(2,time_fit075->GetParameter("p3"));
				eff_075->FixParameter(3,time_fit075->GetParameter("tau_gdls"));
		//	h_time_efficiency_1us_075m->Fit("eff_075","R");
			h_time_efficiency_1us_075m->Write();
			cout << "For 0.75 m... " << endl << "\t400us: " << efficiency(400,time_fit075->GetParameter("p1"),time_fit075->GetParameter("tau_ls"),time_fit075->GetParameter("p3"),time_fit075->GetParameter("tau_gdls")) << endl << "\t600us: " << efficiency(600,time_fit075->GetParameter("p1"),time_fit075->GetParameter("tau_ls"),time_fit075->GetParameter("p3"),time_fit075->GetParameter("tau_gdls")) << endl <<  "\t800us: " << efficiency(800,time_fit075->GetParameter("p1"),time_fit075->GetParameter("tau_ls"),time_fit075->GetParameter("p3"),time_fit075->GetParameter("tau_gdls")) << endl;

			h_time_efficiency_1us_05m->SetStats(0);
			h_time_efficiency_1us_05m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_05m->GetYaxis()->SetTitle("Efficiency");
				TF1* eff_05=new TF1("eff_05","efficiency(x,[0],[1],[2],[3])",0,2000);
				eff_05->FixParameter(0,time_fit05->GetParameter("p1"));
				eff_05->FixParameter(1,time_fit05->GetParameter("tau_ls"));
				eff_05->FixParameter(2,time_fit05->GetParameter("p3"));
				eff_05->FixParameter(3,time_fit05->GetParameter("tau_gdls"));
		//	h_time_efficiency_1us_05m->Fit("eff_05","R");
			h_time_efficiency_1us_05m->Write();

			cout << "For 0.5 m... " << endl << "\t400us: " << efficiency(400,time_fit05->GetParameter("p1"),time_fit05->GetParameter("tau_ls"),time_fit05->GetParameter("p3"),time_fit05->GetParameter("tau_gdls")) << endl << "\t600us: " << efficiency(600,time_fit05->GetParameter("p1"),time_fit05->GetParameter("tau_ls"),time_fit05->GetParameter("p3"),time_fit05->GetParameter("tau_gdls")) << endl <<  "\t800us: " << efficiency(800,time_fit05->GetParameter("p1"),time_fit05->GetParameter("tau_ls"),time_fit05->GetParameter("p3"),time_fit05->GetParameter("tau_gdls")) << endl;
*/

/*			cout << endl << endl << "no dist cut (LS):" << endl;
			//h_coincidence_time_1us_ls->SetStats(0);
			h_coincidence_time_1us_ls->SetMinimum(0);
			h_coincidence_time_1us_ls->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_ls->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_ls->Fit("time_fit_ls","R");
			h_coincidence_time_1us_ls->Write();
				xBins = h_coincidence_time_1us_ls->GetNbinsX();
				acc_const = time_fit_ls->GetParameter("constant");
				eff_denom = h_coincidence_time_1us_ls->Integral(0,xBins);
				eff_denom = eff_denom - acc_const*xBins;
				for(int iBin = 0; iBin < xBins+1; iBin++){
					eff_numer = 0;
					eff_numer = h_coincidence_time_1us_ls->Integral(0,iBin);
					eff_numer = eff_numer - acc_const*iBin;
					h_time_efficiency_1us_ls->SetBinContent(iBin, eff_numer/eff_denom);
				}
			eff_numer = h_coincidence_time_1us_ls->Integral(0,20);
			eff_numer = eff_numer - acc_const*20;
			cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;		


			cout << endl << endl << "1m cut (LS):" << endl;
			//h_coincidence_time_1us_ls_1m->SetStats(0);
			h_coincidence_time_1us_ls_1m->SetMinimum(0);
			h_coincidence_time_1us_ls_1m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_ls_1m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_ls_1m->Fit("time_fit1_ls","R");
			h_coincidence_time_1us_ls_1m->Write();
				xBins = h_coincidence_time_1us_ls_1m->GetNbinsX();
				acc_const = time_fit1_ls->GetParameter("constant");
				eff_denom = h_coincidence_time_1us_ls_1m->Integral(0,xBins);
				eff_denom = eff_denom - acc_const*xBins;
				for(int iBin = 0; iBin < xBins+1; iBin++){
					eff_numer = 0;
					eff_numer = h_coincidence_time_1us_ls_1m->Integral(0,iBin);
					eff_numer = eff_numer - acc_const*iBin;
					h_time_efficiency_1us_ls_1m->SetBinContent(iBin, eff_numer/eff_denom);
				}
			eff_numer = h_coincidence_time_1us_ls_1m->Integral(0,20);
			eff_numer = eff_numer - acc_const*20;
			cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;



			cout << endl << endl << "0.75m cut (LS):" << endl;
			//h_coincidence_time_1us_ls_075m->SetStats(0);
			h_coincidence_time_1us_ls_075m->SetMinimum(0);
			h_coincidence_time_1us_ls_075m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_ls_075m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_ls_075m->Fit("time_fit075_ls","R");
			h_coincidence_time_1us_ls_075m->Write();
				xBins = h_coincidence_time_1us_ls_075m->GetNbinsX();
				acc_const = time_fit075_ls->GetParameter("constant");
				eff_denom = h_coincidence_time_1us_ls_075m->Integral(0,xBins);
				eff_denom = eff_denom - acc_const*xBins;
				for(int iBin = 0; iBin < xBins+1; iBin++){
					eff_numer = 0;
					eff_numer = h_coincidence_time_1us_ls_075m->Integral(0,iBin);
					eff_numer = eff_numer - acc_const*iBin;
					h_time_efficiency_1us_ls_075m->SetBinContent(iBin, eff_numer/eff_denom);
				}
			eff_numer = h_coincidence_time_1us_ls_075m->Integral(0,20);
			eff_numer = eff_numer - acc_const*20;
			cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;


			cout << endl << endl << "0.5m cut (LS):" << endl;
			//h_coincidence_time_1us_ls_05m->SetStats(0);
			h_coincidence_time_1us_ls_05m->SetMinimum(0);
			h_coincidence_time_1us_ls_05m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_ls_05m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_ls_05m->Fit("time_fit05_ls","R");
			h_coincidence_time_1us_ls_05m->Write();
				xBins = h_coincidence_time_1us_ls_05m->GetNbinsX();
				acc_const = time_fit05_ls->GetParameter("constant");
				eff_denom = h_coincidence_time_1us_ls_05m->Integral(0,xBins);
				eff_denom = eff_denom - acc_const*xBins;
				for(int iBin = 0; iBin < xBins+1; iBin++){
					eff_numer = 0;
					eff_numer = h_coincidence_time_1us_ls_05m->Integral(0,iBin);
					eff_numer = eff_numer - acc_const*iBin;
					h_time_efficiency_1us_ls_05m->SetBinContent(iBin, eff_numer/eff_denom);
				}
			eff_numer = h_coincidence_time_1us_ls_05m->Integral(0,20);
			eff_numer = eff_numer - acc_const*20;
			cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;


			h_time_efficiency_1us_ls->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_ls->GetYaxis()->SetTitle("Efficiency");
			h_time_efficiency_1us_ls->Write();

			h_time_efficiency_1us_ls_1m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_ls_1m->GetYaxis()->SetTitle("Efficiency");
			h_time_efficiency_1us_ls_1m->Write();

			h_time_efficiency_1us_ls_075m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_ls_075m->GetYaxis()->SetTitle("Efficiency");
			h_time_efficiency_1us_ls_075m->Write();

			h_time_efficiency_1us_ls_05m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_ls_05m->GetYaxis()->SetTitle("Efficiency");
			h_time_efficiency_1us_ls_05m->Write();



			cout << endl << endl << "no dist cut (GDLS):" << endl;
			//h_coincidence_time_1us_gdls->SetStats(0);
			h_coincidence_time_1us_gdls->SetMinimum(0);
			h_coincidence_time_1us_gdls->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_gdls->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_gdls->Fit("time_fit_gdls","R");
			h_coincidence_time_1us_gdls->Write();
				xBins = h_coincidence_time_1us_gdls->GetNbinsX();
				acc_const = time_fit_gdls->GetParameter("constant");
				eff_denom = h_coincidence_time_1us_gdls->Integral(0,xBins);
				eff_denom = eff_denom - acc_const*xBins;
				for(int iBin = 0; iBin < xBins+1; iBin++){
					eff_numer = 0;
					eff_numer = h_coincidence_time_1us_gdls->Integral(0,iBin);
					eff_numer = eff_numer - acc_const*iBin;
					h_time_efficiency_1us_gdls->SetBinContent(iBin, eff_numer/eff_denom);
				}
			eff_numer = h_coincidence_time_1us_gdls->Integral(0,20);
			eff_numer = eff_numer - acc_const*20;
			cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;



			cout << endl << endl << "1m cut (GDLS):" << endl;
			//h_coincidence_time_1us_gdls_1m->SetStats(0);
			h_coincidence_time_1us_gdls_1m->SetMinimum(0);
			h_coincidence_time_1us_gdls_1m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_gdls_1m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_gdls_1m->Fit("time_fit1_gdls","R");
			h_coincidence_time_1us_gdls_1m->Write();
				xBins = h_coincidence_time_1us_gdls_1m->GetNbinsX();
				acc_const = time_fit1_gdls->GetParameter("constant");
				eff_denom = h_coincidence_time_1us_gdls_1m->Integral(0,xBins);
				eff_denom = eff_denom - acc_const*xBins;
				for(int iBin = 0; iBin < xBins+1; iBin++){
					eff_numer = 0;
					eff_numer = h_coincidence_time_1us_gdls_1m->Integral(0,iBin);
					eff_numer = eff_numer - acc_const*iBin;
					h_time_efficiency_1us_gdls_1m->SetBinContent(iBin, eff_numer/eff_denom);
				}
			eff_numer = h_coincidence_time_1us_gdls_1m->Integral(0,20);
			eff_numer = eff_numer - acc_const*20;
			cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;



			cout << endl << endl << "0.75m cut (GDLS):" << endl;
			//h_coincidence_time_1us_gdls_075m->SetStats(0);
			h_coincidence_time_1us_gdls_075m->SetMinimum(0);
			h_coincidence_time_1us_gdls_075m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_gdls_075m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_gdls_075m->Fit("time_fit075_gdls","R");
			h_coincidence_time_1us_gdls_075m->Write();
				xBins = h_coincidence_time_1us_gdls_075m->GetNbinsX();
				acc_const = time_fit075_gdls->GetParameter("constant");
				eff_denom = h_coincidence_time_1us_gdls_075m->Integral(0,xBins);
				eff_denom = eff_denom - acc_const*xBins;
				for(int iBin = 0; iBin < xBins+1; iBin++){
					eff_numer = 0;
					eff_numer = h_coincidence_time_1us_gdls_075m->Integral(0,iBin);
					eff_numer = eff_numer - acc_const*iBin;
					h_time_efficiency_1us_gdls_075m->SetBinContent(iBin, eff_numer/eff_denom);
				}
			eff_numer = h_coincidence_time_1us_gdls_075m->Integral(0,20);
			eff_numer = eff_numer - acc_const*20;
			cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;


			cout << endl << endl << "0.5m cut (GDLS):" << endl;
			//h_coincidence_time_1us_gdls_05m->SetStats(0);
			h_coincidence_time_1us_gdls_05m->SetMinimum(0);
			h_coincidence_time_1us_gdls_05m->GetXaxis()->SetTitle("Delta Time [us]");
			h_coincidence_time_1us_gdls_05m->GetYaxis()->SetTitle("Counts");
			h_coincidence_time_1us_gdls_05m->Fit("time_fit05_gdls","R");
			h_coincidence_time_1us_gdls_05m->Write();
				xBins = h_coincidence_time_1us_gdls_05m->GetNbinsX();
				acc_const = time_fit05_gdls->GetParameter("constant");
				eff_denom = h_coincidence_time_1us_gdls_05m->Integral(0,xBins);
				eff_denom = eff_denom - acc_const*xBins;
				for(int iBin = 0; iBin < xBins+1; iBin++){
					eff_numer = 0;
					eff_numer = h_coincidence_time_1us_gdls_05m->Integral(0,iBin);
					eff_numer = eff_numer - acc_const*iBin;
					h_time_efficiency_1us_gdls_05m->SetBinContent(iBin, eff_numer/eff_denom);
				}
			eff_numer = h_coincidence_time_1us_gdls_05m->Integral(0,20);
			eff_numer = eff_numer - acc_const*20;
			cout << "Efficiency at 400 us: " << eff_numer/eff_denom << endl;


			h_time_efficiency_1us_gdls->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_gdls->GetYaxis()->SetTitle("Efficiency");
			h_time_efficiency_1us_gdls->Write();

			h_time_efficiency_1us_gdls_1m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_gdls_1m->GetYaxis()->SetTitle("Efficiency");
			h_time_efficiency_1us_gdls_1m->Write();

			h_time_efficiency_1us_gdls_075m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_gdls_075m->GetYaxis()->SetTitle("Efficiency");
			h_time_efficiency_1us_gdls_075m->Write();

			h_time_efficiency_1us_gdls_05m->GetXaxis()->SetTitle("Delta Time [us]");
			h_time_efficiency_1us_gdls_05m->GetYaxis()->SetTitle("Efficiency");
			h_time_efficiency_1us_gdls_05m->Write();

*/



	TCanvas *c1 = new TCanvas("c1","Distance Plots");
	c1->cd();
	h_IBDDistance->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
	h_IBDDistance->GetYaxis()->SetTitle("Counts");
	h_IBDDistance->SetLineColor(kBlack);
	h_IBDDistance->SetLineWidth(3);
	h_accDistance_norm->SetLineColor(kBlue);
	h_accDistance_norm->SetLineWidth(3);
	h_acc_scaled->SetLineColor(kRed);
	h_acc_scaled->SetLineWidth(3);
	h_IBDDistance->Draw();
	h_accDistance_norm->Draw("same");
	h_acc_scaled->Draw("same");


	//400 subset
	TCanvas *c4 = new TCanvas("c4","Distance Plots_400us");
	c4->cd();
	h_IBD_distance_400->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
	h_IBD_distance_400->GetYaxis()->SetTitle("Counts");
	h_IBD_distance_400->SetLineColor(kBlack);
	h_IBD_distance_400->SetLineWidth(3);
	h_acc_norm_400->SetLineColor(kBlue);
	h_acc_norm_400->SetLineWidth(3);
	h_acc_scaled_400->SetLineColor(kRed);
	h_acc_scaled_400->SetLineWidth(3);
	h_IBD_distance_400->Draw();
	h_acc_norm_400->Draw("same");
	h_acc_scaled_400->Draw("same");

	//600 subset
	TCanvas *c6 = new TCanvas("c6","Distance Plots_600us");
	c6->cd();
	h_IBD_distance_600->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
	h_IBD_distance_600->GetYaxis()->SetTitle("Counts");
	h_IBD_distance_600->SetLineColor(kBlack);
	h_IBD_distance_600->SetLineWidth(3);
	h_acc_norm_600->SetLineColor(kBlue);
	h_acc_norm_600->SetLineWidth(3);
	h_acc_scaled_600->SetLineColor(kRed);
	h_acc_scaled_600->SetLineWidth(3);
	h_IBD_distance_600->Draw();
	h_acc_norm_600->Draw("same");
	h_acc_scaled_600->Draw("same");

	//800 subset
	TCanvas *c8 = new TCanvas("c8","Distance Plots_800us");
	c8->cd();
	h_IBD_distance_800->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
	h_IBD_distance_800->GetYaxis()->SetTitle("Counts");
	h_IBD_distance_800->SetLineColor(kBlack);
	h_IBD_distance_800->SetLineWidth(3);
	h_acc_norm_800->SetLineColor(kBlue);
	h_acc_norm_800->SetLineWidth(3);
	h_acc_scaled_800->SetLineColor(kRed);
	h_acc_scaled_800->SetLineWidth(3);
	h_IBD_distance_800->Draw();
	h_acc_norm_800->Draw("same");
	h_acc_scaled_800->Draw("same");

        char title[80];

	//DT distribution
	TCanvas *cDT = new TCanvas("cDT","DT plots");
	cDT->cd();
	h_ibd_DT->SetStats(0);
	h_ibd_DT->SetLineColor(kBlack);
	h_acc_DT_DTnorm->SetLineColor(kRed);
	h_ibd_DT->SetLineWidth(3);
	h_acc_DT_DTnorm->SetLineWidth(3);
	h_ibd_DT->Draw();
	h_acc_DT_DTnorm->Draw("same");

	cDT->BuildLegend();
	sprintf(title,"../nH_files/DTplots_EH%dAD%d.png",hall_num,ad_num);
	cDT->Print(title);

	TCanvas *cDT_Ep = new TCanvas("cDT_Ep","DT plots_Ep");
	cDT_Ep->cd();
	h_ibd_DT_Ep35->SetStats(0);
	h_ibd_DT_Ep35->SetLineColor(kBlack);
	h_acc_DT_Ep35_DTnorm->SetLineColor(kRed);
	h_ibd_DT_Ep35->SetLineWidth(3);
	h_acc_DT_Ep35_DTnorm->SetLineWidth(3);
	h_ibd_DT_Ep35->Draw();
	h_acc_DT_Ep35_DTnorm->Draw("same");

	cDT_Ep->BuildLegend();
	sprintf(title,"../nH_files/DTplots_Ep_EH%dAD%d.png",hall_num,ad_num);
	cDT_Ep->Print(title);

	TCanvas *cPrompt = new TCanvas("cPrompt","Prompt");
	cPrompt->cd();
	h_IBD_Eprompt->SetStats(0);
	h_IBD_Eprompt->SetLineColor(kBlack);
	h_acc_Eprompt_norm->SetLineColor(kRed);
	h_IBD_Eprompt->SetLineWidth(3);
	h_acc_Eprompt_norm->SetLineWidth(3);
	h_IBD_Eprompt->Draw();
	h_acc_Eprompt_norm->Draw("same");

	cPrompt->BuildLegend();

	//After these checks, get rid of the following parts:

	TCanvas *cEnergy = new TCanvas("cEnergy","Energy");
	cEnergy->cd();
		//Prompt subtracted spectra w/o DT cut - Prompt subtracted spectra w/ DT cut
	h_Eprompt_sub_DT800_norm->Rebin(2);
	h_Eprompt_sub_norm->Rebin(2);
	h_Eprompt_sub_DT800_norm->Scale((h_Eprompt_sub_norm->Integral(h_Eprompt_sub_norm->FindBin(3.5),h_Eprompt_sub_norm->FindBin(12.)))/(h_Eprompt_sub_DT800_norm->Integral(h_Eprompt_sub_DT800_norm->FindBin(3.5),h_Eprompt_sub_DT800_norm->FindBin(12.))));
	h_Eprompt_sub_norm->SetStats(0);
	h_Eprompt_sub_norm->GetXaxis()->SetTitle("Energy [MeV]");
	h_Eprompt_sub_norm->GetYaxis()->SetTitle("Counts");
	h_Eprompt_sub_norm->SetLineColor(kBlack);
	h_Eprompt_sub_DT800_norm->SetLineColor(kRed);
	h_Eprompt_sub_norm->SetLineWidth(2);
	h_Eprompt_sub_DT800_norm->SetLineWidth(2);
	h_Eprompt_sub_norm->Draw();
	h_Eprompt_sub_DT800_norm->Draw("same");

	//cEnergy->BuildLegend();
	sprintf(title,"../nH_files/PromptSub_compNoDTwDT800_EH%dAD%d.png",hall_num,ad_num);
	cEnergy->Print(title);




		//Prompt to Delayed Ratio:
	TCanvas *cpVSd = new TCanvas("cpVSd","cpVSd");
	cpVSd->cd();
	h_acc_Edelayed_norm->Rebin(10);

		double p_counts = 0;
		double d_counts = 0;
		double p_binCenter = 0;
		double d_binCenter = 0;
		for(int i = 1; i <= (h_acc_Edelayed_norm->GetNbinsX()); i++){
			if(h_acc_Edelayed_norm->GetBinCenter(i) <1.5) continue;
			d_counts = h_acc_Edelayed_norm->GetBinContent(i);
			//cout << h_acc_Edelayed_norm->GetBinCenter(i) << "\t" << d_counts << endl;
			d_binCenter = h_acc_Edelayed_norm->GetBinCenter(i);
			p_binCenter = h_acc_Eprompt_norm->GetBinCenter(h_acc_Eprompt_norm->FindBin(d_binCenter));
			//if(d_binCenter != p_binCenter) cout << "Bin Centers do not match!!!" << endl;
			if(d_counts == 0) continue;
			p_counts = h_acc_Eprompt_norm->GetBinContent(h_acc_Edelayed_norm->FindBin(h_acc_Edelayed_norm->GetBinCenter(i)));
			h_EpromptToDelayed_norm->SetBinContent(h_acc_Edelayed_norm->FindBin(h_acc_Edelayed_norm->GetBinCenter(i)), p_counts/d_counts);
		}

	h_EpromptToDelayed_norm->SetStats(0);
	h_EpromptToDelayed_norm->SetLineWidth(1);
	h_EpromptToDelayed_norm->GetXaxis()->SetRangeUser(1.5,3.);
	h_EpromptToDelayed_norm->GetXaxis()->SetTitle("Energy [MeV]");
	h_EpromptToDelayed_norm->GetYaxis()->SetTitle("Prompt/Delayed");
	h_EpromptToDelayed_norm->Draw();
	//h_acc_Edelayed_norm->Draw();
	//h_acc_Eprompt_norm->Draw("same");

	//cpVSd->BuildLegend();
	sprintf(title,"../nH_files/promptTOdelayed_EH%dAD%d.png",hall_num,ad_num);
	cpVSd->Print(title);

	TCanvas *d_rate = new TCanvas("d_rate","d_rate");
	d_rate->cd();
	h_Edelayed_sub->Draw();
	sprintf(title,"../nH_files/Edelayed_rate_fit_EH%dAD%d.png",hall_num,ad_num);
	d_rate->Print(title);

	TCanvas *d_norm = new TCanvas("d_norm","d_norm");
	d_norm->cd();
	h_Edelayed_sub_norm->Draw();
	sprintf(title,"../nH_files/Edelayed_norm_fit_EH%dAD%d.png",hall_num,ad_num);
	d_norm->Print(title);

	TCanvas *d_DTnorm = new TCanvas("d_DTnorm","d_DTnorm");
	d_DTnorm->cd();
	h_Edelayed_sub_DTnorm->Draw();
	sprintf(title,"../nH_files/Edelayed_DTnorm_fit_EH%dAD%d.png",hall_num,ad_num);
	d_DTnorm->Print(title);

	TCanvas *d_DT800_rate = new TCanvas("d_DT800_rate","d_DT800_rate");
	d_DT800_rate->cd();
//	h_Edelayed_sub_DT800->Rebin(2);
//	h_Edelayed_sub_DT800->GetXaxis()->SetRangeUser(1.5,3);
	h_Edelayed_sub_DT800->Draw();
	sprintf(title,"../nH_files/Edelayed_DT800_rate_fit_EH%dAD%d.png",hall_num,ad_num);
	d_DT800_rate->Print(title);

	TCanvas *d_DT800_norm = new TCanvas("d_DT800_norm","d_DT800_norm");
	d_DT800_norm->cd();
	h_Edelayed_sub_DT800_norm->Draw();
	sprintf(title,"../nH_files/Edelayed_DT800_norm_fit_EH%dAD%d.png",hall_num,ad_num);
	d_DT800_norm->Print(title);

	TCanvas *d_DT800_DTnorm = new TCanvas("d_DT800_DTnorm","d_DT800_DTnorm");
	d_DT800_DTnorm->cd();
	h_Edelayed_sub_DT800_DTnorm->Draw();
//	delayedFit_DT800_DTnorm->Draw("same");
//	delayedFit_DT800_DTnorm_Sam->Draw("same");
	sprintf(title,"../nH_files/Edelayed_DT800_DTnorm_fit_EH%dAD%d.png",hall_num,ad_num);
	d_DT800_DTnorm->Print(title);


/*	TCanvas *d1 = new TCanvas("d1","Ratio of Distance");
	d1->cd();
	h_dist_ratio_rateTOnorm->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
	h_dist_ratio_rateTOnorm->GetYaxis()->SetTitle("Ratio of Rate to Normalized");
	h_dist_ratio_rateTOnorm->Draw("HIST");

	TCanvas *d2 = new TCanvas("d2","Sub of Acc Distance"); //Come back to this part
	d2->cd();
	h_dist_sub_rateTOnorm->GetXaxis()->SetTitle("Distance Between Prompt and Delayed [m]");
	h_dist_sub_rateTOnorm->GetYaxis()->SetTitle("Rate minus Normalized");
	h_dist_sub_rateTOnorm->Draw("HIST");


	TCanvas *c4 = new TCanvas("c4","Prompt Ratios and Subtractions for Distance Slices",3000,1000);
	c4->Divide(6,3);
		c4->cd(1);
		h_ibd_Penergy_dist[0]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ibd_Penergy_dist[0]->GetYaxis()->SetTitle("Counts");
		h_ibd_Penergy_dist[0]->Draw();
		h_acc_Penergy_dist[0]->SetLineColor(kRed);
		h_acc_Penergy_dist[0]->Draw("same");
		h_acc_Penergy_dist_norm[0]->SetLineColor(kBlue);
		h_acc_Penergy_dist_norm[0]->Draw("same");

		c4->cd(6+1);
		h_ratio_prompt_energy_scaled_dist[0]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaled_dist[0]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_prompt_energy_scaled_dist[0]->SetLineColor(kRed);
		h_ratio_prompt_energy_scaled_dist[0]->Draw();
		h_ratio_prompt_energy_norm_dist[0]->SetLineColor(kBlue);
		h_ratio_prompt_energy_norm_dist[0]->Draw("same");

		c4->cd(6+6+1);
		h_sub_prompt_energy_scaled_dist[0]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_sub_prompt_energy_scaled_dist[0]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_prompt_energy_scaled_dist[0]->SetLineColor(kRed);
		h_sub_prompt_energy_norm_dist[0]->SetLineColor(kBlue);
		h_sub_prompt_energy_norm_dist[0]->Draw();
		h_sub_prompt_energy_scaled_dist[0]->Draw("same");

		c4->cd(2);
		h_ibd_Penergy_dist[1]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ibd_Penergy_dist[1]->GetYaxis()->SetTitle("Counts");
		h_ibd_Penergy_dist[1]->Draw();
		h_acc_Penergy_dist[1]->SetLineColor(kRed);
		h_acc_Penergy_dist[1]->Draw("same");
		h_acc_Penergy_dist_norm[1]->SetLineColor(kBlue);
		h_acc_Penergy_dist_norm[1]->Draw("same");

		c4->cd(6+2);
		h_ratio_prompt_energy_scaled_dist[1]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaled_dist[1]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_prompt_energy_scaled_dist[1]->SetLineColor(kRed);
		h_ratio_prompt_energy_scaled_dist[1]->Draw();
		h_ratio_prompt_energy_norm_dist[1]->SetLineColor(kBlue);
		h_ratio_prompt_energy_norm_dist[1]->Draw("same");

		c4->cd(6+6+2);
		h_sub_prompt_energy_scaled_dist[1]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_sub_prompt_energy_scaled_dist[1]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_prompt_energy_scaled_dist[1]->SetLineColor(kRed);
		h_sub_prompt_energy_scaled_dist[1]->Draw();
		h_sub_prompt_energy_norm_dist[1]->SetLineColor(kBlue);
		h_sub_prompt_energy_norm_dist[1]->Draw("same");


		c4->cd(3);
		h_ibd_Penergy_dist[2]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ibd_Penergy_dist[2]->GetYaxis()->SetTitle("Counts");
		h_ibd_Penergy_dist[2]->Draw();
		h_acc_Penergy_dist[2]->SetLineColor(kRed);
		h_acc_Penergy_dist[2]->Draw("same");
		h_acc_Penergy_dist_norm[2]->SetLineColor(kBlue);
		h_acc_Penergy_dist_norm[2]->Draw("same");

		c4->cd(6+3);
		h_ratio_prompt_energy_scaled_dist[2]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaled_dist[2]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_prompt_energy_scaled_dist[2]->SetLineColor(kRed);
		h_ratio_prompt_energy_scaled_dist[2]->Draw();
		h_ratio_prompt_energy_norm_dist[2]->SetLineColor(kBlue);
		h_ratio_prompt_energy_norm_dist[2]->Draw("same");

		c4->cd(6+6+3);
		h_sub_prompt_energy_scaled_dist[2]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_sub_prompt_energy_scaled_dist[2]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_prompt_energy_scaled_dist[2]->SetLineColor(kRed);
		h_sub_prompt_energy_scaled_dist[2]->Draw();
		h_sub_prompt_energy_norm_dist[2]->SetLineColor(kBlue);
		h_sub_prompt_energy_norm_dist[2]->Draw("same");

		c4->cd(4);
		h_ibd_Penergy_dist[3]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ibd_Penergy_dist[3]->GetYaxis()->SetTitle("Counts");
		h_ibd_Penergy_dist[3]->Draw();
		h_acc_Penergy_dist[3]->SetLineColor(kRed);
		h_acc_Penergy_dist[3]->Draw("same");
		h_acc_Penergy_dist_norm[3]->SetLineColor(kBlue);
		h_acc_Penergy_dist_norm[3]->Draw("same");

		c4->cd(6+4);
		h_ratio_prompt_energy_scaled_dist[3]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaled_dist[3]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_prompt_energy_scaled_dist[3]->SetLineColor(kRed);
		h_ratio_prompt_energy_scaled_dist[3]->Draw();
		h_ratio_prompt_energy_norm_dist[3]->SetLineColor(kBlue);
		h_ratio_prompt_energy_norm_dist[3]->Draw("same");

		c4->cd(6+6+4);
		h_sub_prompt_energy_scaled_dist[3]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_sub_prompt_energy_scaled_dist[3]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_prompt_energy_scaled_dist[3]->SetLineColor(kRed);
		h_sub_prompt_energy_scaled_dist[3]->Draw();
		h_sub_prompt_energy_norm_dist[3]->SetLineColor(kBlue);
		h_sub_prompt_energy_norm_dist[3]->Draw("same");


		c4->cd(5);
		h_ibd_Penergy_dist[4]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ibd_Penergy_dist[4]->GetYaxis()->SetTitle("Counts");
		h_ibd_Penergy_dist[4]->Draw();
		h_acc_Penergy_dist[4]->SetLineColor(kRed);
		h_acc_Penergy_dist[4]->Draw("same");
		h_acc_Penergy_dist_norm[4]->SetLineColor(kBlue);
		h_acc_Penergy_dist_norm[4]->Draw("same");

		c4->cd(6+5);
		h_ratio_prompt_energy_scaled_dist[4]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaled_dist[4]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_prompt_energy_scaled_dist[4]->SetLineColor(kRed);
		h_ratio_prompt_energy_scaled_dist[4]->Draw();
		h_ratio_prompt_energy_norm_dist[4]->SetLineColor(kBlue);
		h_ratio_prompt_energy_norm_dist[4]->Draw("same");

		c4->cd(6+6+5);
		h_sub_prompt_energy_scaled_dist[4]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_sub_prompt_energy_scaled_dist[4]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_prompt_energy_scaled_dist[4]->SetLineColor(kRed);
		h_sub_prompt_energy_scaled_dist[4]->Draw();
		h_sub_prompt_energy_norm_dist[4]->SetLineColor(kBlue);
		h_sub_prompt_energy_norm_dist[4]->Draw("same");

		c4->cd(6);
		h_ibd_Penergy_dist[5]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ibd_Penergy_dist[5]->GetYaxis()->SetTitle("Counts");
		h_ibd_Penergy_dist[5]->Draw();
		h_acc_Penergy_dist[5]->SetLineColor(kRed);
		h_acc_Penergy_dist[5]->Draw("same");
		h_acc_Penergy_dist_norm[5]->SetLineColor(kBlue);
		h_acc_Penergy_dist_norm[5]->Draw("same");

		c4->cd(6+6);
		h_ratio_prompt_energy_scaled_dist[5]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaled_dist[5]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_prompt_energy_scaled_dist[5]->SetLineColor(kRed);
		h_ratio_prompt_energy_scaled_dist[5]->Draw();
		h_ratio_prompt_energy_norm_dist[5]->SetLineColor(kBlue);
		h_ratio_prompt_energy_norm_dist[5]->Draw("same");

		c4->cd(6+6+6);
		h_sub_prompt_energy_scaled_dist[5]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_sub_prompt_energy_scaled_dist[5]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_prompt_energy_scaled_dist[5]->SetLineColor(kRed);
		h_sub_prompt_energy_scaled_dist[5]->Draw();
		h_sub_prompt_energy_norm_dist[5]->SetLineColor(kBlue);
		h_sub_prompt_energy_norm_dist[5]->Draw("same");



	TCanvas *c5 = new TCanvas("c5","Prompt and Delayed Ratios (Rate to Norm) for Distance Slices",2500,500);
	c5->Divide(6,2);
		c5->cd(1);
		h_ratio_prompt_energy_scaledTOnorm_dist[0]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaledTOnorm_dist[0]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_prompt_energy_scaledTOnorm_dist[0]->Draw();

		c5->cd(6+1);
		h_ratio_delayed_energy_scaledTOnorm_dist[0]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaledTOnorm_dist[0]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_delayed_energy_scaledTOnorm_dist[0]->Draw();


		c5->cd(2);
		h_ratio_prompt_energy_scaledTOnorm_dist[1]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaledTOnorm_dist[1]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_prompt_energy_scaledTOnorm_dist[1]->Draw();

		c5->cd(6+2);
		h_ratio_delayed_energy_scaledTOnorm_dist[1]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaledTOnorm_dist[1]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_delayed_energy_scaledTOnorm_dist[1]->Draw();

		c5->cd(3);
		h_ratio_prompt_energy_scaledTOnorm_dist[2]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaledTOnorm_dist[2]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_prompt_energy_scaledTOnorm_dist[2]->Draw();

		c5->cd(6+3);
		h_ratio_delayed_energy_scaledTOnorm_dist[2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaledTOnorm_dist[2]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_delayed_energy_scaledTOnorm_dist[2]->Draw();

		c5->cd(4);
		h_ratio_prompt_energy_scaledTOnorm_dist[3]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaledTOnorm_dist[3]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_prompt_energy_scaledTOnorm_dist[3]->Draw();

		c5->cd(6+4);
		h_ratio_delayed_energy_scaledTOnorm_dist[3]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaledTOnorm_dist[3]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_delayed_energy_scaledTOnorm_dist[3]->Draw();

		c5->cd(5);
		h_ratio_prompt_energy_scaledTOnorm_dist[4]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaledTOnorm_dist[4]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_prompt_energy_scaledTOnorm_dist[4]->Draw();

		c5->cd(6+5);
		h_ratio_delayed_energy_scaledTOnorm_dist[4]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaledTOnorm_dist[4]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_delayed_energy_scaledTOnorm_dist[4]->Draw();

		c5->cd(6);
		h_ratio_prompt_energy_scaledTOnorm_dist[5]->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		h_ratio_prompt_energy_scaledTOnorm_dist[5]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_prompt_energy_scaledTOnorm_dist[5]->Draw();

		c5->cd(6+6);
		h_ratio_delayed_energy_scaledTOnorm_dist[5]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaledTOnorm_dist[5]->GetYaxis()->SetTitle("Scaled/Norm");
		h_ratio_delayed_energy_scaledTOnorm_dist[5]->Draw();


	TCanvas *c6 = new TCanvas("c6","Delayed Ratios and Subtractions for Distance Slices",3000,1000);
	c6->Divide(6,3);
		c6->cd(1);
		h_ibd_Denergy_dist[0]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ibd_Denergy_dist[0]->GetYaxis()->SetTitle("Counts");
		h_ibd_Denergy_dist[0]->Draw();
		h_acc_Denergy_dist[0]->SetLineColor(kRed);
		h_acc_Denergy_dist[0]->Draw("same");
		h_acc_Denergy_dist_norm[0]->SetLineColor(kBlue);
		h_acc_Denergy_dist_norm[0]->Draw("same");

		c6->cd(6+1);
		h_ratio_delayed_energy_scaled_dist[0]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaled_dist[0]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_delayed_energy_scaled_dist[0]->SetLineColor(kRed);
		h_ratio_delayed_energy_scaled_dist[0]->Draw();
		h_ratio_delayed_energy_norm_dist[0]->SetLineColor(kBlue);
		h_ratio_delayed_energy_norm_dist[0]->Draw("same");

		c6->cd(6+6+1);
		h_sub_delayed_energy_scaled_dist[0]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_sub_delayed_energy_scaled_dist[0]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_delayed_energy_scaled_dist[0]->SetLineColor(kRed);
		h_sub_delayed_energy_norm_dist[0]->SetLineColor(kBlue);
		h_sub_delayed_energy_norm_dist[0]->Draw();
		h_sub_delayed_energy_scaled_dist[0]->Draw("same");

		c6->cd(2);
		h_ibd_Denergy_dist[1]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ibd_Denergy_dist[1]->GetYaxis()->SetTitle("Counts");
		h_ibd_Denergy_dist[1]->Draw();
		h_acc_Denergy_dist[1]->SetLineColor(kRed);
		h_acc_Denergy_dist[1]->Draw("same");
		h_acc_Denergy_dist_norm[1]->SetLineColor(kBlue);
		h_acc_Denergy_dist_norm[1]->Draw("same");

		c6->cd(6+2);
		h_ratio_delayed_energy_scaled_dist[1]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaled_dist[1]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_delayed_energy_scaled_dist[1]->SetLineColor(kRed);
		h_ratio_delayed_energy_scaled_dist[1]->Draw();
		h_ratio_delayed_energy_norm_dist[1]->SetLineColor(kBlue);
		h_ratio_delayed_energy_norm_dist[1]->Draw("same");


		c6->cd(6+6+2);
		h_sub_delayed_energy_scaled_dist[1]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_sub_delayed_energy_scaled_dist[1]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_delayed_energy_scaled_dist[1]->SetLineColor(kRed);
		h_sub_delayed_energy_scaled_dist[1]->Draw();
		h_sub_delayed_energy_norm_dist[1]->SetLineColor(kBlue);
		h_sub_delayed_energy_norm_dist[1]->Draw("same");

		c6->cd(3);
		h_ibd_Denergy_dist[2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ibd_Denergy_dist[2]->GetYaxis()->SetTitle("Counts");
		h_ibd_Denergy_dist[2]->Draw();
		h_acc_Denergy_dist[2]->SetLineColor(kRed);
		h_acc_Denergy_dist[2]->Draw("same");
		h_acc_Denergy_dist_norm[2]->SetLineColor(kBlue);
		h_acc_Denergy_dist_norm[2]->Draw("same");

		c6->cd(6+3);
		h_ratio_delayed_energy_scaled_dist[2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaled_dist[2]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_delayed_energy_scaled_dist[2]->SetLineColor(kRed);
		h_ratio_delayed_energy_scaled_dist[2]->Draw();
		h_ratio_delayed_energy_norm_dist[2]->SetLineColor(kBlue);
		h_ratio_delayed_energy_norm_dist[2]->Draw("same");

		c6->cd(6+6+3);
		h_sub_delayed_energy_scaled_dist[2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_sub_delayed_energy_scaled_dist[2]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_delayed_energy_scaled_dist[2]->SetLineColor(kRed);
		h_sub_delayed_energy_scaled_dist[2]->Draw();
		h_sub_delayed_energy_norm_dist[2]->SetLineColor(kBlue);
		h_sub_delayed_energy_norm_dist[2]->Draw("same");

		c6->cd(4);
		h_ibd_Denergy_dist[3]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ibd_Denergy_dist[3]->GetYaxis()->SetTitle("Counts");
		h_ibd_Denergy_dist[3]->Draw();
		h_acc_Denergy_dist[3]->SetLineColor(kRed);
		h_acc_Denergy_dist[3]->Draw("same");
		h_acc_Denergy_dist_norm[3]->SetLineColor(kBlue);
		h_acc_Denergy_dist_norm[3]->Draw("same");

		c6->cd(6+4);
		h_ratio_delayed_energy_scaled_dist[3]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaled_dist[3]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_delayed_energy_scaled_dist[3]->SetLineColor(kRed);
		h_ratio_delayed_energy_scaled_dist[3]->Draw();
		h_ratio_delayed_energy_norm_dist[3]->SetLineColor(kBlue);
		h_ratio_delayed_energy_norm_dist[3]->Draw("same");

		c6->cd(6+6+4);
		h_sub_delayed_energy_scaled_dist[3]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_sub_delayed_energy_scaled_dist[3]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_delayed_energy_scaled_dist[3]->SetLineColor(kRed);
		h_sub_delayed_energy_scaled_dist[3]->Draw();
		h_sub_delayed_energy_norm_dist[3]->SetLineColor(kBlue);
		h_sub_delayed_energy_norm_dist[3]->Draw("same");

		c6->cd(5);
		h_ibd_Denergy_dist[4]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ibd_Denergy_dist[4]->GetYaxis()->SetTitle("Counts");
		h_ibd_Denergy_dist[4]->Draw();
		h_acc_Denergy_dist[4]->SetLineColor(kRed);
		h_acc_Denergy_dist[4]->Draw("same");
		h_acc_Denergy_dist_norm[4]->SetLineColor(kBlue);
		h_acc_Denergy_dist_norm[4]->Draw("same");

		c6->cd(6+5);
		h_ratio_delayed_energy_scaled_dist[4]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaled_dist[4]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_delayed_energy_scaled_dist[4]->SetLineColor(kRed);
		h_ratio_delayed_energy_scaled_dist[4]->Draw();
		h_ratio_delayed_energy_norm_dist[4]->SetLineColor(kBlue);
		h_ratio_delayed_energy_norm_dist[4]->Draw("same");

		c6->cd(6+6+5);
		h_sub_delayed_energy_scaled_dist[4]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_sub_delayed_energy_scaled_dist[4]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_delayed_energy_scaled_dist[4]->SetLineColor(kRed);
		h_sub_delayed_energy_scaled_dist[4]->Draw();
		h_sub_delayed_energy_norm_dist[4]->SetLineColor(kBlue);
		h_sub_delayed_energy_norm_dist[4]->Draw("same");

		c6->cd(6);
		h_ibd_Denergy_dist[5]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ibd_Denergy_dist[5]->GetYaxis()->SetTitle("Counts");
		h_ibd_Denergy_dist[5]->Draw();
		h_acc_Denergy_dist[5]->SetLineColor(kRed);
		h_acc_Denergy_dist[5]->Draw("same");
		h_acc_Denergy_dist_norm[5]->SetLineColor(kBlue);
		h_acc_Denergy_dist_norm[5]->Draw("same");

		c6->cd(6+6);
		h_ratio_delayed_energy_scaled_dist[5]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_ratio_delayed_energy_scaled_dist[5]->GetYaxis()->SetTitle("IBD/Acc");
		h_ratio_delayed_energy_scaled_dist[5]->SetLineColor(kRed);
		h_ratio_delayed_energy_scaled_dist[5]->Draw();
		h_ratio_delayed_energy_norm_dist[5]->SetLineColor(kBlue);
		h_ratio_delayed_energy_norm_dist[5]->Draw("same");

		c6->cd(6+6+6);
		h_sub_delayed_energy_scaled_dist[5]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_sub_delayed_energy_scaled_dist[5]->GetYaxis()->SetTitle("IBD-Acc");
		h_sub_delayed_energy_scaled_dist[5]->SetLineColor(kRed);
		h_sub_delayed_energy_scaled_dist[5]->Draw();
		h_sub_delayed_energy_norm_dist[5]->SetLineColor(kBlue);
		h_sub_delayed_energy_norm_dist[5]->Draw("same");
*/


/*	TCanvas *s1 = new TCanvas("s1","Scales per run");
	s1->cd();
		h_rateScale->SetStats(0);
		h_normScale->SetStats(0);
		h_rateScale->GetXaxis()->SetTitle("Number of Runs (Since Start of P17B)");
		h_rateScale->GetYaxis()->SetTitle("Scale");
		h_rateScale->SetLineColor(kRed);
		h_normScale->SetLineColor(kBlue);
		h_normScale->Draw();
		h_rateScale->Draw("same");*/
/*
	TCanvas *p1 = new TCanvas("p1","Prompt 1 m cut vs no cut");
	p1->cd();
		h_Eprompt_sub->Scale(1/(h_Eprompt_sub->Integral()));
		h_Eprompt_sub->GetXaxis()->SetTitle("Energy [MeV]");
		h_Eprompt_sub->GetYaxis()->SetTitle("Counts");
		h_Eprompt_sub->SetLineColor(kBlack);
		h_Eprompt_sub->Draw("HIST");
		h_sub_prompt_energy_scaled_dist[0]->Scale(1/(h_sub_prompt_energy_scaled_dist[0]->Integral()));
		h_sub_prompt_energy_scaled_dist[0]->SetLineColor(kRed);
		h_sub_prompt_energy_scaled_dist[0]->Draw("same");

	TCanvas *p2 = new TCanvas("p2","Delayed 1 m cut vs no cut");
	p2->cd();
		h_Edelayed_sub->Scale(1/(h_Edelayed_sub->Integral()));
		h_Edelayed_sub->GetXaxis()->SetTitle("Energy [MeV]");
		h_Edelayed_sub->GetYaxis()->SetTitle("Counts");
		h_Edelayed_sub->SetLineColor(kBlack);
		h_Edelayed_sub->Draw("HIST");
		h_sub_delayed_energy_scaled_dist[0]->Scale(1/(h_sub_delayed_energy_scaled_dist[0]->Integral()));
		h_sub_delayed_energy_scaled_dist[0]->SetLineColor(kRed);
		h_sub_delayed_energy_scaled_dist[0]->Draw("same");
		
*/

/*	TCanvas *test = new TCanvas("test","test");//BIG FAT TEST
		test->cd();
		h_Edelayed_IBD_fine->SetLineColor(kBlack);
		h_Edelayed_IBD_fine->Draw();
		h_acc_Edelayed_fine_scaled->SetLineColor(kRed);
		h_acc_Edelayed_fine_scaled->Draw("same");*/


cout << "EH" << hall_num << " AD" << ad_num << ":\t" << endl << "Distance" << endl << h_IBD_distance->Integral(h_IBD_distance->FindBin(2), h_IBD_distance->FindBin(5)) << " +/- " << sqrt(h_IBD_distance->Integral(h_IBD_distance->FindBin(2), h_IBD_distance->FindBin(5))) << "\t" <<h_rateSub_distance->Integral(h_rateSub_distance->FindBin(2), h_rateSub_distance->FindBin(5)) << "\t" << 100*(h_rateSub_distance->Integral(h_rateSub_distance->FindBin(2), h_rateSub_distance->FindBin(5)))/(h_IBD_distance->Integral(h_IBD_distance->FindBin(2), h_IBD_distance->FindBin(5))) << " %" << endl;

cout << "Distance (3sig)" << endl << h_IBD_distance_3sig->Integral(h_IBD_distance_3sig->FindBin(2), h_IBD_distance_3sig->FindBin(5)) << " +/- " << sqrt(h_IBD_distance_3sig->Integral(h_IBD_distance_3sig->FindBin(2), h_IBD_distance_3sig->FindBin(5))) << "\t" <<h_rateSub_distance_3sig->Integral(h_rateSub_distance_3sig->FindBin(2), h_rateSub_distance_3sig->FindBin(5)) << "\t" << 100*(h_rateSub_distance_3sig->Integral(h_rateSub_distance_3sig->FindBin(2), h_rateSub_distance_3sig->FindBin(5)))/(h_IBD_distance_3sig->Integral(h_IBD_distance_3sig->FindBin(2), h_IBD_distance_3sig->FindBin(5))) << " %" << endl;

cout << "DT (3sig)" << endl<< h_ibd_DT_3sig->Integral(h_ibd_DT_3sig->FindBin(3), h_ibd_DT_3sig->FindBin(10)) << " +/- " << sqrt(h_ibd_DT_3sig->Integral(h_ibd_DT_3sig->FindBin(3), h_ibd_DT_3sig->FindBin(10))) << "\t" <<h_sub_DT_3sig_rate->Integral(h_sub_DT_3sig_rate->FindBin(3), h_sub_DT_3sig_rate->FindBin(10)) << "\t" << 100*(h_sub_DT_3sig_rate->Integral(h_sub_DT_3sig_rate->FindBin(3), h_sub_DT_3sig_rate->FindBin(10)))/(h_ibd_DT_3sig->Integral(h_ibd_DT_3sig->FindBin(3), h_ibd_DT_3sig->FindBin(10))) << " %" << endl;


}

void all(int pd_window_microsec){
	int EH[8] = {1,1,2,2,3,3,3,3};
	int AD[8] = {1,2,1,2,1,2,3,4};

	for(int i = 0; i<8; i++){
		subtract(EH[i],AD[i],pd_window_microsec);
	}

}
