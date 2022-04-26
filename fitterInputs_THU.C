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
#include <iomanip>

using namespace std;

int EH[8] = {1,1,2,2,3,3,3,3};
int AD[8] = {1,2,1,2,1,2,3,4};
char prompt_label[64]="THU Selection"; //Label for the prompt spectra
char bgCount_label[64]="THU Selection"; //Label for the background counts
char accSpectra_label[64]="accidentals, Olivia 3/15/2022 NU"; //Label for the accidentals background spectra
char bgSpectra_label[64]="nH backgrounds, Olivia 3/15/2022 NU"; //Label for the non-accidentals background spectra
char muon_label[64]="THU Selection"; //Label for the livetime/muon table

//THU values:
const double livetime[8] = {980.7314, 1103.5306, 1229.7542, 1096.6410, 1679.7422, 1679.6887, 1679.3460, 1498.7921}; //livetime in days
const double tarpRelError = 0.0037; //target proton relative error

const double prompt_counts[8] = {473336, 543502, 541324, 470483, 158991, 162657, 162030,143988};
const double acc_counts_inputs[8] = {52995.79, 58763.43, 60355.65, 51542.52, 76802.72, 81325.65, 80594.85, 71219.39};
const double acc_counts_error[8] = {15.90, 17.63, 18.11, 15.46, 23.04, 24.40, 24.18, 21.37};

//Rates:
const double muon_rate[8] = {199.5163, 199.3905, 150.0165, 149.9443, 15.0450, 15.0450, 15.0449, 14.9246}; //rate in Hz
const double muon_eff[8] = {0.6382, 0.6351, 0.7063, 0.7057, 0.9659, 0.9659, 0.9657, 0.9661};
const double fastn_rate[3] = {2.22,1.69,0.16}; //rate of fast neutrons in 1/day NU
const double fastn_rate_err[3] = {0.27,0.19,0.02}; //error of rate of fast neutrons in 1/day NU
const double li9_rate[3] = {1.61,2.18,0.22}; //rate of li9 in 1/day NU
const double li9_rate_err[3] = {0.80,0.94,0.08}; //error of rate of li9 in 1/day NU
const double radn_rate[1] = {0.151}; //rate of rad n in 1/day
const double radn_rate_err[1] = {0.040}; //error of rate of rad n in 1/day
const double amc_rate[8] = {0.04,0.04,0.04,0.03,0.02,0.01,0.01,0.01}; //rate of amc in 1/day
const double amc_rate_err[8] = {0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01}; //error of rate of amc in 1/day


//nH Paper values
/*const double livetime[8] = {565.436*0.7949, 565.436*0.7920, 568.019*0.8334, 378.407*0.8333, 562.414*0.9814, 562.414*0.9814, 562.414*0.9812, 372.685*0.9814}; //livetime in days
const double tarpRelError = 0.005175; //target proton relative error

const double prompt_counts[8] = {217613, 219721, 208606, 136718, 56880, 56106, 59230, 38037};
const double acc_counts_inputs[8] = {26240, 25721, 25422, 16365, 29920, 30065, 32179, 20427};
const double acc_counts_error[8] = {49,49,43,29,19,20,21,15};

//Rates:
const double muon_rate[8] = {200.32, 200.32, 150.08, 149.80, 15.748, 15.748, 15.748, 15.747}; //rate in Hz
const double muon_eff[8] = {0.7949, 0.7920, 0.8334, 0.8333, 0.9814, 0.9814, 0.9812, 0.9814};
const double fastn_rate[3] = {2.11,1.81,0.16}; //rate of fast neutrons in 1/day NU
const double fastn_rate_err[3] = {0.18,0.17,0.03}; //error of rate of fast neutrons in 1/day NU
const double li9_rate[3] = {2.36,1.73,0.19}; //rate of li9 in 1/day NU
const double li9_rate_err[3] = {1.02,0.75,0.09}; //error of rate of li9 in 1/day NU
const double radn_rate[1] = {0.0000000001}; //rate of rad n in 1/day
const double radn_rate_err[1] = {0.0000000001}; //error of rate of rad n in 1/day
const double amc_rate[8] = {0.07,0.07,0.07,007,0.03,0.03,0.03,0.02}; //rate of amc in 1/day
const double amc_rate_err[8] = {0.04,0.04,0.03,0.03,0.02,0.02,0.02,0.01}; //error of rate of amc in 1/day*/

//nGd Paper values
/*const double livetime[8] = {1536.621,1737.616,1741.235,1554.044,1739.611,1739.611,1739.611,1551.945}; //livetime in days
const double tarpRelError = 0.0003; //target proton relative error

const double prompt_counts[8] = {830036,964381,889171,784736,127107,127726,126666,113922};
const double acc_counts_inputs[8] = {1536.621*8.27,1737.616*8.12,1741.235*6.00,1554.044*5.86,1739.611*1.06,1739.611*1.00,1739.611*1.03,1551.945*0.86};
const double acc_counts_error[8] = {1536.621*0.08,1737.616*0.08,1741.235*0.06,1554.044*0.06,1739.611*0.01,1739.611*0.01,1739.611*0.01,1551.945*0.01};

//Rates:
const double muon_rate[8] = {200.32, 200.32, 150.08, 149.80, 15.748, 15.748, 15.748, 15.747}; //rate in Hz
const double muon_eff[8] = {0.8261,0.8221,0.8577,0.8568,0.9833,0.9833,0.9832,0.9833};
const double fastn_rate[3] = {0.79,0.57,0.05}; //rate of fast neutrons in 1/day NU
const double fastn_rate_err[3] = {0.10,0.07,0.01}; //error of rate of fast neutrons in 1/day NU
const double li9_rate[3] = {2.38,1.59,0.19}; //rate of li9 in 1/day NU
const double li9_rate_err[3] = {0.66,0.49,0.08}; //error of rate of li9 in 1/day NU
const double radn_rate[8] = {0.08,0.06,0.04,0.06,0.04,0.04,0.04,0.04}; //rate of rad n in 1/day
const double radn_rate_err[8] = {0.04,0.03,0.02,0.03,0.02,0.02,0.02,0.02}; //error of rate of rad n in 1/day
const double amc_rate[8] = {0.17,0.15,0.14,0.13,0.06,0.05,0.05,0.04}; //rate of amc in 1/day
const double amc_rate_err[8] = {0.07,0.07,0.06,0.06,0.03,0.02,0.02,0.02}; //error of rate of amc in 1/day
*/

void muon(){

	for(int iad=0; iad<8; iad++){
		cout << "INSERT INTO muon_rates VALUES (" << EH[iad] << "," << AD[iad] << ",\"" << muon_label << "\"," << muon_rate[iad]*livetime[iad]*24*60*60 << "," << livetime[iad]*24*60*60*1.e9 << "," << muon_rate[iad] << "," << muon_eff[iad] << ");" << endl;
	}
}

void prompt(){

	for(int iad=0; iad<8; iad++){

		TFile *in_file = new TFile(Form("../nH_files/spectra_for_Olivia_EH%d.root",EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get("h_ibd_cand");


		char name[64];
		const int numBins = 34;
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		sprintf(name, "h_ibd_eprompt_inclusive_eh%d_ad%d", EH[iad], AD[iad]);
		TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);
		
		sprintf(name, "h_ibd_eprompt_fine_inclusive_eh%d_ad%d", EH[iad], AD[iad]);
		TH1F* h_fine_spectrum=new TH1F(name,name,120,0,12.);


		for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
			h_fine_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		}

		double scale = prompt_counts[iad]/h_rebinned_spectrum->Integral();
		h_rebinned_spectrum->Scale(scale);
		h_fine_spectrum->Scale(scale);


		cout << "INSERT INTO num_coincidences VALUES (" << EH[iad] << "," << AD[iad] << ",1,\"[";
		for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
			cout << (h_rebinned_spectrum->GetBinContent(iBin));
			if(iBin != h_rebinned_spectrum->GetNbinsX()) cout << ",";
		}
		cout << "]\",\"" << prompt_label << "\");" << endl;
		
//		cout << "\tEH" << EH[iad] << " AD" << AD[iad] << "\t" << h_rebinned_spectrum->Integral() << "\t+-\t" << sqrt(h_rebinned_spectrum->Integral())<< endl;

/*		TCanvas *test = new TCanvas("test","test");
		test->Divide(2,1);
		test->cd(1);
		in_hist->Draw();
		test->cd(2);
		h_rebinned_spectrum->Draw();*/
		
		in_file->Close();
	}

}

void acc_spectra(){

	for(int iad=0; iad<8; iad++){
		TFile *in_file = new TFile(Form("../nH_files/spectra_for_Olivia_EH%d.root",EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get("h_acc");

		char name[64];
		const int numBins = 34;
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		sprintf(name, "h_accidental_eprompt_inclusive_eh%d_ad%d", EH[iad], AD[iad]);
		TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);

		sprintf(name, "h_accidental_eprompt_fine_inclusive_eh%d_ad%d", EH[iad], AD[iad]);
		TH1F* h_fine_spectrum=new TH1F(name,name,120,0,12.);

		for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
			h_fine_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		} 


		for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
			cout << "INSERT INTO accidentals_spectrum VALUES (\"" << accSpectra_label << "\"," << EH[iad] << "," << AD[iad] << ",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
		}


		in_file->Close();
	
	}
}


void acc_counts(){

	for(int iad=0; iad<8; iad++){
		cout << "INSERT INTO bg_counts VALUES(\"" << bgCount_label << "\", " << EH[iad] << ", " << AD[iad] << ", \"accidental\", " << acc_counts_inputs[iad] << ", " << acc_counts_error[iad] << ");" << endl;
	}

}

void fastN_spectrum(int eh){

	TFile *in_file = new TFile(Form("../nH_files/spectra_for_Olivia_EH%d.root",eh));
	TH1F *in_hist = (TH1F*)in_file->Get("h_fn");

	char name[64];
	const int numBins = 34;
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	sprintf(name, "h_rebinned_spectrum");
	TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
		if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
		h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
	} 


	for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO fast_neutron_spectrum VALUES (\"" << bgSpectra_label << "\",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
	}
	
	in_file->Close();
}


void fastN_counts(){
	for(int iad=0; iad<8; iad++){
		double temp_livetime_day = livetime[2*(EH[iad]-1)+AD[iad]-1]; //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"fast-neutron\"," << fastn_rate[EH[iad]-1]*temp_livetime_day << "," << fastn_rate_err[EH[iad]-1]*temp_livetime_day << ");" << endl;
	}
}

void li9_spectrum(int eh){

	TFile *in_file = new TFile(Form("../nH_files/spectra_for_Olivia_EH%d.root",eh));
	TH1F *in_hist = (TH1F*)in_file->Get(Form("h_li9"));

        char name[64];
	const int numBins = 34;
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	sprintf(name, "h_rebinned_spectrum");
	TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
		if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
		h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
	} 


	for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO li9_spectrum VALUES (\"" << bgSpectra_label << "\",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
	}

	in_file->Close();
}

void li9_counts(){
	for(int iad = 0; iad < 8; iad++){
		double temp_livetime_day = livetime[iad]; //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"li9\"," << li9_rate[EH[iad]-1]*temp_livetime_day << "," << li9_rate_err[EH[iad]-1]*temp_livetime_day << ");" << endl;
	}
}

void amc_spectrum(int eh){

	TFile *in_file = new TFile(Form("../nH_files/spectra_for_Olivia_EH%d.root",eh));
	TH1F *in_hist = (TH1F*)in_file->Get(Form("h_amc"));

        char name[64];
	const int numBins = 34;
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	sprintf(name, "h_rebinned_spectrum");
	TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
		if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
		h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
	} 


	for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO amc_spectrum VALUES (\"" << bgSpectra_label << "\",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
	}

	in_file->Close();
}

void amc_counts(){
	for(int iad = 0; iad < 8; iad++){
		double temp_livetime_day = livetime[iad]; //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"amc\"," << amc_rate[iad]*temp_livetime_day << "," << amc_rate_err[iad]*temp_livetime_day << ");" << endl;
	}
}

void radN_spectrum(int eh){

	TFile *in_file = new TFile(Form("../nH_files/spectra_for_Olivia_EH%d.root",eh));
	TH1F *in_hist = (TH1F*)in_file->Get(Form("h_radn"));

        char name[64];
	const int numBins = 34;
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	sprintf(name, "h_rebinned_spectrum");
	TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
		if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
		h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
	} 


	for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO rad_n_spectrum VALUES (\"" << bgSpectra_label << "\",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
	}

	in_file->Close();
}

void radN_counts(){
	for(int iad = 0; iad < 8; iad++){
		double temp_livetime_day = livetime[iad]; //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"rad-n\"," << radn_rate[iad]*temp_livetime_day << "," << radn_rate_err[iad]*temp_livetime_day << ");" << endl;
	}
}


void allBg_counts(){
	acc_counts();
	fastN_counts();
	li9_counts();
	amc_counts();
	radN_counts();
}



void effError(double DTeffError_rel, double delayedEffError_rel){

	cout << "Detection Efficiency Error: " << sqrt(pow(DTeffError_rel,2) + pow(delayedEffError_rel,2) + pow(tarpRelError,2)) << endl;

}

void errorOutput(double t13_best, double t13_low, double t13_high){

	cout << "Best fit for sin^2(2t13): " << pow(sin(2*t13_best),2) << endl;
	cout << "Uncertainty: +" << pow(sin(2*t13_high),2)-pow(sin(2*t13_best),2) << " -" << pow(sin(2*t13_best),2)-pow(sin(2*t13_low),2) << endl;

}
