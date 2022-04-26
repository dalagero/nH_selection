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
char prompt_label[64]="toyMC"; //Label for the prompt spectra
char bgCount_label[64]="toyMC"; //Label for the background counts
char accSpectra_label[64]="accidentals, toyMC"; //Label for the accidentals background spectra
char bgSpectra_label[64]="nH backgrounds, toyMC"; //Label for the non-accidentals background spectra
char muon_label[64]="toyMC"; //Label for the livetime/muon table

//toyMC values
const double livetime[8] = {836.1607,940.9862,1090.3267,972.2450,1660.1904,1660.0705,1659.8232,1481.3565}; //livetime in days
const double tarpRelError = 0.0037; //target proton relative error

const double prompt_counts[8] = {406692,467425,485880,423255,163315,163117,165539,144359};
const double acc_counts_inputs[8] = {55.14*836.1607,54.49*940.9862,50.95*1090.3267,49.27*972.2450,48.74*1660.1904,49.38*1660.0705,50.63*1659.8232,48.47*1481.3565};
const double acc_counts_error[8] = {0.04*836.1607,0.04*940.9862,0.03*1090.3267,0.04*972.2450,0.03*1660.1904,0.03*1660.0705,0.03*1659.8232,0.03*1481.3565};

//Rates:
const double muon_rate[8] = {200.32, 200.32, 150.08, 149.80, 15.748, 15.748, 15.748, 15.747}; //rate in Hz
const double muon_eff[8] = {0.5442,0.5415,0.6262,0.6256,0.9543,0.9543,0.9541,0.9545};
const double fastn_rate[3] = {2.15,1.59,0.14}; //rate of fast neutrons in 1/day NU
const double fastn_rate_err[3] = {0.13,0.08,0.02}; //error of rate of fast neutrons in 1/day NU
const double li9_rate[3] = {1.95,1.91,0.15}; //rate of li9 in 1/day NU
const double li9_rate_err[3] = {0.83,0.82,0.05}; //error of rate of li9 in 1/day NU
const double radn_rate[8] = {0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15}; //rate of rad n in 1/day
const double radn_rate_err[8] = {0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04}; //error of rate of rad n in 1/day
const double amc_rate[8] = {0.04,0.04,0.04,0.03,0.02,0.01,0.01,0.01}; //rate of amc in 1/day
const double amc_rate_err[8] = {0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01}; //error of rate of amc in 1/day




void muon(){

	for(int iad=0; iad<8; iad++){
		cout << "INSERT INTO muon_rates VALUES (" << EH[iad] << "," << AD[iad] << ",\"" << muon_label << "\"," << muon_rate[iad]*livetime[iad]*24*60*60 << "," << livetime[iad]*24*60*60*1.e9 << "," << muon_rate[iad] << "," << muon_eff[iad] << ");" << endl;
	}
}

void prompt(){

	TFile *in_file = new TFile("../nH_files/toyMC_spectra/toys.root");
	auto* tree= (TTree*) in_file->Get("tr");
	TH1F* spectrum=0;

	for(int iad=0; iad<8; iad++){
		tree->SetBranchAddress(Form("h_stage2_ad%d",iad+1),&spectrum);
		tree->GetEntry(0);
//		spectrum->Draw("hist");


		cout << "INSERT INTO num_coincidences VALUES (" << EH[iad] << "," << AD[iad] << ",1,\"[";
		for(int iBin=1; iBin < (spectrum->GetNbinsX())+1; iBin++){
			cout << (spectrum->GetBinContent(iBin));
			if(iBin != spectrum->GetNbinsX()) cout << ",";
		}
		cout << "]\",\"" << prompt_label << "\");" << endl;
		
//		cout << "\tEH" << EH[iad] << " AD" << AD[iad] << "\t" << h_rebinned_spectrum->Integral() << "\t+-\t" << sqrt(h_rebinned_spectrum->Integral())<< endl;

/*		TCanvas *test = new TCanvas("test","test");
		test->Divide(2,1);
		test->cd(1);
		in_hist->Draw();
		test->cd(2);
		h_rebinned_spectrum->Draw();*/
		
	}

		in_file->Close();

}

void acc_spectra(){

	TFile *in_file = new TFile("../nH_files/toyMC_spectra/accidental_eprompt_shapes_8ad.root");

	for(int iad=0; iad<8; iad++){
		TH1F *in_hist = (TH1F*)in_file->Get(Form("h_accidental_eprompt_inclusive_eh%d_ad%d",EH[iad],AD[iad]));



		for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
			cout << "INSERT INTO accidentals_spectrum VALUES (\"" << accSpectra_label << "\"," << EH[iad] << "," << AD[iad] << ",1," << iBin-1 << "," << (in_hist->GetBinContent(iBin))/(in_hist->Integral()) << ");" << endl;
		}

	
	}
	in_file->Close();
}


void acc_counts(){

	for(int iad=0; iad<8; iad++){
		cout << "INSERT INTO bg_counts VALUES(\"" << bgCount_label << "\", " << EH[iad] << ", " << AD[iad] << ", \"accidental\", " << acc_counts_inputs[iad] << ", " << acc_counts_error[iad] << ");" << endl;
	}

}

void fastN_spectrum(){

	TFile *in_file = new TFile("../nH_files/toyMC_spectra/P15A_fn_spectrum.root");
	TH1F *in_hist = (TH1F*)in_file->Get("h_1AD_fn");


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO fast_neutron_spectrum VALUES (\"" << bgSpectra_label << "\",1," << iBin-1 << "," << (in_hist->GetBinContent(iBin))/(in_hist->Integral()) << ");" << endl;
	}
	
	in_file->Close();
}


void fastN_counts(){
	for(int iad=0; iad<8; iad++){
		double temp_livetime_day = livetime[2*(EH[iad]-1)+AD[iad]-1]; //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"fast-neutron\"," << fastn_rate[EH[iad]-1]*temp_livetime_day << "," << fastn_rate_err[EH[iad]-1]*temp_livetime_day << ");" << endl;
	}
}

void li9_spectrum(){

	TFile *in_file = new TFile("../nH_files/toyMC_spectra/8he9li_nominal_spectrum.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("h_nominal"));

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

void amc_spectrum(){

	TFile *in_file = new TFile("../nH_files/toyMC_spectra/amc_spectrum.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("h_toy"));

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

void radN_spectrum(){

	TFile *in_file = new TFile("../nH_files/toyMC_spectra/result-DocDB9667.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("AD1"));


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO rad_n_spectrum VALUES (\"" << bgSpectra_label << "\",1," << iBin-1 << "," << (in_hist->GetBinContent(iBin))/(in_hist->Integral()) << ");" << endl;
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
