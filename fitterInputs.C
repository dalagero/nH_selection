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

int EH[8] = {1,1,2,2,3,3,3,3};
int AD[8] = {1,2,1,2,1,2,3,4};
char prompt_label[64]="Olivia 3/15/2022 NU"; //Label for the prompt spectra
char bgCount_label[64]="Olivia 3/15/2022 NU"; //Label for the background counts
char accSpectra_label[64]="accidentals, Olivia 3/15/2022 NU"; //Label for the accidentals background spectra
char bgSpectra_label[64]="nH backgrounds, Olivia 3/15/2022 NU"; //Label for the non-accidentals background spectra

const double livetime[8] = {72283163195734976, 81341302266285248, 94242900368444256, 84031181570224704, 143447962585023856, 143434576886712080, 143413538529066832, 127989198069118256}; //livetime in ns
const double tarpRelError = 0.0037; //target proton relative error

//Rates:
const double fastn_rate[3] = {2.16,1.61,0.14}; //rate of fast neutrons in 1/day NU
const double fastn_rate_err[3] = {0.14,0.08,0.02}; //error of rate of fast neutrons in 1/day NU
const double li9_rate[3] = {1.95,1.91,0.14}; //rate of li9 in 1/day NU
const double li9_rate_err[3] = {0.83,0.82,0.05}; //error of rate of li9 in 1/day NU
const double radn_rate[1] = {0.18}; //rate of rad n in 1/day
const double radn_rate_err[1] = {0.070}; //error of rate of rad n in 1/day
const double amc_rate[8] = {0.04,0.04,0.04,0.03,0.02,0.01,0.01,0.01}; //rate of amc in 1/day
const double amc_rate_err[8] = {0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01}; //error of rate of amc in 1/day


void prompt(){

	for(int iad=0; iad<8; iad++){

		TFile *in_file = new TFile(Form("../nH_files/TotaledPlots_NU_EH%d_1500.root",EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_DT800_3sig_ad%d",AD[iad]));


		char name[64];
		const int numBins = 34;
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		sprintf(name, "h_rebinned_spectrum");
		TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);


		for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
			if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
			h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		} 


		cout << "INSERT INTO num_coincidences VALUES (" << EH[iad] << "," << AD[iad] << ",1,\"[";
		for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
			cout << (h_rebinned_spectrum->GetBinContent(iBin));
			if(iBin != h_rebinned_spectrum->GetNbinsX()) cout << ",";
		}
		cout << "]\",\"" << prompt_label << "\");" << endl;
		
		in_file->Close();
	}


}

void acc_spectra(){

	for(int iad=0; iad<8; iad++){
		TFile *in_file = new TFile(Form("../nH_files/TotaledSingles_NU_1500_EH%d.root",EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_DT800_3sig_scaled_ad%d",AD[iad]));

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
			cout << "INSERT INTO accidentals_spectrum VALUES (\"" << accSpectra_label << "\"," << EH[iad] << "," << AD[iad] << ",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
		}

		in_file->Close();
	
	}

}


void acc_counts(){

	for(int iad=0; iad<8; iad++){
		TFile *in_file = new TFile(Form("../nH_files/TotaledSingles_NU_1500_EH%d.root",EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_DT800_3sig_scaled_ad%d",AD[iad]));

		double counts = 0;
		counts = in_hist->Integral();
		double scale = 0;
		TH1F *scaled_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_scaled_ad%d",AD[iad]));
		TH1F *before_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_before_ad%d",AD[iad]));
		scale = (scaled_hist->Integral())/(before_hist->Integral());

		cout << "INSERT INTO bg_counts VALUES(\"" << bgCount_label << "\", " << EH[iad] << ", " << AD[iad] << ", \"accidental\", " << counts << ", " << sqrt(counts/scale) * scale << ");" << endl;
		
		in_file->Close();
	
	}

}

void fastN_spectrum(int eh){

	TFile *in_file = new TFile("./FastNshape.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("OWP_EH%d_nH_Ep1.5_DT0.8m_Ep_EpLt12MeV",eh));

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
		double temp_livetime_day = livetime[2*(EH[iad]-1)+AD[iad]-1] * 1.e-9/(60*60*24); //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"fast-neutron\"," << fastn_rate[EH[iad]-1]*temp_livetime_day << "," << fastn_rate_err[EH[iad]-1]*temp_livetime_day << ");" << endl;
	}
}

void li9_spectrum(){

	TFile *in_file = new TFile("./Li9He8shape.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("Li9 prompt energy"));

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
		double temp_livetime_day = livetime[iad] * 1.e-9/(60*60*24); //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"li9\"," << li9_rate[EH[iad]-1]*temp_livetime_day << "," << li9_rate_err[EH[iad]-1]*temp_livetime_day << ");" << endl;
	}
}

void amc_counts(){
	for(int iad = 0; iad < 8; iad++){
		double temp_livetime_day = livetime[iad] * 1.e-9/(60*60*24); //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"amc\"," << amc_rate[iad]*temp_livetime_day << "," << amc_rate_err[iad]*temp_livetime_day << ");" << endl;
	}
}

void radN_counts(){
	for(int iad = 0; iad < 8; iad++){
		double temp_livetime_day = livetime[iad] * 1.e-9/(60*60*24); //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"rad-n\"," << radn_rate[0]*temp_livetime_day << "," << radn_rate_err[0]*temp_livetime_day << ");" << endl;
	}
}


void allBg_counts(){
	acc_counts();
	fastN_counts();
	li9_counts();
	amc_counts();
	radN_counts();
}



void effError(double DTeff, double DTeffError, double delayedEff, double delayedEffError){

	cout << "Detection Efficiency Error: " << sqrt(pow(DTeffError/DTeff,2) + pow(delayedEffError/delayedEff,2) + pow(tarpRelError,2)) << endl;

}

void errorOutput(double t13_best, double t13_low, double t13_high){

	cout << "Best fit for sin^2(2t13): " << pow(sin(2*t13_best),2) << endl;
	cout << "Uncertainty: +" << pow(sin(2*t13_high),2)-pow(sin(2*t13_best),2) << " -" << pow(sin(2*t13_best),2)-pow(sin(2*t13_low),2) << endl;

}
