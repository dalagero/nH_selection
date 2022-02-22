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

void fastN(int EH, int AD){

	TFile *in_file = new TFile("./FastN.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("OWP_EH%d_nH_Ep1.5_DT0.8m_Ep_All",EH));

        char name[64];
	const int numBins = 34;
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	sprintf(name, "h_rebinned_spectrum");
	TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);

//	double temp_binCenter = 0;
//	double temp_binContent = 0;

	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
//		temp_binCenter = in_hist->GetBinCenter(iBin);
		if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
		h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
	} 

	TCanvas *test = new TCanvas("test","test");
	test->Divide(2,1);
	test->cd(1);
	in_hist->Draw();
	test->cd(2);
	h_rebinned_spectrum->Draw();


	for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO fast_neutron_spectrum VALUES (\"nH backgrounds, Olivia 1/21/2022\",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
	}



}


void accidentals(int EH, int AD){

	TFile *in_file = new TFile(Form("./TotaledSingles_1500_EH%d.root",EH));
	TH1F *in_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_DT800_3sig_scaled_ad%d",AD));

        char name[64];
	const int numBins = 34;
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	sprintf(name, "h_rebinned_spectrum");
	TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);

//	double temp_binCenter = 0;
//	double temp_binContent = 0;

	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
//		temp_binCenter = in_hist->GetBinCenter(iBin);
		if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
		h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
	} 

	TCanvas *test = new TCanvas("test","test");
	test->Divide(2,1);
	test->cd(1);
	in_hist->Draw();
	test->cd(2);
	h_rebinned_spectrum->Draw();


/*	for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO fast_neutron_spectrum VALUES (\"nH backgrounds, Olivia 1/21/2022\",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
	}*/

	double counts = 0;
	counts = in_hist->Integral();
	double scale = 0;
	TH1F *scaled_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_scaled_ad%d",AD));
	TH1F *before_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_before_ad%d",AD));
	scale = (scaled_hist->Integral())/(before_hist->Integral());

	cout << "INSERT INTO bg_counts VALUES(\"Olivia 1/21/2022\", " << EH << ", " << AD << ", \"accidental\", " << counts << ", " << sqrt(counts/scale) * scale << ");" << endl;

}


void errorOutput(double t13_best, double t13_low, double t13_high){

	cout << "Best fit for sin^2(2t13): " << pow(sin(2*t13_best),2) << endl;
	cout << "Uncertainty: +" << pow(sin(2*t13_high),2)-pow(sin(2*t13_best),2) << " -" << pow(sin(2*t13_best),2)-pow(sin(2*t13_low),2) << endl;

}


