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
#include "TLatex.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include <stdio.h>
#include <stdlib.h>

void plot_fits_rateOnly(){

	double input=0;
	double fit_mc=0;
	double chi_mc=0;
	double fit_nomc=0;
	double chi_nomc=0;

	TH1F* h_fit_results_mc = new TH1F("h_fit_results_mc","",21, -0.0025, 0.1025);
	TH1F* h_chi2_mc = new TH1F("h_chi2_mc","",21, -0.0025, 0.1025);
	TH1F* h_fit_results_nomc = new TH1F("h_fit_results_nomc","",21, -0.0025, 0.1025);
	TH1F* h_chi2_nomc = new TH1F("h_chi2_nomc","",21, -0.0025, 0.1025);

	ifstream myfile ("toyScan.txt"); //File with fit results
	if (myfile.is_open()){
		while (myfile >> input >> fit_mc >> chi_mc >> fit_nomc >> chi_nomc){ //Go line-by-line through the file
			//Fill the plot
			if(input == 0) continue;
			h_fit_results_mc->Fill(input*1.e-3, (fit_mc-(input*1.e-3))/(input*1.e-3));
			if(input == 0) continue;
			h_chi2_mc->Fill(input*1.e-3, chi_mc);
			if(input == 0) continue;
			h_fit_results_nomc->Fill(input*1.e-3, (fit_nomc-(input*1.e-3))/(input*1.e-3));
			if(input == 0) continue;
			h_chi2_nomc->Fill(input*1.e-3, chi_nomc);

		}
	myfile.close();
	}

	for(int ibin=1; ibin<22; ibin++){
		h_fit_results_mc->SetBinError(ibin, 0);
		h_chi2_mc->SetBinError(ibin, 0);
		h_fit_results_nomc->SetBinError(ibin, 0);
		h_chi2_nomc->SetBinError(ibin, 0);
	}

	TLine* line0 = new TLine(0.0025,0,0.1025,0);
	line0->SetLineStyle(2);

	TCanvas* c = new TCanvas("c", "c", 1200, 800);
	c->Divide(1,2);
	c->cd(1);
		h_fit_results_mc->GetXaxis()->SetTitle("Toy MC Input Value");
		h_fit_results_mc->GetYaxis()->SetTitle("(Fit-Input)/Input");
		h_fit_results_mc->GetYaxis()->SetRangeUser(-0.002,0.02);
		h_fit_results_mc->GetXaxis()->SetRangeUser(0.0025,0.1025);
		h_fit_results_mc->SetMarkerStyle(20);
		h_fit_results_mc->SetStats(0);
		h_fit_results_mc->GetXaxis()->SetLabelSize(0.06);
		h_fit_results_mc->GetXaxis()->SetTitleSize(0.06);
		h_fit_results_mc->GetYaxis()->SetLabelSize(0.06);
		h_fit_results_mc->GetYaxis()->SetTitleSize(0.06);
		h_fit_results_mc->Draw("P0");
		h_fit_results_nomc->SetLineColor(kRed);
		h_fit_results_nomc->SetMarkerColor(kRed);
		h_fit_results_nomc->SetMarkerStyle(21);
		h_fit_results_nomc->Draw("same P0");
		line0->Draw();
	c->cd(2);
		h_chi2_mc->GetXaxis()->SetTitle("Toy MC Input Value");
		h_chi2_mc->GetYaxis()->SetTitle("Min chi^{2}");
		h_chi2_mc->GetYaxis()->SetRangeUser(0,2.1*1.e-6);
		h_chi2_mc->GetXaxis()->SetRangeUser(0.0025,0.1025);
		h_chi2_mc->SetMarkerStyle(20);
		h_chi2_mc->SetStats(0);
		h_chi2_mc->GetXaxis()->SetLabelSize(0.06);
		h_chi2_mc->GetXaxis()->SetTitleSize(0.06);
		h_chi2_mc->GetYaxis()->SetLabelSize(0.06);
		h_chi2_mc->GetYaxis()->SetTitleSize(0.06);
		h_chi2_mc->Draw("P0");
		h_chi2_nomc->SetLineColor(kRed);
		h_chi2_nomc->SetMarkerColor(kRed);
		h_chi2_nomc->SetMarkerStyle(21);
		h_chi2_nomc->Draw("same P0");

}
