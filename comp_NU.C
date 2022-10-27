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
#include "TBox.h"
#include "TLine.h"
#include <typeinfo>

using namespace std;

const int maxAD = 8;
const int EH[8] = {1,1,2,2,3,3,3,3};
const int AD[8] = {1,2,1,2,1,2,3,4};

void comp(){

	TH1F *hist_prompt_preNU[maxAD];
	TH1F *hist_prompt_postNU[maxAD];
	TCanvas *comp_prompt[maxAD];
	
	TH1F *hist_delayed_preNU[maxAD];
	TH1F *hist_delayed_postNU[maxAD];
	TCanvas *comp_delayed[maxAD];

	TH1F *hist_ratio_prompt[maxAD];
	TH1F *hist_ratio_delayed[maxAD];
	TCanvas *ratio_prompt[maxAD];
	TCanvas *ratio_delayed[maxAD];
	
        char title[80];

	for(int iad = 0; iad < maxAD; iad++){
		//Before NU correction file
		TFile *file_preNU = new TFile(Form("../nH_files/SubtractedAccidentals_1500_EH%dAD%d.root",EH[iad],AD[iad]));
		//After NU correction file
		TFile *file_postNU = new TFile(Form("../nH_files/SubtractedAccidentals_NU_1500_EH%dAD%d.root",EH[iad],AD[iad]));

		//prompt histograms
		hist_prompt_preNU[iad] = (TH1F*)file_preNU->Get(Form("h_Eprompt_subtract_DT800_ad%d",AD[iad]));
		hist_prompt_postNU[iad] = (TH1F*)file_postNU->Get(Form("h_Eprompt_subtract_DT800_ad%d",AD[iad]));
			comp_prompt[iad] = new TCanvas(Form("comp_prompt_ad%d",iad+1),Form("comp_prompt_ad%d",iad+1));
			comp_prompt[iad]->cd();
			hist_prompt_preNU[iad]->SetLineColor(kBlack);
			hist_prompt_preNU[iad]->Draw();
			hist_prompt_postNU[iad]->SetLineColor(kRed);
			hist_prompt_postNU[iad]->Draw("same");
			sprintf(title,"../nH_files/comp_prompt_EH%dAD%d.png",EH[iad],AD[iad]);
			comp_prompt[iad]->Print(title);
		
		hist_ratio_prompt[iad] = (TH1F*) hist_prompt_postNU[iad]->Clone();
		hist_ratio_prompt[iad]->Divide((TH1F*) hist_prompt_preNU[iad]);
			ratio_prompt[iad] = new TCanvas(Form("ratio_prompt_ad%d",iad+1),Form("ratio_prompt_ad%d",iad+1));
			ratio_prompt[iad]->cd();
			hist_ratio_prompt[iad]->GetXaxis()->SetTitle("Energy [MeV]");
			hist_ratio_prompt[iad]->GetYaxis()->SetTitle("Ratio (After NU/Before)");
			hist_ratio_prompt[iad]->GetYaxis()->SetRangeUser(0.5,1.5);
			hist_ratio_prompt[iad]->SetLineColor(kBlack);
			hist_ratio_prompt[iad]->Draw();
			sprintf(title,"../nH_files/NU_ratio_prompt_EH%dAD%d.png",EH[iad],AD[iad]);
			ratio_prompt[iad]->Print(title);
		
		//delayed histograms
		hist_delayed_preNU[iad] = (TH1F*)file_preNU->Get(Form("h_Edelayed_subtract_DT800_ad%d",AD[iad]));
		hist_delayed_postNU[iad] = (TH1F*)file_postNU->Get(Form("h_Edelayed_subtract_DT800_ad%d",AD[iad]));	
			comp_delayed[iad] = new TCanvas(Form("comp_delayed_ad%d",iad+1),Form("comp_delayed_ad%d",iad+1));
			comp_delayed[iad]->cd();
			hist_delayed_preNU[iad]->SetLineColor(kBlack);
			hist_delayed_preNU[iad]->Draw();
			hist_delayed_postNU[iad]->SetLineColor(kRed);
			hist_delayed_postNU[iad]->Draw("same");
			sprintf(title,"../nH_files/comp_delayed_EH%dAD%d.png",EH[iad],AD[iad]);
			comp_delayed[iad]->Print(title);
		hist_ratio_delayed[iad] = (TH1F*) hist_delayed_postNU[iad]->Clone();
		hist_ratio_delayed[iad]->Divide((TH1F*)hist_delayed_preNU[iad]);
			ratio_delayed[iad] = new TCanvas(Form("ratio_delayed_ad%d",iad+1),Form("ratio_delayed_ad%d",iad+1));
			ratio_delayed[iad]->cd();
			hist_ratio_delayed[iad]->GetXaxis()->SetTitle("Energy [MeV]");
			hist_ratio_delayed[iad]->GetYaxis()->SetTitle("Ratio (After NU/Before)");
			hist_ratio_delayed[iad]->GetYaxis()->SetRangeUser(0.5,1.5);
			hist_ratio_delayed[iad]->SetLineColor(kBlack);
			hist_ratio_delayed[iad]->Draw();
			sprintf(title,"../nH_files/NU_ratio_delayed_EH%dAD%d.png",EH[iad],AD[iad]);
			ratio_delayed[iad]->Print(title);
		
		if(iad==0) cout << "Before:\t\tAfter:\t\tPercent Difference [(after-before)/before*100]" << endl;
		cout << hist_delayed_preNU[iad]->Integral() << "\t\t" << hist_delayed_postNU[iad]->Integral() << "\t\t" << 100*(hist_delayed_postNU[iad]->Integral()-hist_delayed_preNU[iad]->Integral())/hist_delayed_preNU[iad]->Integral() << endl;
	}



}


void comp_thu_one(int eh, int ad){

	int run_num=0;
	if(eh==1) run_num=21221;
	if(eh==2) run_num=51052;
	if(eh==3) run_num=65844;
        char title[80];

	//Files
	TFile *file_none = new TFile(Form("../temp/summary_1500_%d.root",run_num));
	TFile *file_NU = new TFile(Form("../temp/summary_NU_1500_%d.root",run_num));
	TFile *file_thu = new TFile(Form("../temp/summary_thuNU_1500_%d.root",run_num));

	//Get the histogram
	TH1F *hist_none = (TH1F*)file_none->Get(Form("h_delayed_energy_DT800_ad%d",ad));
	TH1F *hist_NU = (TH1F*)file_NU->Get(Form("h_delayed_energy_DT800_ad%d",ad));
	TH1F *hist_thu = (TH1F*)file_thu->Get(Form("h_delayed_energy_DT800_ad%d",ad));
	
	//Rename
	hist_none->SetTitle("No Correction");
	hist_NU->SetTitle("NU Correction");
	hist_thu->SetTitle("THU Correction");
	
	//Rebin
	hist_none->Rebin(5);
	hist_NU->Rebin(5);
	hist_thu->Rebin(5);

	//Draw the histogram
	TCanvas *c1 = new TCanvas("c1","");
	c1->cd();
	hist_none->SetLineWidth(2);
	hist_none->SetLineColor(kBlack);
	hist_none->GetXaxis()->SetRangeUser(1.5,3);
	hist_none->Draw();
	hist_NU->SetLineWidth(2);
	hist_NU->SetLineColor(kBlue);
	hist_NU->Draw("same");
	hist_thu->SetLineWidth(2);
	hist_thu->SetLineColor(kRed);
	hist_thu->Draw("same");
	c1->BuildLegend();


}

void comp_thu_all(int eh, int ad){


	//Files
	TFile *file_none = new TFile(Form("../temp/TotaledPlots_EH%d_1500_1ad.root",eh));
	TFile *file_NU = new TFile(Form("../temp/TotaledPlots_NU_EH%d_1500_1ad.root",eh));
	TFile *file_thu = new TFile(Form("../temp/TotaledPlots_thuNU_EH%d_1500_1ad.root",eh));

	//Get the histogram
	TH1F *hist_none = (TH1F*)file_none->Get(Form("h_total_delayed_energy_DT800_ad%d",ad));
	TH1F *hist_NU = (TH1F*)file_NU->Get(Form("h_total_delayed_energy_DT800_ad%d",ad));
	TH1F *hist_thu = (TH1F*)file_thu->Get(Form("h_total_delayed_energy_DT800_ad%d",ad));
	
	//Rename
	hist_none->SetTitle("No Correction");
	hist_NU->SetTitle("NU Correction");
	hist_thu->SetTitle("THU Correction");
	
	//Rebin
	hist_none->Rebin(2);
	hist_NU->Rebin(2);
	hist_thu->Rebin(2);

	//Ratio plots
	TH1F *ratio_NU = (TH1F*)hist_NU->Clone("ratio_NU");
	TH1F *ratio_thu = (TH1F*)hist_thu->Clone("ratio_thu");
	ratio_NU->Divide(hist_none);
	ratio_thu->Divide(hist_none);
		ratio_NU->SetTitle("NU/raw ratio");
		ratio_thu->SetTitle("THU/raw ratio");

	//Draw the histogram
	TCanvas *c1 = new TCanvas("c1","");
	c1->cd();
	hist_none->SetLineWidth(2);
	hist_none->SetLineColor(kBlack);
	hist_none->GetXaxis()->SetRangeUser(1.5,3);
	hist_none->Draw();
	hist_NU->SetLineWidth(2);
	hist_NU->SetLineColor(kBlue);
	hist_NU->Draw("same");
	hist_thu->SetLineWidth(2);
	hist_thu->SetLineColor(kRed);
	hist_thu->Draw("same");
	c1->BuildLegend();
	
	TCanvas *c2 = new TCanvas("c2","");
	c2->cd();
	ratio_NU->SetLineWidth(2);
	ratio_NU->SetLineColor(kBlue);
	ratio_NU->GetXaxis()->SetRangeUser(1.5,3);
	ratio_NU->GetYaxis()->SetTitle("Ratio to No Correction");
	ratio_NU->Draw();
	ratio_thu->SetLineWidth(2);
	ratio_thu->SetLineColor(kRed);
	ratio_thu->Draw("same");
	c2->BuildLegend();


}

void comp_thu_acc(int eh, int ad){


	//Files
	TFile *file_none = new TFile(Form("../temp/TotaledSingles_1500_EH%d.root",eh));
	TFile *file_NU = new TFile(Form("../temp/TotaledSingles_NU_1500_EH%d.root",eh));
	TFile *file_thu = new TFile(Form("../temp/TotaledSingles_thuNU_1500_EH%d.root",eh));

	//Get the histogram
	TH1F *hist_none = (TH1F*)file_none->Get(Form("h_total_delayed_energy_DT800_ad%d",ad));
	TH1F *hist_NU = (TH1F*)file_NU->Get(Form("h_total_delayed_energy_DT800_ad%d",ad));
	TH1F *hist_thu = (TH1F*)file_thu->Get(Form("h_total_delayed_energy_DT800_ad%d",ad));
	
	//Rename
	hist_none->SetTitle("No Correction");
	hist_NU->SetTitle("NU Correction");
	hist_thu->SetTitle("THU Correction");
	
	//Rebin
	hist_none->Rebin(2);
	hist_NU->Rebin(2);
	hist_thu->Rebin(2);

	//Ratio plots
	TH1F *ratio_NU = (TH1F*)hist_NU->Clone("ratio_NU");
	TH1F *ratio_thu = (TH1F*)hist_thu->Clone("ratio_thu");
	ratio_NU->Divide(hist_none);
	ratio_thu->Divide(hist_none);
		ratio_NU->SetTitle("NU/raw ratio");
		ratio_thu->SetTitle("THU/raw ratio");

	//Draw the histogram
	TCanvas *c1 = new TCanvas("c1","");
	c1->cd();
	hist_none->SetStats(0);
	hist_none->SetLineWidth(2);
	hist_none->SetLineColor(kBlack);
	hist_none->GetXaxis()->SetRangeUser(1.5,3);
	hist_none->Draw();
	hist_NU->SetLineWidth(2);
	hist_NU->SetLineColor(kBlue);
	hist_NU->Draw("same");
	hist_thu->SetLineWidth(2);
	hist_thu->SetLineColor(kRed);
	hist_thu->Draw("same");
	c1->BuildLegend();
	
	TCanvas *c2 = new TCanvas("c2","");
	c2->cd();
	ratio_NU->SetStats(0);
	ratio_NU->SetLineWidth(2);
	ratio_NU->SetLineColor(kBlue);
	ratio_NU->GetXaxis()->SetRangeUser(1.5,3);
	ratio_NU->GetYaxis()->SetTitle("Ratio to No Correction");
	ratio_NU->Draw();
	ratio_thu->SetLineWidth(2);
	ratio_thu->SetLineColor(kRed);
	ratio_thu->Draw("same");
	c2->BuildLegend();


}


void comp_thu_delayedFits(){

	//Files
	TFile *file_NU = new TFile("../nH_files/delayedEnergy.root");
	TFile *file_thu = new TFile("../nH_files/delayedEnergy_thuNU.root");
	TFile *file_yasu = new TFile("../nH_files/delayedEnergy_yasu.root");

	//Get the histogram
	TH1F *peak_NU = (TH1F*)file_NU->Get("h_peak_rate_ADs");
	TH1F *peak_thu = (TH1F*)file_thu->Get("h_peak_rate_ADs");
		TH1F *peak_yasu = (TH1F*)file_yasu->Get("h_peak_rate_ADs");
	TH1F *sigma_NU = (TH1F*)file_NU->Get("h_sigma_rate_ADs");
	TH1F *sigma_thu = (TH1F*)file_thu->Get("h_sigma_rate_ADs");
	TH1F *sigma_yasu = (TH1F*)file_yasu->Get("h_sigma_rate_ADs");
	
	//Rename
	peak_NU->SetTitle("NU Correction");
	peak_thu->SetTitle("THU Correction");
	sigma_NU->SetTitle("NU Correction");
	sigma_thu->SetTitle("THU Correction");
	peak_yasu->SetTitle("w/ Yasu's factors");
	sigma_yasu->SetTitle("w/ Yasu's factors");

	//Draw the histogram
	TCanvas *c1 = new TCanvas("c1","");
	c1->cd();
	peak_NU->SetLineWidth(2);
	peak_NU->SetLineColor(kBlue);
	peak_NU->SetMarkerColor(kBlue);
	peak_NU->GetYaxis()->SetRangeUser(2.235,2.29);
	peak_NU->Draw();
	peak_thu->SetLineWidth(2);
	peak_thu->SetLineColor(kRed);
	peak_thu->Draw("same");
	peak_yasu->SetLineWidth(2);
	peak_yasu->SetMarkerSize(3);
	peak_yasu->SetMarkerStyle(33);
	peak_yasu->SetLineColor(kGreen+2);
	peak_yasu->SetMarkerColor(kGreen+2);
	peak_yasu->Draw("P same");
	c1->BuildLegend();
	
	//Draw the histogram
	TCanvas *c2 = new TCanvas("c2","");
	c2->cd();
	sigma_NU->SetLineWidth(2);
	sigma_NU->SetLineColor(kBlue);
	sigma_NU->SetMarkerColor(kBlue);
	sigma_NU->Draw();
	sigma_thu->SetLineWidth(2);
	sigma_thu->SetLineColor(kRed);
	sigma_thu->Draw("same");
	sigma_yasu->SetLineWidth(2);
	sigma_yasu->SetMarkerSize(3);
	sigma_yasu->SetMarkerStyle(33);
	sigma_yasu->SetLineColor(kGreen+2);
	sigma_yasu->SetMarkerColor(kGreen+2);
	sigma_yasu->Draw("P same");
	c2->BuildLegend();
}

void sub_prompt(){

	const int EH[8]={1,1,2,2,3,3,3,3};
	const int AD[8]={1,2,1,2,1,2,3,4};

	for(int iad=0; iad<8; iad++){
		//Files
		TFile *file_ibd = new TFile(Form("../nH_files/mcSwap/THUfiles/toy_s7_m24.root"));
		TFile *file_acc = new TFile(Form("../nH_files/mcSwap/THUfiles/toy_s7_m24.root"));

		//Get the histogram
		TH1F *hist_ibd_raw = (TH1F*)file_ibd->Get(Form("h_numCoincs_eh%dad%d",EH[iad],AD[iad]));
		TH1F *hist_acc_raw = (TH1F*)file_acc->Get(Form("h_acc_eh%dad%d",EH[iad],AD[iad]));
		
		//Rebin
		const int numBins=34;
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		TH1F *hist_ibd = new TH1F("hist_ibd", "hist_ibd", numBins, binEdges);
		TH1F *hist_acc = new TH1F("hist_acc", "hist_acc", numBins, binEdges);
		for(int iBin=1; iBin < (hist_ibd_raw->GetNbinsX())+1; iBin++){
			if((hist_ibd_raw->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((hist_ibd_raw->GetBinCenter(iBin))<binEdges[0]) continue;
			hist_ibd->Fill(hist_ibd_raw->GetBinCenter(iBin),hist_ibd_raw->GetBinContent(iBin));
		}
		for(int iBin=1; iBin < (hist_acc_raw->GetNbinsX())+1; iBin++){
			if((hist_acc_raw->GetBinCenter(iBin))>binEdges[numBins]) break;
			if((hist_acc_raw->GetBinCenter(iBin))<binEdges[0]) continue;
			hist_acc->Fill(hist_acc_raw->GetBinCenter(iBin),hist_acc_raw->GetBinContent(iBin));
		}
		TH1F *hist_sub = (TH1F*)hist_ibd->Clone();
		hist_sub->Add(hist_acc,-1);
		
		//Rename
		hist_ibd->SetTitle(Form("IBD Candidates EH%dAD%d",EH[iad],AD[iad]));
		hist_acc->SetTitle(Form("Accidentals EH%dAD%d",EH[iad],AD[iad]));
		hist_sub->SetTitle(Form("Subtracted EH%dAD%d",EH[iad],AD[iad]));
		hist_ibd->SetName(Form("IBD Candidates EH%dAD%d",EH[iad],AD[iad]));
		hist_acc->SetName(Form("Accidentals EH%dAD%d",EH[iad],AD[iad]));
		hist_sub->SetName(Form("Subtracted EH%dAD%d",EH[iad],AD[iad]));
		
		cout << "EH" << EH[iad] << "-AD" << AD[iad] << "\tCandidates: " << hist_ibd->Integral() << "\tAccidentals: " << hist_acc->Integral() << endl;

		//Draw the histogram
		TCanvas *c1 = new TCanvas("c1","");
		c1->cd();
		hist_ibd->SetStats(0);
		hist_ibd->SetLineWidth(2);
		hist_ibd->SetLineColor(kBlue);
	//	hist_ibd->SetMarkerColor(kBlue);
		hist_ibd->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		hist_ibd->GetYaxis()->SetTitle("Counts");
		hist_ibd->Draw("hist");
		hist_acc->SetLineWidth(2);
		hist_acc->SetLineColor(kRed);
		hist_acc->Draw("hist same");
		c1->BuildLegend();
		
		//Draw the histogram
		TCanvas *c2 = new TCanvas("c2","");
		c2->cd();
		hist_sub->SetStats(0);
		hist_sub->SetLineWidth(2);
		hist_sub->SetLineColor(kBlue);
	//	hist_sub->SetMarkerColor(kBlue);
		hist_sub->GetXaxis()->SetTitle("Prompt Energy [MeV]");
		hist_sub->GetYaxis()->SetTitle("Counts");
		hist_sub->Draw("hist");
		//c2->BuildLegend();
	}	
}

void plot_results(){

	const int EH[8]={1,1,2,2,3,3,3,3};
	const int AD[8]={1,2,1,2,1,2,3,4};
	const int numPoints=5;
	const double offset=0.2;
	
	//THU reported, raw, NU, thuNU, Yasu
	const char *label[numPoints] = {"Nominal NU", "UCI NU", "THU NU", "Yasu+THU", "Yasu+SDU"};
	const double s22t13_rate_UCI[numPoints] = {0.0764613,0.0754086,0.0777537,0.0774256,0.0747688};
	const double s22t13_rate_THU[numPoints] = {0.0766904,0.0756202,0.0779831,0.0776509,0.0749872};
	const double s22t13_shape_UCI[numPoints] = {0.0808689,0.0784123,0.0770893,0.0783671,0.0768928};
	const double s22t13_shape_THU[numPoints] = {0.0811874,0.0781192,0.0774397,0.0786693,0.0775780};
	const double dm2ee_shape_UCI[numPoints] = {0.00265988,0.00270006,0.00280843,0.00276302,0.00277238};
	const double dm2ee_shape_THU[numPoints] = {0.00262679,0.00267232,0.00276514,0.00272438,0.00273508};

//THU Technote results
		const double s22t13_rate_val = 0.0771;
		const double s22t13_shape_val = 0.0766;
		const double dm2ee_shape_val = 0.0027317;
		const double s22t13_rate_error = 0.0070;
		const double s22t13_shape_error = 0.0057;
		const double dm2ee_shape_error = 0.00014;

//2022 nGd results
	//	const double nGd_s22t13_rate_val = ;
		const double nGd_s22t13_shape_val = 0.0851;
		const double nGd_dm2ee_shape_val = 0.0025177;
	//	const double nGd_s22t13_rate_error = ;
		const double nGd_s22t13_shape_error = 0.0024;
		const double nGd_dm2ee_shape_error = 0.00006;
		
		TBox *THU_rate_s22t13_range = new TBox(0.5,s22t13_rate_val-s22t13_rate_error,numPoints+offset+0.5,s22t13_rate_val+s22t13_rate_error);
		TBox *THU_shape_s22t13_range = new TBox(0.5,s22t13_shape_val-s22t13_shape_error,numPoints+offset+0.5,s22t13_shape_val+s22t13_shape_error);
		TBox *THU_shape_dm2ee_range = new TBox(0.5,dm2ee_shape_val-dm2ee_shape_error,numPoints+offset+0.5,dm2ee_shape_val+dm2ee_shape_error);
		TLine *THU_rate_s22t13_value = new TLine(0.5,s22t13_rate_val,numPoints+offset+0.5,s22t13_rate_val);
		TLine *THU_shape_s22t13_value = new TLine(0.5,s22t13_shape_val,numPoints+offset+0.5,s22t13_shape_val);
		TLine *THU_shape_dm2ee_value = new TLine(0.5,dm2ee_shape_val,numPoints+offset+0.5,dm2ee_shape_val);

		TLine *nGd_shape_s22t13_value = new TLine(0.5,nGd_s22t13_shape_val,numPoints+offset+0.5,nGd_s22t13_shape_val);
		TLine *nGd_shape_dm2ee_value = new TLine(0.5,nGd_dm2ee_shape_val,numPoints+offset+0.5,nGd_dm2ee_shape_val);
		TBox *nGd_shape_s22t13_range = new TBox(0.5,nGd_s22t13_shape_val-nGd_s22t13_shape_error,numPoints+offset+0.5,nGd_s22t13_shape_val+nGd_s22t13_shape_error);
		TBox *nGd_shape_dm2ee_range = new TBox(0.5,nGd_dm2ee_shape_val-nGd_dm2ee_shape_error,numPoints+offset+0.5,nGd_dm2ee_shape_val+nGd_dm2ee_shape_error);

	TH1F *h_frame_s22t13_rate = new TH1F("h_frame_s22t13_rate", "h_frame_s22t13_rate", numPoints, 0.5, numPoints+offset+0.5);
	TH1F *h_frame_s22t13_shape = new TH1F("h_frame_s22t13_shape", "h_frame_s22t13_shape", numPoints, 0.5, numPoints+offset+0.5);
	TH1F *h_frame_dm2ee = new TH1F("h_frame_dm2ee", "h_frame_dm2ee", numPoints, 0.5, numPoints+offset+0.5);

	TH1F *h_rate_UCI = new TH1F("h_rate_UCI", "h_rate_UCI", numPoints, 0.5, numPoints+0.5);
	TH1F *h_rate_THU = new TH1F("h_rate_THU", "h_rate_THU", numPoints, 0.5+offset, numPoints+offset+0.5);
	TH1F *h_shape_s22t13_UCI = new TH1F("h_shape_s22t13_UCI", "h_shape_s22t13_UCI", numPoints, 0.5, numPoints+0.5);
	TH1F *h_shape_s22t13_THU = new TH1F("h_shape_s22t13_THU", "h_shape_s22t13_THU", numPoints, 0.5+offset, numPoints+offset+0.5);
	TH1F *h_shape_dm2ee_UCI = new TH1F("h_shape_dm2ee_UCI", "h_shape_dm2ee_UCI", numPoints, 0.5, numPoints+0.5);
	TH1F *h_shape_dm2ee_THU = new TH1F("h_shape_dm2ee_THU", "h_shape_dm2ee_THU", numPoints, 0.5+offset, numPoints+offset+0.5);

	for(int ipoint=0; ipoint<numPoints; ipoint++){
		//Fill the histograms:
		h_rate_UCI->Fill(ipoint+1, s22t13_rate_UCI[ipoint]);
		h_rate_THU->Fill(ipoint+1, s22t13_rate_THU[ipoint]);
		h_shape_s22t13_UCI->Fill(ipoint+1, s22t13_shape_UCI[ipoint]);
		h_shape_s22t13_THU->Fill(ipoint+1, s22t13_shape_THU[ipoint]);
		h_shape_dm2ee_UCI->Fill(ipoint+1, dm2ee_shape_UCI[ipoint]);
		h_shape_dm2ee_THU->Fill(ipoint+1, dm2ee_shape_THU[ipoint]);
		
		//Change axis labels:
		h_frame_s22t13_rate->GetXaxis()->ChangeLabel(ipoint+1,-1,-1,-1,-1,-1,label[ipoint]);
		h_frame_s22t13_shape->GetXaxis()->ChangeLabel(ipoint+1,-1,-1,-1,-1,-1,label[ipoint]);
		h_frame_dm2ee->GetXaxis()->ChangeLabel(ipoint+1,-1,-1,-1,-1,-1,label[ipoint]);

	}


	//Draw the histogram
	TCanvas *c1 = new TCanvas("c1","");
	c1->cd();
	int numSigmas=2;
	h_frame_s22t13_rate->SetStats(0);
//	h_frame_s22t13_rate->GetYaxis()->SetRangeUser(s22t13_rate_THU[0]-numSigmas*s22t13_rate_error,s22t13_rate_THU[0]+numSigmas*s22t13_rate_error);
	h_frame_s22t13_rate->GetYaxis()->SetRangeUser(0.065,0.09);
	h_frame_s22t13_rate->GetYaxis()->SetTitle("s22t13 result");
	h_frame_s22t13_rate->Draw();
	THU_rate_s22t13_range->SetFillColor(kGray+2);
	THU_rate_s22t13_range->SetFillStyle(3002);
	THU_rate_s22t13_range->Draw("same");
	THU_rate_s22t13_value->SetLineStyle(2);
	THU_rate_s22t13_value->Draw("same");
	nGd_shape_s22t13_range->SetFillColor(kGreen-3);
	nGd_shape_s22t13_range->SetFillStyle(3001);
	nGd_shape_s22t13_range->Draw("same");
	nGd_shape_s22t13_value->SetLineColor(kGreen+3);
	nGd_shape_s22t13_value->Draw("same");
	h_rate_UCI->SetLineWidth(0);
	h_rate_UCI->SetLineColor(kBlue);
	h_rate_UCI->SetMarkerStyle(20);
	h_rate_UCI->SetMarkerSize(2);
	h_rate_UCI->SetMarkerColor(kBlue);
	h_rate_UCI->Draw("P same");
		h_rate_THU->SetLineWidth(0);
		h_rate_THU->SetLineColor(kRed);
		h_rate_THU->SetMarkerStyle(21);
		h_rate_THU->SetMarkerSize(2);
		h_rate_THU->SetMarkerColor(kRed);
		h_rate_THU->Draw("P same");

	TCanvas *c2 = new TCanvas("c2","");
	c2->cd();
	h_frame_s22t13_shape->SetStats(0);
//	h_frame_s22t13_shape->GetYaxis()->SetRangeUser(s22t13_shape_THU[0]-numSigmas*s22t13_shape_error,s22t13_shape_THU[0]+numSigmas*s22t13_shape_error);
	h_frame_s22t13_shape->GetYaxis()->SetRangeUser(0.065,0.09);
	h_frame_s22t13_shape->GetYaxis()->SetTitle("s22t13 result");
	h_frame_s22t13_shape->Draw();
	THU_shape_s22t13_range->SetFillColor(kGray+2);
	THU_shape_s22t13_range->SetFillStyle(3002);
	THU_shape_s22t13_range->Draw("same");
	THU_shape_s22t13_value->SetLineStyle(2);
	THU_shape_s22t13_value->Draw("same");
	nGd_shape_s22t13_range->SetFillColor(kGreen-3);
	nGd_shape_s22t13_range->SetFillStyle(3001);
	nGd_shape_s22t13_range->Draw("same");
	nGd_shape_s22t13_value->SetLineColor(kGreen+3);
	nGd_shape_s22t13_value->Draw("same");
	h_shape_s22t13_UCI->SetLineWidth(0);
	h_shape_s22t13_UCI->SetLineColor(kBlue);
	h_shape_s22t13_UCI->SetMarkerStyle(20);
	h_shape_s22t13_UCI->SetMarkerSize(2);
	h_shape_s22t13_UCI->SetMarkerColor(kBlue);
	h_shape_s22t13_UCI->Draw("P same");
		h_shape_s22t13_THU->SetLineWidth(0);
		h_shape_s22t13_THU->SetLineColor(kRed);
		h_shape_s22t13_THU->SetMarkerStyle(21);
		h_shape_s22t13_THU->SetMarkerSize(2);
		h_shape_s22t13_THU->SetMarkerColor(kRed);
		h_shape_s22t13_THU->Draw("P same");

	TCanvas *c3 = new TCanvas("c3","");
	c3->cd();
	h_frame_dm2ee->SetStats(0);
//	h_frame_dm2ee->GetYaxis()->SetRangeUser(dm2ee_shape_THU[0]-numSigmas*dm2ee_shape_error,dm2ee_shape_THU[0]+numSigmas*dm2ee_shape_error);
	h_frame_dm2ee->GetYaxis()->SetRangeUser(0.0024,0.0029);
	h_frame_dm2ee->GetYaxis()->SetTitle("dm2ee result");
	h_frame_dm2ee->Draw();
	THU_shape_dm2ee_range->SetFillColor(kGray+2);
	THU_shape_dm2ee_range->SetFillStyle(3002);
	THU_shape_dm2ee_range->Draw("same");
	THU_shape_dm2ee_value->SetLineStyle(2);
	THU_shape_dm2ee_value->Draw("same");
	nGd_shape_dm2ee_range->SetFillColor(kGreen-3);
	nGd_shape_dm2ee_range->SetFillStyle(3001);
	nGd_shape_dm2ee_range->Draw("same");
	nGd_shape_dm2ee_value->SetLineColor(kGreen+3);
	nGd_shape_dm2ee_value->Draw("same");
	h_shape_dm2ee_UCI->SetLineWidth(0);
	h_shape_dm2ee_UCI->SetLineColor(kBlue);
	h_shape_dm2ee_UCI->SetMarkerStyle(20);
	h_shape_dm2ee_UCI->SetMarkerSize(2);
	h_shape_dm2ee_UCI->SetMarkerColor(kBlue);
	h_shape_dm2ee_UCI->Draw("P same");
		h_shape_dm2ee_THU->SetLineWidth(0);
		h_shape_dm2ee_THU->SetLineColor(kRed);
		h_shape_dm2ee_THU->SetMarkerStyle(21);
		h_shape_dm2ee_THU->SetMarkerSize(2);
		h_shape_dm2ee_THU->SetMarkerColor(kRed);
		h_shape_dm2ee_THU->Draw("P same");

}


