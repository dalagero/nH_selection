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
			hist_ratio_prompt[iad]->GetYaxis()->SetRangeUser(0,2);
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
			hist_ratio_delayed[iad]->GetYaxis()->SetRangeUser(0,2);
			hist_ratio_delayed[iad]->SetLineColor(kBlack);
			hist_ratio_delayed[iad]->Draw();
			sprintf(title,"../nH_files/NU_ratio_delayed_EH%dAD%d.png",EH[iad],AD[iad]);
			ratio_delayed[iad]->Print(title);
		
		if(iad==0) cout << "Before:\t\tAfter:\t\tPercent Difference [(after-before)/before*100]" << endl;
		cout << hist_delayed_preNU[iad]->Integral() << "\t\t" << hist_delayed_postNU[iad]->Integral() << "\t\t" << 100*(hist_delayed_postNU[iad]->Integral()-hist_delayed_preNU[iad]->Integral())/hist_delayed_preNU[iad]->Integral() << endl;
	}



}
