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
char prompt_label[64]="Olivia 4/22/2022 NU DT2000"; //Label for the prompt spectra
char bgCount_label[64]="Olivia 4/22/2022 NU DT300"; //Label for the background counts
char accSpectra_label[64]="accidentals, Olivia 4/22/2022 NU DT2000"; //Label for the accidentals background spectra
char bgSpectra_label[64]="nH backgrounds, Olivia 3/15/2022 NU"; //Label for the non-accidentals background spectra
char muon_label[64]="Olivia 3/15/2022 NU"; //Label for the livetime/muon table

//const double livetime[8] = {72283163195734976, 81341302266285248, 94242900368444256, 84031181570224704, 143447962585023856, 143434576886712080, 143413538529066832, 127989198069118256}; //livetime in ns
const double livetime[8] = {72244284936199984, 81301206745935696, 94204225977453632, 84001966789146096, 143440456179300048, 143430087896206080, 143408721398819648, 127989198069118256}; //livetime in ns for NU

const double tarpRelError = 0.0037; //target proton relative error

//Rates:
const double muon_rate[8] = {200.32, 200.32, 150.08, 149.80, 15.748, 15.748, 15.748, 15.747}; //rate in Hz
const double muon_eff[8] = {0.5442,0.5415,0.6262,0.6256,0.9543,0.9543,0.9541,0.9545};
const double fastn_rate[3] = {2.16,1.61,0.14}; //rate of fast neutrons in 1/day NU
const double fastn_rate_err[3] = {0.14,0.08,0.02}; //error of rate of fast neutrons in 1/day NU
const double li9_rate[3] = {1.95,1.91,0.14}; //rate of li9 in 1/day NU
const double li9_rate_err[3] = {0.83,0.82,0.05}; //error of rate of li9 in 1/day NU
const double radn_rate[1] = {0.151}; //rate of rad n in 1/day
const double radn_rate_err[1] = {0.040}; //error of rate of rad n in 1/day
const double amc_rate[8] = {0.04,0.04,0.04,0.03,0.02,0.01,0.01,0.01}; //rate of amc in 1/day
const double amc_rate_err[8] = {0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01}; //error of rate of amc in 1/day

void muon(){

	for(int iad=0; iad<8; iad++){
		cout << "INSERT INTO muon_rates VALUES (" << EH[iad] << "," << AD[iad] << ",\"" << muon_label << "\"," << muon_rate[iad]*livetime[iad]*1.e-9 << "," << livetime[iad] << "," << muon_rate[iad] << "," << muon_eff[iad] << ");" << endl;
	}
}


void prompt(){

	TFile* outfile=new TFile("../nH_files/toyMC_spectra/ibd_eprompt_shapes_8ad.root", "RECREATE");

	for(int iad=0; iad<8; iad++){

		TFile *in_file = new TFile(Form("../nH_files/TotaledPlots_NU_EH%d_1500.root",EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_DT2000_3sig_ad%d",AD[iad]));
		
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


		cout << "INSERT INTO num_coincidences VALUES (" << EH[iad] << "," << AD[iad] << ",1,\"[";
		for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
			cout << (h_rebinned_spectrum->GetBinContent(iBin));
			if(iBin != h_rebinned_spectrum->GetNbinsX()) cout << ",";
		}
		cout << "]\",\"" << prompt_label << "\");" << endl;
		
//		cout << "\tEH" << EH[iad] << " AD" << AD[iad] << "\t" << h_rebinned_spectrum->Integral() << "\t+-\t" << sqrt(h_rebinned_spectrum->Integral())<< endl;

		outfile->cd();
			h_rebinned_spectrum->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_rebinned_spectrum->GetYaxis()->SetTitle("Counts");
			h_rebinned_spectrum->Write();
			
			h_fine_spectrum->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_fine_spectrum->GetYaxis()->SetTitle("Counts");
			h_fine_spectrum->Write();

		in_file->Close();
	}
	outfile->Close();

}

void acc_spectra(){

	TFile* outfile=new TFile("../nH_files/toyMC_spectra/accidental_eprompt_shapes_8ad.root", "RECREATE");

	for(int iad=0; iad<8; iad++){
		TFile *in_file = new TFile(Form("../nH_files/TotaledSingles_NU_1500_EH%d.root",EH[iad]));
		TH1F *in_hist = (TH1F*)in_file->Get(Form("h_total_prompt_energy_DT2000_3sig_scaled_ad%d",AD[iad]));

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

		outfile->cd();
			h_rebinned_spectrum->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_rebinned_spectrum->GetYaxis()->SetTitle("Counts");
			h_rebinned_spectrum->Write();
			
			h_fine_spectrum->GetXaxis()->SetTitle("Prompt Energy [MeV]");
			h_fine_spectrum->GetYaxis()->SetTitle("Counts");
			h_fine_spectrum->Write();

		in_file->Close();
	cout << endl;
	}
	outfile->Close();
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

	TFile* outfile=new TFile("../nH_files/toyMC_spectra/P15A_fn_spectrum.root", "RECREATE");

	TFile *in_file = new TFile("./FastNshape.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("OWP_EH%d_nH_Ep1.5_DT0.8m_Ep_EpLt12MeV",eh));

	char name[64];
	const int numBins = 34;
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	sprintf(name, "h_rebinned_spectrum");
	TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);
	
	TH1F* h_fn[8];
	TH1F* h_fn_fine[8];
	for(int i=0; i<8; i++){
		sprintf(name, "h_%dAD_fn",i+1);
		h_fn[i]=new TH1F(name,name,numBins, binEdges);
		sprintf(name, "h_%dAD_fn_fine",i+1);
		h_fn_fine[i]=new TH1F(name,name,120,0,12.);
	}


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
		if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
		h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		for(int i=0; i<8; i++){
			if((in_hist->GetBinCenter(iBin))>1.5) h_fn[i]->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
			if((in_hist->GetBinCenter(iBin))>1.5) h_fn_fine[i]->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		}
	} 


	for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO fast_neutron_spectrum VALUES (\"" << bgSpectra_label << "\",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
	}
	
	outfile->cd();
		for(int i=0; i<8; i++){
			h_fn_fine[i]->GetXaxis()->SetTitle("Energy [MeV]");
			h_fn_fine[i]->GetYaxis()->SetTitle("Counts");
			h_fn_fine[i]->Write();

			h_fn[i]->GetXaxis()->SetTitle("Energy [MeV]");
			h_fn[i]->GetYaxis()->SetTitle("Counts");
			h_fn[i]->Write();
		}
	
	in_file->Close();
	outfile->Close();
}


void fastN_counts(){
	for(int iad=0; iad<8; iad++){
		double temp_livetime_day = livetime[2*(EH[iad]-1)+AD[iad]-1] * 1.e-9/(60*60*24); //converting the livetime from ns to day
		cout << "INSERT INTO bg_counts VALUES (\"" << bgCount_label << "\"," << EH[iad] << "," << AD[iad] << ",\"fast-neutron\"," << fastn_rate[EH[iad]-1]*temp_livetime_day << "," << fastn_rate_err[EH[iad]-1]*temp_livetime_day << ");" << endl;
	}
}

void li9_spectrum(){

	TFile* outfile=new TFile("../nH_files/toyMC_spectra/8he9li_nominal_spectrum.root", "RECREATE");

	TFile *in_file = new TFile("./Li9He8shape.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("Li9 prompt energy"));

        char name[64];
	const int numBins = 34;
	double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
	sprintf(name, "h_rebinned_spectrum");
	TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);

	TH1F* h_fine_spectrum=new TH1F("h_nominal","h_nominal",120,0,12.);


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
		if((in_hist->GetBinCenter(iBin))>binEdges[numBins]) break;
		h_rebinned_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
		if((in_hist->GetBinCenter(iBin))>1.5)h_fine_spectrum->Fill(in_hist->GetBinCenter(iBin),in_hist->GetBinContent(iBin));
	} 


	for(int iBin=1; iBin < (h_rebinned_spectrum->GetNbinsX())+1; iBin++){
		cout << "INSERT INTO li9_spectrum VALUES (\"" << bgSpectra_label << "\",1," << iBin-1 << "," << (h_rebinned_spectrum->GetBinContent(iBin))/(h_rebinned_spectrum->Integral()) << ");" << endl;
	}
	
	outfile->cd();
		h_fine_spectrum->GetXaxis()->SetTitle("Energy [MeV]");
		h_fine_spectrum->GetYaxis()->SetTitle("Counts");
		h_fine_spectrum->Write();

	in_file->Close();

	outfile->Close();
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




void toyMC_li9(){

	TFile* outfile=new TFile("../nH_files/toyMC_spectra/8he9li_nominal_spectrum.root", "RECREATE");

	TFile *in_file = new TFile("../nH_files/toyMC_spectra/8he9li_nominal_spectrum_original.root");
	TH1F *in_hist = (TH1F*)in_file->Get(Form("h_nominal"));


	for(int iBin=1; iBin < (in_hist->GetNbinsX())+2; iBin++){
		if((in_hist->GetBinCenter(iBin))<1.5) in_hist->SetBinContent(iBin,0);
	} 
	
	outfile->cd();
		in_hist->GetXaxis()->SetTitle("Energy [MeV]");
		in_hist->GetYaxis()->SetTitle("Counts");
		in_hist->Write();

	in_file->Close();

	outfile->Close();
}

void toyMC_amc(){

	TFile* outfile=new TFile("../nH_files/toyMC_spectra/amc_spectrum.root", "RECREATE");

	TFile *in_file = new TFile("../nH_files/toyMC_spectra/amc_spectrum_original.root");
	TH1D *in_hist_toy = (TH1D*)in_file->Get(Form("h_toy"));
	TH1D *in_hist_rebin = (TH1D*)in_file->Get(Form("h_rebin"));
	TH1D *in_hist_fit = (TH1D*)in_file->Get(Form("h_sum_fit"));
	TH1F *in_hist_corrAmc = (TH1F*)in_file->Get(Form("hCorrAmCPromptSpec"));
	TF1 *expo_fit = (TF1*)in_file->Get(Form("expo"));



	for(int iBin=1; iBin < (in_hist_toy->GetNbinsX())+2; iBin++){
		if((in_hist_toy->GetBinCenter(iBin))<1.5) in_hist_toy->SetBinContent(iBin,0);
	}
	for(int iBin=1; iBin < (in_hist_rebin->GetNbinsX())+2; iBin++){
		if((in_hist_rebin->GetBinCenter(iBin))<1.5) in_hist_rebin->SetBinContent(iBin,0);
	}
	for(int iBin=1; iBin < (in_hist_fit->GetNbinsX())+2; iBin++){
		if((in_hist_fit->GetBinCenter(iBin))<1.5) in_hist_fit->SetBinContent(iBin,0);
	}
	for(int iBin=1; iBin < (in_hist_corrAmc->GetNbinsX())+2; iBin++){
		if((in_hist_corrAmc->GetBinCenter(iBin))<1.5) in_hist_corrAmc->SetBinContent(iBin,0);
	}



	
	outfile->cd();
		in_hist_toy->GetXaxis()->SetTitle("Energy [MeV]");
		in_hist_toy->GetYaxis()->SetTitle("Counts");
		in_hist_toy->Write();

		in_hist_rebin->GetXaxis()->SetTitle("Energy [MeV]");
		in_hist_rebin->GetYaxis()->SetTitle("Counts");
		in_hist_rebin->Write();

		in_hist_fit->GetXaxis()->SetTitle("Energy [MeV]");
		in_hist_fit->GetYaxis()->SetTitle("Counts");
		in_hist_fit->Write();

		in_hist_corrAmc->GetXaxis()->SetTitle("Energy [MeV]");
		in_hist_corrAmc->GetYaxis()->SetTitle("Counts");
		in_hist_corrAmc->Write();

		expo_fit->GetXaxis()->SetTitle("Energy [MeV]");
		expo_fit->GetYaxis()->SetTitle("Counts");
		expo_fit->Write();

	in_file->Close();

	outfile->Close();
}


void toyMC_alphaN(){

	TFile* outfile=new TFile("../nH_files/toyMC_spectra/result-DocDB9667.root", "RECREATE");

//	TFile *in_file = new TFile("../nH_files/toyMC_spectra/result-DocDB9667_alpha-n_original.root");
	
	for(int i=0; i<8; i++){
	
		char name[64];
		const int numBins = 34;
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		sprintf(name, "AD%d",i+1);
		TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);

		double binContent[numBins] = {0.183738673305789,0.1473001973115,0.131236675163554,0.210434166080378,0.173701502451188,0.042442682095793,0.020343358438646,0.01580270431453,0.013037100196197,0.009911021529084,0.007997333203122,0.007721738444439,0.006984576579051,0.005077747513876,0.003918246185711,0.005410842172025,0.004757312973606,0.003341925152022,0.00122860310898,0.000552871399041,0.000491441243592,0.000245720621796,0.000491441243592,0.000307150777245,0.000737161865388,0.000184290466347,0.000184290466347,0.000122860310898,0.000122860310898,0.000122860310898,0.000637751189138,0.000307150777245,0.000245720621796,0.000860022176286};
		
		for(int iBin=0; iBin<numBins; iBin++){
			cout << h_rebinned_spectrum->GetXaxis()->FindBin(binEdges[iBin+1]-(binEdges[iBin+1]-binEdges[iBin])/2.) << "\t" << binEdges[iBin+1]-(binEdges[iBin+1]-binEdges[iBin])/2. << endl;
			h_rebinned_spectrum->SetBinContent(h_rebinned_spectrum->GetXaxis()->FindBin(binEdges[iBin+1]-(binEdges[iBin+1]-binEdges[iBin])/2.),binContent[iBin]);
		}
		
		outfile->cd();
			h_rebinned_spectrum->GetXaxis()->SetTitle("Energy [MeV]");
			h_rebinned_spectrum->GetYaxis()->SetTitle("Counts");
			h_rebinned_spectrum->Write();
	}

//	in_file->Close();

	outfile->Close();
}

/*void toyMC_radN(){

	TFile* outfile=new TFile("../nH_files/toyMC_spectra/MuonDecaySpec.root", "RECREATE");
	
	for(int i=0; i<3; i++){

		char name[64];
		const int numBins = 34;
		double binEdges[numBins+1] = {1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7,5.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,12.};
		sprintf(name, "MdSpec_EH%d",i+1);
		TH1F* h_rebinned_spectrum=new TH1F(name,name,numBins, binEdges);

		double binContent[numBins] = {0.183738673305789,0.1473001973115,0.131236675163554,0.210434166080378,0.173701502451188,0.042442682095793,0.020343358438646,0.01580270431453,0.013037100196197,0.009911021529084,0.007997333203122,0.007721738444439,0.006984576579051,0.005077747513876,0.003918246185711,0.005410842172025,0.004757312973606,0.003341925152022,0.00122860310898,0.000552871399041,0.000491441243592,0.000245720621796,0.000491441243592,0.000307150777245,0.000737161865388,0.000184290466347,0.000184290466347,0.000122860310898,0.000122860310898,0.000122860310898,0.000637751189138,0.000307150777245,0.000245720621796,0.000860022176286};
		
		for(int iBin=0; iBin<numBins; iBin++){
			cout << h_rebinned_spectrum->GetXaxis()->FindBin(binEdges[iBin+1]-(binEdges[iBin+1]-binEdges[iBin])/2.) << "\t" << binEdges[iBin+1]-(binEdges[iBin+1]-binEdges[iBin])/2. << endl;
			h_rebinned_spectrum->SetBinContent(h_rebinned_spectrum->GetXaxis()->FindBin(binEdges[iBin+1]-(binEdges[iBin+1]-binEdges[iBin])/2.),binContent[iBin]);
		}

		
		outfile->cd();
			h_rebinned_spectrum->GetXaxis()->SetTitle("Energy [MeV]");
			h_rebinned_spectrum->GetYaxis()->SetTitle("Counts");
			h_rebinned_spectrum->GetXaxis()->SetRangeUser(0,12.);
			h_rebinned_spectrum->Write();
	}

	outfile->Close();
}*/

