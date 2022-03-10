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
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLine.h"
#include <typeinfo>

using namespace std;

void calcEff_tD(){
	int num_methods = 3;
	int num_timeCuts = 4; //2000, 800, 600, 400
	int num_distCuts = 5; //2.0, 1.5, 1., 0.75, 0.5
	int maxAD = 8;
	int num_scales = 2;//1 for normalized, 2 for rate corrected

//	int timeCuts[2] = {2000, 400};
//	double distCuts[2] = {2.0, 0.5};

	int timeCuts[4] = {2000, 800, 600, 400};
	double distCuts[5] = {2.0, 1.5, 1., 0.75, 0.5};

	double timeEff[num_timeCuts][num_distCuts][maxAD][num_scales];
	double distEff[num_timeCuts][num_distCuts][maxAD][num_scales];
	double combinedEff[num_methods][num_timeCuts*num_distCuts][maxAD][num_scales];

	double avg_timeEff[num_timeCuts][num_distCuts][num_scales];
	double avg_distEff[num_timeCuts][num_distCuts][num_scales];
	double avg_combinedEff[num_methods][num_timeCuts*num_distCuts][num_scales];
	double avg_combinedEff_farAvg[num_methods][num_timeCuts*num_distCuts][num_scales];

	double time_uncertainty[num_methods][num_timeCuts*num_distCuts][num_scales];
	double dist_uncertainty[num_methods][num_timeCuts*num_distCuts][num_scales];
	double eff_uncertainty[num_methods][num_timeCuts*num_distCuts][num_scales];

	for(int iTime=0; iTime<num_timeCuts; iTime++){
		for(int iDist=0; iDist<num_distCuts; iDist++){
			for(int iScale=0; iScale<num_scales; iScale++){
				avg_timeEff[iTime][iDist][iScale]=0;
				avg_distEff[iTime][iDist][iScale]=0;
				avg_combinedEff[0][num_distCuts*iTime+iDist][iScale]=0;
				avg_combinedEff[1][num_distCuts*iTime+iDist][iScale]=0;
				avg_combinedEff[2][num_distCuts*iTime+iDist][iScale]=0;
				avg_combinedEff_farAvg[0][num_distCuts*iTime+iDist][iScale]=0;
				avg_combinedEff_farAvg[1][num_distCuts*iTime+iDist][iScale]=0;
				avg_combinedEff_farAvg[2][num_distCuts*iTime+iDist][iScale]=0;

				time_uncertainty[0][num_distCuts*iTime+iDist][iScale]=0;
				time_uncertainty[1][num_distCuts*iTime+iDist][iScale]=0;
				dist_uncertainty[0][num_distCuts*iTime+iDist][iScale]=0;
				dist_uncertainty[1][num_distCuts*iTime+iDist][iScale]=0;
				eff_uncertainty[0][num_distCuts*iTime+iDist][iScale]=0;
				eff_uncertainty[1][num_distCuts*iTime+iDist][iScale]=0;
				eff_uncertainty[2][num_distCuts*iTime+iDist][iScale]=0;

				for(int iad=0; iad<maxAD; iad++){
					timeEff[iTime][iDist][iad][iScale]=0;
					distEff[iTime][iDist][iad][iScale]=0;
					combinedEff[0][num_distCuts*iTime+iDist][iad][iScale]=0;
					combinedEff[1][num_distCuts*iTime+iDist][iad][iScale]=0;
					combinedEff[2][num_distCuts*iTime+iDist][iad][iScale]=0;
				}
			}
		}
	}


//MAKING THE PLOTS
	char name[64];
	TMultiGraph* ADs_compared[num_timeCuts*num_distCuts][num_scales];
	TGraph* ad_comp[num_methods][num_timeCuts*num_distCuts][num_scales];
	TMultiGraph* ADs_compared_farAvg[num_timeCuts*num_distCuts][num_scales];
	TGraph* ad_comp_farAvg[num_methods][num_timeCuts*num_distCuts][num_scales];
	for(int iTime = 0; iTime < num_timeCuts; iTime++){
		for(int iDist = 0; iDist < num_distCuts; iDist++){
			for(int iScale=0; iScale < num_scales; iScale++){
				for(int iMethod = 0; iMethod<num_methods; iMethod++){
					ad_comp[iMethod][num_distCuts*iTime+iDist][iScale] = new TGraph();
					if(iScale==0) sprintf(name, "method%d_t%d_d%.2f_norm", iMethod+1,timeCuts[iTime],distCuts[iDist]);
					if(iScale==1) sprintf(name, "method%d_t%d_d%.2f_rate", iMethod+1,timeCuts[iTime],distCuts[iDist]);
					ad_comp[iMethod][num_distCuts*iTime+iDist][iScale]->SetName(name);

					ad_comp_farAvg[iMethod][num_distCuts*iTime+iDist][iScale] = new TGraph();
					if(iScale==0) sprintf(name, "method%d_t%d_d%.2f_norm", iMethod+1,timeCuts[iTime],distCuts[iDist]);
					if(iScale==1) sprintf(name, "method%d_t%d_d%.2f_rate", iMethod+1,timeCuts[iTime],distCuts[iDist]);
					ad_comp_farAvg[iMethod][num_distCuts*iTime+iDist][iScale]->SetName(name);
				}
				if(iScale==0) sprintf(name, "ADs_compared_t%d_d%.2f_norm",timeCuts[iTime],distCuts[iDist]);
				if(iScale==1) sprintf(name, "ADs_compared_t%d_d%.2f_rate",timeCuts[iTime],distCuts[iDist]);
				ADs_compared[num_distCuts*iTime+iDist][iScale] = new TMultiGraph(name,name);

				if(iScale==0) sprintf(name, "ADs_compared_farAvg_t%d_d%.2f_norm",timeCuts[iTime],distCuts[iDist]);
				if(iScale==1) sprintf(name, "ADs_compared_farAvg_t%d_d%.2f_rate",timeCuts[iTime],distCuts[iDist]);
				ADs_compared_farAvg[num_distCuts*iTime+iDist][iScale] = new TMultiGraph(name,name);
			}
		}
	}





//STARTING IN ON THE CALCULATIONS:

		char inputname[64];
		sprintf(inputname,"./combinedEfficiency/distEff2000_results.root");
		TFile *distFile = new TFile(inputname);

		sprintf(inputname,"./combinedEfficiency/TimeCutEffs.root");
		TFile *timeFile = new TFile(inputname);

		//Getting the time efficiencies using various distance cuts

	cout << "Getting the time cut efficiencies" << endl;

		char timeGraphName[64];
		double x = 0;
		double temp_timeEff = 0;
		char distTimeGraphName[64];

		for(int iScale=0; iScale<num_scales; iScale++){
			for(int iDist = 0; iDist < num_distCuts; iDist++){
				//AD time efficiency:
				for(int iad = 0; iad<maxAD; iad++){
					if(distCuts[iDist]==2) sprintf(timeGraphName,"eff_data2_ad%d",iad+1);
					if(distCuts[iDist]==1.5) sprintf(timeGraphName,"eff_data15_ad%d",iad+1);
					if(distCuts[iDist]==1.0) sprintf(timeGraphName,"eff_data1_ad%d",iad+1);
					if(distCuts[iDist]==0.75) sprintf(timeGraphName,"eff_data075_ad%d",iad+1);
					if(distCuts[iDist]==0.5) sprintf(timeGraphName,"eff_data05_ad%d",iad+1);

					TGraph* timeGraph = (TGraph*)timeFile->Get(timeGraphName);

					for(int iTime = 0; iTime < num_timeCuts; iTime++){
						for(int iPoint = 0; iPoint < 2000; iPoint++){
							timeGraph->GetPoint(iPoint,x,temp_timeEff);
							if(x>=timeCuts[iTime]) break;
						}


						timeEff[iTime][iDist][iad][iScale] = temp_timeEff;
					}
				}//end of AD gathering
				

				//average time efficiency:
				if(distCuts[iDist]==2) sprintf(timeGraphName,"avgEff_data2");
				if(distCuts[iDist]==1.5) sprintf(timeGraphName,"avgEff_data15");
				if(distCuts[iDist]==1.0) sprintf(timeGraphName,"avgEff_data1");
				if(distCuts[iDist]==0.75) sprintf(timeGraphName,"avgEff_data075");
				if(distCuts[iDist]==0.5) sprintf(timeGraphName,"avgEff_data05");

				for(int iTime = 0; iTime < num_timeCuts; iTime++){
					TGraph* timeGraph = (TGraph*)timeFile->Get(timeGraphName);
					for(int iPoint = 0; iPoint < 2000; iPoint++){
						timeGraph->GetPoint(iPoint,x,temp_timeEff);
						if(x>=timeCuts[iTime]) break;
					}


					avg_timeEff[iTime][iDist][iScale] = temp_timeEff;
				}//end of average gathering
				

			} //end of gathering time cut efficiencies
		}


		//Getting the distance efficiencies using various time cuts

	cout << "Getting the distance cut efficiencies" << endl;

		char distGraphName[64];
		double temp_distEff = 0;

		for(int iScale=0; iScale<num_scales; iScale++){
			for(int iTime = 0; iTime < num_timeCuts; iTime++){
				//AD distance efficiency:
				for(int iad = 0; iad<maxAD; iad++){
					if(iScale == 0){
						sprintf(distGraphName,"eff_norm_%d_ad%d",timeCuts[iTime],iad+1);
						if(timeCuts[iTime] == 2000) sprintf(distGraphName,"eff_norm_ad%d",iad+1);
					}
					else{
						sprintf(distGraphName,"eff_rate_%d_ad%d",timeCuts[iTime],iad+1);
						if(timeCuts[iTime] == 2000) sprintf(distGraphName,"eff_rate_ad%d",iad+1);
					}

					TGraph* distGraph = (TGraph*)distFile->Get(distGraphName);

					for(int iDist = 0; iDist < num_distCuts; iDist++){
						for(int iPoint = 0; iPoint < 700; iPoint++){
							distGraph->GetPoint(iPoint,x,temp_distEff);
							if(x>=distCuts[iDist]){
								break;
							}
						}


						distEff[iTime][iDist][iad][iScale] = temp_distEff;
					}
				}//end of AD gathering
				

				//average distance efficiency:
				if(iScale == 0){
					sprintf(distGraphName,"eff_norm_%d_avg",timeCuts[iTime]);
					if(timeCuts[iTime] == 2000) sprintf(distGraphName,"eff_norm_avg");
				}
				else{
					sprintf(distGraphName,"eff_rate_%d_avg",timeCuts[iTime]);
					if(timeCuts[iTime] == 2000) sprintf(distGraphName,"eff_rate_avg");
				}
	//	cout << "Getting graph: " << distGraphName << endl;
				for(int iDist = 0; iDist < num_distCuts; iDist++){
					TGraph* distGraph = (TGraph*)distFile->Get(distGraphName);
					for(int iPoint = 0; iPoint < 700; iPoint++){
						distGraph->GetPoint(iPoint,x,temp_distEff);
						if(x>=distCuts[iDist]) break;
					}


					avg_distEff[iTime][iDist][iScale] = temp_distEff;
	//				cout << iTime << "\t" << iDist << "\t" << timeCuts[iTime] << "\t" << x << "\t" << temp_distEff << "\t" << avg_distEff[iTime][iDist] << endl;	
				}//end of average gathering
				

			} //end of gathering time cut efficiencies
		}


	cout << endl << "Combining the efficiencies" << endl << endl;

	double temp_dist_uncertainty = 0;
	double temp_time_uncertainty = 0;
	double temp_dt_uncertainty = 0;
	double nXbins = 0;
	double nYbins = 0;
	double totXbins = 0;
	double totYbins = 0;
	double temp_x_center = 0;
	double temp_y_center = 0;

		//Combining them:
		for(int iMethod=1; iMethod<num_methods+1; iMethod++){
			for(int iScale=0; iScale<num_scales; iScale++){

				cout << "SCALE NUMBER: " << iScale+1 << endl;
				if(iScale==0) cout << "NORMALIZED!!" << endl;
				else cout << "RATE-CORRECTED!!" << endl;

				if(iMethod == 1){ //Method 1: distance cut first, then time cut second
					for(int iDist = 0; iDist < num_distCuts; iDist++){
						for(int iTime = 0; iTime < num_timeCuts; iTime++){
							for(int iad = 0; iad < maxAD; iad++){
								combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale] = distEff[0][iDist][iad][iScale]*timeEff[iTime][iDist][iad][iScale];
								ad_comp[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+0.98,100*combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale]);
								if(iad < 4) ad_comp_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+0.98,100*combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale]);
							} //end of AD loop

							for(int iad = 4; iad < 8; iad++){ //Far hall averaging
								avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] + combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale];
							}
							avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]/4.;
							for(int iad = 4; iad < 8; iad++){
								ad_comp_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+0.98,100*avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]);
							} //end of far hall averaging

							avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_distEff[0][iDist][iScale]*avg_timeEff[iTime][iDist][iScale];

						//Calculate the uncertainty:
							dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = (distEff[0][iDist][4][iScale]+distEff[0][iDist][5][iScale]+distEff[0][iDist][6][iScale]+distEff[0][iDist][7][iScale])/4-avg_distEff[0][iDist][iScale];
							if(dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] < 0.) dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = -1.*dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];
							temp_dist_uncertainty = dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];

							for(int iad = 0; iad < 4; iad++){
								temp_dist_uncertainty = distEff[0][iDist][iad][iScale]-avg_distEff[0][iDist][iScale];
								if(temp_dist_uncertainty < 0.) temp_dist_uncertainty = -1.*temp_dist_uncertainty;
								if(temp_dist_uncertainty > dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale]) dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = temp_dist_uncertainty;
							}

							time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = (timeEff[iTime][iDist][4][iScale]+timeEff[iTime][iDist][5][iScale]+timeEff[iTime][iDist][6][iScale]+timeEff[iTime][iDist][7][iScale])/4-avg_timeEff[iTime][iDist][iScale];
							if(time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] < 0.) time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = -1.*time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];
							temp_time_uncertainty = time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];

							for(int iad = 0; iad < 4; iad++){
								temp_time_uncertainty = timeEff[iTime][iDist][iad][iScale]-avg_timeEff[iTime][iDist][iScale];
								if(temp_time_uncertainty < 0.) temp_time_uncertainty = -1.*temp_time_uncertainty;
								if(temp_time_uncertainty > time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale]) time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = temp_time_uncertainty;
							}

							eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale]*sqrt(pow((dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale])/avg_distEff[0][iDist][iScale],2)+pow(time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale]/avg_timeEff[iTime][iDist][iScale],2));
						//end of calculating the uncertainty for method 1


							

							cout << timeCuts[iTime] << "\t" << distCuts[iDist] << "\t" << 100*avg_timeEff[iTime][iDist][iScale] << "\t" << 100*avg_distEff[0][iDist][iScale] << "\t" << 100*avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale] << "\t" << 100*time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] << "\t" << 100*dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] << "\t" << 100*eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] << endl;
						}
					}
				}//end of method 1

			cout << endl << endl;

				if(iMethod == 2){ //Method 2: time cut first, then distance cut second
					for(int iDist = 0; iDist < num_distCuts; iDist++){
						for(int iTime = 0; iTime < num_timeCuts; iTime++){
							for(int iad = 0; iad < maxAD; iad++){
								combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale] = timeEff[iTime][0][iad][iScale]*distEff[iTime][iDist][iad][iScale];
								ad_comp[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+1.00,100*combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale]);
								if(iad < 4) ad_comp_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+1.02,100*combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale]);
							} //end of AD loop

							for(int iad = 4; iad < 8; iad++){ //Far hall averaging
								avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] + combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale];
							}
							avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]/4.;
							for(int iad = 4; iad < 8; iad++){
								ad_comp_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+1.00,100*avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]);
							} //end of far hall averaging

							avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_timeEff[iTime][0][iScale]*avg_distEff[iTime][iDist][iScale];
						//Calculate the uncertainty here
							dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = (distEff[iTime][iDist][4][iScale]+distEff[iTime][iDist][5][iScale]+distEff[iTime][iDist][6][iScale]+distEff[iTime][iDist][7][iScale])/4-avg_distEff[iTime][iDist][iScale];
							if(dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] < 0.) dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = -1.*dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];
							temp_dist_uncertainty = dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];

							for(int iad = 0; iad < 4; iad++){
								temp_dist_uncertainty = distEff[iTime][iDist][iad][iScale]-avg_distEff[iTime][iDist][iScale];
								if(temp_dist_uncertainty < 0.) temp_dist_uncertainty = -1.*temp_dist_uncertainty;
								if(temp_dist_uncertainty > dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale]) dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = temp_dist_uncertainty;
							}

							time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = (timeEff[iTime][0][4][iScale]+timeEff[iTime][0][5][iScale]+timeEff[iTime][0][6][iScale]+timeEff[iTime][0][7][iScale])/4-avg_timeEff[iTime][0][iScale];
							if(time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] < 0.) time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = -1.*time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];
							temp_time_uncertainty = time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];

							for(int iad = 0; iad < 4; iad++){
								temp_time_uncertainty = timeEff[iTime][0][iad][iScale]-avg_timeEff[iTime][0][iScale];
								if(temp_time_uncertainty < 0.) temp_time_uncertainty = -1.*temp_time_uncertainty;
								if(temp_time_uncertainty > time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale]) time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = temp_time_uncertainty;
							}

							eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale]*sqrt(pow((dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale])/avg_distEff[iTime][iDist][iScale],2)+pow(time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale]/avg_timeEff[iTime][0][iScale],2));
						//End of calculating the uncertainty for method 2


							cout << timeCuts[iTime] << "\t" << distCuts[iDist] << "\t" << 100*avg_timeEff[iTime][0][iScale] << "\t" << 100*avg_distEff[iTime][iDist][iScale] << "\t" << 100*avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale]  << "\t" << 100*time_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] << "\t" << 100*dist_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] << "\t" << 100*eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] << endl;

						}
					}
				}//end of method 2

				if(iMethod == 3){
					for(int iad = 0; iad < maxAD; iad++){
						if(iScale == 0) sprintf(distTimeGraphName,"h_sub_distVStime_norm_ad%d",iad+1);
						else sprintf(distTimeGraphName,"h_sub_distVStime_rate_ad%d",iad+1);
						TH2F* distTimeGraph = (TH2F*)distFile->Get(distTimeGraphName);

						//finding the x or y bins to go to for 100% efficiency
						for(int iX = 0; iX < 2000; iX++){
							temp_x_center = distTimeGraph->GetXaxis()->GetBinCenter(iX);
							if(temp_x_center > timeCuts[0]) break;
							totXbins = iX;
						}
						for(int iY = 0; iY < 700; iY++){
							temp_y_center = distTimeGraph->GetYaxis()->GetBinCenter(iY);
							if(temp_y_center > distCuts[0]) break;
							totYbins = iY;
						}

						for(int iDist = 0; iDist < num_distCuts; iDist++){
							for(int iTime = 0; iTime < num_timeCuts; iTime++){
								temp_x_center = 0;
								temp_y_center = 0;
								//Getting the x or y bins to get cut limits
								for(int iX = 0; iX < 2000; iX++){
									temp_x_center = distTimeGraph->GetXaxis()->GetBinCenter(iX);
									if(temp_x_center > timeCuts[iTime]) break;
									nXbins = iX;
								}
								for(int iY = 0; iY < 700; iY++){
									temp_y_center = distTimeGraph->GetYaxis()->GetBinCenter(iY);
									if(temp_y_center > distCuts[iDist]) break;
									nYbins = iY;
								}

								//calculating the efficiency:
								combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale] = (distTimeGraph->Integral(0,nXbins,0,nYbins))/(distTimeGraph->Integral(0,totXbins,0,totYbins));
								ad_comp[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+1.02,100*combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale]);
								if(iad < 4) ad_comp_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+1.02,100*combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale]);
							}
						}
					}//end of AD loop

					for(int iDist = 0; iDist < num_distCuts; iDist++){
						for(int iTime = 0; iTime < num_timeCuts; iTime++){
							for(int iad = 0; iad < maxAD; iad++){
								avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale] + combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale];
							}
							//calculating the average of the ADs:
							avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale]/8.;

							for(int iad = 4; iad < 8; iad++){ //Averaging the far hall
								avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] + combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale];
							}
							avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale] = avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]/4.;
							for(int iad = 4; iad < 8; iad++){
								ad_comp_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]->SetPoint(iad,iad+1.02,100*avg_combinedEff_farAvg[iMethod-1][num_distCuts*iTime+iDist][iScale]);
							} //end of far hall averaging

							//calculating the uncertainty:
							eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = (combinedEff[iMethod-1][num_distCuts*iTime+iDist][4][iScale]+combinedEff[iMethod-1][num_distCuts*iTime+iDist][5][iScale]+combinedEff[iMethod-1][num_distCuts*iTime+iDist][6][iScale]+combinedEff[iMethod-1][num_distCuts*iTime+iDist][7][iScale])/4-avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale];
							if(eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] < 0.) eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = -1.*eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];
							temp_dt_uncertainty = eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale];

							for(int iad = 0; iad < 4; iad++){
								temp_dt_uncertainty = combinedEff[iMethod-1][num_distCuts*iTime+iDist][iad][iScale]-avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale];
								if(temp_dt_uncertainty < 0.) temp_dt_uncertainty = -1.*temp_dt_uncertainty;
								if(temp_dt_uncertainty > eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale]) eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] = temp_dt_uncertainty;
							}

							cout << timeCuts[iTime] << "\t" << distCuts[iDist] << "\t" << 100*avg_combinedEff[iMethod-1][num_distCuts*iTime+iDist][iScale]  << "\t" << 100*eff_uncertainty[iMethod-1][num_distCuts*iTime+iDist][iScale] << endl;
						}
					}
					

				}// end of method 3

			}// end of scale loop
		}//end of methods 1 and 2 to combine them


	for(int iScale=0; iScale<num_scales; iScale++){
		cout << endl << endl;
		if(iScale==0) cout << "NORMALIZED!!" << endl;
		else cout << "RATE-CORRECTED!!" << endl;

		for(int iTime = 0; iTime < num_timeCuts; iTime++){
			for(int iDist = 0; iDist < num_distCuts; iDist++){
//				cout << timeCuts[iTime] << "\t" << distCuts[iDist] << "\t" << 100*avg_combinedEff[0][num_distCuts*iTime+iDist][iScale] << " +/- " << 100*eff_uncertainty[0][num_distCuts*iTime+iDist][iScale] << "\t\t" << 100*avg_combinedEff[1][num_distCuts*iTime+iDist][iScale] << " +/- " << 100*eff_uncertainty[1][num_distCuts*iTime+iDist][iScale] << "\t\t" << 100*(avg_combinedEff[1][num_distCuts*iTime+iDist][iScale]-avg_combinedEff[0][num_distCuts*iTime+iDist][iScale]) << endl;

				cout << timeCuts[iTime] << "\t" << distCuts[iDist] << "\t" << 100*avg_combinedEff[0][num_distCuts*iTime+iDist][iScale] << " +/- " << 100*eff_uncertainty[0][num_distCuts*iTime+iDist][iScale] << "\t\t" << 100*avg_combinedEff[1][num_distCuts*iTime+iDist][iScale] << " +/- " << 100*eff_uncertainty[1][num_distCuts*iTime+iDist][iScale] << "\t\t" << 100*avg_combinedEff[2][num_distCuts*iTime+iDist][iScale] << " +/- " << 100*eff_uncertainty[2][num_distCuts*iTime+iDist][iScale] << endl;
			}
		}

	}


	int colors[3] = {kBlack, kRed, kBlue};

	//MAKE PLOTS AND SAVE THEM
        char outputname[64];
        sprintf(outputname,"./combinedEfficiency/eff_tD_results.root");
	TFile* outfile=new TFile(outputname, "RECREATE");
		outfile->cd();

	double temp_low = 0;
	double temp_high = 0;
	double y = 0;
		for(int iTime = 0; iTime < num_timeCuts; iTime++){
			for(int iDist = 0; iDist < num_distCuts; iDist++){
				for(int iScale=0; iScale < num_scales; iScale++){
					temp_low = 0;
					temp_high = 0;
					ADs_compared[num_distCuts*iTime+iDist][iScale]->GetXaxis()->SetTitle("AD Number");
					ADs_compared[num_distCuts*iTime+iDist][iScale]->GetYaxis()->SetTitle("Combined Efficiency [%]");
					ADs_compared[num_distCuts*iTime+iDist][iScale]->GetXaxis()->SetLimits(0,9);

					for(int iMethod = 0; iMethod<num_methods; iMethod++){
						ad_comp[iMethod][num_distCuts*iTime+iDist][iScale]->SetMarkerColor(colors[iMethod]);
						ad_comp[iMethod][num_distCuts*iTime+iDist][iScale]->SetMarkerStyle(20);
						ad_comp[iMethod][num_distCuts*iTime+iDist][iScale]->SetLineColor(colors[iMethod]);
						ADs_compared[num_distCuts*iTime+iDist][iScale]->Add(ad_comp[iMethod][num_distCuts*iTime+iDist][iScale]);

						for(int iPoint = 0; iPoint < maxAD; iPoint++){
							ad_comp[iMethod][num_distCuts*iTime+iDist][iScale]->GetPoint(iPoint,x,y);
							if(iPoint == 0 && iMethod == 0){
								temp_low = y;
								temp_high = y;
							}
							else if(y < temp_low) temp_low = y;
							else if(y > temp_high) temp_high = y;
						}
					}
//					ADs_compared[num_distCuts*iTime+iDist][iScale]->GetYaxis()->SetRangeUser(temp_low-0.5,temp_high+0.5);
					ADs_compared[num_distCuts*iTime+iDist][iScale]->GetYaxis()->SetRangeUser((temp_low+temp_high)/2-3.,(temp_low+temp_high)/2+3.);
					ADs_compared[num_distCuts*iTime+iDist][iScale]->Write();
				}
			}
		}

	//Far averaged
		for(int iTime = 0; iTime < num_timeCuts; iTime++){
			for(int iDist = 0; iDist < num_distCuts; iDist++){
				for(int iScale=0; iScale < num_scales; iScale++){
					temp_low = 0;
					temp_high = 0;
					ADs_compared_farAvg[num_distCuts*iTime+iDist][iScale]->GetXaxis()->SetTitle("AD Number");
					ADs_compared_farAvg[num_distCuts*iTime+iDist][iScale]->GetYaxis()->SetTitle("Combined Efficiency [%]");
					ADs_compared_farAvg[num_distCuts*iTime+iDist][iScale]->GetXaxis()->SetLimits(0,9);

					for(int iMethod = 0; iMethod<num_methods; iMethod++){
						ad_comp_farAvg[iMethod][num_distCuts*iTime+iDist][iScale]->SetMarkerColor(colors[iMethod]);
						ad_comp_farAvg[iMethod][num_distCuts*iTime+iDist][iScale]->SetMarkerStyle(20);
						ad_comp_farAvg[iMethod][num_distCuts*iTime+iDist][iScale]->SetLineColor(colors[iMethod]);
						ADs_compared_farAvg[num_distCuts*iTime+iDist][iScale]->Add(ad_comp_farAvg[iMethod][num_distCuts*iTime+iDist][iScale]);
						for(int iPoint = 0; iPoint < maxAD; iPoint++){
							ad_comp_farAvg[iMethod][num_distCuts*iTime+iDist][iScale]->GetPoint(iPoint,x,y);
							if(iPoint == 0 && iMethod == 0){
								temp_low = y;
								temp_high = y;
							}
							else if(y < temp_low) temp_low = y;
							else if(y > temp_high) temp_high = y;
						}
					}
//					ADs_compared_farAvg[num_distCuts*iTime+iDist][iScale]->GetYaxis()->SetRangeUser(temp_low-0.5,temp_high+0.5);
					ADs_compared_farAvg[num_distCuts*iTime+iDist][iScale]->GetYaxis()->SetRangeUser((temp_low+temp_high)/2-3.,(temp_low+temp_high)/2+3.);
					ADs_compared_farAvg[num_distCuts*iTime+iDist][iScale]->Write();
				}
			}
		}


}


void DTcut(int Ep_cut, int weighted, int rateMod){

	char name[64];
	char distTimeGraphName[64];
	char IBDdistTimeGraphName[64];

	//Constants/variables needed
	int maxAD = 8;
	int num_scales = 2;//1 for normalized, 2 for rate corrected
		double DTcut = 800; //units of mm
	//	int Ep_cut = 0; //Is the 3.5 MeV Ep cut implemented?
	double yScale = 4;
	double yScale_farAvg = 3;
	double max_x_value = 2000;
	double max_y_value = 2.0;

	double time_center = 0;
	double dist_center = 0;
	double temp_value = 0;
	double combinedEfficiency[maxAD][num_scales];
	double statError[maxAD][num_scales];
	double errDenom[maxAD][num_scales];
	double sig1[maxAD][num_scales];
	double sig2[maxAD][num_scales];
	double avgFar[num_scales];
	double Sam_avgFar[num_scales];
	double Err_avgFar[num_scales];
	double SamErr_avgFar[num_scales];
	double average[num_scales];
	double Sam_average[num_scales];
	double spread[num_scales];
	double Sam_spread[num_scales];

	double maxSpread = 0;
		for(int iScale=0; iScale<num_scales; iScale++){
			avgFar[iScale]=0;
			Sam_avgFar[iScale]=0;
			Err_avgFar[iScale]=0;
			SamErr_avgFar[iScale]=0;
			average[iScale]=0;
			Sam_average[iScale]=0;
			spread[iScale]=0;
			Sam_spread[iScale]=0;
		}
	double temp_x_center = 0;
	double temp_y_center = 0;
	int nXbins = 0;
	int nYbins = 0;

	//no Ep cut and 800 mm
	double SamEfficiency_norm[8] = {0.693637956601778, 0.691846312216152, 0.694274060062294, 0.692659059055059, 0.689525923187573, 0.687340515173108, 0.67264957907054, 0.680786973870171};
		double SamError_norm[8] = {0.000591709314995, 0.000551299779889, 0.000553218015263, 0.000592877402991, 0.001338277458643, 0.001345896523897, 0.001345890830089, 0.001424837129587};
	double SamEfficiency_rate[8] = {0.704253886550915, 0.690380051863331, 0.691772813550098, 0.69256361135291, 0.709932624524099, 0.661602778381197, 0.689727048087611, 0.63319259891031};
		double SamError_rate[8] = {0.000590285353639, 0.000551435974537, 0.000553470663535, 0.000592887566979, 0.00133192058038, 0.001347645515483, 0.001343659008019, 0.001420343859735};

	//no Ep cut and 300 mm DT cut
/*		double SamEfficiency_norm[8] = {0.159968720088428, 0.158578519219699, 0.160896733935615, 0.159935578263097, 0.158580190397924, 0.157312720964749, 0.155371915857043, 0.15794665393424};
		double SamError_norm[8] = {0.000470532654049, 0.000436142386101, 0.000441210684922, 0.00047100310391, 0.001056547489991, 0.00105707315159, 0.001039029717735, 0.001114664890123};
		double SamEfficiency_rate[8] = {0,0,0,0,0,0,0,0};
		double SamError_rate[8] = {0,0,0,0,0,0,0,0};*/

	//no Ep cut and 500 mm DT cut
/*		double SamEfficiency_norm[8] = {0.406734952067528, 0.404598626937501, 0.407527885055175, 0.405870825030518, 0.404432216537415, 0.401825285866211, 0.394581950773671, 0.398460452867768};
		double SamError_norm[8] = {0.000630965238331, 0.00058590210777, 0.000590444977653, 0.000630924194384, 0.001421152629578, 0.001424377202279, 0.00140153245937, 0.001497519673178};
		double SamEfficiency_rate[8] = {0,0,0,0,0,0,0,0};
		double SamError_rate[8] = {0,0,0,0,0,0,0,0};*/

	//no Ep cut and 1000 mm DT cut
/*		double SamEfficiency_norm[8] = {0.811129573027829, 0.80929033977812, 0.81110078642876, 0.810979173149049, 0.808020380864248, 0.80585882531767, 0.787589324651952, 0.802100985321194};
		double SamError_norm[8] = {0.000502413795086, 0.000469072064629, 0.00047002219087, 0.000503103644241, 0.001139192136498, 0.001151102261117, 0.001172033285516, 0.00122011834718};
		double SamEfficiency_rate[8] = {0,0,0,0,0,0,0,0};
		double SamError_rate[8] = {0,0,0,0,0,0,0,0};*/

	//Ep cut of 3.5 MeV and 800 mm DT cut
/*		double SamEfficiency_norm[8] = {0.700335725872836, 0.699220738470729, 0.701489767580107, 0.700830813136867, 0.702664297257298, 0.704290468956811, 0.700494372740315, 0.697714902237776};
		double SamError_norm[8] = {0.000881043931247, 0.000820273136478, 0.00082260071475, 0.000880709552356, 0.00198450768137, 0.002001466974471, 0.001992591241917, 0.002121760610875};
		double SamEfficiency_rate[8] = {0,0,0,0,0,0,0,0};
		double SamError_rate[8] = {0,0,0,0,0,0,0,0};*/

	//Ep cut and 300 mm DT cut
/*		double SamEfficiency_norm[8] = {0.166810951222475, 0.166630959621351, 0.168261258174192, 0.167006431896991, 0.164867226665123, 0.166808670338921, 0.168820661760003, 0.16924080434732};
		double SamError_norm[8] = {0.000716986526508, 0.000666536479089, 0.000672487680147, 0.000717389647228, 0.001611017433542, 0.001635014252313, 0.001629572049498, 0.001732363933211};
		double SamEfficiency_rate[8] = {0,0,0,0,0,0,0,0};
		double SamError_rate[8] = {0,0,0,0,0,0,0,0};*/

//	cout << SamEfficiency_rate[3] << endl;

	//Correcting for relative error?

/*	cout << "Sam's normalized error before: " << SamError_norm[0] << endl;
	for(int iad = 0; iad < maxAD; iad++){
		SamError_norm[iad] = SamError_norm[iad]/SamEfficiency_norm[iad];
		SamError_rate[iad] = SamError_rate[iad]/SamEfficiency_rate[iad];
	}
	cout << "Sam's normalized error after: " << SamError_norm[0] << endl;*/

	//Make histograms/graphs
	TH1D* h_eff_DT_norm[maxAD];
	TH1D* h_eff_DT_rate[maxAD];
	TH1D* h_eff_DT_DTnorm[maxAD];

	TH1D* h_effDiff_DT_norm[maxAD];
	TH1D* h_effDiff_DT_rate[maxAD];
	TH1D* h_effDiff_DT_DTnorm[maxAD];

	TH1D* h_rel_effDiff_DT_norm[maxAD];
	TH1D* h_rel_effDiff_DT_rate[maxAD];
	TH1D* h_rel_effDiff_DT_DTnorm[maxAD];
	for(int iad=0; iad<maxAD; ++iad){
		sprintf(name,"h_eff_DT_norm_ad%d",iad+1);
		h_eff_DT_norm[iad]=new TH1D(name,name,500,0.01,10.01);

		sprintf(name,"h_eff_DT_rate_ad%d",iad+1);
		h_eff_DT_rate[iad]=new TH1D(name,name,500,0.01,10.01);

		sprintf(name,"h_eff_DT_DTnorm_ad%d",iad+1);
		h_eff_DT_DTnorm[iad]=new TH1D(name,name,500,0.01,10.01);

		sprintf(name,"h_effDiff_DT_norm_ad%d",iad+1);
		h_effDiff_DT_norm[iad]=new TH1D(name,name,500,0.01,10.01);

		sprintf(name,"h_effDiff_DT_rate_ad%d",iad+1);
		h_effDiff_DT_rate[iad]=new TH1D(name,name,500,0.01,10.01);

		sprintf(name,"h_effDiff_DT_DTnorm_ad%d",iad+1);
		h_effDiff_DT_DTnorm[iad]=new TH1D(name,name,500,0.01,10.01);

		sprintf(name,"h_rel_effDiff_DT_norm_ad%d",iad+1);
		h_rel_effDiff_DT_norm[iad]=new TH1D(name,name,500,0.01,10.01);

		sprintf(name,"h_rel_effDiff_DT_rate_ad%d",iad+1);
		h_rel_effDiff_DT_rate[iad]=new TH1D(name,name,500,0.01,10.01);

		sprintf(name,"h_rel_effDiff_DT_DTnorm_ad%d",iad+1);
		h_rel_effDiff_DT_DTnorm[iad]=new TH1D(name,name,500,0.01,10.01);
	}

	TH1D* h_effDiff_DT_norm_farAvg=new TH1D("h_effDiff_DT_norm_farAvg","h_effDiff_DT_norm_farAvg",500,0.01,10.01);
	TH1D* h_effDiff_DT_rate_farAvg=new TH1D("h_effDiff_DT_rate_farAvg","h_effDiff_DT_rate_farAvg",500,0.01,10.01);
	TH1D* h_effDiff_DT_DTnorm_farAvg=new TH1D("h_effDiff_DT_DTnorm_farAvg","h_effDiff_DT_DTnorm_farAvg",500,0.01,10.01);

	TH1D* h_effDiff_DT_norm_nearAvg=new TH1D("h_effDiff_DT_norm_nearAvg","h_effDiff_DT_norm_nearAvg",500,0.01,10.01);
	TH1D* h_effDiff_DT_rate_nearAvg=new TH1D("h_effDiff_DT_rate_nearAvg","h_effDiff_DT_rate_nearAvg",500,0.01,10.01);
	TH1D* h_effDiff_DT_DTnorm_nearAvg=new TH1D("h_effDiff_DT_DTnorm_nearAvg","h_effDiff_DT_DTnorm_nearAvg",500,0.01,10.01);

	TH1D* h_rel_effDiff_DT_norm_farAvg=new TH1D("h_rel_effDiff_DT_norm_farAvg","h_rel_effDiff_DT_norm_farAvg",500,0.01,10.01);
	TH1D* h_rel_effDiff_DT_rate_farAvg=new TH1D("h_rel_effDiff_DT_rate_farAvg","h_rel_effDiff_DT_rate_farAvg",500,0.01,10.01);
	TH1D* h_rel_effDiff_DT_DTnorm_farAvg=new TH1D("h_rel_effDiff_DT_DTnorm_farAvg","h_rel_effDiff_DT_DTnorm_farAvg",500,0.01,10.01);

	TH1D* h_allDT_norm_weighted=new TH1D("h_allDT_norm_weighted","h_allDT_norm_weighted",500,0.,10.);
	TH1D* h_allDT_rate_weighted=new TH1D("h_allDT_rate_weighted","h_allDT_rate_weighted",500,0.,10.);
	TH1D* h_allDT_DTnorm_weighted=new TH1D("h_allDT_DTnorm_weighted","h_allDT_DTnorm_weighted",500,0.,10.);

	TH1D* h_avgEff_norm_weighted=new TH1D("h_avgEff_norm_weighted","h_avgEff_norm_weighted",500,0.01,10.01);
	TH1D* h_avgEff_rate_weighted=new TH1D("h_avgEff_rate_weighted","h_avgEff_rate_weighted",500,0.01,10.01);
	TH1D* h_avgEff_DTnorm_weighted=new TH1D("h_avgEff_DTnorm_weighted","h_avgEff_DTnorm_weighted",500,0.01,10.01);

	TH2F* h_cut_distVStime[maxAD][num_scales];
	TH2F* h_cut_IBD_distVStime[maxAD];
	for(int iad=0; iad<maxAD; ++iad){
		sprintf(name,"h_cut_IBD_distVStime_ad%d",iad+1);
		h_cut_IBD_distVStime[iad]=new TH2F(name,name,1999,1,2000,700,0,7.);

		for(int iScale=0; iScale<num_scales; iScale++){
			if(iScale == 0) sprintf(name,"h_cut_distVStime_norm_ad%d",iad+1);
			if(iScale == 1) sprintf(name,"h_cut_distVStime_rate_ad%d",iad+1);
			h_cut_distVStime[iad][iScale]=new TH2F(name,name,1999,1,2000,700,0,7.);

			combinedEfficiency[iad][iScale] = 0;
			errDenom[iad][iScale] = 0;
			sig1[iad][iScale] = 0;
			sig2[iad][iScale] = 0;
			statError[iad][iScale] = 0;
		}
	}

	TH1F* h_ADeff[num_scales];
	TH1F* h_Sam_ADeff[num_scales];

	TH1F* h_ADeff_norm;
		h_ADeff_norm = new TH1F("h_ADeff_norm","h_ADeff_norm",8,0.5,8.5);
	TH1F* h_ADeff_rate;
		h_ADeff_rate = new TH1F("h_ADeff_rate","h_ADeff_rate",8,0.4,8.4);
	TH1F* h_ADeff_DTnorm;
		h_ADeff_DTnorm = new TH1F("h_ADeff_DTnorm","h_ADeff_DTnorm",8,0.6,8.6);

	TH1F* h_ADeff_farAvg[num_scales];
	TH1F* h_Sam_ADeff_farAvg[num_scales];
	for(int iScale=0; iScale<num_scales; iScale++){
		if(iScale == 0) sprintf(name,"h_ADeff_norm");
		if(iScale == 1) sprintf(name,"h_ADeff_rate");
		h_ADeff[iScale]=new TH1F(name,name,8,0.5,8.5);

		if(iScale == 0) sprintf(name,"h_Sam_ADeff_norm");
		if(iScale == 1) sprintf(name,"h_Sam_ADeff_rate");
		h_Sam_ADeff[iScale]=new TH1F(name,name,8,0.5,8.5);

		if(iScale == 0) sprintf(name,"h_ADeff_norm_farAvg");
		if(iScale == 1) sprintf(name,"h_ADeff_rate_farAvg");
		h_ADeff_farAvg[iScale]=new TH1F(name,name,8,0.5,8.5);

		if(iScale == 0) sprintf(name,"h_Sam_ADeff_norm_farAvg");
		if(iScale == 1) sprintf(name,"h_Sam_ADeff_rate_farAvg");
		h_Sam_ADeff_farAvg[iScale]=new TH1F(name,name,8,0.5,8.5);
	}

	TMultiGraph* ADs_compared;
		ADs_compared = new TMultiGraph("ADs_compared","ADs_compared");
	TGraph* ad_comp[num_scales];
	TGraph* Sam_ad_comp[num_scales];

	TMultiGraph* ADs_compared_farAvg;
		ADs_compared_farAvg = new TMultiGraph("ADs_compared_farAvg","ADs_compared_farAvg");
	TGraph* ad_comp_farAvg[num_scales];
	TGraph* Sam_ad_comp_farAvg[num_scales];

	for(int iScale=0; iScale < num_scales; iScale++){
		ad_comp[iScale] = new TGraph();
		if(iScale==0) ad_comp[iScale]->SetName("DT_norm");
		if(iScale==1) ad_comp[iScale]->SetName("DT_rate");

		Sam_ad_comp[iScale] = new TGraphErrors();
		if(iScale==0) Sam_ad_comp[iScale]->SetName("Sam_DT_norm");
		if(iScale==1) Sam_ad_comp[iScale]->SetName("Sam_DT_rate");

		ad_comp_farAvg[iScale] = new TGraph();
		if(iScale==0) ad_comp_farAvg[iScale]->SetName("DT_norm_farAvg");
		if(iScale==1) ad_comp_farAvg[iScale]->SetName("DT_rate_farAvg");

		Sam_ad_comp_farAvg[iScale] = new TGraph();
		if(iScale==0) Sam_ad_comp_farAvg[iScale]->SetName("Sam_DT_norm_farAvg");
		if(iScale==1) Sam_ad_comp_farAvg[iScale]->SetName("Sam_DT_rate_farAvg");
	}



	//Determining the bin number for doing the integral:
	for(int iX = 0; iX < 2000; iX++){
		temp_x_center = h_cut_distVStime[0][0]->GetXaxis()->GetBinCenter(iX);
		if(temp_x_center > max_x_value) break;
		nXbins = iX;
	}
	for(int iY = 0; iY < 700; iY++){
		temp_y_center = h_cut_distVStime[0][0]->GetYaxis()->GetBinCenter(iY);
		if(temp_y_center > max_y_value) break;
		nYbins = iY;
	}


	//Get files & histograms
	char inputname[64];
	sprintf(inputname,"./AdSimpleNL/3sig/distEff2000_results.root");
	TFile *distFile = new TFile(inputname);

	for(int iScale = 0; iScale < num_scales; iScale++){

		cout << endl << endl << "iScale: " << iScale << endl;
		cout << "AD\tDT Efficiency" << endl;
		for(int iad = 0; iad < maxAD; iad++){

			if(iScale == 0) sprintf(distTimeGraphName,"h_sub_distVStime_norm_ad%d",iad+1);
			else sprintf(distTimeGraphName,"h_sub_distVStime_rate_ad%d",iad+1);

			sprintf(IBDdistTimeGraphName,"h_IBD_distVStime_ad%d",iad+1);

			if(Ep_cut == 1){
				if(iScale == 0) sprintf(distTimeGraphName,"h_sub_distVStime_Ep35_norm_ad%d",iad+1);
				else sprintf(distTimeGraphName,"h_sub_distVStime_Ep35_rate_ad%d",iad+1);

				sprintf(IBDdistTimeGraphName,"h_IBD_distVStime_Ep35_ad%d",iad+1);
			}

	//		cout << distTimeGraphName << endl;
			TH2F* distTimeGraph = (TH2F*)distFile->Get(distTimeGraphName);
			TH2F* IBDdistTimeGraph = (TH2F*)distFile->Get(IBDdistTimeGraphName);

			//Determine if the bin is within the requirements of the DT cut and add it to histogram
			for(int iX = 0; iX < 2000; iX++){
				for(int iY = 0; iY < 700; iY++){
					time_center = distTimeGraph->GetXaxis()->GetBinCenter(iX);
					dist_center = distTimeGraph->GetYaxis()->GetBinCenter(iY);
					temp_value = distTimeGraph->GetBinContent(iX,iY);
					if(dist_center*1.e3 + 1.e3/600.*time_center > DTcut) continue;
					h_cut_distVStime[iad][iScale]->Fill(time_center,dist_center,temp_value);
					if(iScale == 0) h_cut_IBD_distVStime[iad]->Fill(time_center, dist_center, IBDdistTimeGraph->GetBinContent(iX,iY));
				}
			}

			//Calculating the efficiency
			combinedEfficiency[iad][iScale]=(h_cut_distVStime[iad][iScale]->Integral(0, nXbins, 0, nYbins))/(distTimeGraph->Integral(0, nXbins, 0, nYbins));

			//Calculating the statistical error
			errDenom[iad][iScale] = distTimeGraph->Integral(0, nXbins, 0, nYbins);
			sig1[iad][iScale] = (1/errDenom[iad][iScale] - (h_cut_distVStime[iad][iScale]->Integral(0,nXbins,0,nYbins))/(pow(errDenom[iad][iScale],2)))*sqrt(h_cut_IBD_distVStime[iad]->Integral(0,nXbins,0,nYbins));
			sig2[iad][iScale] = ((h_cut_distVStime[iad][iScale]->Integral(0,nXbins,0,nYbins))/(pow(errDenom[iad][iScale],2)))*sqrt((IBDdistTimeGraph->Integral(0,nXbins,0,nYbins))-(h_cut_IBD_distVStime[iad]->Integral(0,nXbins,0,nYbins)));
			statError[iad][iScale] = sqrt(pow(sig1[iad][iScale],2)+pow(sig2[iad][iScale],2));

//			statError[iad][iScale] = sqrt(combinedEfficiency[iad][iScale]*(1-combinedEfficiency[iad][iScale])/(distTimeGraph->Integral(0, nXbins, 0, nYbins)));

			cout << iad+1 << "\t" << combinedEfficiency[iad][iScale] << "\t+/-\t" << statError[iad][iScale] << endl;

			//Plotting the efficiencies
			if(iScale == 0){
				Sam_ad_comp[iScale]->SetPoint(iad, iad+1, SamEfficiency_norm[iad]*100);
				h_Sam_ADeff[iScale]->Fill(iad+1, SamEfficiency_norm[iad]*100);
				h_Sam_ADeff[iScale]->SetBinError(iad+1, SamError_norm[iad]*100);
				if(iad < 4){
					h_Sam_ADeff_farAvg[iScale]->Fill(iad+1, SamEfficiency_norm[iad]*100);
					h_Sam_ADeff_farAvg[iScale]->SetBinError(iad+1, SamError_norm[iad]*100);
				}
			}
			if(iScale == 1){
				Sam_ad_comp[iScale]->SetPoint(iad, iad+1, SamEfficiency_rate[iad]*100);
				h_Sam_ADeff[iScale]->Fill(iad+1, SamEfficiency_rate[iad]*100);
				h_Sam_ADeff[iScale]->SetBinError(iad+1, SamError_rate[iad]*100);
				if(iad < 4){
					h_Sam_ADeff_farAvg[iScale]->Fill(iad+1, SamEfficiency_rate[iad]*100);
					h_Sam_ADeff_farAvg[iScale]->SetBinError(iad+1, SamError_rate[iad]*100);
				}
			}
			ad_comp[iScale]->SetPoint(iad, iad+1, combinedEfficiency[iad][iScale]*100);
			h_ADeff[iScale]->Fill(iad+1, combinedEfficiency[iad][iScale]*100);
			h_ADeff[iScale]->SetBinError(iad+1, statError[iad][iScale]*100);
		}
	}

	//average far ADs:
	for(int iScale=0; iScale<num_scales; iScale++){
		for(int iad=0; iad<maxAD; iad++){
			average[iScale] = average[iScale]+combinedEfficiency[iad][iScale];
			if(iScale == 0) Sam_average[iScale] = Sam_average[iScale]+SamEfficiency_norm[iad];
			if(iScale == 1) Sam_average[iScale] = Sam_average[iScale]+SamEfficiency_rate[iad];


			if(iad < 4){
				if(iScale == 0) Sam_ad_comp_farAvg[iScale]->SetPoint(iad, iad+1, SamEfficiency_norm[iad]*100);
				if(iScale == 1) Sam_ad_comp_farAvg[iScale]->SetPoint(iad, iad+1, SamEfficiency_rate[iad]*100);
				ad_comp_farAvg[iScale]->SetPoint(iad, iad+1, combinedEfficiency[iad][iScale]*100);
				h_ADeff_farAvg[iScale]->Fill(iad+1, combinedEfficiency[iad][iScale]*100);
				h_ADeff_farAvg[iScale]->SetBinError(iad+1, statError[iad][iScale]*100);
			}
			else{//calculate the average
				avgFar[iScale] = avgFar[iScale] + combinedEfficiency[iad][iScale];
				Err_avgFar[iScale] = sqrt(pow(Err_avgFar[iScale],2) + pow(statError[iad][iScale],2));
				if(iScale == 0){
					Sam_avgFar[iScale] = Sam_avgFar[iScale] + SamEfficiency_norm[iad];
					SamErr_avgFar[iScale] = sqrt(pow(SamErr_avgFar[iScale],2) + pow(SamError_norm[iad],2));
				}
				if(iScale == 1){
					Sam_avgFar[iScale] = Sam_avgFar[iScale] + SamEfficiency_rate[iad];
					SamErr_avgFar[iScale] = sqrt(pow(SamErr_avgFar[iScale],2) + pow(SamError_rate[iad],2));
				}
			}
		}
		avgFar[iScale] = avgFar[iScale]/4.;
		Err_avgFar[iScale] = Err_avgFar[iScale]/4.;
		Sam_avgFar[iScale] = Sam_avgFar[iScale]/4.;
		SamErr_avgFar[iScale] = SamErr_avgFar[iScale]/4.;

		for(int iad=4; iad<maxAD; iad++){
			Sam_ad_comp_farAvg[iScale]->SetPoint(iad, iad+1, Sam_avgFar[iScale]*100);
			ad_comp_farAvg[iScale]->SetPoint(iad, iad+1, avgFar[iScale]*100);

			h_Sam_ADeff_farAvg[iScale]->Fill(iad+1, Sam_avgFar[iScale]*100);
			h_Sam_ADeff_farAvg[iScale]->SetBinError(iad+1, SamErr_avgFar[iScale]*100);
			h_ADeff_farAvg[iScale]->Fill(iad+1, avgFar[iScale]*100);
			h_ADeff_farAvg[iScale]->SetBinError(iad+1, Err_avgFar[iScale]*100);
		}

		//Total average
		average[iScale] = average[iScale]/8.;
		Sam_average[iScale] = Sam_average[iScale]/8.;

		//Finding the spread
		spread[iScale] = avgFar[iScale] - average[iScale];
		Sam_spread[iScale] = Sam_avgFar[iScale] - Sam_average[iScale];

		if(spread[iScale] < 0) spread[iScale] = -1 * spread[iScale];
		if(Sam_spread[iScale] < 0) Sam_spread[iScale] = -1 * Sam_spread[iScale];

//		cout << "current spread: " << spread[iScale] << endl;

		for(int iad = 0; iad < 4; iad++){
				if(combinedEfficiency[iad][iScale]-average[iScale] >  0 && combinedEfficiency[iad][iScale]-average[iScale] > spread[iScale]) spread[iScale] = combinedEfficiency[iad][iScale]-average[iScale];
				if(combinedEfficiency[iad][iScale]-average[iScale] <  0 && -1*(combinedEfficiency[iad][iScale]-average[iScale]) > spread[iScale]) spread[iScale] = -1.* (combinedEfficiency[iad][iScale]-average[iScale]);

//		cout << "current difference: " << combinedEfficiency[iad][iScale]-average[iScale] << "\tcurrent spread: " << spread[iScale] << endl;

				if(iScale == 0){
					if(SamEfficiency_norm[iad] - Sam_average[iScale] >  0 && SamEfficiency_norm[iad] - Sam_average[iScale] > Sam_spread[iScale]) Sam_spread[iScale] = SamEfficiency_norm[iad] - Sam_average[iScale];
					if(SamEfficiency_norm[iad] - Sam_average[iScale] <  0 && -1*(SamEfficiency_norm[iad] - Sam_average[iScale]) > Sam_spread[iScale]) Sam_spread[iScale] = -1.* (SamEfficiency_norm[iad] - Sam_average[iScale]);
				}
				if(iScale == 1){
					if(SamEfficiency_rate[iad] - Sam_average[iScale] >  0 && SamEfficiency_rate[iad] - Sam_average[iScale] > Sam_spread[iScale]) Sam_spread[iScale] = SamEfficiency_rate[iad] - Sam_average[iScale];
					if(SamEfficiency_rate[iad] - Sam_average[iScale] <  0 && -1*(SamEfficiency_rate[iad] - Sam_average[iScale]) > Sam_spread[iScale]) Sam_spread[iScale] = -1.* (SamEfficiency_rate[iad] - Sam_average[iScale]);
				}
		}

	cout << endl << "For Scale: " << iScale << endl << "My average:\t" << 100*average[iScale] << "\t+/-\t" << 100*spread[iScale] << endl << "Sam's average:\t" << 100*Sam_average[iScale] << "\t+/-\t" << 100*Sam_spread[iScale] << endl << endl;


	}

	//Finding the max spread:
	maxSpread = spread[0];
	if(maxSpread < spread[1]) maxSpread = spread[1];
	if(maxSpread < Sam_spread[0]) maxSpread = Sam_spread[0];
	if(maxSpread < Sam_spread[1]) maxSpread = Sam_spread[1];



	//********************Using the other DT plots to calculate:************************
	int maxBin = 0;
	float temp_binCenter = 0;
	float DTvalue = .8; // in mm
	int DTbin = 0;
	double BedaError_norm[8];
	double sigma1_norm[8];
	double sigma2_norm[8];
	double sigmaK_norm[8];
	double denom_norm[8];

	double BedaError_rate[8];
	double sigma1_rate[8];
	double sigma2_rate[8];
	double sigmaK_rate[8];
	double denom_rate[8];

	double BedaError_DTnorm[8];
	double sigma1_DTnorm[8];
	double sigma2_DTnorm[8];
	double sigmaK_DTnorm[8];
	double denom_DTnorm[8];
		for(int i=0; i<8; i++){
			BedaError_norm[i]=0;
			sigma1_norm[i]=0;
			sigma2_norm[i]=0;
			sigmaK_norm[i]=0;
			denom_norm[i]=0;

			BedaError_rate[i]=0;
			sigma1_rate[i]=0;
			sigma2_rate[i]=0;
			sigmaK_rate[i]=0;
			denom_rate[i]=0;

			BedaError_DTnorm[i]=0;
			sigma1_DTnorm[i]=0;
			sigma2_DTnorm[i]=0;
			sigmaK_DTnorm[i]=0;
			denom_DTnorm[i]=0;
		}

	for(int iad = 0; iad<maxAD; iad++){

		sprintf(name,"h_sub_DT_norm_ad%d",iad+1);
			if(Ep_cut == 1) sprintf(name,"h_sub_DT_Ep35_norm_ad%d",iad+1);
		TH1D* DTsub_norm = (TH1D*)distFile->Get(name);

		sprintf(name,"h_sub_DT_rate_ad%d",iad+1);
			if(Ep_cut == 1) sprintf(name,"h_sub_DT_Ep35_rate_ad%d",iad+1);
			if(rateMod == 1) sprintf(name,"h_sub_DT_rate_mod_ad%d",iad+1);
			if(Ep_cut == 1 && rateMod == 1) sprintf(name,"h_sub_DT_Ep35_rate_mod_ad%d",iad+1);
		TH1D* DTsub_rate = (TH1D*)distFile->Get(name); //HEREHEREHEREHEREHERE

		sprintf(name,"h_sub_DT_DTnorm_ad%d",iad+1);
			if(Ep_cut == 1) sprintf(name,"h_sub_DT_Ep35_DTnorm_ad%d",iad+1);
		TH1D* DTsub_DTnorm = (TH1D*)distFile->Get(name);

		sprintf(name,"h_ibd_DT_ad%d",iad+1);
			if(Ep_cut == 1) sprintf(name,"h_ibd_DT_Ep35_ad%d",iad+1);
		TH1D* DTibd = (TH1D*)distFile->Get(name);

		sprintf(name,"h_acc_DT_norm_ad%d",iad+1);
			if(Ep_cut == 1) sprintf(name,"h_acc_DT_Ep35_norm_ad%d",iad+1);
		TH1D* DTacc_norm = (TH1D*)distFile->Get(name);

		sprintf(name,"h_acc_DT_rate_ad%d",iad+1);
			if(Ep_cut == 1) sprintf(name,"h_acc_DT_Ep35_rate_ad%d",iad+1);
		TH1D* DTacc_rate = (TH1D*)distFile->Get(name);

		sprintf(name,"h_acc_DT_DTnorm_ad%d",iad+1);
			if(Ep_cut == 1) sprintf(name,"h_acc_DT_Ep35_DTnorm_ad%d",iad+1);
		TH1D* DTacc_DTnorm = (TH1D*)distFile->Get(name);

		h_allDT_norm_weighted->Add(DTsub_norm);
		h_allDT_rate_weighted->Add(DTsub_rate);
		h_allDT_DTnorm_weighted->Add(DTsub_DTnorm);


		for(int iBin = 0; iBin < 502; iBin++){
			temp_binCenter = h_eff_DT_norm[iad]->GetBinCenter(iBin);
			if(temp_binCenter > 3.) break; //here's the place to set 100 percent point
			maxBin = iBin;
		}

		for(int iBin = 0; iBin < 502; iBin++){
			temp_binCenter = h_eff_DT_norm[iad]->GetBinCenter(iBin);
			if(temp_binCenter > DTvalue) break;
			DTbin = iBin;
		}
	//	cout << "DTbin is: " << DTbin << endl;


		for(int iBin = 1; iBin < 502; iBin++){
			h_avgEff_norm_weighted->SetBinContent(iBin,(h_allDT_norm_weighted->Integral(0,iBin))/(h_allDT_norm_weighted->Integral(0,maxBin)));
			h_avgEff_rate_weighted->SetBinContent(iBin,(h_allDT_rate_weighted->Integral(0,iBin))/(h_allDT_rate_weighted->Integral(0,maxBin)));
			h_avgEff_DTnorm_weighted->SetBinContent(iBin,(h_allDT_DTnorm_weighted->Integral(0,iBin))/(h_allDT_DTnorm_weighted->Integral(0,maxBin)));

			h_eff_DT_norm[iad]->SetBinContent(iBin,(DTsub_norm->Integral(0,iBin))/(DTsub_norm->Integral(0,maxBin)));
			h_eff_DT_rate[iad]->SetBinContent(iBin,(DTsub_rate->Integral(0,iBin))/(DTsub_rate->Integral(0,maxBin)));
			h_eff_DT_DTnorm[iad]->SetBinContent(iBin,(DTsub_DTnorm->Integral(0,iBin))/(DTsub_DTnorm->Integral(0,maxBin)));

			denom_norm[iad] = pow(DTsub_norm->Integral(0,maxBin),2);
			sigma1_norm[iad] = (DTsub_norm->Integral(iBin,maxBin))*sqrt(DTibd->Integral(0,iBin))/denom_norm[iad];
			sigma2_norm[iad] = (DTsub_norm->Integral(0,iBin))*sqrt(DTibd->Integral(iBin,maxBin))/denom_norm[iad];
			sigmaK_norm[iad] = ((DTibd->Integral(0,iBin)*DTacc_norm->Integral(iBin,maxBin))-(DTibd->Integral(iBin,maxBin)*DTacc_norm->Integral(0,iBin)))/(denom_norm[iad]*sqrt(DTibd->Integral(maxBin,500)));
			BedaError_norm[iad] = sqrt(pow(sigma1_norm[iad],2)+pow(sigma2_norm[iad],2)+pow(sigmaK_norm[iad],2));

				h_eff_DT_norm[iad]->SetBinError(iBin,BedaError_norm[iad]);

			denom_rate[iad] = pow(DTsub_rate->Integral(0,maxBin),2);
			sigma1_rate[iad] = (DTsub_rate->Integral(iBin,maxBin))*sqrt(DTibd->Integral(0,iBin))/denom_rate[iad];
			sigma2_rate[iad] = (DTsub_rate->Integral(0,iBin))*sqrt(DTibd->Integral(iBin,maxBin))/denom_rate[iad];
			//sigmaK_rate[iad] = ((DTibd->Integral(0,iBin)*DTacc_rate->Integral(iBin,maxBin))-(DTibd->Integral(iBin,maxBin)*DTacc_rate->Integral(0,iBin)))/(denom_rate[iad]*sqrt(DTibd->Integral(maxBin,500)));
			BedaError_rate[iad] = sqrt(pow(sigma1_rate[iad],2)+pow(sigma2_rate[iad],2)+pow(sigmaK_rate[iad],2));

				h_eff_DT_rate[iad]->SetBinError(iBin,BedaError_rate[iad]);

			denom_DTnorm[iad] = pow(DTsub_DTnorm->Integral(0,DTsub_DTnorm->FindBin(3)),2);
			sigma1_DTnorm[iad] = (DTsub_DTnorm->Integral(iBin,DTsub_DTnorm->FindBin(3)))*sqrt(DTibd->Integral(0,iBin))/denom_DTnorm[iad];
			sigma2_DTnorm[iad] = (DTsub_DTnorm->Integral(0,iBin))*sqrt(DTibd->Integral(iBin,DTsub_DTnorm->FindBin(3)))/denom_DTnorm[iad];
			sigmaK_DTnorm[iad] = ((DTibd->Integral(0,iBin)*DTacc_DTnorm->Integral(iBin,DTsub_DTnorm->FindBin(3)))-(DTibd->Integral(iBin,DTsub_DTnorm->FindBin(3))*DTacc_DTnorm->Integral(0,iBin)))/(denom_DTnorm[iad]*sqrt(DTibd->Integral(DTsub_DTnorm->FindBin(3),500)));
			BedaError_DTnorm[iad] = sqrt(pow(sigma1_DTnorm[iad],2)+pow(sigma2_DTnorm[iad],2)+pow(sigmaK_DTnorm[iad],2));

				h_eff_DT_DTnorm[iad]->SetBinError(iBin,BedaError_DTnorm[iad]);


			if(iBin == DTbin){
				h_ADeff_norm->SetBinContent(iad+1,(DTsub_norm->Integral(0,iBin))/(DTsub_norm->Integral(0,maxBin)));
				h_ADeff_rate->SetBinContent(iad+1,(DTsub_rate->Integral(0,iBin))/(DTsub_rate->Integral(0,maxBin)));
				h_ADeff_DTnorm->SetBinContent(iad+1,(DTsub_DTnorm->Integral(0,iBin))/(DTsub_DTnorm->Integral(0,maxBin)));

				h_ADeff_norm->SetBinError(iad+1,BedaError_norm[iad]);
				h_ADeff_rate->SetBinError(iad+1,BedaError_rate[iad]);
				h_ADeff_DTnorm->SetBinError(iad+1,BedaError_DTnorm[iad]);
			}
		}

/*			errDenom[iad][iScale] = distTimeGraph->Integral(0, nXbins, 0, nYbins);
			sig1[iad][iScale] = (1/errDenom[iad][iScale] - (h_cut_distVStime[iad][iScale]->Integral(0,nXbins,0,nYbins))/(pow(errDenom[iad][iScale],2)))*sqrt(h_cut_IBD_distVStime[iad]->Integral(0,nXbins,0,nYbins));
			sig2[iad][iScale] = ((h_cut_distVStime[iad][iScale]->Integral(0,nXbins,0,nYbins))/(pow(errDenom[iad][iScale],2)))*sqrt((IBDdistTimeGraph->Integral(0,nXbins,0,nYbins))-(h_cut_IBD_distVStime[iad]->Integral(0,nXbins,0,nYbins)));
			statError[iad][iScale] = sqrt(pow(sig1[iad][iScale],2)+pow(sig2[iad][iScale],2));

//			statError[iad][iScale] = sqrt(combinedEfficiency[iad][iScale]*(1-combinedEfficiency[iad][iScale])/(distTimeGraph->Integral(0, nXbins, 0, nYbins)));

			cout << iad+1 << "\t" << combinedEfficiency[iad][iScale] << "\t+/-\t" << statError[iad][iScale] << endl;*/

	}

		//Averaging them:
		double avgDT_norm = 0;
		double avgDT_rate = 0;
		double avgDT_DTnorm = 0;

		double avgDT_norm_farAvg = 0;
		double avgDT_rate_farAvg = 0;
		double avgDT_DTnorm_farAvg = 0;

		double avgDT_norm_nearAvg = 0;
		double avgDT_rate_nearAvg = 0;
		double avgDT_DTnorm_nearAvg = 0;

		double avgDT_err_norm_farAvg = 0;
		double avgDT_err_rate_farAvg = 0;
		double avgDT_err_DTnorm_farAvg = 0;
		for(int iBin = 0; iBin <502; iBin++){
			avgDT_norm = 0;
			avgDT_rate = 0;
			avgDT_DTnorm = 0;
			avgDT_norm_farAvg = 0;
			avgDT_rate_farAvg = 0;
			avgDT_DTnorm_farAvg = 0;
			avgDT_norm_nearAvg = 0;
			avgDT_rate_nearAvg = 0;
			avgDT_DTnorm_nearAvg = 0;
			avgDT_err_norm_farAvg = 0;
			avgDT_err_rate_farAvg = 0;
			avgDT_err_DTnorm_farAvg = 0;
			for(int iad = 0; iad<maxAD; iad++){
				avgDT_norm = avgDT_norm + h_eff_DT_norm[iad]->GetBinContent(iBin);
				avgDT_rate = avgDT_rate + h_eff_DT_rate[iad]->GetBinContent(iBin);
				avgDT_DTnorm = avgDT_DTnorm + h_eff_DT_DTnorm[iad]->GetBinContent(iBin);
				if(iad<4){
					avgDT_norm_nearAvg = avgDT_norm_nearAvg + h_eff_DT_norm[iad]->GetBinContent(iBin);
					avgDT_rate_nearAvg = avgDT_rate_nearAvg + h_eff_DT_rate[iad]->GetBinContent(iBin);
					avgDT_DTnorm_nearAvg = avgDT_DTnorm_nearAvg + h_eff_DT_DTnorm[iad]->GetBinContent(iBin);
				}

				if(iad>3){
					avgDT_norm_farAvg = avgDT_norm_farAvg + h_eff_DT_norm[iad]->GetBinContent(iBin);
					avgDT_rate_farAvg = avgDT_rate_farAvg + h_eff_DT_rate[iad]->GetBinContent(iBin);
					avgDT_DTnorm_farAvg = avgDT_DTnorm_farAvg + h_eff_DT_DTnorm[iad]->GetBinContent(iBin);

					avgDT_err_norm_farAvg = sqrt(pow(avgDT_err_norm_farAvg,2) + pow(h_eff_DT_norm[iad]->GetBinError(iBin),2));
					avgDT_err_rate_farAvg = sqrt(pow(avgDT_err_rate_farAvg,2) + pow(h_eff_DT_rate[iad]->GetBinError(iBin),2));
					avgDT_err_DTnorm_farAvg = sqrt(pow(avgDT_err_DTnorm_farAvg,2) + pow(h_eff_DT_DTnorm[iad]->GetBinError(iBin),2));
				}
			}
			avgDT_norm = avgDT_norm/8.;
			avgDT_rate = avgDT_rate/8.;
			avgDT_DTnorm = avgDT_DTnorm/8.;
				if(weighted == 1){
					avgDT_norm = h_avgEff_norm_weighted->GetBinContent(iBin);
					avgDT_rate = h_avgEff_rate_weighted->GetBinContent(iBin);
					avgDT_DTnorm = h_avgEff_DTnorm_weighted->GetBinContent(iBin);
				}

			avgDT_norm_nearAvg = avgDT_norm_nearAvg/4.;
			avgDT_rate_nearAvg = avgDT_rate_nearAvg/4.;
			avgDT_DTnorm_nearAvg = avgDT_DTnorm_nearAvg/4.;

			avgDT_norm_farAvg = avgDT_norm_farAvg/4.;
			avgDT_rate_farAvg = avgDT_rate_farAvg/4.;
			avgDT_DTnorm_farAvg = avgDT_DTnorm_farAvg/4.;

			avgDT_err_norm_farAvg = avgDT_err_norm_farAvg/4.;
			avgDT_err_rate_farAvg = avgDT_err_rate_farAvg/4.;
			avgDT_err_DTnorm_farAvg = avgDT_err_DTnorm_farAvg/4.;

			h_effDiff_DT_norm_farAvg->SetBinContent(iBin, avgDT_norm_farAvg-avgDT_norm);
			h_effDiff_DT_rate_farAvg->SetBinContent(iBin, avgDT_rate_farAvg-avgDT_rate);
			h_effDiff_DT_DTnorm_farAvg->SetBinContent(iBin, avgDT_DTnorm_farAvg-avgDT_DTnorm);

			h_effDiff_DT_norm_nearAvg->SetBinContent(iBin, avgDT_norm_nearAvg-avgDT_norm);
			h_effDiff_DT_rate_nearAvg->SetBinContent(iBin, avgDT_rate_nearAvg-avgDT_rate);
			h_effDiff_DT_DTnorm_nearAvg->SetBinContent(iBin, avgDT_DTnorm_nearAvg-avgDT_DTnorm);

			h_rel_effDiff_DT_norm_farAvg->SetBinContent(iBin, (avgDT_norm_farAvg-avgDT_norm)/(2*avgDT_norm));
			h_rel_effDiff_DT_rate_farAvg->SetBinContent(iBin, (avgDT_rate_farAvg-avgDT_rate)/(2*avgDT_rate));
			h_rel_effDiff_DT_DTnorm_farAvg->SetBinContent(iBin, (avgDT_DTnorm_farAvg-avgDT_DTnorm)/(2*avgDT_DTnorm));

			h_effDiff_DT_norm_farAvg->SetBinError(iBin, avgDT_err_norm_farAvg);
			h_effDiff_DT_rate_farAvg->SetBinError(iBin, avgDT_err_rate_farAvg);
			h_effDiff_DT_DTnorm_farAvg->SetBinError(iBin, avgDT_err_DTnorm_farAvg);

			h_rel_effDiff_DT_norm_farAvg->SetBinError(iBin, avgDT_err_norm_farAvg/(2.*avgDT_norm));
			h_rel_effDiff_DT_rate_farAvg->SetBinError(iBin, avgDT_err_rate_farAvg/(2.*avgDT_rate));
			h_rel_effDiff_DT_DTnorm_farAvg->SetBinError(iBin, avgDT_err_DTnorm_farAvg/(2.*avgDT_DTnorm));

			for(int iad = 0; iad<maxAD; iad++){
				h_effDiff_DT_norm[iad]->SetBinContent(iBin, (h_eff_DT_norm[iad]->GetBinContent(iBin))-avgDT_norm);
				h_effDiff_DT_rate[iad]->SetBinContent(iBin, (h_eff_DT_rate[iad]->GetBinContent(iBin))-avgDT_rate);
				h_effDiff_DT_DTnorm[iad]->SetBinContent(iBin, (h_eff_DT_DTnorm[iad]->GetBinContent(iBin))-avgDT_DTnorm);

				h_effDiff_DT_norm[iad]->SetBinError(iBin, (h_eff_DT_norm[iad]->GetBinError(iBin)));
				h_effDiff_DT_rate[iad]->SetBinError(iBin, (h_eff_DT_rate[iad]->GetBinError(iBin)));
				h_effDiff_DT_DTnorm[iad]->SetBinError(iBin, (h_eff_DT_DTnorm[iad]->GetBinError(iBin)));

				h_rel_effDiff_DT_norm[iad]->SetBinContent(iBin, ((h_eff_DT_norm[iad]->GetBinContent(iBin))-avgDT_norm)/(2*avgDT_norm));
				h_rel_effDiff_DT_rate[iad]->SetBinContent(iBin, ((h_eff_DT_rate[iad]->GetBinContent(iBin))-avgDT_rate)/(2*avgDT_rate));
				h_rel_effDiff_DT_DTnorm[iad]->SetBinContent(iBin, ((h_eff_DT_DTnorm[iad]->GetBinContent(iBin))-avgDT_DTnorm)/(2*avgDT_DTnorm));

				h_rel_effDiff_DT_norm[iad]->SetBinError(iBin, 1/2.*(h_eff_DT_norm[iad]->GetBinError(iBin))); //these are still wrong currently
				h_rel_effDiff_DT_rate[iad]->SetBinError(iBin, 1/2.*(h_eff_DT_rate[iad]->GetBinError(iBin)));
				h_rel_effDiff_DT_DTnorm[iad]->SetBinError(iBin, 1/2.*(h_eff_DT_DTnorm[iad]->GetBinError(iBin)));
			}

		}


//	int colors[4] = {kBlue, kRed, kAzure+5, kMagenta}; //my normalized, my rate, Sam's normalized, Sam's rate

	int colors[8] = {kBlack, kRed, kBlue, kMagenta, kOrange, kSpring-8, kTeal, kGray};

	gStyle->SetPalette(kTemperatureMap);

	//Save histogram
        char outputname[64];
        char title[80];
        sprintf(outputname,"../nH_files/DT_results.root");
	TFile* outfile=new TFile(outputname, "RECREATE");

	for(int iad = 0; iad < maxAD; iad++){
		for(int iScale = 0; iScale < num_scales; iScale++){
			h_cut_distVStime[iad][iScale]->SetStats(0);
			h_cut_distVStime[iad][iScale]->SetOption("COLZ");
			h_cut_distVStime[iad][iScale]->GetXaxis()->SetTitle("Delta Time [us]");
			h_cut_distVStime[iad][iScale]->GetYaxis()->SetTitle("Distance [m]");
			h_cut_distVStime[iad][iScale]->Write();
		}

	//		h_eff_DT_norm[iad]->SetStats(0);
			h_eff_DT_norm[iad]->GetXaxis()->SetTitle("DT [m]");
			h_eff_DT_norm[iad]->GetYaxis()->SetTitle("Efficiency");
			h_eff_DT_norm[iad]->Write();

	//		h_eff_DT_rate[iad]->SetStats(0);
			h_eff_DT_rate[iad]->GetXaxis()->SetTitle("DT [m]");
			h_eff_DT_rate[iad]->GetYaxis()->SetTitle("Efficiency");
			h_eff_DT_rate[iad]->Write();

			h_eff_DT_DTnorm[iad]->SetStats(0);
			h_eff_DT_DTnorm[iad]->GetXaxis()->SetTitle("DT [m]");
			h_eff_DT_DTnorm[iad]->GetYaxis()->SetTitle("Efficiency");
			h_eff_DT_DTnorm[iad]->Write();

	//		h_effDiff_DT_norm[iad]->SetStats(0);
			h_effDiff_DT_norm[iad]->GetXaxis()->SetTitle("DT [m]");
			h_effDiff_DT_norm[iad]->GetYaxis()->SetTitle("Absolute Difference in Efficiency");
			h_effDiff_DT_norm[iad]->Write();

	//		h_effDiff_DT_rate[iad]->SetStats(0);
			h_effDiff_DT_rate[iad]->GetXaxis()->SetTitle("DT [m]");
			h_effDiff_DT_rate[iad]->GetYaxis()->SetTitle("Absolute Difference in Efficiency");
			h_effDiff_DT_rate[iad]->Write();

	//		h_effDiff_DT_DTnorm[iad]->SetStats(0);
			h_effDiff_DT_DTnorm[iad]->GetXaxis()->SetTitle("DT [m]");
			h_effDiff_DT_DTnorm[iad]->GetYaxis()->SetTitle("Absolute Difference in Efficiency");
			h_effDiff_DT_DTnorm[iad]->Write();

	//		h_rel_effDiff_DT_norm[iad]->SetStats(0);
			h_rel_effDiff_DT_norm[iad]->GetXaxis()->SetTitle("DT [m]");
			h_rel_effDiff_DT_norm[iad]->GetYaxis()->SetTitle("Relative Difference in Efficiency");
			h_rel_effDiff_DT_norm[iad]->Write();

	//		h_rel_effDiff_DT_rate[iad]->SetStats(0);
			h_rel_effDiff_DT_rate[iad]->GetXaxis()->SetTitle("DT [m]");
			h_rel_effDiff_DT_rate[iad]->GetYaxis()->SetTitle("Relative Difference in Efficiency");
			h_rel_effDiff_DT_rate[iad]->Write();

	//		h_rel_effDiff_DT_DTnorm[iad]->SetStats(0);
			h_rel_effDiff_DT_DTnorm[iad]->GetXaxis()->SetTitle("DT [m]");
			h_rel_effDiff_DT_DTnorm[iad]->GetYaxis()->SetTitle("Relative Difference in Efficiency");
			h_rel_effDiff_DT_DTnorm[iad]->Write();

	}

	//		h_effDiff_DT_norm_farAvg->SetStats(0);
			h_effDiff_DT_norm_farAvg->GetXaxis()->SetTitle("DT [m]");
			h_effDiff_DT_norm_farAvg->GetYaxis()->SetTitle("Absolute Difference in Efficiency");
			h_effDiff_DT_norm_farAvg->Write();

	//		h_effDiff_DT_rate_farAvg->SetStats(0);
			h_effDiff_DT_rate_farAvg->GetXaxis()->SetTitle("DT [m]");
			h_effDiff_DT_rate_farAvg->GetYaxis()->SetTitle("Absolute Difference in Efficiency");
			h_effDiff_DT_rate_farAvg->Write();

	//		h_effDiff_DT_DTnorm_farAvg->SetStats(0);
			h_effDiff_DT_DTnorm_farAvg->GetXaxis()->SetTitle("DT [m]");
			h_effDiff_DT_DTnorm_farAvg->GetYaxis()->SetTitle("Absolute Difference in Efficiency");
			h_effDiff_DT_DTnorm_farAvg->Write();


	//		h_rel_effDiff_DT_norm_farAvg->SetStats(0);
			h_rel_effDiff_DT_norm_farAvg->GetXaxis()->SetTitle("DT [m]");
			h_rel_effDiff_DT_norm_farAvg->GetYaxis()->SetTitle("Relative Difference in Efficiency");
			h_rel_effDiff_DT_norm_farAvg->Write();

	//		h_rel_effDiff_DT_rate_farAvg->SetStats(0);
			h_rel_effDiff_DT_rate_farAvg->GetXaxis()->SetTitle("DT [m]");
			h_rel_effDiff_DT_rate_farAvg->GetYaxis()->SetTitle("Relative Difference in Efficiency");
			h_rel_effDiff_DT_rate_farAvg->Write();

	//		h_rel_effDiff_DT_DTnorm_farAvg->SetStats(0);
			h_rel_effDiff_DT_DTnorm_farAvg->GetXaxis()->SetTitle("DT [m]");
			h_rel_effDiff_DT_DTnorm_farAvg->GetYaxis()->SetTitle("Relative Difference in Efficiency");
			h_rel_effDiff_DT_DTnorm_farAvg->Write();


	ADs_compared->GetXaxis()->SetTitle("AD Number");
	ADs_compared->GetYaxis()->SetTitle("DT Efficiency");

	for(int iScale = 0; iScale < num_scales; iScale++){
		Sam_ad_comp[iScale]->GetXaxis()->SetTitle("AD Number");
		Sam_ad_comp[iScale]->GetYaxis()->SetTitle("DT Efficiency");
		Sam_ad_comp[iScale]->SetMarkerColor(colors[0]);
		Sam_ad_comp[iScale]->SetLineColor(colors[0]);
		if(iScale == 1) Sam_ad_comp[iScale]->SetLineStyle(7);
		Sam_ad_comp[iScale]->SetMarkerStyle(20);
		if(iScale == 1) Sam_ad_comp[iScale]->SetMarkerStyle(33);

		ad_comp[iScale]->GetXaxis()->SetTitle("AD Number");
		ad_comp[iScale]->GetYaxis()->SetTitle("DT Efficiency");
		ad_comp[iScale]->SetMarkerColor(colors[1]);
		ad_comp[iScale]->SetLineColor(colors[1]);
		if(iScale == 1) ad_comp[iScale]->SetLineStyle(7);
		ad_comp[iScale]->SetMarkerStyle(20);
		if(iScale == 1) ad_comp[iScale]->SetMarkerStyle(33);

		ADs_compared->Add(ad_comp[iScale]);
		if(Ep_cut != 1 || iScale == 0) ADs_compared->Add(Sam_ad_comp[iScale]);
	}

	ADs_compared->GetXaxis()->SetLimits(0.25,8.25);
	ADs_compared->GetYaxis()->SetRangeUser(62,72);
	ADs_compared->Write();

	ADs_compared_farAvg->GetXaxis()->SetTitle("AD Number");
	ADs_compared_farAvg->GetYaxis()->SetTitle("DT Efficiency");

	for(int iScale = 0; iScale < num_scales; iScale++){
		Sam_ad_comp_farAvg[iScale]->GetXaxis()->SetTitle("AD Number");
		Sam_ad_comp_farAvg[iScale]->GetYaxis()->SetTitle("DT Efficiency");
		Sam_ad_comp_farAvg[iScale]->SetMarkerColor(colors[0]);
		Sam_ad_comp_farAvg[iScale]->SetLineColor(colors[0]);
		if(iScale == 1) Sam_ad_comp_farAvg[iScale]->SetLineStyle(7);
		Sam_ad_comp_farAvg[iScale]->SetMarkerStyle(20);
		if(iScale == 1) Sam_ad_comp_farAvg[iScale]->SetMarkerStyle(33);

		ad_comp_farAvg[iScale]->GetXaxis()->SetTitle("AD Number");
		ad_comp_farAvg[iScale]->GetYaxis()->SetTitle("DT Efficiency");
		ad_comp_farAvg[iScale]->SetMarkerColor(colors[1]);
		ad_comp_farAvg[iScale]->SetLineColor(colors[1]);
		if(iScale == 1) ad_comp_farAvg[iScale]->SetLineStyle(7);
		ad_comp_farAvg[iScale]->SetMarkerStyle(20);
		if(iScale == 1) ad_comp_farAvg[iScale]->SetMarkerStyle(33);

		ADs_compared_farAvg->Add(ad_comp_farAvg[iScale]);
		if(Ep_cut != 1 || iScale == 0) ADs_compared_farAvg->Add(Sam_ad_comp_farAvg[iScale]);
	}

	ADs_compared_farAvg->GetXaxis()->SetLimits(0.25,8.25);
	ADs_compared_farAvg->GetYaxis()->SetRangeUser(62,72);
	ADs_compared_farAvg->Write();

/*	TCanvas *adEffs = new TCanvas("adEffs","AD Efficiencies");
	adEffs->cd();
		h_ADeff[0]->SetStats(0);
		adEffs->SetTitle("DT Efficiencies");
		h_ADeff[0]->GetXaxis()->SetTitle("AD Number");
		h_ADeff[0]->GetYaxis()->SetTitle("DT Efficiency");
		h_ADeff[0]->GetYaxis()->SetRangeUser(100*average[0]-yScale, 100*average[0]+yScale);

		h_ADeff[0]->SetLineColor(colors[1]);
		h_ADeff[1]->SetLineColor(colors[1]);
		h_Sam_ADeff[0]->SetLineColor(colors[0]);
		h_Sam_ADeff[1]->SetLineColor(colors[0]);
			h_ADeff[0]->SetMarkerStyle(20);
			h_ADeff[1]->SetMarkerStyle(25);
			h_ADeff[1]->SetLineStyle(7);
			h_Sam_ADeff[0]->SetMarkerStyle(20);
			h_Sam_ADeff[1]->SetMarkerStyle(25);
			h_Sam_ADeff[1]->SetLineStyle(7);
		h_ADeff[0]->SetMarkerColor(colors[1]);
		h_ADeff[1]->SetMarkerColor(colors[1]);
		h_Sam_ADeff[0]->SetMarkerColor(colors[0]);
		h_Sam_ADeff[1]->SetMarkerColor(colors[0]);

		h_ADeff[0]->Draw("e1x0");
	//	h_ADeff[1]->Draw("e1x0same");
		h_Sam_ADeff[0]->Draw("e1x0same");
	//	if(SamEfficiency_rate[0] != 0) h_Sam_ADeff[1]->Draw("e1x0same");

	adEffs->BuildLegend();
//	adEffs->SetGridx();
//	adEffs->SetGridy();
	if(Ep_cut == 0) sprintf(title,"./combinedEfficiency/DTeff_plots/adEfficiencies_DT%.0f_t%.0f_d%.0f.png",DTcut,max_x_value,max_y_value);
	if(Ep_cut == 1) sprintf(title,"./combinedEfficiency/DTeff_plots/adEfficiencies_DT%.0f_t%.0f_d%.0f_Ep.png",DTcut,max_x_value,max_y_value);
	adEffs->Print(title);

	TCanvas *adEffs_farAvg = new TCanvas("adEffs_farAvg","AD Efficiencies");
	adEffs_farAvg->cd();
		h_ADeff_farAvg[0]->SetStats(0);
		adEffs_farAvg->SetTitle("DT Efficiencies - Far ADs averaged");
		h_ADeff_farAvg[0]->GetXaxis()->SetTitle("AD Number");
		h_ADeff_farAvg[0]->GetYaxis()->SetTitle("DT Efficiency");
		h_ADeff_farAvg[0]->GetYaxis()->SetRangeUser(100*average[0]-yScale_farAvg, 100*average[0]+yScale_farAvg);

		h_ADeff_farAvg[0]->SetLineColor(colors[1]);
		h_ADeff_farAvg[1]->SetLineColor(colors[1]);
		h_Sam_ADeff_farAvg[0]->SetLineColor(colors[0]);
		h_Sam_ADeff_farAvg[1]->SetLineColor(colors[0]);
			h_ADeff_farAvg[0]->SetMarkerStyle(20);
			h_ADeff_farAvg[1]->SetMarkerStyle(25);
			h_ADeff_farAvg[1]->SetLineStyle(7);
			h_Sam_ADeff_farAvg[0]->SetMarkerStyle(20);
			h_Sam_ADeff_farAvg[1]->SetMarkerStyle(25);
			h_Sam_ADeff_farAvg[1]->SetLineStyle(7);
		h_ADeff_farAvg[0]->SetMarkerColor(colors[1]);
		h_ADeff_farAvg[1]->SetMarkerColor(colors[1]);
		h_Sam_ADeff_farAvg[0]->SetMarkerColor(colors[0]);
		h_Sam_ADeff_farAvg[1]->SetMarkerColor(colors[0]);

		h_ADeff_farAvg[0]->Draw("e1x0");
	//	h_ADeff_farAvg[1]->Draw("e1x0same");
		h_Sam_ADeff_farAvg[0]->Draw("e1x0same");
	//	if(SamEfficiency_rate[0] != 0) h_Sam_ADeff_farAvg[1]->Draw("e1x0same");

	adEffs_farAvg->BuildLegend();
//	adEffs_farAvg->SetGridx();
//	adEffs_farAvg->SetGridy();
	if(Ep_cut == 0) sprintf(title,"./combinedEfficiency/DTeff_plots/adEfficiencies_farAvg_DT%.0f_t%.0f_d%.0f.png",DTcut,max_x_value,max_y_value);
	if(Ep_cut == 1) sprintf(title,"./combinedEfficiency/DTeff_plots/adEfficiencies_farAvg_DT%.0f_t%.0f_d%.0f_Ep.png",DTcut,max_x_value,max_y_value);
	adEffs_farAvg->Print(title);*/

	TCanvas *eff_norm = new TCanvas("eff_norm","eff_norm");
	eff_norm->cd();
		for(int iad = 0; iad < maxAD; iad++){
			h_eff_DT_norm[iad]->SetStats(0);
			h_eff_DT_norm[iad]->SetLineColor(colors[iad]);
			if(iad == 0) h_eff_DT_norm[iad]->Draw("hist");
			else h_eff_DT_norm[iad]->Draw("hist same");
		}
//		h_eff_DT_norm[0]->GetYaxis()->SetRangeUser(-0.1,0.1);
		h_eff_DT_norm[0]->GetXaxis()->SetRangeUser(0.01,3);
		eff_norm->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/eff_DT_norm.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/eff_DT_Ep35_norm.png");
		eff_norm->Print(title);

	TCanvas *eff_rate = new TCanvas("eff_rate","eff_rate");
	eff_rate->cd();
		for(int iad = 0; iad < maxAD; iad++){
			h_eff_DT_rate[iad]->SetStats(0);
			h_eff_DT_rate[iad]->SetLineColor(colors[iad]);
			if(iad == 0) h_eff_DT_rate[iad]->Draw("hist");
			else h_eff_DT_rate[iad]->Draw("hist same");
		}
//		h_eff_DT_rate[0]->GetYaxis()->SetRangeUser(-0.1,0.1);
		h_eff_DT_rate[0]->GetXaxis()->SetRangeUser(0.01,3);
		eff_rate->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/eff_DT_rate.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/eff_DT_Ep35_rate.png");
		eff_rate->Print(title);

	TCanvas *eff_DTnorm = new TCanvas("eff_DTnorm","eff_DTnorm");
	eff_DTnorm->cd();
		for(int iad = 0; iad < maxAD; iad++){
			h_eff_DT_DTnorm[iad]->SetStats(0);
			h_eff_DT_DTnorm[iad]->SetLineColor(colors[iad]);
			if(iad == 0) h_eff_DT_DTnorm[iad]->Draw("hist");
			else h_eff_DT_DTnorm[iad]->Draw("hist same");
		}
//		h_eff_DT_DTnorm[0]->GetYaxis()->SetRangeUser(-0.05,0.05);
		h_eff_DT_DTnorm[0]->GetXaxis()->SetRangeUser(0.01,3);
		eff_DTnorm->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/eff_DT_DTnorm.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/eff_DT_Ep35_DTnorm.png");
		eff_DTnorm->Print(title);

	TCanvas *effDiff_norm = new TCanvas("effDiff_norm","effDiff_norm");
	effDiff_norm->cd();
		for(int iad = 0; iad < maxAD; iad++){
			h_effDiff_DT_norm[iad]->SetStats(0);
			h_effDiff_DT_norm[iad]->SetLineColor(colors[iad]);
//			if(iad == 0) h_effDiff_DT_norm[iad]->Draw("hist");
//			else h_effDiff_DT_norm[iad]->Draw("hist same");
		}
		for(int iad = 4; iad < maxAD; iad++){
			if(iad == 4) h_effDiff_DT_norm[iad]->Draw();
			else h_effDiff_DT_norm[iad]->Draw("hist same");
		}
		for(int iad = 0; iad < 4; iad++){
			h_effDiff_DT_norm[iad]->Draw("hist same");
		}
		h_effDiff_DT_norm[4]->GetYaxis()->SetRangeUser(-0.05,0.05);
		h_effDiff_DT_norm[4]->GetXaxis()->SetRangeUser(0.01,3);
		effDiff_norm->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_norm.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_norm.png");
		effDiff_norm->Print(title);

	TCanvas *effDiff_rate = new TCanvas("effDiff_rate","effDiff_rate"); //do the thing here
	effDiff_rate->cd();
		for(int iad = 0; iad < maxAD; iad++){
			h_effDiff_DT_rate[iad]->SetStats(0);
			h_effDiff_DT_rate[iad]->SetLineColor(colors[iad]);
		}
		for(int iad = 4; iad < maxAD; iad++){
			if(iad == 4) h_effDiff_DT_rate[iad]->Draw();
			else h_effDiff_DT_rate[iad]->Draw("hist same");
		}
		for(int iad = 0; iad < 4; iad++){
			h_effDiff_DT_rate[iad]->Draw("hist same");
		}
		h_effDiff_DT_rate[4]->GetYaxis()->SetRangeUser(-0.05,0.05);
		if(Ep_cut == 1) h_effDiff_DT_rate[4]->GetYaxis()->SetRangeUser(-0.006,0.006);
		h_effDiff_DT_rate[4]->GetXaxis()->SetRangeUser(0.01,3);
		effDiff_rate->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_rate.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_rate.png");
		effDiff_rate->Print(title);

	TCanvas *effDiff_DTnorm = new TCanvas("effDiff_DTnorm","effDiff_DTnorm");
	effDiff_DTnorm->cd();
		for(int iad = 0; iad < maxAD; iad++){
			h_effDiff_DT_DTnorm[iad]->SetStats(0);
			h_effDiff_DT_DTnorm[iad]->SetLineColor(colors[iad]);
		}
		for(int iad = 4; iad < maxAD; iad++){
			if(iad == 4) h_effDiff_DT_DTnorm[iad]->Draw();
			else h_effDiff_DT_DTnorm[iad]->Draw("hist same");
		}
		for(int iad = 0; iad < 4; iad++){
			h_effDiff_DT_DTnorm[iad]->Draw("hist same");
		}
		h_effDiff_DT_DTnorm[4]->GetYaxis()->SetRangeUser(-0.05,0.05);
		h_effDiff_DT_DTnorm[4]->GetXaxis()->SetRangeUser(0.01,3);
		effDiff_DTnorm->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_DTnorm.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_DTnorm.png");
		effDiff_DTnorm->Print(title);


	//*****Far ADs Averaged*****
	TCanvas *effDiff_norm_farAvg = new TCanvas("effDiff_norm_farAvg","effDiff_norm_farAvg");
	effDiff_norm_farAvg->cd();
			h_effDiff_DT_norm_farAvg->SetStats(0);
			h_effDiff_DT_norm_farAvg->SetLineColor(colors[5]);
			h_effDiff_DT_norm_farAvg->Draw();
		for(int iad = 0; iad < 4; iad++){
			h_effDiff_DT_norm[iad]->SetStats(0);
			h_effDiff_DT_norm[iad]->SetLineColor(colors[iad]);
			h_effDiff_DT_norm[iad]->Draw("hist same");
		}
		h_effDiff_DT_norm_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_effDiff_DT_norm_farAvg->GetYaxis()->SetRangeUser(-0.03,0.03);
		effDiff_norm_farAvg->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_norm_farAvg.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_norm_farAvg.png");
		effDiff_norm_farAvg->Print(title);

	TCanvas *effDiff_rate_farAvg = new TCanvas("effDiff_rate_farAvg","effDiff_rate_farAvg");
	effDiff_rate_farAvg->cd();
			h_effDiff_DT_rate_farAvg->SetStats(0);
			h_effDiff_DT_rate_farAvg->SetLineColor(colors[5]);
			h_effDiff_DT_rate_farAvg->Draw();
		for(int iad = 0; iad < 4; iad++){
			h_effDiff_DT_rate[iad]->SetStats(0);
			h_effDiff_DT_rate[iad]->SetLineColor(colors[iad]);
			h_effDiff_DT_rate[iad]->Draw("hist same");
		}
		h_effDiff_DT_rate_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_effDiff_DT_rate_farAvg->GetYaxis()->SetRangeUser(-0.03,0.03);
		effDiff_rate_farAvg->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_rate_farAvg.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_rate_farAvg.png");
		effDiff_rate_farAvg->Print(title);

	TCanvas *effDiff_DTnorm_farAvg = new TCanvas("effDiff_DTnorm_farAvg","effDiff_DTnorm_farAvg");
	effDiff_DTnorm_farAvg->cd();
			h_effDiff_DT_DTnorm_farAvg->SetStats(0);
			h_effDiff_DT_DTnorm_farAvg->SetLineColor(colors[5]);
			h_effDiff_DT_DTnorm_farAvg->Draw();
		for(int iad = 0; iad < 4; iad++){
			h_effDiff_DT_DTnorm[iad]->SetStats(0);
			h_effDiff_DT_DTnorm[iad]->SetLineColor(colors[iad]);
			h_effDiff_DT_DTnorm[iad]->Draw("hist same");
		}
		h_effDiff_DT_DTnorm_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_effDiff_DT_DTnorm_farAvg->GetYaxis()->SetRangeUser(-0.03,0.03);
		effDiff_DTnorm_farAvg->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_DTnorm_farAvg.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_DTnorm_farAvg.png");
		effDiff_DTnorm_farAvg->Print(title);

	//*****Near VS Far*****
	TCanvas *effDiff_norm_nearVSfar = new TCanvas("effDiff_norm_nearVSfar","effDiff_norm_nearVSfar");
	effDiff_norm_nearVSfar->cd();
			h_effDiff_DT_norm_farAvg->SetStats(0);
			h_effDiff_DT_norm_farAvg->SetLineColor(colors[5]);
			h_effDiff_DT_norm_farAvg->Draw("hist");
			h_effDiff_DT_norm_nearAvg->SetStats(0);
			h_effDiff_DT_norm_nearAvg->SetLineColor(colors[0]);
			h_effDiff_DT_norm_nearAvg->Draw("hist same");
		h_effDiff_DT_norm_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_effDiff_DT_norm_farAvg->GetYaxis()->SetRangeUser(-0.03,0.03);
		effDiff_norm_nearVSfar->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_norm_nearVSfar.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_norm_nearVsFar.png");
		effDiff_norm_nearVSfar->Print(title);

	TCanvas *effDiff_rate_nearVSfar = new TCanvas("effDiff_rate_nearVSfar","effDiff_rate_nearVSfar");
	effDiff_rate_nearVSfar->cd();
			h_effDiff_DT_rate_farAvg->SetStats(0);
			h_effDiff_DT_rate_farAvg->SetLineColor(colors[5]);
			h_effDiff_DT_rate_farAvg->Draw("hist");
			h_effDiff_DT_rate_nearAvg->SetStats(0);
			h_effDiff_DT_rate_nearAvg->SetLineColor(colors[0]);
			h_effDiff_DT_rate_nearAvg->Draw("hist same");
		h_effDiff_DT_rate_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_effDiff_DT_rate_farAvg->GetYaxis()->SetRangeUser(-0.03,0.03);
		effDiff_rate_nearVSfar->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_rate_nearVSfar.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_rate_nearVSfar.png");
		effDiff_rate_nearVSfar->Print(title);

	TCanvas *effDiff_DTnorm_nearVSfar = new TCanvas("effDiff_DTnorm_nearVSfar","effDiff_DTnorm_nearVSfar");
	effDiff_DTnorm_nearVSfar->cd();
			h_effDiff_DT_DTnorm_farAvg->SetStats(0);
			h_effDiff_DT_DTnorm_farAvg->SetLineColor(colors[5]);
			h_effDiff_DT_DTnorm_farAvg->Draw("hist");
			h_effDiff_DT_DTnorm_nearAvg->SetStats(0);
			h_effDiff_DT_DTnorm_nearAvg->SetLineColor(colors[0]);
			h_effDiff_DT_DTnorm_nearAvg->Draw("hist same");
		h_effDiff_DT_DTnorm_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_effDiff_DT_DTnorm_farAvg->GetYaxis()->SetRangeUser(-0.03,0.03);
		effDiff_DTnorm_nearVSfar->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/effDiff_DTnorm_nearVSfar.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/effDiff_DT_Ep35_DTnorm_nearVSfar.png");
		effDiff_DTnorm_nearVSfar->Print(title);


	//**************************************Relative plots********************************************
/*	TCanvas *rel_effDiff_norm = new TCanvas("rel_effDiff_norm","rel_effDiff_norm");
	rel_effDiff_norm->cd();
		for(int iad = 0; iad < maxAD; iad++){
			h_rel_effDiff_DT_norm[iad]->SetStats(0);
			h_rel_effDiff_DT_norm[iad]->SetLineColor(colors[iad]);
			if(iad == 0) h_rel_effDiff_DT_norm[iad]->Draw();
			else h_rel_effDiff_DT_norm[iad]->Draw("same");
		}
		h_rel_effDiff_DT_norm[0]->GetYaxis()->SetRangeUser(-0.1,0.1);
		h_rel_effDiff_DT_norm[0]->GetXaxis()->SetRangeUser(0.01,3);
		rel_effDiff_norm->BuildLegend();
		sprintf(title,"./combinedEfficiency/DTeff_plots/rel_effDiff_DT_norm.png");
		rel_effDiff_norm->Print(title);

	TCanvas *rel_effDiff_rate = new TCanvas("rel_effDiff_rate","rel_effDiff_rate");
	rel_effDiff_rate->cd();
		for(int iad = 0; iad < maxAD; iad++){
			h_rel_effDiff_DT_rate[iad]->SetStats(0);
			h_rel_effDiff_DT_rate[iad]->SetLineColor(colors[iad]);
			if(iad == 0) h_rel_effDiff_DT_rate[iad]->Draw();
			else h_rel_effDiff_DT_rate[iad]->Draw("same");
		}
		h_rel_effDiff_DT_rate[0]->GetYaxis()->SetRangeUser(-0.1,0.1);
		h_rel_effDiff_DT_rate[0]->GetXaxis()->SetRangeUser(0.01,3);
		rel_effDiff_rate->BuildLegend();
		sprintf(title,"./combinedEfficiency/DTeff_plots/rel_effDiff_DT_rate.png");
		rel_effDiff_rate->Print(title);*/

	TCanvas *rel_effDiff_DTnorm = new TCanvas("rel_effDiff_DTnorm","rel_effDiff_DTnorm");
	rel_effDiff_DTnorm->cd();
		h_rel_effDiff_DT_DTnorm[4]->Draw();
		for(int iad = 0; iad < maxAD; iad++){
			h_rel_effDiff_DT_DTnorm[iad]->SetStats(0);
			h_rel_effDiff_DT_DTnorm[iad]->SetLineColor(colors[iad]);
//			h_rel_effDiff_DT_DTnorm[iad]->Draw("same");
			if(iad != 4) h_rel_effDiff_DT_DTnorm[iad]->Draw("hist same");
		}
		h_rel_effDiff_DT_DTnorm[4]->GetYaxis()->SetRangeUser(-0.02,0.02);
		h_rel_effDiff_DT_DTnorm[4]->GetXaxis()->SetRangeUser(0.01,3);
		rel_effDiff_DTnorm->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/rel_effDiff_DT_DTnorm.png");
		rel_effDiff_DTnorm->Print(title);


	//*****Far ADs Averaged*****
	TCanvas *rel_effDiff_norm_farAvg = new TCanvas("rel_effDiff_norm_farAvg","rel_effDiff_norm_farAvg");
	rel_effDiff_norm_farAvg->cd();
			h_rel_effDiff_DT_norm_farAvg->SetStats(0);
			h_rel_effDiff_DT_norm_farAvg->SetLineColor(colors[5]);
			h_rel_effDiff_DT_norm_farAvg->Draw("hist");
		for(int iad = 0; iad < 4; iad++){
			h_rel_effDiff_DT_norm[iad]->SetStats(0);
			h_rel_effDiff_DT_norm[iad]->SetLineColor(colors[iad]);
			h_rel_effDiff_DT_norm[iad]->Draw("hist same");
		}
		h_rel_effDiff_DT_norm_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_rel_effDiff_DT_norm_farAvg->GetYaxis()->SetRangeUser(-0.02,0.02);
		rel_effDiff_norm_farAvg->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/rel_effDiff_DT_norm_farAvg.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/rel_effDiff_DT_Ep35_norm_farAvg.png");
		rel_effDiff_norm_farAvg->Print(title);

	TCanvas *rel_effDiff_rate_farAvg = new TCanvas("rel_effDiff_rate_farAvg","rel_effDiff_rate_farAvg");
	rel_effDiff_rate_farAvg->cd();
			h_rel_effDiff_DT_rate_farAvg->SetStats(0);
			h_rel_effDiff_DT_rate_farAvg->SetLineColor(colors[5]);
			h_rel_effDiff_DT_rate_farAvg->Draw("hist");
		for(int iad = 0; iad < 4; iad++){
			h_rel_effDiff_DT_rate[iad]->SetStats(0);
			h_rel_effDiff_DT_rate[iad]->SetLineColor(colors[iad]);
			h_rel_effDiff_DT_rate[iad]->Draw("hist same");
		}
		h_rel_effDiff_DT_rate_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_rel_effDiff_DT_rate_farAvg->GetYaxis()->SetRangeUser(-0.02,0.02);
		rel_effDiff_rate_farAvg->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/rel_effDiff_DT_rate_farAvg.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/rel_effDiff_DT_Ep35_rate_farAvg.png");
		rel_effDiff_rate_farAvg->Print(title);

	TCanvas *rel_effDiff_DTnorm_farAvg = new TCanvas("rel_effDiff_DTnorm_farAvg","rel_effDiff_DTnorm_farAvg");
	rel_effDiff_DTnorm_farAvg->cd();
			h_rel_effDiff_DT_DTnorm_farAvg->SetStats(0);
			h_rel_effDiff_DT_DTnorm_farAvg->SetLineColor(colors[5]);
		//	h_rel_effDiff_DT_DTnorm_farAvg->SetFillColor(colors[5]-10);
		//	h_rel_effDiff_DT_DTnorm_farAvg->SetFillStyle(3003);
			h_rel_effDiff_DT_DTnorm_farAvg->Draw();
		for(int iad = 0; iad < 4; iad++){
			h_rel_effDiff_DT_DTnorm[iad]->SetStats(0);
			h_rel_effDiff_DT_DTnorm[iad]->SetLineColor(colors[iad]);
			h_rel_effDiff_DT_DTnorm[iad]->Draw("hist same");
		}
		h_rel_effDiff_DT_DTnorm_farAvg->GetXaxis()->SetRangeUser(0.01,3);
		h_rel_effDiff_DT_DTnorm_farAvg->GetYaxis()->SetRangeUser(-0.02,0.02);
//		h_rel_effDiff_DT_DTnorm_farAvg->GetYaxis()->SetRangeUser(-0.005,0.005);
		rel_effDiff_DTnorm_farAvg->BuildLegend();
		sprintf(title,"./AdSimpleNL/3sig/rel_effDiff_DT_DTnorm_farAvg.png");
		if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/rel_effDiff_DT_Ep35_DTnorm_farAvg.png");
		rel_effDiff_DTnorm_farAvg->Print(title);

	TCanvas *adEffs = new TCanvas("adEffs","AD Efficiencies");
	adEffs->cd();
		h_ADeff_DTnorm->SetStats(0);
		h_ADeff_norm->GetXaxis()->SetTitle("AD Number");
		h_ADeff_norm->GetYaxis()->SetTitle("DT Efficiency");
//		h_ADeff_norm->GetYaxis()->SetRangeUser(100*average[0]-yScale, 100*average[0]+yScale);

		h_ADeff_norm->SetLineColor(kBlue);
		h_ADeff_rate->SetLineColor(kRed);
		h_ADeff_DTnorm->SetLineColor(kBlack);
			h_ADeff_norm->SetMarkerStyle(20);
			h_ADeff_rate->SetMarkerStyle(21);
			h_ADeff_DTnorm->SetMarkerStyle(22);
		h_ADeff_norm->SetMarkerColor(kBlue);
		h_ADeff_rate->SetMarkerColor(kRed);
		h_ADeff_DTnorm->SetMarkerColor(kBlack);

		h_ADeff_DTnorm->Draw("e1x0");
		h_ADeff_norm->Draw("p hist same");
		h_ADeff_rate->Draw("p hist same");

	adEffs->BuildLegend();
//	adEffs->SetGridx();
//	adEffs->SetGridy();
	sprintf(title,"./AdSimpleNL/3sig/adEfficiencies_DT_BedaErrors.png");
	if(Ep_cut == 1) sprintf(title,"./AdSimpleNL/3sig/adEfficiencies_DT_Ep35_BedaErrors.png");
	adEffs->Print(title);


}

void delayed(int DT, int a){ //Elow = 1.5, Ehigh = 2.8 for Sam ... used to have this included: (, double Elow, double Ehigh)
	char name[64];
	int EH[8] = {1,1,2,2,3,3,3,3};
	int AD[8] = {1,2,1,2,1,2,3,4};

	const int NzBins = 20;
	const int NzPoints = 5;

	const int Nr2Bins = 20;
	const int Nr2Points = 5;

	double peak_rate[8];
	double sigma_rate[8];
	double alpha_rate[8];
	double lambda_rate[8];
	double efficiency_rate[8];
	double peak_norm[8];
	double sigma_norm[8];
	double alpha_norm[8];
	double lambda_norm[8];
	double efficiency_norm[8];
	double peak_DTnorm[8];
	double sigma_DTnorm[8];
	double alpha_DTnorm[8];
	double lambda_DTnorm[8];
	double efficiency_DTnorm[8];

	double peak_rate_Sam[8];
	double sigma_rate_Sam[8];
	double alpha_rate_Sam[8];
	double lambda_rate_Sam[8];
	double efficiency_rate_Sam[8];
	double peak_norm_Sam[8];
	double sigma_norm_Sam[8];
	double alpha_norm_Sam[8];
	double lambda_norm_Sam[8];
	double efficiency_norm_Sam[8];
	double peak_DTnorm_Sam[8];
	double sigma_DTnorm_Sam[8];
	double alpha_DTnorm_Sam[8];
	double lambda_DTnorm_Sam[8];
	double efficiency_DTnorm_Sam[8];

	double peakError_rate[8];
	double sigmaError_rate[8];
	double alphaError_rate[8];
	double lambdaError_rate[8];
	double efficiencyError_rate[8];
	double peakError_norm[8];
	double sigmaError_norm[8];
	double alphaError_norm[8];
	double lambdaError_norm[8];
	double efficiencyError_norm[8];
	double peakError_DTnorm[8];
	double sigmaError_DTnorm[8];
	double alphaError_DTnorm[8];
	double lambdaError_DTnorm[8];
	double efficiencyError_DTnorm[8];

	double peakError_rate_Sam[8];
	double sigmaError_rate_Sam[8];
	double alphaError_rate_Sam[8];
	double lambdaError_rate_Sam[8];
	double efficiencyError_rate_Sam[8];
	double peakError_norm_Sam[8];
	double sigmaError_norm_Sam[8];
	double alphaError_norm_Sam[8];
	double lambdaError_norm_Sam[8];
	double efficiencyError_norm_Sam[8];
	double peakError_DTnorm_Sam[8];
	double sigmaError_DTnorm_Sam[8];
	double alphaError_DTnorm_Sam[8];
	double lambdaError_DTnorm_Sam[8];
	double efficiencyError_DTnorm_Sam[8];

	double Ehigh_rate[8];
	double Elow_rate[8];
	double Ehigh_norm[8];
	double Elow_norm[8];
	double Ehigh_DTnorm[8];
	double Elow_DTnorm[8];

	double N_nom_rate[8];
	double N_nom_norm[8];
	double N_nom_DTnorm[8];
	double N_nom_rate_Sam[8];
	double N_nom_norm_Sam[8];
	double N_nom_DTnorm_Sam[8];
	double N_ext_rate[8];
	double N_ext_norm[8];
	double N_ext_DTnorm[8];
	double N_ext_rate_Sam[8];
	double N_ext_norm_Sam[8];
	double N_ext_DTnorm_Sam[8];

	double N_nom_rate_z[8][NzBins];
	double N_ext_rate_z[8][NzBins];
	double efficiency_rate_z[8][NzBins];
	double efficiencyError_rate_z[8][NzBins];

	double N_nom_norm_z[8][NzBins];
	double N_ext_norm_z[8][NzBins];
	double efficiency_norm_z[8][NzBins];
	double efficiencyError_norm_z[8][NzBins];

	double N_nom_DTnorm_z[8][NzBins];
	double N_ext_DTnorm_z[8][NzBins];
	double efficiency_DTnorm_z[8][NzBins];
	double efficiencyError_DTnorm_z[8][NzBins];

	double N_nom_rate_r2[8][Nr2Bins];
	double N_ext_rate_r2[8][Nr2Bins];
	double efficiency_rate_r2[8][Nr2Bins];
	double efficiencyError_rate_r2[8][Nr2Bins];

	double N_nom_norm_r2[8][Nr2Bins];
	double N_ext_norm_r2[8][Nr2Bins];
	double efficiency_norm_r2[8][Nr2Bins];
	double efficiencyError_norm_r2[8][Nr2Bins];

	double N_nom_DTnorm_r2[8][Nr2Bins];
	double N_ext_DTnorm_r2[8][Nr2Bins];
	double efficiency_DTnorm_r2[8][Nr2Bins];
	double efficiencyError_DTnorm_r2[8][Nr2Bins];

	double N_nom_rate_zVSr2[8][Nr2Bins][NzBins];
	double N_ext_rate_zVSr2[8][Nr2Bins][NzBins];
	double efficiency_rate_zVSr2[8][Nr2Bins][NzBins];
	double efficiencyError_rate_zVSr2[8][Nr2Bins][NzBins];

	double N_nom_norm_zVSr2[8][Nr2Bins][NzBins];
	double N_ext_norm_zVSr2[8][Nr2Bins][NzBins];
	double efficiency_norm_zVSr2[8][Nr2Bins][NzBins];
	double efficiencyError_norm_zVSr2[8][Nr2Bins][NzBins];

	double N_nom_DTnorm_zVSr2[8][Nr2Bins][NzBins];
	double N_ext_DTnorm_zVSr2[8][Nr2Bins][NzBins];
	double efficiency_DTnorm_zVSr2[8][Nr2Bins][NzBins];
	double efficiencyError_DTnorm_zVSr2[8][Nr2Bins][NzBins];

		double N_nom_rate_zVSr2_near[Nr2Bins][NzBins];
		double N_ext_rate_zVSr2_near[Nr2Bins][NzBins];
		double efficiency_rate_zVSr2_near[Nr2Bins][NzBins];
		double efficiencyError_rate_zVSr2_near[Nr2Bins][NzBins];

		double N_nom_norm_zVSr2_near[Nr2Bins][NzBins];
		double N_ext_norm_zVSr2_near[Nr2Bins][NzBins];
		double efficiency_norm_zVSr2_near[Nr2Bins][NzBins];
		double efficiencyError_norm_zVSr2_near[Nr2Bins][NzBins];

		double N_nom_DTnorm_zVSr2_near[Nr2Bins][NzBins];
		double N_ext_DTnorm_zVSr2_near[Nr2Bins][NzBins];
		double efficiency_DTnorm_zVSr2_near[Nr2Bins][NzBins];
		double efficiencyError_DTnorm_zVSr2_near[Nr2Bins][NzBins];

		double N_nom_rate_zVSr2_far[Nr2Bins][NzBins];
		double N_ext_rate_zVSr2_far[Nr2Bins][NzBins];
		double efficiency_rate_zVSr2_far[Nr2Bins][NzBins];
		double efficiencyError_rate_zVSr2_far[Nr2Bins][NzBins];

		double N_nom_norm_zVSr2_far[Nr2Bins][NzBins];
		double N_ext_norm_zVSr2_far[Nr2Bins][NzBins];
		double efficiency_norm_zVSr2_far[Nr2Bins][NzBins];
		double efficiencyError_norm_zVSr2_far[Nr2Bins][NzBins];

		double N_nom_DTnorm_zVSr2_far[Nr2Bins][NzBins];
		double N_ext_DTnorm_zVSr2_far[Nr2Bins][NzBins];
		double efficiency_DTnorm_zVSr2_far[Nr2Bins][NzBins];
		double efficiencyError_DTnorm_zVSr2_far[Nr2Bins][NzBins];

	for(int i=0; i<8; i++){
		peak_rate[i] = 0;
		sigma_rate[i] = 0;
		alpha_rate[i] = 0;
		lambda_rate[i] = 0;
		efficiency_rate[i]=0;
		peak_norm[i] = 0;
		sigma_norm[i] = 0;
		alpha_norm[i] = 0;
		lambda_norm[i] = 0;
		efficiency_norm[i]=0;
		peak_DTnorm[i] = 0;
		sigma_DTnorm[i] = 0;
		alpha_DTnorm[i] = 0;
		lambda_DTnorm[i] = 0;
		efficiency_DTnorm[i]=0;

		peak_rate_Sam[i] = 0;
		sigma_rate_Sam[i] = 0;
		alpha_rate_Sam[i] = 0;
		lambda_rate_Sam[i] = 0;
		efficiency_rate_Sam[i]=0;
		peak_norm_Sam[i] = 0;
		sigma_norm_Sam[i] = 0;
		alpha_norm_Sam[i] = 0;
		lambda_norm_Sam[i] = 0;
		efficiency_norm_Sam[i]=0;
		peak_DTnorm_Sam[i] = 0;
		sigma_DTnorm_Sam[i] = 0;
		alpha_DTnorm_Sam[i] = 0;
		lambda_DTnorm_Sam[i] = 0;
		efficiency_DTnorm_Sam[i]=0;

		peakError_rate[i] = 0;
		sigmaError_rate[i] = 0;
		alphaError_rate[i] = 0;
		lambdaError_rate[i] = 0;
		efficiencyError_rate[i]=0;
		peakError_norm[i] = 0;
		sigmaError_norm[i] = 0;
		alphaError_norm[i] = 0;
		lambdaError_norm[i] = 0;
		efficiencyError_norm[i]=0;
		peakError_DTnorm[i] = 0;
		sigmaError_DTnorm[i] = 0;
		alphaError_DTnorm[i] = 0;
		lambdaError_DTnorm[i] = 0;
		efficiencyError_DTnorm[i]=0;

		peakError_rate_Sam[i] = 0;
		sigmaError_rate_Sam[i] = 0;
		alphaError_rate_Sam[i] = 0;
		lambdaError_rate_Sam[i] = 0;
		efficiencyError_rate_Sam[i]=0;
		peakError_norm_Sam[i] = 0;
		sigmaError_norm_Sam[i] = 0;
		alphaError_norm_Sam[i] = 0;
		lambdaError_norm_Sam[i] = 0;
		efficiencyError_norm_Sam[i]=0;
		peakError_DTnorm_Sam[i] = 0;
		sigmaError_DTnorm_Sam[i] = 0;
		alphaError_DTnorm_Sam[i] = 0;
		lambdaError_DTnorm_Sam[i] = 0;
		efficiencyError_DTnorm_Sam[i]=0;

		Ehigh_rate[i] = 0;
		Elow_rate[i] = 0;
		Ehigh_norm[i] = 0;
		Elow_norm[i] = 0;
		Ehigh_DTnorm[i] = 0;
		Elow_DTnorm[i] = 0;

		N_nom_rate[i]=0;
		N_nom_norm[i]=0;
		N_nom_DTnorm[i]=0;
		N_nom_rate_Sam[i]=0;
		N_nom_norm_Sam[i]=0;
		N_nom_DTnorm_Sam[i]=0;
		N_ext_rate[i]=0;
		N_ext_norm[i]=0;
		N_ext_DTnorm[i]=0;
		N_ext_rate_Sam[i]=0;
		N_ext_norm_Sam[i]=0;
		N_ext_DTnorm_Sam[i]=0;

		for(int iz = 0; iz < NzBins; iz++){
			N_nom_rate_z[i][iz]=0;
			N_ext_rate_z[i][iz]=0;
			efficiency_rate_z[i][iz]=0;
			efficiencyError_rate_z[i][iz]=0;

			N_nom_norm_z[i][iz]=0;
			N_ext_norm_z[i][iz]=0;
			efficiency_norm_z[i][iz]=0;
			efficiencyError_norm_z[i][iz]=0;

			N_nom_DTnorm_z[i][iz]=0;
			N_ext_DTnorm_z[i][iz]=0;
			efficiency_DTnorm_z[i][iz]=0;
			efficiencyError_DTnorm_z[i][iz]=0;
		}

		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			N_nom_rate_r2[i][ir2]=0;
			N_ext_rate_r2[i][ir2]=0;
			efficiency_rate_r2[i][ir2]=0;
			efficiencyError_rate_r2[i][ir2]=0;

			N_nom_norm_r2[i][ir2]=0;
			N_ext_norm_r2[i][ir2]=0;
			efficiency_norm_r2[i][ir2]=0;
			efficiencyError_norm_r2[i][ir2]=0;

			N_nom_DTnorm_r2[i][ir2]=0;
			N_ext_DTnorm_r2[i][ir2]=0;
			efficiency_DTnorm_r2[i][ir2]=0;
			efficiencyError_DTnorm_r2[i][ir2]=0;
		}

		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			for(int iz = 0; iz < NzBins; iz++){
				N_nom_rate_zVSr2[i][ir2][iz]=0;
				N_ext_rate_zVSr2[i][ir2][iz]=0;
				efficiency_rate_zVSr2[i][ir2][iz]=0;
				efficiencyError_rate_zVSr2[i][ir2][iz]=0;

				N_nom_norm_zVSr2[i][ir2][iz]=0;
				N_ext_norm_zVSr2[i][ir2][iz]=0;
				efficiency_norm_zVSr2[i][ir2][iz]=0;
				efficiencyError_norm_zVSr2[i][ir2][iz]=0;

				N_nom_DTnorm_zVSr2[i][ir2][iz]=0;
				N_ext_DTnorm_zVSr2[i][ir2][iz]=0;
				efficiency_DTnorm_zVSr2[i][ir2][iz]=0;
				efficiencyError_DTnorm_zVSr2[i][ir2][iz]=0;

				N_nom_rate_zVSr2_near[ir2][iz]=0;
				N_ext_rate_zVSr2_near[ir2][iz]=0;
				efficiency_rate_zVSr2_near[ir2][iz]=0;
				efficiencyError_rate_zVSr2_near[ir2][iz]=0;

				N_nom_norm_zVSr2_near[ir2][iz]=0;
				N_ext_norm_zVSr2_near[ir2][iz]=0;
				efficiency_norm_zVSr2_near[ir2][iz]=0;
				efficiencyError_norm_zVSr2_near[ir2][iz]=0;

				N_nom_DTnorm_zVSr2_near[ir2][iz]=0;
				N_ext_DTnorm_zVSr2_near[ir2][iz]=0;
				efficiency_DTnorm_zVSr2_near[ir2][iz]=0;
				efficiencyError_DTnorm_zVSr2_near[ir2][iz]=0;

				N_nom_rate_zVSr2_far[ir2][iz]=0;
				N_ext_rate_zVSr2_far[ir2][iz]=0;
				efficiency_rate_zVSr2_far[ir2][iz]=0;
				efficiencyError_rate_zVSr2_far[ir2][iz]=0;

				N_nom_norm_zVSr2_far[ir2][iz]=0;
				N_ext_norm_zVSr2_far[ir2][iz]=0;
				efficiency_norm_zVSr2_far[ir2][iz]=0;
				efficiencyError_norm_zVSr2_far[ir2][iz]=0;

				N_nom_DTnorm_zVSr2_far[ir2][iz]=0;
				N_ext_DTnorm_zVSr2_far[ir2][iz]=0;
				efficiency_DTnorm_zVSr2_far[ir2][iz]=0;
				efficiencyError_DTnorm_zVSr2_far[ir2][iz]=0;
			}
		}
	}


//Making output file:
	char outputname[64];
	char title[64];


        sprintf(outputname,"../nH_files/delayedEnergy.root");
	TFile* outfile=new TFile(outputname, "RECREATE");

		gStyle->SetOptFit(1111);

	TCanvas *z_slice_Ed_perAD = new TCanvas("z_slice_Ed_perAD","z_slice_Ed_perAD");
	z_slice_Ed_perAD->Divide(4,2);

	TCanvas *r2_slice_Ed_perAD = new TCanvas("r2_slice_Ed_perAD","r2_slice_Ed_perAD");
	r2_slice_Ed_perAD->Divide(4,2);

		TH1D* h_efficiency_sigma_rate[8];
		TH1D* h_efficiency_sigma_norm[8];
		TH1D* h_efficiency_sigma_DTnorm[8];

		TH1F* h_Edelayed_sub_rate_z[8][NzPoints];
		TH1F* h_Edelayed_sub_norm_z[8][NzPoints];
		TH1F* h_Edelayed_sub_DTnorm_z[8][NzPoints];

		TH1F* h_Edelayed_sub_rate_r2[8][Nr2Points];
		TH1F* h_Edelayed_sub_norm_r2[8][Nr2Points];
		TH1F* h_Edelayed_sub_DTnorm_r2[8][Nr2Points];

		TH1F* h_Edelayed_sub_rate_zVSr2[8][Nr2Points][NzPoints];
		TH1F* h_Edelayed_sub_norm_zVSr2[8][Nr2Points][NzPoints];
		TH1F* h_Edelayed_sub_DTnorm_zVSr2[8][Nr2Points][NzPoints];
		for(int iad = 0; iad < 8; iad++){
			h_efficiency_sigma_rate[iad] = new TH1D(Form("h_effiency_sigma_rate_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_sigma_rate_eh%dad%d",EH[iad],AD[iad]),500,0,5);
			h_efficiency_sigma_norm[iad] = new TH1D(Form("h_effiency_sigma_norm_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_sigma_norm_eh%dad%d",EH[iad],AD[iad]),500,0,5);
			h_efficiency_sigma_DTnorm[iad] = new TH1D(Form("h_effiency_sigma_DTnorm_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_sigma_DTnorm_eh%dad%d",EH[iad],AD[iad]),500,0.005,5.005);
		}

	for(int iad = 0; iad < 8; iad++){
		//Getting files:
			sprintf(name,"../nH_files/SubtractedAccidentals_1500_EH%dAD%d.root",EH[iad],AD[iad]);
			TFile *subFile = new TFile(name);
			
			sprintf(name, "../nH_files/TotaledPlots_EH%d_1500.root",EH[iad]);
			TFile *ibdFile = new TFile(name);

			sprintf(name, "../nH_files/TotaledSingles_1500_EH%d.root",EH[iad]);
			TFile *accFile = new TFile(name);

		//Getting the histograms:
			sprintf(name, "h_total_delayed_energy_DT800_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_Ep35_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_total_delayed_energy_DT800_ad%d", AD[iad]);
			TH1F* delayed_rate = (TH1F*)ibdFile->Get(name);
			delayed_rate->SetName(Form("h_Edelayed_subtract_DT800_eh%dad%d",EH[iad],AD[iad]));
			TH1F* delayed_norm = (TH1F*)ibdFile->Get(name);
			delayed_norm->SetName(Form("h_Edelayed_subtract_DT800_norm_eh%dad%d",EH[iad],AD[iad]));
			TH1F* delayed_DTnorm = (TH1F*)ibdFile->Get(name);
			delayed_DTnorm->SetName(Form("h_Edelayed_subtract_DT800_DTnorm_eh%dad%d",EH[iad],AD[iad]));
			TH1F* delayed_IBD = (TH1F*)ibdFile->Get(name);
			delayed_IBD->SetName(Form("h_Edelayed_IBD_DT800_eh%dad%d",EH[iad],AD[iad]));
			
			sprintf(name, "h_total_delayed_energy_fine_DT800_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_Ep35_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_total_delayed_energy_fine_DT800_ad%d", AD[iad]);
			TH1F* delayed_fine_rate = (TH1F*)ibdFile->Get(name);
			delayed_fine_rate->SetName(Form("h_Edelayed_subtract_fine_DT800_eh%dad%d",EH[iad],AD[iad]));
			TH1F* delayed_fine_norm = (TH1F*)ibdFile->Get(name);
			delayed_fine_norm->SetName(Form("h_Edelayed_subtract_fine_DT800_norm_eh%dad%d",EH[iad],AD[iad]));
			TH1F* delayed_fine_DTnorm = (TH1F*)ibdFile->Get(name);
			delayed_fine_DTnorm->SetName(Form("h_Edelayed_subtract_fine_DT800_DTnorm_eh%dad%d",EH[iad],AD[iad]));
			TH1F* delayed_fine_IBD = (TH1F*)ibdFile->Get(name);
			delayed_fine_IBD->SetName(Form("h_Edelayed_IBD_fine_DT800_eh%dad%d",EH[iad],AD[iad]));
			
			sprintf(name, "h_Edelayed_subtract_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_Ep35_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_total_delayed_energy_DT800_scaled_ad%d", AD[iad]);
			TH1F* acc_rate = (TH1F*)accFile->Get(name);
			delayed_rate->Add(acc_rate, -1);
			
			sprintf(name, "h_Edelayed_subtract_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_Ep35_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_total_delayed_energy_DT800_norm_ad%d", AD[iad]);
			TH1F* acc_norm = (TH1F*)accFile->Get(name);
			delayed_norm->Add(acc_norm, -1);

			sprintf(name, "h_Edelayed_subtract_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_Ep35_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_total_delayed_energy_DT800_DTnorm_ad%d", AD[iad]);
			TH1F* acc_DTnorm = (TH1F*)accFile->Get(name);
			delayed_DTnorm->Add(acc_DTnorm, -1);
			

			sprintf(name, "h_total_delayed_energy_fine_DT800_scaled_1_ad%d", AD[iad]);
			TH1F* acc_rate_fine_1 = (TH1F*)accFile->Get(name);
			delayed_fine_rate->Add(acc_rate_fine_1, -1);
			sprintf(name, "h_total_delayed_energy_fine_DT800_scaled_2_ad%d", AD[iad]);
			TH1F* acc_rate_fine_2 = (TH1F*)accFile->Get(name);
			delayed_fine_rate->Add(acc_rate_fine_2, -1);

			sprintf(name, "h_total_delayed_energy_fine_DT800_norm_1_ad%d", AD[iad]);
			TH1F* acc_norm_fine_1 = (TH1F*)accFile->Get(name);
			delayed_fine_norm->Add(acc_norm_fine_1, -1);
			sprintf(name, "h_total_delayed_energy_fine_DT800_norm_2_ad%d", AD[iad]);
			TH1F* acc_norm_fine_2 = (TH1F*)accFile->Get(name);
			delayed_fine_norm->Add(acc_norm_fine_2, -1);

			sprintf(name, "h_total_delayed_energy_fine_DT800_DTnorm_1_ad%d", AD[iad]);
			TH1F* acc_DTnorm_fine_1 = (TH1F*)accFile->Get(name);
			delayed_fine_DTnorm->Add(acc_DTnorm_fine_1, -1);
			sprintf(name, "h_total_delayed_energy_fine_DT800_DTnorm_2_ad%d", AD[iad]);
			TH1F* acc_DTnorm_fine_2 = (TH1F*)accFile->Get(name);
			delayed_fine_DTnorm->Add(acc_DTnorm_fine_2, -1);
			
		//		sprintf(name, "h_Edelayed_IBD_ad%d", AD[iad]); //IBD histogram for error calculation
		//		if(DT == 1) sprintf(name, "h_Edelayed_IBD_DT800_Ep35_ad%d", AD[iad]);
	//			if(DT == 1) sprintf(name, "h_Edelayed_IBD_fine_DT800_Ep35_ad%d", AD[iad]);
		//		if(DT == 1) sprintf(name, "h_Edelayed_IBD_fine_DT800_ad%d", AD[iad]);
		//		TH1F* delayed_IBD = (TH1F*)subFile->Get(name);




	/*		sprintf(name, "h_Edelayed_subtract_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_Ep35_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_ad%d", AD[iad]);
			TH1F* delayed_rate = (TH1F*)subFile->Get(name);
			delayed_rate->SetName(Form("h_Edelayed_subtract_fine_DT800_eh%dad%d",EH[iad],AD[iad]));

			sprintf(name, "h_Edelayed_subtract_norm_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_Ep35_norm_ad%d", AD[iad]);
//			if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_norm_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_norm_ad%d", AD[iad]);
			TH1F* delayed_norm = (TH1F*)subFile->Get(name);

			sprintf(name, "h_Edelayed_subtract_DTnorm_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_Ep35_DTnorm_ad%d", AD[iad]);
//			if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_DTnorm_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_DTnorm_ad%d", AD[iad]);
			TH1F* delayed_DTnorm = (TH1F*)subFile->Get(name);*/



/*			sprintf(name, "h_Edelayed_subtract_fine_ad%d", AD[iad]);
	//		if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_ad%d", AD[iad]);
			TH1F* delayed_fine_rate = (TH1F*)subFile->Get(name);
//			delayed_fine_rate->SetName(Form("h_Edelayed_subtract_fine_DT800_Ep35_eh%dad%d",EH[iad],AD[iad]));
			delayed_fine_rate->SetName(Form("h_Edelayed_subtract_fine_DT800_eh%dad%d",EH[iad],AD[iad]));

			sprintf(name, "h_Edelayed_subtract_fine_norm_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_norm_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_norm_ad%d", AD[iad]);
			TH1F* delayed_fine_norm = (TH1F*)subFile->Get(name);
			delayed_fine_norm->SetName(Form("h_Edelayed_subtract_fine_DT800_norm_eh%dad%d",EH[iad],AD[iad]));

			sprintf(name, "h_Edelayed_subtract_fine_DTnorm_ad%d", AD[iad]);
		//	if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_Ep35_DTnorm_ad%d", AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_fine_DT800_DTnorm_ad%d", AD[iad]);
			TH1F* delayed_fine_DTnorm = (TH1F*)subFile->Get(name);
			delayed_fine_DTnorm->SetName(Form("h_Edelayed_subtract_fine_DT800_DTnorm_eh%dad%d",EH[iad],AD[iad]));

				sprintf(name, "h_Edelayed_IBD_fine_ad%d", AD[iad]); //IBD histogram for error calculation
		//		if(DT == 1) sprintf(name, "h_Edelayed_IBD_fine_DT800_Ep35_ad%d", AD[iad]);
				if(DT == 1) sprintf(name, "h_Edelayed_IBD_fine_DT800_ad%d", AD[iad]);
				TH1F* delayed_fine_IBD = (TH1F*)subFile->Get(name);
				delayed_fine_IBD->SetName(Form("h_Edelayed_IBD_fine_DT800_eh%dad%d",EH[iad],AD[iad]));*/

			cout << "AD"<< iad+1 << endl<< "Rate:\t" << delayed_rate->Integral(delayed_rate->FindBin(2.7),delayed_rate->FindBin(2.8)) << "\twhich is:\t" << 100*delayed_rate->Integral(delayed_rate->FindBin(2.7),delayed_rate->FindBin(2.8))/delayed_rate->Integral(delayed_rate->FindBin(1.5),delayed_rate->FindBin(2.8)) << endl << "Norm:\t" << delayed_norm->Integral(delayed_norm->FindBin(2.7),delayed_norm->FindBin(2.8)) << "\twhich is:\t" << 100*delayed_norm->Integral(delayed_norm->FindBin(2.7),delayed_norm->FindBin(2.8))/delayed_norm->Integral(delayed_norm->FindBin(1.5),delayed_norm->FindBin(2.8)) << endl << "DTnorm:\t" << delayed_DTnorm->Integral(delayed_DTnorm->FindBin(2.7),delayed_DTnorm->FindBin(2.8)) << "\twhich is:\t" << 100*delayed_DTnorm->Integral(delayed_DTnorm->FindBin(2.7),delayed_DTnorm->FindBin(2.8))/delayed_DTnorm->Integral(delayed_DTnorm->FindBin(1.5),delayed_DTnorm->FindBin(2.8)) << endl;

/*			//Getting histograms of the z-dependence
		TH1F* h_Edelayed_ibd_z[NzBins];
		TH1F* h_Edelayed_sub_rate_z[NzBins];
		TH1F* h_Edelayed_ibd_z_points[NzPoints];
		TH1F* h_Edelayed_sub_rate_z_points[NzPoints];
		TH1F* h_Edelayed_sub_norm_z[NzBins];
		TH1F* h_Edelayed_sub_norm_z_points[NzPoints];
		TH1F* h_Edelayed_sub_DTnorm_z[NzBins];
		TH1F* h_Edelayed_sub_DTnorm_z_points[NzPoints];
		for(int iz = 0; iz < NzBins; iz++){
			sprintf(name, "h_Edelayed_ibd_z_ad%d_iz%d", AD[iad], iz+1);
			if(DT == 1) sprintf(name, "h_Edelayed_ibd_DT800_z_ad%d_iz%d", AD[iad], iz+1);
			h_Edelayed_ibd_z[iz] = (TH1F*)subFile->Get(name);

			sprintf(name, "h_Edelayed_subtract_rate_z_ad%d_iz%d", AD[iad], iz+1);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_rate_DT800_z_ad%d_iz%d", AD[iad], iz+1);
			h_Edelayed_sub_rate_z[iz] = (TH1F*)subFile->Get(name);

			sprintf(name, "h_Edelayed_subtract_norm_z_ad%d_iz%d", AD[iad], iz+1);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_norm_DT800_z_ad%d_iz%d", AD[iad], iz+1);
			h_Edelayed_sub_norm_z[iz] = (TH1F*)subFile->Get(name);

			sprintf(name, "h_Edelayed_subtract_DTnorm_z_ad%d_iz%d", AD[iad], iz+1);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_DTnorm_DT800_z_ad%d_iz%d", AD[iad], iz+1);
			h_Edelayed_sub_DTnorm_z[iz] = (TH1F*)subFile->Get(name);

			for(int iBin = 0; iBin < (h_Edelayed_ibd_z[iz]->GetNbinsX())+1; iBin++){
				h_Edelayed_ibd_z[iz]->SetBinError(iBin, sqrt(h_Edelayed_ibd_z[iz]->GetBinContent(iBin)));
				h_Edelayed_sub_rate_z[iz]->SetBinError(iBin, sqrt(h_Edelayed_ibd_z[iz]->GetBinContent(iBin)));
				h_Edelayed_sub_norm_z[iz]->SetBinError(iBin, sqrt(h_Edelayed_ibd_z[iz]->GetBinContent(iBin)));
				h_Edelayed_sub_DTnorm_z[iz]->SetBinError(iBin, sqrt(h_Edelayed_ibd_z[iz]->GetBinContent(iBin)));
			}
		}*/

/*		for(int iz = 0; iz < NzBins; iz++){
			if((iz*NzPoints)%NzBins == 0) h_Edelayed_sub_rate_z_points[int(iz*NzPoints/NzBins)]=(TH1F*)h_Edelayed_sub_rate_z[iz]->Clone();
			else h_Edelayed_sub_rate_z_points[int(iz*NzPoints/NzBins)]->Add(h_Edelayed_sub_rate_z[iz]);

			if((iz*NzPoints)%NzBins == 0) h_Edelayed_ibd_z_points[int(iz*NzPoints/NzBins)]=(TH1F*)h_Edelayed_ibd_z[iz]->Clone();
			else h_Edelayed_ibd_z_points[int(iz*NzPoints/NzBins)]->Add(h_Edelayed_ibd_z[iz]);

			if((iz*NzPoints)%NzBins == 0){
				sprintf(name, "h_Edelayed_sub_rate_z_points_eh%dad%d_iz%d", EH[iad], AD[iad], int(iz*NzPoints/NzBins));
				h_Edelayed_sub_rate_z_points[int(iz*NzPoints/NzBins)]->SetName(name);
			}

			if((iz*NzPoints)%NzBins == 0) h_Edelayed_sub_norm_z_points[int(iz*NzPoints/NzBins)]=(TH1F*)h_Edelayed_sub_norm_z[iz]->Clone();
			else h_Edelayed_sub_norm_z_points[int(iz*NzPoints/NzBins)]->Add(h_Edelayed_sub_norm_z[iz]);

			if((iz*NzPoints)%NzBins == 0){
				sprintf(name, "h_Edelayed_sub_norm_z_points_eh%dad%d_iz%d", EH[iad], AD[iad], int(iz*NzPoints/NzBins));
				h_Edelayed_sub_norm_z_points[int(iz*NzPoints/NzBins)]->SetName(name);
			}

			if((iz*NzPoints)%NzBins == 0) h_Edelayed_sub_DTnorm_z_points[int(iz*NzPoints/NzBins)]=(TH1F*)h_Edelayed_sub_DTnorm_z[iz]->Clone();
			else h_Edelayed_sub_DTnorm_z_points[int(iz*NzPoints/NzBins)]->Add(h_Edelayed_sub_DTnorm_z[iz]);

			if((iz*NzPoints)%NzBins == 0){
				sprintf(name, "h_Edelayed_sub_DTnorm_z_points_eh%dad%d_iz%d", EH[iad], AD[iad], int(iz*NzPoints/NzBins));
				h_Edelayed_sub_DTnorm_z_points[int(iz*NzPoints/NzBins)]->SetName(name);
			}
		}

			//Getting histograms of the r2-dependence
		TH1F* h_Edelayed_ibd_r2[Nr2Bins];
		TH1F* h_Edelayed_sub_rate_r2[Nr2Bins];
		TH1F* h_Edelayed_ibd_r2_points[Nr2Points];
		TH1F* h_Edelayed_sub_rate_r2_points[Nr2Points];
		TH1F* h_Edelayed_sub_norm_r2[Nr2Bins];
		TH1F* h_Edelayed_sub_norm_r2_points[Nr2Points];
		TH1F* h_Edelayed_sub_DTnorm_r2[Nr2Bins];
		TH1F* h_Edelayed_sub_DTnorm_r2_points[Nr2Points];
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			sprintf(name, "h_Edelayed_ibd_r2_ad%d_ir2%d", AD[iad], ir2+1);
			if(DT == 1) sprintf(name, "h_Edelayed_ibd_DT800_r2_ad%d_ir2%d", AD[iad], ir2+1);
			h_Edelayed_ibd_r2[ir2] = (TH1F*)subFile->Get(name);

			sprintf(name, "h_Edelayed_subtract_rate_r2_ad%d_ir2%d", AD[iad], ir2+1);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_rate_DT800_r2_ad%d_ir2%d", AD[iad], ir2+1);
			h_Edelayed_sub_rate_r2[ir2] = (TH1F*)subFile->Get(name);

			sprintf(name, "h_Edelayed_subtract_norm_r2_ad%d_ir2%d", AD[iad], ir2+1);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_norm_DT800_r2_ad%d_ir2%d", AD[iad], ir2+1);
			h_Edelayed_sub_norm_r2[ir2] = (TH1F*)subFile->Get(name);

			sprintf(name, "h_Edelayed_subtract_DTnorm_r2_ad%d_ir2%d", AD[iad], ir2+1);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_DTnorm_DT800_r2_ad%d_ir2%d", AD[iad], ir2+1);
			h_Edelayed_sub_DTnorm_r2[ir2] = (TH1F*)subFile->Get(name);

			for(int iBin = 0; iBin < (h_Edelayed_ibd_r2[ir2]->GetNbinsX())+1; iBin++){
				h_Edelayed_ibd_r2[ir2]->SetBinError(iBin, sqrt(h_Edelayed_ibd_r2[ir2]->GetBinContent(iBin)));
				h_Edelayed_sub_rate_r2[ir2]->SetBinError(iBin, sqrt(h_Edelayed_ibd_r2[ir2]->GetBinContent(iBin)));
				h_Edelayed_sub_norm_r2[ir2]->SetBinError(iBin, sqrt(h_Edelayed_ibd_r2[ir2]->GetBinContent(iBin)));
				h_Edelayed_sub_DTnorm_r2[ir2]->SetBinError(iBin, sqrt(h_Edelayed_ibd_r2[ir2]->GetBinContent(iBin)));
			}
		}


		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			if((ir2*Nr2Points)%Nr2Bins == 0) h_Edelayed_sub_rate_r2_points[int(ir2*Nr2Points/Nr2Bins)]=(TH1F*)h_Edelayed_sub_rate_r2[ir2]->Clone();
			else h_Edelayed_sub_rate_r2_points[int(ir2*Nr2Points/Nr2Bins)]->Add(h_Edelayed_sub_rate_r2[ir2]);

			if((ir2*Nr2Points)%Nr2Bins == 0) h_Edelayed_ibd_r2_points[int(ir2*Nr2Points/Nr2Bins)]=(TH1F*)h_Edelayed_ibd_r2[ir2]->Clone();
			else h_Edelayed_ibd_r2_points[int(ir2*Nr2Points/Nr2Bins)]->Add(h_Edelayed_ibd_r2[ir2]);

			if((ir2*Nr2Points)%Nr2Bins == 0){
				sprintf(name, "h_Edelayed_sub_rate_r2_points_eh%dad%d_ir2%d", EH[iad], AD[iad], int(ir2*Nr2Points/Nr2Bins));
				h_Edelayed_sub_rate_r2_points[int(ir2*Nr2Points/Nr2Bins)]->SetName(name);
			}

			if((ir2*Nr2Points)%Nr2Bins == 0) h_Edelayed_sub_norm_r2_points[int(ir2*Nr2Points/Nr2Bins)]=(TH1F*)h_Edelayed_sub_norm_r2[ir2]->Clone();
			else h_Edelayed_sub_norm_r2_points[int(ir2*Nr2Points/Nr2Bins)]->Add(h_Edelayed_sub_norm_r2[ir2]);

			if((ir2*Nr2Points)%Nr2Bins == 0){
				sprintf(name, "h_Edelayed_sub_norm_r2_points_eh%dad%d_ir2%d", EH[iad], AD[iad], int(ir2*Nr2Points/Nr2Bins));
				h_Edelayed_sub_norm_r2_points[int(ir2*Nr2Points/Nr2Bins)]->SetName(name);
			}

			if((ir2*Nr2Points)%Nr2Bins == 0) h_Edelayed_sub_DTnorm_r2_points[int(ir2*Nr2Points/Nr2Bins)]=(TH1F*)h_Edelayed_sub_DTnorm_r2[ir2]->Clone();
			else h_Edelayed_sub_DTnorm_r2_points[int(ir2*Nr2Points/Nr2Bins)]->Add(h_Edelayed_sub_DTnorm_r2[ir2]);

			if((ir2*Nr2Points)%Nr2Bins == 0){
				sprintf(name, "h_Edelayed_sub_DTnorm_r2_points_eh%dad%d_ir2%d", EH[iad], AD[iad], int(ir2*Nr2Points/Nr2Bins));
				h_Edelayed_sub_DTnorm_r2_points[int(ir2*Nr2Points/Nr2Bins)]->SetName(name);
			}
		}

			//Getting histograms of the zVSr2-dependence
		TH1F* h_Edelayed_ibd_zVSr2[Nr2Bins][NzBins];
		TH1F* h_Edelayed_sub_rate_zVSr2[Nr2Bins][NzBins];
		TH1F* h_Edelayed_ibd_zVSr2_points[Nr2Points][NzPoints];
		TH1F* h_Edelayed_sub_rate_zVSr2_points[Nr2Points][NzPoints];
		TH1F* h_Edelayed_sub_norm_zVSr2[Nr2Bins][NzBins];
		TH1F* h_Edelayed_sub_norm_zVSr2_points[Nr2Points][NzPoints];
		TH1F* h_Edelayed_sub_DTnorm_zVSr2[Nr2Bins][NzBins];
		TH1F* h_Edelayed_sub_DTnorm_zVSr2_points[Nr2Points][NzPoints];
		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			for(int iz = 0; iz < NzBins; iz++){
				sprintf(name, "h_Edelayed_ibd_zVSr2_ad%d_ir2%d_iz%d", AD[iad], ir2+1, iz+1);
				if(DT == 1) sprintf(name, "h_Edelayed_ibd_DT800_zVSr2_ad%d_ir2%d_iz%d", AD[iad], ir2+1, iz+1);
				h_Edelayed_ibd_zVSr2[ir2][iz] = (TH1F*)subFile->Get(name);

				sprintf(name, "h_Edelayed_subtract_rate_zVSr2_ad%d_ir2%d_iz%d", AD[iad], ir2+1, iz+1);
				if(DT == 1) sprintf(name, "h_Edelayed_subtract_rate_DT800_zVSr2_ad%d_ir2%d_iz%d", AD[iad], ir2+1, iz+1);
				h_Edelayed_sub_rate_zVSr2[ir2][iz] = (TH1F*)subFile->Get(name);

				sprintf(name, "h_Edelayed_subtract_norm_zVSr2_ad%d_ir2%d_iz%d", AD[iad], ir2+1, iz+1);
				if(DT == 1) sprintf(name, "h_Edelayed_subtract_norm_DT800_zVSr2_ad%d_ir2%d_iz%d", AD[iad], ir2+1, iz+1);
				h_Edelayed_sub_norm_zVSr2[ir2][iz] = (TH1F*)subFile->Get(name);

				sprintf(name, "h_Edelayed_subtract_DTnorm_zVSr2_ad%d_ir2%d_iz%d", AD[iad], ir2+1, iz+1);
				if(DT == 1) sprintf(name, "h_Edelayed_subtract_DTnorm_DT800_zVSr2_ad%d_ir2%d_iz%d", AD[iad], ir2+1, iz+1);
				h_Edelayed_sub_DTnorm_zVSr2[ir2][iz] = (TH1F*)subFile->Get(name);

				for(int iBin = 0; iBin < (h_Edelayed_ibd_zVSr2[ir2][iz]->GetNbinsX())+1; iBin++){
					h_Edelayed_ibd_zVSr2[ir2][iz]->SetBinError(iBin, sqrt(h_Edelayed_ibd_zVSr2[ir2][iz]->GetBinContent(iBin)));
					h_Edelayed_sub_rate_zVSr2[ir2][iz]->SetBinError(iBin, sqrt(h_Edelayed_ibd_zVSr2[ir2][iz]->GetBinContent(iBin)));
					h_Edelayed_sub_norm_zVSr2[ir2][iz]->SetBinError(iBin, sqrt(h_Edelayed_ibd_zVSr2[ir2][iz]->GetBinContent(iBin)));
					h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]->SetBinError(iBin, sqrt(h_Edelayed_ibd_zVSr2[ir2][iz]->GetBinContent(iBin)));
				}
			}
		}

		for(int ir2 = 0; ir2 < Nr2Bins; ir2++){
			for(int iz = 0; iz < NzBins; iz++){
				if(((iz*NzPoints)%NzBins == 0) && ((ir2*Nr2Points)%Nr2Bins == 0)) h_Edelayed_ibd_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]=(TH1F*)h_Edelayed_ibd_zVSr2[ir2][iz]->Clone();
				else h_Edelayed_ibd_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]->Add(h_Edelayed_ibd_zVSr2[ir2][iz]);

				if(((iz*NzPoints)%NzBins == 0) && ((ir2*Nr2Points)%Nr2Bins == 0)) h_Edelayed_sub_rate_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]=(TH1F*)h_Edelayed_sub_rate_zVSr2[ir2][iz]->Clone();
				else h_Edelayed_sub_rate_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]->Add(h_Edelayed_sub_rate_zVSr2[ir2][iz]);

				if(((iz*NzPoints)%NzBins == 0) && ((ir2*Nr2Points)%Nr2Bins == 0)){
					sprintf(name, "h_Edelayed_sub_rate_zVSr2_points_eh%dad%d_ir2%d_iz%d", EH[iad], AD[iad], int(ir2*Nr2Points/Nr2Bins), int(iz*NzPoints/NzBins));
					h_Edelayed_sub_rate_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]->SetName(name);
				}

				if(((iz*NzPoints)%NzBins == 0) && ((ir2*Nr2Points)%Nr2Bins == 0)) h_Edelayed_sub_norm_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]=(TH1F*)h_Edelayed_sub_norm_zVSr2[ir2][iz]->Clone();
				else h_Edelayed_sub_norm_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]->Add(h_Edelayed_sub_norm_zVSr2[ir2][iz]);

				if(((iz*NzPoints)%NzBins == 0) && ((ir2*Nr2Points)%Nr2Bins == 0)){
					sprintf(name, "h_Edelayed_sub_norm_zVSr2_points_eh%dad%d_ir2%d_iz%d", EH[iad], AD[iad], int(ir2*Nr2Points/Nr2Bins), int(iz*NzPoints/NzBins));
					h_Edelayed_sub_norm_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]->SetName(name);
				}

				if(((iz*NzPoints)%NzBins == 0) && ((ir2*Nr2Points)%Nr2Bins == 0)) h_Edelayed_sub_DTnorm_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]=(TH1F*)h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]->Clone();
				else h_Edelayed_sub_DTnorm_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]->Add(h_Edelayed_sub_DTnorm_zVSr2[ir2][iz]);

				if(((iz*NzPoints)%NzBins == 0) && ((ir2*Nr2Points)%Nr2Bins == 0)){
					sprintf(name, "h_Edelayed_sub_DTnorm_zVSr2_points_eh%dad%d_ir2%d_iz%d", EH[iad], AD[iad], int(ir2*Nr2Points/Nr2Bins), int(iz*NzPoints/NzBins));
					h_Edelayed_sub_DTnorm_zVSr2_points[int(ir2*Nr2Points/Nr2Bins)][int(iz*NzPoints/NzBins)]->SetName(name);
				}
			}
		}
*/
			//renaming
			sprintf(name, "h_Edelayed_subtract_rate_eh%dad%d", EH[iad], AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[iad], AD[iad]);
			delayed_rate->SetName(name);

			sprintf(name, "h_Edelayed_subtract_norm_eh%dad%d", EH[iad], AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_norm_eh%dad%d", EH[iad], AD[iad]);
			delayed_norm->SetName(name);

			sprintf(name, "h_Edelayed_subtract_DTnorm_eh%dad%d", EH[iad], AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_DT800_DTnorm_eh%dad%d", EH[iad], AD[iad]);
			delayed_DTnorm->SetName(name);

			TH1F* delayed_rate_Sam = (TH1F*)delayed_rate->Clone();
			TH1F* delayed_norm_Sam = (TH1F*)delayed_norm->Clone();
			TH1F* delayed_DTnorm_Sam = (TH1F*)delayed_DTnorm->Clone();

			sprintf(name, "h_Edelayed_subtract_Sam_rate_eh%dad%d", EH[iad], AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_Sam_DT800_rate_eh%dad%d", EH[iad], AD[iad]);
			delayed_rate_Sam->SetName(name);

			sprintf(name, "h_Edelayed_subtract_Sam_norm_eh%dad%d", EH[iad], AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_Sam_DT800_norm_eh%dad%d", EH[iad], AD[iad]);
			delayed_norm_Sam->SetName(name);

			sprintf(name, "h_Edelayed_subtract_Sam_DTnorm_eh%dad%d", EH[iad], AD[iad]);
			if(DT == 1) sprintf(name, "h_Edelayed_subtract_Sam_DT800_DTnorm_eh%dad%d", EH[iad], AD[iad]);
			delayed_DTnorm_Sam->SetName(name);



		//Fit function:
			TF1* delayedFit_rate = new TF1("delayedFit_rate", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
			TF1* delayedFit_norm = new TF1("delayedFit_norm", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
			TF1* delayedFit_DTnorm = new TF1("delayedFit_DTnorm", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);

			TF1* delayedFit_rate_Sam = new TF1("delayedFit_rate_Sam", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.5,2.8);
			TF1* delayedFit_norm_Sam = new TF1("delayedFit_norm_Sam", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.5,2.8);
			TF1* delayedFit_DTnorm_Sam = new TF1("delayedFit_DTnorm_Sam", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.5,2.8);

/*			TF1* delayedFit_rate_Sam = new TF1("delayedFit_rate_Sam", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1]) * [4]/((exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]) * exp(2*[4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
			TF1* delayedFit_norm_Sam = new TF1("delayedFit_norm_Sam", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1]) * [4]/((exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]) * exp(2*[4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
			TF1* delayedFit_DTnorm_Sam = new TF1("delayedFit_DTnorm_Sam", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1]) * [4]/((exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]) * exp(2*[4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);*/


		//Saving to output file and fitting:
			outfile->cd();

			delayed_rate->GetXaxis()->SetTitle("Energy [MeV]");
			delayed_rate->GetYaxis()->SetTitle("Counts");
			delayed_rate->GetXaxis()->SetRangeUser(1.5,3);
				delayedFit_rate->SetParameter(0,(delayed_rate->GetBinContent(delayed_rate->FindBin(2.3)))/3.); //rate
				delayedFit_rate->SetParameter(1,0.8); //alpha
				delayedFit_rate->SetParameter(2,2.3); //mu
				delayedFit_rate->SetParameter(3,0.135); //sigma
				delayedFit_rate->SetParameter(4,3); //lambda

				delayedFit_rate->SetParName(0,"N"); //normalization
				delayedFit_rate->SetParName(1,"alpha"); //alpha
				delayedFit_rate->SetParName(2,"mu"); //mu
				delayedFit_rate->SetParName(3,"sigma"); //sigma
				delayedFit_rate->SetParName(4,"lambda"); //lambda
			delayed_rate->Fit("delayedFit_rate", "R");
			delayed_rate->SetStats(11111);
			delayed_rate->Write();
			delayed_fine_rate->Write();
			delayed_rate_Sam->GetXaxis()->SetTitle("Energy [MeV]");
			delayed_rate_Sam->GetYaxis()->SetTitle("Counts");
			delayed_rate_Sam->GetXaxis()->SetRangeUser(1.5,3);
				delayedFit_rate_Sam->SetParameter(0,(delayed_rate_Sam->GetBinContent(delayed_rate_Sam->FindBin(2.3)))/3.); //rate
				delayedFit_rate_Sam->SetParameter(1,0.8); //alpha
				delayedFit_rate_Sam->SetParameter(2,2.3); //mu
				delayedFit_rate_Sam->SetParameter(3,0.136); //sigma
				delayedFit_rate_Sam->SetParameter(4,3); //lambda

				delayedFit_rate_Sam->SetParName(0,"N"); //normalization
				delayedFit_rate_Sam->SetParName(1,"alpha"); //alpha
				delayedFit_rate_Sam->SetParName(2,"mu"); //mu
				delayedFit_rate_Sam->SetParName(3,"sigma"); //sigma
				delayedFit_rate_Sam->SetParName(4,"lambda"); //lambda
			delayed_rate_Sam->Fit("delayedFit_rate_Sam", "R");
			delayed_rate_Sam->Write();

			delayed_norm->GetXaxis()->SetTitle("Energy [MeV]");
			delayed_norm->GetYaxis()->SetTitle("Counts");
			delayed_norm->GetXaxis()->SetRangeUser(1.5,3);
				delayedFit_norm->SetParameter(0,(delayed_norm->GetBinContent(delayed_norm->FindBin(2.3)))/3.); //normalization
				delayedFit_norm->SetParameter(1,0.8); //alpha
				delayedFit_norm->SetParameter(2,2.3); //mu
				delayedFit_norm->SetParameter(3,0.135); //sigma
				delayedFit_norm->SetParameter(4,3); //lambda

				delayedFit_norm->SetParName(0,"N"); //normalization
				delayedFit_norm->SetParName(1,"alpha"); //alpha
				delayedFit_norm->SetParName(2,"mu"); //mu
				delayedFit_norm->SetParName(3,"sigma"); //sigma
				delayedFit_norm->SetParName(4,"lambda"); //lambda
			delayed_norm->Fit("delayedFit_norm", "R");
			delayed_norm->Write();
			delayed_fine_norm->Write();
			delayed_norm_Sam->GetXaxis()->SetTitle("Energy [MeV]");
			delayed_norm_Sam->GetYaxis()->SetTitle("Counts");
			delayed_norm_Sam->GetXaxis()->SetRangeUser(1.5,3);
				delayedFit_norm_Sam->SetParameter(0,(delayed_norm_Sam->GetBinContent(delayed_norm_Sam->FindBin(2.3)))/3.); //normalization
				delayedFit_norm_Sam->SetParameter(1,0.8); //alpha
				delayedFit_norm_Sam->SetParameter(2,2.3); //mu
					if(iad == 6) delayedFit_norm_Sam->SetParameter(2,2.255); //mu
				delayedFit_norm_Sam->SetParameter(3,0.135); //sigma
					if(iad == 6) delayedFit_norm_Sam->SetParameter(3,0.136); //sigma
				delayedFit_norm_Sam->SetParameter(4,3); //lambda

				delayedFit_norm_Sam->SetParName(0,"N"); //normalization
				delayedFit_norm_Sam->SetParName(1,"alpha"); //alpha
				delayedFit_norm_Sam->SetParName(2,"mu"); //mu
				delayedFit_norm_Sam->SetParName(3,"sigma"); //sigma
				delayedFit_norm_Sam->SetParName(4,"lambda"); //lambda
			delayed_norm_Sam->Fit("delayedFit_norm_Sam", "R");
			delayed_norm_Sam->Write();

			delayed_DTnorm->GetXaxis()->SetTitle("Energy [MeV]");
			delayed_DTnorm->GetYaxis()->SetTitle("Counts");
			delayed_DTnorm->GetXaxis()->SetRangeUser(1.5,3);
				delayedFit_DTnorm->SetParameter(0,(delayed_DTnorm->GetBinContent(delayed_DTnorm->FindBin(2.3)))/3.); //DT normalization
				delayedFit_DTnorm->SetParameter(1,0.8); //alpha
				delayedFit_DTnorm->SetParameter(2,2.3); //mu
				delayedFit_DTnorm->SetParameter(3,0.135); //sigma
				delayedFit_DTnorm->SetParameter(4,3); //lambda

				delayedFit_DTnorm->SetParName(0,"N"); //normalization
				delayedFit_DTnorm->SetParName(1,"alpha"); //alpha
				delayedFit_DTnorm->SetParName(2,"mu"); //mu
				delayedFit_DTnorm->SetParName(3,"sigma"); //sigma
				delayedFit_DTnorm->SetParName(4,"lambda"); //lambda
			delayed_DTnorm->Fit("delayedFit_DTnorm", "R");
			delayed_DTnorm->Write();
			delayed_fine_DTnorm->Write();
			delayed_DTnorm_Sam->GetXaxis()->SetTitle("Energy [MeV]");
			delayed_DTnorm_Sam->GetYaxis()->SetTitle("Counts");
			delayed_DTnorm_Sam->GetXaxis()->SetRangeUser(1.5,3);
				delayedFit_DTnorm_Sam->SetParameter(0,(delayed_DTnorm_Sam->GetBinContent(delayed_DTnorm_Sam->FindBin(2.3)))/3.); //DT normalization
				delayedFit_DTnorm_Sam->SetParameter(1,0.8); //alpha
				delayedFit_DTnorm_Sam->SetParameter(2,2.3); //mu
				delayedFit_DTnorm_Sam->SetParameter(3,0.135); //sigma
					if(iad == 6) delayedFit_DTnorm_Sam->SetParameter(3,0.136); //sigma
				delayedFit_DTnorm_Sam->SetParameter(4,3); //lambda

				delayedFit_DTnorm_Sam->SetParName(0,"N"); //normalization
				delayedFit_DTnorm_Sam->SetParName(1,"alpha"); //alpha
				delayedFit_DTnorm_Sam->SetParName(2,"mu"); //mu
				delayedFit_DTnorm_Sam->SetParName(3,"sigma"); //sigma
				delayedFit_DTnorm_Sam->SetParName(4,"lambda"); //lambda
			delayed_DTnorm_Sam->Fit("delayedFit_DTnorm_Sam", "R");
			delayed_DTnorm_Sam->Write();

			peak_rate[iad] = delayedFit_rate->GetParameter("mu");
			peak_norm[iad] = delayedFit_norm->GetParameter("mu");
			peak_DTnorm[iad] = delayedFit_DTnorm->GetParameter("mu");
			peak_rate_Sam[iad] = delayedFit_rate_Sam->GetParameter("mu");
			peak_norm_Sam[iad] = delayedFit_norm_Sam->GetParameter("mu");
			peak_DTnorm_Sam[iad] = delayedFit_DTnorm_Sam->GetParameter("mu");
				peakError_rate[iad] = delayedFit_rate->GetParError(2);
				peakError_norm[iad] = delayedFit_norm->GetParError(2);
				peakError_DTnorm[iad] = delayedFit_DTnorm->GetParError(2);
				peakError_rate_Sam[iad] = delayedFit_rate_Sam->GetParError(2);
				peakError_norm_Sam[iad] = delayedFit_norm_Sam->GetParError(2);
				peakError_DTnorm_Sam[iad] = delayedFit_DTnorm_Sam->GetParError(2);

			sigma_rate[iad] = delayedFit_rate->GetParameter("sigma");
			sigma_norm[iad] = delayedFit_norm->GetParameter("sigma");
			sigma_DTnorm[iad] = delayedFit_DTnorm->GetParameter("sigma");
			sigma_rate_Sam[iad] = delayedFit_rate_Sam->GetParameter("sigma");
			sigma_norm_Sam[iad] = delayedFit_norm_Sam->GetParameter("sigma");
			sigma_DTnorm_Sam[iad] = delayedFit_DTnorm_Sam->GetParameter("sigma");
				sigmaError_rate[iad] = delayedFit_rate->GetParError(3);
				sigmaError_norm[iad] = delayedFit_norm->GetParError(3);
				sigmaError_DTnorm[iad] = delayedFit_DTnorm->GetParError(3);
				sigmaError_rate_Sam[iad] = delayedFit_rate_Sam->GetParError(3);
				sigmaError_norm_Sam[iad] = delayedFit_norm_Sam->GetParError(3);
				sigmaError_DTnorm_Sam[iad] = delayedFit_DTnorm_Sam->GetParError(3);

		//if(iad > 3) sigma_rate[iad] = .136; //forcing the far ADs to have the near ADs' sigma values

			alpha_rate[iad] = delayedFit_rate->GetParameter("alpha");
			alpha_norm[iad] = delayedFit_norm->GetParameter("alpha");
			alpha_DTnorm[iad] = delayedFit_DTnorm->GetParameter("alpha");
			alpha_rate_Sam[iad] = delayedFit_rate_Sam->GetParameter("alpha");
			alpha_norm_Sam[iad] = delayedFit_norm_Sam->GetParameter("alpha");
			alpha_DTnorm_Sam[iad] = delayedFit_DTnorm_Sam->GetParameter("alpha");
				alphaError_rate[iad] = delayedFit_rate->GetParError(1);
				alphaError_norm[iad] = delayedFit_norm->GetParError(1);
				alphaError_DTnorm[iad] = delayedFit_DTnorm->GetParError(1);
				alphaError_rate_Sam[iad] = delayedFit_rate_Sam->GetParError(1);
				alphaError_norm_Sam[iad] = delayedFit_norm_Sam->GetParError(1);
				alphaError_DTnorm_Sam[iad] = delayedFit_DTnorm_Sam->GetParError(1);

			lambda_rate[iad] = delayedFit_rate->GetParameter("lambda");
			lambda_norm[iad] = delayedFit_norm->GetParameter("lambda");
			lambda_DTnorm[iad] = delayedFit_DTnorm->GetParameter("lambda");
			lambda_rate_Sam[iad] = delayedFit_rate_Sam->GetParameter("lambda");
			lambda_norm_Sam[iad] = delayedFit_norm_Sam->GetParameter("lambda");
			lambda_DTnorm_Sam[iad] = delayedFit_DTnorm_Sam->GetParameter(4);
				lambdaError_rate[iad] = delayedFit_rate->GetParError(4);
				lambdaError_norm[iad] = delayedFit_norm->GetParError(4);
				lambdaError_DTnorm[iad] = delayedFit_DTnorm->GetParError(4);
				lambdaError_rate_Sam[iad] = delayedFit_rate_Sam->GetParError(4);
				lambdaError_norm_Sam[iad] = delayedFit_norm_Sam->GetParError(4);
				lambdaError_DTnorm_Sam[iad] = delayedFit_DTnorm_Sam->GetParError(4);

/*		Ehigh_rate[iad] = peak_rate[iad]+5*sigma_rate[iad];
		Elow_rate[iad] = peak_rate[iad]-5*sigma_rate[iad];
		Ehigh_norm[iad] = peak_norm[iad]+5*sigma_norm[iad];
		Elow_norm[iad] = peak_norm[iad]-5*sigma_norm[iad];
		Ehigh_DTnorm[iad] = peak_DTnorm[iad]+5*sigma_DTnorm[iad];
		Elow_DTnorm[iad] = peak_DTnorm[iad]-5*sigma_DTnorm[iad];*/

		Ehigh_rate[iad] = 2.8;
		Elow_rate[iad] = 1.5;
		Ehigh_norm[iad] = 2.8;
		Elow_norm[iad] = 1.5;
		Ehigh_DTnorm[iad] = 2.8;
		Elow_DTnorm[iad] = 1.5;

			N_nom_rate[iad] = delayed_fine_rate->Integral(delayed_fine_rate->FindBin(peak_rate[iad]-3*sigma_rate[iad]), delayed_fine_rate->FindBin(peak_rate[iad]+3*sigma_rate[iad]));
			N_ext_rate[iad] = delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad]));
			efficiency_rate[iad] = N_nom_rate[iad]/N_ext_rate[iad];

			N_nom_norm[iad] = delayed_fine_norm->Integral(delayed_fine_norm->FindBin(peak_norm[iad]-3*sigma_norm[iad]), delayed_fine_norm->FindBin(peak_norm[iad]+3*sigma_norm[iad]));
			N_ext_norm[iad] = delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad]));
			efficiency_norm[iad] = N_nom_norm[iad]/N_ext_norm[iad];

			N_nom_DTnorm[iad] = delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), delayed_fine_DTnorm->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad]));
			N_ext_DTnorm[iad] = delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad]));
			efficiency_DTnorm[iad] = N_nom_DTnorm[iad]/N_ext_DTnorm[iad];

			N_nom_rate_Sam[iad] = delayed_fine_rate->Integral(delayed_fine_rate->FindBin(peak_rate_Sam[iad]-3*sigma_rate_Sam[iad]), delayed_fine_rate->FindBin(peak_rate_Sam[iad]+3*sigma_rate_Sam[iad]));
			N_ext_rate_Sam[iad] = delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad]));
			efficiency_rate_Sam[iad] = N_nom_rate_Sam[iad]/N_ext_rate_Sam[iad];

			N_nom_norm_Sam[iad] = delayed_fine_norm->Integral(delayed_fine_norm->FindBin(peak_norm_Sam[iad]-3*sigma_norm_Sam[iad]), delayed_fine_norm->FindBin(peak_norm_Sam[iad]+3*sigma_norm_Sam[iad]));
			N_ext_norm_Sam[iad] = delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad]));
			efficiency_norm_Sam[iad] = N_nom_norm_Sam[iad]/N_ext_norm_Sam[iad];

			N_nom_DTnorm_Sam[iad] = delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(peak_DTnorm_Sam[iad]-3*sigma_DTnorm_Sam[iad]), delayed_fine_DTnorm->FindBin(peak_DTnorm_Sam[iad]+3*sigma_DTnorm_Sam[iad]));
			N_ext_DTnorm_Sam[iad] = delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad]));
			efficiency_DTnorm_Sam[iad] = N_nom_DTnorm_Sam[iad]/N_ext_DTnorm_Sam[iad];

		//Error for rate
			efficiencyError_rate[iad] = sqrt(pow((delayed_fine_rate->Integral(delayed_fine_rate->FindBin(peak_rate[iad]-3*sigma_rate[iad]), delayed_fine_rate->FindBin(peak_rate[iad]+3*sigma_rate[iad])))/pow(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad])),2),2)*(delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(Elow_rate[iad]), delayed_fine_IBD->FindBin(peak_rate[iad]-3*sigma_rate[iad]))+delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_rate[iad]-3*sigma_rate[iad]), delayed_fine_IBD->FindBin(Ehigh_rate[iad]))) + pow(1./(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad])))-(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(peak_rate[iad]-3*sigma_rate[iad]), delayed_fine_rate->FindBin(peak_rate[iad]+3*sigma_rate[iad])))/pow(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad])),2),2)*delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_rate[iad]-3*sigma_rate[iad]), delayed_fine_IBD->FindBin(peak_rate[iad]+ 3*sigma_rate[iad])));
		//Error for norm
			efficiencyError_norm[iad] = sqrt(pow((delayed_fine_norm->Integral(delayed_fine_norm->FindBin(peak_norm[iad]-3*sigma_norm[iad]), delayed_fine_norm->FindBin(peak_norm[iad]+3*sigma_norm[iad])))/pow(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad])),2),2)*(delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(Elow_norm[iad]), delayed_fine_IBD->FindBin(peak_norm[iad]-3*sigma_norm[iad]))+delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_norm[iad]-3*sigma_norm[iad]), delayed_fine_IBD->FindBin(Ehigh_norm[iad]))) + pow(1./(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad])))-(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(peak_norm[iad]-3*sigma_norm[iad]), delayed_fine_norm->FindBin(peak_norm[iad]+3*sigma_norm[iad])))/pow(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad])),2),2)*delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_norm[iad]-3*sigma_norm[iad]), delayed_fine_IBD->FindBin(peak_norm[iad]+ 3*sigma_norm[iad])));
		//Error for DTnorm
			efficiencyError_DTnorm[iad] = sqrt(pow((delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), delayed_fine_DTnorm->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad])))/pow(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad])),2),2)*(delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(Elow_DTnorm[iad]), delayed_fine_IBD->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]))+delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), delayed_fine_IBD->FindBin(Ehigh_DTnorm[iad]))) + pow(1./(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad])))-(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), delayed_fine_DTnorm->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad])))/pow(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad])),2),2)*delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), delayed_fine_IBD->FindBin(peak_DTnorm[iad]+ 3*sigma_DTnorm[iad])));
		//Error for rate_Sam
			efficiencyError_rate_Sam[iad] = sqrt(pow((delayed_fine_rate->Integral(delayed_fine_rate->FindBin(peak_rate_Sam[iad]-3*sigma_rate_Sam[iad]), delayed_fine_rate->FindBin(peak_rate_Sam[iad]+3*sigma_rate_Sam[iad])))/pow(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad])),2),2)*(delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(Elow_rate[iad]), delayed_fine_IBD->FindBin(peak_rate_Sam[iad]-3*sigma_rate_Sam[iad]))+delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_rate_Sam[iad]-3*sigma_rate_Sam[iad]), delayed_fine_IBD->FindBin(Ehigh_rate[iad]))) + pow(1./(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad])))-(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(peak_rate_Sam[iad]-3*sigma_rate_Sam[iad]), delayed_fine_rate->FindBin(peak_rate_Sam[iad]+3*sigma_rate_Sam[iad])))/pow(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad])),2),2)*delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_rate_Sam[iad]-3*sigma_rate_Sam[iad]), delayed_fine_IBD->FindBin(peak_rate_Sam[iad]+ 3*sigma_rate_Sam[iad])));
		//Error for norm_Sam
			efficiencyError_norm_Sam[iad] = sqrt(pow((delayed_fine_norm->Integral(delayed_fine_norm->FindBin(peak_norm_Sam[iad]-3*sigma_norm_Sam[iad]), delayed_fine_norm->FindBin(peak_norm_Sam[iad]+3*sigma_norm_Sam[iad])))/pow(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad])),2),2)*(delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(Elow_norm[iad]), delayed_fine_IBD->FindBin(peak_norm_Sam[iad]-3*sigma_norm_Sam[iad]))+delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_norm_Sam[iad]-3*sigma_norm_Sam[iad]), delayed_fine_IBD->FindBin(Ehigh_norm[iad]))) + pow(1./(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad])))-(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(peak_norm_Sam[iad]-3*sigma_norm_Sam[iad]), delayed_fine_norm->FindBin(peak_norm_Sam[iad]+3*sigma_norm_Sam[iad])))/pow(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad])),2),2)*delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_norm_Sam[iad]-3*sigma_norm_Sam[iad]), delayed_fine_IBD->FindBin(peak_norm_Sam[iad]+ 3*sigma_norm_Sam[iad])));
		//Error for DTnorm_Sam
			efficiencyError_DTnorm_Sam[iad] = sqrt(pow((delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(peak_DTnorm_Sam[iad]-3*sigma_DTnorm_Sam[iad]), delayed_fine_DTnorm->FindBin(peak_DTnorm_Sam[iad]+3*sigma_DTnorm_Sam[iad])))/pow(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad])),2),2)*(delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(Elow_DTnorm[iad]), delayed_fine_IBD->FindBin(peak_DTnorm_Sam[iad]-3*sigma_DTnorm_Sam[iad]))+delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_DTnorm_Sam[iad]-3*sigma_DTnorm_Sam[iad]), delayed_fine_IBD->FindBin(Ehigh_DTnorm[iad]))) + pow(1./(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad])))-(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(peak_DTnorm_Sam[iad]-3*sigma_DTnorm_Sam[iad]), delayed_fine_DTnorm->FindBin(peak_DTnorm_Sam[iad]+3*sigma_DTnorm_Sam[iad])))/pow(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad])),2),2)*delayed_fine_IBD->Integral(delayed_fine_IBD->FindBin(peak_DTnorm_Sam[iad]-3*sigma_DTnorm_Sam[iad]), delayed_fine_IBD->FindBin(peak_DTnorm_Sam[iad]+ 3*sigma_DTnorm_Sam[iad])));


//For the efficiency as a function of sigma:
		double fract_sigma = 0; //bin center is the sigma coefficient
		for(int iBin = 1; iBin < 501; iBin++){
			fract_sigma = h_efficiency_sigma_rate[iad]->GetBinCenter(iBin);
			h_efficiency_sigma_rate[iad]->SetBinContent(iBin, (delayed_fine_rate->Integral(delayed_fine_rate->FindBin(peak_rate[iad]-fract_sigma*sigma_rate[iad]), delayed_fine_rate->FindBin(peak_rate[iad]+fract_sigma*sigma_rate[iad])))/(delayed_fine_rate->Integral(delayed_fine_rate->FindBin(Elow_rate[iad]), delayed_fine_rate->FindBin(Ehigh_rate[iad]))));
			h_efficiency_sigma_norm[iad]->SetBinContent(iBin, (delayed_fine_norm->Integral(delayed_fine_norm->FindBin(peak_norm[iad]-fract_sigma*sigma_norm[iad]), delayed_fine_norm->FindBin(peak_norm[iad]+fract_sigma*sigma_norm[iad])))/(delayed_fine_norm->Integral(delayed_fine_norm->FindBin(Elow_norm[iad]), delayed_fine_norm->FindBin(Ehigh_norm[iad]))));
			h_efficiency_sigma_DTnorm[iad]->SetBinContent(iBin, (delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(peak_DTnorm[iad]-fract_sigma*sigma_DTnorm[iad]), delayed_fine_DTnorm->FindBin(peak_DTnorm[iad]+fract_sigma*sigma_DTnorm[iad])))/(delayed_fine_DTnorm->Integral(delayed_fine_DTnorm->FindBin(Elow_DTnorm[iad]), delayed_fine_DTnorm->FindBin(Ehigh_DTnorm[iad]))));
		} 

/*		//for z-dependence:
		for(int iz = 0; iz<NzPoints; iz++){
			N_nom_rate_z[iad][iz] = h_Edelayed_sub_rate_z_points[iz]->Integral(h_Edelayed_sub_rate_z_points[iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_z_points[iz]->FindBin(peak_rate[iad]+3*sigma_rate[iad]));
			N_ext_rate_z[iad][iz] = h_Edelayed_sub_rate_z_points[iz]->Integral(h_Edelayed_sub_rate_z_points[iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_z_points[iz]->FindBin(Ehigh_rate[iad]));
			efficiency_rate_z[iad][iz] = N_nom_rate_z[iad][iz]/N_ext_rate_z[iad][iz];

			efficiencyError_rate_z[iad][iz] = sqrt(pow((h_Edelayed_sub_rate_z_points[iz]->Integral(h_Edelayed_sub_rate_z_points[iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_z_points[iz]->FindBin(peak_rate[iad]+3*sigma_rate[iad])))/pow(h_Edelayed_sub_rate_z_points[iz]->Integral(h_Edelayed_sub_rate_z_points[iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_z_points[iz]->FindBin(Ehigh_rate[iad])),2),2)*(h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(Elow_rate[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]))+h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(Ehigh_rate[iad]))) + pow(1./(h_Edelayed_sub_rate_z_points[iz]->Integral(h_Edelayed_sub_rate_z_points[iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_z_points[iz]->FindBin(Ehigh_rate[iad])))-(h_Edelayed_sub_rate_z_points[iz]->Integral(h_Edelayed_sub_rate_z_points[iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_z_points[iz]->FindBin(peak_rate[iad]+3*sigma_rate[iad])))/pow(h_Edelayed_sub_rate_z_points[iz]->Integral(h_Edelayed_sub_rate_z_points[iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_z_points[iz]->FindBin(Ehigh_rate[iad])),2),2)*h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(peak_rate[iad]+ 3*sigma_rate[iad])));

			N_nom_norm_z[iad][iz] = h_Edelayed_sub_norm_z_points[iz]->Integral(h_Edelayed_sub_norm_z_points[iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_z_points[iz]->FindBin(peak_norm[iad]+3*sigma_norm[iad]));
			N_ext_norm_z[iad][iz] = h_Edelayed_sub_norm_z_points[iz]->Integral(h_Edelayed_sub_norm_z_points[iz]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_z_points[iz]->FindBin(Ehigh_norm[iad]));
			efficiency_norm_z[iad][iz] = N_nom_norm_z[iad][iz]/N_ext_norm_z[iad][iz];

			efficiencyError_norm_z[iad][iz] = sqrt(pow((h_Edelayed_sub_norm_z_points[iz]->Integral(h_Edelayed_sub_norm_z_points[iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_z_points[iz]->FindBin(peak_norm[iad]+3*sigma_norm[iad])))/pow(h_Edelayed_sub_norm_z_points[iz]->Integral(h_Edelayed_sub_norm_z_points[iz]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_z_points[iz]->FindBin(Ehigh_norm[iad])),2),2)*(h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(Elow_norm[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]))+h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(Ehigh_norm[iad]))) + pow(1./(h_Edelayed_sub_norm_z_points[iz]->Integral(h_Edelayed_sub_norm_z_points[iz]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_z_points[iz]->FindBin(Ehigh_norm[iad])))-(h_Edelayed_sub_norm_z_points[iz]->Integral(h_Edelayed_sub_norm_z_points[iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_z_points[iz]->FindBin(peak_norm[iad]+3*sigma_norm[iad])))/pow(h_Edelayed_sub_norm_z_points[iz]->Integral(h_Edelayed_sub_norm_z_points[iz]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_z_points[iz]->FindBin(Ehigh_norm[iad])),2),2)*h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(peak_norm[iad]+ 3*sigma_norm[iad])));

			N_nom_DTnorm_z[iad][iz] = h_Edelayed_sub_DTnorm_z_points[iz]->Integral(h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad]));
			N_ext_DTnorm_z[iad][iz] = h_Edelayed_sub_DTnorm_z_points[iz]->Integral(h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(Ehigh_DTnorm[iad]));
			efficiency_DTnorm_z[iad][iz] = N_nom_DTnorm_z[iad][iz]/N_ext_DTnorm_z[iad][iz];

			efficiencyError_DTnorm_z[iad][iz] = sqrt(pow((h_Edelayed_sub_DTnorm_z_points[iz]->Integral(h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad])))/pow(h_Edelayed_sub_DTnorm_z_points[iz]->Integral(h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(Ehigh_DTnorm[iad])),2),2)*(h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]))+h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(Ehigh_DTnorm[iad]))) + pow(1./(h_Edelayed_sub_DTnorm_z_points[iz]->Integral(h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(Ehigh_DTnorm[iad])))-(h_Edelayed_sub_DTnorm_z_points[iz]->Integral(h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad])))/pow(h_Edelayed_sub_DTnorm_z_points[iz]->Integral(h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_z_points[iz]->FindBin(Ehigh_DTnorm[iad])),2),2)*h_Edelayed_ibd_z_points[iz]->Integral(h_Edelayed_ibd_z_points[iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_ibd_z_points[iz]->FindBin(peak_DTnorm[iad]+ 3*sigma_DTnorm[iad])));
		}



		//for r2-dependence:
		for(int ir2 = 0; ir2<Nr2Points; ir2++){
			N_nom_rate_r2[iad][ir2] = h_Edelayed_sub_rate_r2_points[ir2]->Integral(h_Edelayed_sub_rate_r2_points[ir2]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_r2_points[ir2]->FindBin(peak_rate[iad]+3*sigma_rate[iad]));
			N_ext_rate_r2[iad][ir2] = h_Edelayed_sub_rate_r2_points[ir2]->Integral(h_Edelayed_sub_rate_r2_points[ir2]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_r2_points[ir2]->FindBin(Ehigh_rate[iad]));
			efficiency_rate_r2[iad][ir2] = N_nom_rate_r2[iad][ir2]/N_ext_rate_r2[iad][ir2];

			efficiencyError_rate_r2[iad][ir2] = sqrt(pow((h_Edelayed_sub_rate_r2_points[ir2]->Integral(h_Edelayed_sub_rate_r2_points[ir2]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_r2_points[ir2]->FindBin(peak_rate[iad]+3*sigma_rate[iad])))/pow(h_Edelayed_sub_rate_r2_points[ir2]->Integral(h_Edelayed_sub_rate_r2_points[ir2]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_r2_points[ir2]->FindBin(Ehigh_rate[iad])),2),2)*(h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(Elow_rate[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_rate[iad]-3*sigma_rate[iad]))+h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(Ehigh_rate[iad]))) + pow(1./(h_Edelayed_sub_rate_r2_points[ir2]->Integral(h_Edelayed_sub_rate_r2_points[ir2]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_r2_points[ir2]->FindBin(Ehigh_rate[iad])))-(h_Edelayed_sub_rate_r2_points[ir2]->Integral(h_Edelayed_sub_rate_r2_points[ir2]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_r2_points[ir2]->FindBin(peak_rate[iad]+3*sigma_rate[iad])))/pow(h_Edelayed_sub_rate_r2_points[ir2]->Integral(h_Edelayed_sub_rate_r2_points[ir2]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_r2_points[ir2]->FindBin(Ehigh_rate[iad])),2),2)*h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_rate[iad]+ 3*sigma_rate[iad])));

			N_nom_norm_r2[iad][ir2] = h_Edelayed_sub_norm_r2_points[ir2]->Integral(h_Edelayed_sub_norm_r2_points[ir2]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_r2_points[ir2]->FindBin(peak_norm[iad]+3*sigma_norm[iad]));
			N_ext_norm_r2[iad][ir2] = h_Edelayed_sub_norm_r2_points[ir2]->Integral(h_Edelayed_sub_norm_r2_points[ir2]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_r2_points[ir2]->FindBin(Ehigh_norm[iad]));
			efficiency_norm_r2[iad][ir2] = N_nom_norm_r2[iad][ir2]/N_ext_norm_r2[iad][ir2];

			efficiencyError_norm_r2[iad][ir2] = sqrt(pow((h_Edelayed_sub_norm_r2_points[ir2]->Integral(h_Edelayed_sub_norm_r2_points[ir2]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_r2_points[ir2]->FindBin(peak_norm[iad]+3*sigma_norm[iad])))/pow(h_Edelayed_sub_norm_r2_points[ir2]->Integral(h_Edelayed_sub_norm_r2_points[ir2]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_r2_points[ir2]->FindBin(Ehigh_norm[iad])),2),2)*(h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(Elow_norm[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_norm[iad]-3*sigma_norm[iad]))+h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(Ehigh_norm[iad]))) + pow(1./(h_Edelayed_sub_norm_r2_points[ir2]->Integral(h_Edelayed_sub_norm_r2_points[ir2]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_r2_points[ir2]->FindBin(Ehigh_norm[iad])))-(h_Edelayed_sub_norm_r2_points[ir2]->Integral(h_Edelayed_sub_norm_r2_points[ir2]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_r2_points[ir2]->FindBin(peak_norm[iad]+3*sigma_norm[iad])))/pow(h_Edelayed_sub_norm_r2_points[ir2]->Integral(h_Edelayed_sub_norm_r2_points[ir2]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_r2_points[ir2]->FindBin(Ehigh_norm[iad])),2),2)*h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_norm[iad]+ 3*sigma_norm[iad])));

			N_nom_DTnorm_r2[iad][ir2] = h_Edelayed_sub_DTnorm_r2_points[ir2]->Integral(h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad]));
			N_ext_DTnorm_r2[iad][ir2] = h_Edelayed_sub_DTnorm_r2_points[ir2]->Integral(h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(Ehigh_DTnorm[iad]));
			efficiency_DTnorm_r2[iad][ir2] = N_nom_DTnorm_r2[iad][ir2]/N_ext_DTnorm_r2[iad][ir2];

			efficiencyError_DTnorm_r2[iad][ir2] = sqrt(pow((h_Edelayed_sub_DTnorm_r2_points[ir2]->Integral(h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad])))/pow(h_Edelayed_sub_DTnorm_r2_points[ir2]->Integral(h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(Ehigh_DTnorm[iad])),2),2)*(h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(Elow_DTnorm[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]))+h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(Ehigh_DTnorm[iad]))) + pow(1./(h_Edelayed_sub_DTnorm_r2_points[ir2]->Integral(h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(Ehigh_DTnorm[iad])))-(h_Edelayed_sub_DTnorm_r2_points[ir2]->Integral(h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad])))/pow(h_Edelayed_sub_DTnorm_r2_points[ir2]->Integral(h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_r2_points[ir2]->FindBin(Ehigh_DTnorm[iad])),2),2)*h_Edelayed_ibd_r2_points[ir2]->Integral(h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_ibd_r2_points[ir2]->FindBin(peak_DTnorm[iad]+ 3*sigma_DTnorm[iad])));
		}

		//for zVSr2-dependence:
		for(int ir2 = 0; ir2<Nr2Points; ir2++){
			for(int iz = 0; iz<NzPoints; iz++){
				N_nom_rate_zVSr2[iad][ir2][iz] = h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]+3*sigma_rate[iad]));
				N_ext_rate_zVSr2[iad][ir2][iz] = h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Ehigh_rate[iad]));
				efficiency_rate_zVSr2[iad][ir2][iz] = N_nom_rate_zVSr2[iad][ir2][iz]/N_ext_rate_zVSr2[iad][ir2][iz];

				efficiencyError_rate_zVSr2[iad][ir2][iz] = sqrt(pow((h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]+3*sigma_rate[iad])))/pow(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Ehigh_rate[iad])),2),2)*(h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(Elow_rate[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]))+h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(Ehigh_rate[iad]))) + pow(1./(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Ehigh_rate[iad])))-(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]+3*sigma_rate[iad])))/pow(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Ehigh_rate[iad])),2),2)*h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]+ 3*sigma_rate[iad])));

		//		double N_outer = h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Elow_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]-3*sigma_rate[iad])) + h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(peak_rate[iad]+3*sigma_rate[iad]), h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->FindBin(Ehigh_rate[iad]));
		//		efficiencyError_rate_zVSr2[iad][ir2][iz] = pow(efficiency_rate_zVSr2[iad][ir2][iz],2)*N_outer/N_nom_rate_zVSr2[iad][ir2][iz]*sqrt(pow(sqrt(N_outer)/N_outer,2)+pow(sqrt(N_nom_rate_zVSr2[iad][ir2][iz])/N_nom_rate_zVSr2[iad][ir2][iz],2));

		//		N_nom_norm_zVSr2[iad][ir2][iz] = h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]+3*sigma_norm[iad]));
		//		N_ext_norm_zVSr2[iad][ir2][iz] = h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(Ehigh_norm[iad]));
		//		efficiency_norm_zVSr2[iad][ir2][iz] = N_nom_norm_zVSr2[iad][ir2][iz]/N_ext_norm_zVSr2[iad][ir2][iz];

				efficiencyError_norm_zVSr2[iad][ir2][iz] = sqrt(pow((h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]+3*sigma_norm[iad])))/pow(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(Ehigh_norm[iad])),2),2)*(h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(Elow_norm[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]))+h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(Ehigh_norm[iad]))) + pow(1./(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(Ehigh_norm[iad])))-(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]+3*sigma_norm[iad])))/pow(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(Elow_norm[iad]), h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->FindBin(Ehigh_norm[iad])),2),2)*h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]-3*sigma_norm[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_norm[iad]+ 3*sigma_norm[iad])));

				N_nom_DTnorm_zVSr2[iad][ir2][iz] = h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad]));
				N_ext_DTnorm_zVSr2[iad][ir2][iz] = h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(Ehigh_DTnorm[iad]));
				efficiency_DTnorm_zVSr2[iad][ir2][iz] = N_nom_DTnorm_zVSr2[iad][ir2][iz]/N_ext_DTnorm_zVSr2[iad][ir2][iz];

				efficiencyError_DTnorm_zVSr2[iad][ir2][iz] = sqrt(pow((h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad])))/pow(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(Ehigh_DTnorm[iad])),2),2)*(h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]))+h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(Ehigh_DTnorm[iad]))) + pow(1./(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(Ehigh_DTnorm[iad])))-(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]+3*sigma_DTnorm[iad])))/pow(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->Integral(h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(Elow_DTnorm[iad]), h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->FindBin(Ehigh_DTnorm[iad])),2),2)*h_Edelayed_ibd_zVSr2_points[ir2][iz]->Integral(h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]-3*sigma_DTnorm[iad]), h_Edelayed_ibd_zVSr2_points[ir2][iz]->FindBin(peak_DTnorm[iad]+ 3*sigma_DTnorm[iad])));
			}
		}*/

		delayed_IBD->SetName(Form("h_Edelayed_ibd_eh%dad%d",EH[iad],AD[iad]));
		delayed_IBD->Write();
		delayed_fine_IBD->Write();

	int colors_zSlices[8] = {kBlack, kRed, kGreen, kBlue, kMagenta, kOrange, kCyan, kGray+1};
/*		for(int iz = 0; iz < NzPoints; iz++){
			h_Edelayed_sub_rate_z_points[iz]->SetStats(0);
			h_Edelayed_sub_rate_z_points[iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_rate_z_points[iz]->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_rate_z_points[iz]->SetLineColor(colors_zSlices[iz]);
			h_Edelayed_sub_rate_z_points[iz]->SetMarkerColor(colors_zSlices[iz]);
	//		h_Edelayed_sub_rate_z_points[iz]->Write();
		}

		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			h_Edelayed_sub_rate_r2_points[ir2]->SetStats(0);
			h_Edelayed_sub_rate_r2_points[ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_rate_r2_points[ir2]->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_rate_r2_points[ir2]->SetLineColor(colors_zSlices[ir2]);
			h_Edelayed_sub_rate_r2_points[ir2]->SetMarkerColor(colors_zSlices[ir2]);
	//		h_Edelayed_sub_rate_r2_points[ir2]->Write();
		}

		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			for(int iz = 0; iz < NzPoints; iz++){
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->SetStats(0);
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->SetLineColor(colors_zSlices[iz]);
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->SetMarkerColor(colors_zSlices[iz]);
	//			h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Write();

				h_Edelayed_ibd_zVSr2_points[ir2][iz]->SetName(Form("h_Edelayed_ibd_zVSr2_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2,iz));
				h_Edelayed_ibd_zVSr2_points[ir2][iz]->Write();
			}	
		}

//	TCanvas *z_slice_Ed_perAD = new TCanvas(Form("z_slice_Ed_eh%dad%d",EH[iad],AD[iad]),Form("z_slice_Ed_eh%dad%d",EH[iad],AD[iad]));
	z_slice_Ed_perAD->cd(iad);
		h_Edelayed_sub_rate_z_points[2]->GetXaxis()->SetRangeUser(1.5,3);
		h_Edelayed_sub_rate_z_points[2]->Draw();
		for(int iz = 0; iz < NzPoints; iz++){
			if(iz == 2) continue;
			h_Edelayed_sub_rate_z_points[iz]->Draw("same");
		}
		TLine *linLOW = new TLine(peak_rate[iad] - 3*sigma_rate[iad],-50,peak_rate[iad] - 3*sigma_rate[iad],500);
		TLine *linHIGH = new TLine(peak_rate[iad] + 3*sigma_rate[iad],-50,peak_rate[iad] + 3*sigma_rate[iad],500);
		TLine *linPEAK = new TLine(peak_rate[iad],-50,peak_rate[iad],10000);
		linLOW->SetLineStyle(2);
		linHIGH->SetLineStyle(2);
		linPEAK->SetLineStyle(2);
		linLOW->SetLineWidth(2);
		linHIGH->SetLineWidth(2);
		linPEAK->SetLineWidth(2);
		linPEAK->SetLineColor(kCyan);
		linLOW->Draw("same");
		linHIGH->Draw("same");
		linPEAK->Draw("same");
		z_slice_Ed_perAD->SetGridy();
		//z_slice_Ed_perAD->BuildLegend();
		z_slice_Ed_perAD->Print(Form("../nH_files/z_slices_EH%dAD%d.png",EH[iad],AD[iad]));

//	TCanvas *r2_slice_Ed_perAD = new TCanvas(Form("r2_slice_Ed_eh%dad%d",EH[iad],AD[iad]),Form("r2_slice_Ed_eh%dad%d",EH[iad],AD[iad]));
	r2_slice_Ed_perAD->cd(iad);
		h_Edelayed_sub_rate_r2_points[2]->GetXaxis()->SetRangeUser(1.5,3);
		h_Edelayed_sub_rate_r2_points[2]->Draw();
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			if(ir2 == 2) continue;
			h_Edelayed_sub_rate_r2_points[ir2]->Draw("same");
		}
		TLine *linLOWr = new TLine(peak_rate[iad] - 3*sigma_rate[iad],-50,peak_rate[iad] - 3*sigma_rate[iad],500);
		TLine *linHIGHr = new TLine(peak_rate[iad] + 3*sigma_rate[iad],-50,peak_rate[iad] + 3*sigma_rate[iad],500);
		TLine *linPEAKr = new TLine(peak_rate[iad],-50,peak_rate[iad],10000);
		linLOWr->SetLineStyle(2);
		linHIGHr->SetLineStyle(2);
		linPEAKr->SetLineStyle(2);
		linLOWr->SetLineWidth(2);
		linHIGHr->SetLineWidth(2);
		linPEAKr->SetLineWidth(2);
		linPEAKr->SetLineColor(kCyan);
		linLOWr->Draw("same");
		linHIGHr->Draw("same");
		linPEAKr->Draw("same");
		r2_slice_Ed_perAD->SetGridy();
		//r2_slice_Ed_perAD->BuildLegend();
		r2_slice_Ed_perAD->Print(Form("../nH_files/r2_slices_EH%dAD%d.png",EH[iad],AD[iad]));

	TCanvas *zVSr2_slice_Ed_perAD = new TCanvas(Form("zVSr2_slice_Ed_eh%dad%d",EH[iad],AD[iad]),Form("zVSr2_slice_Ed_eh%dad%d",EH[iad],AD[iad]));
	zVSr2_slice_Ed_perAD->cd();  //consider various marker styles
		h_Edelayed_sub_rate_zVSr2_points[2][2]->GetXaxis()->SetRangeUser(1.5,3);
		h_Edelayed_sub_rate_zVSr2_points[2][2]->Draw();
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			for(int iz = 0; iz < NzPoints; iz++){
				if(ir2 == 2 && iz == 2) continue;
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Draw("same");
			}
		}
		TLine *linLOWzVSr2 = new TLine(peak_rate[iad] - 3*sigma_rate[iad],-50,peak_rate[iad] - 3*sigma_rate[iad],500);
		TLine *linHIGHzVSr2 = new TLine(peak_rate[iad] + 3*sigma_rate[iad],-50,peak_rate[iad] + 3*sigma_rate[iad],500);
		TLine *linPEAKzVSr2 = new TLine(peak_rate[iad],-50,peak_rate[iad],10000);
		linLOWzVSr2->SetLineStyle(2);
		linHIGHzVSr2->SetLineStyle(2);
		linPEAKzVSr2->SetLineStyle(2);
		linLOWzVSr2->SetLineWidth(2);
		linHIGHzVSr2->SetLineWidth(2);
		linPEAKzVSr2->SetLineWidth(2);
		linPEAKzVSr2->SetLineColor(kCyan);
		linLOWzVSr2->Draw("same");
		linHIGHzVSr2->Draw("same");
		linPEAKzVSr2->Draw("same");
		zVSr2_slice_Ed_perAD->SetGridy();
		//zVSr2_slice_Ed_perAD->BuildLegend();
		zVSr2_slice_Ed_perAD->Print(Form("../nH_files/zVSr2_slices_EH%dAD%d.png",EH[iad],AD[iad]));
*/
		outfile->cd();
/*		for(int iz = 0; iz < NzPoints; iz++){
			h_Edelayed_sub_rate_z_points[iz]->SetStats(0);
			h_Edelayed_sub_rate_z_points[iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_rate_z_points[iz]->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_rate_z_points[iz]->Write();

			h_Edelayed_sub_norm_z_points[iz]->SetStats(0);
			h_Edelayed_sub_norm_z_points[iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_norm_z_points[iz]->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_norm_z_points[iz]->Write();

			h_Edelayed_sub_DTnorm_z_points[iz]->SetStats(0);
			h_Edelayed_sub_DTnorm_z_points[iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_DTnorm_z_points[iz]->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_DTnorm_z_points[iz]->Write();
		}

		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			h_Edelayed_sub_rate_r2_points[ir2]->SetStats(0);
			h_Edelayed_sub_rate_r2_points[ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_rate_r2_points[ir2]->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_rate_r2_points[ir2]->Write();

			h_Edelayed_sub_norm_r2_points[ir2]->SetStats(0);
			h_Edelayed_sub_norm_r2_points[ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_norm_r2_points[ir2]->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_norm_r2_points[ir2]->Write();

			h_Edelayed_sub_DTnorm_r2_points[ir2]->SetStats(0);
			h_Edelayed_sub_DTnorm_r2_points[ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_DTnorm_r2_points[ir2]->GetYaxis()->SetTitle("Counts");
			h_Edelayed_sub_DTnorm_r2_points[ir2]->Write();
		}

		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			for(int iz = 0; iz < NzPoints; iz++){
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->SetStats(0);
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_rate_zVSr2_points[ir2][iz]->Write();

				h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->SetStats(0);
				h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_norm_zVSr2_points[ir2][iz]->Write();

				h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->SetStats(0);
				h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
				h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->GetYaxis()->SetTitle("Counts");
				h_Edelayed_sub_DTnorm_zVSr2_points[ir2][iz]->Write();
			}
		}*/
	//	subFile->Close();

	}


/*	for(int iad = 0; iad < 8; iad++){
		for(int iz = 0; iz < NzPoints; iz++){
			h_Edelayed_sub_rate_z[iad][iz] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_rate_z_points_eh%dad%d_iz%d", EH[iad], AD[iad],iz));
			h_Edelayed_sub_norm_z[iad][iz] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_norm_z_points_eh%dad%d_iz%d", EH[iad], AD[iad],iz));
			h_Edelayed_sub_DTnorm_z[iad][iz] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_DTnorm_z_points_eh%dad%d_iz%d", EH[iad], AD[iad],iz));

			h_Edelayed_sub_rate_z[iad][iz]->SetTitle(Form("h_Edelayed_sub_rate_z_eh%dad%d_iz%d",EH[iad],AD[iad],iz+1));
			h_Edelayed_sub_norm_z[iad][iz]->SetTitle(Form("h_Edelayed_sub_norm_z_eh%dad%d_iz%d",EH[iad],AD[iad],iz+1));
			h_Edelayed_sub_DTnorm_z[iad][iz]->SetTitle(Form("h_Edelayed_sub_DTnorm_z_eh%dad%d_iz%d",EH[iad],AD[iad],iz+1));

			h_Edelayed_sub_rate_z[iad][iz]->SetName(Form("h_Edelayed_sub_rate_z_eh%dad%d_iz%d",EH[iad],AD[iad],iz+1));
			h_Edelayed_sub_norm_z[iad][iz]->SetName(Form("h_Edelayed_sub_norm_z_eh%dad%d_iz%d",EH[iad],AD[iad],iz+1));
			h_Edelayed_sub_DTnorm_z[iad][iz]->SetName(Form("h_Edelayed_sub_DTnorm_z_eh%dad%d_iz%d",EH[iad],AD[iad],iz+1));
		}
	}

	for(int iad = 0; iad < 8; iad++){
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			h_Edelayed_sub_rate_r2[iad][ir2] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_rate_r2_points_eh%dad%d_ir2%d", EH[iad], AD[iad],ir2));
			h_Edelayed_sub_norm_r2[iad][ir2] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_norm_r2_points_eh%dad%d_ir2%d", EH[iad], AD[iad],ir2));
			h_Edelayed_sub_DTnorm_r2[iad][ir2] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_DTnorm_r2_points_eh%dad%d_ir2%d", EH[iad], AD[iad],ir2));

			h_Edelayed_sub_rate_r2[iad][ir2]->SetTitle(Form("h_Edelayed_sub_rate_r2_eh%dad%d_ir2%d",EH[iad],AD[iad],ir2+1));
			h_Edelayed_sub_norm_r2[iad][ir2]->SetTitle(Form("h_Edelayed_sub_norm_r2_eh%dad%d_ir2%d",EH[iad],AD[iad],ir2+1));
			h_Edelayed_sub_DTnorm_r2[iad][ir2]->SetTitle(Form("h_Edelayed_sub_DTnorm_r2_eh%dad%d_ir2%d",EH[iad],AD[iad],ir2+1));

			h_Edelayed_sub_rate_r2[iad][ir2]->SetName(Form("h_Edelayed_sub_rate_r2_eh%dad%d_ir2%d",EH[iad],AD[iad],ir2+1));
			h_Edelayed_sub_norm_r2[iad][ir2]->SetName(Form("h_Edelayed_sub_norm_r2_eh%dad%d_ir2%d",EH[iad],AD[iad],ir2+1));
			h_Edelayed_sub_DTnorm_r2[iad][ir2]->SetName(Form("h_Edelayed_sub_DTnorm_r2_eh%dad%d_ir2%d",EH[iad],AD[iad],ir2+1));
		}
	}

	TH1F* h_Edelayed_sub_rate_zVSr2_near[Nr2Points][NzPoints];
	TH1F* h_Edelayed_sub_norm_zVSr2_near[Nr2Points][NzPoints];
	TH1F* h_Edelayed_sub_DTnorm_zVSr2_near[Nr2Points][NzPoints];

	TH1F* h_Edelayed_sub_rate_zVSr2_far[Nr2Points][NzPoints];
	TH1F* h_Edelayed_sub_norm_zVSr2_far[Nr2Points][NzPoints];
	TH1F* h_Edelayed_sub_DTnorm_zVSr2_far[Nr2Points][NzPoints];

	for(int iad = 0; iad < 8; iad++){
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			for(int iz = 0; iz < NzPoints; iz++){
				h_Edelayed_sub_rate_zVSr2[iad][ir2][iz] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_rate_zVSr2_points_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2, iz));
				h_Edelayed_sub_norm_zVSr2[iad][ir2][iz] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_norm_zVSr2_points_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2, iz));
				h_Edelayed_sub_DTnorm_zVSr2[iad][ir2][iz] = (TH1F*)outfile->Get(Form("h_Edelayed_sub_DTnorm_zVSr2_points_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2, iz));

				h_Edelayed_sub_rate_zVSr2[iad][ir2][iz]->SetTitle(Form("h_Ed_sub_rate_zVSr2_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2+1, iz+1));
				h_Edelayed_sub_norm_zVSr2[iad][ir2][iz]->SetTitle(Form("h_Edelayed_sub_norm_zVSr2_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2+1, iz+1));
				h_Edelayed_sub_DTnorm_zVSr2[iad][ir2][iz]->SetTitle(Form("h_Edelayed_sub_DTnorm_zVSr2_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2+1, iz+1));

				h_Edelayed_sub_rate_zVSr2[iad][ir2][iz]->SetName(Form("h_Edelayed_sub_rate_zVSr2_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2+1, iz+1));
				h_Edelayed_sub_norm_zVSr2[iad][ir2][iz]->SetName(Form("h_Edelayed_sub_norm_zVSr2_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2+1, iz+1));
				h_Edelayed_sub_DTnorm_zVSr2[iad][ir2][iz]->SetName(Form("h_Edelayed_sub_DTnorm_zVSr2_eh%dad%d_ir2%d_iz%d",EH[iad],AD[iad],ir2+1, iz+1));

	//			if(iad == 0){
	//				h_Edelayed_sub_rate_zVSr2_near[ir2][iz] = (TH1F*)h_Edelayed_sub_rate_zVSr2[iad][ir2][iz]->Clone();
	//				h_Edelayed_sub_norm_zVSr2_near[ir2][iz] = (TH1F*)h_Edelayed_sub_norm_zVSr2[iad][ir2][iz]->Clone();
	//				h_Edelayed_sub_DTnorm_zVSr2_near[ir2][iz] = (TH1F*)h_Edelayed_sub_DTnorm_zVSr2[iad][ir2][iz]->Clone();
	//			}
	//			else if(iad < 5){
	//				h_Edelayed_sub_rate_zVSr2_near[ir2][iz]->Add((TH1F*)h_Edelayed_sub_rate_zVSr2[iad][ir2][iz]);
	//				h_Edelayed_sub_norm_zVSr2_near[ir2][iz]->Add((TH1F*)h_Edelayed_sub_norm_zVSr2[iad][ir2][iz]);
	//				h_Edelayed_sub_DTnorm_zVSr2_near[ir2][iz]->Add((TH1F*)h_Edelayed_sub_DTnorm_zVSr2[iad][ir2][iz]);
	//			}

	//			if(iad == 5){
	//				h_Edelayed_sub_rate_zVSr2_far[ir2][iz] = (TH1F*)h_Edelayed_sub_rate_zVSr2[iad][ir2][iz]->Clone();
	//				h_Edelayed_sub_norm_zVSr2_far[ir2][iz] = (TH1F*)h_Edelayed_sub_norm_zVSr2[iad][ir2][iz]->Clone();
	//				h_Edelayed_sub_DTnorm_zVSr2_far[ir2][iz] = (TH1F*)h_Edelayed_sub_DTnorm_zVSr2[iad][ir2][iz]->Clone();
	//			}
	//			else{
	//				h_Edelayed_sub_rate_zVSr2_far[ir2][iz]->Add((TH1F*)h_Edelayed_sub_rate_zVSr2[iad][ir2][iz]);
	//				h_Edelayed_sub_norm_zVSr2_far[ir2][iz]->Add((TH1F*)h_Edelayed_sub_norm_zVSr2[iad][ir2][iz]);
	//				h_Edelayed_sub_DTnorm_zVSr2_far[ir2][iz]->Add((TH1F*)h_Edelayed_sub_DTnorm_zVSr2[iad][ir2][iz]);
	//			}
			}
		}
	}

	for(int ir2 = 0; ir2 < Nr2Points; ir2++){
		for(int iz = 0; iz < NzPoints; iz++){
			N_nom_rate_zVSr2_near[ir2][iz] = 0;
			N_nom_rate_zVSr2_far[ir2][iz] = 0; 
			N_nom_norm_zVSr2_near[ir2][iz] = 0;
			N_nom_norm_zVSr2_far[ir2][iz] = 0; 
			N_nom_DTnorm_zVSr2_near[ir2][iz] = 0;
			N_nom_DTnorm_zVSr2_far[ir2][iz] = 0;

			N_ext_rate_zVSr2_near[ir2][iz] = 0;
			N_ext_rate_zVSr2_far[ir2][iz] = 0; 
			N_ext_norm_zVSr2_near[ir2][iz] = 0;
			N_ext_norm_zVSr2_far[ir2][iz] = 0; 
			N_ext_DTnorm_zVSr2_near[ir2][iz] = 0;
			N_ext_DTnorm_zVSr2_far[ir2][iz] = 0; 
			for(int iad = 0; iad < 4; iad++){
				N_nom_rate_zVSr2_near[ir2][iz] += N_nom_rate_zVSr2[iad][ir2][iz];
				N_nom_norm_zVSr2_near[ir2][iz] += N_nom_norm_zVSr2[iad][ir2][iz];
				N_nom_norm_zVSr2_near[ir2][iz] += N_nom_norm_zVSr2[iad][ir2][iz];

				N_nom_rate_zVSr2_far[ir2][iz] += N_nom_rate_zVSr2[iad+4][ir2][iz];
				N_nom_norm_zVSr2_far[ir2][iz] += N_nom_norm_zVSr2[iad+4][ir2][iz];
				N_nom_norm_zVSr2_far[ir2][iz] += N_nom_norm_zVSr2[iad+4][ir2][iz];

				N_ext_rate_zVSr2_near[ir2][iz] += N_ext_rate_zVSr2[iad][ir2][iz];
				N_ext_norm_zVSr2_near[ir2][iz] += N_ext_norm_zVSr2[iad][ir2][iz];
				N_ext_norm_zVSr2_near[ir2][iz] += N_ext_norm_zVSr2[iad][ir2][iz];

				N_ext_rate_zVSr2_far[ir2][iz] += N_ext_rate_zVSr2[iad+4][ir2][iz];
				N_ext_norm_zVSr2_far[ir2][iz] += N_ext_norm_zVSr2[iad+4][ir2][iz];
				N_ext_norm_zVSr2_far[ir2][iz] += N_ext_norm_zVSr2[iad+4][ir2][iz];
			}
			efficiency_rate_zVSr2_near[ir2][iz] = N_nom_rate_zVSr2_near[ir2][iz]/N_ext_rate_zVSr2_near[ir2][iz];
			efficiency_norm_zVSr2_near[ir2][iz] = N_nom_norm_zVSr2_near[ir2][iz]/N_ext_norm_zVSr2_near[ir2][iz];
			efficiency_DTnorm_zVSr2_near[ir2][iz] = N_nom_DTnorm_zVSr2_near[ir2][iz]/N_ext_DTnorm_zVSr2_near[ir2][iz];

			efficiency_rate_zVSr2_far[ir2][iz] = N_nom_rate_zVSr2_far[ir2][iz]/N_ext_rate_zVSr2_far[ir2][iz];
			efficiency_norm_zVSr2_far[ir2][iz] = N_nom_norm_zVSr2_far[ir2][iz]/N_ext_norm_zVSr2_far[ir2][iz];
			efficiency_DTnorm_zVSr2_far[ir2][iz] = N_nom_DTnorm_zVSr2_far[ir2][iz]/N_ext_DTnorm_zVSr2_far[ir2][iz];
		}
	}
*/

	//histograms - values vs AD:
	double peakAvg = 0;
	TH1D* h_peak_rate_ADs=new TH1D("h_peak_rate_ADs","h_peak_rate_ADs",8,0.4,8.4);
	TH1D* h_peak_norm_ADs=new TH1D("h_peak_norm_ADs","h_peak_norm_ADs",8,0.5,8.5);
	TH1D* h_peak_DTnorm_ADs=new TH1D("h_peak_DTnorm_ADs","h_peak_DTnorm_ADs",8,0.6,8.6);
	double peakAvg_Sam = 0;
	TH1D* h_peak_rate_ADs_Sam=new TH1D("h_peak_rate_ADs_Sam","h_peak_rate_ADs_1.5",8,0.45,8.45);
	TH1D* h_peak_norm_ADs_Sam=new TH1D("h_peak_norm_ADs_Sam","h_peak_norm_ADs_1.5",8,0.55,8.55);
	TH1D* h_peak_DTnorm_ADs_Sam=new TH1D("h_peak_DTnorm_ADs_Sam","h_peak_DTnorm_ADs_1.5",8,0.65,8.65);

	double sigmaAvg = 0;
	TH1D* h_sigma_rate_ADs=new TH1D("h_sigma_rate_ADs","h_sigma_rate_ADs",8,0.4,8.4);
	TH1D* h_sigma_norm_ADs=new TH1D("h_sigma_norm_ADs","h_sigma_norm_ADs",8,0.5,8.5);
	TH1D* h_sigma_DTnorm_ADs=new TH1D("h_sigma_DTnorm_ADs","h_sigma_DTnorm_ADs",8,0.6,8.6);
	double sigmaAvg_Sam = 0;
	TH1D* h_sigma_rate_ADs_Sam=new TH1D("h_sigma_rate_ADs_Sam","h_sigma_rate_ADs_1.5",8,0.45,8.45);
	TH1D* h_sigma_norm_ADs_Sam=new TH1D("h_sigma_norm_ADs_Sam","h_sigma_norm_ADs_1.5",8,0.55,8.55);
	TH1D* h_sigma_DTnorm_ADs_Sam=new TH1D("h_sigma_DTnorm_ADs_Sam","h_sigma_DTnorm_ADs_1.5",8,0.65,8.65);

	double alphaAvg = 0;
	TH1D* h_alpha_rate_ADs=new TH1D("h_alpha_rate_ADs","h_alpha_rate_ADs",8,0.4,8.4);
	TH1D* h_alpha_norm_ADs=new TH1D("h_alpha_norm_ADs","h_alpha_norm_ADs",8,0.5,8.5);
	TH1D* h_alpha_DTnorm_ADs=new TH1D("h_alpha_DTnorm_ADs","h_alpha_DTnorm_ADs",8,0.6,8.6);
	double alphaAvg_Sam = 0;
	TH1D* h_alpha_rate_ADs_Sam=new TH1D("h_alpha_rate_ADs_Sam","h_alpha_rate_ADs_1.5",8,0.45,8.45);
	TH1D* h_alpha_norm_ADs_Sam=new TH1D("h_alpha_norm_ADs_Sam","h_alpha_norm_ADs_1.5",8,0.55,8.55);
	TH1D* h_alpha_DTnorm_ADs_Sam=new TH1D("h_alpha_DTnorm_ADs_Sam","h_alpha_DTnorm_ADs_1.5",8,0.65,8.65);

	double lambdaAvg = 0;
	TH1D* h_lambda_rate_ADs=new TH1D("h_lambda_rate_ADs","h_lambda_rate_ADs",8,0.4,8.4);
	TH1D* h_lambda_norm_ADs=new TH1D("h_lambda_norm_ADs","h_lambda_norm_ADs",8,0.5,8.5);
	TH1D* h_lambda_DTnorm_ADs=new TH1D("h_lambda_DTnorm_ADs","h_lambda_DTnorm_ADs",8,0.6,8.6);
	double lambdaAvg_Sam = 0;
	TH1D* h_lambda_rate_ADs_Sam=new TH1D("h_lambda_rate_ADs_Sam","h_lambda_rate_ADs_1.5",8,0.45,8.45);
	TH1D* h_lambda_norm_ADs_Sam=new TH1D("h_lambda_norm_ADs_Sam","h_lambda_norm_ADs_1.5",8,0.55,8.55);
	TH1D* h_lambda_DTnorm_ADs_Sam=new TH1D("h_lambda_DTnorm_ADs_Sam","h_lambda_DTnorm_ADs_1.5",8,0.65,8.65);

	double efficiencyAvg = 0;
	TH1D* h_efficiency_rate_ADs=new TH1D("h_efficiency_rate_ADs","h_efficiency_rate_ADs",8,0.4,8.4);
	TH1D* h_efficiency_norm_ADs=new TH1D("h_efficiency_norm_ADs","h_efficiency_norm_ADs",8,0.5,8.5);
	TH1D* h_efficiency_DTnorm_ADs=new TH1D("h_efficiency_DTnorm_ADs","h_efficiency_DTnorm_ADs",8,0.6,8.6);
	double efficiencyAvg_Sam = 0;
	TH1D* h_efficiency_rate_ADs_Sam=new TH1D("h_efficiency_rate_ADs_Sam","h_efficiency_rate_ADs_1.5",8,0.45,8.45);
	TH1D* h_efficiency_norm_ADs_Sam=new TH1D("h_efficiency_norm_ADs_Sam","h_efficiency_norm_ADs_1.5",8,0.55,8.55);
	TH1D* h_efficiency_DTnorm_ADs_Sam=new TH1D("h_efficiency_DTnorm_ADs_Sam","h_efficiency_DTnorm_ADs_1.5",8,0.65,8.65);


	for(int iad = 0; iad < 8; iad++){
		h_peak_rate_ADs->Fill(iad+1, peak_rate[iad]);
		h_peak_norm_ADs->Fill(iad+1, peak_norm[iad]);
		h_peak_DTnorm_ADs->Fill(iad+1, peak_DTnorm[iad]);
		peakAvg = peakAvg + (peak_rate[iad] + peak_norm[iad] + peak_DTnorm[iad])/3.;
		h_peak_rate_ADs_Sam->Fill(iad+1, peak_rate_Sam[iad]);
		h_peak_norm_ADs_Sam->Fill(iad+1, peak_norm_Sam[iad]);
		h_peak_DTnorm_ADs_Sam->Fill(iad+1, peak_DTnorm_Sam[iad]);
		peakAvg_Sam = peakAvg_Sam + (peak_rate_Sam[iad] + peak_norm_Sam[iad] + peak_DTnorm_Sam[iad])/3.;
			h_peak_rate_ADs->SetBinError(h_peak_rate_ADs->FindBin(iad+1), peakError_rate[iad]);
			h_peak_norm_ADs->SetBinError(h_peak_norm_ADs->FindBin(iad+1), peakError_norm[iad]);
			h_peak_DTnorm_ADs->SetBinError(h_peak_DTnorm_ADs->FindBin(iad+1), peakError_DTnorm[iad]);
			h_peak_rate_ADs_Sam->SetBinError(h_peak_rate_ADs_Sam->FindBin(iad+1), peakError_rate_Sam[iad]);
			h_peak_norm_ADs_Sam->SetBinError(h_peak_norm_ADs_Sam->FindBin(iad+1), peakError_norm_Sam[iad]);
			h_peak_DTnorm_ADs_Sam->SetBinError(h_peak_DTnorm_ADs_Sam->FindBin(iad+1), peakError_DTnorm_Sam[iad]);

		h_sigma_rate_ADs->Fill(iad+1, sigma_rate[iad]);
		h_sigma_norm_ADs->Fill(iad+1, sigma_norm[iad]);
		h_sigma_DTnorm_ADs->Fill(iad+1, sigma_DTnorm[iad]);
		sigmaAvg = sigmaAvg + (sigma_rate[iad] + sigma_norm[iad] + sigma_DTnorm[iad])/3.;
		h_sigma_rate_ADs_Sam->Fill(iad+1, sigma_rate_Sam[iad]);
		h_sigma_norm_ADs_Sam->Fill(iad+1, sigma_norm_Sam[iad]);
		h_sigma_DTnorm_ADs_Sam->Fill(iad+1, sigma_DTnorm_Sam[iad]);
		sigmaAvg_Sam = sigmaAvg_Sam + (sigma_rate_Sam[iad] + sigma_norm_Sam[iad] + sigma_DTnorm_Sam[iad])/3.;
			h_sigma_rate_ADs->SetBinError(h_sigma_rate_ADs->FindBin(iad+1), sigmaError_rate[iad]);
			h_sigma_norm_ADs->SetBinError(h_sigma_norm_ADs->FindBin(iad+1), sigmaError_norm[iad]);
			h_sigma_DTnorm_ADs->SetBinError(h_sigma_DTnorm_ADs->FindBin(iad+1), sigmaError_DTnorm[iad]);
			h_sigma_rate_ADs_Sam->SetBinError(h_sigma_rate_ADs_Sam->FindBin(iad+1), sigmaError_rate_Sam[iad]);
			h_sigma_norm_ADs_Sam->SetBinError(h_sigma_norm_ADs_Sam->FindBin(iad+1), sigmaError_norm_Sam[iad]);
			h_sigma_DTnorm_ADs_Sam->SetBinError(h_sigma_DTnorm_ADs_Sam->FindBin(iad+1), sigmaError_DTnorm_Sam[iad]);

		h_alpha_rate_ADs->Fill(iad+1, alpha_rate[iad]);
		h_alpha_norm_ADs->Fill(iad+1, alpha_norm[iad]);
		h_alpha_DTnorm_ADs->Fill(iad+1, alpha_DTnorm[iad]);
		alphaAvg = alphaAvg + (alpha_rate[iad] + alpha_norm[iad] + alpha_DTnorm[iad])/3.;
		h_alpha_rate_ADs_Sam->Fill(iad+1, alpha_rate_Sam[iad]);
		h_alpha_norm_ADs_Sam->Fill(iad+1, alpha_norm_Sam[iad]);
		h_alpha_DTnorm_ADs_Sam->Fill(iad+1, alpha_DTnorm_Sam[iad]);
		alphaAvg_Sam = alphaAvg_Sam + (alpha_rate_Sam[iad] + alpha_norm_Sam[iad] + alpha_DTnorm_Sam[iad])/3.;
			h_alpha_rate_ADs->SetBinError(h_alpha_rate_ADs->FindBin(iad+1), alphaError_rate[iad]);
			h_alpha_norm_ADs->SetBinError(h_alpha_norm_ADs->FindBin(iad+1), alphaError_norm[iad]);
			h_alpha_DTnorm_ADs->SetBinError(h_alpha_DTnorm_ADs->FindBin(iad+1), alphaError_DTnorm[iad]);
			h_alpha_rate_ADs_Sam->SetBinError(h_alpha_rate_ADs_Sam->FindBin(iad+1), alphaError_rate_Sam[iad]);
			h_alpha_norm_ADs_Sam->SetBinError(h_alpha_norm_ADs_Sam->FindBin(iad+1), alphaError_norm_Sam[iad]);
			h_alpha_DTnorm_ADs_Sam->SetBinError(h_alpha_DTnorm_ADs_Sam->FindBin(iad+1), alphaError_DTnorm_Sam[iad]);

		h_lambda_rate_ADs->Fill(iad+1, lambda_rate[iad]);
		h_lambda_norm_ADs->Fill(iad+1, lambda_norm[iad]);
		h_lambda_DTnorm_ADs->Fill(iad+1, lambda_DTnorm[iad]);
		lambdaAvg = lambdaAvg + (lambda_rate[iad] + lambda_norm[iad] + lambda_DTnorm[iad])/3.;
		h_lambda_rate_ADs_Sam->Fill(iad+1, lambda_rate_Sam[iad]);
		h_lambda_norm_ADs_Sam->Fill(iad+1, lambda_norm_Sam[iad]);
		h_lambda_DTnorm_ADs_Sam->Fill(iad+1, lambda_DTnorm_Sam[iad]);
		lambdaAvg_Sam = lambdaAvg_Sam + (lambda_rate_Sam[iad] + lambda_norm_Sam[iad] + lambda_DTnorm_Sam[iad])/3.;
			h_lambda_rate_ADs->SetBinError(h_lambda_rate_ADs->FindBin(iad+1), lambdaError_rate[iad]);
			h_lambda_norm_ADs->SetBinError(h_lambda_norm_ADs->FindBin(iad+1), lambdaError_norm[iad]);
			h_lambda_DTnorm_ADs->SetBinError(h_lambda_DTnorm_ADs->FindBin(iad+1), lambdaError_DTnorm[iad]);
			h_lambda_rate_ADs_Sam->SetBinError(h_lambda_rate_ADs_Sam->FindBin(iad+1), lambdaError_rate_Sam[iad]);
			h_lambda_norm_ADs_Sam->SetBinError(h_lambda_norm_ADs_Sam->FindBin(iad+1), lambdaError_norm_Sam[iad]);
			h_lambda_DTnorm_ADs_Sam->SetBinError(h_lambda_DTnorm_ADs_Sam->FindBin(iad+1), lambdaError_DTnorm_Sam[iad]);

		h_efficiency_rate_ADs->Fill(iad+1, efficiency_rate[iad]);
		h_efficiency_norm_ADs->Fill(iad+1, efficiency_norm[iad]);
		h_efficiency_DTnorm_ADs->Fill(iad+1, efficiency_DTnorm[iad]);
		efficiencyAvg = efficiencyAvg + (efficiency_rate[iad] + efficiency_norm[iad] + efficiency_DTnorm[iad])/3.;
		h_efficiency_rate_ADs_Sam->Fill(iad+1, efficiency_rate_Sam[iad]);
		h_efficiency_norm_ADs_Sam->Fill(iad+1, efficiency_norm_Sam[iad]);
		h_efficiency_DTnorm_ADs_Sam->Fill(iad+1, efficiency_DTnorm_Sam[iad]);
		efficiencyAvg_Sam = efficiencyAvg_Sam + (efficiency_rate_Sam[iad] + efficiency_norm_Sam[iad] + efficiency_DTnorm_Sam[iad])/3.;
			h_efficiency_rate_ADs->SetBinError(h_efficiency_rate_ADs->FindBin(iad+1), efficiencyError_rate[iad]);
			h_efficiency_norm_ADs->SetBinError(h_efficiency_norm_ADs->FindBin(iad+1), efficiencyError_norm[iad]);
			h_efficiency_DTnorm_ADs->SetBinError(h_efficiency_DTnorm_ADs->FindBin(iad+1), efficiencyError_DTnorm[iad]);
			h_efficiency_rate_ADs_Sam->SetBinError(h_efficiency_rate_ADs_Sam->FindBin(iad+1), efficiencyError_rate_Sam[iad]);
			h_efficiency_norm_ADs_Sam->SetBinError(h_efficiency_norm_ADs_Sam->FindBin(iad+1), efficiencyError_norm_Sam[iad]);
			h_efficiency_DTnorm_ADs_Sam->SetBinError(h_efficiency_DTnorm_ADs_Sam->FindBin(iad+1), efficiencyError_DTnorm_Sam[iad]);

		cout << "AD" << iad+1 << "\t" << peak_rate[iad] << "\t" << sigma_rate[iad] << endl;
	}

		peakAvg = peakAvg/8.;
		sigmaAvg = sigmaAvg/8.;
		alphaAvg = alphaAvg/8.;
		lambdaAvg = lambdaAvg/8.;
		efficiencyAvg = efficiencyAvg/8.;

		peakAvg_Sam = peakAvg_Sam/8.;
		sigmaAvg_Sam = sigmaAvg_Sam/8.;
		alphaAvg_Sam = alphaAvg_Sam/8.;
		lambdaAvg_Sam = lambdaAvg_Sam/8.;
		efficiencyAvg_Sam = efficiencyAvg_Sam/8.;

/*
//Plots of Beda's Fit
	TCanvas *peakVsAD = new TCanvas("peakVsAD","peakVsAD");
		peakVsAD->cd();
		h_peak_rate_ADs->SetStats(0);
		h_peak_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_peak_rate_ADs->GetYaxis()->SetTitle("Peak Energy [MeV]");
		h_peak_rate_ADs->SetLineWidth(0);
		h_peak_norm_ADs->SetLineWidth(0);
		h_peak_DTnorm_ADs->SetLineWidth(0);
		h_peak_rate_ADs->SetMarkerStyle(21);
		h_peak_norm_ADs->SetMarkerStyle(20);
		h_peak_DTnorm_ADs->SetMarkerStyle(22);
		h_peak_rate_ADs->SetMarkerColor(kRed);
		h_peak_norm_ADs->SetMarkerColor(kBlue);
		h_peak_DTnorm_ADs->SetMarkerColor(kBlack);
		h_peak_rate_ADs->SetMarkerSize(1.5);
		h_peak_norm_ADs->SetMarkerSize(1.5);
		h_peak_DTnorm_ADs->SetMarkerSize(1.5);

		h_peak_rate_ADs->GetYaxis()->SetRangeUser(peakAvg - peakAvg*0.005, peakAvg + peakAvg*0.005);

		h_peak_rate_ADs->Draw();
		h_peak_norm_ADs->Draw("same");
		h_peak_DTnorm_ADs->Draw("same");
	peakVsAD->BuildLegend();

	TCanvas *sigmaVsAD = new TCanvas("sigmaVsAD","sigmaVsAD");
		sigmaVsAD->cd();
		h_sigma_rate_ADs->SetStats(0);
		h_sigma_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_sigma_rate_ADs->GetYaxis()->SetTitle("Sigma Value");
		h_sigma_rate_ADs->SetLineWidth(0);
		h_sigma_norm_ADs->SetLineWidth(0);
		h_sigma_DTnorm_ADs->SetLineWidth(0);
		h_sigma_rate_ADs->SetMarkerStyle(21);
		h_sigma_norm_ADs->SetMarkerStyle(20);
		h_sigma_DTnorm_ADs->SetMarkerStyle(22);
		h_sigma_rate_ADs->SetMarkerColor(kRed);
		h_sigma_norm_ADs->SetMarkerColor(kBlue);
		h_sigma_DTnorm_ADs->SetMarkerColor(kBlack);
		h_sigma_rate_ADs->SetMarkerSize(1.5);
		h_sigma_norm_ADs->SetMarkerSize(1.5);
		h_sigma_DTnorm_ADs->SetMarkerSize(1.5);

		h_sigma_rate_ADs->GetYaxis()->SetRangeUser(sigmaAvg - 0.005, sigmaAvg + 0.005);

		h_sigma_rate_ADs->Draw();
		h_sigma_norm_ADs->Draw("same");
		h_sigma_DTnorm_ADs->Draw("same");
	sigmaVsAD->BuildLegend();

	TCanvas *efficiencyVsAD = new TCanvas("efficiencyVsAD","efficiencyVsAD");
		efficiencyVsAD->cd();
		h_efficiency_rate_ADs->SetStats(0);
		h_efficiency_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_efficiency_rate_ADs->GetYaxis()->SetTitle("Delayed Cut Efficiency");
		h_efficiency_rate_ADs->SetLineWidth(0);
		h_efficiency_norm_ADs->SetLineWidth(0);
		h_efficiency_DTnorm_ADs->SetLineWidth(0);
		h_efficiency_rate_ADs->SetMarkerStyle(21);
		h_efficiency_norm_ADs->SetMarkerStyle(20);
		h_efficiency_DTnorm_ADs->SetMarkerStyle(22);
		h_efficiency_rate_ADs->SetMarkerColor(kRed);
		h_efficiency_norm_ADs->SetMarkerColor(kBlue);
		h_efficiency_DTnorm_ADs->SetMarkerColor(kBlack);
		h_efficiency_rate_ADs->SetMarkerSize(1.5);
		h_efficiency_norm_ADs->SetMarkerSize(1.5);
		h_efficiency_DTnorm_ADs->SetMarkerSize(1.5);

		h_efficiency_rate_ADs->GetYaxis()->SetRangeUser(efficiencyAvg - 0.015, efficiencyAvg + 0.015);

		h_efficiency_rate_ADs->Draw();
		h_efficiency_norm_ADs->Draw("same");
		h_efficiency_DTnorm_ADs->Draw("same");
	efficiencyVsAD->BuildLegend();

//Plots of Sam's Fit
	TCanvas *peakVsAD_Sam = new TCanvas("peakVsAD_Sam","peakVsAD_Sam");
		peakVsAD_Sam->cd();
		h_peak_rate_ADs_Sam->SetStats(0);
		h_peak_rate_ADs_Sam->GetXaxis()->SetTitle("AD Number");
		h_peak_rate_ADs_Sam->GetYaxis()->SetTitle("Peak Energy [MeV]");
		h_peak_rate_ADs_Sam->SetLineWidth(0);
		h_peak_norm_ADs_Sam->SetLineWidth(0);
		h_peak_DTnorm_ADs_Sam->SetLineWidth(0);
		h_peak_rate_ADs_Sam->SetMarkerStyle(21);
		h_peak_norm_ADs_Sam->SetMarkerStyle(20);
		h_peak_DTnorm_ADs_Sam->SetMarkerStyle(22);
		h_peak_rate_ADs_Sam->SetMarkerColor(kRed);
		h_peak_norm_ADs_Sam->SetMarkerColor(kBlue);
		h_peak_DTnorm_ADs_Sam->SetMarkerColor(kBlack);
		h_peak_rate_ADs_Sam->SetMarkerSize(1.5);
		h_peak_norm_ADs_Sam->SetMarkerSize(1.5);
		h_peak_DTnorm_ADs_Sam->SetMarkerSize(1.5);

		h_peak_rate_ADs_Sam->GetYaxis()->SetRangeUser(peakAvg_Sam - peakAvg_Sam*0.005, peakAvg_Sam + peakAvg_Sam*0.005);

		h_peak_rate_ADs_Sam->Draw();
		h_peak_norm_ADs_Sam->Draw("same");
		h_peak_DTnorm_ADs_Sam->Draw("same");
	peakVsAD_Sam->BuildLegend();

	TCanvas *sigmaVsAD_Sam = new TCanvas("sigmaVsAD_Sam","sigmaVsAD_Sam");
		sigmaVsAD_Sam->cd();
		h_sigma_rate_ADs_Sam->SetStats(0);
		h_sigma_rate_ADs_Sam->GetXaxis()->SetTitle("AD Number");
		h_sigma_rate_ADs_Sam->GetYaxis()->SetTitle("Sigma Value");
		h_sigma_rate_ADs_Sam->SetLineWidth(0);
		h_sigma_norm_ADs_Sam->SetLineWidth(0);
		h_sigma_DTnorm_ADs_Sam->SetLineWidth(0);
		h_sigma_rate_ADs_Sam->SetMarkerStyle(21);
		h_sigma_norm_ADs_Sam->SetMarkerStyle(20);
		h_sigma_DTnorm_ADs_Sam->SetMarkerStyle(22);
		h_sigma_rate_ADs_Sam->SetMarkerColor(kRed);
		h_sigma_norm_ADs_Sam->SetMarkerColor(kBlue);
		h_sigma_DTnorm_ADs_Sam->SetMarkerColor(kBlack);
		h_sigma_rate_ADs_Sam->SetMarkerSize(1.5);
		h_sigma_norm_ADs_Sam->SetMarkerSize(1.5);
		h_sigma_DTnorm_ADs_Sam->SetMarkerSize(1.5);

		h_sigma_rate_ADs_Sam->GetYaxis()->SetRangeUser(sigmaAvg_Sam - 0.005, sigmaAvg_Sam + 0.005);

		h_sigma_rate_ADs_Sam->Draw();
		h_sigma_norm_ADs_Sam->Draw("same");
		h_sigma_DTnorm_ADs_Sam->Draw("same");
	sigmaVsAD_Sam->BuildLegend();

	TCanvas *efficiencyVsAD_Sam = new TCanvas("efficiencyVsAD_Sam","efficiencyVsAD_Sam");
		efficiencyVsAD_Sam->cd();
		h_efficiency_rate_ADs_Sam->SetStats(0);
		h_efficiency_rate_ADs_Sam->GetXaxis()->SetTitle("AD Number");
		h_efficiency_rate_ADs_Sam->GetYaxis()->SetTitle("Delayed Cut Efficiency");
		h_efficiency_rate_ADs_Sam->SetLineWidth(0);
		h_efficiency_norm_ADs_Sam->SetLineWidth(0);
		h_efficiency_DTnorm_ADs_Sam->SetLineWidth(0);
		h_efficiency_rate_ADs_Sam->SetMarkerStyle(21);
		h_efficiency_norm_ADs_Sam->SetMarkerStyle(20);
		h_efficiency_DTnorm_ADs_Sam->SetMarkerStyle(22);
		h_efficiency_rate_ADs_Sam->SetMarkerColor(kRed);
		h_efficiency_norm_ADs_Sam->SetMarkerColor(kBlue);
		h_efficiency_DTnorm_ADs_Sam->SetMarkerColor(kBlack);
		h_efficiency_rate_ADs_Sam->SetMarkerSize(1.5);
		h_efficiency_norm_ADs_Sam->SetMarkerSize(1.5);
		h_efficiency_DTnorm_ADs_Sam->SetMarkerSize(1.5);

		h_efficiency_rate_ADs_Sam->GetYaxis()->SetRangeUser(efficiencyAvg_Sam - 0.015, efficiencyAvg_Sam + 0.015);

		h_efficiency_rate_ADs_Sam->Draw();
		h_efficiency_norm_ADs_Sam->Draw("same");
		h_efficiency_DTnorm_ADs_Sam->Draw("same");
	efficiencyVsAD_Sam->BuildLegend();
*/

//*****************************Plots of all Fits:*******************************************************
	TCanvas *peakVsAD_AllFits = new TCanvas("peakVsAD_AllFits","peakVsAD_AllFits");
		peakVsAD_AllFits->cd();
		h_peak_rate_ADs->SetStats(0);
		h_peak_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_peak_rate_ADs->GetYaxis()->SetTitle("Peak Energy [MeV]");
	//	h_peak_rate_ADs->SetLineWidth(0);
	//	h_peak_norm_ADs->SetLineWidth(0);
	//	h_peak_DTnorm_ADs->SetLineWidth(0);
	//	h_peak_rate_ADs_Sam->SetLineWidth(0);
	//	h_peak_norm_ADs_Sam->SetLineWidth(0);
	//	h_peak_DTnorm_ADs_Sam->SetLineWidth(0);
		h_peak_rate_ADs->SetMarkerStyle(21);
		h_peak_norm_ADs->SetMarkerStyle(20);
		h_peak_DTnorm_ADs->SetMarkerStyle(22);
		h_peak_rate_ADs->SetMarkerColor(kRed);
		h_peak_norm_ADs->SetMarkerColor(kBlue);
		h_peak_DTnorm_ADs->SetMarkerColor(kBlack);
		h_peak_rate_ADs->SetLineColor(kRed);
		h_peak_norm_ADs->SetLineColor(kBlue);
		h_peak_DTnorm_ADs->SetLineColor(kBlack);
		h_peak_rate_ADs->SetMarkerSize(1.5);
		h_peak_norm_ADs->SetMarkerSize(1.5);
		h_peak_DTnorm_ADs->SetMarkerSize(1.5);
		h_peak_rate_ADs_Sam->SetMarkerSize(1.5);
		h_peak_norm_ADs_Sam->SetMarkerSize(1.5);
		h_peak_DTnorm_ADs_Sam->SetMarkerSize(1.5);
		h_peak_rate_ADs_Sam->SetMarkerStyle(25);
		h_peak_norm_ADs_Sam->SetMarkerStyle(24);
		h_peak_DTnorm_ADs_Sam->SetMarkerStyle(26);
		h_peak_rate_ADs_Sam->SetMarkerColor(kMagenta);
		h_peak_norm_ADs_Sam->SetMarkerColor(kCyan);
		h_peak_DTnorm_ADs_Sam->SetMarkerColor(kGray+1);
		h_peak_rate_ADs_Sam->SetLineColor(kMagenta);
		h_peak_norm_ADs_Sam->SetLineColor(kCyan);
		h_peak_DTnorm_ADs_Sam->SetLineColor(kGray+1);

		h_peak_rate_ADs->GetYaxis()->SetRangeUser(peakAvg - peakAvg*0.005, peakAvg + peakAvg*0.005);

		h_peak_rate_ADs->Draw("e1x0");
		h_peak_norm_ADs->Draw("e1x0 same");
		h_peak_DTnorm_ADs->Draw("e1x0 same");
	//	h_peak_rate_ADs_Sam->Draw("e1x0 same");
	//	h_peak_norm_ADs_Sam->Draw("e1x0 same");
	//	h_peak_DTnorm_ADs_Sam->Draw("e1x0 same");
	peakVsAD_AllFits->BuildLegend();

	TCanvas *sigmaVsAD_AllFits = new TCanvas("sigmaVsAD_AllFits","sigmaVsAD_AllFits");
		sigmaVsAD_AllFits->cd();
		h_sigma_rate_ADs->SetStats(0);
		h_sigma_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_sigma_rate_ADs->GetYaxis()->SetTitle("Sigma Value");
	//	h_sigma_rate_ADs->SetLineWidth(0);
	//	h_sigma_norm_ADs->SetLineWidth(0);
	//	h_sigma_DTnorm_ADs->SetLineWidth(0);
	//	h_sigma_rate_ADs_Sam->SetLineWidth(0);
	//	h_sigma_norm_ADs_Sam->SetLineWidth(0);
	//	h_sigma_DTnorm_ADs_Sam->SetLineWidth(0);
		h_sigma_rate_ADs->SetMarkerStyle(21);
		h_sigma_norm_ADs->SetMarkerStyle(20);
		h_sigma_DTnorm_ADs->SetMarkerStyle(22);
		h_sigma_rate_ADs->SetMarkerColor(kRed);
		h_sigma_norm_ADs->SetMarkerColor(kBlue);
		h_sigma_DTnorm_ADs->SetMarkerColor(kBlack);
		h_sigma_rate_ADs->SetMarkerSize(1.5);
		h_sigma_norm_ADs->SetMarkerSize(1.5);
		h_sigma_DTnorm_ADs->SetMarkerSize(1.5);
		h_sigma_rate_ADs_Sam->SetMarkerSize(1.5);
		h_sigma_norm_ADs_Sam->SetMarkerSize(1.5);
		h_sigma_DTnorm_ADs_Sam->SetMarkerSize(1.5);
		h_sigma_rate_ADs_Sam->SetMarkerStyle(25);
		h_sigma_norm_ADs_Sam->SetMarkerStyle(24);
		h_sigma_DTnorm_ADs_Sam->SetMarkerStyle(26);
		h_sigma_rate_ADs_Sam->SetMarkerColor(kMagenta);
		h_sigma_norm_ADs_Sam->SetMarkerColor(kCyan);
		h_sigma_DTnorm_ADs_Sam->SetMarkerColor(kGray+1);
		h_sigma_rate_ADs->SetLineColor(kRed);
		h_sigma_norm_ADs->SetLineColor(kBlue);
		h_sigma_DTnorm_ADs->SetLineColor(kBlack);
		h_sigma_rate_ADs_Sam->SetLineColor(kMagenta);
		h_sigma_norm_ADs_Sam->SetLineColor(kCyan);
		h_sigma_DTnorm_ADs_Sam->SetLineColor(kGray+1);

		h_sigma_rate_ADs->GetYaxis()->SetRangeUser(sigmaAvg - 0.005, sigmaAvg + 0.005);

		h_sigma_rate_ADs->Draw("e1x0");
		h_sigma_norm_ADs->Draw("e1x0 same");
		h_sigma_DTnorm_ADs->Draw("e1x0 same");
	//	h_sigma_rate_ADs_Sam->Draw("e1x0 same");
	//	h_sigma_norm_ADs_Sam->Draw("e1x0 same");
	//	h_sigma_DTnorm_ADs_Sam->Draw("e1x0 same");
	sigmaVsAD_AllFits->BuildLegend();

	TCanvas *alphaVsAD_AllFits = new TCanvas("alphaVsAD_AllFits","alphaVsAD_AllFits");
		alphaVsAD_AllFits->cd();
		h_alpha_rate_ADs->SetStats(0);
		h_alpha_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_alpha_rate_ADs->GetYaxis()->SetTitle("Peak Fraction");
	//	h_alpha_rate_ADs->SetLineWidth(0);
	//	h_alpha_norm_ADs->SetLineWidth(0);
	//	h_alpha_DTnorm_ADs->SetLineWidth(0);
	//	h_alpha_rate_ADs_Sam->SetLineWidth(0);
	//	h_alpha_norm_ADs_Sam->SetLineWidth(0);
	//	h_alpha_DTnorm_ADs_Sam->SetLineWidth(0);
		h_alpha_rate_ADs->SetMarkerStyle(21);
		h_alpha_norm_ADs->SetMarkerStyle(20);
		h_alpha_DTnorm_ADs->SetMarkerStyle(22);
		h_alpha_rate_ADs->SetMarkerColor(kRed);
		h_alpha_norm_ADs->SetMarkerColor(kBlue);
		h_alpha_DTnorm_ADs->SetMarkerColor(kBlack);
		h_alpha_rate_ADs->SetMarkerSize(1.5);
		h_alpha_norm_ADs->SetMarkerSize(1.5);
		h_alpha_DTnorm_ADs->SetMarkerSize(1.5);
		h_alpha_rate_ADs_Sam->SetMarkerSize(1.5);
		h_alpha_norm_ADs_Sam->SetMarkerSize(1.5);
		h_alpha_DTnorm_ADs_Sam->SetMarkerSize(1.5);
		h_alpha_rate_ADs_Sam->SetMarkerStyle(25);
		h_alpha_norm_ADs_Sam->SetMarkerStyle(24);
		h_alpha_DTnorm_ADs_Sam->SetMarkerStyle(26);
		h_alpha_rate_ADs_Sam->SetMarkerColor(kMagenta);
		h_alpha_norm_ADs_Sam->SetMarkerColor(kCyan);
		h_alpha_DTnorm_ADs_Sam->SetMarkerColor(kGray+1);
		h_alpha_rate_ADs->SetLineColor(kRed);
		h_alpha_norm_ADs->SetLineColor(kBlue);
		h_alpha_DTnorm_ADs->SetLineColor(kBlack);
		h_alpha_rate_ADs_Sam->SetLineColor(kMagenta);
		h_alpha_norm_ADs_Sam->SetLineColor(kCyan);
		h_alpha_DTnorm_ADs_Sam->SetLineColor(kGray+1);

	//	h_alpha_rate_ADs->GetYaxis()->SetRangeUser(alphaAvg - alphaAvg*0.005, alphaAvg + alphaAvg*0.005);

		h_alpha_rate_ADs->Draw("e1x0");
		h_alpha_norm_ADs->Draw("e1x0 same");
		h_alpha_DTnorm_ADs->Draw("e1x0 same");
	//	h_alpha_rate_ADs_Sam->Draw("e1x0 same");
	//	h_alpha_norm_ADs_Sam->Draw("e1x0 same");
	//	h_alpha_DTnorm_ADs_Sam->Draw("e1x0 same");
	alphaVsAD_AllFits->BuildLegend();

	TCanvas *lambdaVsAD_AllFits = new TCanvas("lambdaVsAD_AllFits","lambdaVsAD_AllFits");
		lambdaVsAD_AllFits->cd();
		h_lambda_rate_ADs->SetStats(0);
		h_lambda_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_lambda_rate_ADs->GetYaxis()->SetTitle("Lambda Value");
	//	h_lambda_rate_ADs->SetLineWidth(0);
	//	h_lambda_norm_ADs->SetLineWidth(0);
	//	h_lambda_DTnorm_ADs->SetLineWidth(0);
	//	h_lambda_rate_ADs_Sam->SetLineWidth(0);
	//	h_lambda_norm_ADs_Sam->SetLineWidth(0);
	//	h_lambda_DTnorm_ADs_Sam->SetLineWidth(0);
		h_lambda_rate_ADs->SetMarkerStyle(21);
		h_lambda_norm_ADs->SetMarkerStyle(20);
		h_lambda_DTnorm_ADs->SetMarkerStyle(22);
		h_lambda_rate_ADs->SetMarkerColor(kRed);
		h_lambda_norm_ADs->SetMarkerColor(kBlue);
		h_lambda_DTnorm_ADs->SetMarkerColor(kBlack);
		h_lambda_rate_ADs->SetMarkerSize(1.5);
		h_lambda_norm_ADs->SetMarkerSize(1.5);
		h_lambda_DTnorm_ADs->SetMarkerSize(1.5);
		h_lambda_rate_ADs_Sam->SetMarkerSize(1.5);
		h_lambda_norm_ADs_Sam->SetMarkerSize(1.5);
		h_lambda_DTnorm_ADs_Sam->SetMarkerSize(1.5);
		h_lambda_rate_ADs_Sam->SetMarkerStyle(25);
		h_lambda_norm_ADs_Sam->SetMarkerStyle(24);
		h_lambda_DTnorm_ADs_Sam->SetMarkerStyle(26);
		h_lambda_rate_ADs_Sam->SetMarkerColor(kMagenta);
		h_lambda_norm_ADs_Sam->SetMarkerColor(kCyan);
		h_lambda_DTnorm_ADs_Sam->SetMarkerColor(kGray+1);
		h_lambda_rate_ADs->SetLineColor(kRed);
		h_lambda_norm_ADs->SetLineColor(kBlue);
		h_lambda_DTnorm_ADs->SetLineColor(kBlack);
		h_lambda_rate_ADs_Sam->SetLineColor(kMagenta);
		h_lambda_norm_ADs_Sam->SetLineColor(kCyan);
		h_lambda_DTnorm_ADs_Sam->SetLineColor(kGray+1);

	//	h_lambda_rate_ADs->GetYaxis()->SetRangeUser(lambdaAvg - 0.005, lambdaAvg + 0.005);

		h_lambda_rate_ADs->Draw("e1x0");
		h_lambda_norm_ADs->Draw("e1x0 same");
		h_lambda_DTnorm_ADs->Draw("e1x0 same");
	//	h_lambda_rate_ADs_Sam->Draw("e1x0 same");
	//	h_lambda_norm_ADs_Sam->Draw("e1x0 same");
	//	h_lambda_DTnorm_ADs_Sam->Draw("e1x0 same");
	lambdaVsAD_AllFits->BuildLegend();

	TCanvas *efficiencyVsAD_AllFits = new TCanvas("efficiencyVsAD_AllFits","efficiencyVsAD_AllFits");
		efficiencyVsAD_AllFits->cd();
		h_efficiency_rate_ADs->SetStats(0);
		h_efficiency_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_efficiency_rate_ADs->GetYaxis()->SetTitle("Delayed Cut Efficiency");
	//	h_efficiency_rate_ADs->SetLineWidth(0);
	//	h_efficiency_norm_ADs->SetLineWidth(0);
	//	h_efficiency_DTnorm_ADs->SetLineWidth(0);
	//	h_efficiency_rate_ADs_Sam->SetLineWidth(0);
	//	h_efficiency_norm_ADs_Sam->SetLineWidth(0);
	//	h_efficiency_DTnorm_ADs_Sam->SetLineWidth(0);
		h_efficiency_rate_ADs->SetMarkerStyle(21);
		h_efficiency_norm_ADs->SetMarkerStyle(20);
		h_efficiency_DTnorm_ADs->SetMarkerStyle(22);
		h_efficiency_rate_ADs->SetMarkerColor(kRed);
		h_efficiency_norm_ADs->SetMarkerColor(kBlue);
		h_efficiency_DTnorm_ADs->SetMarkerColor(kBlack);
		h_efficiency_rate_ADs->SetMarkerSize(1.5);
		h_efficiency_norm_ADs->SetMarkerSize(1.5);
		h_efficiency_DTnorm_ADs->SetMarkerSize(1.5);
		h_efficiency_rate_ADs_Sam->SetMarkerSize(1.5);
		h_efficiency_norm_ADs_Sam->SetMarkerSize(1.5);
		h_efficiency_DTnorm_ADs_Sam->SetMarkerSize(1.5);
		h_efficiency_rate_ADs_Sam->SetMarkerStyle(25);
		h_efficiency_norm_ADs_Sam->SetMarkerStyle(24);
		h_efficiency_DTnorm_ADs_Sam->SetMarkerStyle(26);
		h_efficiency_rate_ADs_Sam->SetMarkerColor(kMagenta);
		h_efficiency_norm_ADs_Sam->SetMarkerColor(kCyan);
		h_efficiency_DTnorm_ADs_Sam->SetMarkerColor(kGray+1);
		h_efficiency_rate_ADs->SetLineColor(kRed);
		h_efficiency_norm_ADs->SetLineColor(kBlue);
		h_efficiency_DTnorm_ADs->SetLineColor(kBlack);
		h_efficiency_rate_ADs_Sam->SetLineColor(kMagenta);
		h_efficiency_norm_ADs_Sam->SetLineColor(kCyan);
		h_efficiency_DTnorm_ADs_Sam->SetLineColor(kGray+1);

		h_efficiency_rate_ADs->GetYaxis()->SetRangeUser(efficiencyAvg - 0.015, efficiencyAvg + 0.015);

		h_efficiency_rate_ADs->Draw("e1x0 ");
		h_efficiency_norm_ADs->Draw("e1x0 same");
		h_efficiency_DTnorm_ADs->Draw("e1x0 same");
	//	h_efficiency_rate_ADs_Sam->Draw("e1x0 same");
	//	h_efficiency_norm_ADs_Sam->Draw("e1x0 same");
	//	h_efficiency_DTnorm_ADs_Sam->Draw("e1x0 same");
	efficiencyVsAD_AllFits->BuildLegend();

	TF1* nomVsExt_Fit_rate = new TF1("nomVsExt_Fit_rate", "[0]+[1]*x",0,500000);
		if(a == 0) nomVsExt_Fit_rate->FixParameter(0,0);
		nomVsExt_Fit_rate->SetParameter(1,0.95);
	TF1* nomVsExt_Fit_norm = new TF1("nomVsExt_Fit_norm", "[0]+[1]*x",0,500000);
		if(a == 0) nomVsExt_Fit_norm->FixParameter(0,0);
		nomVsExt_Fit_norm->SetParameter(1,0.95);
	TF1* nomVsExt_Fit_DTnorm = new TF1("nomVsExt_Fit_DTnorm", "[0]+[1]*x",0,500000);
		if(a == 0) nomVsExt_Fit_DTnorm->FixParameter(0,0);
		nomVsExt_Fit_DTnorm->SetParameter(1,0.95);

	TCanvas* c1 = new TCanvas("c1","Sam/Jinjing Plot - Rate");
	c1->cd();
	const Int_t n = 8;
	Double_t e_nom_rate[n];
	Double_t e_ext_rate[n];
	for(int i = 0; i<n; i++){
		e_nom_rate[i] = sqrt(N_nom_rate[i]);
		e_ext_rate[i] = sqrt(N_ext_rate[i]);
	}
	auto gr_rate = new TGraphErrors(n,N_ext_rate,N_nom_rate,e_ext_rate,e_nom_rate);
	gr_rate->SetMarkerColor(4);
	gr_rate->SetMarkerStyle(21);
	gr_rate->Draw("AP");
		gr_rate->Fit("nomVsExt_Fit_rate", "R");

	TCanvas* c2 = new TCanvas("c2","Sam/Jinjing Plot - Normalized");
	c2->cd();
	Double_t e_nom_norm[n];
	Double_t e_ext_norm[n];
	for(int i = 0; i<n; i++){
		e_nom_norm[i] = sqrt(N_nom_norm[i]);
		e_ext_norm[i] = sqrt(N_ext_norm[i]);
	}
	auto gr_norm = new TGraphErrors(n,N_ext_norm,N_nom_norm,e_ext_norm,e_nom_norm);
	gr_norm->SetMarkerColor(4);
	gr_norm->SetMarkerStyle(21);
	gr_norm->Draw("AP");
		gr_norm->Fit("nomVsExt_Fit_norm", "R");

	TCanvas* c3 = new TCanvas("c3","Sam/Jinjing Plot - DT Normalized");
	c3->cd();
	Double_t e_nom_DTnorm[n];
	Double_t e_ext_DTnorm[n];
	for(int i = 0; i<n; i++){
		e_nom_DTnorm[i] = sqrt(N_nom_DTnorm[i]);
		e_ext_DTnorm[i] = sqrt(N_ext_DTnorm[i]);
	}
	auto gr_DTnorm = new TGraphErrors(n,N_ext_DTnorm,N_nom_DTnorm,e_ext_DTnorm,e_nom_DTnorm);
	gr_DTnorm->SetMarkerColor(4);
	gr_DTnorm->SetMarkerStyle(21);
	gr_DTnorm->Draw("AP");
		gr_DTnorm->Fit("nomVsExt_Fit_DTnorm", "R");

/*	TF1* nomVsExt_Fit_rate_Sam = new TF1("nomVsExt_Fit_rate_Sam", "[0]+[1]*x",0,500000);
		if(a == 0) nomVsExt_Fit_rate_Sam->FixParameter(0,0);
		nomVsExt_Fit_rate_Sam->SetParameter(1,0.95);
	TF1* nomVsExt_Fit_norm_Sam = new TF1("nomVsExt_Fit_norm_Sam", "[0]+[1]*x",0,500000);
		if(a == 0) nomVsExt_Fit_norm_Sam->FixParameter(0,0);
		nomVsExt_Fit_norm_Sam->SetParameter(1,0.95);
	TF1* nomVsExt_Fit_DTnorm_Sam = new TF1("nomVsExt_Fit_DTnorm_Sam", "[0]+[1]*x",0,500000);
		if(a == 0) nomVsExt_Fit_DTnorm_Sam->FixParameter(0,0);
		nomVsExt_Fit_DTnorm_Sam->SetParameter(1,0.95);

	TCanvas* c4 = new TCanvas("c4","Sam/Jinjing Plot - Rate");
	c4->cd();
	Double_t e_nom_rate_Sam[n];
	Double_t e_ext_rate_Sam[n];
	for(int i = 0; i<n; i++){
		e_nom_rate_Sam[i] = sqrt(N_nom_rate_Sam[i]);
		e_ext_rate_Sam[i] = sqrt(N_ext_rate_Sam[i]);
	}
	auto gr_rate_Sam = new TGraphErrors(n,N_ext_rate_Sam,N_nom_rate_Sam,e_ext_rate_Sam,e_nom_rate_Sam);
	gr_rate_Sam->SetMarkerColor(4);
	gr_rate_Sam->SetMarkerStyle(21);
	gr_rate_Sam->Draw("AP");
		gr_rate_Sam->Fit("nomVsExt_Fit_rate_Sam", "R");

	TCanvas* c5 = new TCanvas("c5","Sam/Jinjing Plot - Normalized");
	c5->cd();
	Double_t e_nom_norm_Sam[n];
	Double_t e_ext_norm_Sam[n];
	for(int i = 0; i<n; i++){
		e_nom_norm_Sam[i] = sqrt(N_nom_norm_Sam[i]);
		e_ext_norm_Sam[i] = sqrt(N_ext_norm_Sam[i]);
	}
	auto gr_norm_Sam = new TGraphErrors(n,N_ext_norm_Sam,N_nom_norm_Sam,e_ext_norm_Sam,e_nom_norm_Sam);
	gr_norm_Sam->SetMarkerColor(4);
	gr_norm_Sam->SetMarkerStyle(21);
	gr_norm_Sam->Draw("AP");
		gr_norm_Sam->Fit("nomVsExt_Fit_norm_Sam", "R");

	TCanvas* c6 = new TCanvas("c6","Sam/Jinjing Plot - DT Normalized");
	c6->cd();
	Double_t e_nom_DTnorm_Sam[n];
	Double_t e_ext_DTnorm_Sam[n];
	for(int i = 0; i<n; i++){
		e_nom_DTnorm_Sam[i] = sqrt(N_nom_DTnorm_Sam[i]);
		e_ext_DTnorm_Sam[i] = sqrt(N_ext_DTnorm_Sam[i]);
	}
	auto gr_DTnorm_Sam = new TGraphErrors(n,N_ext_DTnorm_Sam,N_nom_DTnorm_Sam,e_ext_DTnorm_Sam,e_nom_DTnorm_Sam);
	gr_DTnorm_Sam->SetMarkerColor(4);
	gr_DTnorm_Sam->SetMarkerStyle(21);
	gr_DTnorm_Sam->Draw("AP");
		gr_DTnorm_Sam->Fit("nomVsExt_Fit_DTnorm_Sam", "R");*/

	TH1D* h_absDif_rate_ADs=new TH1D("h_absDif_rate_ADs","h_absDif_rate_ADs",8,0.4,8.4);
	TH1D* h_absDif_norm_ADs=new TH1D("h_absDif_norm_ADs","h_absDif_norm_ADs",8,0.5,8.5);
	TH1D* h_absDif_DTnorm_ADs=new TH1D("h_absDif_DTnorm_ADs","h_absDif_DTnorm_ADs",8,0.6,8.6);

	TH1D* h_relDif_rate_ADs=new TH1D("h_relDif_rate_ADs","h_relDif_rate_ADs",8,0.4,8.4);
	TH1D* h_relDif_norm_ADs=new TH1D("h_relDif_norm_ADs","h_relDif_norm_ADs",8,0.5,8.5);
	TH1D* h_relDif_DTnorm_ADs=new TH1D("h_relDif_DTnorm_ADs","h_relDif_DTnorm_ADs",8,0.6,8.6);

/*	TH1D* h_absDif_rate_Sam_ADs=new TH1D("h_absDif_rate_Sam_ADs","h_absDif_rate_Sam_ADs",8,0.45,8.45);
	TH1D* h_absDif_norm_Sam_ADs=new TH1D("h_absDif_norm_Sam_ADs","h_absDif_norm_Sam_ADs",8,0.55,8.55);
	TH1D* h_absDif_DTnorm_Sam_ADs=new TH1D("h_absDif_DTnorm_Sam_ADs","h_absDif_DTnorm_Sam_ADs",8,0.65,8.65);

	TH1D* h_relDif_rate_Sam_ADs=new TH1D("h_relDif_rate_Sam_ADs","h_relDif_rate_Sam_ADs",8,0.45,8.45);
	TH1D* h_relDif_norm_Sam_ADs=new TH1D("h_relDif_norm_Sam_ADs","h_relDif_norm_Sam_ADs",8,0.55,8.55);
	TH1D* h_relDif_DTnorm_Sam_ADs=new TH1D("h_relDif_DTnorm_Sam_ADs","h_relDif_DTnorm_Sam_ADs",8,0.65,8.65);*/

	double N_pred_rate[8];
	double N_pred_norm[8];
	double N_pred_DTnorm[8];
	double relDif_rate[8];
	double relDif_norm[8];
	double relDif_DTnorm[8];

	double N_pred_rate_Sam[8];
	double N_pred_norm_Sam[8];
	double N_pred_DTnorm_Sam[8];
	double relDif_rate_Sam[8];
	double relDif_norm_Sam[8];
	double relDif_DTnorm_Sam[8];
		for(int i=0; i<8; i++){
			N_pred_rate[i] = 0;
			N_pred_norm[i] = 0;
			N_pred_DTnorm[i] = 0;

			relDif_rate[i] = 0;
			relDif_norm[i] = 0;
			relDif_DTnorm[i] = 0;
			N_pred_rate[i] = nomVsExt_Fit_rate->GetParameter(0) + nomVsExt_Fit_rate->GetParameter(1)*N_ext_rate[i];
			N_pred_norm[i] = nomVsExt_Fit_norm->GetParameter(0) + nomVsExt_Fit_norm->GetParameter(1)*N_ext_norm[i];
			N_pred_DTnorm[i] = nomVsExt_Fit_DTnorm->GetParameter(0) + nomVsExt_Fit_DTnorm->GetParameter(1)*N_ext_DTnorm[i];

			relDif_rate[i] = 100*(N_nom_rate[i]-N_pred_rate[i])/N_pred_rate[i];
			relDif_norm[i] = 100*(N_nom_norm[i]-N_pred_norm[i])/N_pred_norm[i];
			relDif_DTnorm[i] = 100*(N_nom_DTnorm[i]-N_pred_DTnorm[i])/N_pred_DTnorm[i];

/*			N_pred_rate_Sam[i] = 0;
			N_pred_norm_Sam[i] = 0;
			N_pred_DTnorm_Sam[i] = 0;

			relDif_rate_Sam[i] = 0;
			relDif_norm_Sam[i] = 0;
			relDif_DTnorm_Sam[i] = 0;
			N_pred_rate_Sam[i] = nomVsExt_Fit_rate_Sam->GetParameter(0) + nomVsExt_Fit_rate_Sam->GetParameter(1)*N_ext_rate_Sam[i];
			N_pred_norm_Sam[i] = nomVsExt_Fit_norm_Sam->GetParameter(0) + nomVsExt_Fit_norm_Sam->GetParameter(1)*N_ext_norm_Sam[i];
			N_pred_DTnorm_Sam[i] = nomVsExt_Fit_DTnorm_Sam->GetParameter(0) + nomVsExt_Fit_DTnorm_Sam->GetParameter(1)*N_ext_DTnorm_Sam[i];

			relDif_rate_Sam[i] = 100*(N_nom_rate_Sam[i]-N_pred_rate_Sam[i])/N_pred_rate_Sam[i];
			relDif_norm_Sam[i] = 100*(N_nom_norm_Sam[i]-N_pred_norm_Sam[i])/N_pred_norm_Sam[i];
			relDif_DTnorm_Sam[i] = 100*(N_nom_DTnorm_Sam[i]-N_pred_DTnorm_Sam[i])/N_pred_DTnorm_Sam[i];*/
		}

	cout << "RelDif (%)\tRate\tNorm\tDTnorm" << endl;
	for(int i=0; i<8; i++){
		//cout << i+1 << "\t" << N_pred_rate[i] << "\t" << N_pred_norm[i] << "\t" << N_pred_DTnorm[i] << endl;
		cout << i+1 << "\t\t" << relDif_rate[i] << "\t" << relDif_norm[i] << "\t" << relDif_DTnorm[i] << endl;

		h_absDif_rate_ADs->Fill(i+1, N_nom_rate[i]-N_pred_rate[i]);
		h_absDif_norm_ADs->Fill(i+1, N_nom_norm[i]-N_pred_norm[i]);
		h_absDif_DTnorm_ADs->Fill(i+1, N_nom_DTnorm[i]-N_pred_DTnorm[i]);

		h_relDif_rate_ADs->Fill(i+1, relDif_rate[i]);
		h_relDif_norm_ADs->Fill(i+1, relDif_norm[i]);
		h_relDif_DTnorm_ADs->Fill(i+1, relDif_DTnorm[i]);

	/*	h_absDif_rate_Sam_ADs->Fill(i+1, N_nom_rate_Sam[i]-N_pred_rate_Sam[i]);
		h_absDif_norm_Sam_ADs->Fill(i+1, N_nom_norm_Sam[i]-N_pred_norm_Sam[i]);
		h_absDif_DTnorm_Sam_ADs->Fill(i+1, N_nom_DTnorm_Sam[i]-N_pred_DTnorm_Sam[i]);

		h_relDif_rate_Sam_ADs->Fill(i+1, relDif_rate_Sam[i]);
		h_relDif_norm_Sam_ADs->Fill(i+1, relDif_norm_Sam[i]);
		h_relDif_DTnorm_Sam_ADs->Fill(i+1, relDif_DTnorm_Sam[i]);*/
	}
	
	double temp_eff_low = 0;
	double temp_eff_high = 0;
	cout << "Delayed energy cut efficiencies:" << endl;
	for(int i=0; i<8; i++){
		if(i==0){
			temp_eff_low = efficiency_rate[i];
			temp_eff_high = efficiency_rate[i];
		}
		cout << i+1 << "\t" << efficiency_rate[i] << endl;
		if(efficiency_rate[i] < temp_eff_low) temp_eff_low = efficiency_rate[i];
		if(efficiency_rate[i] > temp_eff_high) temp_eff_high = efficiency_rate[i];
	}
	cout << "Uncertainty = " << (temp_eff_high - temp_eff_low)/2. << endl;
	

	TCanvas *absDifVsAD_AllFits = new TCanvas("absDifVsAD_AllFits","absDifVsAD_AllFits");
		absDifVsAD_AllFits->cd();
		h_absDif_rate_ADs->SetStats(0);
		h_absDif_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_absDif_rate_ADs->GetYaxis()->SetTitle("Absolute Difference from Fit");
		h_absDif_rate_ADs->SetMarkerStyle(21);
		h_absDif_norm_ADs->SetMarkerStyle(20);
		h_absDif_DTnorm_ADs->SetMarkerStyle(22);
		h_absDif_rate_ADs->SetMarkerColor(kRed);
		h_absDif_norm_ADs->SetMarkerColor(kBlue);
		h_absDif_DTnorm_ADs->SetMarkerColor(kBlack);
		h_absDif_rate_ADs->SetMarkerSize(1.5);
		h_absDif_norm_ADs->SetMarkerSize(1.5);
		h_absDif_DTnorm_ADs->SetMarkerSize(1.5);
	//	h_absDif_rate_Sam_ADs->SetMarkerSize(1.5);
	//	h_absDif_norm_Sam_ADs->SetMarkerSize(1.5);
	//	h_absDif_DTnorm_Sam_ADs->SetMarkerSize(1.5);
	//	h_absDif_rate_Sam_ADs->SetMarkerStyle(25);
	//	h_absDif_norm_Sam_ADs->SetMarkerStyle(24);
	//	h_absDif_DTnorm_Sam_ADs->SetMarkerStyle(26);
	//	h_absDif_rate_Sam_ADs->SetMarkerColor(kMagenta);
	//	h_absDif_norm_Sam_ADs->SetMarkerColor(kCyan);
	//	h_absDif_DTnorm_Sam_ADs->SetMarkerColor(kGray+1);
		h_absDif_rate_ADs->SetLineColor(kRed);
		h_absDif_norm_ADs->SetLineColor(kBlue);
		h_absDif_DTnorm_ADs->SetLineColor(kBlack);
	//	h_absDif_rate_Sam_ADs->SetLineColor(kMagenta);
	//	h_absDif_norm_Sam_ADs->SetLineColor(kCyan);
	//	h_absDif_DTnorm_Sam_ADs->SetLineColor(kGray+1);

		h_absDif_rate_ADs->SetLineWidth(0);
		h_absDif_norm_ADs->SetLineWidth(0);
		h_absDif_DTnorm_ADs->SetLineWidth(0);
	//	h_absDif_rate_Sam_ADs->SetLineWidth(0);
	//	h_absDif_norm_Sam_ADs->SetLineWidth(0);
	//	h_absDif_DTnorm_Sam_ADs->SetLineWidth(0);

		h_absDif_rate_ADs->Draw("e1x0 ");
		h_absDif_norm_ADs->Draw("e1x0 same");
		h_absDif_DTnorm_ADs->Draw("e1x0 same");
	//	h_absDif_rate_Sam_ADs->Draw("e1x0 same");
	//	h_absDif_norm_Sam_ADs->Draw("e1x0 same");
	//	h_absDif_DTnorm_Sam_ADs->Draw("e1x0 same");
	absDifVsAD_AllFits->BuildLegend();

	TCanvas *relDifVsAD_AllFits = new TCanvas("relDifVsAD_AllFits","relDifVsAD_AllFits");
		relDifVsAD_AllFits->cd();
		h_relDif_rate_ADs->SetStats(0);
		h_relDif_rate_ADs->GetXaxis()->SetTitle("AD Number");
		h_relDif_rate_ADs->GetYaxis()->SetTitle("Relative Difference from Fit");
		h_relDif_rate_ADs->SetMarkerStyle(21);
		h_relDif_norm_ADs->SetMarkerStyle(20);
		h_relDif_DTnorm_ADs->SetMarkerStyle(22);
		h_relDif_rate_ADs->SetMarkerColor(kRed);
		h_relDif_norm_ADs->SetMarkerColor(kBlue);
		h_relDif_DTnorm_ADs->SetMarkerColor(kBlack);
		h_relDif_rate_ADs->SetMarkerSize(1.5);
		h_relDif_norm_ADs->SetMarkerSize(1.5);
		h_relDif_DTnorm_ADs->SetMarkerSize(1.5);
//		h_relDif_rate_Sam_ADs->SetMarkerSize(1.5);
//		h_relDif_norm_Sam_ADs->SetMarkerSize(1.5);
//		h_relDif_DTnorm_Sam_ADs->SetMarkerSize(1.5);
//		h_relDif_rate_Sam_ADs->SetMarkerStyle(25);
//		h_relDif_norm_Sam_ADs->SetMarkerStyle(24);
//		h_relDif_DTnorm_Sam_ADs->SetMarkerStyle(26);
//		h_relDif_rate_Sam_ADs->SetMarkerColor(kMagenta);
//		h_relDif_norm_Sam_ADs->SetMarkerColor(kCyan);
//		h_relDif_DTnorm_Sam_ADs->SetMarkerColor(kGray+1);
		h_relDif_rate_ADs->SetLineColor(kRed);
		h_relDif_norm_ADs->SetLineColor(kBlue);
		h_relDif_DTnorm_ADs->SetLineColor(kBlack);
//		h_relDif_rate_Sam_ADs->SetLineColor(kMagenta);
//		h_relDif_norm_Sam_ADs->SetLineColor(kCyan);
//		h_relDif_DTnorm_Sam_ADs->SetLineColor(kGray+1);

		h_relDif_rate_ADs->SetLineWidth(0);
		h_relDif_norm_ADs->SetLineWidth(0);
		h_relDif_DTnorm_ADs->SetLineWidth(0);
//		h_relDif_rate_Sam_ADs->SetLineWidth(0);
//		h_relDif_norm_Sam_ADs->SetLineWidth(0);
//		h_relDif_DTnorm_Sam_ADs->SetLineWidth(0);

		h_relDif_rate_ADs->Draw("e1x0 ");
		h_relDif_norm_ADs->Draw("e1x0 same");
		h_relDif_DTnorm_ADs->Draw("e1x0 same");
	//	h_relDif_rate_Sam_ADs->Draw("e1x0 same");
	//	h_relDif_norm_Sam_ADs->Draw("e1x0 same");
	//	h_relDif_DTnorm_Sam_ADs->Draw("e1x0 same");
	relDifVsAD_AllFits->BuildLegend();

	int colors_AD[8] = {kBlack, kRed, kGreen, kBlue, kMagenta, kOrange, kCyan, kGray+1};

/*	TH1D* h_efficiency_z_rate[8];
	for(int iad = 0; iad < 8; iad++){
		h_efficiency_z_rate[iad]=new TH1D(Form("h_effiency_z_rate_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_z_rate_eh%dad%d",EH[iad],AD[iad]),NzPoints,0.42+iad*.02,0.42+iad*.02+NzPoints);
		for(int iz = 0; iz < NzPoints; iz++){ //filling histograms
			h_efficiency_z_rate[iad]->SetBinContent(h_efficiency_z_rate[iad]->FindBin(iz+1), efficiency_rate_z[iad][iz]);
			h_efficiency_z_rate[iad]->SetBinError(h_efficiency_z_rate[iad]->FindBin(iz+1), efficiencyError_rate_z[iad][iz]);
			//if(iad == 0) cout << efficiencyError_rate_z[iad][iz] << endl;
		}
		h_efficiency_z_rate[iad]->SetMarkerStyle(21);
		h_efficiency_z_rate[iad]->SetLineWidth(1);
		h_efficiency_z_rate[iad]->SetMarkerSize(1.5);
		h_efficiency_z_rate[iad]->SetLineColor(colors_AD[iad]);
		h_efficiency_z_rate[iad]->SetMarkerColor(colors_AD[iad]);
		h_efficiency_z_rate[iad]->GetXaxis()->SetTitle("z slice");
		h_efficiency_z_rate[iad]->GetYaxis()->SetTitle("Efficiency");
		h_efficiency_z_rate[iad]->GetYaxis()->SetRangeUser(0.87,0.99);
		h_efficiency_z_rate[iad]->SetStats(0);
	}

	TCanvas *z_slice_effs_rate = new TCanvas("z_slice_effs_rate","z_slice_effs_rate");
		z_slice_effs_rate->cd();
		h_efficiency_z_rate[0]->Draw("e1x0");
		h_efficiency_z_rate[0]->GetXaxis()->SetTitle("z slice");
		h_efficiency_z_rate[0]->GetYaxis()->SetTitle("Efficiency");
		for(int i = 1; i < 8; i++){
			h_efficiency_z_rate[i]->Draw("e1x0 same");
		}
	z_slice_effs_rate->SetGridy();
	z_slice_effs_rate->BuildLegend();

	TH1D* h_efficiency_z_norm[8];
	for(int iad = 0; iad < 8; iad++){
		h_efficiency_z_norm[iad]=new TH1D(Form("h_effiency_z_norm_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_z_norm_eh%dad%d",EH[iad],AD[iad]),NzPoints,0.42+iad*.02,0.42+iad*.02+NzPoints);
		for(int iz = 0; iz < NzPoints; iz++){ //filling histograms
			h_efficiency_z_norm[iad]->SetBinContent(h_efficiency_z_norm[iad]->FindBin(iz+1), efficiency_norm_z[iad][iz]);
			h_efficiency_z_norm[iad]->SetBinError(h_efficiency_z_norm[iad]->FindBin(iz+1), efficiencyError_norm_z[iad][iz]);
			//if(iad == 0) cout << efficiencyError_norm_z[iad][iz] << endl;
		}
		h_efficiency_z_norm[iad]->SetMarkerStyle(21);
		h_efficiency_z_norm[iad]->SetLineWidth(1);
		h_efficiency_z_norm[iad]->SetMarkerSize(1.5);
		h_efficiency_z_norm[iad]->SetLineColor(colors_AD[iad]);
		h_efficiency_z_norm[iad]->SetMarkerColor(colors_AD[iad]);
		h_efficiency_z_norm[iad]->GetXaxis()->SetTitle("z slice");
		h_efficiency_z_norm[iad]->GetYaxis()->SetTitle("Efficiency");
		h_efficiency_z_norm[iad]->GetYaxis()->SetRangeUser(0.87,0.99);
		h_efficiency_z_norm[iad]->SetStats(0);
	}

	TCanvas *z_slice_effs_norm = new TCanvas("z_slice_effs_norm","z_slice_effs_norm");
		z_slice_effs_norm->cd();
		h_efficiency_z_norm[0]->Draw("e1x0");
		h_efficiency_z_norm[0]->GetXaxis()->SetTitle("z slice");
		h_efficiency_z_norm[0]->GetYaxis()->SetTitle("Efficiency");
		for(int i = 1; i < 8; i++){
			h_efficiency_z_norm[i]->Draw("e1x0 same");
		}
	z_slice_effs_norm->SetGridy();
	z_slice_effs_norm->BuildLegend();

	TH1D* h_efficiency_z_DTnorm[8];
	for(int iad = 0; iad < 8; iad++){
		h_efficiency_z_DTnorm[iad]=new TH1D(Form("h_effiency_z_DTnorm_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_z_DTnorm_eh%dad%d",EH[iad],AD[iad]),NzPoints,0.42+iad*.02,0.42+iad*.02+NzPoints);
		for(int iz = 0; iz < NzPoints; iz++){ //filling histograms
			h_efficiency_z_DTnorm[iad]->SetBinContent(h_efficiency_z_DTnorm[iad]->FindBin(iz+1), efficiency_DTnorm_z[iad][iz]);
			h_efficiency_z_DTnorm[iad]->SetBinError(h_efficiency_z_DTnorm[iad]->FindBin(iz+1), efficiencyError_DTnorm_z[iad][iz]);
			//if(iad == 0) cout << efficiencyError_DTnorm_z[iad][iz] << endl;
		}
		h_efficiency_z_DTnorm[iad]->SetMarkerStyle(21);
		h_efficiency_z_DTnorm[iad]->SetLineWidth(1);
		h_efficiency_z_DTnorm[iad]->SetMarkerSize(1.5);
		h_efficiency_z_DTnorm[iad]->SetLineColor(colors_AD[iad]);
		h_efficiency_z_DTnorm[iad]->SetMarkerColor(colors_AD[iad]);
		h_efficiency_z_DTnorm[iad]->GetXaxis()->SetTitle("z slice");
		h_efficiency_z_DTnorm[iad]->GetYaxis()->SetTitle("Efficiency");
		h_efficiency_z_DTnorm[iad]->GetYaxis()->SetRangeUser(0.87,0.99);
		h_efficiency_z_DTnorm[iad]->SetStats(0);
	}

	TCanvas *z_slice_effs_DTnorm = new TCanvas("z_slice_effs_DTnorm","z_slice_effs_DTnorm");
		z_slice_effs_DTnorm->cd();
		h_efficiency_z_DTnorm[0]->Draw("e1x0");
		h_efficiency_z_DTnorm[0]->GetXaxis()->SetTitle("z slice");
		h_efficiency_z_DTnorm[0]->GetYaxis()->SetTitle("Efficiency");
		for(int i = 1; i < 8; i++){
			h_efficiency_z_DTnorm[i]->Draw("e1x0 same");
		}
	z_slice_effs_DTnorm->SetGridy();
	z_slice_effs_DTnorm->BuildLegend();

	for(int iz = 0; iz < NzPoints; iz++){
		TCanvas *z_slice_compADs_rate = new TCanvas(Form("z_slice_compADs_rate_iz%d",iz+1),Form("z_slice_compADs_rate_iz%d",iz+1));
		z_slice_compADs_rate->cd();
		for(int iad = 0; iad < 8; iad++){
			h_Edelayed_sub_rate_z[iad][iz]->Rebin(2);
			h_Edelayed_sub_rate_z[iad][iz]->SetStats(1111);
			h_Edelayed_sub_rate_z[iad][iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_rate_z[iad][iz]->GetXaxis()->SetRangeUser(1.5,3);
			h_Edelayed_sub_rate_z[iad][iz]->GetYaxis()->SetTitle("Counts");
			//h_Edelayed_sub_rate_z[iad][iz]->GetYaxis()->SetRangeUser(-50,3000);
			h_Edelayed_sub_rate_z[iad][iz]->SetLineColor(colors_AD[iad]);
			h_Edelayed_sub_rate_z[iad][iz]->SetMarkerColor(colors_AD[iad]);
		}
		h_Edelayed_sub_rate_z[2][iz]->Draw();
		for(int iad = 0; iad < 8; iad++){
			if(iad == 2) continue;
			else h_Edelayed_sub_rate_z[iad][iz]->Draw("same");
		}
		z_slice_compADs_rate->SetGridy();
		//z_slice_compADs_rate->BuildLegend();
		z_slice_compADs_rate->Print(Form("../nH_files/z_slice_compADs_rate_iz%d.png",iz+1));
	}

	TH1D* h_efficiency_r2_rate[8];
//	int colors_AD[8] = {kBlack, kRed, kGreen, kBlue, kMagenta, kOrange, kCyan, kGray+1};
	for(int iad = 0; iad < 8; iad++){
		h_efficiency_r2_rate[iad]=new TH1D(Form("h_effiency_r2_rate_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_r2_rate_eh%dad%d",EH[iad],AD[iad]),Nr2Points,0.42+iad*.02,0.42+iad*.02+Nr2Points);
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){ //filling histograms
			h_efficiency_r2_rate[iad]->SetBinContent(h_efficiency_r2_rate[iad]->FindBin(ir2+1), efficiency_rate_r2[iad][ir2]);
			h_efficiency_r2_rate[iad]->SetBinError(h_efficiency_r2_rate[iad]->FindBin(ir2+1), efficiencyError_rate_r2[iad][ir2]);
			//if(iad == 0) cout << efficiencyError_rate_r2[iad][ir2] << endl;
		}
		h_efficiency_r2_rate[iad]->SetMarkerStyle(21);
		h_efficiency_r2_rate[iad]->SetLineWidth(1);
		h_efficiency_r2_rate[iad]->SetMarkerSize(1.5);
		h_efficiency_r2_rate[iad]->SetLineColor(colors_AD[iad]);
		h_efficiency_r2_rate[iad]->SetMarkerColor(colors_AD[iad]);
		h_efficiency_r2_rate[iad]->GetXaxis()->SetTitle("r2 slice");
		h_efficiency_r2_rate[iad]->GetYaxis()->SetTitle("Efficiency");
		h_efficiency_r2_rate[iad]->GetYaxis()->SetRangeUser(0.87,0.985);
		h_efficiency_r2_rate[iad]->SetStats(0);
	}

	TCanvas *r2_slice_effs_rate = new TCanvas("r2_slice_effs_rate","r2_slice_effs_rate");
		r2_slice_effs_rate->cd();
		h_efficiency_r2_rate[0]->Draw("e1x0");
		h_efficiency_r2_rate[0]->GetXaxis()->SetTitle("r2 slice");
		h_efficiency_r2_rate[0]->GetYaxis()->SetTitle("Efficiency");
		for(int i = 1; i < 8; i++){
			h_efficiency_r2_rate[i]->Draw("e1x0 same");
		}
	r2_slice_effs_rate->SetGridy();
	r2_slice_effs_rate->BuildLegend();

	TH1D* h_efficiency_r2_norm[8];
	for(int iad = 0; iad < 8; iad++){
		h_efficiency_r2_norm[iad]=new TH1D(Form("h_effiency_r2_norm_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_r2_norm_eh%dad%d",EH[iad],AD[iad]),Nr2Points,0.42+iad*.02,0.42+iad*.02+Nr2Points);
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){ //filling histograms
			h_efficiency_r2_norm[iad]->SetBinContent(h_efficiency_r2_norm[iad]->FindBin(ir2+1), efficiency_norm_r2[iad][ir2]);
			h_efficiency_r2_norm[iad]->SetBinError(h_efficiency_r2_norm[iad]->FindBin(ir2+1), efficiencyError_norm_r2[iad][ir2]);
			//if(iad == 0) cout << efficiencyError_norm_r2[iad][ir2] << endl;
		}
		h_efficiency_r2_norm[iad]->SetMarkerStyle(21);
		h_efficiency_r2_norm[iad]->SetLineWidth(1);
		h_efficiency_r2_norm[iad]->SetMarkerSize(1.5);
		h_efficiency_r2_norm[iad]->SetLineColor(colors_AD[iad]);
		h_efficiency_r2_norm[iad]->SetMarkerColor(colors_AD[iad]);
		h_efficiency_r2_norm[iad]->GetXaxis()->SetTitle("r2 slice");
		h_efficiency_r2_norm[iad]->GetYaxis()->SetTitle("Efficiency");
		h_efficiency_r2_norm[iad]->GetYaxis()->SetRangeUser(0.87,0.985);
		h_efficiency_r2_norm[iad]->SetStats(0);
	}

	TCanvas *r2_slice_effs_norm = new TCanvas("r2_slice_effs_norm","r2_slice_effs_norm");
		r2_slice_effs_norm->cd();
		h_efficiency_r2_norm[0]->Draw("e1x0");
		h_efficiency_r2_norm[0]->GetXaxis()->SetTitle("r2 slice");
		h_efficiency_r2_norm[0]->GetYaxis()->SetTitle("Efficiency");
		for(int i = 1; i < 8; i++){
			h_efficiency_r2_norm[i]->Draw("e1x0 same");
		}
	r2_slice_effs_norm->SetGridy();
	r2_slice_effs_norm->BuildLegend();

	TH1D* h_efficiency_r2_DTnorm[8];
	for(int iad = 0; iad < 8; iad++){
		h_efficiency_r2_DTnorm[iad]=new TH1D(Form("h_effiency_r2_DTnorm_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_r2_DTnorm_eh%dad%d",EH[iad],AD[iad]),Nr2Points,0.42+iad*.02,0.42+iad*.02+Nr2Points);
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){ //filling histograms
			h_efficiency_r2_DTnorm[iad]->SetBinContent(h_efficiency_r2_DTnorm[iad]->FindBin(ir2+1), efficiency_DTnorm_r2[iad][ir2]);
			h_efficiency_r2_DTnorm[iad]->SetBinError(h_efficiency_r2_DTnorm[iad]->FindBin(ir2+1), efficiencyError_DTnorm_r2[iad][ir2]);
			//if(iad == 0) cout << efficiencyError_DTnorm_r2[iad][ir2] << endl;
		}
		h_efficiency_r2_DTnorm[iad]->SetMarkerStyle(21);
		h_efficiency_r2_DTnorm[iad]->SetLineWidth(1);
		h_efficiency_r2_DTnorm[iad]->SetMarkerSize(1.5);
		h_efficiency_r2_DTnorm[iad]->SetLineColor(colors_AD[iad]);
		h_efficiency_r2_DTnorm[iad]->SetMarkerColor(colors_AD[iad]);
		h_efficiency_r2_DTnorm[iad]->GetXaxis()->SetTitle("r2 slice");
		h_efficiency_r2_DTnorm[iad]->GetYaxis()->SetTitle("Efficiency");
		h_efficiency_r2_DTnorm[iad]->GetYaxis()->SetRangeUser(0.87,0.985);
		h_efficiency_r2_DTnorm[iad]->SetStats(0);
	}

	TCanvas *r2_slice_effs_DTnorm = new TCanvas("r2_slice_effs_DTnorm","r2_slice_effs_DTnorm");
		r2_slice_effs_DTnorm->cd();
		h_efficiency_r2_DTnorm[0]->Draw("e1x0");
		h_efficiency_r2_DTnorm[0]->GetXaxis()->SetTitle("r2 slice");
		h_efficiency_r2_DTnorm[0]->GetYaxis()->SetTitle("Efficiency");
		for(int i = 1; i < 8; i++){
			h_efficiency_r2_DTnorm[i]->Draw("e1x0 same");
		}
	r2_slice_effs_DTnorm->SetGridy();
	r2_slice_effs_DTnorm->BuildLegend();

	for(int ir2 = 0; ir2 < Nr2Points; ir2++){
		TCanvas *r2_slice_compADs_rate = new TCanvas(Form("r2_slice_compADs_rate_ir2%d",ir2+1),Form("r2_slice_compADs_rate_ir2%d",ir2+1));
		r2_slice_compADs_rate->cd();
		for(int iad = 0; iad < 8; iad++){
			h_Edelayed_sub_rate_r2[iad][ir2]->Rebin(2);
			h_Edelayed_sub_rate_r2[iad][ir2]->SetStats(1111);
			h_Edelayed_sub_rate_r2[iad][ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
			h_Edelayed_sub_rate_r2[iad][ir2]->GetXaxis()->SetRangeUser(1.5,3);
			h_Edelayed_sub_rate_r2[iad][ir2]->GetYaxis()->SetTitle("Counts");
			//h_Edelayed_sub_rate_r2[iad][ir2]->GetYaxis()->SetRangeUser(-50,3000);
			h_Edelayed_sub_rate_r2[iad][ir2]->SetLineColor(colors_AD[iad]);
			h_Edelayed_sub_rate_r2[iad][ir2]->SetMarkerColor(colors_AD[iad]);
		}
		h_Edelayed_sub_rate_r2[2][ir2]->Draw();
		for(int iad = 0; iad < 8; iad++){
			if(iad == 2) continue;
			else h_Edelayed_sub_rate_r2[iad][ir2]->Draw("same");
		}
		r2_slice_compADs_rate->SetGridy();
		//r2_slice_compADs_rate->BuildLegend();
		r2_slice_compADs_rate->Print(Form("../nH_files/r2_slice_compADs_rate_ir2%d.png",ir2+1));
	}


//zVSr2 efficiency comparisons
	TH1D* h_efficiency_zVSr2_rate[8];
//	int colors_AD[8] = {kBlack, kRed, kGreen, kBlue, kMagenta, kOrange, kCyan, kGray+1};
	for(int iad = 0; iad < 8; iad++){
		h_efficiency_zVSr2_rate[iad]=new TH1D(Form("h_effiency_zVSr2_rate_eh%dad%d",EH[iad],AD[iad]),Form("h_effiency_zVSr2_rate_eh%dad%d",EH[iad],AD[iad]),NzPoints*Nr2Points,0.42+iad*.02,0.42+iad*.02+Nr2Points*NzPoints);
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){ //filling histograms
			for(int iz = 0; iz < NzPoints; iz++){
				//if(ir2+1 == Nr2Points) continue;
				h_efficiency_zVSr2_rate[iad]->SetBinContent(h_efficiency_zVSr2_rate[iad]->FindBin(iz*Nr2Points + ir2+1), efficiency_rate_zVSr2[iad][ir2][iz]);
				h_efficiency_zVSr2_rate[iad]->SetBinError(h_efficiency_zVSr2_rate[iad]->FindBin(iz*Nr2Points + ir2+1), efficiencyError_rate_zVSr2[iad][ir2][iz]);
				//if(iad == 0) cout << efficiencyError_rate_zVSr2[iad][ir2] << endl;
			}
		}
		h_efficiency_zVSr2_rate[iad]->SetMarkerStyle(21);
		h_efficiency_zVSr2_rate[iad]->SetLineWidth(1);
		h_efficiency_zVSr2_rate[iad]->SetMarkerSize(1.5);
		h_efficiency_zVSr2_rate[iad]->SetLineColor(colors_AD[iad]);
		h_efficiency_zVSr2_rate[iad]->SetMarkerColor(colors_AD[iad]);
		h_efficiency_zVSr2_rate[iad]->GetXaxis()->SetTitle("voxel ID");
		h_efficiency_zVSr2_rate[iad]->GetYaxis()->SetTitle("Efficiency");
		//h_efficiency_zVSr2_rate[iad]->GetYaxis()->SetRangeUser(0.87,0.985);
		h_efficiency_zVSr2_rate[iad]->SetStats(0);
	}

	TCanvas *zVSr2_slice_effs_rate = new TCanvas("zVSr2_slice_effs_rate","zVSr2_slice_effs_rate");
		zVSr2_slice_effs_rate->cd();
		h_efficiency_zVSr2_rate[0]->Draw("e1x0");
		h_efficiency_zVSr2_rate[0]->GetXaxis()->SetTitle("voxel ID");
		h_efficiency_zVSr2_rate[0]->GetYaxis()->SetTitle("Efficiency");
		for(int i = 1; i < 8; i++){
			h_efficiency_zVSr2_rate[i]->Draw("e1x0 same");
		}
	zVSr2_slice_effs_rate->SetGridy();
	zVSr2_slice_effs_rate->BuildLegend();

*/


/*	TCanvas *z_slice_effs_nearVSfar = new TCanvas("z_slice_effs_nearVSfar","z_slice_effs_nearVSfar");
	z_slice_effs_nearVSfar->Divide(2,1);
	z_slice_effs_nearVSfar->cd(1);
		h_efficiency_z_rate[0]->GetYaxis()->SetRangeUser(0.75,1);
		h_efficiency_z_rate[0]->Draw("e1x0");
		for(int i = 1; i< 4; i++){
			h_efficiency_z_rate[i]->Draw("e1x0 same");
		}
	z_slice_effs_nearVSfar->cd(2);
		h_efficiency_z_rate[4]->GetYaxis()->SetRangeUser(0.75,1);
		h_efficiency_z_rate[4]->Draw("e1x0");
		for(int i = 4; i< 8; i++){
			h_efficiency_z_rate[i]->Draw("e1x0 same");
		}
	z_slice_effs_nearVSfar->BuildLegend();*/

	TCanvas *eff_sigma_rate = new TCanvas("eff_sigma_rate","eff_sigma_rate");
	eff_sigma_rate->cd();
		for(int iad = 0; iad < 8; iad++){
			h_efficiency_sigma_rate[iad]->SetStats(0);
			h_efficiency_sigma_rate[iad]->SetLineColor(colors_AD[iad]);
			h_efficiency_sigma_rate[iad]->SetMarkerColor(colors_AD[iad]);
			h_efficiency_sigma_rate[iad]->GetXaxis()->SetTitle("Sigma");
			h_efficiency_sigma_rate[iad]->GetYaxis()->SetTitle("Delayed Cut Efficiency");
			if(iad == 0) h_efficiency_sigma_rate[iad]->Draw();
			else h_efficiency_sigma_rate[iad]->Draw("same");
		}
	eff_sigma_rate->SetGridy();
	eff_sigma_rate->BuildLegend();
	eff_sigma_rate->Print("../nH_files/effSigma_rate.png");

	TCanvas *eff_sigma_norm = new TCanvas("eff_sigma_norm","eff_sigma_norm");
	eff_sigma_norm->cd();
		for(int iad = 0; iad < 8; iad++){
			h_efficiency_sigma_norm[iad]->SetStats(0);
			h_efficiency_sigma_norm[iad]->SetLineColor(colors_AD[iad]);
			h_efficiency_sigma_norm[iad]->SetMarkerColor(colors_AD[iad]);
			h_efficiency_sigma_norm[iad]->GetXaxis()->SetTitle("Sigma");
			h_efficiency_sigma_norm[iad]->GetYaxis()->SetTitle("Delayed Cut Efficiency");
			if(iad == 0) h_efficiency_sigma_norm[iad]->Draw();
			else h_efficiency_sigma_norm[iad]->Draw("same");
		}
	eff_sigma_norm->SetGridy();
	eff_sigma_norm->BuildLegend();

	TCanvas *eff_sigma_DTnorm = new TCanvas("eff_sigma_DTnorm","eff_sigma_DTnorm");
	eff_sigma_DTnorm->cd();
		for(int iad = 0; iad < 8; iad++){
			h_efficiency_sigma_DTnorm[iad]->SetStats(0);
			h_efficiency_sigma_DTnorm[iad]->SetLineColor(colors_AD[iad]);
			h_efficiency_sigma_DTnorm[iad]->SetMarkerColor(colors_AD[iad]);
			h_efficiency_sigma_DTnorm[iad]->GetXaxis()->SetTitle("Sigma");
			h_efficiency_sigma_DTnorm[iad]->GetYaxis()->SetTitle("Delayed Cut Efficiency");
			if(iad == 0) h_efficiency_sigma_DTnorm[iad]->Draw();
			else h_efficiency_sigma_DTnorm[iad]->Draw("same");
		}
	eff_sigma_DTnorm->SetGridy();
	eff_sigma_DTnorm->BuildLegend();



	TH1D* h_eff_sigma_near_rate = new TH1D("h_eff_sigma_near_rate","h_eff_sigma_near_rate",500,0,5);
	TH1D* h_eff_sigma_near_norm = new TH1D("h_eff_sigma_near_norm","h_eff_sigma_near_norm",500,0,5);
	TH1D* h_eff_sigma_near_DTnorm = new TH1D("h_eff_sigma_near_DTnorm","h_eff_sigma_near_DTnorm",500,0,5);
	TH1D* h_eff_sigma_far_rate = new TH1D("h_eff_sigma_far_rate","h_eff_sigma_far_rate",500,0,5);
	TH1D* h_eff_sigma_far_norm = new TH1D("h_eff_sigma_far_norm","h_eff_sigma_far_norm",500,0,5);
	TH1D* h_eff_sigma_far_DTnorm = new TH1D("h_eff_sigma_far_DTnorm","h_eff_sigma_far_DTnorm",500,0,5);

	TH1D* h_eff_sigma_nearFarDiff_rate = new TH1D("h_eff_sigma_nearFarDiff_rate","h_eff_sigma_nearFarDiff_rate",500,0,5);
	TH1D* h_eff_sigma_nearFarDiff_norm = new TH1D("h_eff_sigma_nearFarDiff_norm","h_eff_sigma_nearFarDiff_norm",500,0,5);
	TH1D* h_eff_sigma_nearFarDiff_DTnorm = new TH1D("h_eff_sigma_nearFarDiff_DTnorm","h_eff_sigma_nearFarDiff_DTnorm",500,0,5);

	for(int iBin = 1; iBin < 501; iBin++){
		h_eff_sigma_near_rate->SetBinContent(iBin, .25*(h_efficiency_sigma_rate[0]->GetBinContent(iBin)+h_efficiency_sigma_rate[1]->GetBinContent(iBin)+h_efficiency_sigma_rate[2]->GetBinContent(iBin)+h_efficiency_sigma_rate[3]->GetBinContent(iBin)));
		h_eff_sigma_far_rate->SetBinContent(iBin, .25*(h_efficiency_sigma_rate[4]->GetBinContent(iBin)+h_efficiency_sigma_rate[5]->GetBinContent(iBin)+h_efficiency_sigma_rate[6]->GetBinContent(iBin)+h_efficiency_sigma_rate[7]->GetBinContent(iBin)));
		h_eff_sigma_nearFarDiff_rate->SetBinContent(iBin, h_eff_sigma_near_rate->GetBinContent(iBin) - h_eff_sigma_far_rate->GetBinContent(iBin));

		h_eff_sigma_near_norm->SetBinContent(iBin, .25*(h_efficiency_sigma_norm[0]->GetBinContent(iBin)+h_efficiency_sigma_norm[1]->GetBinContent(iBin)+h_efficiency_sigma_norm[2]->GetBinContent(iBin)+h_efficiency_sigma_norm[3]->GetBinContent(iBin)));
		h_eff_sigma_far_norm->SetBinContent(iBin, .25*(h_efficiency_sigma_norm[4]->GetBinContent(iBin)+h_efficiency_sigma_norm[5]->GetBinContent(iBin)+h_efficiency_sigma_norm[6]->GetBinContent(iBin)+h_efficiency_sigma_norm[7]->GetBinContent(iBin)));
		h_eff_sigma_nearFarDiff_norm->SetBinContent(iBin, h_eff_sigma_near_norm->GetBinContent(iBin) - h_eff_sigma_far_norm->GetBinContent(iBin));

		h_eff_sigma_near_DTnorm->SetBinContent(iBin, .25*(h_efficiency_sigma_DTnorm[0]->GetBinContent(iBin)+h_efficiency_sigma_DTnorm[1]->GetBinContent(iBin)+h_efficiency_sigma_DTnorm[2]->GetBinContent(iBin)+h_efficiency_sigma_DTnorm[3]->GetBinContent(iBin)));
		h_eff_sigma_far_DTnorm->SetBinContent(iBin, .25*(h_efficiency_sigma_DTnorm[4]->GetBinContent(iBin)+h_efficiency_sigma_DTnorm[5]->GetBinContent(iBin)+h_efficiency_sigma_DTnorm[6]->GetBinContent(iBin)+h_efficiency_sigma_DTnorm[7]->GetBinContent(iBin)));
		h_eff_sigma_nearFarDiff_DTnorm->SetBinContent(iBin, h_eff_sigma_near_DTnorm->GetBinContent(iBin) - h_eff_sigma_far_DTnorm->GetBinContent(iBin));
	}

	TCanvas *eff_sigma_nearFarDiff_rate = new TCanvas("eff_sigma_nearFarDiff_rate","eff_sigma_nearFarDiff_rate");
	eff_sigma_nearFarDiff_rate->cd();
		h_eff_sigma_nearFarDiff_rate->SetStats(0);
		h_eff_sigma_nearFarDiff_rate->GetXaxis()->SetTitle("Sigma");
		h_eff_sigma_nearFarDiff_rate->GetYaxis()->SetTitle("Near Avg - Far Avg Efficiency");
		h_eff_sigma_nearFarDiff_rate->Draw();
	eff_sigma_nearFarDiff_rate->SetGridy();
	eff_sigma_nearFarDiff_rate->Print("../nH_files/effSigma_nearFarDiff_rate.png");

	TCanvas *eff_sigma_nearFarDiff_norm = new TCanvas("eff_sigma_nearFarDiff_norm","eff_sigma_nearFarDiff_norm");
	eff_sigma_nearFarDiff_norm->cd();
		h_eff_sigma_nearFarDiff_norm->SetStats(0);
		h_eff_sigma_nearFarDiff_norm->GetXaxis()->SetTitle("Sigma");
		h_eff_sigma_nearFarDiff_norm->GetYaxis()->SetTitle("Near Avg - Far Avg Efficiency");
		h_eff_sigma_nearFarDiff_norm->Draw();
	eff_sigma_nearFarDiff_norm->SetGridy();

	TCanvas *eff_sigma_nearFarDiff_DTnorm = new TCanvas("eff_sigma_nearFarDiff_DTnorm","eff_sigma_nearFarDiff_DTnorm");
	eff_sigma_nearFarDiff_DTnorm->cd();
		h_eff_sigma_nearFarDiff_DTnorm->SetStats(0);
		h_eff_sigma_nearFarDiff_DTnorm->GetXaxis()->SetTitle("Sigma");
		h_eff_sigma_nearFarDiff_DTnorm->GetYaxis()->SetTitle("Near Avg - Far Avg Efficiency");
		h_eff_sigma_nearFarDiff_DTnorm->Draw();
	eff_sigma_nearFarDiff_DTnorm->SetGridy();



//Fit the z slices separately
/*	TF1* z_Ed_Fit_rate = new TF1("z_Ed_Fit_rate", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.5,2.8);
	double mu_rate_z[8][NzPoints];
	double sigma_rate_z[8][NzPoints];
	double alpha_rate_z[8][NzPoints];
	double lambda_rate_z[8][NzPoints];

	double muErr_rate_z[8][NzPoints];
	double sigmaErr_rate_z[8][NzPoints];
	double alphaErr_rate_z[8][NzPoints];
	double lambdaErr_rate_z[8][NzPoints];
	for(int iad = 0; iad < 8; iad++){
		for(int iz = 0; iz < NzPoints; iz++){
			mu_rate_z[iad][iz] = 0;
			sigma_rate_z[iad][iz] = 0;
			alpha_rate_z[iad][iz] = 0;
			lambda_rate_z[iad][iz] = 0;

			muErr_rate_z[iad][iz] = 0;
			sigmaErr_rate_z[iad][iz] = 0;
			alphaErr_rate_z[iad][iz] = 0;
			lambdaErr_rate_z[iad][iz] = 0;
		}
	}

	TH1D* h_mu_rate_z[8];
	TH1D* h_sigma_rate_z[8];
	TH1D* h_alpha_rate_z[8];
	TH1D* h_lambda_rate_z[8];
	for(int iad = 0; iad < 8; iad++){
		h_mu_rate_z[iad] = new TH1D(Form("h_mu_rate_z_eh%dad%d",EH[iad],AD[iad]),Form("h_mu_rate_z_eh%dad%d",EH[iad],AD[iad]),NzPoints,0.42+iad*.02,0.42+iad*.02+NzPoints);
		h_sigma_rate_z[iad] = new TH1D(Form("h_sigma_rate_z_eh%dad%d",EH[iad],AD[iad]),Form("h_sigma_rate_z_eh%dad%d",EH[iad],AD[iad]),NzPoints,0.42+iad*.02,0.42+iad*.02+NzPoints);
		h_alpha_rate_z[iad] = new TH1D(Form("h_alpha_rate_z_eh%dad%d",EH[iad],AD[iad]),Form("h_alpha_rate_z_eh%dad%d",EH[iad],AD[iad]),NzPoints,0.42+iad*.02,0.42+iad*.02+NzPoints);
		h_lambda_rate_z[iad] = new TH1D(Form("h_lambda_rate_z_eh%dad%d",EH[iad],AD[iad]),Form("h_lambda_rate_z_eh%dad%d",EH[iad],AD[iad]),NzPoints,0.42+iad*.02,0.42+iad*.02+NzPoints);
	}
	for(int iad = 0; iad < 8; iad++){
		for(int iz = 0; iz < NzPoints; iz++){
			cout << "AD " << iad+1 << "\tz slice: " << iz+1 << endl;

				z_Ed_Fit_rate->SetParameter(0,(h_Edelayed_sub_rate_z[iad][iz]->GetBinContent(h_Edelayed_sub_rate_z[iad][iz]->FindBin(2.3)))/3.);
				z_Ed_Fit_rate->SetParameter(1,0.8); //alpha
				z_Ed_Fit_rate->SetParameter(2,2.3); //mu
					if(iad == 5) z_Ed_Fit_rate->SetParameter(2,2.2); //sigma
				z_Ed_Fit_rate->SetParameter(3,0.135); //sigma
					if(iad == 6) z_Ed_Fit_rate->SetParameter(3,0.136); //sigma
					if(iad == 5) z_Ed_Fit_rate->SetParameter(3,0.13); //sigma
				z_Ed_Fit_rate->SetParameter(4,3); //lambda

				z_Ed_Fit_rate->SetParName(0,"N"); //normalization
				z_Ed_Fit_rate->SetParName(1,"alpha"); //alpha
				z_Ed_Fit_rate->SetParName(2,"mu"); //mu
				z_Ed_Fit_rate->SetParName(3,"sigma"); //sigma
				z_Ed_Fit_rate->SetParName(4,"lambda"); //lambda
			h_Edelayed_sub_rate_z[iad][iz]->Fit("z_Ed_Fit_rate", "R");
			h_Edelayed_sub_rate_z[iad][iz]->Write();

			mu_rate_z[iad][iz] = z_Ed_Fit_rate->GetParameter("mu");
			sigma_rate_z[iad][iz] = z_Ed_Fit_rate->GetParameter("sigma");
			alpha_rate_z[iad][iz] = z_Ed_Fit_rate->GetParameter("alpha");
			lambda_rate_z[iad][iz] = z_Ed_Fit_rate->GetParameter("lambda");

			muErr_rate_z[iad][iz] = z_Ed_Fit_rate->GetParError(2);
			sigmaErr_rate_z[iad][iz] = z_Ed_Fit_rate->GetParError(3);
			alphaErr_rate_z[iad][iz] = z_Ed_Fit_rate->GetParError(1);
			lambdaErr_rate_z[iad][iz] = z_Ed_Fit_rate->GetParError(4);

			h_mu_rate_z[iad]->SetBinContent(iz+1, mu_rate_z[iad][iz]);
			h_sigma_rate_z[iad]->SetBinContent(iz+1, sigma_rate_z[iad][iz]);
			h_alpha_rate_z[iad]->SetBinContent(iz+1, alpha_rate_z[iad][iz]);
			h_lambda_rate_z[iad]->SetBinContent(iz+1, lambda_rate_z[iad][iz]);

			h_mu_rate_z[iad]->SetBinError(iz+1, muErr_rate_z[iad][iz]);
			h_sigma_rate_z[iad]->SetBinError(iz+1, sigmaErr_rate_z[iad][iz]);
			h_alpha_rate_z[iad]->SetBinError(iz+1, alphaErr_rate_z[iad][iz]);
			h_lambda_rate_z[iad]->SetBinError(iz+1, lambdaErr_rate_z[iad][iz]);
		}

		h_mu_rate_z[iad]->SetStats(0);
		h_mu_rate_z[iad]->GetXaxis()->SetTitle("z slice");
		h_mu_rate_z[iad]->GetYaxis()->SetTitle("mu fit value");
		h_mu_rate_z[iad]->GetYaxis()->SetRangeUser(2.,2.5);
		h_mu_rate_z[iad]->SetMarkerStyle(21);
		h_mu_rate_z[iad]->SetMarkerColor(colors_AD[iad]);
		h_mu_rate_z[iad]->SetLineColor(colors_AD[iad]);

		h_sigma_rate_z[iad]->SetStats(0);
		h_sigma_rate_z[iad]->GetXaxis()->SetTitle("z slice");
		h_sigma_rate_z[iad]->GetYaxis()->SetTitle("sigma fit value");
		h_sigma_rate_z[iad]->GetYaxis()->SetRangeUser(0.1,0.15);
		h_sigma_rate_z[iad]->SetMarkerStyle(21);
		h_sigma_rate_z[iad]->SetMarkerColor(colors_AD[iad]);
		h_sigma_rate_z[iad]->SetLineColor(colors_AD[iad]);

		h_alpha_rate_z[iad]->SetStats(0);
		h_alpha_rate_z[iad]->GetXaxis()->SetTitle("z slice");
		h_alpha_rate_z[iad]->GetYaxis()->SetTitle("alpha fit value");
		//h_alpha_rate_z[iad]->GetYaxis()->SetRangeUser(2.,2.5);
		h_alpha_rate_z[iad]->SetMarkerStyle(21);
		h_alpha_rate_z[iad]->SetMarkerColor(colors_AD[iad]);
		h_alpha_rate_z[iad]->SetLineColor(colors_AD[iad]);

		h_lambda_rate_z[iad]->SetStats(0);
		h_lambda_rate_z[iad]->GetXaxis()->SetTitle("z slice");
		h_lambda_rate_z[iad]->GetYaxis()->SetTitle("lambda fit value");
		//h_lambda_rate_z[iad]->GetYaxis()->SetRangeUser(0.1,0.15);
		h_lambda_rate_z[iad]->SetMarkerStyle(21);
		h_lambda_rate_z[iad]->SetMarkerColor(colors_AD[iad]);
		h_lambda_rate_z[iad]->SetLineColor(colors_AD[iad]);
	}

	TCanvas *mu_rate_z_ADs = new TCanvas("mu_rate_z_ADs","mu_rate_z_ADs");
	mu_rate_z_ADs->cd();
		h_mu_rate_z[0]->Draw("e1x0");
		for(int iad = 1; iad < 8; iad++){
			h_mu_rate_z[iad]->Draw("e1x0 same");
		}
	mu_rate_z_ADs->BuildLegend();

	TCanvas *sigma_rate_z_ADs = new TCanvas("sigma_rate_z_ADs","sigma_rate_z_ADs");
	sigma_rate_z_ADs->cd();
		h_sigma_rate_z[0]->Draw("e1x0");
		for(int iad = 1; iad < 8; iad++){
			h_sigma_rate_z[iad]->Draw("e1x0 same");
		}
	sigma_rate_z_ADs->BuildLegend();

	TCanvas *alpha_rate_z_ADs = new TCanvas("alpha_rate_z_ADs","alpha_rate_z_ADs");
	alpha_rate_z_ADs->cd();
		h_alpha_rate_z[0]->Draw("e1x0");
		for(int iad = 1; iad < 8; iad++){
			h_alpha_rate_z[iad]->Draw("e1x0 same");
		}
	alpha_rate_z_ADs->BuildLegend();

	TCanvas *lambda_rate_z_ADs = new TCanvas("lambda_rate_z_ADs","lambda_rate_z_ADs");
	lambda_rate_z_ADs->cd();
		h_lambda_rate_z[0]->Draw("e1x0");
		for(int iad = 1; iad < 8; iad++){
			h_lambda_rate_z[iad]->Draw("e1x0 same");
		}
	lambda_rate_z_ADs->BuildLegend();

*/

	int colors_zSlices[8] = {kBlack, kRed, kGreen, kBlue, kMagenta, kOrange, kCyan, kGray+1};

//Do a near-far comparison of shape for the z-slices here
/*	TH1F* h_Edelayed_sub_rate_z_near[NzPoints];
	TH1F* h_Edelayed_sub_rate_z_far[NzPoints];
	TH1F* h_Edelayed_sub_rate_z_nfRatio[NzPoints];
	for(int iz = 0; iz < NzPoints; iz++){
		h_Edelayed_sub_rate_z_near[iz] = (TH1F*)h_Edelayed_sub_rate_z[0][iz]->Clone();
		for(int iad = 1; iad < 4; iad++){
			h_Edelayed_sub_rate_z_near[iz]->Add(h_Edelayed_sub_rate_z[iad][iz]);
		}
		h_Edelayed_sub_rate_z_far[iz] = (TH1F*)h_Edelayed_sub_rate_z[4][iz]->Clone();
		for(int iad = 5; iad < 8; iad++){
			h_Edelayed_sub_rate_z_far[iz]->Add(h_Edelayed_sub_rate_z[iad][iz]);
		}

		h_Edelayed_sub_rate_z_near[iz]->Scale(1./(h_Edelayed_sub_rate_z_near[iz]->Integral(h_Edelayed_sub_rate_z_near[iz]->FindBin(1.5),h_Edelayed_sub_rate_z_near[iz]->FindBin(2.8))));
		h_Edelayed_sub_rate_z_far[iz]->Scale(1./(h_Edelayed_sub_rate_z_far[iz]->Integral(h_Edelayed_sub_rate_z_far[iz]->FindBin(1.5),h_Edelayed_sub_rate_z_far[iz]->FindBin(2.8))));

		h_Edelayed_sub_rate_z_near[iz]->Rebin(2);
		h_Edelayed_sub_rate_z_near[iz]->SetStats(0);
		h_Edelayed_sub_rate_z_near[iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_Edelayed_sub_rate_z_near[iz]->GetXaxis()->SetRangeUser(1.5,3.);
		h_Edelayed_sub_rate_z_near[iz]->GetYaxis()->SetTitle("Counts");

		h_Edelayed_sub_rate_z_far[iz]->Rebin(2);
		h_Edelayed_sub_rate_z_far[iz]->SetStats(0);
		h_Edelayed_sub_rate_z_far[iz]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_Edelayed_sub_rate_z_far[iz]->GetXaxis()->SetRangeUser(1.5,3.);
		h_Edelayed_sub_rate_z_far[iz]->GetYaxis()->SetTitle("Counts");

		h_Edelayed_sub_rate_z_nfRatio[iz] = (TH1F*)h_Edelayed_sub_rate_z_near[iz]->Clone();
		h_Edelayed_sub_rate_z_nfRatio[iz]->Divide(h_Edelayed_sub_rate_z_far[iz]);


	TCanvas *nearFar_EdSub_z = new TCanvas(Form("nearFar_EdSub_z%d",iz+1),Form("nearFar_EdSub_z%d",iz+1));
	nearFar_EdSub_z->cd();
		h_Edelayed_sub_rate_z_near[iz]->Draw();
		h_Edelayed_sub_rate_z_far[iz]->Draw("same");
	//nearFar_EdSub_z->BuildLegend();
		h_Edelayed_sub_rate_z_near[iz]->SetLineColor(colors_zSlices[0]);
		h_Edelayed_sub_rate_z_far[iz]->SetLineColor(colors_zSlices[4]);

	TCanvas *nearFar_EdRatio_z = new TCanvas(Form("nearFar_EdRatio_z%d",iz+1),Form("nearFar_EdRatio_z%d",iz+1));
	nearFar_EdRatio_z->cd();
		h_Edelayed_sub_rate_z_nfRatio[iz]->Draw();

	}
*/
/*	TCanvas *nearADs_EdSub_z = new TCanvas("nearADs_EdSub_z","nearADs_EdSub_z");
	nearADs_EdSub_z->cd();
		h_Edelayed_sub_rate_z_near[0]->Draw();
		for(int iz = 1; iz < NzPoints; iz++){
			h_Edelayed_sub_rate_z_near[iz]->Draw("same");
		}
	//nearADs_EdSub_z->BuildLegend();*/


//Do a near-far comparison of shape for the r2-slices here
/*	TH1F* h_Edelayed_sub_rate_r2_near[Nr2Points];
	TH1F* h_Edelayed_sub_rate_r2_far[Nr2Points];
	TH1F* h_Edelayed_sub_rate_r2_nfRatio[Nr2Points];
	for(int ir2 = 0; ir2 < Nr2Points; ir2++){
		h_Edelayed_sub_rate_r2_near[ir2] = (TH1F*)h_Edelayed_sub_rate_r2[0][ir2]->Clone();
		for(int iad = 1; iad < 4; iad++){
			h_Edelayed_sub_rate_r2_near[ir2]->Add(h_Edelayed_sub_rate_r2[iad][ir2]);
		}
		h_Edelayed_sub_rate_r2_far[ir2] = (TH1F*)h_Edelayed_sub_rate_r2[4][ir2]->Clone();
		for(int iad = 5; iad < 8; iad++){
			h_Edelayed_sub_rate_r2_far[ir2]->Add(h_Edelayed_sub_rate_r2[iad][ir2]);
		}

		h_Edelayed_sub_rate_r2_near[ir2]->Scale(1./(h_Edelayed_sub_rate_r2_near[ir2]->Integral(h_Edelayed_sub_rate_r2_near[ir2]->FindBin(1.5),h_Edelayed_sub_rate_r2_near[ir2]->FindBin(2.8))));
		h_Edelayed_sub_rate_r2_far[ir2]->Scale(1./(h_Edelayed_sub_rate_r2_far[ir2]->Integral(h_Edelayed_sub_rate_r2_far[ir2]->FindBin(1.5),h_Edelayed_sub_rate_r2_far[ir2]->FindBin(2.8))));

		h_Edelayed_sub_rate_r2_near[ir2]->Rebin(2);
		h_Edelayed_sub_rate_r2_near[ir2]->SetStats(0);
		h_Edelayed_sub_rate_r2_near[ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_Edelayed_sub_rate_r2_near[ir2]->GetXaxis()->SetRangeUser(1.5,3.);
		h_Edelayed_sub_rate_r2_near[ir2]->GetYaxis()->SetTitle("Counts");

		h_Edelayed_sub_rate_r2_far[ir2]->Rebin(2);
		h_Edelayed_sub_rate_r2_far[ir2]->SetStats(0);
		h_Edelayed_sub_rate_r2_far[ir2]->GetXaxis()->SetTitle("Delayed Energy [MeV]");
		h_Edelayed_sub_rate_r2_far[ir2]->GetXaxis()->SetRangeUser(1.5,3.);
		h_Edelayed_sub_rate_r2_far[ir2]->GetYaxis()->SetTitle("Counts");

		h_Edelayed_sub_rate_r2_nfRatio[ir2] = (TH1F*)h_Edelayed_sub_rate_r2_near[ir2]->Clone();
		h_Edelayed_sub_rate_r2_nfRatio[ir2]->Divide(h_Edelayed_sub_rate_r2_far[ir2]);


	TCanvas *nearFar_EdSub_r2 = new TCanvas(Form("nearFar_EdSub_r2%d",ir2+1),Form("nearFar_EdSub_r2%d",ir2+1));
	nearFar_EdSub_r2->cd();
		h_Edelayed_sub_rate_r2_near[ir2]->Draw();
		h_Edelayed_sub_rate_r2_far[ir2]->Draw("same");
	//nearFar_EdSub_r2->BuildLegend();
		h_Edelayed_sub_rate_r2_near[ir2]->SetLineColor(colors_zSlices[0]);
		h_Edelayed_sub_rate_r2_far[ir2]->SetLineColor(colors_zSlices[4]);

	TCanvas *nearFar_EdRatio_r2 = new TCanvas(Form("nearFar_EdRatio_r2%d",ir2+1),Form("nearFar_EdRatio_r2%d",ir2+1));
	nearFar_EdRatio_r2->cd();
		h_Edelayed_sub_rate_r2_nfRatio[ir2]->Draw();
	}
*/
/*	TCanvas *nearADs_EdSub_r2 = new TCanvas("nearADs_EdSub_r2","nearADs_EdSub_r2");
	nearADs_EdSub_r2->cd();
		h_Edelayed_sub_rate_r2_near[0]->Draw();
		for(int ir2 = 1; ir2 < Nr2Points; ir2++){
			h_Edelayed_sub_rate_r2_near[ir2]->Draw("same");
		}
	//nearADs_EdSub_r2->BuildLegend();*/

		TH1F* Ed_rate_totAD_near = (TH1F*)outfile->Get(Form("h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[0], AD[0]));
		TH1F* Ed_rate_totAD_far = (TH1F*)outfile->Get(Form("h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[4], AD[4]));

		for(int iad = 1; iad < 4; iad++){
				Ed_rate_totAD_near->Add((TH1F*)outfile->Get(Form("h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[0+iad], AD[0+iad])));
				Ed_rate_totAD_far->Add((TH1F*)outfile->Get(Form("h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[4+iad], AD[4+iad])));
		}

		double peak_near = (peak_rate[0] + peak_rate[1] + peak_rate[2] + peak_rate[3])/4.;
		double peak_far = (peak_rate[4] + peak_rate[5] + peak_rate[6] + peak_rate[7])/4.;

		double sigma_near = (sigma_rate[0] + sigma_rate[1] + sigma_rate[2] + sigma_rate[3])/4.;
		double sigma_far = (sigma_rate[4] + sigma_rate[5] + sigma_rate[6] + sigma_rate[7])/4.;

		Ed_rate_totAD_near->Scale(1./(Ed_rate_totAD_near->Integral(Ed_rate_totAD_near->FindBin(peak_near - 3*sigma_near),Ed_rate_totAD_near->FindBin(peak_near + 3*sigma_far))));
		Ed_rate_totAD_far->Scale(1./(Ed_rate_totAD_far->Integral(Ed_rate_totAD_far->FindBin(peak_far - 3*sigma_far),Ed_rate_totAD_far->FindBin(peak_far + 3*sigma_far))));

	TCanvas *nearFar_EdSub_wholeAD = new TCanvas(Form("nearFar_EdSub_wholeAD"),Form("nearFar_EdSub_wholeAD"));
	nearFar_EdSub_wholeAD->cd();
		Ed_rate_totAD_near->Rebin(20);
		Ed_rate_totAD_far->Rebin(20);
		Ed_rate_totAD_near->GetXaxis()->SetRangeUser(1.5,3);
		Ed_rate_totAD_near->GetYaxis()->SetRangeUser(0,0.06);
		Ed_rate_totAD_near->SetStats(0);
		Ed_rate_totAD_near->Draw();
		Ed_rate_totAD_far->Draw("same");
		Ed_rate_totAD_near->SetLineColor(colors_zSlices[0]);
		Ed_rate_totAD_far->SetLineColor(colors_zSlices[4]);
		Ed_rate_totAD_near->SetLineWidth(2);
		Ed_rate_totAD_far->SetLineWidth(2);






//zVSr2 efficiency comparisons
/*	TH1D* h_efficiency_zVSr2_rate_nearVSfar[2];
	h_efficiency_zVSr2_rate_nearVSfar[0]=new TH1D("h_effiency_zVSr2_rate_near","h_effiency_zVSr2_rate_near",NzPoints*Nr2Points,0.42+0*.02,0.42+0*.02+Nr2Points*NzPoints);
	h_efficiency_zVSr2_rate_nearVSfar[1]=new TH1D("h_effiency_zVSr2_rate_far","h_effiency_zVSr2_rate_far",NzPoints*Nr2Points,0.42+(1+4)*.02,0.42+(1+4)*.02+Nr2Points*NzPoints);
	for(int iad = 0; iad < 2; iad++){
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){ //filling histograms
			for(int iz = 0; iz < NzPoints; iz++){
			//	if(ir2+1 == Nr2Points) continue;
				if(iad == 0){
					h_efficiency_zVSr2_rate_nearVSfar[iad]->SetBinContent(h_efficiency_zVSr2_rate_nearVSfar[iad]->FindBin(iz*Nr2Points + ir2+1), efficiency_rate_zVSr2_near[ir2][iz]);
				//	h_efficiency_zVSr2_rate_nearVSfar[iad]->SetBinError(h_efficiency_zVSr2_rate_nearVSfar[iad]->FindBin(iz*Nr2Points + ir2+1), efficiencyError_rate_zVSr2_near[ir2][iz]);
				}
				else{
					h_efficiency_zVSr2_rate_nearVSfar[iad]->SetBinContent(h_efficiency_zVSr2_rate_nearVSfar[iad]->FindBin(iz*Nr2Points + ir2+1), efficiency_rate_zVSr2_far[ir2][iz]);
				//	h_efficiency_zVSr2_rate_nearVSfar[iad]->SetBinError(h_efficiency_zVSr2_rate_nearVSfar[iad]->FindBin(iz*Nr2Points + ir2+1), efficiencyError_rate_zVSr2_far[ir2][iz]);
				}
			}
		}
		h_efficiency_zVSr2_rate_nearVSfar[iad]->SetMarkerStyle(21);
		h_efficiency_zVSr2_rate_nearVSfar[iad]->SetLineWidth(1);
		h_efficiency_zVSr2_rate_nearVSfar[iad]->SetMarkerSize(1.5);
		h_efficiency_zVSr2_rate_nearVSfar[iad]->SetLineColor(colors_AD[iad]);
		h_efficiency_zVSr2_rate_nearVSfar[iad]->SetMarkerColor(colors_AD[iad]);
		h_efficiency_zVSr2_rate_nearVSfar[iad]->GetXaxis()->SetTitle("voxel ID");
		h_efficiency_zVSr2_rate_nearVSfar[iad]->GetYaxis()->SetTitle("Efficiency");
		//h_efficiency_zVSr2_rate_nearVSfar[iad]->GetYaxis()->SetRangeUser(0.87,0.985);
		h_efficiency_zVSr2_rate_nearVSfar[iad]->SetStats(0);
	}

	TCanvas *zVSr2_slice_effs_rate_nearVSfar = new TCanvas("zVSr2_slice_effs_rate_nearVSfar","zVSr2_slice_effs_rate_nearVSfar");
		zVSr2_slice_effs_rate_nearVSfar->cd();
		h_efficiency_zVSr2_rate_nearVSfar[0]->Draw("e1x0");
		h_efficiency_zVSr2_rate_nearVSfar[0]->GetXaxis()->SetTitle("voxel ID");
		h_efficiency_zVSr2_rate_nearVSfar[0]->GetYaxis()->SetTitle("Efficiency");
		h_efficiency_zVSr2_rate_nearVSfar[1]->Draw("e1x0 same");
	zVSr2_slice_effs_rate_nearVSfar->SetGridy();
	zVSr2_slice_effs_rate_nearVSfar->BuildLegend();
*/


}


void delayed_NvsF(int far_shift, int extRange){ //for extRange: 0 is fixed 1.5-2.8 MeV, 1 is just 1.5 - peak+3sig MeV, 2 is peak-3sig - 2.8 MeV, 3 is +-5sig, 3 is +-4sig
	char name[64];
	int EH[8] = {1,1,2,2,3,3,3,3};
	int AD[8] = {1,2,1,2,1,2,3,4};

	const int NzPoints = 5;
	const int Nr2Points = 5;

//Opening input file:
	char inputname[64];
	char title[64];

cout << "1" << endl;
        sprintf(inputname,"../nH_files/delayedEnergy.root");
	TFile* infile=new TFile(inputname, "READ");
cout << "2" << endl;
	TH1F* Ed_rate_totAD_near = (TH1F*)infile->Get(Form("h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[0], AD[0]));
	TH1F* Ed_rate_totAD_far = (TH1F*)infile->Get(Form("h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[4], AD[4]));
	TH1F* Ed_ibd_totAD_near = (TH1F*)infile->Get(Form("h_Edelayed_ibd_eh%dad%d", EH[0], AD[0]));
	TH1F* Ed_ibd_totAD_far = (TH1F*)infile->Get(Form("h_Edelayed_ibd_eh%dad%d", EH[4], AD[4]));
cout << "3" << endl;
	for(int iad = 1; iad < 4; iad++){
		Ed_rate_totAD_near->Add((TH1F*)infile->Get(Form("h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[0+iad], AD[0+iad])));
		Ed_rate_totAD_far->Add((TH1F*)infile->Get(Form("h_Edelayed_subtract_DT800_rate_eh%dad%d", EH[4+iad], AD[4+iad])));

		Ed_ibd_totAD_near->Add((TH1F*)infile->Get(Form("h_Edelayed_ibd_eh%dad%d", EH[0+iad], AD[0+iad])));
		Ed_ibd_totAD_far->Add((TH1F*)infile->Get(Form("h_Edelayed_ibd_eh%dad%d", EH[4+iad], AD[4+iad])));
	}
cout << "4" << endl;
	TH1F* Ed_sub_near = (TH1F*)infile->Get(Form("h_Edelayed_subtract_fine_DT800_eh%dad%d", EH[0], AD[0]));
	TH1F* Ed_sub_far = (TH1F*)infile->Get(Form("h_Edelayed_subtract_fine_DT800_eh%dad%d", EH[4], AD[4]));
	TH1F* Ed_ibd_near = (TH1F*)infile->Get(Form("h_Edelayed_IBD_fine_DT800_eh%dad%d", EH[0], AD[0]));
	TH1F* Ed_ibd_far = (TH1F*)infile->Get(Form("h_Edelayed_IBD_fine_DT800_eh%dad%d", EH[4], AD[4]));
cout << "5" << endl;
	for(int iad = 1; iad < 4; iad++){
		Ed_sub_near->Add((TH1F*)infile->Get(Form("h_Edelayed_subtract_fine_DT800_eh%dad%d", EH[0+iad], AD[0+iad])));
		Ed_sub_far->Add((TH1F*)infile->Get(Form("h_Edelayed_subtract_fine_DT800_eh%dad%d", EH[4+iad], AD[4+iad])));

		Ed_ibd_near->Add((TH1F*)infile->Get(Form("h_Edelayed_IBD_fine_DT800_eh%dad%d", EH[0+iad], AD[0+iad])));
		Ed_ibd_far->Add((TH1F*)infile->Get(Form("h_Edelayed_IBD_fine_DT800_eh%dad%d", EH[4+iad], AD[4+iad])));
	}

cout << "6" << endl;
	TH1F* Ed_rate_zVSr2[8][Nr2Points][NzPoints];
	TH1F* Ed_ibd_zVSr2[8][Nr2Points][NzPoints];
	TH1F* Ed_rate_zVSr2_near[Nr2Points][NzPoints];
	TH1F* Ed_rate_zVSr2_far[Nr2Points][NzPoints];
	TH1F* Ed_ibd_zVSr2_near[Nr2Points][NzPoints];
	TH1F* Ed_ibd_zVSr2_far[Nr2Points][NzPoints];
	for(int ir2 = 0; ir2 < Nr2Points; ir2++){
		for(int iz = 0; iz < NzPoints; iz++){
			for(int iad = 0; iad < 8; iad++){
				Ed_rate_zVSr2[iad][ir2][iz] = (TH1F*)infile->Get(Form("h_Edelayed_sub_rate_zVSr2_points_eh%dad%d_ir2%d_iz%d", EH[iad], AD[iad], ir2, iz));
				Ed_ibd_zVSr2[iad][ir2][iz] = (TH1F*)infile->Get(Form("h_Edelayed_ibd_zVSr2_eh%dad%d_ir2%d_iz%d", EH[iad], AD[iad], ir2, iz));
			}
			for(int iad = 0; iad < 4; iad++){
				if(iad == 0){
					Ed_rate_zVSr2_near[ir2][iz] = (TH1F*)Ed_rate_zVSr2[iad][ir2][iz]->Clone();
					Ed_rate_zVSr2_far[ir2][iz] = (TH1F*)Ed_rate_zVSr2[iad+4][ir2][iz]->Clone();

					Ed_rate_zVSr2_near[ir2][iz]->SetName(Form("h_Ed_rate_zVSr2_near_ir2%d_iz%d", ir2, iz));
					Ed_rate_zVSr2_near[ir2][iz]->SetName(Form("h_Ed_rate_zVSr2_far_ir2%d_iz%d", ir2, iz));

					Ed_ibd_zVSr2_near[ir2][iz] = (TH1F*)Ed_ibd_zVSr2[iad][ir2][iz]->Clone();
					Ed_ibd_zVSr2_far[ir2][iz] = (TH1F*)Ed_ibd_zVSr2[iad+4][ir2][iz]->Clone();

					Ed_ibd_zVSr2_near[ir2][iz]->SetName(Form("h_Ed_ibd_zVSr2_near_ir2%d_iz%d", ir2, iz));
					Ed_ibd_zVSr2_near[ir2][iz]->SetName(Form("h_Ed_ibd_zVSr2_far_ir2%d_iz%d", ir2, iz));
				}
				else{
					Ed_rate_zVSr2_near[ir2][iz]->Add(Ed_rate_zVSr2[iad][ir2][iz]);
					Ed_rate_zVSr2_far[ir2][iz]->Add(Ed_rate_zVSr2[iad+4][ir2][iz]);

					Ed_ibd_zVSr2_near[ir2][iz]->Add(Ed_ibd_zVSr2[iad][ir2][iz]);
					Ed_ibd_zVSr2_far[ir2][iz]->Add(Ed_ibd_zVSr2[iad+4][ir2][iz]);
				}
			}
		}
	} // end of zVSr2 voxels part 1
cout << "7" << endl;

//Fitting the summed histograms

	TF1* delayedFit_rate_near = new TF1("delayedFit_rate_near", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);
	TF1* delayedFit_rate_far = new TF1("delayedFit_rate_far", "[0]*([1]*exp(-pow(x-[2],2)/(2*[3]*[3])) / ([3]*sqrt(2*TMath::Pi()))+(1.-[1])*[4]/(2*(exp([4]*[2])-1)) * exp([3]*[3]*[4]*[4]/2) * exp([4]*x) * ( TMath::Erf(([2]-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) - TMath::Erf((0-(x+[3]*[3]*[4]))/(sqrt(2)*[3])) ))",1.6,2.8);


		delayedFit_rate_near->SetParameter(0,(Ed_rate_totAD_near->GetBinContent(Ed_rate_totAD_near->FindBin(2.3)))/3.); //rate
		delayedFit_rate_near->SetParameter(1,0.8); //alpha
		delayedFit_rate_near->SetParameter(2,2.3); //mu
		delayedFit_rate_near->SetParameter(3,0.135); //sigma
		delayedFit_rate_near->SetParameter(4,3); //lambda
			delayedFit_rate_near->SetParName(0,"N"); //normalization
			delayedFit_rate_near->SetParName(1,"alpha"); //alpha
			delayedFit_rate_near->SetParName(2,"mu"); //mu
			delayedFit_rate_near->SetParName(3,"sigma"); //sigma
			delayedFit_rate_near->SetParName(4,"lambda"); //lambda
	Ed_rate_totAD_near->Fit("delayedFit_rate_near", "R");

		delayedFit_rate_far->SetParameter(0,(Ed_rate_totAD_far->GetBinContent(Ed_rate_totAD_far->FindBin(2.3)))/3.); //rate
		delayedFit_rate_far->SetParameter(1,0.8); //alpha
		delayedFit_rate_far->SetParameter(2,2.3); //mu
		delayedFit_rate_far->SetParameter(3,0.135); //sigma
		delayedFit_rate_far->SetParameter(4,3); //lambda
			delayedFit_rate_far->SetParName(0,"N"); //normalization
			delayedFit_rate_far->SetParName(1,"alpha"); //alpha
			delayedFit_rate_far->SetParName(2,"mu"); //mu
			delayedFit_rate_far->SetParName(3,"sigma"); //sigma
			delayedFit_rate_far->SetParName(4,"lambda"); //lambda
	Ed_rate_totAD_far->Fit("delayedFit_rate_far", "R");

//Getting the fit values:
	double peak_near = 0;
	double peak_far = 0;
	double sigma_near = 0;
	double sigma_far = 0;
	double peakErr_near = 0;
	double peakErr_far = 0;
	double sigmaErr_near = 0;
	double sigmaErr_far = 0;

	peak_near = delayedFit_rate_near->GetParameter("mu");
	peak_far = delayedFit_rate_far->GetParameter("mu");
	peakErr_near = delayedFit_rate_near->GetParError(2);
	peakErr_far = delayedFit_rate_far->GetParError(2);

	sigma_near = delayedFit_rate_near->GetParameter("sigma");
	sigma_far = delayedFit_rate_far->GetParameter("sigma");
	sigmaErr_near = delayedFit_rate_near->GetParError(3);
	sigmaErr_far = delayedFit_rate_far->GetParError(3);


//Calculate the efficiencies and error
	double Elow = 1.5;
	double Ehigh = 2.8;
	double Elow_far = Elow;
	double Ehigh_far = Ehigh;
	if(extRange == 1){
		Ehigh = peak_near + 3*sigma_near;
		Ehigh_far = peak_far + 3*sigma_far;
	}
	if(extRange == 2){
		Elow = peak_near - 3*sigma_near;
		Elow_far = peak_far - 3*sigma_far;
	}
	if(far_shift == 1 && extRange != 2) Elow_far = Elow_far+peak_far-peak_near;
	if(far_shift == 1 && extRange != 1) Ehigh_far = Ehigh_far+peak_far-peak_near;

	if(extRange == 3){
		Elow = peak_near - 5*sigma_near;
		Elow_far = peak_far - 5*sigma_far;
		Ehigh = peak_near + 5*sigma_near;
		Ehigh_far = peak_far + 5*sigma_far;
	}

	if(extRange == 4){
		Elow = peak_near - 4*sigma_near;
		Elow_far = peak_far - 4*sigma_far;
		Ehigh = peak_near + 4*sigma_near;
		Ehigh_far = peak_far + 4*sigma_far;
	}

	double N_nom_near = 0;
	double N_nom_far = 0;
	double N_ext_near = 0;
	double N_ext_far = 0;
	double efficiency_near = 0;
	double efficiency_far = 0;
	double efficiencyError_near = 0;
	double efficiencyError_far = 0;

	double N_nom_near_zVSr2[Nr2Points][NzPoints];
	double N_nom_far_zVSr2[Nr2Points][NzPoints];
	double N_ext_near_zVSr2[Nr2Points][NzPoints];
	double N_ext_far_zVSr2[Nr2Points][NzPoints];
	double efficiency_near_zVSr2[Nr2Points][NzPoints];
	double efficiency_far_zVSr2[Nr2Points][NzPoints];
	double efficiencyError_near_zVSr2[Nr2Points][NzPoints];
	double efficiencyError_far_zVSr2[Nr2Points][NzPoints];

	for(int ir2 = 0; ir2 < Nr2Points; ir2++){
		for(int iz = 0; iz < NzPoints; iz++){
			N_nom_near_zVSr2[ir2][iz] = 0;
			N_nom_far_zVSr2[ir2][iz] = 0;
			N_ext_near_zVSr2[ir2][iz] = 0;
			N_ext_far_zVSr2[ir2][iz] = 0;
			efficiency_near_zVSr2[ir2][iz] = 0;
			efficiency_far_zVSr2[ir2][iz] = 0;
			efficiencyError_near_zVSr2[ir2][iz] = 0;
			efficiencyError_far_zVSr2[ir2][iz] = 0;
		}
	}


//Near
	N_nom_near = Ed_sub_near->Integral(Ed_sub_near->FindBin(peak_near-3*sigma_near), Ed_sub_near->FindBin(peak_near+3*sigma_near));
	N_ext_near = Ed_sub_near->Integral(Ed_sub_near->FindBin(Elow), Ed_sub_near->FindBin(Ehigh));
	efficiency_near = N_nom_near/N_ext_near;
	//Error
	efficiencyError_near = sqrt(pow((Ed_sub_near->Integral(Ed_sub_near->FindBin(peak_near-3*sigma_near), Ed_sub_near->FindBin(peak_near+3*sigma_near)))/pow(Ed_sub_near->Integral(Ed_sub_near->FindBin(Elow), Ed_sub_near->FindBin(Ehigh)),2),2)*(Ed_ibd_near->Integral(Ed_ibd_near->FindBin(Elow), Ed_ibd_near->FindBin(peak_near-3*sigma_near))+Ed_ibd_near->Integral(Ed_ibd_near->FindBin(peak_near-3*sigma_near), Ed_ibd_near->FindBin(Ehigh))) + pow(1./(Ed_sub_near->Integral(Ed_sub_near->FindBin(Elow), Ed_sub_near->FindBin(Ehigh)))-(Ed_sub_near->Integral(Ed_sub_near->FindBin(peak_near-3*sigma_near), Ed_sub_near->FindBin(peak_near+3*sigma_near)))/pow(Ed_sub_near->Integral(Ed_sub_near->FindBin(Elow), Ed_sub_near->FindBin(Ehigh)),2),2)*Ed_ibd_near->Integral(Ed_ibd_near->FindBin(peak_near-3*sigma_near), Ed_ibd_near->FindBin(peak_near+ 3*sigma_near)));

//Far
	N_nom_far = Ed_sub_far->Integral(Ed_sub_far->FindBin(peak_far-3*sigma_far), Ed_sub_far->FindBin(peak_far+3*sigma_far));
	N_ext_far = Ed_sub_far->Integral(Ed_sub_far->FindBin(Elow_far), Ed_sub_far->FindBin(Ehigh_far));
	efficiency_far = N_nom_far/N_ext_far;
	//Error
	efficiencyError_far = sqrt(pow((Ed_sub_far->Integral(Ed_sub_far->FindBin(peak_far-3*sigma_far), Ed_sub_far->FindBin(peak_far+3*sigma_far)))/pow(Ed_sub_far->Integral(Ed_sub_far->FindBin(Elow_far), Ed_sub_far->FindBin(Ehigh_far)),2),2)*(Ed_ibd_far->Integral(Ed_ibd_far->FindBin(Elow_far), Ed_ibd_far->FindBin(peak_far-3*sigma_far))+Ed_ibd_far->Integral(Ed_ibd_far->FindBin(peak_far-3*sigma_far), Ed_ibd_far->FindBin(Ehigh_far))) + pow(1./(Ed_sub_far->Integral(Ed_sub_far->FindBin(Elow_far), Ed_sub_far->FindBin(Ehigh_far)))-(Ed_sub_far->Integral(Ed_sub_far->FindBin(peak_far-3*sigma_far), Ed_sub_far->FindBin(peak_far+3*sigma_far)))/pow(Ed_sub_far->Integral(Ed_sub_far->FindBin(Elow_far), Ed_sub_far->FindBin(Ehigh_far)),2),2)*Ed_ibd_far->Integral(Ed_ibd_far->FindBin(peak_far-3*sigma_far), Ed_ibd_far->FindBin(peak_far+ 3*sigma_far)));

//Efficiencies and errors for voxels
	for(int ir2 = 0; ir2 < Nr2Points; ir2++){
		for(int iz = 0; iz < NzPoints; iz++){
		//Near
			N_nom_near_zVSr2[ir2][iz] = Ed_rate_zVSr2_near[ir2][iz]->Integral(Ed_rate_zVSr2_near[ir2][iz]->FindBin(peak_near-3*sigma_near), Ed_rate_zVSr2_near[ir2][iz]->FindBin(peak_near+3*sigma_near));
			N_ext_near_zVSr2[ir2][iz] = Ed_rate_zVSr2_near[ir2][iz]->Integral(Ed_rate_zVSr2_near[ir2][iz]->FindBin(Elow), Ed_rate_zVSr2_near[ir2][iz]->FindBin(Ehigh));
			efficiency_near_zVSr2[ir2][iz] = N_nom_near_zVSr2[ir2][iz]/N_ext_near_zVSr2[ir2][iz];
			//Error
			efficiencyError_near_zVSr2[ir2][iz] = sqrt(pow((Ed_rate_zVSr2_near[ir2][iz]->Integral(Ed_rate_zVSr2_near[ir2][iz]->FindBin(peak_near-3*sigma_near), Ed_rate_zVSr2_near[ir2][iz]->FindBin(peak_near+3*sigma_near)))/pow(Ed_rate_zVSr2_near[ir2][iz]->Integral(Ed_rate_zVSr2_near[ir2][iz]->FindBin(Elow), Ed_rate_zVSr2_near[ir2][iz]->FindBin(Ehigh)),2),2)*(Ed_ibd_zVSr2_near[ir2][iz]->Integral(Ed_ibd_zVSr2_near[ir2][iz]->FindBin(Elow), Ed_ibd_zVSr2_near[ir2][iz]->FindBin(peak_near-3*sigma_near))+Ed_ibd_zVSr2_near[ir2][iz]->Integral(Ed_ibd_zVSr2_near[ir2][iz]->FindBin(peak_near-3*sigma_near), Ed_ibd_zVSr2_near[ir2][iz]->FindBin(Ehigh))) + pow(1./(Ed_rate_zVSr2_near[ir2][iz]->Integral(Ed_rate_zVSr2_near[ir2][iz]->FindBin(Elow), Ed_rate_zVSr2_near[ir2][iz]->FindBin(Ehigh)))-(Ed_rate_zVSr2_near[ir2][iz]->Integral(Ed_rate_zVSr2_near[ir2][iz]->FindBin(peak_near-3*sigma_near), Ed_rate_zVSr2_near[ir2][iz]->FindBin(peak_near+3*sigma_near)))/pow(Ed_rate_zVSr2_near[ir2][iz]->Integral(Ed_rate_zVSr2_near[ir2][iz]->FindBin(Elow), Ed_rate_zVSr2_near[ir2][iz]->FindBin(Ehigh)),2),2)*Ed_ibd_zVSr2_near[ir2][iz]->Integral(Ed_ibd_zVSr2_near[ir2][iz]->FindBin(peak_near-3*sigma_near), Ed_ibd_zVSr2_near[ir2][iz]->FindBin(peak_near+ 3*sigma_near)));
		//Far
			N_nom_far_zVSr2[ir2][iz] = Ed_rate_zVSr2_far[ir2][iz]->Integral(Ed_rate_zVSr2_far[ir2][iz]->FindBin(peak_far-3*sigma_far), Ed_rate_zVSr2_far[ir2][iz]->FindBin(peak_far+3*sigma_far));
			N_ext_far_zVSr2[ir2][iz] = Ed_rate_zVSr2_far[ir2][iz]->Integral(Ed_rate_zVSr2_far[ir2][iz]->FindBin(Elow_far), Ed_rate_zVSr2_far[ir2][iz]->FindBin(Ehigh_far));
			efficiency_far_zVSr2[ir2][iz] = N_nom_far_zVSr2[ir2][iz]/N_ext_far_zVSr2[ir2][iz];
			//Error
			efficiencyError_far_zVSr2[ir2][iz] = sqrt(pow((Ed_rate_zVSr2_far[ir2][iz]->Integral(Ed_rate_zVSr2_far[ir2][iz]->FindBin(peak_far-3*sigma_far), Ed_rate_zVSr2_far[ir2][iz]->FindBin(peak_far+3*sigma_far)))/pow(Ed_rate_zVSr2_far[ir2][iz]->Integral(Ed_rate_zVSr2_far[ir2][iz]->FindBin(Elow_far), Ed_rate_zVSr2_far[ir2][iz]->FindBin(Ehigh_far)),2),2)*(Ed_ibd_zVSr2_far[ir2][iz]->Integral(Ed_ibd_zVSr2_far[ir2][iz]->FindBin(Elow_far), Ed_ibd_zVSr2_far[ir2][iz]->FindBin(peak_far-3*sigma_far))+Ed_ibd_zVSr2_far[ir2][iz]->Integral(Ed_ibd_zVSr2_far[ir2][iz]->FindBin(peak_far-3*sigma_far), Ed_ibd_zVSr2_far[ir2][iz]->FindBin(Ehigh_far))) + pow(1./(Ed_rate_zVSr2_far[ir2][iz]->Integral(Ed_rate_zVSr2_far[ir2][iz]->FindBin(Elow_far), Ed_rate_zVSr2_far[ir2][iz]->FindBin(Ehigh_far)))-(Ed_rate_zVSr2_far[ir2][iz]->Integral(Ed_rate_zVSr2_far[ir2][iz]->FindBin(peak_far-3*sigma_far), Ed_rate_zVSr2_far[ir2][iz]->FindBin(peak_far+3*sigma_far)))/pow(Ed_rate_zVSr2_far[ir2][iz]->Integral(Ed_rate_zVSr2_far[ir2][iz]->FindBin(Elow_far), Ed_rate_zVSr2_far[ir2][iz]->FindBin(Ehigh_far)),2),2)*Ed_ibd_zVSr2_far[ir2][iz]->Integral(Ed_ibd_zVSr2_far[ir2][iz]->FindBin(peak_far-3*sigma_far), Ed_ibd_zVSr2_far[ir2][iz]->FindBin(peak_far+ 3*sigma_far)));
		}
	}




//Plot the efficiencies: whole ADs
	TH1F* h_eff_whole = new TH1F("h_eff_whole","h_eff_whole",2,0.5,2.5);
	h_eff_whole->SetBinContent(h_eff_whole->FindBin(1),efficiency_near);
	h_eff_whole->SetBinContent(h_eff_whole->FindBin(2),efficiency_far);
	h_eff_whole->SetBinError(h_eff_whole->FindBin(1),efficiencyError_near);
	h_eff_whole->SetBinError(h_eff_whole->FindBin(2),efficiencyError_far);
	h_eff_whole->GetXaxis()->SetTitle("Near, Far");
	h_eff_whole->GetYaxis()->SetTitle("Delayed Energy Cut Efficiency");

	TCanvas *eff_whole_comp = new TCanvas("eff_whole_comp","eff_whole_comp");
		eff_whole_comp->cd();
		h_eff_whole->SetStats(0);
		h_eff_whole->Draw();
	cout << "Near efficiency for whole detector: " << efficiency_near << endl << "Far efficiency for whole detector: " << efficiency_far << endl;
	cout << "% difference: " << 100*(efficiency_far - efficiency_near)/efficiency_near << endl;

//Plot the efficiencies: voxels
	TH1F* h_eff_voxels_near = new TH1F("h_eff_voxels_near","h_eff_voxels_near",Nr2Points*NzPoints,0.5,Nr2Points*NzPoints+0.5);
	TH1F* h_eff_voxels_far = new TH1F("h_eff_voxels_far","h_eff_voxels_far",Nr2Points*NzPoints,0.5,Nr2Points*NzPoints+0.5);
	for(int ir2 = 0; ir2 < Nr2Points; ir2++){
		for(int iz = 0; iz < NzPoints; iz++){
			if(ir2 == Nr2Points-1) continue;
			h_eff_voxels_near->SetBinContent(h_eff_voxels_near->FindBin(iz*Nr2Points+ir2+1),efficiency_near_zVSr2[ir2][iz]);
			h_eff_voxels_near->SetBinError(h_eff_voxels_near->FindBin(iz*Nr2Points+ir2+1),efficiencyError_near_zVSr2[ir2][iz]);
			h_eff_voxels_far->SetBinContent(h_eff_voxels_far->FindBin(iz*Nr2Points+ir2+1),efficiency_far_zVSr2[ir2][iz]);
			h_eff_voxels_far->SetBinError(h_eff_voxels_far->FindBin(iz*Nr2Points+ir2+1),efficiencyError_far_zVSr2[ir2][iz]);
		}
	}

	h_eff_voxels_near->GetXaxis()->SetTitle("Voxel ID");
	h_eff_voxels_near->GetYaxis()->SetTitle("Delayed Energy Cut Efficiency");
	h_eff_voxels_near->SetLineColor(kBlack);
	h_eff_voxels_far->SetLineColor(kRed);
	h_eff_voxels_near->SetMarkerStyle(22);
	h_eff_voxels_far->SetMarkerStyle(22);
	h_eff_voxels_near->SetMarkerColor(kBlack);
	h_eff_voxels_far->SetMarkerColor(kRed);

	TCanvas *eff_voxels_comp = new TCanvas("eff_voxels_comp","eff_voxels_comp");
		eff_voxels_comp->cd();
		h_eff_voxels_near->SetStats(0);
		h_eff_voxels_near->Draw();
		h_eff_voxels_far->Draw("same");
		eff_voxels_comp->BuildLegend();

	TH1F* h_eff_voxels_diff = new TH1F("h_eff_voxels_diff","h_eff_voxels_diff",Nr2Points*NzPoints,0.5,Nr2Points*NzPoints+0.5);
	for(int ir2 = 0; ir2 < Nr2Points; ir2++){
		for(int iz = 0; iz < NzPoints; iz++){
			if(ir2 == Nr2Points-1) continue;
			h_eff_voxels_diff->SetBinContent(h_eff_voxels_diff->FindBin(iz*Nr2Points+ir2+1),100*(efficiency_far_zVSr2[ir2][iz]-efficiency_near_zVSr2[ir2][iz])/efficiency_near_zVSr2[ir2][iz]);
			h_eff_voxels_diff->SetBinError(h_eff_voxels_diff->FindBin(iz*Nr2Points+ir2+1),-1*h_eff_voxels_diff->GetBinContent(h_eff_voxels_diff->FindBin(iz*Nr2Points+ir2+1))*sqrt(pow((efficiencyError_far_zVSr2[ir2][iz]/efficiency_far_zVSr2[ir2][iz]),2)+pow((efficiencyError_near_zVSr2[ir2][iz]/efficiency_near_zVSr2[ir2][iz]),2)));
		}
	}

	h_eff_voxels_diff->GetXaxis()->SetTitle("Voxel ID");
	h_eff_voxels_diff->GetYaxis()->SetTitle("% Difference");
	h_eff_voxels_diff->SetLineColor(kBlack);
	h_eff_voxels_diff->SetMarkerStyle(22);
	h_eff_voxels_diff->SetMarkerColor(kBlack);


	TCanvas *eff_voxels_diff = new TCanvas("eff_voxels_diff","eff_voxels_diff");
		eff_voxels_diff->cd();
		h_eff_voxels_diff->SetStats(0);
		h_eff_voxels_diff->Draw();

//Plot the shapes of the near and far spectra normalized to +- 3sigma
	TCanvas *Ed_totAD_3sigNormalized = new TCanvas("Ed_totAD_3sigNormalized","Ed_totAD_3sigNormalized");
		Ed_totAD_3sigNormalized->cd();
		Ed_sub_near->Scale(1./(Ed_sub_near->Integral(Ed_sub_near->FindBin(peak_near-3*sigma_near),Ed_sub_near->FindBin(peak_near+3*sigma_near))));
		Ed_sub_far->Scale(1./(Ed_sub_far->Integral(Ed_sub_far->FindBin(peak_far-3*sigma_far),Ed_sub_far->FindBin(peak_far+3*sigma_far))));
		Ed_sub_near->SetLineColor(kBlack);
		Ed_sub_far->SetLineColor(kRed);
		Ed_sub_near->Rebin(20);
		Ed_sub_far->Rebin(20);
		Ed_sub_near->SetLineWidth(3);
		Ed_sub_far->SetLineWidth(3);
		Ed_sub_near->GetXaxis()->SetRangeUser(1.5,3);
		Ed_sub_near->SetStats(0);
		Ed_sub_near->Draw();
		Ed_sub_far->Draw("same");
			TLine *lin_3sigLow_near = new TLine(peak_near - 3*sigma_near,0,peak_near - 3*sigma_near,0.3);
			TLine *lin_3sigHigh_near = new TLine(peak_near + 3*sigma_near,0,peak_near + 3*sigma_near,0.3);
			TLine *lin_3sigLow_far = new TLine(peak_far - 3*sigma_far,0,peak_far - 3*sigma_far,0.3);
			TLine *lin_3sigHigh_far = new TLine(peak_far + 3*sigma_far,0,peak_far + 3*sigma_far,0.3);
			lin_3sigLow_near->SetLineColor(kBlack);
			lin_3sigHigh_near->SetLineColor(kBlack);
			lin_3sigLow_far->SetLineColor(kRed);
			lin_3sigHigh_far->SetLineColor(kRed);
			lin_3sigLow_near->SetLineStyle(2);
			lin_3sigHigh_near->SetLineStyle(2);
			lin_3sigLow_far->SetLineStyle(2);
			lin_3sigHigh_far->SetLineStyle(2);
			lin_3sigLow_near->SetLineWidth(2);
			lin_3sigHigh_near->SetLineWidth(2);
			lin_3sigLow_far->SetLineWidth(2);
			lin_3sigHigh_far->SetLineWidth(2);
			lin_3sigLow_near->Draw("same");
			lin_3sigHigh_near->Draw("same");
			lin_3sigLow_far->Draw("same");
			lin_3sigHigh_far->Draw("same");
		TLine *lin_Elow_near = new TLine(Elow,0,Elow,0.3);
		TLine *lin_Ehigh_near = new TLine(Ehigh,0,Ehigh,0.3);
		TLine *lin_Elow_far = new TLine(Elow_far,0,Elow_far,0.3);
		TLine *lin_Ehigh_far = new TLine(Ehigh_far,0,Ehigh_far,0.3);
		lin_Elow_near->SetLineColor(kBlack);
		lin_Ehigh_near->SetLineColor(kBlack);
		lin_Elow_far->SetLineColor(kRed);
		lin_Ehigh_far->SetLineColor(kRed);
		lin_Elow_near->SetLineStyle(1);
		lin_Ehigh_near->SetLineStyle(1);
		lin_Elow_far->SetLineStyle(1);
		lin_Ehigh_far->SetLineStyle(1);
		lin_Elow_near->SetLineWidth(2);
		lin_Ehigh_near->SetLineWidth(2);
		lin_Elow_far->SetLineWidth(2);
		lin_Ehigh_far->SetLineWidth(2);
		lin_Elow_near->Draw("same");
		lin_Ehigh_near->Draw("same");
		lin_Elow_far->Draw("same");
		lin_Ehigh_far->Draw("same");

	cout << "Shift of peak value: " << peak_far - peak_near << endl;

//Plot the individual voxels in the same way:
	TCanvas *Ed_vox_3sigNorm[Nr2Points*NzPoints];
	for(int iz = 0; iz < NzPoints; iz++){
		for(int ir2 = 0; ir2 < Nr2Points; ir2++){
			Ed_vox_3sigNorm[(iz*Nr2Points)+ir2] = new TCanvas(Form("Ed_vox_3sigNorm_%d",(iz*Nr2Points)+ir2+1),Form("Ed_vox_3sigNorm_%d",(iz*Nr2Points)+ir2+1));
			Ed_vox_3sigNorm[(iz*Nr2Points)+ir2]->cd();

		Ed_rate_zVSr2_near[ir2][iz]->Scale(1./(Ed_rate_zVSr2_near[ir2][iz]->Integral(Ed_rate_zVSr2_near[ir2][iz]->FindBin(peak_near-3*sigma_near),Ed_rate_zVSr2_near[ir2][iz]->FindBin(peak_near+3*sigma_near))));
		Ed_rate_zVSr2_far[ir2][iz]->Scale(1./(Ed_rate_zVSr2_far[ir2][iz]->Integral(Ed_rate_zVSr2_far[ir2][iz]->FindBin(peak_far-3*sigma_far),Ed_rate_zVSr2_far[ir2][iz]->FindBin(peak_far+3*sigma_far))));
		Ed_rate_zVSr2_near[ir2][iz]->SetLineColor(kBlack);
		Ed_rate_zVSr2_far[ir2][iz]->SetLineColor(kRed);
		Ed_rate_zVSr2_near[ir2][iz]->Rebin(10);
		Ed_rate_zVSr2_far[ir2][iz]->Rebin(10);
		Ed_rate_zVSr2_near[ir2][iz]->SetLineWidth(3);
		Ed_rate_zVSr2_far[ir2][iz]->SetLineWidth(3);
		Ed_rate_zVSr2_near[ir2][iz]->GetXaxis()->SetRangeUser(1.5,3);
		Ed_rate_zVSr2_near[ir2][iz]->SetStats(0);
		Ed_rate_zVSr2_near[ir2][iz]->Draw();
		Ed_rate_zVSr2_far[ir2][iz]->Draw("same");

			lin_3sigLow_near->Draw("same");
			lin_3sigHigh_near->Draw("same");
			lin_3sigLow_far->Draw("same");
			lin_3sigHigh_far->Draw("same");
			lin_Elow_near->Draw("same");
			lin_Ehigh_near->Draw("same");
			lin_Elow_far->Draw("same");
			lin_Ehigh_far->Draw("same");
			Ed_vox_3sigNorm[(iz*Nr2Points)+ir2]->SetGridy();
		}
	}

}









