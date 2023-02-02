///////////////////////////////////////////////////////////////////////////////////////
// LOG ////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//  
//  Original program is from T. Kikawa.
//  Arranged by A. Hiramot.
//  
// 2020/07/28 H.Oshima
//  Arranged for run6 analysis.
// 
// 
///////////////////////////////////////////////////////////////////////////////////////


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<iostream>
#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <TPaletteAxis.h>
#include <TList.h>
#include <TRint.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TPave.h>
#include <TAttPad.h> 
#include <TTree.h>
#include <TFile.h>
#include <TArrow.h>
#include <TLegend.h>
#include <TText.h>
#include <TLine.h>
#include <THStack.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include <TRandom3.h>
#include <TMatrixDSym.h>

using namespace std;


// Covariance matrix bin //
#define nbin 40




int main(int argc, char *argv[]){
	
	TCanvas *can_dummy = new TCanvas("can_dummy","",1600,1400);
	
/* ---------------------------------------------------------------------------------------------------------------------------- */
	
	/* Setting of gStyle */
	//gStyle->SetOptStat(11111111);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetPalette(1);
	//gStyle->SetPadGridX(1);
	//gStyle->SetPadGridY(1);
	gStyle->SetStatY(1.00);
	gStyle->SetStatX(1.00);
	gStyle->SetStatW(0.09);
	gStyle->SetStatH(0.10);
	
	
/* ---------------------------------------------------------------------------------------------------------------------------- */
	
	
	// neutrino beam mode //
	int beam = 0;
	bool flg_b = false;
	
	// neutrino component //
	int neutrino = 0;
	bool flg_n = false;
	
	//  input file name
	char root_filename_i1[1024];
	bool flg_i = false;
	
	// output file name
	char root_filename_o1[1024];
	bool flg_o = false;
	
	// output PDF file name
	char pdf_filename_o1[1024];
	bool flg_d = false;
	
	// Module (ND) //
	int nd = 0;
	bool flg_m = false;
	
	// interaction target //
	char target_name[1024];
	bool flg_a = false;
	
	
	
	// option //
	int c = -1;
	
	while((c=getopt(argc, argv, "i:o:f:b:n:m:a:d:")) != -1){
		
		switch(c){
			case 'b':
				beam = atoi(optarg);
				flg_b = true;
				break;
			case 'n':
				neutrino = atoi(optarg);
				flg_n = true;
				break;
			case 'i':
				sprintf(root_filename_i1,"%s",optarg);
				flg_i = true;
				break;
			case 'o':
				sprintf(root_filename_o1,"%s",optarg);
				flg_o = true;
				break;
			case 'm':
				nd = atoi(optarg);
				flg_m = true;
				break;
			case 'a':
				sprintf(target_name,"%s",optarg);
				flg_a = true;
				break;
			case 'd':
				sprintf(pdf_filename_o1,"%s",optarg);
				flg_d = true;
				break;
		}
		
	}
	
	
	
	// Beam Mode //
	if(!flg_b){
		cerr << "Error : -b option is need for neutrino beam mode." << endl;
		cerr << "1 ... Neutrino Beam Mode,  -1 ... Anti-neutrino Beam Mode" << endl;
		exit(1);
	}
	if(beam==1){
		cerr << "Neutrino Beam Mode" << endl;
	}
	else if(beam==-1){
		cerr << "Anti-Neutrino Beam Mode" << endl;
	}
	else{
		cerr << "Error : -b option must be 1/-1." << endl;
		cerr << "1 ... Neutrino Beam Mode,  -1 ... Anti-neutrino Beam Mode" << endl;
		exit(1);
	}
	
	
	// Neutrino component //
	char name_neutrino[128];
	if(!flg_n){
		cerr << "Error : -n option is need for neutrino component." << endl;
		cerr << "1 ... Neutrino component,  -1 ... Anti-neutrino component" << endl;
		exit(1);
	}
	if(neutrino==1){
		cerr << "Neutrino component" << endl;
		sprintf(name_neutrino,"numu");
	}
	else if(neutrino==-1){
		cerr << "Anti-Neutrino component" << endl;
		sprintf(name_neutrino,"numubar");
	}
	else if(neutrino==0){
		cerr << "Neutrino + Anti-Neutrino component" << endl;
		sprintf(name_neutrino,"numu_numubar");
	}
	else{
		cerr << "Error : -n option must be 1/-1." << endl;
		cerr << "1 ... Neutrino component,  -1 ... Anti-neutrino component" << endl;
		exit(1);
	}
	
	
	// EVENT ROOT File //
	
	// intput event root file //
	if(!flg_i){
		cerr << "Error : -i option is need for input event root file." << endl;
		exit(1);
	}
	
	
	// interaction target //
	if(!flg_a){
		cerr << "Error : -a option is need for interaction target material." << endl;
		exit(1);
	}
	
	
	// Neutrino Detector //
	if(!flg_m){
		cerr << "Error : -m option is need for ND." << endl;
		exit(1);
	}
	
	
	TFile *rootfile_i1 = new TFile(root_filename_i1,"read");  // Open the Geant output ROOT file
	if (rootfile_i1->IsZombie()){
		cerr << "Failed to open " << root_filename_i1 << "." << endl;
		return 0;
	}
	
	TTree* tree = (TTree*)rootfile_i1->Get("tree");
	
	// Flux Uncertainties //
	const int nbin_mul              =  10;
	const int nbin_mul_pipm         =  10;
	const int nbin_mul_proton       =  10;
	const int nbin_mul_pipm_ppi     =   5;
	const int nbin_mul_proton_ppi   =   5;
	const int nbin_angle_muon       =  36; //   5 deg.  bin
	const int nbin_angle_pipm       =  18; //  10 deg.  bin
	const int nbin_angle_proton     =  18; //  10 deg.  bin
	const int nbin_angle_cos_muon   = 100; // 0.02      bin
	const int nbin_angle_cos_pipm   =  20; // 0.10      bin
	const int nbin_angle_cos_proton =  20; // 0.10      bin
	const int nbin_p_muon           = 250; // 0.2 GeV/c bin
	const int nbin_p_pipm           = 500; // 0.1 GeV/c bin
	const int nbin_p_proton         = 500; // 0.1 GeV/c bin
	
	double fluxerr_mul_total[nbin_mul]                                  = {0.0};
	double fluxerr_mul_pipm_total[nbin_mul_pipm]                        = {0.0};
	double fluxerr_mul_proton_total[nbin_mul_proton]                    = {0.0};
	double fluxerr_angle_corr_deg_muon_total[nbin_angle_muon]           = {0.0};
	double fluxerr_angle_corr_deg_pipm_total[nbin_angle_pipm]           = {0.0};
	double fluxerr_angle_corr_deg_proton_total[nbin_angle_proton]       = {0.0};
	double fluxerr_angle_corr_cos_muon_total[nbin_angle_cos_muon]       = {0.0};
	double fluxerr_angle_corr_cos_pipm_total[nbin_angle_cos_pipm]       = {0.0};
	double fluxerr_angle_corr_cos_proton_total[nbin_angle_cos_proton]   = {0.0};
	double fluxerr_p_muon_total[nbin_p_muon]                            = {0.0};
	double fluxerr_p_pipm_total[nbin_p_pipm]                            = {0.0};
	double fluxerr_p_proton_total[nbin_p_proton]                        = {0.0};
	
	double fluxerr_mul_proton_Npipm_total[nbin_mul_pipm_ppi][nbin_mul_proton_ppi];
	
	
	
	// Tree Branch //
	tree->SetBranchAddress(Form("fluxerr_mul_%s_total", name_neutrino),                       fluxerr_mul_total);
	tree->SetBranchAddress(Form("fluxerr_mul_pipm_%s_total", name_neutrino),                  fluxerr_mul_pipm_total);
	tree->SetBranchAddress(Form("fluxerr_mul_proton_%s_total", name_neutrino),                fluxerr_mul_proton_total);
	tree->SetBranchAddress(Form("fluxerr_angle_corr_deg_muon_%s_total", name_neutrino),       fluxerr_angle_corr_deg_muon_total);
	tree->SetBranchAddress(Form("fluxerr_angle_corr_deg_pipm_%s_total", name_neutrino),       fluxerr_angle_corr_deg_pipm_total);
	tree->SetBranchAddress(Form("fluxerr_angle_corr_deg_proton_%s_total", name_neutrino),     fluxerr_angle_corr_deg_proton_total);
	tree->SetBranchAddress(Form("fluxerr_angle_corr_cos_muon_%s_total", name_neutrino),       fluxerr_angle_corr_cos_muon_total);
	tree->SetBranchAddress(Form("fluxerr_angle_corr_cos_pipm_%s_total", name_neutrino),       fluxerr_angle_corr_cos_pipm_total);
	tree->SetBranchAddress(Form("fluxerr_angle_corr_cos_proton_%s_total", name_neutrino),     fluxerr_angle_corr_cos_proton_total);
	tree->SetBranchAddress(Form("fluxerr_p_muon_%s_total", name_neutrino),                    fluxerr_p_muon_total);
	tree->SetBranchAddress(Form("fluxerr_p_pipm_%s_total", name_neutrino),                    fluxerr_p_pipm_total);
	tree->SetBranchAddress(Form("fluxerr_p_proton_%s_total", name_neutrino),                  fluxerr_p_proton_total);
	tree->SetBranchAddress(Form("fluxerr_mul_proton_Npipm_%s_total", name_neutrino),          fluxerr_mul_proton_Npipm_total);
	
	
	
	
	// output event root file //
	if(!flg_o){
		cerr << "Error : -o option is need for output event root file." << endl;
		exit(1);
	}
	TFile *root_o1 = new TFile(root_filename_o1,"recreate");  // Open the Geant output ROOT file
	if (root_o1->IsZombie()){
		cerr << "Failed to open " << root_filename_o1 << "." << endl;
		return 0;
	}
	
	
	// multiplicity //
	// total			>= 5	>= 5
	// pion				>= 3	>= 3
	double xbins_mul[]        = {-0.5,0.5,1.5,2.5,3.5,4.5,10.5};
	double xbins_mul_pipm[]   = {-0.5,0.5,1.5,2.5,10.5};
	double xbins_mul_proton[] = {-0.5,0.5,1.5,2.5,10.5};
	
	// proton			>= 3	>= 3
	// proton vs. pi	>= 3,3	>= 3,3
	double xbins_mul_ppi_pipm[]   = {-0.5,0.5,1.5,2.5,10.5};
	double xbins_mul_ppi_proton[] = {-0.5,0.5,1.5,2.5,10.5};
	
	// angle (deg.) //
	// muon				>=  45,  5deg	>=     45,  5deg
	// pion				60-120, 10deg	>= 60-120, 10deg
	// proton			60-120, 10deg	>= 60-120, 10deg
	double xbins_angle_corr_deg_muon[]   = {0.0,5.0,10.0,15.0,20.0,25.0,30.0,90.0,180.0};
	double xbins_angle_corr_deg_pipm[]   = {0.0,20.0,40.0,90.0,140.0,160.0,180.0};
	double xbins_angle_corr_deg_proton[] = {0.0,20.0,40.0,90.0,140.0,160.0,180.0};
	
	// angle (cos.) //
	// muon				< 0.80,    0.04		< 0.80,    0.04
	// pion				-0.5-+0.5, 0.10		-0.5-+0.5, 0.10
	// proton			-0.5-+0.5, 0.10		-0.5-+0.5, 0.10
	double xbins_angle_corr_cos_muon[]   = {-1.00,0.00,0.80,0.84,0.88,0.92,0.96,1.00};
	double xbins_angle_corr_cos_pipm[]   = {-1.00,-0.90,-0.80,-0.20,0.20,0.80,0.90,1.00};
	double xbins_angle_corr_cos_proton[] = {-1.00,-0.90,-0.80,-0.20,0.20,0.80,0.90,1.00};
	
	// momentum //
	// muon				>= 4.0GeV/c, 0.4GeV/c		>= 4.0GeV/c, 0.4GeV/c
	// pion				>= 1.4GeV/c, 0.2GeV/c		>= 1.4GeV/c, 0.2GeV/c
	// proton			>= 1.4GeV/c, 0.2GeV/c		>= 1.4GeV/c, 0.2GeV/c
	double xbins_p_muon[]   = {0.0,0.4,0.8,1.2,1.6,2.0,5.0,50.0};
	double xbins_p_pipm[]   = {0.0,0.2,0.4,0.6,1.0,2.0,50.0};
	double xbins_p_proton[] = {0.0,0.4,0.6,0.8,1.0,2.0,50.0};
	
	int   nbinsx_mul                   = sizeof(xbins_mul)/sizeof(double);
	int   nbinsx_mul_pipm              = sizeof(xbins_mul_pipm)/sizeof(double);
	int   nbinsx_mul_proton            = sizeof(xbins_mul_proton)/sizeof(double);
	int   nbinsx_mul_ppi_pipm          = sizeof(xbins_mul_ppi_pipm)/sizeof(double);
	int   nbinsx_mul_ppi_proton        = sizeof(xbins_mul_ppi_proton)/sizeof(double);
	int   nbinsx_angle_corr_deg_muon   = sizeof(xbins_angle_corr_deg_muon)/sizeof(double);
	int   nbinsx_angle_corr_deg_pipm   = sizeof(xbins_angle_corr_deg_pipm)/sizeof(double);
	int   nbinsx_angle_corr_deg_proton = sizeof(xbins_angle_corr_deg_proton)/sizeof(double);
	int   nbinsx_angle_corr_cos_muon   = sizeof(xbins_angle_corr_cos_muon)/sizeof(double);
	int   nbinsx_angle_corr_cos_pipm   = sizeof(xbins_angle_corr_cos_pipm)/sizeof(double);
	int   nbinsx_angle_corr_cos_proton = sizeof(xbins_angle_corr_cos_proton)/sizeof(double);
	int   nbinsx_p_muon                = sizeof(xbins_p_muon)/sizeof(double);
	int   nbinsx_p_pipm                = sizeof(xbins_p_pipm)/sizeof(double);
	int   nbinsx_p_proton              = sizeof(xbins_p_proton)/sizeof(double);
	
	
	TH1D *h_fluxerr_mul_total        = new TH1D(Form("h_fluxerr_mul_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_mul-1,xbins_mul);
	TH1D *h_fluxerr_mul_pipm_total   = new TH1D(Form("h_fluxerr_mul_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_mul_pipm-1,xbins_mul_pipm);
	TH1D *h_fluxerr_mul_proton_total = new TH1D(Form("h_fluxerr_mul_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_mul_proton-1,xbins_mul_proton);
	
	TH2D *h_fluxerr_mul_proton_pipm_total  = new TH2D(Form("h_fluxerr_mul_proton_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_mul_ppi_proton-1,xbins_mul_ppi_proton,nbinsx_mul_ppi_pipm-1,xbins_mul_ppi_pipm);
	
	TH1D *h_fluxerr_angle_corr_deg_muon_total   = new TH1D(Form("h_fluxerr_angle_corr_deg_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_angle_corr_deg_muon-1,xbins_angle_corr_deg_muon);
	TH1D *h_fluxerr_angle_corr_deg_pipm_total   = new TH1D(Form("h_fluxerr_angle_corr_deg_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_angle_corr_deg_pipm-1,xbins_angle_corr_deg_pipm);
	TH1D *h_fluxerr_angle_corr_deg_proton_total = new TH1D(Form("h_fluxerr_angle_corr_deg_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_angle_corr_deg_proton-1,xbins_angle_corr_deg_proton);
	
	TH1D *h_fluxerr_angle_corr_cos_muon_total   = new TH1D(Form("h_fluxerr_angle_corr_cos_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_angle_corr_cos_muon-1,xbins_angle_corr_cos_muon);
	TH1D *h_fluxerr_angle_corr_cos_pipm_total   = new TH1D(Form("h_fluxerr_angle_corr_cos_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_angle_corr_cos_pipm-1,xbins_angle_corr_cos_pipm);
	TH1D *h_fluxerr_angle_corr_cos_proton_total = new TH1D(Form("h_fluxerr_angle_corr_cos_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_angle_corr_cos_proton-1,xbins_angle_corr_cos_proton);
	
	TH1D *h_fluxerr_p_muon_total   = new TH1D(Form("h_fluxerr_p_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_p_muon-1,xbins_p_muon);
	TH1D *h_fluxerr_p_pipm_total   = new TH1D(Form("h_fluxerr_p_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_p_pipm-1,xbins_p_pipm);
	TH1D *h_fluxerr_p_proton_total = new TH1D(Form("h_fluxerr_p_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbinsx_p_proton-1,xbins_p_proton);
	
	//TH1D *h_fluxerr_mul_total        = new TH1D(Form("h_fluxerr_mul_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_mul,-0.5,9.5);
	//TH1D *h_fluxerr_mul_pipm_total   = new TH1D(Form("h_fluxerr_mul_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_mul_pipm,-0.5,9.5);
	//TH1D *h_fluxerr_mul_proton_total = new TH1D(Form("h_fluxerr_mul_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_mul_proton,-0.5,9.5);
	//
	//TH2D *h_fluxerr_mul_proton_pipm_total  = new TH2D(Form("h_fluxerr_mul_proton_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_mul_proton_ppi,-0.5,4.5,nbin_mul_pipm_ppi,-0.5,4.5);
	//
	//TH1D *h_fluxerr_angle_corr_deg_muon_total   = new TH1D(Form("h_fluxerr_angle_corr_deg_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_angle_muon,0.0,180.0);
	//TH1D *h_fluxerr_angle_corr_deg_pipm_total   = new TH1D(Form("h_fluxerr_angle_corr_deg_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_angle_pipm,0.0,180.0);
	//TH1D *h_fluxerr_angle_corr_deg_proton_total = new TH1D(Form("h_fluxerr_angle_corr_deg_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_angle_proton,0.0,180.0);
	//
	//TH1D *h_fluxerr_angle_corr_cos_muon_total   = new TH1D(Form("h_fluxerr_angle_corr_cos_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_angle_cos_muon,-1.0,1.0);
	//TH1D *h_fluxerr_angle_corr_cos_pipm_total   = new TH1D(Form("h_fluxerr_angle_corr_cos_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_angle_cos_pipm,-1.0,1.0);
	//TH1D *h_fluxerr_angle_corr_cos_proton_total = new TH1D(Form("h_fluxerr_angle_corr_cos_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_angle_cos_proton,-1.0,1.0);
	//
	//TH1D *h_fluxerr_p_muon_total   = new TH1D(Form("h_fluxerr_p_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_p_muon,0.0,50.0);
	//TH1D *h_fluxerr_p_pipm_total   = new TH1D(Form("h_fluxerr_p_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_p_pipm,0.0,50.0);
	//TH1D *h_fluxerr_p_proton_total = new TH1D(Form("h_fluxerr_p_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name),"",nbin_p_proton,0.0,50.0);
	
	
	int NEVENT = tree->GetEntries();
	
	
	// Fit //
	TF1 *f_gaus = new TF1("f_gaus","gaus(0)");
	f_gaus->SetParName(0,"height");
	f_gaus->SetParName(1,"center");
	f_gaus->SetParName(2,"sigma");
	f_gaus->SetParLimits(0,0,1e7);
	f_gaus->SetParLimits(1,0.500,1.500);
	f_gaus->SetParLimits(2,0.000,0.500);
	f_gaus->SetLineColor(kBlack);
	
	
	// multiplicity //
	cerr << "multiplicity ..." << endl;
	for(int ix = 0; ix < nbinsx_mul-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_mul; i++){
			
			double mul = (double)i * 1.0;
			
			if(mul < xbins_mul[ix] || xbins_mul[ix+1] <= mul)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_mul_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		//cerr << endl;
		//cerr << mul << "     " << sigma << endl;
		//cerr << endl;
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_mul_total->Fill((xbins_mul[ix]+xbins_mul[ix+1])/2.0, 0.0);
		else
			h_fluxerr_mul_total->Fill((xbins_mul[ix]+xbins_mul[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	// Number of pions //
	cerr << "number of pions ..." << endl;
	for(int ix = 0; ix < nbinsx_mul_pipm-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_mul_pipm; i++){
			
			double mul = (double)i * 1.0;
			
			if(mul < xbins_mul_pipm[ix] || xbins_mul_pipm[ix+1] <= mul)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_mul_pipm_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_mul_pipm_total->Fill((xbins_mul_pipm[ix]+xbins_mul_pipm[ix+1])/2.0, 0.0);
		else
			h_fluxerr_mul_pipm_total->Fill((xbins_mul_pipm[ix]+xbins_mul_pipm[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	
	// Number of protons //
	cerr << "number of protons ..." << endl;
	for(int ix = 0; ix < nbinsx_mul_proton-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_mul_proton; i++){
			
			double mul = (double)i * 1.0;
			
			if(mul < xbins_mul_proton[ix] || xbins_mul_proton[ix+1] <= mul)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_mul_proton_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_mul_proton_total->Fill((xbins_mul_proton[ix]+xbins_mul_proton[ix+1])/2.0, 0.0);
		else
			h_fluxerr_mul_proton_total->Fill((xbins_mul_proton[ix]+xbins_mul_proton[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	
	// Number of protons vs. pions //
	cerr << "number of protons vs. pions ..." << endl;
	for(int ix = 0; ix < nbinsx_mul_ppi_pipm-1; ix++){
		
		for(int iy = 0; iy < nbinsx_mul_ppi_proton-1; iy++){
			
			TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
			
			for(int i = 0; i < nbin_mul_pipm_ppi; i++){
				
				double mul_pipm   = (double)i * 1.0;
				
				if(mul_pipm < xbins_mul_ppi_pipm[ix] || xbins_mul_ppi_pipm[ix+1] <= mul_pipm)
					continue;
				
				for(int j = 0; j < nbin_mul_proton_ppi; j++){
					
					double mul_proton = (double)j * 1.0;
					
					if(mul_proton < xbins_mul_ppi_proton[iy] || xbins_mul_ppi_proton[iy+1] <= mul_proton)
						continue;
					
					for(int k = 0; k < NEVENT; k++){
						
						if(k%10000==0) fprintf(stderr,"\r Bin:%3d, %3d  %10lld / %10lld", i, j, k, NEVENT);
						
						tree->GetEntry(k);
						
						tmp->Fill(fluxerr_mul_proton_Npipm_total[i][j]);
						
					}
					
				}
				
			}
			
			tmp->Fit(f_gaus,"Q","",0.0,2.0);
			
			double center = f_gaus->GetParameter(1);
			double sigma  = f_gaus->GetParameter(2);
			
			if(tmp->GetMean()<0.01)
				h_fluxerr_mul_proton_pipm_total->Fill((xbins_mul_ppi_proton[iy]+xbins_mul_ppi_proton[iy+1])/2.0, (xbins_mul_ppi_pipm[ix]+xbins_mul_ppi_pipm[ix+1])/2.0, 0.0);
			else
				h_fluxerr_mul_proton_pipm_total->Fill((xbins_mul_ppi_proton[iy]+xbins_mul_ppi_proton[iy+1])/2.0, (xbins_mul_ppi_pipm[ix]+xbins_mul_ppi_pipm[ix+1])/2.0, sigma);
			
			tmp->Delete();
			
		}
		
	}
	cerr << endl;
	
	
	// Angle of muons //
	cerr << "Muon angle ..." << endl;
	for(int ix = 0; ix < nbinsx_angle_corr_deg_muon-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_angle_muon; i++){
			
			double angle = (double)i * 5.0;
			
			if(angle < xbins_angle_corr_deg_muon[ix] || xbins_angle_corr_deg_muon[ix+1] <= angle)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_angle_corr_deg_muon_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_angle_corr_deg_muon_total->Fill((xbins_angle_corr_deg_muon[ix]+xbins_angle_corr_deg_muon[ix+1])/2.0, 0.0);
		else
			h_fluxerr_angle_corr_deg_muon_total->Fill((xbins_angle_corr_deg_muon[ix]+xbins_angle_corr_deg_muon[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	
	// Angle of pions //
	cerr << "Pion angle ..." << endl;
	for(int ix = 0; ix < nbinsx_angle_corr_deg_pipm-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_angle_pipm; i++){
			
			double angle = (double)i * 10.0;
			
			if(angle < xbins_angle_corr_deg_pipm[ix] || xbins_angle_corr_deg_pipm[ix+1] <= angle)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_angle_corr_deg_pipm_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_angle_corr_deg_pipm_total->Fill((xbins_angle_corr_deg_pipm[ix]+xbins_angle_corr_deg_pipm[ix+1])/2.0, 0.0);
		else
			h_fluxerr_angle_corr_deg_pipm_total->Fill((xbins_angle_corr_deg_pipm[ix]+xbins_angle_corr_deg_pipm[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	
	// Angle of protons //
	cerr << "Proton angle ..." << endl;
	for(int ix = 0; ix < nbinsx_angle_corr_deg_proton-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_angle_proton; i++){
			
			double angle = (double)i * 10.0;
			
			if(angle < xbins_angle_corr_deg_proton[ix] || xbins_angle_corr_deg_proton[ix+1] <= angle)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_angle_corr_deg_proton_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_angle_corr_deg_proton_total->Fill((xbins_angle_corr_deg_proton[ix]+xbins_angle_corr_deg_proton[ix+1])/2.0, 0.0);
		else
			h_fluxerr_angle_corr_deg_proton_total->Fill((xbins_angle_corr_deg_proton[ix]+xbins_angle_corr_deg_proton[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	
	// Angle of muons (cos) //
	cerr << "Muon angle (cos) ..." << endl;
	for(int ix = 0; ix < nbinsx_angle_corr_cos_muon-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_angle_cos_muon; i++){
			
			double angle = (double)i * 2.0/(double)nbin_angle_cos_muon - 1.0;
			
			if(angle < xbins_angle_corr_cos_muon[ix] || xbins_angle_corr_cos_muon[ix+1] <= angle)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_angle_corr_cos_muon_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_angle_corr_cos_muon_total->Fill((xbins_angle_corr_cos_muon[ix]+xbins_angle_corr_cos_muon[ix+1])/2.0, 0.0);
		else
			h_fluxerr_angle_corr_cos_muon_total->Fill((xbins_angle_corr_cos_muon[ix]+xbins_angle_corr_cos_muon[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	
	// Angle of pions (cos) //
	cerr << "Pion angle (cos) ..." << endl;
	for(int ix = 0; ix < nbinsx_angle_corr_cos_pipm-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_angle_cos_pipm; i++){
			
			double angle = (double)i * 2.0/(double)nbin_angle_cos_pipm - 1.0;
			
			if(angle < xbins_angle_corr_cos_pipm[ix] || xbins_angle_corr_cos_pipm[ix+1] <= angle)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_angle_corr_cos_pipm_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_angle_corr_cos_pipm_total->Fill((xbins_angle_corr_cos_pipm[ix]+xbins_angle_corr_cos_pipm[ix+1])/2.0, 0.0);
		else
			h_fluxerr_angle_corr_cos_pipm_total->Fill((xbins_angle_corr_cos_pipm[ix]+xbins_angle_corr_cos_pipm[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	
	// Angle of protons (cos) //
	cerr << "Proton angle (cos) ..." << endl;
	for(int ix = 0; ix < nbinsx_angle_corr_cos_proton-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_angle_cos_proton; i++){
			
			double angle = (double)i * 2.0/(double)nbin_angle_cos_proton - 1.0;
			
			if(angle < xbins_angle_corr_cos_proton[ix] || xbins_angle_corr_cos_proton[ix+1] <= angle)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_angle_corr_cos_proton_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_angle_corr_cos_proton_total->Fill((xbins_angle_corr_cos_proton[ix]+xbins_angle_corr_cos_proton[ix+1])/2.0, 0.0);
		else
			h_fluxerr_angle_corr_cos_proton_total->Fill((xbins_angle_corr_cos_proton[ix]+xbins_angle_corr_cos_proton[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		
	}
	cerr << endl;
	
	
	// Momentum of muons //
	cerr << "Muon momentum ..." << endl;
	for(int ix = 0; ix < nbinsx_p_muon-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_p_muon; i++){
			
			double momentum = (double)i * 0.20;
			
			if(momentum < xbins_p_muon[ix] || xbins_p_muon[ix+1] <= momentum)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_p_muon_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_p_muon_total->Fill((xbins_p_muon[ix]+xbins_p_muon[ix+1])/2.0, 0.0);
		else
			h_fluxerr_p_muon_total->Fill((xbins_p_muon[ix]+xbins_p_muon[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		//if(momentum > 5.0) break;
		
		
	}
	cerr << endl;
	
	
	// Momentum of pions //
	cerr << "Pion momentum ..." << endl;
	for(int ix = 0; ix < nbinsx_p_pipm-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_p_pipm; i++){
			
			double momentum = (double)i * 0.10;
			
			if(momentum < xbins_p_pipm[ix] || xbins_p_pipm[ix+1] <= momentum)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_p_pipm_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_p_pipm_total->Fill((xbins_p_pipm[ix]+xbins_p_pipm[ix+1])/2.0, 0.0);
		else
			h_fluxerr_p_pipm_total->Fill((xbins_p_pipm[ix]+xbins_p_pipm[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		//if(momentum > 5.0) break;
		
		
	}
	cerr << endl;
	
	
	// Momentum of protons //
	cerr << "Proton momentum ..." << endl;
	for(int ix = 0; ix < nbinsx_p_muon-1; ix++){
		
		TH1D *tmp = new TH1D("tmp","",200,0.0,2.0);
		
		for(int i = 0; i < nbin_p_proton; i++){
			
			double momentum = (double)i * 0.10;
			
			if(momentum < xbins_p_proton[ix] || xbins_p_proton[ix+1] <= momentum)
				continue;
			
			for(int j = 0; j < NEVENT; j++){
				
				if(j%10000==0) fprintf(stderr,"\r Bin:%3d  %10lld / %10lld", i, j, NEVENT);
				
				tree->GetEntry(j);
				
				tmp->Fill(fluxerr_p_proton_total[i]);
				
			}
			
		}
		
		tmp->Fit(f_gaus,"Q","",0.0,2.0);
		
		double center = f_gaus->GetParameter(1);
		double sigma  = f_gaus->GetParameter(2);
		
		if(tmp->GetMean()<0.01)
			h_fluxerr_p_proton_total->Fill((xbins_p_proton[ix]+xbins_p_proton[ix+1])/2.0, 0.0);
		else
			h_fluxerr_p_proton_total->Fill((xbins_p_proton[ix]+xbins_p_proton[ix+1])/2.0, sigma);
		
		tmp->Delete();
		
		//if(momentum > 5.0) break;
		
		//fprintf(stderr,"\r%10lld / %10lld\n", nthrows, nthrows);
		
	}
	cerr << endl;
	
	
	
	// Draw //
	
	// multiplicity //
	TCanvas *can_mul = new TCanvas("can_mul","",1600,1400);
	can_mul->cd();
	can_mul->SetTicks();
	TH1D *h_frame_mul = new TH1D("h_frame_mul","",10,-0.5,9.5);
	h_frame_mul->Draw();
	h_frame_mul->SetTitle("Multiplicity;multiplicity;Fractional error");
	h_frame_mul->GetXaxis()->SetRangeUser(0.5,9.5);
	h_frame_mul->SetMinimum(0.00);
	h_frame_mul->SetMaximum(0.20);
	h_frame_mul->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_mul_total->Draw("sames");
	h_fluxerr_mul_total->SetLineColor(kBlack);
	h_fluxerr_mul_total->SetLineWidth(2);
	
	// Number of pions //
	TCanvas *can_mul_pipm = new TCanvas("can_mul_pipm","",1600,1400);
	can_mul_pipm->cd();
	can_mul_pipm->SetTicks();
	TH1D *h_frame_mul_pipm = new TH1D("h_frame_mul_pipm","",10,-0.5,9.5);
	h_frame_mul_pipm->Draw();
	h_frame_mul_pipm->SetTitle("Number of charged pions;# of pions;Fractional error");
	h_frame_mul_pipm->GetXaxis()->SetRangeUser(-0.5,9.5);
	h_frame_mul_pipm->SetMinimum(0.00);
	h_frame_mul_pipm->SetMaximum(0.20);
	h_frame_mul_pipm->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_mul_pipm_total->Draw("sames");
	h_fluxerr_mul_pipm_total->SetLineColor(kBlack);
	h_fluxerr_mul_pipm_total->SetLineWidth(2);
	
	// Number of protons //
	TCanvas *can_mul_proton = new TCanvas("can_mul_proton","",1600,1400);
	can_mul_proton->cd();
	can_mul_proton->SetTicks();
	TH1D *h_frame_mul_proton = new TH1D("h_frame_mul_proton","",10,-0.5,9.5);
	h_frame_mul_proton->Draw();
	h_frame_mul_proton->SetTitle("Number of charged protons;# of protons;Fractional error");
	h_frame_mul_proton->GetXaxis()->SetRangeUser(-0.5,9.5);
	h_frame_mul_proton->SetMinimum(0.00);
	h_frame_mul_proton->SetMaximum(0.20);
	h_frame_mul_proton->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_mul_proton_total->Draw("sames");
	h_fluxerr_mul_proton_total->SetLineColor(kBlack);
	h_fluxerr_mul_proton_total->SetLineWidth(2);
	
	// Number of protons vs. pions //
	TCanvas *can_mul_proton_pipm = new TCanvas("can_mul_proton_pipm","",1600,1400);
	can_mul_proton_pipm->cd();
	can_mul_proton_pipm->SetTicks();
	TH2D *h_frame_mul_proton_pipm = new TH2D("h_frame_mul_proton_pipm","",10,-0.5,9.5,10,-0.5,9.5);
	h_frame_mul_proton_pipm->Draw();
	h_frame_mul_proton_pipm->SetTitle("Number of charged protons vs. protons;# of protons;# of pions;Fractional error");
	h_frame_mul_proton_pipm->GetXaxis()->SetTitleOffset(1.3);
	h_frame_mul_proton_pipm->GetYaxis()->SetTitleOffset(1.3);
	h_frame_mul_proton_pipm->GetXaxis()->SetRangeUser(-0.5,4.5);
	h_frame_mul_proton_pipm->GetYaxis()->SetRangeUser(-0.5,4.5);
	h_frame_mul_proton_pipm->SetMinimum(0.00);
	h_frame_mul_proton_pipm->SetMaximum(0.20);
	h_frame_mul_proton_pipm->GetXaxis()->SetTitleOffset(1.33);
	h_frame_mul_proton_pipm->GetYaxis()->SetTitleOffset(1.33);
	h_frame_mul_proton_pipm->GetZaxis()->SetTitleOffset(1.33);
	h_fluxerr_mul_proton_pipm_total->Draw("colz sames");
	h_fluxerr_mul_proton_pipm_total->SetLineColor(kBlack);
	h_fluxerr_mul_proton_pipm_total->SetLineWidth(2);
	
	// Angle of muons //
	TCanvas *can_angle_corr_deg_muon = new TCanvas("can_angle_corr_deg_muon","",1600,1400);
	can_angle_corr_deg_muon->cd();
	can_angle_corr_deg_muon->SetTicks();
	TH1D *h_frame_angle_corr_deg_muon = new TH1D("h_frame_angle_corr_deg_muon","",180,0.0,180.0);
	h_frame_angle_corr_deg_muon->Draw();
	h_frame_angle_corr_deg_muon->SetTitle("Emission angle of muons;#theta_{#mu} (deg.);Fractional error");
	h_frame_angle_corr_deg_muon->GetXaxis()->SetRangeUser(0.0,90.0);
	h_frame_angle_corr_deg_muon->SetMinimum(0.00);
	h_frame_angle_corr_deg_muon->SetMaximum(0.20);
	h_frame_angle_corr_deg_muon->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_angle_corr_deg_muon_total->Draw("sames");
	h_fluxerr_angle_corr_deg_muon_total->SetLineColor(kBlack);
	h_fluxerr_angle_corr_deg_muon_total->SetLineWidth(2);
	
	// Angle of pions //
	TCanvas *can_angle_corr_deg_pipm = new TCanvas("can_angle_corr_deg_pipm","",1600,1400);
	can_angle_corr_deg_pipm->cd();
	can_angle_corr_deg_pipm->SetTicks();
	TH1D *h_frame_angle_corr_deg_pipm = new TH1D("h_frame_angle_corr_deg_pipm","",180,0.0,180.0);
	h_frame_angle_corr_deg_pipm->Draw();
	h_frame_angle_corr_deg_pipm->SetTitle("Emission angle of charged pions;#theta_{#pi^{#pm}} (deg.);Fractional error");
	h_frame_angle_corr_deg_pipm->GetXaxis()->SetRangeUser(0.0,180.0);
	h_frame_angle_corr_deg_pipm->SetMinimum(0.00);
	h_frame_angle_corr_deg_pipm->SetMaximum(0.20);
	h_frame_angle_corr_deg_pipm->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_angle_corr_deg_pipm_total->Draw("sames");
	h_fluxerr_angle_corr_deg_pipm_total->SetLineColor(kBlack);
	h_fluxerr_angle_corr_deg_pipm_total->SetLineWidth(2);
	
	// Angle of protons //
	TCanvas *can_angle_corr_deg_proton = new TCanvas("can_angle_corr_deg_proton","",1600,1400);
	can_angle_corr_deg_proton->cd();
	can_angle_corr_deg_proton->SetTicks();
	TH1D *h_frame_angle_corr_deg_proton = new TH1D("h_frame_angle_corr_deg_proton","",180,0.0,180.0);
	h_frame_angle_corr_deg_proton->Draw();
	h_frame_angle_corr_deg_proton->SetTitle("Emission angle of protons;#theta_{p} (deg.);Fractional error");
	h_frame_angle_corr_deg_proton->GetXaxis()->SetRangeUser(0.0,180.0);
	h_frame_angle_corr_deg_proton->SetMinimum(0.00);
	h_frame_angle_corr_deg_proton->SetMaximum(0.20);
	h_frame_angle_corr_deg_proton->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_angle_corr_deg_proton_total->Draw("sames");
	h_fluxerr_angle_corr_deg_proton_total->SetLineColor(kBlack);
	h_fluxerr_angle_corr_deg_proton_total->SetLineWidth(2);
	
	// Angle of muons (cos) //
	TCanvas *can_angle_corr_cos_muon = new TCanvas("can_angle_corr_cos_muon","",1600,1400);
	can_angle_corr_cos_muon->cd();
	can_angle_corr_cos_muon->SetTicks();
	TH1D *h_frame_angle_corr_cos_muon = new TH1D("h_frame_angle_corr_cos_muon","",200,-1.0,1.0);
	h_frame_angle_corr_cos_muon->Draw();
	h_frame_angle_corr_cos_muon->SetTitle("Emission angle of muons;cos#theta_{#mu};Fractional error");
	h_frame_angle_corr_cos_muon->GetXaxis()->SetRangeUser(0.0,1.0);
	h_frame_angle_corr_cos_muon->SetMinimum(0.00);
	h_frame_angle_corr_cos_muon->SetMaximum(0.20);
	h_frame_angle_corr_cos_muon->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_angle_corr_cos_muon_total->Draw("sames");
	h_fluxerr_angle_corr_cos_muon_total->SetLineColor(kBlack);
	h_fluxerr_angle_corr_cos_muon_total->SetLineWidth(2);
	
	// Angle of pions (cos) //
	TCanvas *can_angle_corr_cos_pipm = new TCanvas("can_angle_corr_cos_pipm","",1600,1400);
	can_angle_corr_cos_pipm->cd();
	can_angle_corr_cos_pipm->SetTicks();
	TH1D *h_frame_angle_corr_cos_pipm = new TH1D("h_frame_angle_corr_cos_pipm","",200,-1.0,1.0);
	h_frame_angle_corr_cos_pipm->Draw();
	h_frame_angle_corr_cos_pipm->SetTitle("Emission angle of charged pions;cos#theta_{#pi^{#pm}};Fractional error");
	h_frame_angle_corr_cos_pipm->GetXaxis()->SetRangeUser(-1.0,1.0);
	h_frame_angle_corr_cos_pipm->SetMinimum(0.00);
	h_frame_angle_corr_cos_pipm->SetMaximum(0.20);
	h_frame_angle_corr_cos_pipm->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_angle_corr_cos_pipm_total->Draw("sames");
	h_fluxerr_angle_corr_cos_pipm_total->SetLineColor(kBlack);
	h_fluxerr_angle_corr_cos_pipm_total->SetLineWidth(2);
	
	// Angle of protons (cos) //
	TCanvas *can_angle_corr_cos_proton = new TCanvas("can_angle_corr_cos_proton","",1600,1400);
	can_angle_corr_cos_proton->cd();
	can_angle_corr_cos_proton->SetTicks();
	TH1D *h_frame_angle_corr_cos_proton = new TH1D("h_frame_angle_corr_cos_proton","",200,-1.0,1.0);
	h_frame_angle_corr_cos_proton->Draw();
	h_frame_angle_corr_cos_proton->SetTitle("Emission angle of protons;cos#theta_{p};Fractional error");
	h_frame_angle_corr_cos_proton->GetXaxis()->SetRangeUser(-1.0,1.0);
	h_frame_angle_corr_cos_proton->SetMinimum(0.00);
	h_frame_angle_corr_cos_proton->SetMaximum(0.20);
	h_frame_angle_corr_cos_proton->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_angle_corr_cos_proton_total->Draw("sames");
	h_fluxerr_angle_corr_cos_proton_total->SetLineColor(kBlack);
	h_fluxerr_angle_corr_cos_proton_total->SetLineWidth(2);
	
	// Momentum of muons //
	TCanvas *can_p_muon = new TCanvas("can_p_muon","",1600,1400);
	can_p_muon->cd();
	can_p_muon->SetTicks();
	TH1D *h_frame_p_muon = new TH1D("h_frame_p_muon","",500,0.0,50.0);
	h_frame_p_muon->Draw();
	h_frame_p_muon->SetTitle("Momentum of muons;P_{#mu} (GeV/c);Fractional error");
	h_frame_p_muon->GetXaxis()->SetRangeUser(0.0,6.0);
	h_frame_p_muon->SetMinimum(0.00);
	h_frame_p_muon->SetMaximum(0.20);
	h_frame_p_muon->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_p_muon_total->Draw("sames");
	h_fluxerr_p_muon_total->SetLineColor(kBlack);
	h_fluxerr_p_muon_total->SetLineWidth(2);
	
	// Momentum of pions //
	TCanvas *can_p_pipm = new TCanvas("can_p_pipm","",1600,1400);
	can_p_pipm->cd();
	can_p_pipm->SetTicks();
	TH1D *h_frame_p_pipm = new TH1D("h_frame_p_pipm","",500,0.0,50.0);
	h_frame_p_pipm->Draw();
	h_frame_p_pipm->SetTitle("Momentum of charged pions;P_{#pi^{#pm}} (GeV/c);Fractional error");
	h_frame_p_pipm->GetXaxis()->SetRangeUser(0.0,3.0);
	h_frame_p_pipm->SetMinimum(0.00);
	h_frame_p_pipm->SetMaximum(0.20);
	h_frame_p_pipm->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_p_pipm_total->Draw("sames");
	h_fluxerr_p_pipm_total->SetLineColor(kBlack);
	h_fluxerr_p_pipm_total->SetLineWidth(2);
	
	// Momentum of protons //
	TCanvas *can_p_proton = new TCanvas("can_p_proton","",1600,1400);
	can_p_proton->cd();
	can_p_proton->SetTicks();
	TH1D *h_frame_p_proton = new TH1D("h_frame_p_proton","",500,0.0,50.0);
	h_frame_p_proton->Draw();
	h_frame_p_proton->SetTitle("Momentum of protons;P_{p} (GeV/c);Fractional error");
	h_frame_p_proton->GetXaxis()->SetRangeUser(0.0,3.0);
	h_frame_p_proton->SetMinimum(0.00);
	h_frame_p_proton->SetMaximum(0.20);
	h_frame_p_proton->GetYaxis()->SetTitleOffset(1.33);
	h_fluxerr_p_proton_total->Draw("sames");
	h_fluxerr_p_proton_total->SetLineColor(kBlack);
	h_fluxerr_p_proton_total->SetLineWidth(2);
	
	
	
	
	
	
	
	//////////////////////
	//   CREATE PLOTS   //
	//////////////////////
	char pdfname[128];
	sprintf(pdfname, "%s", pdf_filename_o1);
	char pdfname_sta[128];
	sprintf(pdfname_sta,"%s(",pdfname);
	char pdfname_end[128];
	sprintf(pdfname_end,"%s)",pdfname);
	
	TCanvas *can_sta = new TCanvas("can_sta","",1600,1200);
	can_sta->Print(pdfname_sta,"pdf Portrait");
	
	
	can_mul->Print(pdfname,"pdf Portrait");
	can_mul_pipm->Print(pdfname,"pdf Portrait");
	can_mul_proton->Print(pdfname,"pdf Portrait");
	can_mul_proton_pipm->Print(pdfname,"pdf Portrait");
	can_angle_corr_deg_muon->Print(pdfname,"pdf Portrait");
	can_angle_corr_deg_pipm->Print(pdfname,"pdf Portrait");
	can_angle_corr_deg_proton->Print(pdfname,"pdf Portrait");
	can_angle_corr_cos_muon->Print(pdfname,"pdf Portrait");
	can_angle_corr_cos_pipm->Print(pdfname,"pdf Portrait");
	can_angle_corr_cos_proton->Print(pdfname,"pdf Portrait");
	can_p_muon->Print(pdfname,"pdf Portrait");
	can_p_pipm->Print(pdfname,"pdf Portrait");
	can_p_proton->Print(pdfname,"pdf Portrait");
	
	
	TCanvas *can_end = new TCanvas("can_end","",1600,1200);
	can_sta->Print(pdfname_end,"pdf Portrait");
	
	
	
	// Write //
	root_o1->cd();
	
	h_fluxerr_mul_total-> Write();
	h_fluxerr_mul_pipm_total-> Write();
	h_fluxerr_mul_proton_total-> Write();
	h_fluxerr_mul_proton_pipm_total-> Write();
	
	h_fluxerr_angle_corr_deg_muon_total-> Write();
	h_fluxerr_angle_corr_deg_pipm_total-> Write();
	h_fluxerr_angle_corr_deg_proton_total-> Write();
	
	h_fluxerr_angle_corr_cos_muon_total-> Write();
	h_fluxerr_angle_corr_cos_pipm_total-> Write();
	h_fluxerr_angle_corr_cos_proton_total-> Write();
	
	h_fluxerr_p_muon_total-> Write();
	h_fluxerr_p_pipm_total-> Write();
	h_fluxerr_p_proton_total-> Write();
	
	root_o1->Close();
	
	
	return 0;
	
	
}


