///////////////////////////////////////////////////////////////////////////////////////
// LOG ////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//  
//  Original program is from T. Kikawa.
//  Arranged by A. Hiramot.
//  
// 2020/07/15 H.Oshima
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


// Number of MC throws //
#define nthrows 100000

// Covariance matrix bin //
#define nbin 40


void draw_matrix(float mat[nbin][nbin]);
void cholcov_conv(float covmat[nbin][nbin], float cholcovmat[nbin][nbin]);


int main(int argc, char *argv[]){
	
	// neutrino beam mode //
	int beam = 0;
	bool flg_b = false;
	
	// neutrino component //
	int neutrino = 0;
	bool flg_n = false;
	
	//  input file name
	char root_filename_i1[1024];
	bool flg_i = false;
	
	//  input covariance file name
	char covariance_filename_i1[1024];
	bool flg_c = false;
	
	// output file name
	char root_filename_o1[1024];
	bool flg_o = false;
	
	// Module (ND) //
	int nd = 0;
	bool flg_m = false;
	
	// interaction target //
	char target_name[1024];
	bool flg_a = false;
	
	
	
	// option //
	int c = -1;
	
	while((c=getopt(argc, argv, "i:c:o:f:b:n:m:a:")) != -1){
		
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
			case 'c':
				sprintf(covariance_filename_i1,"%s",optarg);
				flg_c = true;
				break;
			case 'm':
				nd = atoi(optarg);
				flg_m = true;
				break;
			case 'a':
				sprintf(target_name,"%s",optarg);
				flg_a = true;
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
	else{
		cerr << "Error : -n option must be 1/-1." << endl;
		cerr << "1 ... Neutrino component,  -1 ... Anti-neutrino component" << endl;
		exit(1);
	}
	
	
	//  TFile *fop = TFile::Open("flux_covariance.root");
	//  TH2D *hcov = (TH2D*)fop->Get("hcov");
	//  if(hcov==NULL){
	//    cout << "Cannot get matrix from " << fop->GetPath() << endl;
	//    exit;
	//  }
	
	// Flux covariance matrix file //
	if(!flg_c){
		cerr << "Error : -c option is need for flux covariance root file." << endl;
		exit(1);
	}
	
	TFile *fop = TFile::Open(covariance_filename_i1,"read");
	if (fop->IsZombie()){
		cerr << "Failed to open " << covariance_filename_i1 << "." << endl;
		return 0;
	}
	
	//TMatrixDSym *mat = (TMatrixDSym*)fop->Get("flux_cov");
	TMatrixDSym *mat = (TMatrixDSym*)fop->Get("total_flux_cov");
	if(mat==NULL){
		cout << "Cannot get matrix from " << fop->GetPath() << endl;
		exit(1);
	}
	
	
	
	// Flux Covariance matrix Nbin //
	int bin_offset_beammode  = 0;
	
	// FHC //
	if(beam==1){
		bin_offset_beammode = 0;
	}
	// RHC //
	else if(beam==-1){
		bin_offset_beammode = 40;
	}
	
	
	float cov_mat[nbin][nbin];
	
	for(int i=0; i<nbin; i++){
		for(int j=0; j<nbin; j++){
			cov_mat[i][j] = (*mat)(i+bin_offset_beammode,j+bin_offset_beammode);
			//cout<<i<<" "<<j<<" "<<cov_mat[i][j]<<endl;
		}
	}
	
	
	TCanvas *can_dummy = new TCanvas("can_dummy","",1600,1400);
	
	#ifdef MKPLOT
		draw_matrix(cov_mat);
	#endif
	
	float chol_mat[nbin][nbin];
	cholcov_conv(cov_mat, chol_mat);
	
	#ifdef MKPLOT2
		draw_matrix(chol_mat);
	#endif
	
	
	
	
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
	
	
	TFile *root_i1 = new TFile(root_filename_i1,"read");  // Open the Geant output ROOT file
	if (root_i1->IsZombie()){
		cerr << "Failed to open " << root_filename_i1 << "." << endl;
		return 0;
	}
	
	float weight[nbin],nrand[nbin];
	
	
	TH2D *h_stack_Enu_mul_total                        =   (TH2D*)root_i1->Get(Form("h_stack_Enu_mul_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_mul_pipm_total                   =   (TH2D*)root_i1->Get(Form("h_stack_Enu_mul_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_mul_proton_total                 =   (TH2D*)root_i1->Get(Form("h_stack_Enu_mul_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_angle_corr_deg_muon_total        =   (TH2D*)root_i1->Get(Form("h_stack_Enu_angle_corr_deg_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_angle_corr_deg_pipm_total        =   (TH2D*)root_i1->Get(Form("h_stack_Enu_angle_corr_deg_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_angle_corr_deg_proton_total      =   (TH2D*)root_i1->Get(Form("h_stack_Enu_angle_corr_deg_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_angle_corr_cos_muon_total        =   (TH2D*)root_i1->Get(Form("h_stack_Enu_angle_corr_cos_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_angle_corr_cos_pipm_total        =   (TH2D*)root_i1->Get(Form("h_stack_Enu_angle_corr_cos_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_angle_corr_cos_proton_total      =   (TH2D*)root_i1->Get(Form("h_stack_Enu_angle_corr_cos_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_p_muon_total                     =   (TH2D*)root_i1->Get(Form("h_stack_Enu_p_muon_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_p_proton_total                   =   (TH2D*)root_i1->Get(Form("h_stack_Enu_p_proton_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_p_pipm_total                     =   (TH2D*)root_i1->Get(Form("h_stack_Enu_p_pipm_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	
	TH2D *h_stack_Enu_opang_cos_pp_0pi2p_total         =   (TH2D*)root_i1->Get(Form("h_stack_Enu_opang_cos_pp_0pi2p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_total_p_pp_0pi2p_total           =   (TH2D*)root_i1->Get(Form("h_stack_Enu_total_p_pp_0pi2p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_inv_mass_pp_0pi2p_total          =   (TH2D*)root_i1->Get(Form("h_stack_Enu_inv_mass_pp_0pi2p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_p_low_p_0pi2p_total              =   (TH2D*)root_i1->Get(Form("h_stack_Enu_p_low_p_0pi2p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_p_high_p_0pi2p_total             =   (TH2D*)root_i1->Get(Form("h_stack_Enu_p_high_p_0pi2p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_inv_mass_pip_1pi1p_total         =   (TH2D*)root_i1->Get(Form("h_stack_Enu_inv_mass_pip_1pi1p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_STV_dp_0pi1p_total               =   (TH2D*)root_i1->Get(Form("h_stack_Enu_STV_dp_0pi1p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_STV_dalpha_0pi1p_total           =   (TH2D*)root_i1->Get(Form("h_stack_Enu_STV_dalpha_0pi1p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_STV_dphi_0pi1p_total             =   (TH2D*)root_i1->Get(Form("h_stack_Enu_STV_dphi_0pi1p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_IPK_dp_scalar_sum_0pi1p_total    =   (TH2D*)root_i1->Get(Form("h_stack_Enu_IPK_dp_scalar_sum_0pi1p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_IPK_dtheta_0pi1p_total           =   (TH2D*)root_i1->Get(Form("h_stack_Enu_IPK_dtheta_0pi1p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	TH2D *h_stack_Enu_IPK_dp_vector_sum_0pi1p_total    =   (TH2D*)root_i1->Get(Form("h_stack_Enu_IPK_dp_vector_sum_0pi1p_total_%s_nd%d_%s", name_neutrino, nd, target_name));
	
	
	const int nbin_mul_pipm_ppi   = 5;
	const int nbin_mul_proton_ppi = 5;
	TH2D *h_stack_Enu_mul_proton_Npipm_total[nbin_mul_pipm_ppi];
	
	for(int i = 0; i < nbin_mul_pipm_ppi; i++){
		h_stack_Enu_mul_proton_Npipm_total[i]    = (TH2D*)root_i1->Get(Form("h_stack_Enu_mul_proton_pipm%d_total_%s_nd%d_%s", i, name_neutrino, nd, target_name));
	}
	
	
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
	
	TTree *tree = new TTree("tree","tree");
	
	const int nbin_mul              =  10;
	const int nbin_mul_pipm         =  10;
	const int nbin_mul_proton       =  10;
	const int nbin_angle_muon       =  36; //   5 deg.  bin
	const int nbin_angle_pipm       =  18; //  10 deg.  bin
	const int nbin_angle_proton     =  18; //  10 deg.  bin
	const int nbin_angle_cos_muon   = 100; // 0.02      bin
	const int nbin_angle_cos_pipm   =  20; // 0.10      bin
	const int nbin_angle_cos_proton =  20; // 0.10      bin
	const int nbin_p_muon           = 250; // 0.2 GeV/c bin
	const int nbin_p_pipm           = 500; // 0.1 GeV/c bin
	const int nbin_p_proton         = 500; // 0.1 GeV/c bin
	
	const int nbin_opang_cos_pp_0pi2p           =   10; // 0.2     /bin
	const int nbin_total_p_pp_0pi2p             =  125; // 0.4     /bin
	const int nbin_inv_mass_pp_0pi2p            = 1250; // 0.04    /bin
	const int nbin_p_low_p_0pi2p                =  500; // 0.1     /bin
	const int nbin_p_high_p_0pi2p               =  500; // 0.1     /bin
	const int nbin_low_angle_corr_deg_p_0pi2p   =   18; // 10 deg. /bin
	const int nbin_high_angle_corr_deg_p_0pi2p  =   18; // 10 deg. /bin
	const int nbin_inv_mass_pip_1pi1p           =  500; // 0.1     /bin
	const int nbin_STV_dp_0pi1p                 =  500; // 0.2     /bin
	const int nbin_STV_dalpha_0pi1p             =   36; // 10 deg. /bin
	const int nbin_STV_dphi_0pi1p               =   36; // 10 deg. /bin
	const int nbin_IPK_dp_scalar_sum_0pi1p      =  500; // 0.2     /bin
	const int nbin_IPK_dtheta_0pi1p             =   36; // 10 deg. /bin
	const int nbin_IPK_dp_vector_sum_0pi1p      =  500; // 0.2     /bin
	
	const int rebin_inv_mass_pp_0pi2p = 1;
	
	// [0] ... nominal, [1] ... weighted
	double mul_total[nbin_mul][2];
	double mul_pipm_total[nbin_mul_pipm][2];
	double mul_proton_total[nbin_mul_proton][2];
	double angle_corr_deg_muon_total[nbin_angle_muon][2];
	double angle_corr_deg_pipm_total[nbin_angle_pipm][2];
	double angle_corr_deg_proton_total[nbin_angle_proton][2];
	double angle_corr_cos_muon_total[nbin_angle_cos_muon][2];
	double angle_corr_cos_pipm_total[nbin_angle_cos_pipm][2];
	double angle_corr_cos_proton_total[nbin_angle_cos_proton][2];
	double p_muon_total[nbin_p_muon][2];
	double p_pipm_total[nbin_p_pipm][2];
	double p_proton_total[nbin_p_proton][2];
	
	double mul_proton_Npipm_total[nbin_mul_pipm_ppi][nbin_mul_proton_ppi][2];
	
	double opang_cos_pp_0pi2p_total[nbin_opang_cos_pp_0pi2p][2];
	double total_p_pp_0pi2p_total[nbin_total_p_pp_0pi2p][2];
	double inv_mass_pp_0pi2p_total[nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p][2];
	double p_low_p_0pi2p_total[nbin_p_low_p_0pi2p][2];
	double p_high_p_0pi2p_total[nbin_p_high_p_0pi2p][2];
	double inv_mass_pip_1pi1p_total[nbin_inv_mass_pip_1pi1p][2];
	double STV_dp_0pi1p_total[nbin_STV_dp_0pi1p][2];
	double STV_dalpha_0pi1p_total[nbin_STV_dalpha_0pi1p][2];
	double STV_dphi_0pi1p_total[nbin_STV_dphi_0pi1p][2];
	double IPK_dp_scalar_sum_0pi1p_total[nbin_IPK_dp_scalar_sum_0pi1p][2];
	double IPK_dtheta_0pi1p_total[nbin_IPK_dtheta_0pi1p][2];
	double IPK_dp_vector_sum_0pi1p_total[nbin_IPK_dp_vector_sum_0pi1p][2];
	
	// initialize //
	for(int j = 0; j < nbin_mul; j++){
		mul_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_mul_pipm; j++){
		mul_pipm_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_mul_proton; j++){
		mul_proton_total[j][0] = 0;
	}
	for(int k = 0; k < nbin_mul_pipm_ppi; k++){
		for(int j = 0; j < nbin_mul_proton_ppi; j++){
			mul_proton_Npipm_total[k][j][0] = 0;
		}
	}
	for(int j = 0; j < nbin_angle_muon; j++){
		angle_corr_deg_muon_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_angle_pipm; j++){
		angle_corr_deg_pipm_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_angle_proton; j++){
		angle_corr_deg_proton_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_angle_cos_muon; j++){
		angle_corr_cos_muon_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_angle_cos_pipm; j++){
		angle_corr_cos_pipm_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_angle_cos_proton; j++){
		angle_corr_cos_proton_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_p_muon; j++){
		p_muon_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_p_pipm; j++){
		p_pipm_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_p_proton; j++){
		p_muon_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_opang_cos_pp_0pi2p; j++){
		opang_cos_pp_0pi2p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_total_p_pp_0pi2p; j++){
		total_p_pp_0pi2p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_p_low_p_0pi2p; j++){
		p_low_p_0pi2p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_p_high_p_0pi2p; j++){
		p_high_p_0pi2p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p; j++){
		inv_mass_pp_0pi2p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_inv_mass_pip_1pi1p; j++){
		inv_mass_pip_1pi1p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_STV_dp_0pi1p; j++){
		STV_dp_0pi1p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_STV_dalpha_0pi1p; j++){
		STV_dalpha_0pi1p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_STV_dphi_0pi1p; j++){
		STV_dphi_0pi1p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_IPK_dp_scalar_sum_0pi1p; j++){
		IPK_dp_scalar_sum_0pi1p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_IPK_dtheta_0pi1p; j++){
		IPK_dtheta_0pi1p_total[j][0] = 0;
	}
	for(int j = 0; j < nbin_IPK_dp_vector_sum_0pi1p; j++){
		IPK_dp_vector_sum_0pi1p_total[j][0] = 0;
	}
	
	
	// Error //
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
	
	double fluxerr_mul_proton_Npipm_total[nbin_mul_pipm_ppi][nbin_mul_proton_ppi]    = {};
	
	double fluxerr_opang_cos_pp_0pi2p_total[nbin_opang_cos_pp_0pi2p]           = {0.0};
	double fluxerr_total_p_pp_0pi2p_total[nbin_total_p_pp_0pi2p]               = {0.0};
	double fluxerr_inv_mass_pp_0pi2p_total[nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p]             = {0.0};
	double fluxerr_p_low_p_0pi2p_total[nbin_p_low_p_0pi2p]                     = {0.0};
	double fluxerr_p_high_p_0pi2p_total[nbin_p_high_p_0pi2p]                   = {0.0};
	double fluxerr_inv_mass_pip_1pi1p_total[nbin_inv_mass_pip_1pi1p]           = {0.0};
	double fluxerr_STV_dp_0pi1p_total[nbin_STV_dp_0pi1p]                       = {0.0};
	double fluxerr_STV_dalpha_0pi1p_total[nbin_STV_dalpha_0pi1p]               = {0.0};
	double fluxerr_STV_dphi_0pi1p_total[nbin_STV_dphi_0pi1p]                   = {0.0};
	double fluxerr_IPK_dp_scalar_sum_0pi1p_total[nbin_IPK_dp_scalar_sum_0pi1p] = {0.0};
	double fluxerr_IPK_dtheta_0pi1p_total[nbin_IPK_dtheta_0pi1p]               = {0.0};
	double fluxerr_IPK_dp_vector_sum_0pi1p_total[nbin_IPK_dp_vector_sum_0pi1p] = {0.0};
	
	// Tree Branch //
	tree->Branch(Form("fluxerr_mul_%s_total", name_neutrino),                       fluxerr_mul_total,                       Form("fluxerr_mul_%s_total[%d]/D",name_neutrino,nbin_mul));
	tree->Branch(Form("fluxerr_mul_pipm_%s_total", name_neutrino),                  fluxerr_mul_pipm_total,                  Form("fluxerr_mul_pipm_%s_total[%d]/D",name_neutrino,nbin_mul_pipm));
	tree->Branch(Form("fluxerr_mul_proton_%s_total", name_neutrino),                fluxerr_mul_proton_total,                Form("fluxerr_mul_proton_%s_total[%d]/D",name_neutrino,nbin_mul_proton));
	tree->Branch(Form("fluxerr_angle_corr_deg_muon_%s_total", name_neutrino),       fluxerr_angle_corr_deg_muon_total,       Form("fluxerr_angle_corr_deg_muon_%s_total[%d]/D",name_neutrino,nbin_angle_muon));
	tree->Branch(Form("fluxerr_angle_corr_deg_pipm_%s_total", name_neutrino),       fluxerr_angle_corr_deg_pipm_total,       Form("fluxerr_angle_corr_deg_pipm_%s_total[%d]/D",name_neutrino,nbin_angle_pipm));
	tree->Branch(Form("fluxerr_angle_corr_deg_proton_%s_total", name_neutrino),     fluxerr_angle_corr_deg_proton_total,     Form("fluxerr_angle_corr_deg_proton_%s_total[%d]/D",name_neutrino,nbin_angle_proton));
	tree->Branch(Form("fluxerr_angle_corr_cos_muon_%s_total", name_neutrino),       fluxerr_angle_corr_cos_muon_total,       Form("fluxerr_angle_corr_cos_muon_%s_total[%d]/D",name_neutrino,nbin_angle_cos_muon));
	tree->Branch(Form("fluxerr_angle_corr_cos_pipm_%s_total", name_neutrino),       fluxerr_angle_corr_cos_pipm_total,       Form("fluxerr_angle_corr_cos_pipm_%s_total[%d]/D",name_neutrino,nbin_angle_cos_pipm));
	tree->Branch(Form("fluxerr_angle_corr_cos_proton_%s_total", name_neutrino),     fluxerr_angle_corr_cos_proton_total,     Form("fluxerr_angle_corr_cos_proton_%s_total[%d]/D",name_neutrino,nbin_angle_cos_proton));
	tree->Branch(Form("fluxerr_p_muon_%s_total", name_neutrino),                    fluxerr_p_muon_total,                    Form("fluxerr_p_muon_%s_total[%d]/D",name_neutrino,nbin_p_muon));
	tree->Branch(Form("fluxerr_p_pipm_%s_total", name_neutrino),                    fluxerr_p_pipm_total,                    Form("fluxerr_p_pipm_%s_total[%d]/D",name_neutrino,nbin_p_pipm));
	tree->Branch(Form("fluxerr_p_proton_%s_total", name_neutrino),                  fluxerr_p_proton_total,                  Form("fluxerr_p_proton_%s_total[%d]/D",name_neutrino,nbin_p_proton));
	tree->Branch(Form("fluxerr_mul_proton_Npipm_%s_total", name_neutrino),          fluxerr_mul_proton_Npipm_total,          Form("fluxerr_mul_proton_Npipm_%s_total[%d][%d]/D",name_neutrino,nbin_mul_pipm_ppi,nbin_mul_proton_ppi));
	
	tree->Branch(Form("fluxerr_opang_cos_pp_0pi2p_%s_total", name_neutrino),        fluxerr_opang_cos_pp_0pi2p_total,        Form("fluxerr_opang_cos_pp_0pi2p_%s_total[%d]/D",name_neutrino,nbin_opang_cos_pp_0pi2p));
	tree->Branch(Form("fluxerr_total_p_pp_0pi2p_%s_total", name_neutrino),          fluxerr_total_p_pp_0pi2p_total,          Form("fluxerr_total_p_pp_0pi2p_%s_total[%d]/D",name_neutrino,nbin_total_p_pp_0pi2p));
	tree->Branch(Form("fluxerr_inv_mass_pp_0pi2p_%s_total", name_neutrino),         fluxerr_inv_mass_pp_0pi2p_total,         Form("fluxerr_inv_mass_pp_0pi2p_%s_total[%d]/D",name_neutrino,nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p));
	tree->Branch(Form("fluxerr_p_low_p_0pi2p_%s_total", name_neutrino),             fluxerr_p_low_p_0pi2p_total,             Form("fluxerr_p_low_p_0pi2p_%s_total[%d]/D",name_neutrino,nbin_p_low_p_0pi2p));
	tree->Branch(Form("fluxerr_p_high_p_0pi2p_%s_total", name_neutrino),            fluxerr_p_high_p_0pi2p_total,            Form("fluxerr_p_high_p_0pi2p_%s_total[%d]/D",name_neutrino,nbin_p_high_p_0pi2p));
	tree->Branch(Form("fluxerr_inv_mass_pip_1pi1p_%s_total", name_neutrino),        fluxerr_inv_mass_pip_1pi1p_total,        Form("fluxerr_inv_mass_pip_1pi1p_%s_total[%d]/D",name_neutrino,nbin_inv_mass_pip_1pi1p));
	tree->Branch(Form("fluxerr_STV_dp_0pi1p_%s_total", name_neutrino),              fluxerr_STV_dp_0pi1p_total,              Form("fluxerr_STV_dp_0pi1p_%s_total[%d]/D",name_neutrino,nbin_STV_dp_0pi1p));
	tree->Branch(Form("fluxerr_STV_dalpha_0pi1p_%s_total", name_neutrino),          fluxerr_STV_dalpha_0pi1p_total,          Form("fluxerr_STV_dalpha_0pi1p_%s_total[%d]/D",name_neutrino,nbin_STV_dalpha_0pi1p));
	tree->Branch(Form("fluxerr_STV_dphi_0pi1p_%s_total", name_neutrino),            fluxerr_STV_dphi_0pi1p_total,            Form("fluxerr_STV_dphi_0pi1p_%s_total[%d]/D",name_neutrino,nbin_STV_dphi_0pi1p));
	tree->Branch(Form("fluxerr_IPK_dp_scalar_sum_%s_total", name_neutrino),         fluxerr_IPK_dp_scalar_sum_0pi1p_total,   Form("fluxerr_IPK_dp_scalar_sum_0pi1p_%s_total[%d]/D",name_neutrino,nbin_IPK_dp_scalar_sum_0pi1p));
	tree->Branch(Form("fluxerr_IPK_dtheta_0pi1p_%s_total", name_neutrino),          fluxerr_IPK_dtheta_0pi1p_total,          Form("fluxerr_IPK_dtheta_0pi1p_%s_total[%d]/D",name_neutrino,nbin_IPK_dtheta_0pi1p));
	tree->Branch(Form("fluxerr_IPK_dp_vector_sum_%s_total", name_neutrino),         fluxerr_IPK_dp_vector_sum_0pi1p_total,   Form("fluxerr_IPK_dp_vector_sum_0pi1p_%s_total[%d]/D",name_neutrino,nbin_IPK_dp_vector_sum_0pi1p));
	
	
	for(int i = 0; i < nbin/2; i++){
		
		// multiplicity //
		for(int j = 0; j < nbin_mul; j++){
			mul_total[j][0] += h_stack_Enu_mul_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_mul_pipm; j++){
			mul_pipm_total[j][0] += h_stack_Enu_mul_pipm_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_mul_proton; j++){
			mul_proton_total[j][0] += h_stack_Enu_mul_proton_total->GetBinContent(i+1, j+1);
		}
		
		for(int k = 0; k < nbin_mul_pipm_ppi; k++){
			for(int j = 0; j < nbin_mul_proton_ppi; j++){
				mul_proton_Npipm_total[k][j][0] += h_stack_Enu_mul_proton_Npipm_total[k]->GetBinContent(i+1, j+1);
			}
		}
		
		
		
		// angle //
		for(int j = 0; j < nbin_angle_muon; j++){
			angle_corr_deg_muon_total[j][0] += h_stack_Enu_angle_corr_deg_muon_total->GetBinContent(i+1, j+1);
		}
		
		for(int j = 0; j < nbin_angle_pipm; j++){
			angle_corr_deg_pipm_total[j][0] += h_stack_Enu_angle_corr_deg_pipm_total->GetBinContent(i+1, j+1);
		}
		
		for(int j = 0; j < nbin_angle_proton; j++){
			angle_corr_deg_proton_total[j][0] += h_stack_Enu_angle_corr_deg_proton_total->GetBinContent(i+1, j+1);
		}
		
		
		
		// angle (cos) //
		for(int j = 0; j < nbin_angle_cos_muon; j++){
			angle_corr_cos_muon_total[j][0] += h_stack_Enu_angle_corr_cos_muon_total->GetBinContent(i+1, j+1);
		}
		
		for(int j = 0; j < nbin_angle_cos_pipm; j++){
			angle_corr_cos_pipm_total[j][0] += h_stack_Enu_angle_corr_cos_pipm_total->GetBinContent(i+1, j+1);
		}
		
		for(int j = 0; j < nbin_angle_cos_proton; j++){
			angle_corr_cos_proton_total[j][0] += h_stack_Enu_angle_corr_cos_proton_total->GetBinContent(i+1, j+1);
		}
		
		
		
		// momentum //
		for(int j = 0; j < nbin_p_muon; j++){
			p_muon_total[j][0] += h_stack_Enu_p_muon_total->GetBinContent(i+1, j+1);
		}
		
		for(int j = 0; j < nbin_p_pipm; j++){
			p_pipm_total[j][0] += h_stack_Enu_p_pipm_total->GetBinContent(i+1, j+1);
		}
		
		for(int j = 0; j < nbin_p_proton; j++){
			p_proton_total[j][0] += h_stack_Enu_p_proton_total->GetBinContent(i+1, j+1);
		}
		
		
		
		// CC 0pi 2p //
		for(int j = 0; j < nbin_opang_cos_pp_0pi2p; j++){
			opang_cos_pp_0pi2p_total[j][0] += h_stack_Enu_opang_cos_pp_0pi2p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_total_p_pp_0pi2p; j++){
			total_p_pp_0pi2p_total[j][0] += h_stack_Enu_total_p_pp_0pi2p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p; j++){
			inv_mass_pp_0pi2p_total[j][0] += h_stack_Enu_inv_mass_pp_0pi2p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_p_low_p_0pi2p; j++){
			p_low_p_0pi2p_total[j][0] += h_stack_Enu_p_low_p_0pi2p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_p_high_p_0pi2p; j++){
			p_high_p_0pi2p_total[j][0] += h_stack_Enu_p_high_p_0pi2p_total->GetBinContent(i+1, j+1);
		}
		
		
		
		// CC 1pi 1p //
		for(int j = 0; j < nbin_inv_mass_pip_1pi1p; j++){
			inv_mass_pip_1pi1p_total[j][0] += h_stack_Enu_inv_mass_pip_1pi1p_total->GetBinContent(i+1, j+1);
		}
		
		
		
		// CC 0pi 1p //
		for(int j = 0; j < nbin_STV_dp_0pi1p; j++){
			STV_dp_0pi1p_total[j][0] += h_stack_Enu_STV_dp_0pi1p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_STV_dalpha_0pi1p; j++){
			STV_dalpha_0pi1p_total[j][0] += h_stack_Enu_STV_dalpha_0pi1p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_STV_dphi_0pi1p; j++){
			STV_dphi_0pi1p_total[j][0] += h_stack_Enu_STV_dphi_0pi1p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_IPK_dp_scalar_sum_0pi1p; j++){
			IPK_dp_scalar_sum_0pi1p_total[j][0] += h_stack_Enu_IPK_dp_scalar_sum_0pi1p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_IPK_dtheta_0pi1p; j++){
			IPK_dtheta_0pi1p_total[j][0] += h_stack_Enu_IPK_dtheta_0pi1p_total->GetBinContent(i+1, j+1);
		}
		for(int j = 0; j < nbin_IPK_dp_vector_sum_0pi1p; j++){
			IPK_dp_vector_sum_0pi1p_total[j][0] += h_stack_Enu_IPK_dp_vector_sum_0pi1p_total->GetBinContent(i+1, j+1);
		}
		
	}
	
	
	for(int n = 0; n < nthrows; n++){
		
		if(n%1000==0) fprintf(stderr,"\r%10lld / %10lld", n, nthrows);
		
		// initialize //
		// multiplicity //
		for(int j = 0; j < nbin_mul; j++){
			mul_total[j][1] = 0.0;
		}
		for(int j = 0; j < nbin_mul_pipm; j++){
			mul_pipm_total[j][1] = 0.0;
		}
		for(int j = 0; j < nbin_mul_proton; j++){
			mul_proton_total[j][1] = 0.0;
		}
		
		for(int k = 0; k < nbin_mul_pipm_ppi; k++){
			for(int j = 0; j < nbin_mul_proton_ppi; j++){
				mul_proton_Npipm_total[k][j][1] = 0.0;
			}
		}
		
		// angle //
		for(int j = 0; j < nbin_angle_muon; j++){
			angle_corr_deg_muon_total[j][1] = 0.0;
		}
		
		for(int j = 0; j < nbin_angle_pipm; j++){
			angle_corr_deg_pipm_total[j][1] = 0.0;
		}
		
		for(int j = 0; j < nbin_angle_proton; j++){
			angle_corr_deg_proton_total[j][1] = 0.0;
		}
		
		// angle (cos) //
		for(int j = 0; j < nbin_angle_cos_muon; j++){
			angle_corr_cos_muon_total[j][1] = 0.0;
		}
		
		for(int j = 0; j < nbin_angle_cos_pipm; j++){
			angle_corr_cos_pipm_total[j][1] = 0.0;
		}
		
		for(int j = 0; j < nbin_angle_cos_proton; j++){
			angle_corr_cos_proton_total[j][1] = 0.0;
		}
		
		// momentum //
		for(int j = 0; j < nbin_p_muon; j++){
			p_muon_total[j][1] = 0.0;
		}
		
		for(int j = 0; j < nbin_p_pipm; j++){
			p_pipm_total[j][1] = 0.0;
		}
		
		for(int j = 0; j < nbin_p_proton; j++){
			p_proton_total[j][1] = 0.0;
		}
		
		// CC 0pi 2p //
		for(int j = 0; j < nbin_opang_cos_pp_0pi2p; j++){
			opang_cos_pp_0pi2p_total[j][1] = 0;
		}
		
		for(int j = 0; j < nbin_total_p_pp_0pi2p; j++){
			total_p_pp_0pi2p_total[j][1] = 0;
		}
		
		for(int j = 0; j < nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p; j++){
			inv_mass_pp_0pi2p_total[j][1] = 0;
		}
		for(int j = 0; j < nbin_p_low_p_0pi2p; j++){
			p_low_p_0pi2p_total[j][1] = 0;
		}
		for(int j = 0; j < nbin_p_high_p_0pi2p; j++){
			p_high_p_0pi2p_total[j][1] += 0;
		}
		
		// CC 1pi 1p //
		for(int j = 0; j < nbin_inv_mass_pip_1pi1p; j++){
			inv_mass_pip_1pi1p_total[j][1] = 0;
		}
		
		// CC 0pi 1p //
		for(int j = 0; j < nbin_STV_dp_0pi1p; j++){
			STV_dp_0pi1p_total[j][1] = 0;
		}
		
		for(int j = 0; j < nbin_STV_dalpha_0pi1p; j++){
			STV_dalpha_0pi1p_total[j][1] = 0;
		}
		
		for(int j = 0; j < nbin_STV_dphi_0pi1p; j++){
			STV_dphi_0pi1p_total[j][1] = 0;
		}
		
		for(int j = 0; j < nbin_IPK_dp_scalar_sum_0pi1p; j++){
			IPK_dp_scalar_sum_0pi1p_total[j][1] = 0;
		}
		
		for(int j = 0; j < nbin_IPK_dtheta_0pi1p; j++){
			IPK_dtheta_0pi1p_total[j][1] = 0;
		}
		
		for(int j = 0; j < nbin_IPK_dp_vector_sum_0pi1p; j++){
			IPK_dp_vector_sum_0pi1p_total[j][1] = 0;
		}
		
		
		
		// gauss random number //
		for(int i = 0; i < nbin; i++){
			nrand[i] = gRandom->Gaus();
		}
		
		for(int i = 0; i < nbin; i++){
			
			weight[i] = 1.0;
			
			// weight //
			for(int k = 0; k <= i; k++){
				weight[i] += chol_mat[i][k] * nrand[k];
			}
			
			// integral //
			if(neutrino==1){
				
				if(i < nbin/2){
					
					// multiplicity //
					for(int j = 0; j < nbin_mul; j++){
						mul_total[j][1] += h_stack_Enu_mul_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_mul_pipm; j++){
						mul_pipm_total[j][1] += h_stack_Enu_mul_pipm_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_mul_proton; j++){
						mul_proton_total[j][1] += h_stack_Enu_mul_proton_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					for(int k = 0; k < nbin_mul_pipm_ppi; k++){
						for(int j = 0; j < nbin_mul_proton_ppi; j++){
							mul_proton_Npipm_total[k][j][1] += h_stack_Enu_mul_proton_Npipm_total[k]->GetBinContent(i+1, j+1)    *  weight[i];
						}
					}
					
					
					
					// angle //
					for(int j = 0; j < nbin_angle_muon; j++){
						angle_corr_deg_muon_total[j][1] += h_stack_Enu_angle_corr_deg_muon_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_angle_pipm; j++){
						angle_corr_deg_pipm_total[j][1] += h_stack_Enu_angle_corr_deg_pipm_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_angle_proton; j++){
						angle_corr_deg_proton_total[j][1] += h_stack_Enu_angle_corr_deg_proton_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					
					
					// angle (cos) //
					for(int j = 0; j < nbin_angle_cos_muon; j++){
						angle_corr_cos_muon_total[j][1] += h_stack_Enu_angle_corr_cos_muon_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_angle_cos_pipm; j++){
						angle_corr_cos_pipm_total[j][1] += h_stack_Enu_angle_corr_cos_pipm_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_angle_cos_proton; j++){
						angle_corr_cos_proton_total[j][1] += h_stack_Enu_angle_corr_cos_proton_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					
					
					// momentum //
					for(int j = 0; j < nbin_p_muon; j++){
						p_muon_total[j][1] += h_stack_Enu_p_muon_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_p_pipm; j++){
						p_pipm_total[j][1] += h_stack_Enu_p_pipm_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_p_proton; j++){
						p_proton_total[j][1] += h_stack_Enu_p_proton_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					
					
					// CC 0pi 2p //
					for(int j = 0; j < nbin_opang_cos_pp_0pi2p; j++){
						opang_cos_pp_0pi2p_total[j][1] += h_stack_Enu_opang_cos_pp_0pi2p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_total_p_pp_0pi2p; j++){
						total_p_pp_0pi2p_total[j][1] += h_stack_Enu_total_p_pp_0pi2p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p; j++){
						inv_mass_pp_0pi2p_total[j][1] += h_stack_Enu_inv_mass_pp_0pi2p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_p_low_p_0pi2p; j++){
						p_low_p_0pi2p_total[j][1] += h_stack_Enu_p_low_p_0pi2p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_p_high_p_0pi2p; j++){
						p_high_p_0pi2p_total[j][1] += h_stack_Enu_p_high_p_0pi2p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					
					
					// CC 1pi 1p //
					for(int j = 0; j < nbin_inv_mass_pip_1pi1p; j++){
						inv_mass_pip_1pi1p_total[j][1] += h_stack_Enu_inv_mass_pip_1pi1p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
					
					
					// CC 0pi 1p //
					for(int j = 0; j < nbin_STV_dp_0pi1p; j++){
						STV_dp_0pi1p_total[j][1] += h_stack_Enu_STV_dp_0pi1p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_STV_dalpha_0pi1p; j++){
						STV_dalpha_0pi1p_total[j][1] += h_stack_Enu_STV_dalpha_0pi1p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_STV_dphi_0pi1p; j++){
						STV_dphi_0pi1p_total[j][1] += h_stack_Enu_STV_dphi_0pi1p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_IPK_dp_scalar_sum_0pi1p; j++){
						IPK_dp_scalar_sum_0pi1p_total[j][1] += h_stack_Enu_IPK_dp_scalar_sum_0pi1p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_IPK_dtheta_0pi1p; j++){
						IPK_dtheta_0pi1p_total[j][1] += h_stack_Enu_IPK_dtheta_0pi1p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_IPK_dp_vector_sum_0pi1p; j++){
						IPK_dp_vector_sum_0pi1p_total[j][1] += h_stack_Enu_IPK_dp_vector_sum_0pi1p_total->GetBinContent(i+1, j+1)    *  weight[i];
					}
					
				}
				
			}
			if(neutrino==-1){
				
				if(i >= nbin/2){
					
					// multiplicity //
					for(int j = 0; j < nbin_mul; j++){
						mul_total[j][1] += h_stack_Enu_mul_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_mul_pipm; j++){
						mul_pipm_total[j][1] += h_stack_Enu_mul_pipm_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_mul_proton; j++){
						mul_proton_total[j][1] += h_stack_Enu_mul_proton_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					for(int k = 0; k < nbin_mul_pipm_ppi; k++){
						for(int j = 0; j < nbin_mul_proton_ppi; j++){
							mul_proton_Npipm_total[k][j][1] += h_stack_Enu_mul_proton_Npipm_total[k]->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
						}
					}
					
					
					
					// angle //
					for(int j = 0; j < nbin_angle_muon; j++){
						angle_corr_deg_muon_total[j][1] += h_stack_Enu_angle_corr_deg_muon_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_angle_pipm; j++){
						angle_corr_deg_pipm_total[j][1] += h_stack_Enu_angle_corr_deg_pipm_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_angle_proton; j++){
						angle_corr_deg_proton_total[j][1] += h_stack_Enu_angle_corr_deg_proton_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					
					
					// angle (cos) //
					for(int j = 0; j < nbin_angle_cos_muon; j++){
						angle_corr_cos_muon_total[j][1] += h_stack_Enu_angle_corr_cos_muon_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_angle_cos_pipm; j++){
						angle_corr_cos_pipm_total[j][1] += h_stack_Enu_angle_corr_cos_pipm_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_angle_cos_proton; j++){
						angle_corr_cos_proton_total[j][1] += h_stack_Enu_angle_corr_cos_proton_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					
					
					// momentum //
					for(int j = 0; j < nbin_p_muon; j++){
						p_muon_total[j][1] += h_stack_Enu_p_muon_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_p_pipm; j++){
						p_pipm_total[j][1] += h_stack_Enu_p_pipm_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					for(int j = 0; j < nbin_p_proton; j++){
						p_proton_total[j][1] += h_stack_Enu_p_proton_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					
					
					// CC 0pi 2p //
					for(int j = 0; j < nbin_opang_cos_pp_0pi2p; j++){
						opang_cos_pp_0pi2p_total[j][1] += h_stack_Enu_opang_cos_pp_0pi2p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_total_p_pp_0pi2p; j++){
						total_p_pp_0pi2p_total[j][1] += h_stack_Enu_total_p_pp_0pi2p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p; j++){
						inv_mass_pp_0pi2p_total[j][1] += h_stack_Enu_inv_mass_pp_0pi2p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_p_low_p_0pi2p; j++){
						p_low_p_0pi2p_total[j][1] += h_stack_Enu_p_low_p_0pi2p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_p_high_p_0pi2p; j++){
						p_high_p_0pi2p_total[j][1] += h_stack_Enu_p_high_p_0pi2p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					
					
					// CC 1pi 1p //
					for(int j = 0; j < nbin_inv_mass_pip_1pi1p; j++){
						inv_mass_pip_1pi1p_total[j][1] += h_stack_Enu_inv_mass_pip_1pi1p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
					
					
					// CC 0pi 1p //
					for(int j = 0; j < nbin_STV_dp_0pi1p; j++){
						STV_dp_0pi1p_total[j][1] += h_stack_Enu_STV_dp_0pi1p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_STV_dalpha_0pi1p; j++){
						STV_dalpha_0pi1p_total[j][1] += h_stack_Enu_STV_dalpha_0pi1p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_STV_dphi_0pi1p; j++){
						STV_dphi_0pi1p_total[j][1] += h_stack_Enu_STV_dphi_0pi1p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_IPK_dp_scalar_sum_0pi1p; j++){
						IPK_dp_scalar_sum_0pi1p_total[j][1] += h_stack_Enu_IPK_dp_scalar_sum_0pi1p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_IPK_dtheta_0pi1p; j++){
						IPK_dtheta_0pi1p_total[j][1] += h_stack_Enu_IPK_dtheta_0pi1p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					for(int j = 0; j < nbin_IPK_dp_vector_sum_0pi1p; j++){
						IPK_dp_vector_sum_0pi1p_total[j][1] += h_stack_Enu_IPK_dp_vector_sum_0pi1p_total->GetBinContent(i+1-nbin/2, j+1)    *  weight[i];
					}
					
				}
				
			}
			
			
		}
		
		
		
		// ratio //
		
		// multiplicity //
		for(int j = 0; j < nbin_mul; j++){
			//if(fabs(mul_total[j][0])<0.0001 && fabs(mul_total[j][1])<0.0001)
			//	fluxerr_mul_total[j] = 0.0;
			//else
			//	fluxerr_mul_total[j] =    mul_total[j][1]        /  mul_total[j][0];
			fluxerr_mul_total[j] =    mul_total[j][1]        /  mul_total[j][0];
		}
		for(int j = 0; j < nbin_mul_pipm; j++){
			//if(fabs(mul_pipm_total[j][0])<0.0001 && fabs(mul_pipm_total[j][1])<0.0001)
			//	fluxerr_mul_pipm_total[j] = 0.0;
			//else
			//	fluxerr_mul_pipm_total[j] =    mul_pipm_total[j][1]   /  mul_pipm_total[j][0];
			fluxerr_mul_pipm_total[j] =    mul_pipm_total[j][1]   /  mul_pipm_total[j][0];
		}
		for(int j = 0; j < nbin_mul_proton; j++){
			//if(fabs(mul_proton_total[j][0])<0.0001 && fabs(mul_proton_total[j][1])<0.0001)
			//	fluxerr_mul_proton_total[j] = 0.0;
			//else
			//	fluxerr_mul_proton_total[j] =    mul_proton_total[j][1] /  mul_proton_total[j][0];
			fluxerr_mul_proton_total[j] =    mul_proton_total[j][1] /  mul_proton_total[j][0];
		}
		
		for(int k = 0; k < nbin_mul_pipm_ppi; k++){
			for(int j = 0; j < nbin_mul_proton_ppi; j++){
				//if(fabs(mul_proton_Npipm_total[k][j][0])<0.0001 && fabs(mul_proton_Npipm_total[k][j][1])<0.0001)
				//	fluxerr_mul_proton_Npipm_total[k][j] = 0.0;
				//else
				//	fluxerr_mul_proton_Npipm_total[k][j] =    mul_proton_Npipm_total[k][j][1] /  mul_proton_Npipm_total[k][j][0];
				fluxerr_mul_proton_Npipm_total[k][j] =    mul_proton_Npipm_total[k][j][1] /  mul_proton_Npipm_total[k][j][0];
			}
		}
		
		
		
		// angle //
		for(int j = 0; j < nbin_angle_muon; j++){
			//if(fabs(angle_corr_deg_muon_total[j][0])<0.0001 && fabs(angle_corr_deg_muon_total[j][1])<0.0001)
			//	fluxerr_angle_corr_deg_muon_total[j] = 0.0;
			//else
			//	fluxerr_angle_corr_deg_muon_total[j] =    angle_corr_deg_muon_total[j][1]    /  angle_corr_deg_muon_total[j][0];
			fluxerr_angle_corr_deg_muon_total[j] =    angle_corr_deg_muon_total[j][1]    /  angle_corr_deg_muon_total[j][0];
		}
		
		for(int j = 0; j < nbin_angle_pipm; j++){
			//if(fabs(angle_corr_deg_pipm_total[j][0])<0.0001 && fabs(angle_corr_deg_pipm_total[j][1])<0.0001)
			//	fluxerr_angle_corr_deg_pipm_total[j] = 0.0;
			//else
			//	fluxerr_angle_corr_deg_pipm_total[j] =    angle_corr_deg_pipm_total[j][1]    /  angle_corr_deg_pipm_total[j][0];
			fluxerr_angle_corr_deg_pipm_total[j] =    angle_corr_deg_pipm_total[j][1]    /  angle_corr_deg_pipm_total[j][0];
		}
		
		for(int j = 0; j < nbin_angle_proton; j++){
			//if(fabs(angle_corr_deg_proton_total[j][0])<0.0001 && fabs(angle_corr_deg_proton_total[j][1])<0.0001)
			//	fluxerr_angle_corr_deg_proton_total[j] = 0.0;
			//else
			//	fluxerr_angle_corr_deg_proton_total[j] =    angle_corr_deg_proton_total[j][1]  /  angle_corr_deg_proton_total[j][0];
			fluxerr_angle_corr_deg_proton_total[j] =    angle_corr_deg_proton_total[j][1]  /  angle_corr_deg_proton_total[j][0];
		}
		
		
		
		// angle (cos) //
		for(int j = 0; j < nbin_angle_cos_muon; j++){
			fluxerr_angle_corr_cos_muon_total[j] =    angle_corr_cos_muon_total[j][1]    /  angle_corr_cos_muon_total[j][0];
		}
		
		for(int j = 0; j < nbin_angle_cos_pipm; j++){
			fluxerr_angle_corr_cos_pipm_total[j] =    angle_corr_cos_pipm_total[j][1]    /  angle_corr_cos_pipm_total[j][0];
		}
		
		for(int j = 0; j < nbin_angle_cos_proton; j++){
			fluxerr_angle_corr_cos_proton_total[j] =    angle_corr_cos_proton_total[j][1]  /  angle_corr_cos_proton_total[j][0];
		}
		
		
		
		// momentum //
		for(int j = 0; j < nbin_p_muon; j++){
			//if(fabs(p_muon_total[j][0])<0.0001 && fabs(p_muon_total[j][1])<0.0001)
			//	fluxerr_p_muon_total[j] = 0.0;
			//else
			//	fluxerr_p_muon_total[j] =    p_muon_total[j][1]                 /  p_muon_total[j][0];
			fluxerr_p_muon_total[j] =    p_muon_total[j][1]                 /  p_muon_total[j][0];
		}
		
		for(int j = 0; j < nbin_p_pipm; j++){
			//if(fabs(p_pipm_total[j][0])<0.0001 && fabs(p_pipm_total[j][1])<0.0001)
			//	fluxerr_p_pipm_total[j] = 0.0;
			//else
			//	fluxerr_p_pipm_total[j]                    =    p_pipm_total[j][1]                 /  p_pipm_total[j][0];
			fluxerr_p_pipm_total[j]                    =    p_pipm_total[j][1]                 /  p_pipm_total[j][0];
		}
		
		for(int j = 0; j < nbin_p_proton; j++){
			//if(fabs(p_proton_total[j][0])<0.0001 && fabs(p_proton_total[j][1])<0.0001)
			//	fluxerr_p_proton_total[j] = 0.0;
			//else
			//	fluxerr_p_proton_total[j]                  =    p_proton_total[j][1]               /  p_proton_total[j][0];
			fluxerr_p_proton_total[j]                  =    p_proton_total[j][1]               /  p_proton_total[j][0];
		}
		
		
		
		// CC 0pi 2p //
		for(int j = 0; j < nbin_opang_cos_pp_0pi2p; j++){
			fluxerr_opang_cos_pp_0pi2p_total[j]          =    opang_cos_pp_0pi2p_total[j][1]               /  opang_cos_pp_0pi2p_total[j][0];
		}
		for(int j = 0; j < nbin_total_p_pp_0pi2p; j++){
			fluxerr_total_p_pp_0pi2p_total[j]          =    total_p_pp_0pi2p_total[j][1]               /  total_p_pp_0pi2p_total[j][0];
		}
		for(int j = 0; j < nbin_inv_mass_pp_0pi2p/rebin_inv_mass_pp_0pi2p; j++){
			fluxerr_inv_mass_pp_0pi2p_total[j]          =    inv_mass_pp_0pi2p_total[j][1]               /  inv_mass_pp_0pi2p_total[j][0];
		}
		for(int j = 0; j < nbin_p_low_p_0pi2p; j++){
			fluxerr_p_low_p_0pi2p_total[j]          =    p_low_p_0pi2p_total[j][1]               /  p_low_p_0pi2p_total[j][0];
		}
		for(int j = 0; j < nbin_p_high_p_0pi2p; j++){
			fluxerr_p_high_p_0pi2p_total[j]          =    p_high_p_0pi2p_total[j][1]               /  p_high_p_0pi2p_total[j][0];
		}
		
		
		
		// CC 1pi 1p //
		for(int j = 0; j < nbin_inv_mass_pip_1pi1p; j++){
			fluxerr_inv_mass_pip_1pi1p_total[j]          =    inv_mass_pip_1pi1p_total[j][1]               /  inv_mass_pip_1pi1p_total[j][0];
		}
		
		
		
		// CC 0pi 1p //
		for(int j = 0; j < nbin_STV_dp_0pi1p; j++){
			fluxerr_STV_dp_0pi1p_total[j]          =    STV_dp_0pi1p_total[j][1]               /  STV_dp_0pi1p_total[j][0];
		}
		for(int j = 0; j < nbin_STV_dalpha_0pi1p; j++){
			fluxerr_STV_dalpha_0pi1p_total[j]          =    STV_dalpha_0pi1p_total[j][1]               /  STV_dalpha_0pi1p_total[j][0];
		}
		for(int j = 0; j < nbin_STV_dphi_0pi1p; j++){
			fluxerr_STV_dphi_0pi1p_total[j]          =    STV_dphi_0pi1p_total[j][1]               /  STV_dphi_0pi1p_total[j][0];
		}
		for(int j = 0; j < nbin_IPK_dp_scalar_sum_0pi1p; j++){
			fluxerr_IPK_dp_scalar_sum_0pi1p_total[j]          =    IPK_dp_scalar_sum_0pi1p_total[j][1]               /  IPK_dp_scalar_sum_0pi1p_total[j][0];
		}
		for(int j = 0; j < nbin_IPK_dtheta_0pi1p; j++){
			fluxerr_IPK_dtheta_0pi1p_total[j]          =    IPK_dtheta_0pi1p_total[j][1]               /  IPK_dtheta_0pi1p_total[j][0];
		}
		for(int j = 0; j < nbin_IPK_dp_vector_sum_0pi1p; j++){
			fluxerr_IPK_dp_vector_sum_0pi1p_total[j]          =    IPK_dp_vector_sum_0pi1p_total[j][1]               /  IPK_dp_vector_sum_0pi1p_total[j][0];
		}
		
		
		tree->Fill();
	}
	
	fprintf(stderr,"\r%10lld / %10lld\n", nthrows, nthrows);
	
	
	
	
	
	// Write //
	root_o1->cd();
	
	tree->Write();
	h_stack_Enu_mul_total                      -> Write();
	h_stack_Enu_mul_pipm_total                 -> Write();
	h_stack_Enu_mul_proton_total               -> Write();
	h_stack_Enu_angle_corr_deg_muon_total      -> Write();
	h_stack_Enu_angle_corr_deg_pipm_total      -> Write();
	h_stack_Enu_angle_corr_deg_proton_total    -> Write();
	h_stack_Enu_angle_corr_cos_muon_total      -> Write();
	h_stack_Enu_angle_corr_cos_pipm_total      -> Write();
	h_stack_Enu_angle_corr_cos_proton_total    -> Write();
	h_stack_Enu_p_muon_total                   -> Write();
	h_stack_Enu_p_proton_total                 -> Write();
	h_stack_Enu_p_pipm_total                   -> Write();
	h_stack_Enu_opang_cos_pp_0pi2p_total       -> Write();
	h_stack_Enu_total_p_pp_0pi2p_total         -> Write();
	h_stack_Enu_inv_mass_pp_0pi2p_total        -> Write();
	h_stack_Enu_p_low_p_0pi2p_total            -> Write();
	h_stack_Enu_p_high_p_0pi2p_total           -> Write();
	h_stack_Enu_inv_mass_pip_1pi1p_total       -> Write();
	h_stack_Enu_STV_dp_0pi1p_total             -> Write();
	h_stack_Enu_STV_dalpha_0pi1p_total         -> Write();
	h_stack_Enu_STV_dphi_0pi1p_total           -> Write();
	h_stack_Enu_IPK_dp_scalar_sum_0pi1p_total  -> Write();
	h_stack_Enu_IPK_dtheta_0pi1p_total         -> Write();
	h_stack_Enu_IPK_dp_vector_sum_0pi1p_total  -> Write();
	
	root_o1->Close();
	
	
	return 0;
	
	
}


void cholcov_conv(float covmat[nbin][nbin], float cholcovmat[nbin][nbin]){
	
	memset(cholcovmat,0,sizeof(cholcovmat));
	for ( Int_t j=0; j<nbin; j++ ) { 
		Double_t s = covmat[j][j] ;
		for ( Int_t k=0; k<j; k++ ) {
			s -= cholcovmat[j][k]*cholcovmat[j][k] ;
		}
		if(s<0){
			std::cout << "strange !" << j << " " << s << std::endl ;
			exit(0) ;
		}
		cholcovmat[j][j] = sqrt(s) ;
		
		for ( Int_t i=j+1; i<nbin; i++ ) {
			s = covmat[i][j] ;
			for ( Int_t k=0; k<j; k++ ) {
				s -= cholcovmat[i][k]*cholcovmat[j][k] ;
			}
			if ( TMath::Abs(s)<0.000000000001 )
				cholcovmat[i][j] = 0. ;
			else
				cholcovmat[i][j] = s/cholcovmat[j][j] ;
			
			//cout<<cholcovmat[i][j]<<endl;
		}
	}
	
}

void draw_matrix(float mat[nbin][nbin]){
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
	TH2F *h1=new TH2F("h1","h1",nbin,0,nbin,nbin,0,nbin);
	h1->GetZaxis()->SetRangeUser(-0.045,0.045);
	h1->Reset();
	for(int i=0;i<nbin;i++){
		for(int j=0;j<nbin;j++){
			h1->Fill(i+0.0001,j+0.0001,mat[i][j]);
		}
	}
	h1->Draw("colz");
	c1->Update();
	
	char ans[8];
	fgets(ans,8,stdin);
	
}
