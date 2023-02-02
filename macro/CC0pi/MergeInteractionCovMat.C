TString GetGeneralFileName(int i);
TString GetBinErrFileName(int i);
TString GetDialName(int i);

void MergeInteractionCovMat() {

  int kinematics = 13;

  TString general_filename = GetGeneralFileName(kinematics);
  TString bin_err_filename = GetBinErrFileName(kinematics);

  TString ofilename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/xsecsys/covariance_matrix";
  ofilename += general_filename;
  TFile *ofile = new TFile(ofilename, "recreate");
  TMatrixDSym *mat_tot;

  const Int_t total_multi_bin_size = 7;
  const double total_multi_bins[total_multi_bin_size] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5};
  const Int_t hadron_multi_bin_size = 7;
  const double hadron_multi_bins[hadron_multi_bin_size] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 9.5};
  const Int_t muon_mom_bin_size = 7;
  const double muon_mom_bins[muon_mom_bin_size] = {0., 400., 600., 800., 1000., 1200., 2000.};
  const Int_t hadron_mom_bin_size = 6;
  const double hadron_mom_bins[hadron_mom_bin_size] = {0., 200., 400., 600., 800., 1000.};
  const Int_t muon_deg_bin_size = 7;
  const double muon_deg_bins[muon_deg_bin_size] = {0.0, 10., 20., 30., 40., 50., 90.};
  const Int_t hadron_deg_bin_size = 10;
  const double hadron_deg_bins[hadron_deg_bin_size] = {0.0, 20., 40., 60., 80., 100., 120., 140., 160., 180.};


  TString errfilename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/xsecsys/bin_err";
  errfilename += bin_err_filename;
  TFile *errfile = new TFile(errfilename, "recreate");
  TH1D *hist;
  if ( kinematics == 0 )
    hist = new TH1D("h", ";# of tracks;Fractional error", total_multi_bin_size-1, total_multi_bins);
  else if ( kinematics < 3 )
    hist = new TH1D("h", ";# of tracks;Fractional error", hadron_multi_bin_size-1, hadron_multi_bins);
  else if ( kinematics < 6 )
    hist = new TH1D("h", ";# of tracks;Fractional error", muon_mom_bin_size-1, muon_mom_bins);
  else if ( kinematics < 12 )
    hist = new TH1D("h", ";# of tracks;Fractional error", hadron_mom_bin_size-1, hadron_mom_bins);
  else if ( kinematics < 13 )
    hist = new TH1D("h", ";# of tracks;Fractional error", muon_deg_bin_size-1, muon_deg_bins);
  else if ( kinematics < 15 )
    hist = new TH1D("h", ";# of tracks;Fractional error", hadron_deg_bin_size-1, hadron_deg_bins);

  for ( int i = 0; i < 21; i++ ) {
    TString filename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/xsecsys/covariance_matrix";
    filename += "/";
    filename += GetDialName(i);
    filename += general_filename;

    TFile *file = new TFile(filename, "read");
    TMatrixDSym *mat = (TMatrixDSym*)file->Get("cov_mat");

    if ( i == 0 ) {
      mat_tot = new TMatrixDSym(mat->GetNcols());
    }

    for ( int ibin = 0; ibin < mat_tot->GetNcols(); ibin++ ) 
      for ( int jbin = 0; jbin < mat_tot->GetNcols(); jbin++ ) 
      	(*mat_tot)(ibin, jbin) += (*mat)(ibin, jbin);

  }

  for ( int ibin = 0; ibin < hist->GetNbinsX(); ibin++ ) {
    hist->Fill(hist->GetBinCenter(ibin+1), std::sqrt((*mat_tot)(ibin, ibin)));
  }

  ofile->cd();
  mat_tot->Write("cov_mat");
  ofile->Close();

  errfile->cd();
  hist->Write();
  errfile->Close();

}

TString GetDialName(int i) {
  switch (i) {
  case 0 :
    return "MaCCQE";
  case 1 :
    return "QETwk_HighQ2Weight_1";
  case 2 :
    return "QETwk_HighQ2Weight_2";
  case 3 :
    return "QETwk_HighQ2Weight_3";
  case 4 :
    return "SF_OptPotTwkDial_O16";
  case 5 :
    return "MECTwkDial_Norm_O16";
  case 6 :
    return "MECTwkDial_PDDWeight_O16_NN";
  case 7 :
    return "MECTwkDial_PDDWeight_O16_np";
  case 8 :
    return "MECTwkDial_PNNN_Shape";
  case 9 :
    return "RES_Eb_O_numu";
  case 10 :
    return "BgSclRES";
  case 11 :
    return "CA5RES";
  case 12 :
    return "MaRES";
  case 13 :
    return "PionFSI_AbsProb";
  case 14 :
    return "PionFSI_CExHighMomProb";
  case 15 :
    return "PionFSI_CExLowMomProb";
  case 16 :
    return "PionFSI_InelProb";
  case 17 :
    return "PionFSI_QEHighMomProb";
  case 18 :
    return "PionFSI_QELowMomProb";
  case 19 :
    return "TwkDial_FateNucleonFSI";
  case 20 :
    return "CC_DIS_MultiPi_Norm_Nu";
  default :
    throw std::runtime_error("Invalid dial id : " + std::to_string(i));
  }
}

TString GetGeneralFileName(int i) {
  switch (i) {
  case 0 :
    return "/hist_water_total_multi_covmat.root";
  case 1 :
    return "/hist_water_proton_multi_covmat.root";
  case 2 :
    return "/hist_water_pion_multi_covmat.root";
  case 3 :
    return "/hist_muon_mom_covmat.root";
  case 4 :
    return "/hist_muon_mom_range_covmat.root";
  case 5 :
    return "/hist_muon_mom_mcs_covmat.root";
  case 6 :
    return "/hist_proton_mom_covmat.root";
  case 7 :
    return "/hist_proton_mom_range_covmat.root";
  case 8 :
    return "/hist_proton_mom_mcs_covmat.root";
  case 9 :
    return "/hist_pion_mom_covmat.root";
  case 10 :
    return "/hist_pion_mom_range_covmat.root";
  case 11 :
    return "/hist_pion_mom_mcs_covmat.root";
  case 12 :
    return "/hist_muon_ang_deg_covmat.root";
  case 13 :
    return "/hist_proton_ang_deg_covmat.root";
  case 14 :
    return "/hist_pion_ang_deg_covmat.root";
  default :
    throw std::runtime_error("Invalid kinematics id : " + std::to_string(i));
  }
}

TString GetBinErrFileName(int i) {
  switch (i) {
  case 0 :
    return "/hist_xsec_total_multi.root";
  case 1 :
    return "/hist_xsec_proton_multi.root";
  case 2 :
    return "/hist_xsec_pion_multi.root";
  case 3 :
    return "/hist_xsec_muon_mom_total.root";
  case 4 :
    return "/hist_xsec_muon_mom_range.root";
  case 5 :
    return "/hist_xsec_muon_mom_mcs.root";
  case 6 :
    return "/hist_xsec_proton_mom_total.root";
  case 7 :
    return "/hist_xsec_proton_mom_range.root";
  case 8 :
    return "/hist_xsec_proton_mom_mcs.root";
  case 9 :
    return "/hist_xsec_pion_mom_total.root";
  case 10 :
    return "/hist_xsec_pion_mom_range.root";
  case 11 :
    return "/hist_xsec_pion_mom_mcs.root";
  case 12 :
    return "/hist_xsec_muon_deg.root";
  case 13 :
    return "/hist_xsec_proton_deg.root";
  case 14 :
    return "/hist_xsec_pion_deg.root";
  default :
    throw std::runtime_error("Invalid kinematics id : " + std::to_string(i));     
  }
}
