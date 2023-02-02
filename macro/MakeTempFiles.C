TString GetGeneralFileName(int i);
TString GetBinErrFileName(int i);

void MakeTempFiles() {

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

  // const double temporal_error = 0.01;
  const double temporal_error = 0.05;

  for ( int kinematics = 0; kinematics < 15; kinematics++ ) {
    
    TString errfilename = "/hsm/nu/ninja/pra_tmp/detsys/bin_err/BabyMind";
    //TString errfilename = "/hsm/nu/ninja/pra_tmp/detsys/bin_err/MaterialThickness";
    TString covfilename = "/hsm/nu/ninja/pra_tmp/detsys/covariance_matrix/BabyMind";
    //TString covfilename = "/hsm/nu/ninja/pra_tmp/detsys/covariance_matrix/MaterialThickness";

    errfilename += GetBinErrFileName(kinematics);
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

    covfilename += GetGeneralFileName(kinematics);
    TFile *covfile = new TFile(covfilename, "recreate");
    TMatrixDSym *mat = new TMatrixDSym(hist->GetNbinsX());

    for ( int ibin = 0; ibin < mat->GetNcols(); ibin++ ) {
      hist->SetBinContent(ibin+1, temporal_error);
      for ( int jbin = 0; jbin < mat->GetNcols(); jbin++ ) {
	(*mat)(ibin, jbin) = temporal_error * temporal_error;
      }
    }

    covfile->cd();
    mat->Write("cov_mat");
    covfile->Close();

    errfile->cd();
    errfile->Write();
    errfile->Close();  
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
    return "/hist_det_total_multi.root";
  case 1 :
    return "/hist_det_proton_multi.root";
  case 2 :
    return "/hist_det_pion_multi.root";
  case 3 :
    return "/hist_det_muon_mom_total.root";
  case 4 :
    return "/hist_det_muon_mom_range.root";
  case 5 :
    return "/hist_det_muon_mom_mcs.root";
  case 6 :
    return "/hist_det_proton_mom_total.root";
  case 7 :
    return "/hist_det_proton_mom_range.root";
  case 8 :
    return "/hist_det_proton_mom_mcs.root";
  case 9 :
    return "/hist_det_pion_mom_total.root";
  case 10 :
    return "/hist_det_pion_mom_range.root";
  case 11 :
    return "/hist_det_pion_mom_mcs.root";
  case 12 :
    return "/hist_det_muon_deg.root";
  case 13 :
    return "/hist_det_proton_deg.root";
  case 14 :
    return "/hist_det_pion_deg.root";
  default :
    throw std::runtime_error("Invalid kinematics id : " + std::to_string(i));     
  }
}
