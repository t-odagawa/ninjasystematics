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

TString GetParameterName(int i) {
  switch (i) {
  case 0 :
    return "M_{A}^{QE}";
  case 1 :
    return "Reweight CCQE (0.25 < Q^{2} < 0.5)";
  case 2 :
    return "Reweight CCQE (0.5 < Q^{2} < 1.0)";
  case 3 :
    return "Reweight CCQE (1.0 < Q^{2})";
  case 4 :
    return "SF Optical Potential Correction";
  case 5 :
    return "2p2h Normalization";
  case 6 :
    return "2p2h Shape (delta-like/NN-like for nn)";
  case 7 :
    return "2p2h Shape (delta-like/NN-like for np)";
  case 8 :
    return "2p2h Shape (nn or np)";
  case 9 :
    return "Resonant SPP E_{B}";
  case 10 :
    return "Isospin 1/2 non-resonant bkg.";
  case 11 :
    return "C_{5}^{A}(0)";
  case 12 :
    return "M_{A}^{RES}";
  case 13 : 
    return "Pion FSI Absorption";
  case 14 :
    return "Pion FSI Single Charge Exchange (High E)";
  case 15 :
    return "Pion FSI Single Charge Exchange (Low E)";
  case 16 :
    return "Pion FSI Hadron Production";
  case 17 :
    return "Pion FSI QE Scattering (High E)";
  case 18 :
    return "Pion FSI QE Scattering (Low E)";
  case 19 :
    return "FSI Nucleon Fate";
  case 20 :
    return "DIS/Multi Pion Normalization";
  default :
    throw std::runtime_error("Invalid dial id : " + std::to_string(i));
  }
}

int GetLineColor(int i) {
  switch (i) {
  case 0 :
  case 1 :
  case 2 :
  case 3 :
  case 4 :
    return 424;
  case 5 :
  case 6 :
  case 7 :
  case 8 :
    return 624;
  case 9 :
  case 10 :
  case 11 :
  case 12 :
    return 394;
  case 13 :
  case 14 :
  case 15 :
  case 16 :
  case 17 :
  case 18 :
  case 19 :
    return kRed;
  case 20 :
    return kBlue;
  default :
    throw std::runtime_error("Invalid dial id : " + std::to_string(i));
  }
}

int GetLineStyle(int i) {
  switch (i) {
  case 0 : return 1;
  case 1 : return 2;
  case 2 : return 3;
  case 3 : return 4;
  case 4 : return 5;
  case 5 : return 1;
  case 6 : return 2;
  case 7 : return 3;
  case 8 : return 4;
  case 9 : return 1;
  case 10 : return 2;
  case 11 : return 3;
  case 12 : return 4;
  case 13 : return 1;
  case 14 : return 2;
  case 15 : return 3;
  case 16 : return 4;
  case 17 : return 5;
  case 18 : return 6;
  case 19 : return 7;
  case 20 : return 1;
  default : throw std::runtime_error("Invalid dial id : " + std::to_string(i));
  }
}

void CompareInteractionError() {  

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

  gStyle->SetOptStat(0);

  TLegend *leg = new TLegend(0.1, 0.1, 0.9, 0.9);
  TH1D *hist[21];
  double highest_entry = 0.0;
  int nbins = 0;

  TCanvas *c1 = new TCanvas("c1", "c1", 2000, 800);
  c1->Divide(2, 1);
  c1->cd(1);

  for ( int i = 0; i < 21; i++ ) {
    std::cout << i << std::endl;
    TString filename = "/hsm/nu/ninja/pra_tmp/xsecsys/covariance_matrix";
    filename += "/";
    filename += GetDialName(i);
    // filename += "/hist_pion_mom_covmat.root";
    filename += "/hist_water_proton_multi_covmat.root";
    // filename += "/hist_proton_ang_deg_covmat.root";

    TFile *file = new TFile(filename, "read");
    TMatrixDSym *mat = (TMatrixDSym*)file->Get("cov_mat");
    if ( i == 0 )
      nbins = mat->GetNcols();
    hist[i] = new TH1D(Form("hist%d",i), "", nbins, 0, nbins);
    // hist[i] = new TH1D(Form("hist%d",i), "", hadron_deg_bin_size-1, hadron_deg_bins);
    hist[i]->SetLineColor(GetLineColor(i));
    hist[i]->SetLineStyle(GetLineStyle(i));
    hist[i]->SetLineWidth(3);
    leg->AddEntry(hist[i], GetParameterName(i), "l");
    for ( int ibin = 0; ibin < nbins; ibin++ ) {
      double error = std::sqrt((*mat)(ibin, ibin));
      hist[i]->SetBinContent(ibin + 1, error);
      if ( error > highest_entry ) highest_entry = error;
    }   

    if ( i == 0 )
      hist[i]->Draw("hist");
    else
      hist[i]->Draw("hist same");

  }

  hist[0]->SetTitle(";Bin id;Fractional error");
  // hist[0]->SetTitle(";Momentum [MeV/c];Fractional error");
  // hist[0]->SetTitle(";Angle [deg];Fractional error");
  hist[0]->GetYaxis()->SetRangeUser(0, highest_entry + 0.1);
  hist[0]->GetYaxis()->CenterTitle();
  hist[0]->GetXaxis()->CenterTitle();

  gPad->SetGrid();

  c1->cd(2);
  leg->Draw();

  TFile *ofile = new TFile("~/test_bin_err.root", "recreate");
  
  leg->Write("leg");
  for ( int i = 0; i < 1; i++ ) 
    hist[i]->Write();
  c1->Write();
  ofile->Close();

}
