void FluxReWeight() {

  // TFile *fop = TFile::Open("/home/t2k/odagawa/NinjaSystematics/flux_tuning/tuned_bv21v2.root");
  TFile *fop = TFile::Open("/home/t2k/odagawa/NinjaSystematics/flux_tuning_wall/tuned_bv21v2.root");
  
  auto *numu_nom = (TH1D*)fop->Get("nd2_nom_numu");
  auto *numu_tune = (TH1D*)fop->Get("nd2_tune_numu");

  const Int_t nu_ene_bin_size = 20;
  const Double_t nu_ene_bins[nu_ene_bin_size] = {0.0, 0.1, 0.2, 0.3, 0.4,
						 0.5, 0.6, 0.7, 0.8, 1.0,
						 1.2, 1.5, 2.0, 2.5, 3.0,
						 3.5, 4.0, 5.0, 7.0, 10.0};
  
  auto nd7_nom_numu_rebin = numu_nom->Rebin(nu_ene_bin_size-1, "nd2_nom_numu_rebin", nu_ene_bins);
  auto nd7_tune_numu_rebin = numu_tune->Rebin(nu_ene_bin_size-1, "nd2_tune_numu_rebin", nu_ene_bins);

  nd7_tune_numu_rebin->Divide(nd7_nom_numu_rebin);

  nd7_tune_numu_rebin->Draw("hist e");

  for ( int ibin = 1; ibin < nu_ene_bin_size; ibin++ ) {
    std::cout << "Bin " << ibin << " : " << nd7_tune_numu_rebin->GetBinContent(ibin) << std::endl;
  }

  TString ofilename = "numu_flux_reweight_wall.root";
  TFile *ofile = new TFile(ofilename, "recreate");

  nd7_tune_numu_rebin->Write();
  
  ofile->Close();

}
