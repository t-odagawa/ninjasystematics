void DrawFluxCovariance() {

  gStyle->SetOptStat(0);
  const int nstep = 5;
  double stops[nstep] = {0.0, 0.25, 0.5, 0.75, 1.0};
  double red[nstep] = {0.0, 0.8, 1.0, 1.0, 1.0};
  double green[nstep] = {0.0, 0.8, 1.0, 0.8, 0.0};
  double blue[nstep] = {1.0, 1.0, 1.0, 0.8, 0.0};
  TColor::CreateGradientColorTable(nstep,stops,red,green,blue,100);

  TString ifilename = "/hsm/nu/ninja/pra_tmp/fluxsys/bin_err/hist_flux_total_multi.root";
  TFile *ifile = new TFile(ifilename, "read");
  auto h = (TH1D*)ifile->Get("h");

  TString ofilename = "/hsm/nu/ninja/pra_tmp/fluxsys/covariance_matrix/cov_flux_total_multi.root";
  TFile *ofile = new TFile(ofilename, "recreate");

  TMatrixDSym *mat = new TMatrixDSym(h->GetNbinsX());

  for ( int ibin = 1; ibin <= h->GetNbinsX(); ibin++ ) {    
    for ( int jbin = 1; jbin <= h->GetNbinsX(); jbin++ ) {
      (*mat)(ibin-1, jbin-1) = h->GetBinContent(ibin) * h->GetBinContent(jbin);
    }
  }

  // mat->Draw("col text");

  ofile->cd();
  mat-Write("mat");
  ofile->Close();

}
