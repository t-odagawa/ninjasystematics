void DrawResult() {

  // TString filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode0.root.hist_flux_pi.toymc.root";
  // TString filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode2.root.hist_flux_pi_mcs.toymc.root";
  TString filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode3.root.hist_flux_pi_deg.toymc.root";
  TFile *file = new TFile(filename, "read");
  auto hsignal = (TH2D*)file->Get("hsignal");
  auto tree = (TTree*)file->Get("tree");

  TString ofilename = "/hsm/nu/ninja/pra_tmp/fluxsys/bin_err/hist_flux_pion_deg.root";
  TFile *ofile = new TFile(ofilename, "recreate");
  
  const int nbinsx = hsignal->GetNbinsX();
  double x_bins[nbinsx + 1];
  hsignal->ProjectionX()->GetLowEdge(x_bins);
  x_bins[nbinsx] = x_bins[nbinsx - 1] + hsignal->ProjectionX()->GetBinWidth(nbinsx);

  for ( int i = 0; i <= nbinsx; i++ ) cout << x_bins[i] << endl;

  TH1D* h = new TH1D("h", "", nbinsx, x_bins);

  int n = tree->GetEntries();
  double signal;
  double observable;
  tree->SetBranchAddress("signal", &signal);
  tree->SetBranchAddress("observable", &observable);

  for ( int i = 1; i < nbinsx + 1; i++ ) {
    cout << h->GetBinCenter(i) << endl;
    TH1D *tmp = new TH1D("tmp", "", 100, 0., 2.);
    for ( int j = 0; j < n; j++ ) {
      tree->GetEntry(j);
      if ( observable == h->GetBinCenter(i) ) tmp->Fill(signal);
    }
    
    if ( tmp->GetMean() > 0 ) {
      tmp->Fit("gaus", "Q0", "", 0., 2.);
      auto gaus = (TF1*)gROOT->FindObject("gaus");
      double err = gaus->GetParameter(2);
      h->Fill(h->GetBinCenter(i), err);
    }
    else {
      h->Fill(h->GetBinCenter(i), 0.);
    }
    tmp->Delete();

  }

  h->GetYaxis()->SetRangeUser(0, 0.2);
  h->GetYaxis()->SetTitle("Fractional error");
  h->GetXaxis()->SetTitle(hsignal->GetXaxis()->GetTitle());
  h->SetStats(0);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  // h->Draw("hist");

  ofile->cd();
  h->Write();
  ofile->Close();

}
