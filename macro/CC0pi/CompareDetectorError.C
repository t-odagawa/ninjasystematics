TString GetDialName(int i);
TString GetParameterName(int i);
int GetLineColor(int i);
int GetLineStyle(int i);

void CompareDetectorError() {

  gStyle->SetOptStat(0);

  const Int_t numpar = 12;
  const std::vector<Int_t> pars = {0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11};

  TLegend *leg = new TLegend(0.1, 0.1, 0.9, 0.9);
  TH1D *hist[numpar];
  double highest_entry = 0.0;
  int nbins = 0;

  TCanvas *c1 = new TCanvas("c1", "c1", 2000, 800);
  c1->Divide(2,1);
  c1->cd(1);

  for ( auto i : pars ) {
    // for ( int i = 0; i < numpar; i++ ) {
    std::cout << i << std::endl;
    if ( i == 7 ) continue;
    TString filename = "/hsm/nu/ninja/pra_tmp/CC0pi_20221213/detsys/covariance_matrix";
    filename += "/";
    filename += GetDialName(i);
    filename += "/hist_water_proton_multi_covmat.root";
    // filename += "/hist_water_total_multi_covmat.root";

    TFile *file = new TFile(filename, "read");
    TMatrixDSym *mat = (TMatrixDSym*)file->Get("cov_mat");
    if ( i == 0 )
      nbins = mat->GetNcols();
    hist[i] = new TH1D(Form("hist_%d",i), "", nbins, 0, nbins);
    hist[i]->SetLineColor(GetLineColor(i));
    hist[i]->SetLineStyle(GetLineStyle(i));
    hist[i]->SetLineWidth(3);
    leg->AddEntry(hist[i], GetParameterName(i), "l");
    for ( int ibin = 0; ibin < nbins; ibin++ ) {
      double error = std::sqrt((*mat)(ibin, ibin));
      hist[i]->Fill(ibin + 0.5, error);
      if ( error > highest_entry ) highest_entry = error;
    }

    if ( i == 0 )
      hist[i]->Draw("hist");
    else
      hist[i]->Draw("hist same");
  }

  hist[0]->SetTitle(";Bin id;Fractional error");
  hist[0]->GetYaxis()->SetRangeUser(0, highest_entry + 0.1);
  hist[0]->GetYaxis()->CenterTitle();
  hist[0]->GetXaxis()->CenterTitle();

  gPad->SetGrid();
  
  c1->cd(2);
  leg->Draw();

}

TString GetDialName(int i) {
  switch (i) {
  case 0 : return "BabyMind";
  case 1 : return "NinjaBabyMindDistance";
  case 2 : return "HitThreshold";
  case 3 : return "MPPCNoise";
  case 4 : return "Alignment";
  case 5 : return "AngularResolution";
  case 6 : return "DetectionEfficiency";
  case 7 : return "McsScaling";
  case 8 : return "VphMean";
  case 9 : return "VphSigma";
  case 10 : return "MaterialThickness";    
  case 11 : return "PhysicsList";
  default : throw std::runtime_error("Invalid parameter id : " + std::to_string(i));
  }
}

TString GetParameterName(int i) {
  switch (i) {
  case 0 : return "Baby MIND";
  case 1 : return "Distance b/w NJ & BM";
  case 2 : return "Tracker light yield";
  case 3 : return "Tracker MPPC noise";
  case 4 : return "Tracker scintillator alignment";
  case 5 : return "ECC film angular resolution";
  case 6 : return "ECC film detection efficiency";
  case 7 : return "S paramter in MCS";
  case 8 : return "VPH PDF mean";
  case 9 : return "VPH PDF sigma";
  case 10 : return "ECC detector structure";
  case 11 : return "Geant4 physics list";
  default : throw std::runtime_error("Invalid paraemter id : " + std::to_string(i));
  }
}

int GetLineColor(int i) {
  switch (i) {
  case 0 :
  case 1 :
    return kOrange;
  case 2 :
  case 3 :
  case 4 :
    return kRed;
  case 5 :
  case 6 :
    return kBlue;
  case 7 :
  case 8 :
  case 9 :
    return kMagenta;
  case 10 :
  case 11 :
    return kGray;
  default : throw std::runtime_error("Invalid parameter id : " + std::to_string(i));
  }
}

int GetLineStyle(int i) {
  switch (i) {
  case 0 : return 1;
  case 1 : return 2;
  case 2 : return 1;
  case 3 : return 2;
  case 4 : return 3;
  case 5 : return 1;
  case 6 : return 2;
  case 7 : return 1;
  case 8 : return 2;
  case 9 : return 3;
  case 10 : return 1;
  case 11 : return 2;
  default : throw std::runtime_error("Invalid parameter id : " + std::to_string(i));
  }
}
