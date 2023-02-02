//#define MKPLOT

#define nbin 40

void draw_matrix(double mat[nbin][nbin]) {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 800, 600);
  TH2D *h1 = new TH2D("h1", "h1", nbin, 0, nbin, nbin, 0, nbin);
  h1->GetZaxis()->SetRangeUser(-0.045, 0.045);
  h1->Reset();
  for ( int i = 0; i < nbin; i++ ) {
    for ( int j = 0; j < nbin; j++ ) {
      h1->Fill(i + 0.0001, j + 0.0001, mat[i][j]);
    }
  }
  h1->Draw("colz");
  c1->Update();

  char ans[8];
  fgets(ans, 8, stdin);
};

void cholcov_conv(double covmat[nbin][nbin], double cholcovmat[nbin][nbin]) {

  for ( int j = 0; j < nbin; j++ ) {
    double s = covmat[j][j];
    for ( int k = 0; k < j; k++ ) {
      s -= cholcovmat[j][k] * cholcovmat[j][k];
    }
    if ( s < 0 ) {
      std::cout << "strange !" << j << ", " << s << std::endl;
      std::exit(1);
    }
    cholcovmat[j][j] = TMath::Sqrt(s);

    for ( int i = j + 1; i < nbin; i++ ) {
      s = covmat[i][j];
      for ( int k = 0; k < j; k++ ) {
	s -= cholcovmat[i][k] * cholcovmat[j][k];
      }
      if ( TMath::Abs(s) < 1.e-12 )
	cholcovmat[i][j] = 0.;
      else 
	cholcovmat[i][j] = s / cholcovmat[j][j];
    }
  }
};

void FluxErr() {

  gRandom->SetSeed(time(NULL));

  TRandom3 rand;

  TFile *fop = TFile::Open("~/NinjaSystematics/yasutome-san/flux_covariance_wagasci_2021_total_bv21v2.root");
  TMatrixDSym *mat = (TMatrixDSym*)fop->Get("total_flux_cov");
  if ( mat == NULL ) {
    std::cout << "Cannot get matrix from " << fop->GetPath() << std::endl;
    std::exit(1);
  }

  double cov_mat[nbin][nbin];
  for ( int i = 0; i < nbin; i++ )
    for ( int j = 0; j < nbin; j++ )
      cov_mat[i][j] = (*mat)(i, j);

#ifdef MKPLOT
  draw_matrix(cov_mat);
#endif

  double chol_mat[nbin][nbin];
  cholcov_conv(cov_mat, chol_mat);

#ifdef MKPLOT
  // draw_matrix(chol_mat);
#endif

  int nthrows = 100000;
  double weight[nbin], nrand[nbin];
  TString filename;
  TH2D *hsignal;

  for ( int ifile = 0; ifile < 3; ifile++ ) {
    
    if ( ifile == 0 ) filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode0.root"; // multiplicty
    else if ( ifile == 1 ) filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode2.root"; // momentum
    else if ( ifile == 2 ) filename = "/group/nu/ninja/work/odagawa/20221020-phd-thesis-preliminary/signal/output/output_mode3.root"; // angle

    TFile *f = new TFile(filename, "read");

    int nhist;
    if ( ifile == 0 ) nhist = 3;
    else if ( ifile == 1 ) nhist = 9;
    else if ( ifile == 2 ) nhist = 3;

    for ( int ihist = 0; ihist < nhist; ihist++ ) {
      if ( ifile == 0 ) {
	if ( ihist == 0 ) hsignal = (TH2D*)f->Get("hist_flux");
	else if ( ihist == 1 ) hsignal = (TH2D*)f->Get("hist_flux_p");
	else if ( ihist == 2 ) hsignal = (TH2D*)f->Get("hist_flux_pi");
      }
      else if ( ifile == 1 ) {
	if ( ihist == 0 ) hsignal = (TH2D*)f->Get("hist_flux_mu_total");
	else if ( ihist == 1 ) hsignal = (TH2D*)f->Get("hist_flux_mu_range");
	else if ( ihist == 2 ) hsignal = (TH2D*)f->Get("hist_flux_mu_mcs");
	else if ( ihist == 3 ) hsignal = (TH2D*)f->Get("hist_flux_p_total");
	else if ( ihist == 4 ) hsignal = (TH2D*)f->Get("hist_flux_p_range");
	else if ( ihist == 5 ) hsignal = (TH2D*)f->Get("hist_flux_p_mcs");
	else if ( ihist == 6 ) hsignal = (TH2D*)f->Get("hist_flux_pi_total");
	else if ( ihist == 7 ) hsignal = (TH2D*)f->Get("hist_flux_pi_range");
	else if ( ihist == 8 ) hsignal = (TH2D*)f->Get("hist_flux_pi_mcs");
      }
      else if ( ifile == 2 ) {
	if ( ihist == 0 ) hsignal = (TH2D*)f->Get("hist_flux_mu_deg");
	else if ( ihist == 1 ) hsignal = (TH2D*)f->Get("hist_flux_p_deg");
	else if ( ihist == 2 ) hsignal = (TH2D*)f->Get("hist_flux_pi_deg");
      }

      stringstream fout_ss;
      fout_ss << filename << "." << hsignal->GetName() << ".toymc.root";
      TFile *fout = new TFile(fout_ss.str().c_str(), "recreate"); // ここを自動でかけるようにしておきたい
      TTree *tree = new TTree("tree", "tree");

      double signal, observable;
      
      tree->Branch("signal", &signal, "signal/D");
      tree->Branch("observable", &observable, "observable/D"); 
      
      cout << hsignal->GetNbinsX() << endl;
      
      for ( int l = 0; l < hsignal->GetNbinsX(); l++ ) { // bin 数を自動でとれるようにしたい
	double signal_nomi = 0.;
	observable = hsignal->GetXaxis()->GetBinCenter(l+1); // hsignal から自動でとってくる
	
	for ( int j = 0; j < nbin / 2; j++ ) {
	  signal_nomi += hsignal->GetBinContent(l + 1, j + 1);
	}
	
	for ( int i = 0; i < nthrows; i++ ) {
	  signal = 0.;
	  for ( int k = 0; k < nbin; k++ )
	    nrand[k] = rand.Gaus();
	  
	  for ( int j = 0; j < nbin; j++ ) {
	    weight[j] = 1.;
	    for ( int k = 0; k <= j; k++ )
	      weight[j] += chol_mat[j][k] * nrand[k];
	    if ( j < nbin / 2 )
	      signal += hsignal->GetBinContent(l+1, j+1) * weight[j];
	  }
	  
	  if ( signal_nomi > 0 )
	    signal /= signal_nomi;
	  else 
	    signal = 0.;
	  
	  tree->Fill();
	}
      }
      
      fout->cd();
      hsignal->Write("hsignal");
      tree->Write();
      fout->Close();
    }

    delete f;
    
  }

}
