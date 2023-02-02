void DrawBMSystematics() {

  TString filename = "../data/bmsyst/detector_systematics_ch.root";
  TFile *file = new TFile(filename, "read");

  auto test = (TMatrixT<float>*)file->Get("cov_matrix_syst_p_7");

  // test->Draw("colz");
  TH1D* hist = new TH1D("hist", "hist", 8, 0, 8);
  for ( int i = 0; i < 8; i++ ) {
    hist->Fill(i + 0.5, (*test)(i,i) * (*test)(i,i));
  }

  hist->Draw("hist");

}
