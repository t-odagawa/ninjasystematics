void CheckFadingErr() {

  double sigma_p = 50.;
  double sigma_pi = 10.;
  double mean_p = 200.;
  double mean_pi =30.;

  double nominal_thr = 0.;

  for ( double vph = mean_pi; vph < mean_p; vph += 1.) {
    double l_p = TMath::Gaus(vph, mean_p, sigma_p, kTRUE);
    double l_pi = TMath::Gaus(vph, mean_pi, sigma_pi, kTRUE);

    if ( l_p >= l_pi ) {
      nominal_thr = vph;
      break;
    }
  }

  cout << nominal_thr << endl;

  TGraph *g = new TGraph();
  int ipoint = 0;

  for ( double scale = 0.9; scale < 1.1; scale += 0.005 ) {

    double sigma_p_scale = sigma_p * scale;
    double sigma_pi_scale = sigma_pi * scale;

    double thr = 0.;

    for ( double vph = mean_pi; vph < mean_p; vph += 0.1) {
      double l_p = TMath::Gaus(vph, mean_p, sigma_p_scale, kTRUE);
      double l_pi = TMath::Gaus(vph, mean_pi, sigma_pi_scale, kTRUE);

      if ( l_p >= l_pi ) {
	thr = vph;
	break;
      }
    }

    g->SetPoint(ipoint, scale, thr);
    ipoint++;

  }
  
  g->SetMarkerStyle(20);
  g->Draw("AP");


}
