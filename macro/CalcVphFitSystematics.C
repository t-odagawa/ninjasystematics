void CalcVphFitSystematics() {

  const double fixed_dedx_value = 1.64;

  TString filename = "/home/t2k/odagawa/NinjaMomentumRecon/data/pid/vph_func_param.txt";
  std::ifstream ifs(filename);

  double input_ang_min, input_ang_max;
  double mean_slope, mean_slope_err;
  double mean_inter, mean_inter_err;
  double sigma_inter, sigma_inter_err;
  double sigma_scale, sigma_scale_err;

  double nominal_mu;
  double nominal_sigma;

  while ( ifs >> input_ang_min >> input_ang_max
	  >> mean_slope >> mean_slope_err
	  >> mean_inter >> mean_inter_err
	  >> sigma_inter >> sigma_inter_err
	  >> sigma_scale >> sigma_scale_err ) {
    
    std::cout << "Angle : (" << input_ang_min << ", " << input_ang_max << ")" << std::endl;
    std::cout << "mu (pbeta = 1200MeV/c) = " << mean_slope << " * " << fixed_dedx_value << " + " << mean_inter << std::endl;
    nominal_mu = mean_slope * fixed_dedx_value + mean_inter;
    std::cout << "mu (pbeta = 1200MeV/c) = " << nominal_mu << std::endl;
    std::cout << "b = " << nominal_mu << " - a * " << fixed_dedx_value << std::endl;
    std::cout << "mu = " << "a * (dE/dx - " << fixed_dedx_value << ") + " << nominal_mu << std::endl;
    std::cout << "sigma_a = " << mean_slope_err << std::endl;
    std::cout << "sigma_b = " << mean_inter_err << std::endl;
    std::cout << "sigma_b (a correlated) = " << -1.64 * mean_slope_err << std::endl;

    std::vector<double> dedx_vec;
    std::vector<double> vph, vph_low, vph_high;
    for ( double dedx = 1.; dedx <= 10.; dedx += 0.1 ) {
      dedx_vec.push_back(dedx);
      vph.push_back(mean_slope * (dedx - fixed_dedx_value) + nominal_mu);
      vph_high.push_back((mean_slope + mean_slope_err) * (dedx - fixed_dedx_value) + nominal_mu);
      vph_low.push_back((mean_slope - mean_slope_err) * (dedx - fixed_dedx_value) + nominal_mu);
    }
    for ( int i = 0; i < vph.size(); i++ ) {
      vph_high.at(i) /= vph.at(i);
      vph_low.at(i) /= vph.at(i);
    }

    TGraph *g_high = new TGraph(vph.size(), &dedx_vec[0], &vph_high[0]);
    TGraph *g_low = new TGraph(vph.size(), &dedx_vec[0], &vph_low[0]);

    TCanvas *c = new TCanvas("c", "c");
    c->DrawFrame(0,0.8,10,1.2);

    g_high->SetMarkerColor(kRed);
    g_high->SetMarkerStyle(20);
    g_high->Draw("P");
    g_low->SetMarkerColor(kBlue);
    g_low->SetMarkerStyle(20);
    g_low->Draw("P SAME");

    std::cout << "sigma (pbeta = 1200MeV/c) = " << sigma_scale << " * sqrt(" << nominal_mu << ") + " << sigma_inter << std::endl;
    nominal_sigma = sigma_scale * std::sqrt(nominal_mu) + sigma_inter;
    std::cout << "sigma (pbeta = 1200MeV/c) = " << nominal_sigma << std::endl;
    std::cout << "b = " << nominal_sigma << " - a * " << std::sqrt(nominal_mu) << std::endl;
    std::cout << "sigma = " << "a * (sqrt(VPH) - " << std::sqrt(nominal_mu) << ") + " << nominal_sigma << std::endl;
    std::cout << "sigma_a = " << sigma_scale_err << std::endl;
    std::cout << "sigma_b = " << sigma_inter_err << std::endl;
    std::cout << "sigma_b (a correlated) = " << -std::sqrt(nominal_mu) * sigma_scale_err << std::endl;
  }

}
