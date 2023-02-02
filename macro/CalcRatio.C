void CalcRatio() {

  gStyle->SetOptStat(0);

  // std::string systematics = "mcs_scale_syst";
  // std::string systematics = "angres_syst";
  // std::string systematics = "bt_eff_syst";
  std::string systematics = "md_syst";
  // std::string systematics = "maccqe_syst";

  TString ofilename = "/home/t2k/odagawa/NinjaSystematics/data/ratio/ratio_hist_" + systematics + ".root";
  TFile *ofile = new TFile(ofilename, "recreate");

  for ( int mode = 0; mode < 4; mode++ ) {

    std::vector<TString> hist_vec;

    if ( mode == 0 ) {
      hist_vec.push_back("hist_water_mode_multi_");
      hist_vec.push_back("hist_water_mode_proton_multi_");
      hist_vec.push_back("hist_water_mode_pion_multi_");
    }
    else if ( mode == 1 ) continue;
    else if ( mode == 2 ) {
      hist_vec.push_back("hist_muon_mom_");
      hist_vec.push_back("hist_muon_mom_range_");
      hist_vec.push_back("hist_muon_mom_mcs_");
      hist_vec.push_back("hist_proton_mom_");
      hist_vec.push_back("hist_proton_mom_range_");
      hist_vec.push_back("hist_proton_mom_mcs_");
      hist_vec.push_back("hist_pion_mom_");
      hist_vec.push_back("hist_pion_mom_range_");
      hist_vec.push_back("hist_pion_mom_mcs_");
    }
    else if ( mode == 3 ) {
      hist_vec.push_back("hist_muon_ang_deg_");
      hist_vec.push_back("hist_proton_ang_deg_");
      hist_vec.push_back("hist_pion_ang_deg_");
    }

    std::stringstream filename_nom_ss;
    filename_nom_ss << "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output/output_mode" << mode << ".root";
    std::stringstream filename_plus_ss;
    filename_plus_ss << "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output_" << systematics << "/plus/output_mode" << mode << ".root";
    std::stringstream filename_minus_ss;
    filename_minus_ss << "/hsm/nu/ninja/pra_tmp/mc_tmp_20220620/output_" << systematics << "/minus/output_mode" << mode << ".root";
    TString filename_nom = filename_nom_ss.str();
    TString filename_plus = filename_plus_ss.str();
    TString filename_minus = filename_minus_ss.str();
    
    TFile *file_nom = new TFile(filename_nom, "read");
    TFile *file_plus = new TFile(filename_plus, "read");
    TFile *file_minus = new TFile(filename_minus, "read");

    for ( auto hist_name : hist_vec ) {

      TH1D *hist = (TH1D*)file_nom->Get(Form(hist_name + "%d", 0));

      const int nbins = hist->GetNbinsX();
      double bins[nbins+1];
      for ( int ibin = 0; ibin < nbins; ibin++ ) {
	bins[ibin] = hist->GetBinLowEdge(ibin+1);
      }
      bins[nbins] = hist->GetBinLowEdge(nbins) + hist->GetBinWidth(nbins);
      
      TH1D *hist_nom = new TH1D("hist_nom", "", nbins, bins);
      TH1D *hist_plus = new TH1D("hist_plus", "", nbins, bins);
      TH1D *hist_minus = new TH1D("hist_minus", "", nbins, bins);
      
      TList *list_nom = new TList();
      TList *list_plus = new TList();
      TList *list_minus = new TList();
      
      for ( int i = 0; i < 5; i++ ) {
	list_nom->Add((TH1D*)file_nom->Get(Form(hist_name + "%d", i)));
	list_plus->Add((TH1D*)file_plus->Get(Form(hist_name + "%d", i)));
	list_minus->Add((TH1D*)file_minus->Get(Form(hist_name + "%d", i)));
      }
      
      hist_nom->Merge(list_nom);
      hist_plus->Merge(list_plus);
      hist_minus->Merge(list_minus);
      hist_plus->Divide(hist_nom);
      hist_minus->Divide(hist_nom);
      
      hist_plus->SetLineWidth(2);
      hist_minus->SetLineWidth(2);
      hist_plus->SetLineColor(kRed);
      hist_minus->SetLineColor(kBlue);
      
      TLegend *l = new TLegend(0.15,0.7,0.3,0.89);
      l->AddEntry(hist_plus, "+1#sigma", "l");
      l->AddEntry(hist_minus, "-1#sigma", "l");
      
      TCanvas *c = new TCanvas();
      c->SetName("c" + hist_name);
      
      auto frame = c->DrawFrame(bins[0], 0.9, bins[nbins], 1.1);
      frame->SetTitle(";;Ratio");
      
      
      hist_plus->Draw("hist same"); 
      hist_minus->Draw("hist same");
      
      gPad->SetGrid();
      
      l->Draw();

      ofile->cd();
      c->Write();
      
      std::stringstream ss;
      ss << "/home/t2k/odagawa/NinjaSystematics/data/ratio/" << systematics << "_" << hist_name << "ratio.txt";
      std::ofstream ofs(ss.str());
      
      for ( int i = 0; i < nbins; i++ ) {
	ofs << bins[i] << " " << bins[i+1] << " "
	    << hist_plus->GetBinContent(i+1) << " "
	    << hist_minus->GetBinContent(i+1) << " "
	    << 1 - hist_plus->GetBinContent(i+1) << " "
	    << 1 - hist_minus->GetBinContent(i+1) << std::endl;
      }
      
      ofs.close();

      hist_nom->Delete();
      hist_plus->Delete();
      hist_minus->Delete();
      
    }
    
  }

  ofile->Close();

}
