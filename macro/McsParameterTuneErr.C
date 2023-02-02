void McsParameterTuneErr() {

  TFile *file = new TFile("muon_wo_s_tune_graph.root");
  TGraphErrors *ge = (TGraphErrors*)file->Get("ge");

  TF1 *line = new TF1("line", "[0]", 200, 2000);

  TGraph *g = new TGraph();
  TGraph *g1 = new TGraph();
  int ipoint = 0;
  
  for ( int i = 2; i < 7; i++ ) {
    for ( int j = 16; j < 21; j++ ) {

      double min = i * 100;
      double max = j * 100;

      ge->Fit(line,"Q","",min, max);

      g->SetPoint(ipoint, ipoint, line->GetParameter(0));
      if ( min == 500. && max == 2000. ) g1->SetPoint(0, ipoint, line->GetParameter(0));
      ipoint++;
      
    }
  }

  g->SetMarkerStyle(20);
  g->Draw("AP");
  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(kRed);
  g1->Draw("SAME P");
  
}
