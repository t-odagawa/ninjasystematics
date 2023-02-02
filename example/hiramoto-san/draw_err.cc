void draw_err(){

  TH1D*h = new TH1D("muang_err","",18,0,90);

  int n=tree->GetEntries();
  int ang;
  float signal;

  tree->SetBranchAddress("signal",&signal);
  tree->SetBranchAddress("ang",&ang);
  
  for(int i=0;i<18;i++){
    TH1D*tmp = new TH1D("tmp","",100,0,2);
    for(int j=0;j<n;j++){
      tree->GetEntry(j);
      if(ang==i*5)tmp->Fill(signal);
    }
    tmp->Fit("gaus");
    double err=gaus->GetParameter(2);
    if(i>10)err=0;
    h->Fill(i*5,err);
  }
  h->GetYaxis()->SetRangeUser(0,0.2);
  h->GetXaxis()->SetTitle("Muon angle (deg)");
  h->GetYaxis()->SetTitle("Fractional error");
  h->SetStats(0);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->Draw();
}
