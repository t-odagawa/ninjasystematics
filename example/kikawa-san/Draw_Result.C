float result=0.6305;
float result_qe=0.307;
float bgrat=0.4665;

float GetError(TH1F *h1, bool posi){
  float width = h1->GetBinWidth(1);
  int nbint = h1->GetNbinsX();
  int nbin = nbint/2;
  double total=0;
  for(int i=0;i<nbin;i++){
    if(posi)total+=h1->GetBinContent(nbin+i+1);
    else    total+=h1->GetBinContent(nbin-i);
  }
  double sum=0;
  float err;
  for(int i=0;i<nbin;i++){
    if(posi)sum+=h1->GetBinContent(nbin+i+1);
    else    sum+=h1->GetBinContent(nbin-i);
    if(sum/total>0.682){
      err=i*width;
      break;
    }
  }

  return err;
};


void Draw_Result(){

  TFile *_file0 = TFile::Open("toymc.root");

  float bg_ratio=1-((qedet->GetSumOfWeights())/(qedet->GetSumOfWeights()+nqedet->GetSumOfWeights()));
  bg_ratio=bgrat;

  TCanvas *c1 = new TCanvas("c1","c1",0,0,1300,900);
  gStyle->SetOptTitle(0);

  tree->Draw(Form("((%e-%e*nqedet)/(%e-%e)/flux/qeeff-1)*100>>h1(4000,-200,200)",result, bg_ratio, result, bg_ratio));
  h1->GetXaxis()->SetTitle("Variation of #sigma_{CC-coh.#pi} [%]");
  h1->SetFillColor(2);
  gStyle->SetOptStat(1100);
  float err=h1->GetRMS();

  cout<<err<<endl;

  cout<<GetError(h1,0)<<" "<<GetError(h1,1)<<endl;
}
