//#define MKPLOT
//#define MKPLOT2
#define nbin 40
TRandom3 rand;

void draw_matrix(float mat[nbin][nbin]){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  TH2F *h1=new TH2F("h1","h1",nbin,0,nbin,nbin,0,nbin);
  h1->GetZaxis()->SetRangeUser(-0.045,0.045);
  h1->Reset();
  for(int i=0;i<nbin;i++){
    for(int j=0;j<nbin;j++){
      h1->Fill(i+0.0001,j+0.0001,mat[i][j]);
    }
  }
  h1->Draw("colz");
  c1->Update();

  char ans[8];
  fgets(ans,8,stdin);

};

void cholcov_conv(float covmat[nbin][nbin], float cholcovmat[nbin][nbin])
{
  memset(cholcovmat,0,sizeof(cholcovmat));
  for ( Int_t j=0; j<nbin; j++ ) { 
    Double_t s = covmat[j][j] ;
    for ( Int_t k=0; k<j; k++ ) {
      s -= cholcovmat[j][k]*cholcovmat[j][k] ;
    }
    if(s<0){
      std::cout << "strange !" << j << " " << s << std::endl ;
      exit(0) ;
    }
    cholcovmat[j][j] = sqrt(s) ;

    for ( Int_t i=j+1; i<nbin; i++ ) {
      s = covmat[i][j] ;
      for ( Int_t k=0; k<j; k++ ) {
        s -= cholcovmat[i][k]*cholcovmat[j][k] ;
      }
      if ( TMath::Abs(s)<0.000000000001 )
        cholcovmat[i][j] = 0. ;
      else
        cholcovmat[i][j] = s/cholcovmat[j][j] ;
      //cout<<cholcovmat[i][j]<<endl;
    }
  }

};


void FluxErr_all(){

  TFile *fop = TFile::Open("covariance_tn336.root");
  TH2D *hcov = (TH2D*)fop->Get("hcov");
  if(hcov==NULL){
    cout << "Cannot get matrix from " << fop->GetPath() << endl;
    exit;
  }

  float cov_mat[nbin][nbin];
  for(int i=0; i<nbin; i++){
    for(int j=0; j<nbin; j++){
      cov_mat[i][j] = hcov->GetBinContent(hcov->GetBin(i+481,j+481));
      cout<<i<<" "<<j<<" "<<cov_mat[i][j]<<endl;
    }
  }
#ifdef MKPLOT
  draw_matrix(cov_mat);
#endif

  float chol_mat[nbin][nbin];
  cholcov_conv(cov_mat, chol_mat);

#ifdef MKPLOT2
  draw_matrix(chol_mat);
#endif

  int nthrows = 100000;
  float weight[nbin],nrand[nbin];
  TH1F*  hsignal;

  TFile *f  = new TFile("muang.root");

  hsignal  = (TH1F*)f  -> Get("muang");
  hsignal_numu  = (TH1F*)f  -> Get("muang_numu");
  hsignal_numub  = (TH1F*)f  -> Get("muang_numubar");

  TFile *fout = new TFile("toymc.root","recreate");
  TTree *tree = new TTree("tree","tree");

  float signal, signal_numu, signal_numub;
  int ang;

  tree->Branch("signal",&signal,"signal/F");
  //tree->Branch("signal_numu",&signal_numu,"signal_numu/F");
  //tree->Branch("signal_numub",&signal_numub,"signal_numub/F");
  tree->Branch("ang",&ang,"ang/I");


  for(int l=0;l<18;l++){
    float signal_nomi=0;
    ang = 5*l;
    
    for(int j=0;j<nbin/2;j++){
      signal_nomi   += hsignal_numu   -> GetBinContent(l+1,j+1);
      signal_nomi   += hsignal_numub   -> GetBinContent(l+1,j+1);
    }
    
    
    for(int i=0;i<nthrows;i++){
      if(i%10000==0)cout<<ang<<" "<<i<<endl;
      
      signal=0;
      signal_numu=0;
      signal_numub=0;
      
      for(int k=0;k<nbin;k++){
	nrand[k]=rand.Gaus();
      }
    
      for(int j=0;j<nbin;j++){
	weight[j]=1;
	for(int k=0;k<=j;k++){
	  weight[j]+=chol_mat[j][k]*nrand[k];
 	}
	if(j<nbin/2){
	  signal_numu   += hsignal_numu   -> GetBinContent(l+1,j+1)*weight[j];
	}
	else{
	  signal_numub   += hsignal_numub   -> GetBinContent(l+1,j+1-nbin/2)*weight[j];
	}
      }
      
      signal = (signal_numu+signal_numub)/signal_nomi;
      
      tree->Fill();
    }
  }
  
  tree->Write();
  
  hsignal  ->Write("signal");
  
  fout->Close();
}
