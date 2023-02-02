{

  // const char* plotdir = "/hsm/nu/wagasci/yasuken/wagasci-babymind/plot/20211226";
  const char* plotdir = "~/data/plots/20220822_flux_err";
  //Int_t eachbin = 320;
  Int_t eachbin = 80;
  //TFile *fmult = new TFile("flux_covariance_wagasci_2021_total_bv21v2.root");
  TFile *fmult = new TFile("../yasutome-san/flux_covariance_wagasci_2021_total_bv21v2_reduced_bins.root");

   const int nerrors=21;
   const int ntots=1; // number of "total" fluxes
   TMatrixD *err_cov[nerrors];
   err_cov[1] = (TMatrixD*)fmult->Get("meson_multiplicity_cov");
   err_cov[2] = (TMatrixD*)fmult->Get("pion_rescatter_cov");
   err_cov[3] = (TMatrixD*)fmult->Get("secondary_nucleon_cov");
   err_cov[4] = (TMatrixD*)fmult->Get("interaction_length_cov");

   err_cov[5] = (TMatrixD*)fmult->Get("pbeam_profile_cov");
   err_cov[6] = (TMatrixD*)fmult->Get("horn_current_cov");
   err_cov[7] = (TMatrixD*)fmult->Get("horn_field_asymm_cov");
   err_cov[8] = (TMatrixD*)fmult->Get("horn_misalign_cov");
   err_cov[9] = (TMatrixD*)fmult->Get("target_misalign_cov");
   err_cov[10] = (TMatrixD*)fmult->Get("simulation_materials_cov");
   err_cov[11] = (TMatrixD*)fmult->Get("proton_numb_cov");
   err_cov[12] = (TMatrixD*)fmult->Get("offaxis_angle_cov");
  
   err_cov[0] = (TMatrixD*)fmult->Get("total_flux_cov");//total
   err_cov[13] = new TMatrixD(eachbin,eachbin);//hadron
   err_cov[14] = new TMatrixD(eachbin,eachbin);//beam
   err_cov[15] = new TMatrixD(eachbin,eachbin);//ProtonBeam+Offaxis
   err_cov[16] = new TMatrixD(eachbin,eachbin);//HornCurrent+Field
   err_cov[17] = new TMatrixD(eachbin,eachbin);//Horn+TarGetAlign
   err_cov[18] = new TMatrixD(eachbin,eachbin);//Material
   err_cov[19] = new TMatrixD(eachbin,eachbin);//NumProton
   err_cov[20] = (TMatrixD*)fmult->Get("total_flux_cov");//Total
   
   for(int i=0; i<eachbin; i++) 
     {
     for(int j=0; j<eachbin; j++)
       {
	 (*err_cov[13])(i,j) 
	   = (*err_cov[1])(i,j)+(*err_cov[2])(i,j)
	   +(*err_cov[3])(i,j)+(*err_cov[4])(i,j);
       }
     }

   for(int i=0; i<eachbin; i++) 
     {
     for(int j=0; j<eachbin; j++)
       {
	 (*err_cov[14])(i,j) 
	 = (*err_cov[5])(i,j)+(*err_cov[6])(i,j)
	 +(*err_cov[7])(i,j)+(*err_cov[8])(i,j)
	 +(*err_cov[9])(i,j)+(*err_cov[10])(i,j)
	 +(*err_cov[11])(i,j)+(*err_cov[12])(i,j);
       }
     }

   for(int i=0; i<eachbin; i++) 
     {
     for(int j=0; j<eachbin; j++)
       {
	 (*err_cov[15])(i,j) 
	   = (*err_cov[5])(i,j)+(*err_cov[12])(i,j);
       }
     }

   for(int i=0; i<eachbin; i++) 
     {
     for(int j=0; j<eachbin; j++)
       {
	 (*err_cov[16])(i,j) 
	   = (*err_cov[6])(i,j)+(*err_cov[7])(i,j);
       }
     }

   for(int i=0; i<eachbin; i++) 
     {
     for(int j=0; j<eachbin; j++)
       {
	 (*err_cov[17])(i,j) 
	   = (*err_cov[8])(i,j) + +(*err_cov[9])(i,j);
       }
     }

   for(int i=0; i<eachbin; i++) 
     {
     for(int j=0; j<eachbin; j++)
       {
	 (*err_cov[18])(i,j) 
	   = (*err_cov[10])(i,j);
       }
     }

   for(int i=0; i<eachbin; i++) 
     {
     for(int j=0; j<eachbin; j++)
       {
	 (*err_cov[19])(i,j) 
	   = (*err_cov[11])(i,j);
       }
     }
 
   TH1D *err_hist[nerrors][2][2][4];

   TH1D *flux_hist[2][2][4];
   double flux_rebins[21] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0,
			     1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 30.0};
   double flux_reweights[20] = {1., 1., 1., 1., 1., 1., 1., 1.,
				1/2.,1/2.,
				1/3.,
				1/5., 1/5., 1/5., 1/5., 1/5.,
				1/10., 1/20., 1/30., 1/200.};


   TString nomdir = "/hsm/nu/wagasci/yasuken/flux_uncertainty/bv21v2/3-tuned-p250/";
   TString nomfile = nomdir + "tuned_bv21v2.root";

   TFile *fflux[2][2];
   fflux[0][0] = new TFile(nomfile);
   fflux[1][0] = new TFile(nomfile);
   fflux[0][1] = new TFile(nomfile);
   fflux[1][1] = new TFile(nomfile);
   
   int colors[21] = {2,3,4,5,12,6,7,8,9,15,11,13,14, 2, 3, 4, 5, 6, 8, 9, 1};
   int styles[21] = {1,1,1,1,1,1,1,1,1,1,1,1,1,      1, 1, 1, 1, 1, 1, 1, 1};

   char err_name[21][80] = {"total","meson_multiplicity", "pion_rescatter", "secondary_nucleon","interaction_length","pbeam_profile","horn_current","horn_field_asymm","horn_misalign","target_misalign","simulation_materials","proton_numb","offaxis_angle","hadron","non-hadron","Proton Beam Profile + Off-axis Angle", "Horn Current and Field", "Horn and Target Alignment", "Material Modeling", "Number of Protons", "Total"};

   char err_namet[21][80] = {"Total","meson_multiplicity","pion_rescater","secondary_nucleon","interaction_length","pbeam_profile","horn_current","horn_field_asymmetry","horn_misalignment","target_misalignment","simulation_materials","proton number","Offaxis","Hadron Interaction","Non-Hadron (Beamline)","Proton Beam Profile + Off-axis Angle", "Horn Current and Field", "Horn and Target Alignment", "Material Modeling", "Number of Protons", "Total"};

   char det_name[2][30] = {"nd7","sk"};
   char mode_name[1][30] = {"numode"};
   char flavor_name[4][30] = {"numu","numub","nue","nueb"};

   char det_namet[2][30] = {"WAGASCI", "SK"};
   char mode_namet[1][100] = {"Neutrino Mode"};
   char flavor_namet[4][30] = {"#nu_{#mu}","#bar{#nu}_{#mu}","#nu_{e}", "#bar(#nu)_{e}"};
   
   for(int i=0; i<2; i++)
     for(int j=0; j<1; j++)
       for(int k=0; k<4; k++){
          flux_hist[i][j][k] = (TH1D*)fflux[i][j]->Get(Form("%s_tune_%s",det_name[i],flavor_name[k])); // says numode but depends on file
          for(int l=1; l<=flux_hist[i][j][k]->GetNbinsX(); l++)
             flux_hist[i][j][k]->SetBinContent(l,flux_hist[i][j][k]->GetBinContent(l)*flux_hist[i][j][k]->GetXaxis()->GetBinCenter(l));
	  flux_hist[i][j][k] = (TH1D*)flux_hist[i][j][k]->Rebin(20, "", flux_rebins);
	  for ( int bin = 0; bin < 20; bin++ )
	    flux_hist[i][j][k]->SetBinContent(bin+1, flux_hist[i][j][k]->GetBinContent(bin+1)*flux_reweights[bin]);
          flux_hist[i][j][k]->Scale(0.18/flux_hist[i][j][k]->GetMaximum());
        }


   TCanvas *canv[2][2][4];

   int nbins = 20;
   double bins[21] = {0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0,
                      1.2 , 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 30.0};

   double max[8] = {0.35,0.35,0.35,0.45,0.35,0.35,0.35,0.4};

   TLegend *leg[2][2][4];
   TLegend *leg2[2][2][4];
   TLegend *leg3[2][2][4];
   TLegend *legU[2][2][4];

   for(int id=0; id<1; id++)
     for(int im=0; im<1; im++)
       for(int ifl=0; ifl<1; ifl++){
         printf("canv[id=%d][im=%d][ifl=%d]\n", id, im, ifl);
         canv[id][im][ifl] = new TCanvas(Form("c%d%d%d",id,im,ifl),Form("c%d%d%d",id,im,ifl),700,500);
         gPad->SetLogx();
         gPad->SetRightMargin(0.05);
         leg[id][im][ifl] = new TLegend(0.12,0.63,0.62,0.87);
         //leg2[id][im][ifl] = new TLegend(0.60,0.59,0.92,0.87);
         //leg3[id][im][ifl] = new TLegend(0.15,0.59,0.65,0.63);
         legU[id][im][ifl] = new TLegend(0.62,0.59,0.92,0.87);
         for(int ie=0; ie<nerrors; ie++){
           // if(im==1 && ie==nerrors-1) continue;          
	   if(ie>=1 && ie<=12) continue;
	   //if(ie>=13) continue;
           err_hist[ie][id][im][ifl] = new TH1D(Form("%s_%s_%s_%s",err_name[ie],det_name[id],mode_name[im],flavor_name[ifl]),"",nbins,bins);
           for(int k=1; k<=err_hist[ie][id][im][ifl]->GetNbinsX(); k++){
	     //int entry = (k-1)+id*4+im*2+ifl*1;
             int entry = (k-1)+id*8+im*80+ifl*20;
             /*
             if(ie==nerrors-1) // use different bins for the comparison to old prediction
               entry = (k-1)+id*80+ifl*20;
               */
             err_hist[ie][id][im][ifl]->SetBinContent(k, sqrt( (*(err_cov[ie]))(entry,entry) ));
           } 
           err_hist[ie][id][im][ifl]->SetLineWidth(2);
           err_hist[ie][id][im][ifl]->SetLineColor(colors[ie]);
           err_hist[ie][id][im][ifl]->SetLineStyle(styles[ie]);

           if(ie==0){
            err_hist[ie][id][im][ifl]->SetMaximum(max[ifl+im*4]);
            err_hist[ie][id][im][ifl]->SetMinimum(0.0);
            err_hist[ie][id][im][ifl]->SetStats(false);
            err_hist[ie][id][im][ifl]->GetXaxis()->CenterTitle();
            err_hist[ie][id][im][ifl]->GetYaxis()->CenterTitle();	    
            // err_hist[ie][id][im][ifl]->GetXaxis()->SetTitle("E_{#nu} (GeV)");
            // err_hist[ie][id][im][ifl]->GetYaxis()->SetTitle("Fractional Error");
            err_hist[ie][id][im][ifl]->GetYaxis()->SetNdivisions(505);
            // err_hist[ie][id][im][ifl]->SetTitle(Form("%s: %s, %s",det_namet[id],mode_namet[im],flavor_namet[ifl]));
	    err_hist[ie][id][im][ifl]->SetTitle(";E_{#nu} [GeV];Fractional Error");
            err_hist[ie][id][im][ifl]->Draw("hist");
            flux_hist[id][im][ifl]->SetFillColor(18);
            flux_hist[id][im][ifl]->Draw("hist same");
            //err_hist[ie][id][im][ifl]->Draw("hist same");

            legU[id][im][ifl]->AddEntry(flux_hist[id][im][ifl],"#Phi#timesE_{#nu}, Arb. Norm.","f");
	    //legU[id][im][ifl]->AddEntry(err_hist[ie][id][im][ifl],err_namet[ie],"l");	    
	   
	   } else {
	     err_hist[ie][id][im][ifl]->Draw("hist same");
	     
	     if(ie==14 || ie==18 || ie==19 || ie == 20){
	       legU[id][im][ifl]->AddEntry(err_hist[ie][id][im][ifl],err_namet[ie],"l");
	     }else if(ie==13 || (ie>=15 && ie<=17)){
	       leg[id][im][ifl]->AddEntry(err_hist[ie][id][im][ifl],err_namet[ie],"l");
	     }
	     
	   }
	   	   
	 }
	
        leg[id][im][ifl]->SetFillColor(0);
        leg[id][im][ifl]->SetBorderSize(0);
        leg[id][im][ifl]->Draw("");

        legU[id][im][ifl]->SetFillColor(0);
        legU[id][im][ifl]->SetBorderSize(0);
        legU[id][im][ifl]->Draw("");

        canv[id][im][ifl]->Print(Form("%s/flux_error_t2k_%s_%s_%s.pdf",plotdir,det_name[id],mode_name[im],flavor_name[ifl]));
        canv[id][im][ifl]->Print(Form("%s/flux_error_t2k_%s_%s_%s.png",plotdir,det_name[id],mode_name[im],flavor_name[ifl]));
       }

 
             
}
