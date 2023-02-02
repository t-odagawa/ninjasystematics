#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <TMatrixDSym.h>

#include <B2Enum.hh>

namespace fs = boost::filesystem;
namespace logging = boost::log;

int main(int argc, char* argv[]) {

  if ( argc != 6 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0]
			     << " <nominal file name> <xsec syst momch dir name> <dial name> <hist name> <output file name>";
    return 1;
  }

  TString nominal_filename = argv[1];
  fs::path nominal_path(argv[1]);
  std::cout << nominal_filename << std::endl;
  TFile *nominal_file = new TFile(nominal_filename, "read");
  TString histogram_name = argv[4];
  auto hist_nom = (TH1D*)nominal_file->Get(histogram_name);

  std::string sys_dir_name = std::string(argv[2]);
  std::string dial_name = std::string(argv[3]);

  std::stringstream ss_varied_filename;  
  ss_varied_filename << sys_dir_name
		     << dial_name + "/";
  TString plus_filename = ss_varied_filename.str() + "plus/output/" + nominal_path.leaf().string();
  std::cout << plus_filename << std::endl;
  TFile *plus_file = new TFile(plus_filename, "read");
  auto hist_plus = (TH1D*)plus_file->Get(histogram_name);
  TString minus_filename;
  if ( dial_name == "/SF_OptPotTwkDial_O16" ||
       dial_name == "/PionFSI_InelProb" )
    minus_filename = ss_varied_filename.str() + "plus/output/" + nominal_path.leaf().string();
  else 
    minus_filename = ss_varied_filename.str() + "minus/output/" + nominal_path.leaf().string();
  std::cout << minus_filename << std::endl;
  TFile *minus_file = new TFile(minus_filename, "minus");
  auto hist_minus = (TH1D*)minus_file->Get(histogram_name);

  TString ofilename = argv[5];
  TFile *ofile = new TFile(ofilename, "recreate");
  TMatrixDSym *mat = new TMatrixDSym(hist_nom->GetNbinsX());

  for ( int ibin = 1; ibin <= hist_nom->GetNbinsX(); ibin++ ) {
    double Ni = hist_nom->GetBinContent(ibin);
    double N_plus_i = hist_plus->GetBinContent(ibin);
    double N_minus_i = hist_minus->GetBinContent(ibin);
    for ( int jbin = 1; jbin <= hist_nom->GetNbinsX(); jbin++ ) {
      double Nj = hist_nom->GetBinContent(jbin);
      double N_plus_j = hist_plus->GetBinContent(jbin);
      double N_minus_j = hist_minus->GetBinContent(jbin);

      if (std::fabs(Ni) < 1.e-11 ||
	  std::fabs(Nj) < 1.e-11 ) {
	(*mat)(ibin-1, jbin-1) = 0.;
      }
      else {
	(*mat)(ibin-1, jbin-1) = 0.5 * (Ni - N_plus_i) * (Nj - N_plus_j) / Ni / Nj
	  + 0.5 * (Ni - N_minus_i) * (Nj - N_minus_j) / Ni / Nj;
      }
    }
  }

  
  ofile->cd();
  hist_nom->Write("hist_nom");
  hist_plus->Write("hist_plus");
  hist_minus->Write("hist_minus");
  mat->Write("cov_mat");
  ofile->Close();

  return 0;

}

