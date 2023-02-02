#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <TFile.h>

#include <B2Enum.hh>

#include <WeightTree.hpp>

#include <McsClass.hpp>

namespace logging = boost::log;

int main(int argc, char* argv[]) {

  if ( argc != 6 ) {
    BOOST_LOG_TRIVIAL(error) << "Usage : " << argv[0] << " <input file> <reweight file> <dial> <sigma> <output file>";
    return 1;
  }

  logging::core::get()->set_filter
    (
     logging::trivial::severity >= logging::trivial::info
     );

  try {

    std::string imomchname = std::string(argv[1]);
    BOOST_LOG_TRIVIAL(info) << "Input file name : " << imomchname;
    auto ev_vec = Momentum_recon::ReadEventInformationBin(imomchname);
    BOOST_LOG_TRIVIAL(info) << "Event size : " << ev_vec.size();
    
    TFile rwfile(argv[2], "read");
    wgrw::WeightTree rwtree(rwfile, wgrw::TreeMode::kReadWeights);

    BOOST_LOG_TRIVIAL(info) << "ReWeight file name : " << argv[2];
    // BOOST_LOG_TRIVIAL(debug) << "ReWeight tree size : " << rwtree.GetEntries();

    std::string dial_name = (std::string)argv[3];
    std::string dial_sigma = (std::string)argv[4];
    const auto &dial = rwtree.GetDialList().At(dial_name);

    BOOST_LOG_TRIVIAL(info) << "Creating weight event map...";
    
    for ( auto &ev : ev_vec ) {
      
      int entry = ev.entry_in_daily_file;
      BOOST_LOG_TRIVIAL(debug) << "Groupid : " << ev.groupid << ", "
			       << "Entry : " << entry;

      rwtree.ReadWeightEntry(entry);
      const auto &rw = rwtree.GetWeightEvent();

      double rwfactor = 1.;

      if ( dial_sigma == "plus" ) {
	rwfactor = rw.GetWeight(dial, B2Sigma::kPlusOneSigma);
	BOOST_LOG_TRIVIAL(debug) << "ReWeight factor = " << rwfactor << "(+1sigma)";
      }
      else if ( dial_sigma == "minus" ) {
	rwfactor = rw.GetWeight(dial, B2Sigma::kMinusOneSigma);
	BOOST_LOG_TRIVIAL(debug) << "ReWeight factor = " << rwfactor << "(-1sigma)";
      }

      ev.weight *= rwfactor;
      
    }

    std::string ofilename = argv[5];
    Momentum_recon::WriteEventInformationBin(ofilename, ev_vec);

  } catch (const std::runtime_error &error) {
    BOOST_LOG_TRIVIAL(fatal) << "Runtime error : " << error.what();
    std::exit(1);
  }

  return 0;

}
